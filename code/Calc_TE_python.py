#Calc_TE_python.py
import os
import pandas as pd
import argparse
import numpy as np
import multiprocessing
import time
from scipy import sparse
from multiprocessing import Process, Array
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

from scipy.special import digamma
from sklearn.neighbors import KDTree  # Fallback KDTree
from scipy.spatial import cKDTree as CKDTree  # Faster C-backed KDTree (supports p=inf)
from dataclasses import dataclass
from typing import List, Tuple

# Small utility for fast window views if available
try:
    from numpy.lib.stride_tricks import sliding_window_view as _sliding_window_view
except Exception:  # numpy < 1.20 or unavailable
    _sliding_window_view = None

# Optional: Numba-accelerated counting if available
try:
    from numba import njit
except Exception:
    njit = None

def _neighbors_to_csr(neighbors: List[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    """Convert list-of-arrays neighbors to CSR-like (indptr, indices)."""
    indptr = np.zeros(len(neighbors) + 1, dtype=np.int64)
    total = 0
    for i, arr in enumerate(neighbors):
        total += int(arr.size)
        indptr[i + 1] = total
    indices = np.empty(total, dtype=np.int64)
    pos = 0
    for arr in neighbors:
        n = int(arr.size)
        if n:
            indices[pos:pos + n] = arr
            pos += n
    return indptr, indices

if njit is not None:
    @njit(cache=True)
    def _count_csr(indptr: np.ndarray, indices: np.ndarray, svals: np.ndarray) -> np.ndarray:
        out = np.empty(svals.shape[0], dtype=np.int32)
        for i in range(svals.shape[0]):
            s_i = svals[i]
            cnt = 0
            start = indptr[i]
            end = indptr[i + 1]
            for p in range(start, end):
                j = indices[p]
                if abs(svals[j] - s_i) <= 1.0:
                    cnt += 1
            out[i] = cnt
        return out


@dataclass
class EmpiricalMeasurementDistribution:
    """
    Class to store the empirical distribution of Transfer Entropy under permutation testing.

    Attributes:
        distribution (np.ndarray): Array containing the Transfer Entropy values obtained from permutation testing.
        p_value (float, optional): The p-value of the actual Transfer Entropy, calculated based on the permutation distribution. Defaults to None.
        actual_value (float, optional): The actual Transfer Entropy value calculated from the original data. Defaults to None.
    """
    distribution: np.ndarray
    p_value: float = None
    actual_value: float = None


class KernelEstimatorTransferEntropy:
    """
    Kernel Estimator for Transfer Entropy.
    Handles density estimation and counting within a specified kernel width.
    Optimized for batch processing and vectorized operations.
    Uses KDTree with Chebyshev metric for efficient neighbor queries.

    Attributes:
        normalise (bool): Whether to normalise the data by standard deviation. Default is True.
        kernel_width (float): The kernel width for density estimation. Default is 0.5.
        debug (bool): Enable debug prints. Default is False.
        time_window (int): Window for dynamic correlation exclusion. Default is 0 (disabled).
        force_compare_to_all (bool): Force comparison to all data points (not fully implemented/used). Default is False.
        indices (np.ndarray): Indices of the observations, used for dynamic correlation exclusion.
        kernel_widths (dict): Dictionary storing kernel widths for different dimensions.
        X_past, X_past_source, X_past_next, X_past_source_next (np.ndarray): Scaled observation matrices.
        tree (dict): Dictionary storing KDTree models for different observation sets.
    """
    def __init__(self, normalise=True, kernel_width=0.5, debug=False):
        self.normalise = normalise
        self.kernel_width = kernel_width
        self.debug = debug
        self.time_window = 0  # Dynamic correlation exclusion window
        self.force_compare_to_all = False
        self.indices = None  # To keep track of observation indices for exclusion

    def set_normalise(self, normalise: bool):
        """Set whether to normalise the data."""
        if not isinstance(normalise, bool):
            raise ValueError("normalise must be a boolean.")
        self.normalise = normalise

    def set_kernel_width(self, kernel_width: float):
        """Set the kernel width."""
        if not isinstance(kernel_width, (int, float)) or kernel_width <= 0:
            raise ValueError("kernel_width must be a positive number.")
        self.kernel_width = kernel_width

    def set_dynamic_correlation_exclusion(self, time_window: int):
        """Set the time window for dynamic correlation exclusion."""
        if not isinstance(time_window, int) or time_window < 0:
            raise ValueError("time_window must be a non-negative integer.")
        self.time_window = time_window

    def clear_dynamic_correlation_exclusion(self):
        """Clear the dynamic correlation exclusion window (disable exclusion)."""
        self.time_window = 0

    def set_force_compare_to_all(self, force: bool):
        """Set whether to force comparison to all data points (not fully implemented/used)."""
        if not isinstance(force, bool):
            raise ValueError("force must be a boolean.")
        self.force_compare_to_all = force

    def initialise(self, k: int, kernel_width: float = None):
        """Initialise the estimator with history length k and optional kernel width."""
        if not isinstance(k, int) or k <= 0:
            raise ValueError("k must be a positive integer.")
        self.k = k
        if kernel_width is not None:
            self.set_kernel_width(kernel_width)

    def set_observations(self, dest_past_vectors: np.ndarray, dest_next_values: np.ndarray, source_values: np.ndarray):
        """
        Sets the observations and prepares the data for neighbor searches.
        Scales the data based on standard deviation if normalisation is enabled.
        Builds KDTree models for efficient neighbor queries using Chebyshev metric.

        Args:
            dest_past_vectors (np.ndarray): 2D array of destination past vectors.
            dest_next_values (np.ndarray): 1D array of destination next values.
            source_values (np.ndarray): 1D array of source values.
        """
        if not isinstance(dest_past_vectors, np.ndarray) or dest_past_vectors.ndim != 2:
            raise ValueError("dest_past_vectors must be a 2D numpy array.")
        if not isinstance(dest_next_values, np.ndarray) or dest_next_values.ndim != 1:
            raise ValueError("dest_next_values must be a 1D numpy array.")
        if not isinstance(source_values, np.ndarray) or source_values.ndim != 1:
            raise ValueError("source_values must be a 1D numpy array.")
        if len(dest_past_vectors) != len(dest_next_values) or len(dest_past_vectors) != len(source_values):
            raise ValueError("Input arrays must have the same length.")


        if self.normalise:
            std_dest_past = np.std(dest_past_vectors, axis=0, ddof=1)
            std_source = np.std(source_values, ddof=1)
            std_dest_next = np.std(dest_next_values, ddof=1)

            epsilon = 1e-9 # Small epsilon to avoid division by zero and numerical instability
            # Add epsilon to std to handle cases where std is zero or very close to zero
            kw_dest_past = self.kernel_width * (std_dest_past + epsilon)
            kw_source = self.kernel_width * (std_source + epsilon)
            kw_dest_next = self.kernel_width * (std_dest_next + epsilon)
        else:
            kw_dest_past = np.full(dest_past_vectors.shape[1], self.kernel_width)
            kw_source = self.kernel_width
            kw_dest_next = self.kernel_width

        self.kernel_widths = {
            'dest_past': kw_dest_past,
            'source': kw_source,
            'dest_next': kw_dest_next
        }

        self.X_past, self.X_past_source, self.X_past_next, self.X_past_source_next = self._scale_observations(
            dest_past_vectors, dest_next_values, source_values
        )

        # Build KDTree models for fast radius queries (with Chebyshev metric)
        # Chebyshev metric is used as it's computationally efficient and often suitable for kernel density estimation
        self._fit_kdtrees()

        self.indices = np.arange(len(dest_past_vectors))

    def _scale_observations(self, dest_past, dest_next, source):
        """Scale observations using kernel widths."""
        # Scale observations using kernel widths, handling potential division by zero
        scaled_past = np.divide(dest_past, self.kernel_widths['dest_past'], out=np.zeros_like(dest_past),
                                 where=self.kernel_widths['dest_past'] != 0)
        scaled_source = np.divide(source.reshape(-1, 1), self.kernel_widths['source'], out=np.zeros_like(source.reshape(-1, 1)),
                                   where=self.kernel_widths['source'] != 0)
        scaled_next = np.divide(dest_next.reshape(-1, 1), self.kernel_widths['dest_next'], out=np.zeros_like(dest_next.reshape(-1, 1)),
                                 where=self.kernel_widths['dest_next'] != 0)

        X_past_source = np.hstack((scaled_past, scaled_source))
        X_past_next = np.hstack((scaled_past, scaled_next))
        X_past_source_next = np.hstack((scaled_past, scaled_source, scaled_next))

        return scaled_past, X_past_source, X_past_next, X_past_source_next

    def _fit_kdtrees(self):
        """Build KDTree models for each set of observations."""
        # Build KDTree models for each set of observations using Chebyshev metric
        self.tree = {}
        self.tree['past'] = KDTree(self.X_past, metric='chebyshev') # Consider Euclidean metric as an alternative
        self.tree['past_source'] = KDTree(self.X_past_source, metric='chebyshev')
        self.tree['past_next'] = KDTree(self.X_past_next, metric='chebyshev')
        self.tree['past_source_next'] = KDTree(self.X_past_source_next, metric='chebyshev')

    def get_counts_batch(self, dest_past: np.ndarray, dest_next: np.ndarray, source: np.ndarray, indices: np.ndarray):
        """
        Calculates neighbor counts for a batch of observations.

        Args:
            dest_past (np.ndarray): Batch of destination past vectors.
            dest_next (np.ndarray): Batch of destination next values.
            source (np.ndarray): Batch of source values.
            indices (np.ndarray): Indices of the observations in the batch.

        Returns:
            dict: Dictionary containing neighbor counts for different contexts.
                  Keys: 'countPast', 'countPastSource', 'countNextPast', 'countNextPastSource'
        """
        # Re-scale observations for the batch (if they differ from the ones used for building trees - although they should be the same scaling)
        scaled_past, scaled_past_source, scaled_past_next, scaled_past_source_next = self._scale_observations(dest_past, dest_next, source)

        # For non-dynamic exclusion, use vectorized count queries using KDTree's query_radius
        count_past = self.tree['past'].query_radius(scaled_past, r=1.0, count_only=True)
        count_past_source = self.tree['past_source'].query_radius(scaled_past_source, r=1.0, count_only=True)
        count_past_next = self.tree['past_next'].query_radius(scaled_past_next, r=1.0, count_only=True)
        count_past_source_next = self.tree['past_source_next'].query_radius(scaled_past_source_next, r=1.0, count_only=True)

        # If dynamic correlation exclusion is enabled, adjust counts accordingly
        if self.time_window > 0:
            count_past, count_past_source, count_past_next, count_past_source_next = self._exclude_dynamic_correlations(
                count_past, count_past_source, count_past_next, count_past_source_next,
                indices, scaled_past, scaled_past_source, scaled_past_next, scaled_past_source_next
            )

        if self.debug:
            print(f"Batch Counts: count_past={count_past}, count_past_source={count_past_source}, "
                  f"count_past_next={count_past_next}, count_past_source_next={count_past_source_next}")

        return {
            'countPast': count_past,
            'countPastSource': count_past_source,
            'countNextPast': count_past_next,
            'countNextPastSource': count_past_source_next
        }

    def _exclude_dynamic_correlations(self, count_past, count_past_source, count_past_next, count_past_source_next,
                                      indices, scaled_past, scaled_past_source, scaled_past_next, scaled_past_source_next):
        """
        Excludes dynamically correlated neighbors within the specified time window.
        This is done by adjusting the neighbor counts obtained from KDTree queries.

        Note: This implementation iterates through each observation to exclude neighbors within the time window.
        For very large datasets and/or large time windows, this part could become computationally intensive.
        Profiling and potential optimizations (e.g., modifying KDTree query or alternative data structures)
        might be needed if this becomes a performance bottleneck.
        """
        # Compute exclusion boundaries (note: indices are contiguous in each batch)
        exclusion_start = np.maximum(indices - self.time_window, 0)
        exclusion_end = np.minimum(indices + self.time_window + 1, len(self.indices) if self.indices is not None else len(indices) ) # Use pre-set indices length if available, else batch indices length

        # Batch query: get neighbor indices for all observations at once for each tree
        neighbors_past = self.tree['past'].query_radius(scaled_past, r=1.0, return_distance=False)
        neighbors_past_source = self.tree['past_source'].query_radius(scaled_past_source, r=1.0, return_distance=False)
        neighbors_past_next = self.tree['past_next'].query_radius(scaled_past_next, r=1.0, return_distance=False)
        neighbors_past_source_next = self.tree['past_source_next'].query_radius(scaled_past_source_next, r=1.0, return_distance=False)

        # For each observation, count neighbors that fall into the exclusion range and subtract
        for i in range(len(indices)):
            # Since the exclusion indices form a contiguous range within the batch, we count neighbors that fall within [exclusion_start, exclusion_end)
            count_past[i] -= np.count_nonzero((neighbors_past[i] >= exclusion_start[i]) & (neighbors_past[i] < exclusion_end[i]))
            count_past_source[i] -= np.count_nonzero((neighbors_past_source[i] >= exclusion_start[i]) & (neighbors_past_source[i] < exclusion_end[i]))
            count_past_next[i] -= np.count_nonzero((neighbors_past_next[i] >= exclusion_start[i]) & (neighbors_past_next[i] < exclusion_end[i]))
            count_past_source_next[i] -= np.count_nonzero((neighbors_past_source_next[i] >= exclusion_start[i]) & (neighbors_past_source_next[i] < exclusion_end[i]))

        return count_past, count_past_source, count_past_next, count_past_source_next


class TransferEntropyCalculatorKernel:
    """
    Calculates Transfer Entropy using Kernel Estimation.

    Attributes:
        DEFAULT_KERNEL_WIDTH (float): Default kernel width.
        LOG2 (float): Natural logarithm of 2.
        te_kernel_estimator (KernelEstimatorTransferEntropy): Kernel estimator object.
        normalise (bool): Whether to normalise the data.
        kernel_width (float): Kernel width for estimation.
        k (int): History length.
        dest_past_vectors (np.ndarray): Destination past vectors.
        dest_next_values (np.ndarray): Destination next values.
        source_values (np.ndarray): Source values.
        total_observations (int): Total number of observations.
        debug (bool): Enable debug mode.
        last_average (float): Last computed average Transfer Entropy value.
        vector_of_destination_observations (list): List of destination time series.
        vector_of_source_observations (list): List of source time series.
    """
    DEFAULT_KERNEL_WIDTH = 0.25
    LOG2 = np.log(2)

    def __init__(self):
        self.te_kernel_estimator = KernelEstimatorTransferEntropy()
        self.normalise = True
        self.kernel_width = self.DEFAULT_KERNEL_WIDTH
        self.k = None
        self.dest_past_vectors = None
        self.dest_next_values = None
        self.source_values = None
        self.total_observations = 0
        self.debug = False
        self.last_average = None
        self.vector_of_destination_observations = []
        self.vector_of_source_observations = []

    def set_debug(self, debug: bool):
        """Set debug mode."""
        if not isinstance(debug, bool):
            raise ValueError("debug must be a boolean.")
        self.debug = debug
        self.te_kernel_estimator.debug = debug

    def set_property(self, property_name: str, property_value):
        """Set a specific property of the Transfer Entropy calculator."""
        if not isinstance(property_name, str):
            raise ValueError("property_name must be a string.")
        property_name_upper = property_name.upper()
        if property_name_upper == "KERNEL_WIDTH":
            try:
                self.kernel_width = float(property_value)
                self.te_kernel_estimator.set_kernel_width(self.kernel_width)
            except ValueError:
                raise ValueError("KERNEL_WIDTH must be a float.")
        elif property_name_upper == "NORMALISE":
            if not isinstance(property_value, str):
                raise ValueError("NORMALISE must be a string ('true' or 'false').")
            self.normalise = property_value.lower() == 'true'
            self.te_kernel_estimator.set_normalise(self.normalise)
        elif property_name_upper == "DYN_CORR_EXCL":
            try:
                self.te_kernel_estimator.set_dynamic_correlation_exclusion(int(property_value))
            except ValueError:
                raise ValueError("DYN_CORR_EXCL must be an integer.")
        elif property_name_upper == "FORCE_KERNEL_COMPARE_TO_ALL":
            if not isinstance(property_value, str):
                raise ValueError("FORCE_KERNEL_COMPARE_TO_ALL must be a string ('true' or 'false').")
            self.te_kernel_estimator.set_force_compare_to_all(property_value.lower() == 'true')
        else:
            raise ValueError(f"Unknown property name: {property_name}")


    def initialise(self, k: int, kernel_width: float = None):
        """Initialise the Transfer Entropy calculator with history length k and optional kernel width."""
        if not isinstance(k, int) or k <= 0:
            raise ValueError("k must be a positive integer.")
        self.k = k
        if kernel_width is not None:
            self.set_property("KERNEL_WIDTH", kernel_width) # Use set_property for validation
        self.te_kernel_estimator.initialise(k, self.kernel_width)

    def add_observations(self, source, destination):
        """
        Add source and destination time series observations.

        Args:
            source (list or np.ndarray): Source time series.
            destination (list or np.ndarray): Destination time series.
        """
        # Validate that source and destination are list or np.ndarray
        if not isinstance(source, (list, np.ndarray)):
            raise TypeError("source must be a list or numpy array.")
        if not isinstance(destination, (list, np.ndarray)):
            raise TypeError("destination must be a list or numpy array.")

        # Validate that source and destination have the same length
        if len(source) != len(destination):
            raise ValueError("Source and destination time series must have the same length.")
        self.vector_of_source_observations.append(np.array(source))
        self.vector_of_destination_observations.append(np.array(destination))

    def finalise_add_observations(self):
        """Finalise adding observations and prepare data for TE calculation."""
        # Compute total observations based on history length k
        self.total_observations = sum(len(dest) - self.k for dest in self.vector_of_destination_observations if len(dest) > self.k)
        if self.total_observations == 0:
            raise ValueError("Not enough observations after considering history length k. "
                             "Ensure destination time series are longer than k.")

        self.dest_past_vectors, self.dest_next_values, self.source_values, self.indices = self._prepare_data_for_batch()
        self.te_kernel_estimator.set_observations(self.dest_past_vectors, self.dest_next_values, self.source_values)
        self.vector_of_source_observations = [] # Clear to free memory
        self.vector_of_destination_observations = [] # Clear to free memory

    def _prepare_data_for_batch(self):
        """Prepare data batches from added observations."""
        dest_past_vectors = np.empty((self.total_observations, self.k))
        dest_next_values = np.empty(self.total_observations)
        source_values = np.empty(self.total_observations)
        indices = np.empty(self.total_observations, dtype=int)

        start_observation = 0
        for source, destination in zip(self.vector_of_source_observations, self.vector_of_destination_observations):
            num_observations = len(destination) - self.k
            if num_observations > 0:
                end_observation = start_observation + num_observations
                # Build destination past vectors using list comprehension (vectorized slicing) for efficiency
                dest_past_vectors[start_observation:end_observation] = np.array([destination[i:i+self.k] for i in range(num_observations)])
                dest_next_values[start_observation:end_observation] = destination[self.k:]
                # Align source observations appropriately (shifted by k-1)
                source_values[start_observation:end_observation] = source[self.k - 1:-1]
                indices[start_observation:end_observation] = np.arange(start_observation, end_observation)
                start_observation = end_observation

        return dest_past_vectors, dest_next_values, source_values, indices

    def compute_average_local_of_observations(self):
        """Compute the average local Transfer Entropy for the current observations."""
        counts = self.te_kernel_estimator.get_counts_batch(self.dest_past_vectors, self.dest_next_values, self.source_values, self.indices)
        # Valid indices where counts are positive to avoid division by zero in log
        valid = (counts['countPastSource'] > 0) & (counts['countPast'] > 0) & (counts['countNextPast'] > 0)
        numerator = counts['countNextPastSource'][valid] / counts['countPastSource'][valid]
        denominator = counts['countNextPast'][valid] / counts['countPast'][valid]
        with np.errstate(divide='ignore', invalid='ignore'): # Ignore divide by zero and invalid (NaN) warnings in log
            log_terms = np.log(numerator / denominator)
            log_terms = np.where(np.isfinite(log_terms), log_terms, 0.0) # Replace inf/NaN with 0

        te_sum = np.sum(log_terms)
        self.last_average = te_sum / self.total_observations / self.LOG2
        return self.last_average

    def _compute_permutation_te(self, ordering):
        """
        Worker function to compute TE for a given permutation ordering.
        Each worker resets the estimator's source observations accordingly.
        """
        permuted_source = self.source_values[ordering]
        # Create a local copy of the estimator to avoid state conflicts in multiprocessing
        local_estimator = KernelEstimatorTransferEntropy(
            normalise=self.te_kernel_estimator.normalise,
            kernel_width=self.te_kernel_estimator.kernel_width,
            debug=self.te_kernel_estimator.debug
        )
        local_estimator.set_observations(self.dest_past_vectors, self.dest_next_values, permuted_source)
        return self.compute_average_local_of_observations_with_estimator(local_estimator, source_values=permuted_source)

    def compute_average_local_of_observations_with_estimator(self, estimator, source_values=None):
        """
        Compute average local TE using a provided estimator (used for permutation testing).

        Args:
            estimator: Kernel estimator configured with observations for the desired source.
            source_values (np.ndarray, optional): Source series aligned with destination windows.
                                                  Defaults to the calculator's stored source.
        """
        if source_values is None:
            source_values = self.source_values
        counts = estimator.get_counts_batch(self.dest_past_vectors, self.dest_next_values, source_values, self.indices)
        valid = (counts['countPastSource'] > 0) & (counts['countPast'] > 0) & (counts['countNextPast'] > 0)
        numerator = counts['countNextPastSource'][valid] / counts['countPastSource'][valid]
        denominator = counts['countNextPast'][valid] / counts['countPast'][valid]
        with np.errstate(divide='ignore', invalid='ignore'):
            log_terms = np.log(numerator / denominator)
            log_terms = np.where(np.isfinite(log_terms), log_terms, 0.0) # Replace inf/NaN with 0
        te_sum = np.sum(log_terms)
        return te_sum / self.total_observations / self.LOG2

    def compute_significance(self, num_permutations_to_check: int = None, new_orderings: list = None):
        """
        Compute permutation-based significance of Transfer Entropy.
        Parallelizes permutation testing using multiprocessing for speed.

        Args:
            num_permutations_to_check (int, optional): Number of permutations to generate and check.
                                                       If provided, random permutations are generated.
            new_orderings (list, optional): List of pre-generated permutation orderings.
                                            If provided, these orderings are used instead of generating new ones.

        Returns:
            EmpiricalMeasurementDistribution: Object containing the permutation distribution, p-value, and actual TE value.
        """
        if num_permutations_to_check is not None:
            if not isinstance(num_permutations_to_check, int) or num_permutations_to_check <= 0:
                raise ValueError("num_permutations_to_check must be a positive integer.")
            rng = np.random.default_rng()
            new_orderings = [rng.permutation(self.total_observations) for _ in range(num_permutations_to_check)]
        elif new_orderings is None:
            raise ValueError("Either num_permutations_to_check or new_orderings must be provided for significance testing.")
        elif not isinstance(new_orderings, list):
            raise TypeError("new_orderings must be a list of permutation arrays.")

        # Compute actual TE value using the current (original) source observations
        actual_te = self.compute_average_local_of_observations()
        old_source_values = self.source_values.copy() # Keep a copy to reset after permutation testing
        meas_distribution = EmpiricalMeasurementDistribution(distribution=np.zeros(len(new_orderings)), actual_value=actual_te)

        # Parallel permutation TE computation using process_map for progress bar and multiprocessing
        te_values = process_map(self._compute_permutation_te, new_orderings,
                                max_workers=multiprocessing.cpu_count(), desc="Permutations") # Consider setting max_workers manually if needed

        te_values = np.array(te_values)
        meas_distribution.distribution = te_values
        count_more_significant = np.sum(te_values >= actual_te)
        meas_distribution.p_value = count_more_significant / len(new_orderings)

        # Reset estimator to original source values after permutation testing
        self.te_kernel_estimator.set_observations(self.dest_past_vectors, self.dest_next_values, old_source_values)
        self.last_average = actual_te # Update last average to the actual TE
        return meas_distribution


def compute_te_for_target(
    dest_series: np.ndarray,
    sources: List[np.ndarray],
    k: int,
    kernel_width: float = 0.5,
    normalise: bool = True,
    reuse_dest_neighbors: bool = False,
) -> List[float]:
    """
    Compute Transfer Entropy values from many sources to a single destination
    by reusing destination-specific computations (past and next kd-trees).

    Args:
        dest_series: 1D array of destination time series values.
        sources: List of 1D arrays for each source time series.
        k: History length.
        kernel_width: Kernel width used in estimator.
        normalise: Whether to normalise by standard deviation.

    Returns:
        List of TE values, aligned with 'sources' order.
    """
    # Ensure float dtype to avoid integer division/casting errors
    dest_series = np.asarray(dest_series, dtype=float)
    N = len(dest_series)
    if N <= k:
        return [0.0 for _ in sources]

    # Prepare destination past and next arrays
    total_obs = N - k
    if _sliding_window_view is not None:
        # sliding_window_view returns (N-k+1, k); we need first (N-k) windows to align with next values
        dest_past_vectors = _sliding_window_view(dest_series, window_shape=k)[:-1].astype(float, copy=False)
    else:
        dest_past_vectors = np.empty((total_obs, k), dtype=float)
        for i in range(total_obs):
            dest_past_vectors[i] = dest_series[i:i + k]
    dest_next_values = dest_series[k:]

    # Kernel widths
    if normalise:
        epsilon = 1e-9
        std_dest_past = np.std(dest_past_vectors, axis=0, ddof=1)
        std_dest_next = np.std(dest_next_values, ddof=1)
        kw_dest_past = kernel_width * (std_dest_past + epsilon)
        kw_dest_next = kernel_width * (std_dest_next + epsilon)
    else:
        kw_dest_past = np.full(k, kernel_width, dtype=float)
        kw_dest_next = float(kernel_width)

    # Scale destination parts
    scaled_past = np.divide(dest_past_vectors, kw_dest_past,
                            out=np.zeros_like(dest_past_vectors, dtype=float),
                            where=kw_dest_past != 0)
    scaled_next = np.divide(dest_next_values.reshape(-1, 1), kw_dest_next,
                            out=np.zeros_like(dest_next_values.reshape(-1, 1), dtype=float),
                            where=kw_dest_next != 0)

    # Build kd-trees for destination-only contexts
    tree_past = KDTree(scaled_past, metric='chebyshev')
    X_past_next = np.hstack((scaled_past, scaled_next))
    tree_past_next = KDTree(X_past_next, metric='chebyshev')

    # Precompute counts and, optionally, neighbor indices independent of source
    count_past = tree_past.query_radius(scaled_past, r=1.0, count_only=True)
    count_next_past = tree_past_next.query_radius(X_past_next, r=1.0, count_only=True)

    neighbors_past = None
    neighbors_past_next = None
    if reuse_dest_neighbors:
        # Note: returning neighbor indices once lets us cheaply filter by source dimension later
        neighbors_past = tree_past.query_radius(scaled_past, r=1.0, return_distance=False)
        neighbors_past_next = tree_past_next.query_radius(X_past_next, r=1.0, return_distance=False)

    # Compute TE per source
    te_values: List[float] = []
    LOG2 = np.log(2)
    for src in sources:
        if len(src) != N:
            te_values.append(0.0)
            continue
        source_values = src[k - 1:-1]
        if normalise:
            epsilon = 1e-9
            std_source = np.std(source_values, ddof=1)
            kw_source = kernel_width * (std_source + epsilon)
        else:
            kw_source = float(kernel_width)

        scaled_source_1d = np.divide(
            source_values.astype(float), kw_source,
            out=np.zeros_like(source_values, dtype=float),
            where=(kw_source != 0)
        )

        if reuse_dest_neighbors and (neighbors_past is not None) and (neighbors_past_next is not None):
            # Compute counts by filtering destination neighbor lists with the 1D source constraint
            # Chebyshev radius in (past, source) space is equivalent to: max(past_diff_max, |source_i - source_j|) <= 1
            # Since neighbors_past already satisfy past_diff_max <= 1, remaining filter is |source_i - source_j| <= 1
            if njit is not None:
                indptr_past, indices_past = _neighbors_to_csr(neighbors_past)
                indptr_pn, indices_pn = _neighbors_to_csr(neighbors_past_next)
                count_past_source = _count_csr(indptr_past, indices_past, scaled_source_1d)
                count_next_past_source = _count_csr(indptr_pn, indices_pn, scaled_source_1d)
            else:
                def _count_with_source(nbrs_list, svals):
                    out_counts = np.empty(len(svals), dtype=np.int32)
                    for i, nbrs in enumerate(nbrs_list):
                        if nbrs.size == 0:
                            out_counts[i] = 0
                        else:
                            out_counts[i] = int(np.count_nonzero(np.abs(svals[nbrs] - svals[i]) <= 1.0))
                    return out_counts

                count_past_source = _count_with_source(neighbors_past, scaled_source_1d)
                count_next_past_source = _count_with_source(neighbors_past_next, scaled_source_1d)
        else:
            # Fallback to per-source KDTree (original slower path)
            scaled_source = scaled_source_1d.reshape(-1, 1)
            X_past_source = np.hstack((scaled_past, scaled_source))
            X_past_source_next = np.hstack((scaled_past, scaled_source, scaled_next))
            tree_past_source = KDTree(X_past_source, metric='chebyshev')
            tree_past_source_next = KDTree(X_past_source_next, metric='chebyshev')
            count_past_source = tree_past_source.query_radius(X_past_source, r=1.0, count_only=True)
            count_next_past_source = tree_past_source_next.query_radius(X_past_source_next, r=1.0, count_only=True)

        valid = (count_past_source > 0) & (count_past > 0) & (count_next_past > 0)
        with np.errstate(divide='ignore', invalid='ignore'):
            numerator = count_next_past_source[valid] / count_past_source[valid]
            denominator = count_next_past[valid] / count_past[valid]
            log_terms = np.log(numerator / denominator)
            log_terms = np.where(np.isfinite(log_terms), log_terms, 0.0)

        te_sum = float(np.sum(log_terms))
        te_values.append(te_sum / total_obs / LOG2)

    return te_values


def _radius_counts_ckdtree(tree: CKDTree, points: np.ndarray, r: float, workers: int = 1):
    nbrs = tree.query_ball_point(points, r=r, p=np.inf, workers=workers)
    # Convert list of arrays to counts (int32)
    return np.fromiter((len(v) for v in nbrs), dtype=np.int32, count=len(nbrs)), nbrs


def _radius_counts_sklearn(tree: KDTree, points: np.ndarray, r: float):
    return tree.query_radius(points, r=r, count_only=True), None


def prepare_dest_context(
    dest_series: np.ndarray,
    k: int,
    kernel_width: float = 0.5,
    normalise: bool = True,
    precompute_neighbors: bool = False,
    backend: str = "sklearn",
    workers: int = 1,
):
    """
    Precompute destination-only structures (scaled arrays, KDTree, counts) for reuse
    across many sources for the same target.

    Returns a dict context to be used with compute_te_for_sources_with_context.
    """
    dest_series = np.asarray(dest_series, dtype=float)
    N = len(dest_series)
    if N <= k:
        return {
            'N': N,
            'k': k,
            'kernel_width': kernel_width,
            'normalise': normalise,
            'total_obs': 0,
        }

    total_obs = N - k
    if _sliding_window_view is not None:
        dest_past_vectors = _sliding_window_view(dest_series, window_shape=k)[:-1].astype(float, copy=False)
    else:
        dest_past_vectors = np.empty((total_obs, k), dtype=float)
        for i in range(total_obs):
            dest_past_vectors[i] = dest_series[i:i + k]
    dest_next_values = dest_series[k:]

    if normalise:
        epsilon = 1e-9
        std_dest_past = np.std(dest_past_vectors, axis=0, ddof=1)
        std_dest_next = np.std(dest_next_values, ddof=1)
        kw_dest_past = kernel_width * (std_dest_past + epsilon)
        kw_dest_next = kernel_width * (std_dest_next + epsilon)
    else:
        kw_dest_past = np.full(k, kernel_width, dtype=float)
        kw_dest_next = float(kernel_width)

    scaled_past = np.divide(dest_past_vectors, kw_dest_past,
                            out=np.zeros_like(dest_past_vectors, dtype=float),
                            where=kw_dest_past != 0)
    scaled_next = np.divide(dest_next_values.reshape(-1, 1), kw_dest_next,
                            out=np.zeros_like(dest_next_values.reshape(-1, 1), dtype=float),
                            where=kw_dest_next != 0)

    # Build trees in selected backend
    if backend == "ckdtree":
        tree_past = CKDTree(scaled_past, leafsize=64)
        X_past_next = np.hstack((scaled_past, scaled_next))
        tree_past_next = CKDTree(X_past_next, leafsize=64)
        count_past, neighbors_past = _radius_counts_ckdtree(tree_past, scaled_past, r=1.0, workers=workers)
        count_next_past, neighbors_past_next = _radius_counts_ckdtree(tree_past_next, X_past_next, r=1.0, workers=workers)
        if not precompute_neighbors:
            neighbors_past = None
            neighbors_past_next = None
    else:
        tree_past = KDTree(scaled_past, metric='chebyshev')
        X_past_next = np.hstack((scaled_past, scaled_next))
        tree_past_next = KDTree(X_past_next, metric='chebyshev')
        count_past, _ = _radius_counts_sklearn(tree_past, scaled_past, r=1.0)
        count_next_past, _ = _radius_counts_sklearn(tree_past_next, X_past_next, r=1.0)
        neighbors_past = None
        neighbors_past_next = None

    return {
        'N': N,
        'k': k,
        'kernel_width': kernel_width,
        'normalise': normalise,
        'total_obs': total_obs,
        'scaled_past': scaled_past,
        'scaled_next': scaled_next,
        'tree_past': tree_past,
        'tree_past_next': tree_past_next,
        'count_past': count_past,
        'count_next_past': count_next_past,
        'neighbors_past': neighbors_past,
        'neighbors_past_next': neighbors_past_next,
        'backend': backend,
        'workers': workers,
    }


def compute_te_for_sources_with_context(
    ctx: dict,
    sources: List[np.ndarray],
    reuse_dest_neighbors: bool = False,
    dense_threshold: int = 250,
    return_local: bool = False,
) -> List[float] | tuple[List[float], List[np.ndarray]]:
    """
    Compute TE values for multiple sources using a precomputed destination context.
    """
    if ctx.get('total_obs', 0) <= 0:
        return [0.0 for _ in sources]

    k = ctx['k']
    kernel_width = ctx['kernel_width']
    normalise = ctx['normalise']
    total_obs = ctx['total_obs']
    scaled_past = ctx['scaled_past']
    scaled_next = ctx['scaled_next']
    count_past = ctx['count_past']
    count_next_past = ctx['count_next_past']
    backend = ctx.get('backend', 'ckdtree')
    workers = int(ctx.get('workers', 1))

    neighbors_past = ctx.get('neighbors_past') if reuse_dest_neighbors else None
    neighbors_past_next = ctx.get('neighbors_past_next') if reuse_dest_neighbors else None

    te_values: List[float] = []
    local_values: List[np.ndarray] | None = [] if return_local else None
    LOG2 = np.log(2)
    N = ctx['N']

    # Heuristic: if destination neighborhoods are dense, prefer per-source tree backend
    avg_nbr = max(float(np.mean(count_past)), float(np.mean(count_next_past)))
    prefer_per_source_tree = (avg_nbr >= dense_threshold) or (not reuse_dest_neighbors)

    for src in sources:
        if len(src) != N:
            te_values.append(0.0)
            continue
        source_values = np.asarray(src, dtype=float)[k - 1:-1]
        if normalise:
            epsilon = 1e-9
            std_source = np.std(source_values, ddof=1)
            kw_source = kernel_width * (std_source + epsilon)
        else:
            kw_source = float(kernel_width)

        scaled_source_1d = np.divide(
            source_values.astype(float), kw_source,
            out=np.zeros_like(source_values, dtype=float),
            where=(kw_source != 0)
        )

        if (not prefer_per_source_tree) and reuse_dest_neighbors and (neighbors_past is not None) and (neighbors_past_next is not None):
            if njit is not None:
                indptr_past, indices_past = _neighbors_to_csr(neighbors_past)
                indptr_pn, indices_pn = _neighbors_to_csr(neighbors_past_next)
                count_past_source = _count_csr(indptr_past, indices_past, scaled_source_1d)
                count_next_past_source = _count_csr(indptr_pn, indices_pn, scaled_source_1d)
            else:
                def _count_with_source(nbrs_list, svals):
                    out_counts = np.empty(len(svals), dtype=np.int32)
                    for i, nbrs in enumerate(nbrs_list):
                        if nbrs.size == 0:
                            out_counts[i] = 0
                        else:
                            out_counts[i] = int(np.count_nonzero(np.abs(svals[nbrs] - svals[i]) <= 1.0))
                    return out_counts

                count_past_source = _count_with_source(neighbors_past, scaled_source_1d)
                count_next_past_source = _count_with_source(neighbors_past_next, scaled_source_1d)
        else:
            scaled_source = scaled_source_1d.reshape(-1, 1)
            X_past_source = np.hstack((scaled_past, scaled_source))
            X_past_source_next = np.hstack((scaled_past, scaled_source, scaled_next))
            if backend == 'ckdtree':
                tree_ps = CKDTree(X_past_source, leafsize=64)
                tree_psn = CKDTree(X_past_source_next, leafsize=64)
                count_past_source, _ = _radius_counts_ckdtree(tree_ps, X_past_source, r=1.0, workers=workers)
                count_next_past_source, _ = _radius_counts_ckdtree(tree_psn, X_past_source_next, r=1.0, workers=workers)
            else:
                tree_ps = KDTree(X_past_source, metric='chebyshev')
                tree_psn = KDTree(X_past_source_next, metric='chebyshev')
                count_past_source, _ = _radius_counts_sklearn(tree_ps, X_past_source, r=1.0)
                count_next_past_source, _ = _radius_counts_sklearn(tree_psn, X_past_source_next, r=1.0)

        valid = (count_past_source > 0) & (count_past > 0) & (count_next_past > 0)
        with np.errstate(divide='ignore', invalid='ignore'):
            numerator = count_next_past_source[valid] / count_past_source[valid]
            denominator = count_next_past[valid] / count_past[valid]
            log_terms = np.log(numerator / denominator)
            log_terms = np.where(np.isfinite(log_terms), log_terms, 0.0)

        te_sum = float(np.sum(log_terms))
        te = te_sum / total_obs / LOG2
        te_values.append(te)
        if return_local and local_values is not None:
            local_bits = np.zeros(total_obs, dtype=np.float32)
            if np.any(valid):
                local_bits[valid] = log_terms / LOG2
            local_values.append(local_bits)

    if return_local and local_values is not None:
        return te_values, local_values
    return te_values
