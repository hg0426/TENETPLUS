import numpy as np

class TransferEntropyCalculatorKernel:
    def __init__(self, k=1, kernel_width=0.25, normalize=True, dyn_corr_excl_time=0, force_compare_to_all=False):
        self.k = k
        self.kernel_width = kernel_width
        self.normalize = normalize
        self.dyn_corr_excl_time = dyn_corr_excl_time
        self.force_compare_to_all = force_compare_to_all
        self.te_kernel_estimator = KernelEstimatorTransferEntropy(k, kernel_width, normalize, dyn_corr_excl_time, force_compare_to_all)
        self.dest_past_vectors = None
        self.dest_next_values = None
        self.source_values = None
        self.last_average = None
        self.total_observations = 0
        self.added_more_than_one_observation_set = False

    def set_property(self, property_name, property_value):
        if property_name == "KERNEL_WIDTH" or property_name == "EPSILON":
            self.kernel_width = float(property_value)
            self.te_kernel_estimator.set_kernel_width(self.kernel_width)
        elif property_name == "NORMALISE":
            self.normalize = bool(property_value)
            self.te_kernel_estimator.set_normalize(self.normalize)
        elif property_name == "DYN_CORR_EXCL_TIME":
            self.dyn_corr_excl_time = int(property_value)
            self.dyn_corr_excl = self.dyn_corr_excl_time > 0
            self.te_kernel_estimator.set_dynamic_correlation_exclusion(self.dyn_corr_excl_time)
        elif property_name == "FORCE_KERNEL_COMPARE_TO_ALL":
            self.force_compare_to_all = bool(property_value)
            self.te_kernel_estimator.set_force_compare_to_all(self.force_compare_to_all)
        else:
            # Try superclass:
            super().set_property(property_name, property_value)

    def get_property(self, property_name):
        if property_name == "KERNEL_WIDTH" or property_name == "EPSILON":
            return self.kernel_width
        elif property_name == "NORMALISE":
            return self.normalize
        elif property_name == "DYN_CORR_EXCL_TIME":
            return self.dyn_corr_excl_time
        elif property_name == "FORCE_KERNEL_COMPARE_TO_ALL":
            return self.force_compare_to_all
        else:
            # Try superclass:
            return super().get_property(property_name)

    def initialise(self, k):
        self.k = k
        self.te_kernel_estimator.initialise(k)
        self.dest_past_vectors = None
        self.dest_next_values = None
        self.source_values = None
        self.last_average = None
        self.total_observations = 0
        self.added_more_than_one_observation_set = False

    def finalise_add_observations(self, destination_observations, source_observations):
        self.total_observations = 0
        for dest_obs in destination_observations:
            self.total_observations += len(dest_obs) - self.k
        self.dest_past_vectors = np.zeros((self.total_observations, self.k))
        self.dest_next_values = np.zeros(self.total_observations)
        self.source_values = np.zeros(self.total_observations)

        start_observation = 0
        for dest_obs, source_obs in zip(destination_observations, source_observations):
            past_vectors = self.make_joint_vector_for_past(dest_obs)
            self.dest_past_vectors[start_observation: start_observation+past_vectors.shape[0], :] = past_vectors
            self.dest_next_values[start_observation: start_observation+len(dest_obs)-self.