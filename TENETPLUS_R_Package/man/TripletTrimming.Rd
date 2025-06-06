% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TripletTrimming.R
\name{TripletTrimming}
\alias{TripletTrimming}
\title{Triplet Trimming for Regulatory Network Analysis}
\usage{
TripletTrimming(object, trim_indirect = TRUE, trim_distance = 1e+06)
}
\arguments{
\item{object}{An object (often an S4 or list) containing the necessary input and output data.
It is expected to have:
\itemize{
\item \code{object@input$feature_list$gtf}: A data.frame with gene annotations, including columns
\code{gene}, \code{gene_chr}, \code{gene_start}, and \code{gene_end}.
\item \code{object@output$result$TAR_TG}: A data.frame with peak data. The \code{source} column should
contain strings formatted as "chr-start-end", along with columns such as \code{score} and \code{target}.
\item \code{object@output$result$TF_TG_indirect}: A data.frame of indirect TF-target interactions.
\item \code{object@output$result$TF_TG}: A data.frame of direct TF-target interactions.
\item \code{object@output$result$TF_TAR}: A data.frame of TF-peak interactions, with a \code{TAR} column.
}}

\item{trim_indirect}{Logical; if \code{TRUE}, indirect TF-target interactions (\code{TF_TG_indirect}) will
be used for trimming. If \code{FALSE}, direct interactions (\code{TF_TG}) are used. Default is \code{TRUE}.}

\item{trim_distance}{Numeric or character; specifies the maximum allowed distance between peaks and gene
annotations for a valid interaction. When a numeric value is provided, only peak-gene pairs with a
calculated distance less than \code{trim_distance} are retained. If set to \code{"None"}, no distance-based
trimming is applied. Default is \code{1000000}.}
}
\value{
The input \code{object} with updated results in its \code{@output$result} slot:
\itemize{
\item \code{$TAR_TG}: A data.frame with added columns for peak chromosome (\code{peak_chr}),
start (\code{peak_start}), end (\code{peak_end}), and the calculated \code{distance}
between the peak and gene annotation.
\item \code{$Trimm}: A data.frame representing the trimmed triplet regulatory network. It includes
columns \code{TF}, \code{TG}, \code{TAR}, \code{TF_TG_TE}, \code{TF_TAR_TE}, \code{TAR_TG_TE},
and \code{distance}.
}
}
\description{
This function performs triplet trimming on regulatory network data stored in a given object.
It calculates the genomic distance between peak regions and gene annotations, and then trims
the triplet network based on a specified distance threshold and the choice between direct or
indirect transcription factor (TF) to target gene (TG) interactions.
}
\details{
The function executes the following steps:
\enumerate{
\item Extracts and parses the peak data from the \code{source} column to retrieve chromosome,
start, and end positions.
\item Processes the gene annotation data by ensuring numeric types for start and end positions,
and removes duplicate gene entries.
\item Merges the peak data with gene annotation data and verifies that the chromosome information
matches. An error is raised if mismatches are found.
\item Calculates the distance between each peak and the corresponding gene using a custom
distance function.
\item Depending on the \code{trim_indirect} flag, selects either indirect or direct TF-target
interaction data for further processing.
\item Merges the filtered peak-gene list with TF-peak and TF-target interaction data to form a
trimmed triplet network based on the specified \code{trim_distance}.
}
}
\examples{
\dontrun{
# Assuming 'my_object' is a properly formatted object:
result_object <- TripletTrimming(my_object, trim_indirect = TRUE, trim_distance = 500000)
}

}
