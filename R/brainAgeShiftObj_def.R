################################################################################
# brainAgeShiftObj class constructor.                                          #
################################################################################


# Creates a brainAgeShift object.
# counts needs to be a genome-wide counts matrix with samples in columns and
# genes in rows.
# metadata needs to be a data.frame with at least the variable of interest the
# user wants to make comparisons of transcriptomic age, with the row names
# matching the column names of the counts matrix. If metadata is not given the
# counts still can be normalized and ages can be predicted, but downstream
# analyses won't be performed.
# variable. The variable where the comparisons factors are stored.
# comparisons. either a vector of length = 2, or a list of vectors of length 2.
# The first item in each comparison should be the initial or base state in the
# comparison, and the second one the final or altered state, so in the output
# the log2FC will be positive if in the samples belonging to comparisons[[n]][2]
# are higher, and the other way around.

#' Create a brainAgeShift object.
#'
#' Creates a brainAgeShift object to do brain transcriptomic age predictions,
#' and subsequent statistical analyses to detect transcriptomic age-shifting
#' factors.
#'
#' @param counts A matrix of gene counts. Genes in rows, samples in columns. Row names must be ENSEMBL IDs
#' @param metadata A data.frame with at least the variable of interest for comparing transcriptomic ages. Row names must match column names of counts matrix.
#' @param variable The column in the metadata dataframe where the comparisons are stored.
#' @param comparisons Either a vector of two factors, or a list of multiple two-factor vectors. First item in each comparison will be initial or base state during the difference computation, second the final or altered state.
#'
#' @return A new brainAgeShift object initialized with the required slots and default values.
#' @export
create_brainAgeShiftObj <- function(counts,
                                    metadata = NULL,
                                    variable = NULL,
                                    comparisons = NULL){
        if(!is.null(metadata)){
                if (is.null(variable) | is.null(comparisons)){
                        stop("If metadata is included both variable and comparisons need to be indicated.",
                             call. = F)
                }
                if (ncol(counts) != nrow(metadata)){
                        stop("Number of columns of counts and number of rows of metadata differ.",
                             call. = F)
                }
                if (!all(colnames(counts) == rownames(metadata))){
                        stop("Column names of counts and row mames of metadata differ.",
                             call. = F)
                }
                if (!variable %in% colnames(metadata)){
                        stop("The introduced variable is not included in the metadata.",
                             call. = F)
                }
                if (!is.list(comparisons)){
                        comparisons <- list(comparisons)
                }
                if (any(unlist(lapply(comparisons,
                                      function(x) length(x) != 2)))){
                        stop("Comparisons must be of length 2.",
                             call. = F)
                }
                comps_uniq <- unique(unlist(comparisons))
                if (any(!comps_uniq %in% metadata[, variable])){
                        notInMetDat <- comps_uniq[!comps_uniq %in% metadata[, variable]]
                        paste(notInMetDat, collapse = ", ")
                        stop(sprintf("%s not included in the %s.",
                                     notInMetDat,
                                     variable),
                             call. = F)
                }
                if(any(!metadata[, variable] %in% comps_uniq)){
                        print(sprintf("Filtering out the samples that are not labeled as %s.",
                                      paste(comps_uniq, collapse = ", ")))
                        keep <- metadata[, variable] %in% comps_uniq
                        metadata <- metadata[keep, ]
                        counts <- counts[, keep]
                }
        }else{
                warning("No metadata, variable and comparison indicated. This object only will work for doing predictions with the clock.")
                if (!is.null(variable)){
                        warning("If metadata is not provided variable won't be used.")
                        variable <- NULL
                }
                if (!is.null(comparisons)){
                        warning("If metadata is not provided comparisons won't be used.")
                        comparisons <- NULL
                }
        }
        obj <- list(counts = counts,
                    norm_counts = NULL,
                    metadata = metadata,
                    variable = variable,
                    comparisons = comparisons,
                    stats = NULL,
                    sign_genes = NULL)
        obj$default_slot <- "counts"
        obj$used_for_preds <- NULL
        class(obj) <- "brainAgeShiftObj"
        return(obj)
}
