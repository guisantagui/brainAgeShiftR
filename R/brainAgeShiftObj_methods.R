################################################################################
# brainAgeShiftObj methods.                                                    #
################################################################################

#' Change the default_slot
#'
#' Changes the default_slot, which is the one that will be used by predictAge.
#' to compute transcriptomic age.
#'
#' @param obj A `brainAgeShift` object.
#'
#' @example examples/default_slot_example.R
#'
#' @export
default_slot <- function(obj, ...){
        obj$default_slot
}

#' @rdname default_slot
#' @export
`default_slot<-` <- function(obj, value) {
        UseMethod("default_slot<-")
}

#' @rdname default_slot
#' @method default_slot<- brainAgeShiftObj
#' @export
`default_slot<-.brainAgeShiftObj` <- function(obj, value) {
        if (!value %in% c("counts", "norm_counts", "frozen_SVAed")) {
                stop("default_slot must be either counts, norm_counts or frozen_SVAed")
        }
        obj$default_slot <- value
        obj  # Return the modified object
}

# Performs quantile normalization
quantNorm <- function(m, axis = 2, train_means = NULL) {
        if (axis == 1) {
                m <- t(m)
        }
        # Get ranks of the matrix
        m_rank <- apply(m, 2, rank, ties.method = "average")

        # Sort the matrix
        m_sort <- apply(m, 2, sort)

        # Calculate row means of the sorted matrix
        if(!is.null(train_means)){
                means <- train_means
        }else{
                means <- rowMeans(m_sort)
        }

        # Create a normalized matrix with the same dimensions
        m_norm <- matrix(0, nrow = nrow(m), ncol = ncol(m))

        for (i in 1:ncol(m)) {
                m_norm[, i] <- means[rank(m[, i], ties.method = "average")]
        }

        if (axis == 1) {
                m_norm <- t(m_norm)
        }

        dimnames(m_norm) <- dimnames(m)
        return(m_norm)
}

# Normalize the counts of the brainAgeShift object to make them suitable
# to predict ages with it.

#' Normalize brainAgeShift counts
#'
#' Normalizes the counts in a `brainAgeShift` object to prepare them for age
#' prediction. The normalization process includes filtering genes to retain only
#' those used in training, applying log2 transformation (`log2(counts + 1)`),
#' and performing quantile normalization. After normalization, the
#' `default_slot` of the object is updated to `"norm_counts"`.
#'
#' @param obj A `brainAgeShift` object.
#' @param useTrainMeans Logical. If `TRUE`, uses the quantile means from the
#' training data for nomralization, ensuring distributions are closer to the
#' training set. If `FALSE`, computes quantile means directly from the input
#' data.
#'
#' @return A `brainAgeShift` object with the `norm_counts` slot filled with the
#' normalized counts matrix.
#'
#' @example examples/normalizeCounts_example.R
#'
#' @export
normalizeCounts <- function(obj, ...){
        UseMethod("normalizeCounts")
}

#' @rdname normalizeCounts
#' @method normalizeCounts brainAgeShiftObj
#' @export
normalizeCounts.brainAgeShiftObj <- function(obj, useTrainMeans = F){
        bg_genes <- .brainAgeShiftR_env$bg_genes
        obj$norm_counts <- obj$counts
        trainGenes_notInDat <- bg_genes[!bg_genes %in% rownames(obj$counts)]
        obj$norm_counts <- obj$norm_counts[rownames(obj$norm_counts) %in% bg_genes, ]
        obj$norm_counts <- log2(obj$norm_counts + 1)
        if (useTrainMeans){
                print("Quantile normalizing with training set gene means.")
                if (length(trainGenes_notInDat) != 0){
                        warning(sprintf("%s genes of our training data (%s %%) are not in the counts matrix. The top %s quantiles of the training data won't be considered.",
                                        length(trainGenes_notInDat),
                                        round(length(trainGenes_notInDat)/length(bg_genes) * 100, digits = 3),
                                        length(trainGenes_notInDat)))
                        print("Training genes not in counts data:")
                        print(trainGenes_notInDat)
                }
                train_means <- .brainAgeShiftR_env$train_quant_means[, 1]
                obj$norm_counts <- quantNorm(obj$norm_counts,
                                             train_means = train_means)
        }else{
                print("Quantile normalizing with computed gene means.")
                obj$norm_counts <- quantNorm(obj$norm_counts)
        }
        obj$default_slot <- "norm_counts"
        return(obj)
}

# Apply frozen SVA to the norm_counts matrix. SVs are predicted based on a
# linear model we obtained with our training data, and are regressed out from
# the norm_counts matrix.

#' Apply frozen SVA
#'
#' Applies frozen SVA on the `norm_counts` slot of a `brainAgeShift` object.
#' Predicts SVs based on the ones computed in our training data, and regresses
#' them out. Default slot is changed to `frozen_SVAed` after running this
#' method.
#'
#' @param obj A `brainAgeShift` object.
#'
#' @return A `brainAgeShift` object with the `frozen_SVAed` slot filled with the
#' SVAed normalized counts matrix.
#'
#' @example examples/normalizeCounts_example.R
#'
#' @export
do_frozenSVA <- function(obj, ...){
        UseMethod("do_frozenSVA")
}

#' @rdname do_frozenSVA
#' @method do_frozenSVA brainAgeShiftObj
#' @export
do_frozenSVA.brainAgeShiftObj <- function(obj){
        if(is.null(obj$norm_counts)){
                stop(sprintf("Frozen SVA needs to be applied on norm_counts, and this slot is not available. Run normalizeCounts() on your object before."),
                     call. = F)
        }
        sva_mod_coefs <- .brainAgeShiftR_env$sva_mod_coefs
        sva_mod_genes <- rownames(sva_mod_coefs)
        sva_mod_genes <- sva_mod_genes[sva_mod_genes != "(Intercept)"]
        sva_mod_genes_notInDat <- sva_mod_genes[!sva_mod_genes %in% rownames(obj$norm_counts)]
        if(length(sva_mod_genes_notInDat) > 0){
                stop(sprintf("%s genes are not available in the input data. Frozen SVA cannot be applied.",
                             paste(sva_mod_genes_notInDat, collapse = ", ")),
                     call. = F)
        }
        predMat <- data.frame(t(obj$norm_counts[match(sva_mod_genes,
                                                      rownames(obj$norm_counts)), ]))
        pred <- function(expVec, coefVec){
                predSV <- sum(expVec * coefVec[2:length(coefVec)]) + coefVec[1]
                return(predSV)
        }
        pred_SVs <- apply(sva_mod_coefs,
                          2,
                          function(y) apply(predMat,
                                            1,
                                            function(x) pred(x, y)))
        b <- matrix(nrow = 0, ncol = ncol(pred_SVs),
                    dimnames = list(NULL, colnames(pred_SVs)))
        for(i in seq_along(colnames(predMat))){
                form <- paste(colnames(predMat)[i],
                              paste(colnames(pred_SVs), collapse = " + "),
                              sep = " ~ ")
                lm_fit <- lm(as.formula(form), cbind(predMat, pred_SVs))
                toBind <- coef(lm_fit)[2:length(coef(lm_fit))]
                toBind <- matrix(toBind, nrow = 1,
                                 dimnames = list(colnames(predMat)[i],
                                                 names(toBind)))
                b <- rbind(b, toBind)
        }
        b[is.na(b)] <- 0
        res <- t(predMat) - b %*% t(pred_SVs)
        obj$frozen_SVAed <- res
        default_slot(obj) <- "frozen_SVAed"
        return (obj)
}

# Uses norm_counts of brainAgeShiftObj to predict age with the brain clock.
# Stores result as a new column in the metadata dataframe within the
# brainAgeShiftObj.

#' Predict ages
#'
#' Predicts transcriptomic ages using the brain transcriptomic clock. The method
#' operates based on the default slot, which must be one of `counts`,
#' `norm_counts` or `frozen_SVAed`.
#'
#' @param obj A `brainAgeShift` object.
#'
#' @return A `brainAgeShift` object with a `predicted_age` column added to the
#' metadata slot's dataframe.
#'
#' @example examples/predictAge_example.R
#'
#' @export
predictAge <- function(obj, ...){
        UseMethod("predictAge")
}

#' @rdname predictAge
#' @method predictAge brainAgeShiftObj
#' @export
predictAge.brainAgeShiftObj <- function(obj){
        print(sprintf("Predicting ages on %s slot...",
                      obj$default_slot))
        normed4Pred <- obj[[obj$default_slot]]
        mod_coef <- .brainAgeShiftR_env$mod_coef
        mod_coef <- mod_coef[mod_coef$coefficients != 0, ]
        mod_genes <- mod_coef$names
        mod_genes <- mod_genes[mod_genes != "Intercept"]
        normed4Pred <- normed4Pred[match(mod_genes, rownames(normed4Pred)), ]
        if (!all(mod_genes %in% rownames(normed4Pred))){
                missGenes <- mod_genes[!mod_genes %in% rownames(normed4Pred)]
                warning(sprintf("%s genes in the model don't appear in the input data. Their values will be imputed to zero. The missing genes are:",
                                length(missGenes)))
                print(missGenes)
                normed4Pred[is.na(rownames(normed4Pred)), ] <- 0
                rownames(normed4Pred)[is.na(rownames(normed4Pred))] <- missGenes
        }
        # Helper function for predicting the ages using the coefficient vector.
        # Only usable within this scope, as predictAge orders expVec according
        # to
        pred <- function(expVec, coefVec){
                predAge <- sum(expVec * coefVec[2:length(coefVec)]) + coefVec[1]
                return(predAge)
        }
        predVec <- apply(normed4Pred, 2, pred, coefVec = mod_coef$coefficients)
        if (is.null(obj$metadata)){
                predDF <- data.frame(predicted_age = predVec)
                rownames(predDF) <- colnames(obj$counts)
                obj$metadata <- predDF
        }
        obj$metadata$predicted_age <- predVec
        obj$used_for_preds <- obj$default_slot
        print("Done!")
        return(obj)
}

# Compute statistical significance for predicted age comparisons


#' Compute statistical significance for predicted age comparisons
#'
#' Performs statistical analysis of the transcriptomic ages on the specified
#' comparisons within the metadata slot of a `brainAgeShift` object. The
#' function computes the p-value (t-test), adjusted p-value, log2 fold change
#' and difference for each requested comparison.
#'
#' @param obj A `brainAgeShift` object with a non-empty `metadata`, `variable`
#' and `comparisons` slots.
#'
#' @param adjust_method A character string specifying the the p-value adjustment
#' method. This is passed to `p.adjust`, such as `"BH"` (Benjamini-Hochberg) or
#' `"bonferroni"`. Default is `"BH"`.
#'
#' @return A `brainAgeShift` object with the `stats` slot updated with a
#' dataframe containing the following columns:
#' \itemize{
#'      \item \strong{comparison}: The name of the comparison.
#'      \item \strong{p-value}: The p-value from the t-test
#'      \item \strong{log2FC}: log2 fold change.
#'      \item \strong{difference}: The raw difference between group means.
#'      \item \strong{p_adj}: The adjusted p-value.
#' }
#'
#' @example examples/do_signTest_example.R
#'
#' @export
do_signTest <- function(obj, ...){
        UseMethod("do_signTest")
}

#' @rdname do_signTest
#' @method do_signTest brainAgeShiftObj
#' @export
do_signTest.brainAgeShiftObj <- function(obj, adjust_method = "BH"){
        if (is.null(obj$metadata) | !any(colnames(obj$metadata) != "predicted_age")){
                stop("The object introduced doesn't have a proper metadata slot.",
                     call. = F)
        }
        if (is.null(obj$variable) | is.null(obj$comparisons)){
                stop("The object introduced doesn't have any indicated variable and/or comparisons.",
                     call. = F)
        }
        stats_df <- data.frame(matrix(nrow = 0,
                                      ncol = 4,
                                      dimnames = list(NULL,
                                                      c("comparison",
                                                        "p_value",
                                                        "log2FC",
                                                        "difference"))))
        for (i in seq_along(obj$comparisons)){
                c_1 <- obj$comparisons[[i]][1]
                c_2 <- obj$comparisons[[i]][2]
                comp_name <- sprintf("%s_vs_%s",
                                     c_2, c_1)
                c_1_vec <- obj$metadata$predicted_age[obj$metadata[, obj$variable] == c_1]
                c_2_vec <- obj$metadata$predicted_age[obj$metadata[, obj$variable] == c_2]
                p_val <- t.test(c_1_vec, c_2_vec)$p.value
                log2FC <- log2(mean(c_2_vec)/mean(c_1_vec))
                diff <- mean(c_2_vec) - mean(c_1_vec)
                toBind <- data.frame(comparison = comp_name,
                                     p_value = p_val,
                                     log2FC = log2FC,
                                     difference = diff)
                stats_df <- rbind.data.frame(stats_df, toBind)
        }
        stats_df$p_adj <- p.adjust(stats_df$p_value, method = adjust_method)
        obj$stats <- stats_df
        return(obj)
}

# Given a brainAgeShift object which has metadata slot with computed
# predicted ages, and a comparison (in the format of C2_vs_C1), returns a
# dataframe with geneIDs of the genes the clock uses, mean differences
# (C2 - C1), coefficients and weighted differences (difference * coefficient).
# This function is a helper function for do_permTest.
getWeightDiffDF <- function(obj, comp){
        mod_coef <- .brainAgeShiftR_env$mod_coef
        mod_coef <- mod_coef[mod_coef$coefficients != 0, ]
        c_2 <- gsub("\\_vs_.*", "", comp)
        c_1 <- gsub(".*_vs_", "", comp)
        c_2_samps <- rownames(obj$metadata)[obj$metadata[, obj$variable] == c_2]
        c_1_samps <- rownames(obj$metadata)[obj$metadata[, obj$variable] == c_1]
        c_2_samps_idxs <- which(colnames(obj[[obj$used_for_preds]]) %in% c_2_samps)
        c_1_samps_idxs <- which(colnames(obj[[obj$used_for_preds]]) %in% c_1_samps)
        diffs <- apply(obj[[obj$used_for_preds]],
                       1,
                       function(x) mean(x[c_2_samps_idxs]) - mean(x[c_1_samps_idxs]))
        diffs <- diffs[names(diffs) %in% mod_coef$names]
        diffs_df <- data.frame(geneID = names(diffs),
                               difference = diffs,
                               coefficient = mod_coef$coefficients[match(names(diffs),
                                                                         mod_coef$names)])
        diffs_df$weighted_diff <- diffs_df$difference * diffs_df$coefficient
        return(diffs_df)
}

# Given a brainAgeShift object which has metadata slot with computed
# predicted ages, a comparison (in the format of C2_vs_C1), the number of
# permutations and the p_adjust method, returns a dataframe with geneIDs of the
# genes the clock uses, mean differences (C2 - C1), coefficients, weighted
# differences (difference * coefficient), and p_value and p_adj columns, which
# are computed with a permutation test, obtaining a null distribution by
# randomly permuting the labels of the comparison. This function is a helper
# function for get_signGenes.
do_permTest <- function(obj, comp, n_perms, adjust_method = "BH"){
        real_weight_df <- getWeightDiffDF(obj, comp)
        perms_df <- data.frame(matrix(ncol = 0, nrow = nrow(real_weight_df),
                                      dimnames = list(rownames(real_weight_df),
                                                      NULL)))
        pb <- txtProgressBar(min = 0, max = n_perms, style = 3)
        for(i in 1:n_perms){
                obj_perm <- obj
                obj_perm$metadata[, obj_perm$variable] <- sample(obj$metadata[, obj$variable],
                                                                 size = nrow(obj$metadata),
                                                                 replace = F)
                perm_weight_df <- getWeightDiffDF(obj_perm, comp)
                toBind <- data.frame(matrix(perm_weight_df$weighted_diff,
                                            ncol = 1,
                                            dimnames = list(rownames(perm_weight_df),
                                                            sprintf("perm_%s", i))))
                perms_df <- cbind.data.frame(perms_df, toBind)
                setTxtProgressBar(pb, i)
        }
        close(pb)
        pVals <- rowSums(perms_df >= abs(real_weight_df$weighted_diff)) / n_perms
        real_weight_df$p_value <- pVals
        real_weight_df$p_adj <- p.adjust(pVals, method = adjust_method)
        return(real_weight_df)
}

# Given a brainAgeShift object which has metadata slot with computed predicted
# ages and a stats slot, computes the level of significance of the contribution
# each gene in the brain clock had towards the observed shift in predicted age
# for each comparison that was significantly shifted in age.
# alpha_comparisons is the level of significance for the comparison.
# alpha_genes is the level of significance for the genes.
# n_perms is the number of permutations for obtaining the null distribution
# adjust_method is the method for adjusting the p_value.
# sort_genes is a logical that indicates if significant genes should be sorted
# according to the magitude of their countribution (weighted difference). If the
# mean predicted age difference of the comparison is positive they would be
# sorted in a decreasing manner, and if it's negative the other way around.

#' Identify genes significantly contributing to age shifts
#'
#' For each significant comparison in the `stats` slot of a `brainAgeShift`
#' object, this function identifies the genes that significantly contribute to
#' observed shifts. The weighted mean difference between comparison groups is
#' calculated for each gene and compared against a null distribution generated
#' through random permutations of group labels. P-values are obtained for each
#' gene and adjusted for multiple comparisons.
#'
#' @param obj A `brainAgeShift` object with a non-empty `stats` slot.
#' @param alpha_comparisons Numeric. The significance threshold for selecting
#' significant comparisons (adjusted p-value <= `alpha_comparisons`). Defaults
#' to `0.05`.
#' @param alpha_genes Numeric. The significance threshold for selecting
#' significant genes (adjusted p-value <= `alpha_genes`). Defaults to `0.05`.
#' @param n_perms Integer. The number of permutations to perform for the
#' permutation test. Defaults to `1000`.
#' @param adjust_method Character. The p-value adjustment method to use for
#' genes' p-values (e.g. `"BH"` for Benjamini-Hochberg or `"bonferroni"`).
#' Defaults to `"BH"`.
#' @param sort_genes Logical. If `TRUE`, genes in each `sign_genes` dataframe
#' will be sorted by weighted differences. Sorting is performed such that
#' positive log2 fold change in the comparison implies decreasing order of
#' genes' weighted differences, and viceversa. Defaults to `TRUE`.
#'
#' @return A `brainAgeShift` object with the `sign_genes` slot updated. The
#' `sign_genes` slot contains a list where each element corresponds to a
#' significant comparison and is a dataframe of significant genes. The
#' dataframes include the following columns:
#' \itemize{
#'      \item \strong{geneID}: the ENSEMBL ID of the gene.
#'      \item \strong{difference}: the raw mean difference between groups.
#'      \item \strong{coefficient}: the brain clock coefficient.
#'      \item \strong{weighted_diff}: the weighted mean difference
#'      \item \strong{p_value}: the p-value from the permutation test.
#'      \item \strong{p_adj}: the adjusted p-value.
#' }
#'
#' @example examples/get_signGenes_example.R
#'
#' @export
get_signGenes <- function(obj, ...){
        UseMethod("get_signGenes")
}

#' @rdname get_signGenes
#' @method get_signGenes brainAgeShiftObj
#' @export
get_signGenes.brainAgeShiftObj <- function(obj,
                                           alpha_comparisons = 0.05,
                                           alpha_genes = 0.05,
                                           n_perms = 1000,
                                           adjust_method = "BH",
                                           sort_genes = T){
        mod_coef <- .brainAgeShiftR_env$mod_coef
        mod_coef <- mod_coef[mod_coef$coefficients != 0, ]
        mod_gene_info <- .brainAgeShiftR_env$mod_gene_info
        if(is.null(obj$stats)){
                stop("The object introduced doen't have a stats slot.",
                     call. = F)
        }
        sign_comps <- obj$stats[obj$stats$p_adj <= alpha_comparisons, ]
        signGenesList <- list()
        for (i in seq_along(sign_comps$comparison)){
                comp <- sign_comps$comparison[i]
                print(sprintf("Obtaining significant genes for %s comparison...",
                              gsub("_", " ", comp)))
                signGenesDF <- do_permTest(obj, comp, n_perms, adjust_method)
                signGenesDF <- signGenesDF[signGenesDF$p_adj <= alpha_genes, ]
                if(sort_genes){
                        if(obj$stats$difference[obj$stats$comparison == comp] > 0){
                                signGenesDF <- signGenesDF[order(signGenesDF$weighted_diff,
                                                                 decreasing = T), ]
                        }else{
                                signGenesDF <- signGenesDF[order(signGenesDF$weighted_diff,
                                                                 decreasing = F), ]
                        }
                }
                signGenes_info <- mod_gene_info[match(signGenesDF$geneID,
                                                      mod_gene_info$ensembl_gene_id),
                                                colnames(mod_gene_info) != "ensembl_gene_id"]
                signGenesDF <- cbind.data.frame(signGenesDF, signGenes_info)
                signGenesDF <- signGenesDF[, c("geneID",
                                               "hgnc_symbol",
                                               "entrezgene_id",
                                               "description",
                                               "difference",
                                               "coefficient",
                                               "weighted_diff",
                                               "p_value",
                                               "p_adj")]
                rownames(signGenesDF) <- 1:nrow(signGenesDF)
                colnames(signGenesDF) <- gsub("geneID", "ensembl_gene_id")
                signGenesList[[comp]] <- signGenesDF
        }
        obj$sign_genes <- signGenesList
        return (obj)
}
