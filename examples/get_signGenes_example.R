# Assuming obj is a brainAgeShift object with non-empty stats slot:
obj <- get_signGenes(obj,
                     alpha_comparisons = 0.05, 
                     alpha_genes = 0.01,
                     n_perms = 1000, 
                     adjust_method = "BH",
                     sort_genes = TRUE)
