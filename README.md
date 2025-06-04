# brainAgeShiftR

## Overview
An R package for detecting shifts in transcriptomic age in samples originating from human brain cells or tissues. Starting from a **matrix of gene counts**, and a set of experimental factors and comparisons:
- Performs preprocessing steps.
- Computes **transcriptomic age**.
- Detects **significant transcriptomic age shifts**.
- Identifies the **genes driving these shifts**.

## Installation

### From GitHub

```r
if (!require("devtools", quietly = T))
{
	install.packages("devtools")
}

devtools::install_github("guisantagui/brainAgeShiftR")
```

## Usage
Here's a basic example of how to use **brainAgeShiftR**

### Load the package
```r
library(brainAgeShiftR)
```
### Load example data
The example data consists on transcriptomic data originating from the DLPFC of healthy individuals. The included samples are a selection from the BrainSeq Phase 1 study[1]. They are grouped in young (age <= 40) and old (age >= 70).
```r
data(example_data)
```
### Create a `brainAgeShiftObj`
```r
baso_old_vs_young <- create_brainAgeShiftObj(counts = example_counts,
                                             metadata = example_metadata,
                                             variable = "group",
                                             comparisons = c("young", "old"))
```
### Normalize the counts data
```r
baso_old_vs_young <- normalizeCounts(baso_old_vs_young, useTrainMeans = T)
```
### Compute transcriptomic age
```r
baso_old_vs_young <- predictAge(baso_old_vs_young)
```
### Compute significance of the transcriptomic age shift of the comparison.
```r
# Compute significance.
baso_old_vs_young <- do_signTest(baso_old_vs_young, adjust_method = "BH")

# Access results of the significance test(s).
baso_old_vs_young$stats
```
### Get the genes that are significantly contributing to the shift
In some cases, the observed shifts in transcriptional age are driven by a small number of dominant genes. The function `get_signGenes` tests whether any individual gene's weighted expression difference contributes more than expected by chance. However, in many biological scenarios, the signal is diffuse: no single gene dominates, but many contribute weakly and collectively. In such cases, no individual gene may reach statistical significance. This is what happens in the example dataset. If there's still interest in identifying and ranking genes by their relative contribution to the age shift (even if not statistically significant), the `alpha_gene` parameter can increased to return a sorted dataframe of genes based on their weighted contribution.
```r
# Compute the genes that are contributing to the significant(s) age shift(s).
baso_old_vs_young <- get_signGenes(baso_old_vs_young, alpha_comparisons = .05, alpha_genes = 0.1)

# Access the gene significance matrix for the comparison(s).
print(baso_old_vs_young$sign_genes$old_vs_young)
```
## References
1. Jaffe AE, Straub RE, Shin JH, Tao R, Gao Y, Collado-Torres L, Kam-Thong T, Xi HS, Quan J, Chen Q, Colantuoni C, Ulrich WS, Maher BJ, Deep-Soboslay A, BrainSeq Consortium, Cross AJ, Brandon NJ, Leek JT, Hyde TM, Kleinman JE, Weinberger DR. *Developmental and genetic regulation of the human cortex transcriptome illuminate schizophrenia pathogenesis.* Nat Neurosci. 2018; 21(8):1117-1125. [DOI: 10.1038/s41593-018-0197-y](https://doi.org/10.1038/s41593-018-0197-y)
