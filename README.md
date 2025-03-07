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
The example data consists on transcriptomic data from neurons trans-differentiated from fibroblasts originating from either healthy controls (`HC`) or late-onset Alzheimer's disease (`LOAD`), generated by Sun et al. in 2024 [1].
```r
data(example_data)
```
### Create a `brainAgeShiftObj`
```r
baso_load_vs_hc <- create_brainAgeShiftObj(counts = example_counts,
                                           metadata = example_metadata,
                                           variable = "pheno",
                                           comparisons = c("HC", "LOAD"))
```
### Normalize the counts data
```r
baso_load_vs_hc <- normalizeCounts(baso_load_vs_hc)
```
### Compute transcriptomic age
```r
baso_load_vs_hc <- predictAge(baso_load_vs_hc)
```
### Compute significance of the transcriptomic age shift of the comparison.
```r
# Compute significance.
baso_load_vs_hc <- do_signTest(baso_load_vs_hc, adjust_method = "BH")

# Access results of the significance test(s).
baso_load_vs_hc$stats
```
### Get the genes that are significantly contributing to the shift
```r
# Compute the genes that are contributing to the significant(s) age shift(s).
baso_load_vs_hc <- get_signGenes(baso_load_vs_hc, alpha_comparisons = .05)

# Access the gene significance matrix for the comparison(s).
print(baso_load_vs_hc$sign_genes$LOAD_vs_HC)
```
## References
1. Sun Z, Kwon JL, Ren Y, Chen S, Walker CK, Lu X, Cates K, Karahan H, Sviben S, Fitzpatrick JAJ, Valdez C, Houlden H, Karch CM, Bateman RJ, Sato C, Mennerick SJ, Diamond MI, Kim J, Tanzi RE, Holtzman DM, Yoo AS. *Modeling late-onset Alzheimer's disease neuropathology via direct neuronal reprogramming.* Science. 2024;385(6708):adl2992. [DOI: 10.1126/science.adl2992](https://doi.org/10.1126/science.adl2992).
