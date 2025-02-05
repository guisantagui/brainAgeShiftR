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
