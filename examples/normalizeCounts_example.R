# Normalize counts using computed quantile means
obj <- normalizeCounts(obj, useTrainMeans = FALSE)

# Normalize counts using training quantile means
obj <- normalizeCounts(obj, useTrainMeans = TRUE)
