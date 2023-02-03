# Load the imputeTS library
library(imputeTS)

# Impute the missing values using na_kalman
imputed_matrix <- na_kalman(norm)

write.table(imputed_matrix, file='formicicum_proteomics/results/quant_norm.tsv', sep = "\t", row.names = TRUE,
            col.names = TRUE, quote=FALSE)
