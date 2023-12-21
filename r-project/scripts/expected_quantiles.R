# ! IMPORTANTE
# Este script está diseñado para ser utilizado en otros análisis
# Favor de implementar la función en el análisis como:
#
# source(snp2hapmap.R)

### Generar QQ-Plot - GWAS ###
expected_quantiles <- function(pvalues){
  n = length(pvalues)
  actual_quantiles = sort(pvalues)
  expected_quantiles = seq_along(pvalues)/n 
  data.frame(expected = expected_quantiles, actual = actual_quantiles)
}