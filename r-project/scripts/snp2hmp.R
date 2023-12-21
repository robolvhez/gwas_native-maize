# ! IMPORTANTE
# Este script está diseñado para ser utilizado en otros análisis
# Favor de implementar la función en el análisis como:
#
# source(snp2hapmap.R)

### SNP2HAPMAP ###
snp2hapmap <- function(file){
  snp.report <- file # leer la tabla
  snp.report <- snp.report[snp.report$snp!="",] # Eliminar los renglones vacios
  
  ## Sustituir 
  alleles <- data.frame(
    allele      = gsub("\\d+:(\\w)>(\\w)", "\\1", snp.report$snp),
    sustitution = gsub("\\d+:(\\w)>(\\w)", "\\2", snp.report$snp)
  )
  alleles$snp     <- paste0(alleles$allele,"/",alleles$sustitution)
  
  alleles.results <- cbind(alleles,snp.report[,7:ncol(snp.report)])
  
  # Sustituir valores numericos por alelos
  alleles.results[,4:ncol(alleles.results)] <- ifelse( # Sentencia condicional para cambiar valores.
    alleles.results[,4:ncol(alleles.results)] == 1, # homocigotos
    strrep(alleles.results$sustitution,2),
    ifelse(
      alleles.results[,4:ncol(alleles.results)] == 2, # heterocigotos
      paste0(alleles.results$allele,alleles.results$sustitution),
      ifelse(
        alleles.results[,4:ncol(alleles.results)] == 0,
        strrep(alleles.results$allele,2),
        NA
      )
    )
  )
  
  alleles.results <- alleles.results[,4:ncol(alleles.results)]
  head(alleles.results)
  
  hapmap <- data.frame(
    rs        = snp.report$id,
    alleles   = alleles$snp,          
    chrom     = snp.report$chrom,     # Número de cromosoma (si hay letras hay que quitarlas, o filtrar los contigs)
    pos       = snp.report$chrom_pos, # Posicion en cromosoma
    strand    = NA,
    assembly  = NA,
    center    = NA,
    protLSID  = NA,
    assayLSID = NA,
    panel     = NA,
    QCcode    = NA
  )
  hapmap <- cbind(hapmap,alleles.results)
  return(hapmap)
}