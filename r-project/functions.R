# ! IMPORTANT
# These scripts are coded to run on main.Rmd,
# so please, be sure to add `source("functions.R")` at the beginning
# of the notebook so everything runs smooth as butter.

### SNP2HMP ###
snp2hmp <- function(file){
# Read table
  snp.report <- file 
# Ignore empty rows.
  snp.report <- snp.report[snp.report$snp!="",]
  
# Subtitute
  ## Obtain the information of changes in the allele.
  ## The notation for each change is '21:A>G'
  ## 
  ## So, we create an 'alleles' dataframe containing
  ## the original allele (1) and its mutation (2).

  alleles <- data.frame(
    allele      = gsub("\\d+:(\\w)>(\\w)", "\\1", snp.report$snp),
    sustitution = gsub("\\d+:(\\w)>(\\w)", "\\2", snp.report$snp)
  )

  ## HAPMAP files requiere the notation for changes
  ## as a 'A/G' notation.
  ## 
  ## These new notation is now on the 'allele' dataframe.

  alleles$snp <- paste0(alleles$allele,"/",alleles$sustitution)

  ## Adding samples cygotes in a different table.
  cygotes <- snp.report[,7:ncol(snp.report)]

  ## Make a joined table
  results <- cbind(alleles,cygotes)

# Sustituir valores numericos por alelos
  results[,4:ncol(results)] =
    ifelse(
      ## Homocygotes
      results[,4:ncol(results)] == 1,
      strrep(results$sustitution,2),

        ifelse(
        ## Heterocygotes
        results[,4:ncol(results)] == 2,
        paste0(results$allele,results$sustitution),

          ifelse(
            ## No change
            results[,4:ncol(results)] == 0,
            strrep(results$allele,2),
            ## Empty data
            NA
        )
      )
    )
    
  results <- results[,4:ncol(results)]
  
  hapmap <- data.frame(
      rs        = snp.report$id,
      alleles   = alleles$snp,          
      chrom     = snp.report$chrom,
      pos       = snp.report$chrom_pos,
      strand    = NA,
      assembly  = NA,
      center    = NA,
      protLSID  = NA,
      assayLSID = NA,
      panelLDID = NA,
      QCcode    = NA
    )

    hapmap <- cbind(hapmap,results)

  return(hapmap)
}

### SITESUMMARY ANALYSIS ###

sitesummary_analysis <- function(complete,filter){
  # Juntar ambos resumenes de datos
  ## Summary sin filtrar
  complete <- complete %>% 
    select( # Seleccionar solamente las columnas que nos interesan, no todas son necesarias. 
      Site.Name,Chromosome,Physical.Position,Minor.Allele.Frequency
    ) %>% 
    mutate(summary_type = as.factor("Completo"))
  
  ## Summary filtrado
  filter <- filter %>% 
    select(
      Site.Name,Chromosome,Physical.Position,Minor.Allele.Frequency
    ) %>% 
    mutate(summary_type = as.factor("Filtrado"))
  
  ## Combinar (bind)
  binded <- rbind(complete,filter)
  
  # ¿Cuántos datos perdí en cada cromosoma?
  loss <- binded %>% group_by(summary_type, Chromosome) %>% 
    summarise(snp = n()) %>% 
    spread(key = summary_type, value = snp) %>% 
    mutate(
      SNP.Perdidos     = Completo - Filtrado,
      Prop.Perdidos    = (Completo - Filtrado)/Completo,
      Prop.Conservados = 1 - Prop.Perdidos
    )
  
 return(loss) 
}

### EXPECTED QUANTILES ###
expected_quantiles <- function(pvalues){
  n = length(pvalues)
  actual_quantiles = sort(pvalues)
  expected_quantiles = seq_along(pvalues)/n 
  data.frame(expected = expected_quantiles, actual = actual_quantiles)
}

