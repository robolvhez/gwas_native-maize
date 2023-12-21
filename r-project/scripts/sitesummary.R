
# ! IMPORTANTE
# Este script está diseñado para ser utilizado en otros análisis
# Favor de implementar la función en el análisis como:
#
# source(sitesummary.R)

### GATHER SITESUMMARY ANALYSIS ###
# Juntar ambos resumenes de datos
## Summary sin filtrar
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
