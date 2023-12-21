### CLEAN EVERYTHING AND START ###
cat("\014");rm(list=ls())

### ESSENTIAL PACKAGES ###
require(dplyr)
require(tidyr)
require(ggplot2)
require(cowplot)

### IMPORT FUNCTIONS ###
source("functions.R")

### CREATE 'results/' DIR ###
dir.create( file.path("results/plots"), recursive = TRUE, showWarnings = FALSE )

# DATA INFORMATION ----
list.files("data/")

# PHENOTYPE DATA ----
## Load complete phenotype data
  phenodata_complete <- read.csv("data/phenotype-complete.csv", header=TRUE)

## Mean in 3 last columns
  phenodata_complete$respuesta <- rowMeans(phenodata_complete[,5:7],na.rm=TRUE) 

## Delete those columns
  phenodata_complete <- phenodata_complete[,c(1:3,8)]

## Load phenotype of sequenced individuals 
  phenodata.df <- read.table("data/dmz119_phenotype.txt",header=TRUE,skip=2)

## Histogram: Hydrotropic Response histogram
  PhenoHist <- phenodata_complete %>% filter(respuesta > 0) %>%
    ggplot(aes(x=respuesta))+
      geom_histogram(alpha=0,color="black")+
      labs(
        x="Response (angle)",
        y="Count"
      )+
      theme_minimal()

  ## Save plot as PNG
  ggsave(
    file  = "results/plots/PhenoHist.png",
    plot  = PhenoHist,
    device = "png",
    width = 4, height = 4, dpi = 300
  )


# GENOTYPE DATA ----
# Read sequencing data
  snp_report.path <- "data/genotype_dmz119.txt"
  snp_report.file <- read.table(snp_report.path,sep="\t",header = TRUE,row.names = NULL)

# 5 individuals as an example
  head(snp_report.file[,c(7:12)])
  names(snp_report.file[,1:6])
## CONVERT TO HAPMAP ----
## 'snp2hmp()' imported from 'source("functions.R")'
  # source("scripts/snp2hmp.R")
  hapmap_report <- snp2hmp(file = snp_report.file[,-6])

# 5 individuals as an example.
  head(hapmap_report[,12:17])

# CLEAN DATA ----
## Show unique chromosome data 
  unique(hapmap_report$chrom)[1:20]

## Data proportion ----
## Chromosomes
  perc.chrom <- sum(grepl("_chromosome",snp_report.file$chrom))/length(snp_report.file$chrom)
  perc.chrom <- sprintf("%.2f%%",perc.chrom*100)

## Contigs
  perc.ctg <- sum(grepl("_ctg",snp_report.file$chrom))/length(snp_report.file$chrom)
  perc.ctg <- sprintf("%.2f%%",perc.ctg*100)

## Empty Data
  perc.rest <- length(which(snp_report.file == ""))/length(snp_report.file$chrom)
  perc.rest <- sprintf("%.2f%%",perc.rest*100)

## Total SNPs
  total_snps.all <- scales::comma( nrow(snp_report.file) ) 

## Print results
  cat("Total SNPs: ", total_snps.all, "\n")
  cat("Identified chromosome: ", perc.chrom, "\n")
  cat("Contigs: ", perc.ctg, "\n")
  cat("Empty data: ", perc.rest, "\n")

## Select every string with 'chromosome'
  hapmap_report.clean <- hapmap_report %>% filter( grepl('chromosome',chrom) )

## Delete 'chromosome'
  hapmap_report.clean$chrom <- gsub("_chromosome","",hapmap_report.clean$chrom)

# Assign a chromosomal number code
## TASSEL does not accept strings in 'chr' column.
## 11 = Mitocondrial Chrom. ; 12 = Chloroplast Chrom.
  hapmap_report.clean$chrom <- gsub("Mt", 11,hapmap_report.clean$chrom)
  hapmap_report.clean$chrom <- gsub("Pt", 12,hapmap_report.clean$chrom)

  hapmap_report.clean$chrom <- as.integer(hapmap_report.clean$chrom) 

## Total SNPs
  cat("Total de SNPs: ", scales::comma(nrow(hapmap_report.clean)), "para los cromosomas: [", unique(hapmap_report.clean$chrom), "]\n")

## Save HAPMAP 
  names(hapmap_report.clean)[1] <- "rs#"
  names(hapmap_report.clean)[6] <- "assembly#"

  head(hapmap_report.clean)[1:15]

  write.table(
   hapmap_report.clean,
   "results/dmz119_genotype.hmp.txt",
   sep = "\t", row.names = FALSE, quote = FALSE
  )

# PERFORM GWAS ----
# Site Summary (Complete)
  SiteSum.complete <- read.table("data/dmz119_sitesum.txt",
                                      sep="\t",header=TRUE)

# Site Summary (Filtered)
  SiteSum.filter <- read.table("data/dmz119_filtered-sitesum.txt",
                                      sep="\t",header=TRUE)

  psnps.loss <- ( nrow(SiteSum.complete) - nrow(SiteSum.filter) ) / nrow(SiteSum.complete)

  cat("SNPs before filter:",scales::comma(nrow(SiteSum.complete)),"\n",
      "SNPs after filter:",scales::comma(nrow(SiteSum.filter)),"\n",
      "SNPs lost:", sprintf("%.2f%%",psnps.loss*100),"\n")

# 'sitesummary_analysis()' imported from 'source("function.R")'
  SiteSum.loss <- sitesummary_analysis(
    SiteSum.complete,
    SiteSum.filter
  ) 
  names(SiteSum.loss)[1] <- "chrom"

## "Translate" again mitocondria and chloroplast codes
  SiteSum.loss$chrom <- gsub(11,"Mt",SiteSum.loss$chrom)
  SiteSum.loss$chrom <- gsub(12,"Pt",SiteSum.loss$chrom)

## Arrange chromosome order
  SiteSum.loss$chrom <- factor(SiteSum.loss$chrom,
                               levels = c(as.character(1:10),"Mt") )
  names(SiteSum.loss)[5:6] <- c("Perdidos","Conservados")

## Print summary table
  SiteSum.loss

# DATA VISUALIZATION ----
## Histograms: Site Summary SNP Count ----
## BEFORE applying MAF < 10% filter
  SiteSumComplete.hist <- SiteSum.complete %>% select(Chromosome,Minor.Allele.Frequency) %>% 
    ggplot( aes(x=Minor.Allele.Frequency) )+
    geom_histogram(color="black",fill="white",bins = 25)+
    scale_x_continuous(limits = c(0,0.5), breaks = seq(0,0.5,0.05))+
    labs(x="")+
    theme_minimal()+
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

## AFTER applying MAF < 10% filter
  SiteSumFiltered.hist <- SiteSum.filter %>% select(Chromosome,Minor.Allele.Frequency) %>% 
    ggplot( aes(x=Minor.Allele.Frequency) )+
    geom_histogram(color="black",fill="white", bins = 25)+
    scale_x_continuous(limits = c(0,0.5), breaks = seq(0,0.5,0.05))+
    labs(x="MAF")+
    scale_x_continuous(limits = c(0,0.5), breaks = seq(0,0.5,0.05))+
    theme_minimal()+
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

## Join both plots in a grid
  SiteSumMAF <- plot_grid(
    SiteSumComplete.hist,
    SiteSumFiltered.hist,
    ncol=1,nrow=2,labels=c('A','B')
  )

# Save plot as a PNG 
  ggsave(
    file = "results/plots/SiteSumLoss.png",
    plot     = SiteSumMAF,
    device = "png",
    width = 4, height = 4, dpi = 300
  )

## Bar Chart: SNPs loss per chromosome afer MAF filter ----
  SiteSumChrom <- SiteSum.loss %>% filter(chrom != "NA") %>% 
  select(chrom, Perdidos, Conservados) %>% 
  gather("Status", "Proportion", -chrom) %>% 
    ggplot( aes(x=chrom, y=Proportion, fill=factor(Status,levels=c('Perdidos','Conservados')) )  )+
      geom_bar(stat='identity')+
      geom_text(
      aes( label = sprintf("%.1f%%",Proportion*100) ),
        position = position_stack(vjust = 0.5), color="white", size = 2, angle = 90
      )+
      labs(
      fill = "Status",
      x = "Cromosoma"
      )+
      scale_y_continuous( labels = scales::percent, expand = c(0,0) )+
      theme_minimal()+
      theme(
        panel.grid.major.x = element_blank()
      )

# Save SiteSumChrom plot as PNG
  ggsave(
    file = "results/plots/SiteSumChrom.png",
    plot     = SiteSumChrom,
    device = "png",
    width = 4, height = 4, dpi = 300
  )

# GWAS ANALYSIS ----
  GWAS_results <- read.table("data/snp-dmz119_gwas-mlm.txt",
                             sep="\t",header=TRUE)

## "Translate" back chromosome numbers top
## mitocondrial and chloroplast tags

  GWAS_results$Chr <- gsub(11,"Mt",GWAS_results$Chr)

## Note: MAF filter did not return Chrom12 SNPs,
## that's why its omited.

  # GWAS_results$Chr <- gsub(12,"Pt",GWAS_results$Chr)

## Select Top 20 SNPs according to
## logarithmic ( -log10(p) ) transformation

  GWAS_Best <- GWAS_results %>% 
    filter(Chr != "NA" & Trait == "RH") %>% 
    top_n(20, -log10(p))

## Write table with best SNPs
  write.table(
    GWAS_Best,
    "results/GWASBestResults.txt",
    sep = "\t", row.names = FALSE
  )

## Arrange chromosomes
  GWAS_results$Chr <- as.factor(GWAS_results$Chr)
  levels(GWAS_results$Chr) <- c( as.character(1:11) )

## Bar Chart: Total SNPs per chromosome
  ChromBar.plot <- GWAS_results %>% group_by(Chr) %>% 
  summarise( SNPs = n() ) %>% 
  mutate( Chr = as.factor(Chr) ) %>%
    ggplot( aes(x=Chr,y=SNPs) )+
    geom_bar( stat='identity', fill="#000066", alpha=0.8 )+
    geom_text( aes( label = scales::comma(SNPs) ),
      angle = 90, color = "white", position = position_stack(vjust = 0.9)
    )+
    scale_y_continuous(expand=c(0,0),
    labels = scales::comma)+
    labs(
      x = "Cromosoma",
      y = "SNPs"
    )+
    theme_minimal()+
    theme(
      panel.grid = element_blank()
    )

## Save plot as PNG
  ggsave(
    plot = ChromBar.plot,
    file = "results/plots/ChromBar.png",
    device="png",
    width=6, height=3, dpi = 300, units="in" 
  )

# GWAS RESULTS ---- 
## 'expected_quantiles()' imported from 'source("functions.R")'
qq <- expected_quantiles( na.omit(GWAS_results$p) )

# QQ-Plot ----
  QQPlot.plot <- qq %>% 
  ggplot(
    aes(x = -log10(expected), y = -log10(actual) )
  )+
  geom_point(color="red",size=1,alpha=0.45)+
  geom_text(
    label = "-log10(p) = 2.48",
    x = 1.5, y = 4.5, color ="blue"
  )+
  geom_abline(slope=1)+
  geom_vline(xintercept = 2.48,color="blue",linetype="dashed",size=1)+
  scale_y_continuous(
    limits = c(0,5),
    breaks = seq(0,5,0.5)
  )+
  scale_x_continuous(
    limits = c(0,5),
    breaks = seq(0,5,0.5)
  )+
  labs(
    x = "-log10(p-value) esperado",
    y = "-log10(p-value) real"
  )+
  theme_minimal()+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90)
  )

## Save as PNG
  ggsave(
    plot = QQPlot.plot,
    file = "results/plots/QQPlot.png",
    device="png",
    width = 4, height = 3, dpi = 150, units = "in" 
  )

  GWAS_results$Chr <- as.factor(GWAS_results$Chr)
  levels(GWAS_results$Chr) <- as.character(1:11)

  ManhattanPlot.plot <- 
    GWAS_results %>% filter(Chr != "NA" & Trait == "rh") %>% 
    ggplot( aes( x=Pos, y=-log10(p), color=Chr )  )+
      geom_point(alpha=0.45, size = 0.8)+
      # Top 20 SNPs
      geom_point( 
        data = . %>%
        # group_by(Chr) %>%
        top_n(10, -log10(p)),
        aes(x=Pos,y=-log10(p) ),
        color = "red", size = 2
      )+
      geom_hline(
        yintercept = 2.48,
        color="red",
        linetype="dashed"
      )+
      facet_grid( cols = vars(Chr), scales="free_x", switch = "both" )+
      guides( color='none' )+
      scale_y_continuous(
        limits = c(0,4),
        expand = c(0,0)
      )+
      theme_minimal()+
      labs(
        x = "Cromosoma"
      )+
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing.x = unit(0,"lines")
      )

## Save plot as PNG
  ggsave(
    file   = "results/plots/ManhattanPlot.png",
    plot   = ManhattanPlot.plot,
    device = "png",
    width = 8, height = 4, units = "in", dpi = 300
  )
