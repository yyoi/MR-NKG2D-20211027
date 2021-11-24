#! /usr/bin/R

library(devtools)
library(TwoSampleMR)
library(tidyverse)
library(data.table)

setwd("./")
genes=c("MICA","MICB")
VALs=c("1E51.0Z.CHB","5A02.0.Graves")

# Analysis 
## read exposure
for(j in seq(1, length(genes))){ 
gene <- genes[j]
  exp_dat <- read_exposure_data(
    filename = paste0("./eQTLs/instruments.",gene,".forMR.txt"),
    phenotype_col = "Phenotype.e",
    snp_col = "ID",
    beta_col = "BETA.e",
    se_col = "SE.e",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    eaf_col = "FREQ.e",
    pval_col = "P.e",
    samplesize_col = "N.e",
    units_col = "units.e") 
  
for(i in seq(1, length(VALs))){ 
  val <- VALs[i]  
  ## Outocome
  outcome_dat <- read_outcome_data(
    filename = paste0("./Outcomes/filtered_outcome_",gene,"_",val,".un.txt"),
    phenotype_col = "Phenotype",
    snp_col = "ID",
    beta_col = "BETA.o",
    se_col = "SE.o",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    eaf_col = "FREQ.o",
    pval_col = "P.o",
    samplesize_col = "N.o",
    units_col = "units.o"
    )

  ## Harmonize 
  dat <- harmonise_data(
    exposure_dat = exp_dat,
    outcome_dat = outcome_dat)
  dat <- dat
  datsd <- standardise_units(dat)

  ## find outliners
  pre <- run_mr_presso(datsd, NbDistribution = 5000, SignifThreshold = 0.05)
  if (is.na(pre[[1]][1]$`Main MR results`$`P-value`[2])){
    datsd1 <- datsd
    print(datsd1)
  } else {
    rmins <- c(pre[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
    rmsnp <- dat$SNP[rmins]
    datsd1 <- datsd %>% filter(!(SNP %in% rmsnp))
    print(datsd %>% filter((SNP %in% rmsnp)))
  }
  
  ## MR
  param <- default_parameters()
  param$nboot <- 5000
  res <- mr(datsd1, parameters = param)
  
  print(res2 <- generate_odds_ratios(res))
  filename=paste0("./mr/",gene,"_",val,".mr.or.sd.rmsnp.txt")
  write.table(res2, filename, row.names = F, quote = F)
  
  print(o <- mr_pleiotropy_test(datsd1))
  filename=paste0("./mr/",gene,"_",val,".pleio.sd.rmsnp.txt")
  write.table(o, filename, row.names = F, quote = F)
  
}}
