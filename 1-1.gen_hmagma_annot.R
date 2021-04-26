library(dplyr)
library(data.table)
library(tidyverse)
library(optparse)
library(ggplot2)
library(ggrepel)
library(qqman)
library(parallel)
library(data.table)
options(scipen = 99)

### 0. Read in GWAS and eQTL files
# gwas.path = "/z/Comp/lu_group/Members/khuang82/Hi-C_GWAS/results/GWAS_sumstat_H-MAGMA_compatible/Huang_ASD_2019.1KG_rsid_update.add_Z.txt"
# gwas = fread(gwas.path)
# gwas.chr = lapply(1:22, FUN = function(chr) gwas %>% filter(CHR == chr)) 

### 1. Prepare the sumstat for plotting
res = 10000
tss = fread("./condor_submit_txt/gene_chr.location.txt") %>% 
  `colnames<-`(c("Chr", "Gene", "Start", "End", "Direction")) %>% 
  mutate(promoter = ifelse(Direction == "-", End + 2000, Start - 2000)) %>% 
  mutate(bin = (2 * floor(promoter/res) + 1) * res/2)

info = tss %>% 
  mutate(promoter = ifelse(Direction == "-", End + 2000, Start - 2000)) %>% 
  mutate(bin = (2 * floor(promoter/res) + 1) * res/2) %>% 
  mutate(key = 0) %>% 
  left_join(data.frame(region = c("FBD", "FBP"), key = 0), by = "key") %>% 
  mutate(hic = paste0("../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fit_hi-c_unzipped/", region,"/", region, "_chr", Chr, "_FitHiC.spline_pass2.res", res, ".significances.txt")) %>% 
  mutate(out = paste0("../data/H-MAGMA/input/annot_KL/", region, "/chr", Chr, "/", Gene, ".prenatal.", region, "res", res, ".annot.txt"))


gene.body = fread("/z/Comp/lu_group/Members/khuang82/Resources/Public_database/Gencode/human/release_26lift37/gencode.v26lift37.annotation.gene_exon_transcript.curated.txt") 
### 2. Stack the summary statistics together
source("./Rfiles/organize_GWAS_hiC.R")

hmagma_annot = function(idx, upstream = 1, downstream = 1, upstream.hmagma = 0.035, downstream.hmagma = 0.01){
  #FIXME
  TSS.info = info[idx, ]
  upstream = 1
  downstream = 1
  upstream.hmagma = 0.035
  downstream.hmagma = 0.01
  # ---1. Organize Hi-C data
  bin.gene = fread(TSS.info$hic) %>% 
    as.data.frame %>% 
    select(chr1, fragmentMid1, fragmentMid2, contactCount:`q-value`) %>% 
    filter(fragmentMid1 == TSS.info$bin | fragmentMid2 == TSS.info$bin) %>% 
    mutate(fragment = case_when(
      fragmentMid1 == TSS.info$bin ~ fragmentMid2,
      fragmentMid2 == TSS.info$bin ~ fragmentMid1), 
      direction = TSS.info$Direction) %>% 
    mutate(include = case_when(
      (direction == "+") & (fragment <= TSS.info$bin + downstream * 1e6) & (fragment >= TSS.info$bin - upstream * 1e6) ~ T,
      (direction == "-") & (fragment <= TSS.info$bin + upstream * 1e6) & (fragment >= TSS.info$bin - downstream * 1e6) ~ T,
      TRUE ~ FALSE
    )) %>% 
    filter(include) %>% 
    select(fragment, contactCount:`q-value`) 
  
  # ---2. Organize exon information
  
  exon = gene.body %>% filter(name == TSS.info$Gene) %>% 
    filter(Category == "exon")
  
  # ---3. Organize GWAS summary statistics
  
  gwas.gene = gwas.chr[[TSS.info$Chr]] %>% 
    mutate(include = case_when(
      TSS.info$Direction == "+" & (BP <= TSS.info$promoter + downstream * 1e6) & (BP >= TSS.info$promoter - upstream * 1e6) ~ T,
      TSS.info$Direction == "-" & (BP <= TSS.info$promoter + upstream * 1e6) & (BP >= TSS.info$promoter - downstream * 1e6) ~ T,
      TRUE ~ F
    )) %>% 
    filter(include) %>% 
    mutate(Gene = TSS.info$Gene, promoter = TSS.info$promoter, promoter.bin = TSS.info$bin) %>% 
    mutate(bin = (2 * floor(BP/res) + 1) * res/2) %>% 
    select(Gene, CHR:ID, bin, promoter, promoter.bin, everything()) %>% 
    
    # annotate H-MAGMA SNPs
    
    mutate(hmagma = case_when(
      TSS.info$Direction == "+" & (BP <= TSS.info$promoter + downstream.hmagma * 1e6) & (BP >= TSS.info$promoter - upstream.hmagma * 1e6) ~ 1,
      TSS.info$Direction == "-" & (BP <= TSS.info$promoter + upstream.hmagma * 1e6) & (BP >= TSS.info$promoter - downstream.hmagma * 1e6) ~ 1,
      TRUE ~ 0
    )) %>% 
    
    # add enhancer information
    
    left_join(bin.gene, by = c("bin" = "fragment")) %>%
    mutate(enhancer = ifelse(`q-value` > 0.01 | is.na(`q-value`), 0, 1))
    
  # ---4. Add exon information
    
   exon.snps = gwas.gene %>% 
     left_join(exon, by = c("Gene" = "name")) %>% 
     filter(BP >= start & BP <= end)
   
   gwas.gene.exon = gwas.gene %>% 
     mutate(exon = ifelse(SNP %in% exon.snps$SNP, 1, 0)) %>% 
     select(Gene:promoter.bin, A1, A2, Z, P, hmagma, `q-value`:exon)
   
      
   hmagma.annot = gwas.gene.exon %>% 
     mutate(promoter.region = ifelse(abs(BP -  promoter.bin) <= 5000, 1, 0)) %>% 
     filter(hmagma == 1 & (enhancer == 1 | exon == 1 | promoter == 1)) %>% 
     select(Gene, SNP) 
   
   if(nrow(hmagma.annot) != 0){
     hmagma.annot.final = c(TSS.info$Gene, hmagma.annot$SNP)
     write.table(t(hmagma.annot.final), file = TSS.info$out, col.names = F, row.names = F, 
                 sep = "\t", quote = F)
   }
     
   
   
  
}

mclapply(1:nrow(info), hmagma_annot, mc.cores = 64)

# mclapply(1:10, hmagma_annot, mc.cores = 60)




