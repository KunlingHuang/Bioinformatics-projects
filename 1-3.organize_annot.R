library(dplyr)
library(tidyverse)
library(optparse)
library(parallel)
library(data.table)
options(stringsAsFactors=FALSE)
### 1. Prepare the sumstat for plotting
res = 40000
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
  mutate(out.kl = paste0("../data/H-MAGMA/input/annot_kl/", region, "/chr", Chr, "/", Gene, ".prenatal.", region, ".res", res, ".annot.txt"),
         out.hmagma = paste0("../data/H-MAGMA/input/annot_hmagma/", region, "/chr", Chr, "/", Gene, ".prenatal.", region, ".res", res, ".annot.txt"))


file.create("../data/H-MAGMA/input/annot_kl/full_table/annot.kl.full.FBD.res40000.txt")
file.create("../data/H-MAGMA/input/annot_kl/full_table/annot.kl.full.FBP.res40000.txt")
file.create("../data/H-MAGMA/input/annot_hmagma/full_table/annot.hmagma.full.FBD.res40000.txt")
file.create("../data/H-MAGMA/input/annot_hmagma/full_table/annot.hmagma.full.FBP.res40000.txt")


annot1 = list()
annot2 = list()

stack_annot = function(idx){
  
  TSS.info = info[idx]
  
  if(file.exists(TSS.info$out.kl)){
    annot_kl.gene = fread(TSS.info$out.kl, header = F) %>% 
      mutate(gene.info = paste0(TSS.info$Chr, ":", TSS.info$Start, ":", TSS.info$End)) %>% 
      select(V1, gene.info, everything()) 
    # write.table(annot_kl[[idx]], 
    #             file = paste0("../data/H-MAGMA/input/annot_kl/full_table/annot.kl.full.", TSS.info$region, ".res40000.txt"), 
    #             col.names = F, row.names = F, sep = "\t", quote = F, append = TRUE)
  } else{
    annot_kl.gene = NA
    
  }
  
  if(file.exists(TSS.info$out.hmagma)){
    annot_hmagma.gene = fread(TSS.info$out.hmagma, header = F) %>%
      mutate(gene.info = paste0(TSS.info$Chr, ":", TSS.info$Start, ":", TSS.info$End)) %>%
      select(V1, gene.info, everything())
    # write.table(annot_hmagma[[idx]],
    #             file = paste0("../data/H-MAGMA/input/annot_hmagma/full_table/annot.hmagma.full.", TSS.info$region, ".res40000.txt"),
    #             col.names = F, row.names = F, sep = "\t", quote = F, append = TRUE)
    #
  } else{
    annot_hmagma.gene = NA
  }
   
  # annot1[[idx]] = annot_kl.gene
  # annot2[[idx]] = annot_hmagma.gene
  # print(annot1)
  return()
}


# lapply(1:nrow(info), stack_annot)

tmp = mclapply(1:3, stack_annot, mc.cores = 20)


tmp
# lapply(1:nrow(info), stack_annot)

# mclapply(1:10, hmagma_annot, mc.cores = 60)




