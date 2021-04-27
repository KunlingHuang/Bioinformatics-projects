organize.gwas_hic = function(Chr, Gene, hic, gwas.path, gwas.colnames, Region, upstream, downstream, q.threshold, gwas.threshold, TSS.info){
  ### 1. Read in Hi-C and GWAS data
  
  # 1-2. Extract the Hi-C interactions interacting with the promoter region
  bin.gene = fread(hic) %>% 
    as.data.frame %>% 
    select(chr1, fragmentMid1, fragmentMid2, contactCount:`q-value`) %>% 
    filter(fragmentMid1 == TSS.info$bin | fragmentMid2 == TSS.info$bin) %>% 
    mutate(fragment = case_when(
      fragmentMid1 == TSS.info$bin ~ fragmentMid2,
      fragmentMid2 == TSS.info$bin ~ fragmentMid1
    ), direction = TSS$direction) %>% 
    mutate(include = case_when(
      (direction == "+") & (fragment <= TSS.info$bin + upstream * 1e6) & (fragment >= TSS.info$bin - downstream * 1e6) ~ T,
      (direction == "-") & (fragment <= TSS.info$bin + downstream * 1e6) & (fragment >= TSS.info$bin - upstream * 1e6) ~ T,
      TRUE ~ FALSE
    )) %>% 
    select(fragment, contactCount:`q-value`) 
  
  if(nrow(bin.gene) == 0){
    gwas = data.frame()
    return(gwas)
  }
  
  # 1-3. Read in GWAS summary statistics
  gwas = fread(gwas.path) %>% 
    as.data.frame %>% 
    `colnames<-`(strsplit(gwas.colnames, ",")[[1]]) %>% 
    mutate(Z = BETA/SE, gwas.bin = (2 * floor(as.numeric(bp)/10000) + 1) * 5000, direction = TSS.info$direction) %>% 
    mutate(include = case_when(
      (direction == "+") & (gwas.bin <= TSS.info$bin + upstream * 1e6) & (gwas.bin >= TSS.info$bin - downstream * 1e6) ~ T,
      (direction == "-") & (gwas.bin <= TSS.info$bin + downstream * 1e6) & (gwas.bin >= TSS.info$bin - upstream * 1e6) ~ T,
      TRUE ~ FALSE
    )) %>% 
    filter(include == T) %>% 
    left_join(bin.gene, by = c("gwas.bin" = "fragment")) %>% 
    mutate(contactCount = tidyr::replace_na(contactCount, 0), 
           `p-value` = tidyr::replace_na(`p-value`, 1),
           `q-value` = tidyr::replace_na(`q-value`, 1)) %>% 
    mutate(contactCount = ifelse(gwas.bin == TSS.info$bin, max(.$contactCount), contactCount), 
           `p-value` = ifelse(gwas.bin == TSS.info$bin, min(.$`p-value`), `p-value`),
           `q-value` = ifelse(gwas.bin == TSS.info$bin, min(.$`q-value`), `q-value`)) %>% 
    mutate(contact.binary = ifelse(`q-value` <= q.threshold, 1, 0)) %>% 
    select(chr, bp, snp, A1, A2, Z, P, gwas.bin:contact.binary)
  
  return(gwas)
  
}

plot.hic_gwas = function(gwas, col.x, col.y, title){
  gwas = gwas %>% 
    select(as.name(col.x), as.name(col.y)) %>% 
    `colnames<-`(c("V1", "V2")) 
  gg = ggplot(gwas, aes(x = factor(V1), y = V2)) +
    geom_violin(fill = "#f3e6e3", draw_quantiles = c(0.25, 0.5, 0.75)) +
    # geom_smooth(method = "lm", col = "#1e5f74") +
    geom_text(stat="count", aes(label = paste0("N = ",..count..)), y = 1.1*max(gwas$V2)) +
    ggtitle(title) +
    ylim(c(min(gwas$V2), 1.15*max(gwas$V2))) +
    xlab(col.x) +
    ylab(col.y) +
    theme_bw() + 
    theme(panel.grid.major.y = element_blank(), 
          panel.grid.minor=element_blank(), 
          axis.line.y = element_line(color = "#2b2b2b", size = 0.15),
          plot.subtitle=element_text(margin = margin(b = 10)),
          plot.title=element_text(size = 8), 
          plot.caption=element_text(size = 8, margin = margin(t=10)),
          legend.position="bottom"
    )
  return(gg)
}


plot.hic_gwas.distance = function(gwas, col.x, col.y, title, min.abs_distance.promoter){
  gg = gwas %>% 
    select(as.name(col.x), as.name(col.y), `q-value`) %>% 
    `colnames<-`(c("V1", "V2", "q")) %>% 
    mutate(q = -log(q, 10)) %>% 
    arrange(q) %>% 
    ggplot(., aes(x = V1, y = V2)) +
    geom_point(aes(fill = q), shape = 21, aplha = 0.4) +
    geom_vline(xintercept = min.abs_distance.promoter) +
    scale_fill_continuous(type = "viridis") + 
    scale_x_continuous(breaks=seq(-1.01 * 1e6, 1.01 * 1e6, 1e4)) + 
    ggtitle(title) +
    xlab(col.x) +
    ylab(col.y) +
    theme_bw() + 
    theme(panel.grid.major.y = element_blank(), 
          panel.grid.minor=element_blank(), 
          axis.text.x = element_blank(),
          axis.line.y = element_line(color = "#2b2b2b", size = 0.15),
          plot.subtitle=element_text(margin = margin(b = 10)),
          plot.title=element_text(face = "bold"), 
          plot.caption=element_text(size = 8, margin = margin(t=10)),
          legend.position="bottom"
    )
  return(gg)
}
