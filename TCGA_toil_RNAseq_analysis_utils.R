#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TCGA survival analysis using UCSC TOIL dastaset  #####
# 2021-10-13 Hwang Seok Young 
# 
# Re-write functions using UCSC XENA Toil-RNAseq 
#
# 
# Functions
# 1. ge_tcga : get gene expression profiles of target gene 
# 2. ge_survival : survival analysis of target_gene, target_cancer pair at given expression threshold
# 3. ge_survival_allcancer : survival analysis of target gene for all cancers, output : list 
# 4. ge_survival_allcancer_genelists : survival analysis of genelist for all cancers, output : list 
# 5. survivalgraph2ppt : survival graphs to pptx (input : fit ,output : pptx file )
# 6. fit2dataframe : fit data to ordered dataframe by log rank p value 
#
# TO-DO LIST
# 1. Cox p value calculate
# 2. Auto cut threshold 
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# load libraries
library(tidyverse)
library(dplyr)
library(data.table)
library(survival)
library(survminer)
library(progress)
library(stringr)
source('https://raw.githubusercontent.com/syhwang-snu/tcga-utils/main/graph_ppt_utils.R')


# path for TOIL dataset 
toil_path <- '~/Documents/TCGA-survival-ASC/toil-data'

# Load TCGA dataset 

tcga_names <- strsplit('ACC	Adrenocortical Cancer
BLCA	Bladder Urothelial Carcinoma
BRCA	Breast invasive carcinoma
CESC	Cervical & Endocervical Cancer
CHOL	Cholangiocarcinoma
COAD	Colon adenocarcinoma
DLBC	Diffuse Large B-Cell Lymphoma
ESCA	Esophageal carcinoma
GBM	Glioblastoma multiforme
HNSC	Head & Neck Squamous Cell Carcinoma
KICH	Kidney Chromophobe
KIRC	Kidney Clear Cell Carcinoma
KIRP	Kidney Papillary Cell Carcinoma
LAML	Acute Myeloid Leukemia
LGG	Brain Lower Grade Glioma
LIHC	Liver hepatocellular carcinoma
LUAD	Lung adenocarcinoma
LUSC	Lung squamous cell carcinoma
MESO	Mesothelioma
OV	Ovarian serous cystadenocarcinoma
PAAD	Pancreatic adenocarcinoma
PCPG	Pheochromocytoma & Paraganglioma
PRAD	Prostate adenocarcinoma
READ	Rectum adenocarcinoma
SARC	Sarcoma
SKCM	Skin Cutaneous Melanoma
STAD	Stomach adenocarcinoma
TGCT	Testicular Germ Cell Tumor
THCA	Thyroid carcinoma
THYM	Thymoma
UCEC	Uterine Corpus Endometrioid Carcinoma
UCS	Uterine Carcinosarcoma
UVM	Uveal Melanoma',split = "\n")[[1]]

tcga_symbols <- trimws(sapply(strsplit(tcga_names, split = "\t"), "[[", 1))
tcga_fullnames <- sapply(strsplit(tcga_names, split = "\t"), "[[", 2)
tcga_fullnames <- str_to_title(tcga_fullnames)

tcga_types <- cbind(tcga_symbols, tcga_fullnames) %>% as_tibble()
tcga_types$tcga_symbols <- as.factor(tcga_types$tcga_symbols)
annotation <- read.csv('~/Documents/TCGA-survival-ASC/toil-data/gencode.v23.reorder.csv')
tcga_genes <- annotation$gene


##### preprocessing  codes #########################################################################################
# # reordering gencodes
# gencode <- fread('~/Documents/TCGA-survival-ASC/toil-data/probeMap%2Fgencode.v23.annotation.gene.probemap') %>% as.data.frame()
# sample_exps_order <- fread('~/Documents/TCGA-survival-ASC/toil-data/tcga_RSEM_gene_tpm.gz', sep='\t',select = 1)
# all(gencode$id %in% sample_exps_order$sample)
# all(gencode$id == sample_exps_order$sample)
# gencode <- gencode[order(match(gencode$id,sample_exps_order$sample)),]
# row.names(gencode) <- 1:nrow(gencode)
# write.csv(gencode, file = '~/Documents/TCGA-survival-ASC/toil-data/gencode.v23.reorder.csv', row.names = FALSE)


#### Load datasets #####

# After get data and tcga_type and annotation
# Data from USCS XENA GDC TCGA Hiseq & survival Dataset

sample_exps_all <- fread(file.path(toil_path, 'tcga_RSEM_gene_tpm.gz'))
sample_category <- fread(file.path(toil_path, 'TcgaTargetGTEX_phenotype.txt.gz'))


#### Functions ####
# 1. ge_tcga : get gene expression profiles of target gene 
# 2. ge_survival : survival analysis of target_gene, target_cancer pair at given expression threshold
# 3. ge_survival_allcancer : survival analysis of target gene for all cancers, output : list 
# 4. ge_survival_allcancer_genelists : survival analysis of genelist for all cancers, output : list 
# 5. survivalgraph2ppt : survival graphs to pptx (input : fit ,output : pptx file )
# 6. fit2dataframe : fit data to ordered dataframe by log rank p value 


ge_tcga <- function(target_gene, target_cancers = 'all', loading = TRUE) {
  ge_data <- NULL
  
  # get index using gene symbol
  ind_gene <- which(target_gene == annotation$gene)
  
  # use full names of cancer
  if(target_cancers == 'all'){
    target_cancers <- tcga_types$tcga_fullnames
  } else{
    target_cancers <- tcga_types$tcga_fullnames[tcga_types$tcga_symbols %in% target_cancers]
    
  }
  
  if(loading == FALSE){
    # sample phenotype file -  sample detailed_category primary disease or tissue _primary_site  _sample_type _gender _study
    sample_category <- fread('~/Documents/TCGA-survival-ASC/toil-data/TcgaTargetGTEX_phenotype.txt.gz')
    sample_exps_ids <- fread('~/Documents/TCGA-survival-ASC/toil-data/tcga_RSEM_gene_tpm.gz', nrows = 1, sep='\t') %>% colnames()
    sample_exps_ids <- sample_exps_ids[-1]
    sample_exps_raw <- fread('~/Documents/TCGA-survival-ASC/toil-data/tcga_RSEM_gene_tpm.gz', skip = ind_gene ,nrows = 1, sep='\t', header = FALSE)
    sample_exps <- as.numeric(sample_exps_raw[1,-1])
    sample_exps <- data.frame(sample = sample_exps_ids, exp = sample_exps)
    
  }else{
      sample_exps_ids <- colnames(sample_exps_all)[-1]
      sample_exps <- as.numeric(sample_exps_all[ind_gene, -1])
      sample_exps <- data.frame(sample = sample_exps_ids, exp = sample_exps)
    
  }
  
  # add sample info 
  sample_ids <- sample_category %>% dplyr::filter(`_study` == 'TCGA') %>%  dplyr::filter(detailed_category %in% target_cancers)
  sample_exps <- sample_exps %>%  inner_join(sample_ids, by=c('sample'='sample'))
  
  # Add tcga cancer type and symbol 
  sample_exps <- sample_exps %>% 
    left_join(tcga_types, by = c('detailed_category'='tcga_fullnames'))
  
  # Add tumor or normal info ( tumor = 0 , normal = 1 )
  sample_exps <- sample_exps %>% mutate(tumor_sample_or_not = factor(ifelse(grepl('^(\\S*)-(\\S*)-(\\S*)-0(\\S*)$', sample), 'tumor', 'normal'), levels = c('tumor', 'normal')))
  
  # Load to container
  ge <- list()
  ge$data <- sample_exps
  expplot <- ggstripchart(sample_exps, x = "tcga_symbols", y = "exp",
                          color = "tumor_sample_or_not",
                          xlab = "TCGA Cancer Types",
                          ylab = expression(paste(log[2], "(TPM + 0.001)",sep=""))) +
    geom_hline(aes(yintercept = median(exp)), color = alpha("blue", 0.5)) +
    labs(title = paste0(target_gene, " expression")) +
    #scale_y_continuous(expand = c(0, 0)) +
    geom_boxplot(alpha = 0, outlier.color = NA, width = 0.5) +
    stat_summary(fun = median, geom = "point", shape = 18,
                 size = 2.5, color = "red") +
    rotate_x_text(45) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))
  ge$expplot <- expplot
  print(expplot)
  return(ge)
}




ge_survival <- function(target_gene,target_cancer, thrsld = 0.2, loading=TRUE, regression = FALSE, return_plot = FALSE,legend_dodge = FALSE, risk_table = FALSE) {
  
  # get index using gene symbol
  ind_gene <- which(target_gene == annotation$gene)
  
  # use full names of cancer
  target_cancer <- tcga_types$tcga_fullnames[tcga_types$tcga_symbols == target_cancer]
  
  # get tcga sample ids of target cancer 
  
  sample_ids <- sample_category %>% dplyr::filter(`_study` == 'TCGA') %>%  dplyr::filter(detailed_category == target_cancer)
  
  if(loading == FALSE){
    sample_exps_ids <- fread('~/Documents/TCGA-survival-ASC/toil-data/tcga_RSEM_gene_tpm.gz', nrows = 1, sep='\t') %>% colnames()
    sample_exps_ids <- sample_exps_ids[-1]
    sample_exps_raw <- fread('~/Documents/TCGA-survival-ASC/toil-data/tcga_RSEM_gene_tpm.gz', skip = ind_gene ,nrows = 1, sep='\t', header = FALSE)
    sample_exps <- as.numeric(sample_exps_raw[1,-1])
    sample_exps <- data.frame(sample = sample_exps_ids, exp = sample_exps)
  }
  else{
    
    sample_exps_ids <- colnames(sample_exps_all)[-1]
    sample_exps <- as.numeric(sample_exps_all[ind_gene, -1])
    sample_exps <- data.frame(sample = sample_exps_ids, exp = sample_exps)
  }
  
  # add sample info 
  sample_exps <- sample_exps %>%  inner_join(sample_ids, by=c('sample'='sample'))
  
  
  # Add tcga cancer type and symbol 
  sample_exps <- sample_exps %>% 
    left_join(tcga_types, by = c('detailed_category'='tcga_fullnames'))
  
  # Add tumor or normal info ( tumor = 0 , normal = 1 )
  sample_exps <- sample_exps %>% mutate(tumor_sample_or_not = as.factor(ifelse(grepl('^(\\S*)-(\\S*)-(\\S*)-0(\\S*)$', sample), 0, 1)))
  # filter tumor sample only
  sample_exps <- sample_exps %>% dplyr::filter(tumor_sample_or_not == 0)
  exp <- sample_exps$exp
  
  cutoff <- quantile(exp, probs = c(thrsld, 1 - thrsld), na.rm = TRUE)
  
  high_exp <- sample_exps$exp[sample_exps$exp >= cutoff[2]]
  low_exp <- sample_exps$exp[sample_exps$exp <= cutoff[1]]
  
  survival <- fread("~/Documents/TCGA-survival-ASC/toil-data/TCGA_survival_data")
  sample_exps <- sample_exps %>% left_join(survival, by = 'sample')
  sample_exps <- sample_exps %>% dplyr::mutate(group = ifelse(exp >= cutoff[2], paste0("High ", target_gene), ifelse(exp <= cutoff[1], paste0("Low ", target_gene), NA)))
  
  
  high_n <- nrow(sample_exps[which(sample_exps$'group' == paste0("High ", target_gene)), ])
  low_n <- nrow(sample_exps[which(sample_exps$'group' == paste0("Low ", target_gene)), ])
  sample_exps$group <- as.factor(sample_exps$group)
  sample_exps$OS.time <- sample_exps$OS.time / 30
  
  fit <- NULL
  try({
    fit <- survfit(Surv(OS.time, OS) ~ group, data = sample_exps[!is.na(sample_exps$group),])
    fit$pval <- surv_pvalue(fit, data = sample_exps[!is.na(sample_exps$group),])$pval
  }, silent = TRUE)
  
  
  if(!is.null(fit)){
    fit$target_cancer <- target_cancer
    fit$high_exp <- high_exp
    fit$low_exp <- low_exp
    fit$exp_fc <- 2^(mean(high_exp) - mean(low_exp))
    fit$exp_pval <- t.test(high_exp, low_exp)$p.val
    fit$high_n <- high_n
    fit$low_n <- low_n
    fit$gene <- target_gene
    fit$sample_exps <- sample_exps
    
  }
  else{
    fit$target_cancer <- target_cancer
    fit$pval <- NA
    fit$high_exp <- NA
    fit$low_exp <- NA
    fit$exp_fc <- NA
    fit$exp_pval <- NA
    fit$high_n <- NA
    fit$low_n <- NA
    fit$gene <- target_gene
    fit$sample_exps <- sample_exps
    
  }
  
  
  if(regression == TRUE){
    fit$day <- NA
    try(
      {   reg <- survreg(Surv(OS.time, OS) ~ group, data = sample_exps[!is.na(sample_exps$group),])
      pct <- c(0.5, 0.7)
      ptime <- predict(reg, newdata=data.frame(group=names(table(sample_exps[!is.na(sample_exps$group),]$group))), type='quantile', p=pct)
      day <- round(as.numeric(ptime) * 30)
      names(day) <- c("high 50% surv", "low 50% surv", "high 70% surv", "low 70% surv")
      fit$day <- day
      }, silent = TRUE)
    
  }
  if(return_plot == TRUE) {
    try({
      
      fit$survplot <- ggsurvplot(fit, data = sample_exps[!is.na(sample_exps$group),] ,xlab = "Follow-Up Time (Months)", ylab = "Percent Survival",
                                 palette = c("red", "blue"), 
                                 risk.table = FALSE, conf.int = TRUE, )
      plab <- paste0("Log-rank P ",
                     ifelse(fit$pval < 0.001, "< 0.001",
                            paste0("= ", round(fit$pval, 4))), '\n', 'high : ',high_n, ' low : ', low_n)
      suppressMessages(
        fit$survplot$plot <- fit$survplot$plot +
          scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
          scale_x_continuous(expand = c(0, 0)) +
          labs(title = paste0(target_gene, " expression ~ overall survival in ", target_cancer),
               fill = plab, color = plab) +
          guides(color = guide_legend(title.position = ifelse(legend_dodge, "top", "bottom"))) +
          theme(aspect.ratio = 1,
                plot.title = element_text(hjust = 0.5, size = 12), 
                legend.title.align = 1,
                legend.title = element_text(size = 12, hjust = 1),
                legend.text = element_text(size = 12), 
                legend.position = c(0.98, ifelse(legend_dodge, 0.12, 0.93)),
                legend.justification = "right",
                legend.background = element_blank())
      )
      fit$survplot$plot <- ggpar(fit$survplot$plot, 
                                 font.x = c(12),
                                 font.y = c(12),
                                 font.tickslab = c(12))
    }, silent = TRUE)
    
  }
  
  return(fit)
}




ge_survival_allcancer <- function(target_gene, thrsld, return_plot = TRUE, regression = TRUE){
  # need tcga_types vector
  survp <- list()
  for(i in 1:length(tcga_types$tcga_symbols)){
    
    tcga_type <- tcga_types$tcga_symbols[i]
    try(
      
      {survp[[i]] <- ge_survival(target_gene = target_gene, target_cancer = tcga_type, thrsld = thrsld,return_plot = return_plot, regression = regression)}, silent = TRUE
      
    )
    
  }
  
  return(survp)
  
}

ge_survival_allcancer_genelists <- function(target_genes, thrsld = 0.2, return_plot = TRUE, regression = TRUE){
  survp <- list()
  for(i in 1:length(target_genes)){
    target_gene <- target_genes[i]
    # cat(i, ' / ', length(target_genes), '\n')
    plts <- ge_survival_allcancer(target_gene, thrsld = thrsld, return_plot = return_plot, regression = regression)
    survp <- append(survp, plts)
  }
  
  return(survp)
  
}

# get graph 


survivalgraph2ppt <- function(fit, path = '~/Documents/ER-UPR/210923-ER-UPR-all-survival_plot-thrsld-02.pptx'){
  # get graphs for fit objects 
  l_surv <- list()
  list_idxs <- as.numeric(fit$fitdf$list_idx)
  for(i in 1:length(fit)){
    idx <- list_idxs[i]
    try(l_surv[[i]] <- fit[[idx]]$survplot$plot)
    
  }
  
  graph2ppt(l_surv, path = path, width = 7, height = 7)
  
  
}

fit2dataframe <- function(fit){
  
  fitdf <- lapply(fit, x <- function(x){list(x$gene,x$target_cancer,x$pval, 
                                             median(x$high_exp, na.rm = TRUE), median(x$low_exp, na.rm = TRUE),
                                             x$high_n, x$low_n)}) %>%
    rbindlist() %>% rownames_to_column('list_idx')%>% 
    `colnames<-`(c('list_idx','gene', 'cancer', 'pval', 'high_exp', 'low_exp', 'high_n', 'low_n')) %>% 
    arrange(pval)
  
  fit$fitdf <- fitdf 
  
  return(fit)
}



ge_survival_parallel <- function(target_cancer, target_genes, thrsld = 0.2, return_plot = FALSE, regression = TRUE, path = '~/Documents/TCGA_survival_analysis/ALL/results/thrs50'){
  # survival analysis for given gene lists for one cancer type. 
  # further parallel for cancer types 
  survp <- data.frame(gene = character(),
                      target_cancer = character(),
                      pval = numeric(),
                      high_exp = numeric(),
                      low_exp = numeric(),
                      high_n = numeric(),
                      low_n = numeric()
    
  )
  name_rds <- paste0('_thrs_', thrsld*100, '.RDS')
  name_csv <- paste0('_thrs_', thrsld*100, '.csv')
  n<- 1
  if(file.exists(file.path(path, paste0('ALL_survival_', target_cancer, name_rds)) )){
    survp <- readRDS(file.path(path, paste0('ALL_survival_', target_cancer, name_rds)))
    print(paste0('load survp ',nrow(survp)))

    n <- nrow(survp) + 1
    
  }

  
  
  for(i in n:length(target_genes)){
    target_gene <- target_genes[i]
    x <- NULL
    try({
        x <- ge_survival(target_gene = target_gene, target_cancer = target_cancer, thrsld = thrsld,return_plot = return_plot, regression = regression)

      }, silent = TRUE)
    
    if(is.null(x)){x <- NA}
    if(!is.na(x[1])){
      survp <- rbindlist(list(survp, list(x$gene,x$target_cancer,x$pval, 
                                            median(x$high_exp, na.rm = TRUE), median(x$low_exp, na.rm = TRUE),
                                            x$high_n, x$low_n)), fill=TRUE)
      
    }else{
      survp <- rbindlist(list(survp, as.list(c(target_gene,tcga_types$tcga_fullnames[tcga_types$tcga_symbols == target_cancer],
                                               rep(NA,5)))), fill=TRUE)
      
    }
    
    if(i %% 1000 == 12|i == length(target_genes)){
      
      survp <- survp %>% arrange(pval)
      saveRDS(survp, file = file.path(path, paste0('ALL_survival_', target_cancer, name_rds)))
    }
    
  }
  survp <- survp %>% arrange(pval)
  survp$pval <- as.numeric(survp$pval)
  survp$high_exp <- as.numeric(survp$high_exp)
  survp$low_exp <- as.numeric(survp$high_low)
  write.csv(survp, file.path(path,paste0('ALL_survival_', target_cancer, name_csv)),row.names = FALSE)
  saveRDS(survp, file = file.path(path, paste0('ALL_survival_', target_cancer, name_rds)))
  survp <- NULL
  return(survp)
  
}



gene_tcga_eval <- function(cancer, thrsld = 0.2){
  if(!dir.exists(paste0('results-thrsld-',thrsld))){dir.create(paste0('results-thrsld-',thrsld))}
  screen <- data.frame(gene = character(),
                       pval = numeric(),
                       exp_pval = numeric(),
                       exp_fc = numeric(),
                       high_n = numeric(),
                       low_n = numeric()
                       
  )
  name <- paste0('result-',cancer)
  n <- 1
  if(file.exists(file.path(paste0('results-thrsld-',thrsld), paste0('result-',cancer,'.RData') ))){
    load(file.path(paste0('results-thrsld-',thrsld), paste0('result-',cancer,'.RData')))
    print(paste0('load screen ',nrow(screen)))
    n <- nrow(screen) + 1
    
  }
  
  samples <- fread(paste0("./data/TCGA-", cancer, ".htseq_fpkm.tsv.gz"),
                   stringsAsFactors = FALSE, sep = "\t", skip = 0, nrows = 1, header = FALSE)
  samples <- as.character(samples[, -1])
  
  for(i in n:length(genes)){
    
    gene <- genes[i]
    fit <- ge_survival_loop(idx = i ,target_gene = gene, target_cancer = cancer, samples = samples,
                            thrsld = thrsld, return_plot = FALSE, regression = FALSE)
    if(is.null(fit)){fit <- NA}
    if(!is.na(fit[1])){
      screen <- rbindlist(list(screen, as.list(c(gene, fit$pval,fit$exp_pval,fit$exp_fc,fit$high_n,fit$low_n))), fill=TRUE)
      
    }else{
      screen <- rbindlist(list(screen, as.list(c(gene,rep(NA,5)))), fill=TRUE)
      
    }
    if(i %% 1000 == 10|i == length(genes)){
      save(screen, file = file.path(paste0('results-thrsld-',thrsld), paste0('result-',cancer,'.RData')))
    }
  }
  
  write.csv(screen, file = file.path(paste0('results-thrsld-',thrsld), paste0('result-',cancer,'.csv')),row.names = FALSE)
  save(screen, file = file.path(paste0('results-thrsld-',thrsld), paste0('result-',cancer,'.RData')))
  screen <- NULL
  return(screen)
  
}




cat(" TCGA-survival functions loaded !!! 

      # 1. ge_tcga : get gene expression profiles of target gene 
      # 2. ge_survival : survival analysis of target_gene, target_cancer pair at given expression threshold
      # 3. ge_survival_allcancer : survival analysis of target gene for all cancers, output : list 
      # 4. ge_survival_allcancer_genelists : survival analysis of genelist for all cancers, output : list 
      # 5. survivalgraph2ppt : survival graphs to pptx (input : fit ,output : pptx file )
      # 6. fit2dataframe : fit data to ordered dataframe by log rank p value 
      
      
      ")


