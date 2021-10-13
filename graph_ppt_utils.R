#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 210827 HSY
# Gene to Graph and PPTX utils
# 
# gene_to_graph
# genelist_to_graph
# graph2ppt
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(ggpubr)


gene_to_graph <- function(gene, probes, legend = 'none'){
  # need probes dataset

  probes_ids <- probes$`Probe Set ID`[probes$`Gene symbol` == gene]
  myplots <- list()
  
  mypal = pal_npg()(10)
  
  for(i in 1:length(probes_ids)){
    probe_id <- probes_ids[i]
    genename <- probes$Genename[probes$`Probe Set ID` == probe_id]
    gene_exp <- nci_exp[nci_exp$ID_REF == probe_id, ][-1] %>% t(.)
    gene_exp <- as.vector(gene_exp)

    graph_val <- pheno
    graph_val$broad_group <- factor(graph_val$broad_group, levels = unique(pheno$broad_group))
    graph_val$exp <- gene_exp
    
    n <- length(unique(graph_val$broad_group))
    set.seed(005)
    pal <- sample(colorRampPalette(mypal)(n))
    plt <- ggbarplot(graph_val, x = 'cell', y = 'exp', add = c('mean_se','jitter'), 
                     color = 'broad_group',fill = 'broad_group', palette = pal,
                     title = paste(genename, paste(gene, probe_id, sep = ' : '), sep= '\n'), 
                     ylab = 'expression' , xlab = FALSE, label = FALSE, legend = legend, lab.size = 3) +
      rotate_x_text(angle=65)
    
    myplots[[i]] <- plt
  }
  
  return(myplots)
}



genelist_to_graph <- function(geneset, probes){
  l <- list()
  for(i in 1:length(geneset)){
    plts <- NULL
    genei <- geneset[i]
    try(
      plts <- gene_to_graph(genei, probes), silent = T)
    
    l <- append(l, plts)
    cat(genei, i, '/', length(geneset), '\n')
  }
  return(l)
}



graph2ppt <- function(l, path, left = 0, top = 1, width = 14, height = 4.95){
  dmls <- list()
  for(i in 1:length(l)){
    try({
    img <- l[[i]]
    dml <- rvg::dml(ggobj = img)
    dmls[[i]] <- dml
    })
  }
  
  if (!file.exists(path)) {
    out <- officer::read_pptx(path = '~/Documents/ER-UPR/NCI60_template.pptx')
  }
  # if file exist, append slides to exisiting file ----
  else {
    out <- officer::read_pptx(path)
  }
  
  free_loc <- officer::ph_location(
    left = left, top = top, 
    width = width, height = height)
  
  for(i in 1:length(dmls)){
    
    try({
      dml <- dmls[[i]]
      out <- out %>% 
      officer::add_slide() %>% 
      officer::ph_with(value = dml, location = free_loc) })
  }
  print(out,target = path)
}









