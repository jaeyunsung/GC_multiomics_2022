library("BiocManager")

tutorial_pkgs <- c('devtools','doBy','dplyr','phyloseq','vegan','DESeq2','pheatmap','scales','RColorBrewer','pcaExplorer','AnnotationDbi','org.Hs.eg.db','gage','gageData','PMA','data.table','patchwork','xCell','grid','gridExtra','reshape','ggcorrplot')

for (pkg in tutorial_pkgs) {
  if( !is.element(pkg, .packages(all.available = TRUE)) ) {
    message(sprintf("installing %s ...",pkg))
    BiocManager::install(pkgs=pkg)
    message("Done.")
  }
  message(sprintf("loading %s ...",pkg))
  library(pkg,character.only = TRUE)
}