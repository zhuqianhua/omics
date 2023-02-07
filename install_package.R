#!/ifswh1/BC_PUB/biosoft/pipeline/Package/R-3.2.1/bin/Rscript

#--- install packages
path = '/ifswh1/BC_PUB/biosoft/pipeline/RNA/Omics_approach/Omics_approach_2017a/R_package'

.libPaths(path)
rpos = 'http://cloud.r-project.org'
if (!file.exists(path))
  dir.create(path, recursive = TRUE)

pCran = c("getopt", "hash", "RColorBrewer", "ape", "plyr", "pheatmap", "plotrix", "gplots")
pBioc = c("WGCNA", "KEGG.db", "GOstats", "GSEABase", "AnnotationDbi", "limma",
          "org.Hs.eg.db", "org.Sc.sgd.db", "org.Ag.eg.db", "org.Pt.eg.db",
          "org.Rn.eg.db", "org.Ss.eg.db", "org.At.tair.db", "org.Bt.eg.db",
          "org.Ce.eg.db", "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db",
          "org.Gg.eg.db", "org.Mm.eg.db", "org.Mmu.eg.db", "org.Xl.eg.db")

#--- install packages from CRAN
invisible(lapply(pCran, 
                 function(x) {
                   if (x %in% installed.packages()[, 1] == FALSE)
                     install.packages(x, lib = path, repos = rpos)
                 }))
#--- install packages from Bioconductor
source("http://bioconductor.org/biocLite.R")
invisible(lapply(pBioc, 
                 function(x) { 
                   if (x %in% installed.packages()[, 1] == FALSE)
                     biocLite(x, lib = path)
                 }))
