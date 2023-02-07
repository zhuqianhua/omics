#--- generate nodes and edge files for network

prepare_network <- function(path)
{
  library(RColorBrewer)
  if (!dir.exists(paste(path, "/Enrich_overlap", sep = "")))
    return(0)
  if (file.exists(paste(path, "/Network", sep = "")))
    unlink(paste(path, "/Network", sep = ""), recursive = TRUE)
  invisible(lapply(Sys.glob(file.path(paste(path, "/Enrich_overlap", sep = ""), '*')), node_edge))
}

node_edge <- function(indir)
{
  if (file.exists(paste(indir, "/enrich.txt", sep = "")) == FALSE)
    return(0)
  otdir = gsub("Enrich_overlap", "Network", indir)
  if (!dir.exists(otdir)) 
	dir.create(otdir, recursive = TRUE)
  enrTb = as.vector(unlist(read.table(paste(indir, "/enrich.txt", sep = ""))))
  if (length(enrTb) < 8)
    return(0)
  intTb = read.table(enrTb[5], header = TRUE, row.names = 1)
  nshp = c('sphere', 'star', 'crectangle', 'triangle')
  ncol = brewer.pal(9, "OrRd")
  onod = paste(otdir, "/node.txt", sep = "")
  oedg = paste(otdir, "/edge.txt", sep = "")
  write.table(t(c('name', 'shape', 'label.cex', "color", "frame.color")), onod, 
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(t(c('from', 'to')), oedg, sep = "\t", quote = FALSE, 
              col.names = FALSE, row.names = FALSE)
  lapply(7:length(enrTb), 
         function(x) {
           mshp = ifelse((x - 6) %% 4 == 0, 4, (x - 6) %% 4)
           expTb = read.table(enrTb[x], header = TRUE, row.names = 1)
           link = paste(dirname(dirname(indir)), "/Trans_", 
                         gsub(".txt$", "", basename(enrTb[x])), "/link.xls", sep = "")
           if (file.exists(link) == TRUE) {
             lnkTb = read.table(link, sep = "\t")
             write.table(lnkTb[which(lnkTb$V2 %in% rownames(intTb)),], oedg, sep = "\t",
                         quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
           }
           expOd = order(rowMeans(expTb))
           lapply(1:nrow(expTb), 
                  function(y) {
                    mmen = round(expOd[y] * 9 / nrow(expTb))
                    mcol = ifelse(mmen == 0, 1, mmen)
                    write.table(t(c(rownames(expTb[y,]), nshp[mshp], 0.8, ncol[mcol], ncol[mcol])),
                                onod, sep = "\t", quote = FALSE, col.names = FALSE, 
                                row.names = FALSE, append = TRUE)
                  })
           
         })
  nodCk = read.csv(onod, header = TRUE, sep = "\t")
  write.table(unique(nodCk), onod, sep = "\t", quote = FALSE, row.names = FALSE)
  edgCk = read.table(oedg, header = TRUE, sep = "\t")
  if (nrow(edgCk) < 2) {
    unlink(otdir, recursive = TRUE)
    return(0)
  }
  write.table(edgCk[which(edgCk$from %in% nodCk$name & edgCk$to %in% nodCk$name),],
              oedg, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(t(c("[options]", "degree = 6")), paste(otdir, "/config.txt", sep = ""),
              sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(t(c(paste(gsub("00.Intermediate/Enrich_overlap", "04.Network", indir), "/", basename(indir),
                        sep = ""), onod, oedg, paste(otdir, "/config.txt", sep = ""))), 
              paste(otdir, "/network.txt", sep = ""), sep = "\n", row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
}

opts <- commandArgs(trailingOnly = T)
if (length(opts) != 2) {
  cat("./Rscript s3.omics_network.R <indir> <rlib>\n")
  quit(save = "no", status = 1, runLast = TRUE)
}
.libPaths(opts[2])
prepare_network(opts[1])
