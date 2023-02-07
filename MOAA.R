#!/ifs4/BC_PUB/biosoft/pipeline/Package/R-3.3.1/bin/Rscript

#--- Multiple Omics Analysis Approache
rbin = '/ifs4/BC_PUB/biosoft/pipeline/RNA/Omics_approache/Omics_approache_2017a'
rscr = '/ifs4/BC_PUB/biosoft/pipeline/Package/R-3.3.1/bin/Rscript'
qsub = '/ifs4/BC_PUB/biosoft/pipeline/Package/pymonitor-1.1/monitor.cpu'
xbio = '/ifs4/BC_PUB/biosoft/bin/upload.sh.x'
pyth = 'python'
netw = '/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_ceRNA/ceRNA_2017a/ceRNA_Network/network.py'

rlib = paste(rbin, "/R_package", sep = "")
.libPaths(rlib)

OMAA_usage <- function() {
  cat(paste("Program: Multiple Omics Analysis Approache",
            "Version: 1.0",
            "Contact: Qianhua ZHU <zhuqianhua@bgi.com>  Chichuan LIU <liuchichuan@bgi.com>",
            "Usage  : ./MOAA.R [options]\n",
            "Options: -c, --config  *<s>  configuer file",
            "         -o, --outdir   <s>  outdir, default \"./\"",
      paste("         -r, --rscript  <s>  path of Rscript, default ", rscr, sep = ""),
            "         -h, --help     <b>  print this information\n\n",
            sep = "\n"))
  quit(status = 1)
}

monitor_shell <- function(shell, command) {
  write.table(c("#!/bin/bash", 
                "echo ==========start at: `date` ========== && \\", command,
                "echo ==========end at: `date` ========== && \\",
                "echo Still_waters_run_deep 1>&2 && \\",
          paste("echo Still_waters_run_deep > ", shell, ".sign", sep = "")), shell,
                quote = FALSE, sep = "\n", col.names = FALSE, row.names = FALSE)
}

OMAA_main <- function(opt) {
  library(hash)
  datPara = hash('queue' = 'bc.q', 'project' = 'Project', 'subpro' = '', 'mail' = 'user@genomics.cn',
                 'cdtsuser' = 'user', 'cdtspasswd' = 'passwd', 'cdtsloc' = '',
                 'slicenum' = 20, 'interenri' = 0)
  datFile = file(opt$config, "r")
  datLine = readLines(datFile, n = 1, warn = FALSE)
  while (length(datLine) != 0) {
    lineVec = unlist(strsplit(datLine, '(\\s+)?=(\\s+)?'))
    if (length(lineVec) == 2 & length(grep('#', as.character(lineVec[1]))) == 0)
      if (has.key(tolower(lineVec[1]), datPara) == TRUE)
        datPara[tolower(lineVec[1])] = lineVec[2]
    datLine = readLines(datFile, n = 1, warn = FALSE)
  }
  close(datFile)
  shDir = paste(opt$outdir, "/shell", sep = "")
  reDir = paste(opt$outdir, "/Omics_analysis", sep = "")
  if (!file.exists(shDir))
    dir.create(shDir, recursive = TRUE)
  monitor_shell(paste(shDir, "/omics_overlap.sh", sep = ""),
                paste(rscr, paste(rbin, "/omics_bins/s1.omics_overlap.R", sep = ""), 
                opt$config, rlib, reDir, sep = " "))
  dep = paste(shDir, "/omics_overlap.sh:1.5G", sep = "")
  if (datPara[['interenri']] != 0) {
    monitor_shell(paste(shDir, "/split_enrichment.sh", sep = ""),
                  paste(rscr, paste(rbin, "/omics_bins/s2.split_script.R", sep = ""),
                    paste(opt$outdir, "/Omics_analysis/00.Intermediate", sep = ""), "*/*/enrich.txt",
                    paste(shDir, "/enrich/enrich", sep = ""), datPara[['slicenum']], rscr, rlib, 
                    paste(rbin, "/omics_bins/s3.omics_enrichment.R", sep = ""), sep = " "))
    if (!file.exists(paste(shDir, "/enrich", sep = "")))
      dir.create(paste(shDir, "/enrich", sep = ""), recursive = TRUE)
    monitor_shell(paste(shDir, "/omics_venn.sh", sep = ""),
                  paste(rscr, paste(rbin, "/omics_bins/s4.omics_enrich_venn.R", sep = ""),
                    paste(opt$outdir, "/Omics_analysis/00.Intermediate", sep = ""), rlib, sep = " "))
    dep = c(dep, paste(paste(shDir, "/omics_overlap.sh:1.5G", sep = ""),
                       paste(shDir, "/split_enrichment.sh:0.1G", sep = ""), sep = "\t"))
    invisible(lapply(1:datPara[['slicenum']],
                function(x) {
                  ensh = paste(shDir, "/enrich/enrich_", x, ".sh:1G", sep = "")
                  dep <<- c(dep, paste(paste(shDir, "/split_enrichment.sh:0.1G", sep = ""), ensh, sep = "\t"))
                  dep <<- c(dep, paste(ensh, paste(shDir, "/omics_venn.sh:0.5G", sep = ""), sep = "\t"))
                }))
  }
  monitor_shell(paste(shDir, "/omics_network.sh", sep = ""),
                paste(rscr, paste(rbin, "/omics_bins/s3.omics_network.R", sep = ""),
                  paste(opt$outdir, "/Omics_analysis/00.Intermediate", sep = ""), rlib, sep = " "))
  dep = c(dep, paste(paste(shDir, "/omics_overlap.sh:1.5G", sep = ""), 
                       paste(shDir, "/omics_network.sh:0.1G", sep = ""), sep = "\t")) 
  monitor_shell(paste(shDir, "/split_network.sh", sep = ""),
                paste(rscr, paste(rbin, "/omics_bins/s2.split_script.R", sep = ""),
                  paste(opt$outdir, "/Omics_analysis/00.Intermediate", sep = ""), "*/*/network.txt",
                  paste(shDir, "/network/network", sep = ""), datPara[['slicenum']], pyth, 'NA', netw, sep = " "))
  if (!file.exists(paste(shDir, "/network", sep = "")))
    dir.create(paste(shDir, "/network", sep = ""), recursive = TRUE)
  dep = c(dep, paste(paste(shDir, "/omics_network.sh:0.1G", sep = ""),
                       paste(shDir, "/split_network.sh", sep = ""), sep = "\t"))
  invisible(lapply(1:datPara[['slicenum']], 
              function(x) {dep <<- c(dep, paste(paste(shDir, "/split_network.sh:0.1G", sep = ""),
                                     paste(shDir, "/network/network_", x, ".sh:0.5G", sep = ""), sep = "\t"))
              }))
  write.table(dep, paste(shDir, "/dependence.txt", sep = ""), sep = "\n", 
              quote = FALSE, col.names = FALSE, row.names = FALSE)
  # s1.qsub.sh
  write.table(paste(qsub, "taskmonitor -q ", datPara[['queue']], "-P", datPara[['subpro']],
                    "-p", datPara[['project']], "-i", paste(shDir, "/dependence.txt", sep = ""),
                    sep = " "), paste(opt$outdir, "/s1.qsub.sh", sep = ""),
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
  # s2.report.sh
  write.table(paste(paste(rscr, paste(rbin, "/omics_bins/s5.omics_report.R", sep = ""), 
                          paste(opt$outdir, "/Omics_analysis", sep = ""), datPara[['project']], 
                          paste(rbin, "/omics_bins/omics_report", sep = ""), sep = " "), " && \\\ncd ",
                    opt$outdir, " && \\\ntar zcf Omics_report.tar.gz arf resource Omics_analysis/0[^0].* && \\\n",
                    paste(xbio, "-usr", datPara[['cdtsuser']], "-pwd",  datPara[['cdtspasswd']], 
                          "-area",  datPara[['cdtsloc']], paste("-f Omics_report.tar.gz -p /", 
                          datPara[['project']], "/", datPara[['subpro']], sep = ""), "-mail", 
                          datPara[['mail']], "-id", datPara[['subpro']], "-n", "'Omics Report'", sep = " "), sep = ""), 
              paste(opt$outdir, "/s2.report.sh", sep = ""),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

library(getopt)
opts = matrix(c(
  'config', 'c', 1, 'character',
  'outdir', 'o', 2, 'character',
  'rscript', 'r', 2, 'character',
  'help', 'h', 0, 'logical'
  ), byrow = TRUE, ncol = 4);
opt = getopt(opts)


if (is.null(opt$config) || !is.null(opt$help))
  OMAA_usage()
if (is.null(opt$outdir))
  opt$outdir = './'
opt$config = normalizePath(opt$config)
if (!file.exists(opt$outdir))
  dir.create(opt$outdir, recursive = TRUE)
opt$outdir = normalizePath(opt$outdir)
OMAA_main(opt)
