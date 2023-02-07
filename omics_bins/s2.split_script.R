split_script <- function(indir, suffix, prefix, piece, rscr, rlib, prog) {
  dat = Sys.glob(file.path(indir, suffix))
  if (length(dat) < piece)
    dat[(length(dat) + 1):piece] = rep("echo do nothing", (piece - length(dat)))
  fac = vector()
  #--- Non-uniform distribution, modified on JAN 22th, 2018
  #num = ifelse(length(dat) %% piece == 0, length(dat) / piece, floor(length(dat) / piece)+1)
  #num = floor(length(dat) / piece)
  #invisible(lapply(1:length(dat), function(x) {
  #  tfc = ifelse(x %% num == 0, x / num, floor(x / num)+1)
  #  if (tfc > piece)
  #    tfc = piece
  #  fac <<- c(fac, tfc)}))
  mak = 1
  invisible(lapply(1:length(dat), function(x) {
    if (mak > piece)
      mak = 1
    fac <<- c(fac, mak)
    mak <<- mak + 1}))
  spt = split(dat, factor(fac))
  invisible(lapply(1:length(spt), function(x) {
    key = names(spt[x])
    val = spt[[x]]
    cmd = NULL
    lapply(1:length(val), function(y) {
      tsk = paste(rscr, prog, val[y], rlib, " && \\", sep = " ")
      if (rlib == "NA" & val[y] != "echo do nothing") {
        net = unlist(read.table(val[y]))
        tsk = paste(rscr, prog, "--nodes", net[2], "--edges", net[3], 
                    "--config", net[4], "--output_prefix", net[1], " && \\", sep = " ")
      }
      cmd <<- c(cmd, ifelse(grepl("echo do nothing", val[y]) == TRUE, 
                            paste(val[y], " && \\", sep = ""), tsk))
    })
    if (!file.exists(dirname(prefix)))
      dir.create(dirname(prefix), recursive = TRUE)
    write.table(c("#!/bin/bash",
                  "echo ==========start at: `date` ========== && \\", cmd,
                  "echo ==========end at: `date` ========== && \\",
                  "echo Still_waters_run_deep 1>&2 && \\",
                  paste("echo Still_waters_run_deep > ", paste(prefix, "_", key, ".sh", sep = ""),
                    ".sign", sep = "")), paste(prefix, "_", key, ".sh", sep = ""), 
                  quote = FALSE, sep = "\n", col.names = FALSE, row.names = FALSE)
  }))
}

opts <- commandArgs(trailingOnly = T)
if (length(opts) != 7) {
  cat("./Rscript s2.split_script.R <indir> <suffix> <prefix> <piece> <rscript> <rlib> <prog>\n")
  quit(save = "no", status = 1, runLast = TRUE);
}
split_script(opts[1], opts[2], opts[3], as.numeric(opts[4]), opts[5], opts[6], opts[7])
