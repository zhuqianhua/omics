omics_report <- function(indir, proid, rbin) {
  arf = paste(indir, "/../arf", sep = "")
  if (!file.exists(arf))
    dir.create(arf, recursive = TRUE)
  file.copy(from = paste(rbin, "/resource", sep = ""), to = paste(indir, "/..", sep = ""), 
            copy.mode = TRUE, recursive = TRUE)
  rep_overlap(proid, indir)  
  arf_output(paste(arf, "/results.arf", sep = ""), paste(rbin, "/r1.overlap.arf", sep = ""),
             paste(indir, "/overlap.xls", sep = ""), FALSE)
  file.copy(paste(rbin, "/m1.overlap.arf", sep = ""), paste(arf, "/methods.arf", sep = ""), overwrite = TRUE)
  file.copy(paste(rbin, "/h1.overlap.arf", sep = ""), paste(arf, "/help.arf", sep = ""), overwrite = TRUE)
  file.copy(paste(rbin, "/faq.arf", sep = ""), paste(arf, "/FAQs.arf", sep = ""), overwrite = TRUE)
  if (file.exists(paste(indir, "/03.Enrichment_pattern", sep = ""))) {
    rep_enrich(indir) 
    arf_output(paste(arf, "/results.arf", sep = ""), paste(rbin, "/r2.enrich.arf", sep = ""),
               paste(indir, "/enrich.xls", sep = ""), TRUE)
    file.append(paste(arf, "/methods.arf", sep = ""), paste(rbin, "/m2.enrich.arf", sep = ""))
    file.append(paste(arf, "/help.arf", sep = ""), paste(rbin, "/h2.enrich.arf", sep = ""))
  }
  rep_network(indir)
  arf_output(paste(arf, "/results.arf", sep = ""), paste(rbin, "/r3.network.arf", sep = ""),
             paste(indir, "/network.xls", sep = ""), TRUE)
  file.append(paste(arf, "/methods.arf", sep = ""), paste(rbin, "/m3.network.arf", sep = ""))
  file.append(paste(arf, "/help.arf", sep = ""), paste(rbin, "/h3.network.arf", sep = ""))
}

arf_output <- function(out, mod, val, app) {
  dmod = readLines(mod)
  dval = read.table(val, sep = "\t")
  invisible(lapply(1:nrow(dval), function(x) {
                                   dmod <<- gsub(dval[x,]$V1, dval[x,]$V2, dmod)}))
  write.table(dmod, out, quote = FALSE, row.names = FALSE, 
              append = app, col.names = FALSE)
  file.remove(val)
}

conf_output <- function(outfile, flag, value, app) {
  write.table(t(c(flag, value)), outfile, sep = "\t", append = app, 
              col.names = FALSE, row.names = FALSE, quote = FALSE)
}

rep_fig <- function(path, pattern, suffix) {
  name = vector()
  invisible(lapply(Sys.glob(file.path(path, pattern)), 
                   function(x) {
                     key = gsub(suffix, "", basename(x))
                     fig = sub(".*[^Omics_analysis]+(Omics_analysis*)", "\\1", x)
                     name <<- c(paste("file = <url = ", fig, "; label = \"", key, "\">", sep = ""), name)}))
  return(name)
}

rep_tab <- function(path, pattern, mak) {
  name = Sys.glob(file.path(path, pattern))
  lmak = 10
  if (length(name) < lmak)
    lmak = length(name)
  count = vector()
  invisible(lapply(name[1:lmak],
                   function(x) {
                     key = basename(x)
                     tab = sub(".*[^Omics_analysis]+(Omics_analysis*)", "\\1", x)
                     if (mak == 1) {
                       count <<- c(paste("@table\ntitle = \"", key, "\"\nfile = <url = ", tab, ">\nview = 0\n", sep = ""), count)
                     } else {
                       htmlf = gsub(".html", "_files", x)
                       vec = vector()
                       lapply(c(Sys.glob(file.path(htmlf, "*/*")), Sys.glob(file.path(htmlf, "*/*/*/*"))), 
                              function(y) {
                                nam = sub(".*[^Omics_analysis]+(Omics_analysis*)", "\\1", y)
                                vec <<- c(paste("\\\\href{", nam, "}{} ", sep = ""), vec) })
                       count <<- c(paste("\\\\href{", tab, "}{", key, "} ", 
                                         paste(vec, collapse = " "), sep = ""), count)
                     }}))
  return(count)
}

rep_network <- function(indir) {
  #--- remove empty directory
  invisible(lapply(Sys.glob(file.path(paste(indir, "/04.Network", sep = ""),  "*")), 
                   function(x) { #print(paste(x, ": ", file.info(x)$size, sep = ""))
                     if (file.info(x)$size == 0)
                       file.remove(x, recursive = TRUE) }))
  conf = paste(indir, "/network.xls", sep = "")
  name = rep_fig(paste(indir, "/04.Network/*", sep = ""), "*.pdf", ".pdf")
  conf_output(conf, 'ARZ1', paste("\'", paste(name, collapse = "\n"), "\'", sep = ""), FALSE)
  name = rep_tab(paste(indir, "/04.Network/*", sep = ""), "*.html", 2)
  conf_output(conf, 'ARZ2', paste("\'", paste(name, collapse = "\n"), "\'", sep = ""), TRUE)
}

rep_enrich <- function(indir) {
  conf = paste(indir, "/enrich.xls", sep = "")
  name = rep_fig(paste(indir, "/03.Enrichment_pattern", sep = ""), "*_heatmap.pdf", "_heatmap.pdf")
  conf_output(conf, 'ARZ1', paste("\'", paste(name, collapse = "\n"), "\'", sep = ""), FALSE)
  name = rep_tab(paste(indir, "/03.Enrichment_pattern", sep = ""), "*.xls", 1)
  conf_output(conf, 'ARZ2', paste("\'", paste(name, collapse = "\n"), "\'", sep = ""), TRUE)
  name = rep_fig(paste(indir, "/03.Enrichment_pattern", sep = ""), "*_enrichment.pdf", "_enrichment.pdf")
  conf_output(conf, 'ARZ3', paste("\'", paste(name, collapse = "\n"), "\'", sep = ""), TRUE)
  name = rep_fig(paste(indir, "/03.Enrichment_overlap", sep = ""), "*_trend.pdf", "_trend.pdf")
  conf_output(conf, 'ARZ4', paste("\'", paste(name, collapse = "\n"), "\'", sep = ""), TRUE)
  name = rep_fig(paste(indir, "/03.Enrichment_overlap", sep = ""), "*_venn.pdf", "_venn.pdf")
  conf_output(conf, 'ARZ5', paste("\'", paste(name, collapse = "\n"), "\'", sep = ""), TRUE)
  name = rep_tab(paste(indir, "/03.Enrichment_overlap", sep = ""), "*.xls", 1)
  conf_output(conf, 'ARZ6', paste("\'", paste(name, collapse = "\n"), "\'", sep = ""), TRUE)
}

rep_overlap <- function(proid, indir) {
  conf = paste(indir, "/overlap.xls", sep = "")
  conf_output(conf, 'ARZ1', proid, FALSE) 
  name = count = vector()
  invisible(lapply(Sys.glob(file.path(paste(indir, "/01.Expression_pattern", sep = ""), "*")), 
                   function(x) {
                     name <<- c(basename(x), name)
                     count <<- c(length(Sys.glob(file.path(x, "*.txt"))), count)
                   }))
  conf_output(conf, 'ARZ2', length(name), TRUE)
  conf_output(conf, 'ARZ3', paste(name, collapse = ","), TRUE)
  conf_output(conf, 'ARZ4', paste(count, collapse = "个、"), TRUE)
  conf_output(conf, 'ARZ4', paste(count, collapse = "个、"), TRUE)
  name = vector()
  invisible(lapply(Sys.glob(file.path(paste(indir, "/02.Module_overlap", sep = ""), "*.xls")),
                   function(x) {
                     modc = gsub(".xls", "", basename(x))
                     modt = read.table(x, sep = "\t", header = TRUE)
                     name <<- c(paste(modc, nrow(modt), sep = ":"), name)
                   }))
  conf_output(conf, 'ARZ5', paste(name, collapse = "种、"), TRUE)
  name = c('pattern', 'eigengene', 'heatmap')
  invisible(lapply(1:3, 
            function(x) {
              count = rep_fig(paste(indir, "/01.Expression_pattern", sep = ""), paste("*/*_", name[x], ".pdf", sep = ""),
                              paste("_", name[x], ".pdf", sep = ""))
              conf_output(conf, paste('ARZ', x+5, sep = ""), 
                          paste("\'", paste(count, collapse = "\n"), "\'", sep = ""), TRUE)
             }))
  count = rep_tab(paste(indir, "/01.Expression_pattern", sep = ""), "*/*.txt", 1)
  conf_output(conf, 'ARZ9', paste("\'", paste(count, collapse = "\n"), "\'", sep = ""), TRUE)
  name = Sys.glob(file.path(paste(indir, "/02.Module_overlap", sep = ""), "*.xls"))
  conf_output(conf, 'ARZ0', sub(".*[^Omics_analysis]+(Omics_analysis*)", "\\1", name[1]), TRUE)
}

opts <- commandArgs(trailingOnly = T)
if (length(opts) != 3) {
  cat("./Rscript s5.omics_report.R <indir> <proid> <rbin>\n")
  quit(save = "no", status = 1, runLast = TRUE)
}
omics_report(opts[1], opts[2], opts[3])
