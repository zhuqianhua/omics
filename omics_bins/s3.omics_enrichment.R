omics_enrichment <- function(cfg)
{
  ctb = as.vector(unlist(read.table(cfg)))
  speciesName = ctb[1]
  goAnn = ctb[2]
  prefix = ctb[3]
  pvlCut = as.numeric(ctb[4])
  geneList = rownames(read.table(ctb[5], header = TRUE, row.names = 1))
  if (!file.exists(dirname(prefix)))
    dir.create(dirname(prefix), recursive = TRUE)
  library(GOstats, quietly = TRUE)
  #--- avaiable database in Bioconductor
  if (has.key(tolower(speciesName), annSpe) == TRUE) {
    library(annSpe[[speciesName]], character.only = TRUE, quietly = TRUE)
    go_db = gsub('.db', 'GO', annSpe[[speciesName]])
    go_table = toTable(get(go_db))
    #--- perform kegg enrichment
    library(KEGG.db, quietly = TRUE)
    ko_hgp = new("KEGGHyperGParams", geneIds = geneList, universeGeneIds = NULL,
                 annotation = annSpe[[speciesName]], categoryName = "KEGG",
                 pvalueCutoff = pvlCut, testDirection = "over")
    ko_hpt = summary(hyperGTest(ko_hgp))
    if (nrow(ko_hpt) > 2) {
      write.table(ko_hpt, paste(prefix, "_KEGG.xls", sep = ""),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      if (length(ctb) > 7)
        draw_pattern(paste(prefix, "_KEGG.xls", sep = ""), ctb[7:length(ctb)])
    }
  } else { #--- customized annotation file, GO: gene_id, go_id, Evidence
    go_table = read.table(goAnn, sep = "\t", header = T)
  } 
  go_frame = data.frame(go_table$go_id, go_table$Evidence, go_table$gene_id)
  #--- create GOAllFrame
  library(AnnotationDbi)
  library(RColorBrewer)
  go_af = GOAllFrame(GOFrame(go_frame, organism = speciesName))
  #--- create gene set
  library(GSEABase)
  go_gsc <- GeneSetCollection(go_af, setType = GOCollection())
  if (length(ctb) == 7)
    pattern_pheatmap(ctb[7], ctb[3])
  ogi = paste(dirname(dirname(prefix)), "/00.Intermediate/GO_gene", sep = "")
  if (!dir.exists(ogi))
    dir.create(ogi, recursive = TRUE)
  #--- perform GO enrichment analysis
  for (i in c("BP", "CC", "MF")) {
    go_hgp = GSEAGOHyperGParams(name=paste(speciesName, 'GO', sep = " "),
                                geneSetCollection = go_gsc, geneIds = geneList,
                                universeGeneIds = NULL, ontology = i, 
                                pvalueCutoff = pvlCut, conditional = FALSE, 
                                testDirection = 'over')
    go_hpt = hyperGTest(go_hgp)
    oen = paste(prefix, "_", i, ".xls", sep = "")
    if (nrow(summary(go_hpt)) < 2)
      next
    write.table(summary(go_hpt), oen, sep = "\t", row.names = FALSE, quote = FALSE)
    if (length(ctb) > 7) {  #--- enrimenth for overlap
      draw_pattern(oen, ctb[7:length(ctb)])
    } else {
      etb = read.table(oen, header = TRUE, sep = "\t")
      ptb = read.table(ctb[7], header = TRUE, row.names = 1, sep = "\t")
      if (file.exists(ctb[6]) == TRUE) 
        ltb = read.table(ctb[6], sep = "\t")
      gCut = 10
      mCut = 0
      if (nrow(etb) < gCut)
        gCut = nrow(etb)
      pdf(paste(prefix, "_", i, "_enrichment.pdf", sep = ""), height = 1.25 * gCut)
      tnr = 10;
      if (nrow(etb) < tnr)
        tnr = nrow(etb)
      layout(matrix(1:tnr, ncol = 2))
      lapply(1:nrow(etb), 
             function(x) {
               if (mCut >= gCut)
                 return(0)
               goe = as.vector(as.matrix(etb[x,]))
               gog = geneList[which(geneList %in% go_hpt@goDag@nodeData@data[[goe[1]]]$geneIds)]
               ftb = ptb[which(rownames(ptb) %in% gog),]
               if (file.exists(ctb[6]) == TRUE) 
                 ftb = ptb[which(rownames(ptb) %in% unique(ltb$V1[which(ltb$V2 %in% gog)])),]
               write.table(ftb, paste(ogi, "/", basename(prefix), "_", i, ".xls", sep = ""), sep = "\t")
               if (nrow(ftb) < 3)
                 return(0)
               mCut <<- mCut + 1
               dhcl = hclust(dist(ftb))
               par(mar = c(4, 2, 4, 2), las = 2)
               dcol = 1 / (ncol(ftb) - 1)
               if (nchar(goe[7]) > 35)
                 goe[7] = paste(substr(goe[7], 1, 35), "...", sep = "")
               image(t(as.matrix(ftb[dhcl$order,])), col = brewer.pal(9, "BuGn"), 
                     main = paste(goe[7], paste("Pvalue: ", goe[2], sep = ""), sep = "\n"),
                     axes = FALSE)
               axis(1, seq(from = 1e-5, to = 1 + dcol, by = dcol), labels = names(ftb))
             })
      dev.off()
    }
  }
}

pattern_pheatmap <- function(pat, pfx)
{
  library(pheatmap)
  dat = read.table(pat, header = TRUE, row.names = 1, sep = "\t")
  vsd = vector()
  invisible(lapply(1:nrow(dat), function(x) { if (sd(dat[x,]) == 0) vsd <<- c(vsd, -x) }))
  pdf(paste(pfx, "_heatmap.pdf", sep = ""), width = 8, height = 6, onefile = FALSE)
  if (length(vsd) != 0)
    dat = dat[c(vsd),]
  pheatmap(dat, scale = 'row', cluster_cols = FALSE, show_rownames = F, fontsize_row = .8,
           main = paste("Heatmap of ", basename(pfx), sep = ""))
  dev.off()
}

draw_pattern <- function(ann, pat)
{
  patFig = gsub(".xls$", "_trend.pdf", ann)
  patNum = length(pat)
  width = 10
  if (patNum == 2)
    width = 12
  pdf(patFig, width = width, height = 7)
  library(limma)
  library(RColorBrewer)
  layout(matrix(c(seq(1:patNum), rep(patNum+1, patNum)), patNum, 2, byrow = FALSE))
  lapply(1:patNum, 
         function(x){
           tab = read.table(as.character(pat[x]), sep = "\t", header = TRUE, row.names = 1)
           #dat = as.data.frame(voom(tab)$E)
           dat = data.frame(matrix(rank(tab, ties.method = "first"), nrow(tab), ncol(tab), byrow = FALSE),
                            row.names = row.names(tab))
           names(dat) = names(tab)
           dat = round(dat / round(max(dat) / 10 + 0.4999999, 0) + 0.4999999, 0)
           tit = gsub("M[^M]+$|.txt", "", basename(pat[x]))
           par(mar=c(5,5,3,1), las = 2)
           switch(
             ifelse(x %% 3 == 0, 3, x %% 3),
             barplot(colMeans(dat), main = tit, col = brewer.pal(8, "Set1"),
                     ylab = "Mean of expression", border = NA, cex.axis = 0.8),
             my_heatmap(dat, tit),
             boxplot(dat, range = 5, col = brewer.pal(8, "Dark2"), main = tit,
                     pars = list(boxwex = 0.2, staplewex = 0.5, outwex = 0.5),
                     ylab = "Expression")
           )})
  my_hist(ann)
  dev.off()
}

my_heatmap <- function(data, title)
{
  dhcl = hclust(dist(data))
  dcol = 1 / (ncol(data) - 1)
  image(t(as.matrix(data[dhcl$order,])), col = brewer.pal(9, "Purples"), 
        main = title, axes = F, ylab = "Expression")
  axis(1, seq(from = 1e-5, to = 1 + dcol, by = dcol), labels = names(data))
}

my_hist <- function(enri)
{
  dat = read.table(enri, header = TRUE, sep = "\t")
  num = 30
  if (nrow(dat) < 30)
    num = nrow(dat)
  dat = na.omit(dat[1:num,])
  par(mar = c(1, 1, 3, 25))
  myc = colorRamp(brewer.pal(9, "YlOrRd"), space = "Lab")
  enc = rgb(myc((max(dat$Pvalue) - dat$Pvalue) / max(dat$Pvalue)), max = 255)
  barplot(-dat$Count / max(dat$Count), axes = FALSE, col = enc, 
          main = "Enrichment", adj = 1, space = 0.5, width = 5, 
          las = 2, horiz = T, xlim = c(-1, 0), border = NA)
  myl = lapply(as.vector(dat$Term), 
               function(x){
                 if (nchar(x) > 50)
                   print(paste(substr(x, 1, 50), "...", sep = ""))
                 else
                   print(x)})
  axis(4, seq(5, 7.45 * num, length = num), labels = myl, las = 2, cex = 1)
}

options <- commandArgs(trailingOnly = T)
if (length(options) != 2) {
  cat("./Rscript s1.omics_enrichment.R <enrich> <rlib>\n")
  quit(save = "no", status = 1, runLast = TRUE);
}
.libPaths(options[2])
#--- GOdb information
library(hash)
#--- available database for the model/common species
annKey = c("human", "hsa", "homo sapiens", "yeast", "sce", "saccharomyces cerevisiae",
           "anopheles", "aga", "anopheles gambiae", "chimpanzee", "ptr", "pan troglodytes",
           "rat", "rno", "rattus norvegicus",  "pig", "ssc", "sus scrofa",
           "arabidopsis", "ath", "arabidopsis thaliana", "cow", "bta", "bos taurus",
           "nematode", "cel", "caenorhabditis elegans", "dog", "cfa", "canis familiaris",
           "fruit fly", "dme", "drosophila melanogaster", "zebrafish", "dre", "danio rerio",
           "chicken", "gga", "gallus gallus", "mouse", "mmu", "mus musculus",
           "rhesus monkey", "mcc", "macaca mulatta", "african clawed frog", "xla", "xenopus laevis"
)
annVal = c("org.Hs.eg.db", "org.Sc.sgd.db", "org.Ag.eg.db", "org.Pt.eg.db", 
           "org.Rn.eg.db", "org.Ss.eg.db", "org.At.tair.db", "org.Bt.eg.db", 
           "org.Ce.eg.db", "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db",
           "org.Gg.eg.db", "org.Mm.eg.db", "org.Mmu.eg.db", "org.Xl.eg.db")
annSpe = hash()
.set(annSpe, keys = annKey, values = rep(annVal, each = 3))

omics_enrichment(options[1])

