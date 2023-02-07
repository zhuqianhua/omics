expr_pattern <- function(expr, prefix, missRate = 0.2,
                         topVariant = 5000, betaThreshold = 0.9, minModuleSize = 20)
{
  #--- load expr
  dat = read.table(expr, header = T, row.names = 1, sep = "\t")
  #--- remove samples which expression estimates with counts in less than 20% of cases
  dat = dat[apply(dat, 1, function(x) sum(x == 0)) < ncol(dat) * (1 - missRate),]
  #--- figure output
  if (!file.exists(dirname(prefix)))
    dir.create(dirname(prefix), recursive = TRUE)
  pdf(paste(prefix, "_threshold.pdf", sep = ""), 12, 8)
  
  #--- transform count data to log2-counts per million (logCPM)
  #library(limma)
  #dat_voom = voom(dat)$E
  #--- choose top 5000 most variant genes by MAD
  #if (nrow(dat) < topVariant)
  #  topVariant = nrow(dat)
  #dat_matrix = t(dat_voom[order(apply(dat_voom, 1, mad), decreasing = T)[1:topVariant],])
  #dat_matrix = t(dat[order(apply(dat, 1, mad), decreasing = T)[1:topVariant],])
  dat_matrix = na.omit(t(dat))
  #--- similarity measure between gene profiles (cor or bicor)
  library(WGCNA)
  dat_cor = abs(cor(dat_matrix, use = 'pairwise.complete.obs'))
  #--- assess the scale free topology of the network
  powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
  sft = pickSoftThreshold(dat_matrix, powerVector = powers, verbose = 0)
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab='Soft Threshold (power)', ylab='Scale Free Topology Model Fit,signed R^2',
       type='n', main = paste('Scale independence'));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
       labels = powers, cex = 1, col = 'red'); 
  abline(h = betaThreshold, col = 'red')
  dev.off()
  #--- determine the beta value of the lowest power
  #--- 1) R^2 >= 0.9
  betaValue = 0
  firstSlope = 0
  for (i in 1:length(powers)) {
    if (firstSlope == 0 & sft$fitIndices[i,3] < 0)
      firstSlope = sft$fitIndices[i,1]
    if (-sign(sft$fitIndices[i,3])*sft$fitIndices[i,2] >= betaThreshold & betaValue == 0)
      betaValue = sft$fitIndices[i,1]
  }
  #--- 2) inflection points
  if (betaValue == 0) {
    dat_loess = loess(-sign(sft$fitIndices[,3])*sft$fitIndices[,2]~sft$fitIndices[,1])
    dat_predict = predict(dat_loess, powers)
    dat_diff = c(FALSE, diff(diff(dat_predict)>0)!=0)
    if (length(powers[dat_diff]) > 0)
      betaValue = powers[dat_diff][1]
  }
  #--- 3) worse betaValue
  if (betaValue == 0) {
    #--- 3.1) the first point which slope is less than 0
    if (firstSlope != 0)
      betaValue = firstSlope
    #--- 3.2) the biggest power while the curve is not saturated
    else 
      betaValue = powers[length(powers)]
  }
  print(paste("Using beta value: ", betaValue, sep=""))
  
  #--- calculation of adjacency matrix
  dat_adjacency = dat_cor^betaValue
  rm(dat_cor)
  gc()
  #--- dissimilarity measure 
  dat_diss = 1 - dat_adjacency
  dat_diss[is.na(dat_diss)] = 0
  rm(dat_adjacency) 
  gc()
  #--- create gene tree by average linkage hierarchical clustering 
  dat_tree = hclust(as.dist(dat_diss), method = "average")
  
  #--- module identification using dynamic tree cut algorithm
  dat_module = cutreeDynamic(dendro = dat_tree, distM = dat_diss, deepSplit = 0,
                             pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
  #--- assign module colours
  dat_colour = labels2colors(dat_module)
  #--- plot the dendrogram and corresponding colour bars underneath
  pdf(paste(prefix, "_pattern.pdf", sep = ""), 12, 8)
  plotDendroAndColors(dat_tree, dat_colour, 'Module colours', dendroLabels = FALSE,
                      hang = 0.03, addGuide = TRUE, guideHang = 0.05, main='')
  dev.off()
  #--- calculate eigengenes
  library(ape)
  dat_eigengene = moduleEigengenes(dat_matrix, colors = dat_colour, excludeGrey = FALSE)$eigengenes
  #--- calculate dissimilarity of module eigengenes
  dat_mdiss = 1 - cor(dat_eigengene, use = 'pairwise.complete.obs')
  #--- cluster module eigengenes
  dat_mdiss[is.na(dat_mdiss)] = 0
  dat_mtree = hclust(as.dist(dat_mdiss), method = 'average');
  pdf(paste(prefix, "_eigengene.pdf", sep = ""), 12, 8)
  par(mar=c(2,2,2,2))
  plot.phylo(as.phylo(dat_mtree), type = 'fan', show.tip.label = FALSE, 
             main='Clustering of module eigengenes')
  tiplabels(frame = 'circle', col='black', text=rep('', length(unique(dat_module))), 
            bg = levels(as.factor(dat_colour)), cex=0.9)
  dev.off()
  #--- module eigengenes heatmap
  pdf(paste(prefix, "_heatmap.pdf", sep = ""), 12, 8)
  TOMplot(dat_mdiss, dat_mtree, as.character(table(dat_colour)),
          main="Heatmap of module eigengenes")
  dev.off()
  
  #--- output module genes
  dat_mname = substring(names(dat_eigengene), 3)
  for (i in 1:length(dat_mname)) {
    #--- skip for those modules with less than 2 genes
    if (length(names(as.data.frame(dat_matrix))[dat_colour == dat_mname[i]]) < 2)
      next
    outFile = paste(prefix, "M", dat_mname[i], ".txt", sep="")
    #write.table(names(as.data.frame(dat_matrix))[dat_colour == dat_mname[i]], outFile, quote=F, row.names=F)
    write.table(t(as.data.frame(dat_matrix)[dat_colour == dat_mname[i]]), 
                outFile, sep = "\t", quote = FALSE, row.names = TRUE)
  }
  gc()
}

mod_cycle <- function(modVec, modList)
{
  if (length(modList) == 1) {
    return(paste(rep(modVec, each = length(unlist(modList[1]))), 
                 rep(unlist(modList[1]), times = length(modVec)), sep = "_sModName_"))
  } else {
    mod_cycle(mod_cycle(modVec, modList[1]), modList[2:length(modList)])
  }
}

my_comb <- function(n)
{
  return(unlist(apply(as.matrix(2:n), 1, 
                      function(x) apply(combn(1:n, x), 2, 
                                        function(x) paste(x, collapse = "_sModComb_")))))
}

my_inter <- function(modVec, modNam)
{
  if (length(modNam) == 1) {
    mod_gene = read.table(modNam[1], header = T)
    return(intersect(modVec, unlist(mod_gene)))
  } else {
    return(my_inter(my_inter(modVec, modNam[1]), modNam[2:length(modNam)]))
  }
}

mod_phyper <- function(module, bgGene, outFile, minInterSize, pvlCut, 
                       speNam, goAnn, ePvlCut, interEn)
{
  modComb = unlist(strsplit(module, '_sModName_'))
  modGa = read.table(modComb[1], header = T, row.names = 1)
  modGa = unlist(rownames(modGa))
  modGb = read.table(modComb[2], header = T, row.names = 1)
  interVec = modGa
  lapply(2:length(modComb),
         function(x) {
           mod_gene = read.table(modComb[x], header = TRUE, row.names = 1)
           interVec <<- intersect(interVec, unlist(rownames(mod_gene)))
         })
  #--- threshold: count of inter genes
  if (length(interVec) < minInterSize)
    return(0)
  if (length(modComb) == 2) {
    unionVec = length(unlist(rownames(modGb)))
  } else {
    modAll= unique(unlist(apply(as.matrix(modComb[3:length(modComb)]), 1,
                                function(x) {
                                  mod_gene = read.table(x, header = TRUE, row.names = 1)
                                  unlist(rownames(mod_gene))
                                })))
    unionVec = length(union(modAll, unlist(rownames(modGb))))
    rm(modAll)
    rm(modGb)
    gc()
  }
  interLst = paste(apply(as.matrix(modComb), 1, basename), collapse = "")
  interLst = gsub(".txt$", "", interLst)
  interLst = gsub(".txt", "_", interLst)
  #--- hypergeometric distribution 
  pvalue = phyper(length(interVec), unionVec, bgGene - unionVec, length(unlist(modGa)))
  #--- threshold: p-value
  if (pvalue > pvlCut)
    return(0)
  write.table(t(c(interLst, pvalue, length(interVec))), 
              sep = "\t", outFile, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  #--- perform go enrichment analysis
  if (interEn != 0) {
    library(plyr)
    if (goAnn == "")
      goAnn = 'na'
    ldply(modComb, pre_enrich, ing = interVec, inm = interLst, spe = speNam, 
           ann = goAnn, pct = ePvlCut)
  }
}

pre_enrich <- function(mod, ing, inm, spe, ann, pct) 
{
  pfx = gsub(".txt", "", basename(mod))
  key = gsub("M[^M]+.txt$", '', basename(mod))
  out = paste(dirname(dirname(dirname(mod))), "/00.Intermediate", sep = "")
  opt = paste(out, "/Enrich_pattern/", pfx, sep = "")
  oit = paste(out, "/Enrich_overlap/", inm, sep = "")
  ovn = paste(out, "/Enrich_venn/", inm, sep = "")
  if (!file.exists(opt))
    dir.create(opt, recursive = TRUE)
  if (!file.exists(oit))
    dir.create(oit, recursive = TRUE)
  if (!file.exists(ovn))
    dir.create(ovn, recursive = TRUE)
  fpt = paste(opt, "/enrich.txt", sep = "")
  fit = paste(oit, "/enrich.txt", sep = "")
  fvn = paste(ovn, "/venn.txt", sep = "")
  if (!file.exists(fit)) {
    write.table(c("gene", ing), paste(oit, "/inter.txt", sep = ""),
                           quote = FALSE, col.names = FALSE, row.names = FALSE)
    write.table(c(spe, ann, paste(dirname(out), "/03.Enrichment_overlap/", inm, sep = ""),
                  pct, paste(oit, "/inter.txt", sep = ""), 'na'), fit, quote = FALSE, 
                sep = "\n", col.names = FALSE, row.names = FALSE)
  }
  if (!file.exists(fvn)) 
    write.table(paste(dirname(out), "/03.Enrichment_overlap/", inm, sep = ""), fvn,
                quote = FALSE, sep = "\n", col.names = FALSE, row.names = FALSE)
  lnk = paste(out, "/Trans_", key, "/link.xls", sep = "")
  if (file.exists(lnk)) {
    tpt = paste(dirname(out), "/01.Expression_pattern/", key, "/", basename(mod), sep = "")
    dpt = read.table(tpt, header = TRUE, row.names = 1)
    dlk = read.table(paste(out, "/Trans_", key, "/link.xls", sep = ""), sep =  "\t")
    write.table(dpt[which(rownames(dpt) %in% unique(dlk$V1[which(dlk$V2 %in% ing)])),],
                paste(oit, "/", key, ".txt", sep = ""), sep = "\t", quote = FALSE)
  } else {
    dpt = read.table(mod, header = TRUE, row.names = 1)
    write.table(dpt[which(rownames(dpt) %in% ing),], 
                paste(oit, "/", key, ".txt", sep = ""), sep = "\t", quote = FALSE)
    tpt = mod
    lnk = 'na'
  }
  write.table(paste(oit, "/", key, ".txt", sep = ""), fit, sep = "\n", append = TRUE,
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(paste(dirname(out), "/03.Enrichment_pattern/", pfx, sep = ""), fvn, append = TRUE,
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  if (!file.exists(paste(opt, "/enrich.txt", sep = "")))
    write.table(c(spe, ann, paste(dirname(out), "/03.Enrichment_pattern/", pfx, sep = ""), 
                  pct, mod, lnk, tpt), paste(opt, "/enrich.txt", sep = ""), quote = FALSE, 
                sep = "\n", col.names = FALSE, row.names = FALSE)
}

mod_overlap <- function(modPath, outdir, speNam, goAnn, ePvlCut, interEn, 
                        minInterSize = 10, pvlCut = 0.05)
{
  if (!file.exists(outdir))
    dir.create(outdir, recursive = TRUE)
  #--- get background genes
  bgLst = vector(length = 0)
  bgLst = unlist(apply(as.matrix(seq(2, length(modPath), 2)), 1, 
                       function(x) append(bgLst, as.vector(unlist(modPath[x])))))
  modGA = read.table(bgLst[1], header = T, row.names = 1)
  modAll= unique(unlist(apply(as.matrix(bgLst[2:length(bgLst)]), 1, 
                              function(x) {
                                mod_gene = read.table(x, header = TRUE, row.names = 1)
                                unlist(rownames(mod_gene))
                              })))
  bgGene = length(union(modAll, unlist(rownames(modGA))))
  rm(modAll)
  rm(modGA)
  gc()
  combVec = my_comb(length(modPath) / 2)
  #--- omics combination
  for (i in 1:length(combVec)) {
    splitVec = as.numeric(unlist(strsplit(combVec[i], "_sModComb_")))
    modVec = unlist(modPath[splitVec[1] * 2])
    modLst = as.list(modPath[splitVec[2:length(splitVec)] * 2])
    modNam = paste(modPath[splitVec * 2 - 1], collapse = "_")
    modCmbR= mod_cycle(modVec, modLst)
    #--- module combination
    outFile= paste(outdir, "/", modNam, ".xls", sep="")
    write.table(t(c("Module name", "P-value", "Inter count")), 
                sep = "\t", outFile, quote = FALSE, row.names = FALSE, col.names = FALSE, append = FALSE)
    apply(as.matrix(modCmbR), 1, mod_phyper, bgGene = bgGene, minInterSize = minInterSize,
          speNam = speNam, ePvlCut = ePvlCut, goAnn = goAnn, interEn = interEn,
          pvlCut = pvlCut, outFile = outFile)
    datt = read.table(outFile, sep = "\t")
    if (nrow(datt) == 1) 
      file.remove(outFile)
  }
}

my_trans <- function(path, suffix, link, maxLink = 30)
{
  datLink = read.table(link, sep = "\t")
  #--- choose top 30 target
  library(plyr)
  if (maxLink != 0)
    datLink = ldply(unique(datLink$V1),function(x) head(datLink[datLink$V1 == x,], maxLink))
  datMod = Sys.glob(file.path(path, paste(suffix, "M*.txt", sep = "")))
  outdir = paste(dirname(dirname(path)), "/00.Intermediate/Trans_", suffix, sep = "")
  if (!file.exists(outdir))
    dir.create(outdir, recursive = TRUE)
  file.copy(link, paste(outdir, "/link.xls", sep = ""), overwrite = TRUE)
  apply(as.matrix(datMod), 1, 
        function(x) {
          datLst = read.table(x, header = TRUE, sep = "\t")
          outFile = paste(outdir, "/", basename(x), sep = "")
          write.table(unique(datLink$V2[which(datLink$V1 %in% row.names(datLst))]),
                      row.names = FALSE, outFile, quote = FALSE)
        })
  
}

#--- perform omics analysis
omics_overlap <- function(config, outdir)
{
  #--- reading config file
  library(hash)
  datPara = hash('speciesabb' = 'unk', 'outdir' = outdir, 'missrate' = 0.2,
                 'topvariant' = 500000, 'betathreshold' = 0.8, 'minmodulesize' = 20,
                 'promodule' = 1,
                 'minintersize' = 10, 'mpvlcut' = 0.05, 'maxlink' = 30,
                 'interenri' = 0, 'goprefix' = '', 'epvlcut' = 0.5,
                 'queue' = 'bc.q', 'project' = 'Project', 'subpro' = '', 'mail' = 'user@genomics.cn',
                 'cdtsuser' = 'user', 'cdtspasswd' = 'passwd', 'cdtsloc' = '',
                 'slicenum' = 20
  )
  datExpr = vector()
  datFile = file(config, "r")
  datLine = readLines(datFile, n = 1, warn = FALSE)
  while (length(datLine) != 0) {
    lineVec = unlist(strsplit(datLine, '(\\s+)?=(\\s+)?'))
    if (length(lineVec) == 2 & length(grep('#', as.character(lineVec[1]))) == 0) {
      if (has.key(tolower(lineVec[1]), datPara) == TRUE)
        datPara[tolower(lineVec[1])] = lineVec[2]
      else
        datExpr = append(datExpr, c(lineVec[1], gsub('\\s+$', '', lineVec[2])))
    }
    datLine = readLines(datFile, n = 1, warn = FALSE)
  }
  close(datFile)
  #--- perform expression pattern analysis
  if (file.exists(paste(datPara[['outdir']], "/00.Intermediate", sep = "")))
    unlink(paste(datPara[['outdir']], "/00.Intermediate", sep = ""), recursive = TRUE)
  datOmics = list()
  for (x in seq(1, length(datExpr)-1, 2)) {
    exprVec = unlist(strsplit(datExpr[x+1], ','))
    outPfx = paste(datPara[['outdir']], "/01.Expression_pattern/", datExpr[x], "/", datExpr[x], sep = "")
    if (as.numeric(datPara[['promodule']]) == 1)
      expr_pattern(exprVec[1], outPfx, as.numeric(datPara[['missrate']]), 
                   as.numeric(datPara[['topvariant']]), as.numeric(datPara[['betathreshold']]),
                   as.numeric(datPara[['minmodulesize']]))
    modPfx = paste(datExpr[x], "M*.txt", sep = "")
    if (length(exprVec) == 2) {
      #--- transfer to gene id
      my_trans(dirname(outPfx), datExpr[x], exprVec[2], as.numeric(datPara[['maxlink']]))
      outPfx = paste(datPara[['outdir']], "/00.Intermediate/Trans_", datExpr[x], "/", datExpr[x], sep = "")
    }
    datOmics = c(datOmics, list(datExpr[x], Sys.glob(file.path(dirname(outPfx), modPfx))))
  }
  #--- caculate overlap
  mod_overlap(datOmics, paste(datPara[['outdir']], "/02.Module_overlap", sep = ""),
              datPara[['speciesabb']], datPara[['goprefix']], 
              as.numeric(datPara[['epvlcut']]), as.numeric(datPara[['interenri']]),
              as.numeric(datPara[['minintersize']]), as.numeric(datPara[['mpvlcut']]))
}

opts <- commandArgs(trailingOnly = T)
if (length(opts) != 3) {
  cat("./Rscript s1.omics_overlap.R <config> <rlib> <outdir>\n")
  quit(save = "no", status = 1, runLast = TRUE);
}
.libPaths(opts[2])
omics_overlap(opts[1], opts[3])

