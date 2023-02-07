draw_venn <- function(type, vnTab)
{
  FileTab = vnTab[2:length(vnTab)]
  outPfx = paste(vnTab[1], "_", type, "_venn.pdf", sep = "")
  
  ############Draw a Venn Diagram with Two Sets
  FileTab[1] = paste(as.character(FileTab[1]), "_", type, ".xls", sep = "")
  FileTab[2] = paste(as.character(FileTab[2]), "_", type, ".xls", sep = "")
  
  if(length(FileTab) == 2)
  {
    if (!file.exists(FileTab[1]) || !file.exists(FileTab[2]))
      return(0)

	###Read GO table and Extract Venn labels name
    GOlistA <- read.table(FileTab[1], head=T, sep="\t")
    VennlabelA <- basename(as.character(FileTab[1]))
    VennlabelA <- gsub("_.*", "", VennlabelA)

    GOlistB <- read.table(FileTab[2], head=T, sep="\t")
    VennlabelB <- basename(as.character(FileTab[2]))
    VennlabelB <- gsub("_.*", "", VennlabelB)
    
    ###Extract GO term
    GOtermA <- GOlistA$Term
    GOtermB <- GOlistB$Term
    
    
    ###Contents in every section of Venn diagram 
    commonAB <- intersect(GOtermA, GOtermB) 
    uniqA <- setdiff(GOtermA, commonAB)
    uniqB <- setdiff(GOtermB, commonAB)
    
	
    ###Number in every section of Venn diagram
    commonABnum <- length(commonAB)
    if ( commonABnum == 0 )
	  return(0)		
	
	uniqAnum <- length(uniqA)
    uniqBnum <- length(uniqB)
    

    
    ###Color of Venn diagram
    color <- c("cyan3", "palevioletred1") 
    color_ellipse <- adjustcolor(color, alpha.f = 0.4) 
    color_border <- adjustcolor(color, alpha.f = 1)
    
    pdf(file=outPfx, width=10, height=8)
    par(mar=c(13,6,3,6)+0.1,xpd=TRUE)
    plot(c(-18,18), c(-18,18), type="n", axes=F, xlab="", ylab="", main="")
    draw.ellipse(c(-5,5), c(0,0), c(8,8), c(9,9), col = color_ellipse, border=color_border, lty=1, lwd = 2)
    
    ###Color of numbers in Venn diagram
    num_x <- c(-7,0,7)
    num_y   <- c(0,0,0)
    num <- c(uniqAnum, commonABnum, uniqBnum)
    num_color <- c("turquoise3", "darkorchid2", "palevioletred3")
    text(num_x, num_y, num, cex=1.5, col=num_color)
    
    ###Color of labers in Venn diagram
    label_x <- c(-9,9)
    label_y <- c(-11,-11)
    label_name <- c(VennlabelA, VennlabelB)
    label_color <- c("turquoise3", "palevioletred3")
    text(label_x, label_y, label_name, cex=1.2, col=label_color)
    
    ###Color of first five GO term
    #uniqA
    if (length(uniqA) > 0)
    {
	if (length(uniqA) >= 5)
    	{
      		uniqA <- uniqA[1:5]
    	}
    	uniqA <- lapply(uniqA,
                    function(x){
                      if(nchar(x) > 30){
                        print(paste(substr(x, 1, 30), "...", sep = ""))
                      }
                      else{
                        print(x)
                      }
                    })
    	legend(-23, -16, uniqA, fill="turquoise3", border="turquoise3", cex=1, text.col="black", bty="n")
    }
    
    #commonAB
    if (length(commonAB) > 0)
    {
	if (length(commonAB) >= 5)
    	{
      		commonAB <- commonAB[1:5]
    	}
    	commonAB <- lapply(commonAB,
                       function(x){
                         if(nchar(x) > 30){
                           print(paste(substr(x, 1, 30), "...", sep = ""))
                         }
                         else{
                           print(x)
                         }
                       })
    	legend(-7, -16, commonAB, fill="darkorchid2", border="darkorchid2", cex=1, text.col="black", bty="n")
    }
    
    #uniqB
    if (length(uniqB) > 0)
    {
	if (length(uniqB) >= 5)
    	{
      		uniqB <- uniqB[1:5]
    	}
    	uniqB <- lapply(uniqB,
                    function(x){
                      if(nchar(x) > 30){
                        print(paste(substr(x, 1, 30), "...", sep = ""))
                      }
                      else{
                        print(x)
                      }
                    })
    	legend(9, -16, uniqB, fill="palevioletred3", border="palevioletred3", cex=1, text.col="black", bty="n")
    }

    dev.off()
  }
  
  ############Draw a Venn Diagram with Three Sets
  if(length(FileTab) == 3)
  {
    FileTab[3] = paste(as.character(FileTab[3]), "_", type, ".xls", sep = "")
    if (!file.exists(FileTab[1]) || !file.exists(FileTab[2]) || !file.exists(FileTab[3]))
      return(0)
    ###Read GO table and Extract Venn labels name
    GOlistA <- read.table(FileTab[1], head=T, sep="\t")
    VennlabelA <- basename(as.character(FileTab[1]))
    VennlabelA <- gsub("_.*", "", VennlabelA)
    
    GOlistB <- read.table(FileTab[2], head=T, sep="\t")
    VennlabelB <- basename(as.character(FileTab[2]))
    VennlabelB <- gsub("_.*", "", VennlabelB)
    
    GOlistC <- read.table(FileTab[3], head=T, sep="\t")
    VennlabelC <- basename(as.character(FileTab[3]))
    VennlabelC <- gsub("_.*", "", VennlabelC)
    
    
    ###Extract GO term
    GOtermA <- GOlistA$Term
    GOtermB <- GOlistB$Term
    GOtermC <- GOlistC$Term
    
    
    ###Contents in every section of Venn diagram 
    commonAB <- intersect(GOtermA, GOtermB) 
    commonAC <- intersect(GOtermA, GOtermC) 
    commonBC <- intersect(GOtermB, GOtermC) 
    commonABAC <- intersect(commonAB, commonAC)
    commonABC <- intersect(commonABAC, commonBC)
    
    uniqAB <- setdiff(commonAB, commonABC)
    uniqAC <- setdiff(commonAC, commonABC)
    uniqBC <- setdiff(commonBC, commonABC)
    
    remainA <- setdiff(GOtermA, commonAB)
    uniqA <- setdiff(remainA, uniqAC)
    remainB <- setdiff(GOtermB, commonAB)
    uniqB <- setdiff(remainB, uniqBC)
    remainC <- setdiff(GOtermC, commonAC)
    uniqC <- setdiff(remainC, uniqBC)
    
    ###Number in every section of Venn diagram
    uniqAnum <- length(uniqA)
    uniqBnum <- length(uniqB)
    uniqCnum <- length(uniqC)
    uniqABnum <- length(uniqAB)
    uniqACnum <- length(uniqAC)
    uniqBCnum <- length(uniqBC)
	commonABCnum <- length(commonABC)
    
	if ( uniqABnum == 0 && uniqACnum == 0 && uniqBCnum == 0 && commonABCnum == 0)
		return(0)


    ###Color of Venn diagram
    color <- c("cyan3", "palevioletred1", "goldenrod1") 
    color_ellipse <- adjustcolor(color, alpha.f = 0.4) 
    color_border <- adjustcolor(color, alpha.f = 1)
    
    #outputname <- paste(VennlabelA, "_", VennlabelB, "_", VennlabelC, ".pdf")
    pdf(file=outPfx, width=10, height=8)
    par(mar=c(20,12,0,12)+0.1,xpd=TRUE)
    plot(c(-18,18), c(-18,18), type="n", axes=F, xlab="", ylab="", main="")
    draw.ellipse(c(0,4,-4), c(3.02,-3.912,-3.912), c(10,10,10), c(10,10,10), col = color_ellipse, border=color_border,  angle=c(60,120,0), lty=1, lwd = 2)
    
    ###Color of numbers in Venn diagram
    num_x <- c(0, -7, -8, 0, 0, 7, 8)
    num_y   <- c(10, 2, -8, -1, -10, 2, -8)
    num <- c(uniqAnum, uniqABnum, uniqBnum, commonABCnum, uniqBCnum, uniqACnum, uniqCnum)
    num_color <- c("turquoise3", "green3", "goldenrod2", "red3", "darkorange2", "purple4", "palevioletred3")
    text(num_x, num_y, num, cex=1.5, col=num_color)
    
    ###Color of labers in Venn diagram
    label_x <- c(0, -15, 15)
    label_y <- c(15, -16, -16)
    label_name <- c(VennlabelA, VennlabelB, VennlabelC)
    label_color <- c("turquoise3", "goldenrod2", "palevioletred3")
    text(label_x, label_y, label_name, cex=1.2, col=label_color)
    
    ###Color of first five GO term
    #uniqA
    if (length(uniqA) > 0)
    {
	if (length(uniqA) >= 5)
        {
      		uniqA <- uniqA[1:5]
    	}
    	uniqA <- lapply(uniqA,
                    function(x){
                      if(nchar(x) > 30){
                        print(paste(substr(x, 1, 30), "...", sep = ""))
                      }
                      else{
                        print(x)
                      }
                    })
    	legend(-35, -20, uniqA, fill="turquoise3", border="turquoise3", cex=1, text.col="black", bty="n")
    }
    
	#uniqAB
    if (length(uniqAB) > 0)
    {	
	if (length(uniqAB) >= 5)
    	{
      		uniqAB <- uniqAB[1:5]
    	}
    	uniqAB <- lapply(uniqAB,
                     function(x){
                       if(nchar(x) > 30){
                         print(paste(substr(x, 1, 30), "...", sep = ""))
                       }
                       else{
                         print(x)
                       }
                     })
    	legend(-35, -32, uniqAB, fill="green3", border="green3", cex=1, text.col="black", bty="n")
    }
    
	#uniqB
    if (length(uniqB) > 0)
    {
	if (length(uniqB) >= 5)
    	{
      		uniqB <- uniqB[1:5]
    	}
    	uniqB <- lapply(uniqB,
                    function(x){
                      if(nchar(x) > 30){
                        print(paste(substr(x, 1, 30), "...", sep = ""))
                      }
                      else{
                        print(x)
                      }
                    })
    	legend(-35, -44, uniqB, fill="goldenrod2", border="goldenrod2", cex=1, text.col="black", bty="n")
    }
    
	#commonABC
    if (length(commonABC) > 0)
    {
	if (length(commonABC) >= 5)
    	{
      		commonABC <- commonABC[1:5]
    	}
    	commonABC <- lapply(commonABC,
                        function(x){
                          if(nchar(x) > 30){
                            print(paste(substr(x, 1, 30), "...", sep = ""))
                          }
                          else{
                            print(x)
                          }
                        })
    	legend(-12, -32, commonABC, fill="red3", border="red3", cex=1, text.col="black", bty="n")
    }

    #uniqBC
    if (length(uniqBC) > 0)
    {
	if (length(uniqBC) >= 5)
    	{
      		uniqBC <- uniqBC[1:5]
    	}
    	uniqBC <- lapply(uniqBC,
                     function(x){
                       if(nchar(x) > 30){
                         print(paste(substr(x, 1, 30), "...", sep = ""))
                       }
                       else{
                         print(x)
                       }
                     })
    	legend(11, -20, uniqBC, fill="darkorange2", border="darkorange2", cex=1, text.col="black", bty="n")
    }

    #uniqAC
    if (length(uniqAC) > 0)
    {
	if (length(uniqAC) >= 5)
    	{
      		uniqAC <- uniqAC[1:5]
    	}
    	uniqAC <- lapply(uniqAC,
                     function(x){
                       if(nchar(x) > 30){
                         print(paste(substr(x, 1, 30), "...", sep = ""))
                       }
                       else{
                         print(x)
                       }
                     })
    	legend(11, -32, uniqAC, fill="purple4", border="purple4", cex=1, text.col="black", bty="n")
    }
    
    #uniqC
    if (length(uniqC) > 0)
    {
	if (length(uniqC) >= 5)
    	{
      		uniqC <- uniqC[1:5]
    	}
    	uniqC <- lapply(uniqC,
                    function(x){
                      if(nchar(x) > 30){
                        print(paste(substr(x, 1, 30), "...", sep = ""))
                      }
                      else{
                        print(x)
                      }
                    })
    	legend(11, -44, uniqC, fill="palevioletred3", border="palevioletred3", cex=1, text.col="black", bty="n")
    }
    dev.off()
  }
}

omics_enrich_venn <- function(config)
{
  vTb = as.vector(unlist(read.table(config)))
  if (length(vTb) < 3 || length(vTb) > 4)
    return(0)
  invisible(lapply(c('BP', 'CC', 'MF'), draw_venn, vTb))
}

opts <- commandArgs(trailingOnly = T)
if (length(opts) != 2) {
  cat("./Rscript s4.omics_enrich_venn.R <indir> <rlib>\n")
  quit(save = "no", status = 1, runLast = TRUE)
}
.libPaths(opts[2])
library(plotrix)
library(gplots)
invisible(lapply(Sys.glob(file.path(opts[1], "*/*/venn.txt")), omics_enrich_venn))
