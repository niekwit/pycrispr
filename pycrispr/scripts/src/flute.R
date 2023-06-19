library(MAGeCKFlute)
library(clusterProfiler)
library(ggplot2)

#path to gene summary file 
file1 = snakemake@input[[1]]

#path to sgRNA summary file 
file2 = snakemake@input[[2]]

#load data
gdata <- ReadRRA(file1)
gdata$LogFDR <- -log10(gdata$FDR)

#volcano plot
pdf(snakemake@output[[1]])
VolcanoView(gdata, x = "Score", y = "FDR", Label = "id")
dev.off()

#dot plots
gdata$RandomIndex <- sample(1:nrow(gdata), nrow(gdata))
gdata <- gdata[order(-gdata$Score), ]

gg <- gdata[gdata$Score > 0, ]
pdf(snakemake@output[[2]])
ScatterView(gg, x = "RandomIndex", y = "Score", label = "id",
                 y_cut = CutoffCalling(gdata$Score,2), 
                 groups = "top", top = 10, ylab = "Log2FC")
dev.off()

gg <- gdata[gdata$Score < 0, ]
pdf(snakemake@output[[3]])
ScatterView(gg, x = "RandomIndex", y = "Score", label = "id",
                 y_cut = CutoffCalling(gdata$Score,2), 
                 groups = "bottom", top = 10, ylab = "Log2FC")
dev.off()

#sgRNA rank plot
sdata <- ReadsgRRA(file2)

pdf(snakemake@output[[4]])
sgRankView(sdata, top = 5, bottom = 5)
dev.off()

#enrichment analysis
tryCatch( #analysis can fail if a small gRNA library is used, so skip if any error occurs
 {
    geneList <- gdata$Score
    names(geneList) <- gdata$id
    
    enrich_pos <- EnrichAnalyzer(geneList = geneList[geneList > 0.5], 
                                 method = "GSEA", 
                                 type = "KEGG",
                                 organism = snakemake@params[[1]],
                                 filter = TRUE)
    
    pdf(snakemake@output[[5]])
    EnrichedView(enrich_pos, mode = 1, top = 10, bottom = 0)
    dev.off()
    
    enrich_neg <- EnrichAnalyzer(geneList = geneList[geneList < -0.5], 
                                 method = "GSEA", 
                                 type = "KEGG",
                                 organism = snakemake@params[[1]],
                                 filter = TRUE)
    
    pdf(snakemake@output[[6]])
    EnrichedView(enrich_neg, mode = 1, top = 0, bottom = 10)
    dev.off()
  },
  error = function(error_message) {
    message("Skipping enrichment analysis: error was raised by R (most likely caused by using smaller gRNA libraries)")
    message("This is the R error message:")
    message(error_message)
    return(NA)
  }
)






