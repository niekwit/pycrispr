###BASED ON BOX 2 OF HANNA AND DOENCH (2020 Nature Gen.)

#redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

#load required libraries
library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(viridis)
library(yaml)

#open count table
counts <- read.delim(snakemake@input[[1]])

#get replicate information
replicates <- snakemake@params[["reps"]]

#get control condition
control <- snakemake@params[["control"]]

#get sgrna_summary.txt files
sgrna_summaries <- Sys.glob("mageck/*/*.sgrna_summary.txt")

#create output dir
dir.create("qc/replicate-correlation/", showWarnings = FALSE)

#create scatter plot and calculate correlation (Spearman) between each replicate pair
for (r in replicates){
  
  #get each replicate sample
  r1 <- str_split(r,"-")[[1]][1]
  r2 <- str_split(r,"-")[[1]][2]
  
  #get MAGeCK sgrna_summary.txt for each sample versus control
  r1_control <- sgrna_summaries[grepl(paste0(r1,"_vs_",control), sgrna_summaries)]
  r1_control <- r1_control[grep("-", r1_control, invert = TRUE)] #remove file(s) with replicate
  
  r2_control <- sgrna_summaries[grepl(paste0(r2,"_vs_",control), sgrna_summaries)]
  r2_control <- r2_control[grep("-", r2_control, invert = TRUE)] #remove file(s) with replicate
  
  #samples <- c(r1,r2)
  
  #load sgrna data
  column1 <- paste0(r1," LFC from control")
  df1 <- read.delim(r1_control) %>%
    select(c("sgrna",LFC)) %>%
    rename(!!quo_name(column1) := "LFC")
    
  column2 <- paste0(r2," LFC from control")
  df2 <- read.delim(r2_control) %>%
    select(c("sgrna",LFC)) %>%
    rename(!!quo_name(column2) := "LFC")
  
  #create one df with lfc of both samples
  df <- full_join(df1,df2, by="sgrna")
  
  #calculate density colours
  df$dens <- densCols(df, colramp = colorRampPalette(viridis(15)))
  
  #create lm for plotting regression line
  df.lm <- lm(df[[column2]] ~ df[[column1]], df)
  
  #calculate correlation
  c <- format(round(cor(df[[column1]],df[[column2]], method = "spearman"),3), nsmall = 3)
  
 
  
  p <- ggplot(df, aes(x=.data[[column2]], 
                      y=.data[[column1]],
                      col = dens)) +
    geom_point() +
    scale_color_identity() +
    theme_bw(base_size = 18) +
    annotate("text", 
             x = 0.75 * max(df[[column2]]),
             y = 0.85 * min(df[[column1]]),
             label = paste0("r = ",c),
             size = 5) +
    geom_abline(slope = coef(df.lm)[[2]],
                intercept = coef(df.lm)[[1]])
  
  ggsave(paste0("qc/replicate-correlation/",r,".pdf"), p)
  
}

sink(log, type = "output")
sink(log, type = "message")

