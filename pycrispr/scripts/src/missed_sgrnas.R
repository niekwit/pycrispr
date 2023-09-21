#redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(ggplot2)
library(reshape2)
library(tibble)
library(yaml)

#load count table
counts <- read.delim(snakemake@input[[1]]) %>%
  dplyr::select(-c("sgRNA","gene"))

#get number of sgRNAs with zero counts
df <- colSums(counts==0) %>%
  melt() %>%
  rownames_to_column(var = "sample") 

#load plot settings
settings <- read_yaml("envs/plot_settings.yaml",readLines.warn=FALSE)
theme <- paste0(settings["ggplot2_theme"][[1]],"(base_size = ",settings["font_size"][[1]],")")
fill <- settings["bar_graph_fill_colour"][[1]]
colour <- settings["bar_graph_line_colour"][[1]]

p <- ggplot(data=df, aes(x=sample, y=value)) + 
  geom_bar(stat="identity",
           fill=fill,
           colour=colour) +
  eval(parse(text=as.character(theme))) +
  xlab(NULL) +
  ylab("Missed sgRNAs") +
  scale_x_discrete(guide = guide_axis(angle = 45))

#save plot
ggsave(snakemake@output[[1]], p)

sink(log, type = "output")
sink(log, type = "message")




