#redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(ggplot2)
library(ineq)
library(yaml)


#load count table
counts <- read.delim(snakemake@input[[1]])

#create df to store Gini indices
df <- as.data.frame(matrix(data = NA,
                           ncol = 2,
                           nrow = ncol(counts) - 2))
names(df) <- c("sample","Gini_index")
df$sample <- names(counts)[3:ncol(counts)]

#calculate Gini index for each sample
for (i in seq_len(nrow(df))){
  
  df[i,"Gini_index"] <- Gini(as.vector(counts[i+2])[[1]])
}

#load plot settings
settings <- read_yaml("envs/plot_settings.yaml",readLines.warn=FALSE)
theme <- paste0(settings["ggplot2_theme"][[1]],"(base_size = ",settings["font_size"][[1]],")")
fill <- settings["bar_graph_fill_colour"][[1]]
colour <- settings["bar_graph_line_colour"][[1]]

#plot Gini index bar graph
p <- ggplot(data=df, aes(x=sample, y=Gini_index)) + 
  geom_bar(stat="identity",
           fill=fill,
           colour=colour) +
  eval(parse(text=as.character(theme))) +
  ylab("Gini index") +
  xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Evenness of sgRNAs") +
  scale_x_discrete(guide = guide_axis(angle = 45))

#save plot
ggsave(snakemake@output[[1]], p)


sink(log, type = "output")
sink(log, type = "message")

