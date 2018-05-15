library(dplyr)
library(ggplot2)
library(reshape2)

files <- list.files(pattern = "summary\\.csv$")
dlist <- lapply(files, read.csv)
names <- lapply(files, function(x) gsub("happy_", "", gsub(".summary.csv", "", x)))
dnamed <- mapply(cbind, dlist, "Name"=names, SIMPLIFY=F)
merged <- Reduce(function(...) merge(..., all=T), dnamed)
melted <- melt(merged, id.vars=c("Name", "Filter", "Type"))

metrics <- subset(melted, variable%in%c("METRIC.Recall", "METRIC.F1_Score", "METRIC.Precision"))
ggplot(metrics, aes(x=Name, y=value, color=Filter)) +
  geom_point(stat="identity", position = position_jitter(w = 0.06, h = 0)) +
  geom_text(aes(label=ifelse((Filter=="PASS"), round(value,3), "")), color="black", size=2.5, hjust=-0.4, vjust=0.5) +
  facet_grid( variable ~ Type, scales="free_y" ) +
  theme(axis.text.x=element_text(angle=30, hjust = 1))
ggsave(filename = 'metrics.png')

counts <- subset(melted, variable%in%c("TRUTH.TP", "TRUTH.FN", "QUERY.FP"))
ggplot(counts, aes(x=Name, y=value, color=Filter)) +
  geom_point(stat="identity", position = position_jitter(w = 0.06, h = 0))  +
  geom_text(aes(label=ifelse((Filter=="PASS"), value, "")), color="black", size=2.5, hjust=-0.4, vjust=0.5) +
  facet_grid( variable ~ Type, scales="free_y" ) +
  theme(axis.text.x=element_text(angle=30, hjust = 1))
ggsave(filename = 'counts.png')
