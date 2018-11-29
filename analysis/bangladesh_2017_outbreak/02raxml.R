library(ape)
library(ggplot2)
library(ggtree)
library(lubridate)
library(magrittr)
library(countrycode)
library(dplyr)
library(seqinr)

output <- new.env()

output$mltree_raw <- mltree_raw <- read.tree("../../data/outbreak_e1_aln.fasta.raxml.bestTree")

mltree_df <-
  data.frame(seq=mltree_raw$tip.label) %>%
  mutate(country=sub(".*_ *(.*?) *_.*", "\\1", seq)) %>%
  mutate(continent=countrycode(country, origin='country.name', destination="continent")) %>%
  mutate(date=as.Date(sub(".*_", "", seq))) %>%
  mutate(chrf=grepl("CHRF", seq))
rownames(mltree_df) <- mltree_df$seq
output$mltree_df <- mltree_df

mltree_plot <- ggtree(mltree_raw) + theme_tree2() + xlim(0, 0.05)

mltree_plot <- mltree_plot %<+%
  mltree_df +
  geom_tiplab(aes(colour=chrf), align=TRUE) +
  scale_color_manual(values=c("gray50", "orangered1"), labels=c("Others", "CHRF isolates")) + 
  guides(colour = guide_legend(reverse=TRUE))
mltree_plot %<>%
  gheatmap(data.frame(Continent=mltree_df[, "continent"], row.names=rownames(mltree_df)), 
           offset=.01, width=0.1, colnames=FALSE) %>%
  scale_x_ggtree()
output$mltree_plot <- mltree_plot <-
  mltree_plot +
  theme(legend.title=element_text(size=12)) +
  scale_fill_discrete("Continent") +
  scale_color_manual("CHRF", values=c("gray50", "orangered1"), labels=c("Others", "CHRF isolates"))
ggsave("../../figures/mltree_outbreak_plot.pdf", mltree_plot,
       width=10, height=30)

save(list=names(output), envir = output, file="parse_raxml.RData")
