library(ape)
library(ggplot2)
library(ggtree)
library(lubridate)
library(magrittr)
library(countrycode)
library(dplyr)
library(seqinr)

dataenv <- new.env()
load("parse_raxml_output.RData", envir=dataenv)

output <- new.env()


# summarize output --------------------------------------------------------

lapply(c("log", "trees"), function (x) {
  paste0("/Applications/BEAST\\ 2.5.0/bin/logcombiner -log ../../beast/outbreak_e1_aln_1.", x, 
         " -log ../../beast/outbreak_e1_aln_2.", x, 
         " -burnin 95 -o ../../beast/outbreak_e1_aln_combined.", x) %>%
    sapply(system)
})
paste0("/Applications/BEAST\\ 2.5.0/bin/treeannotator -burnin 95 -heights 'median' ../../beast/outbreak_e1_aln_combined.trees ../../beast/outbreak_e1_aln_combined.mcc.tree") %>%
  sapply(system)


# Trees -------------------------------------------------------------------

output$beast_tree <- beast_tree <- 
  ggtree::read.beast("../../beast/outbreak_e1_aln_combined.mcc.tree")

output$beast_tree_df <-
  beast_tree_df <-
  data.frame(seq=beast_tree@phylo$tip.label) %>%
  mutate(country=sub(".*_ *(.*?) *_.*", "\\1", seq)) %>%
  mutate(continent=countrycode(country, origin='country.name', destination="continent")) %>%
  mutate(date=as.Date(sub(".*_", "", seq))) %>%
  mutate(chrf=grepl("CHRF", seq))

beast_tree_plot <-
  ggtree(beast_tree, right=TRUE, mrsd=as.character(max(beast_tree_df$date))) +
  theme_tree2()
beast_tree_plot <- beast_tree_plot +
  xlim(min(beast_tree_plot$data$x)-3, 2018) +
  geom_nodepoint(aes(alpha=posterior))
beast_tree_plot <- beast_tree_plot %<+% beast_tree_df +
  geom_tiplab(aes(color=chrf), align=TRUE)
output$gheatmap_df <- gheatmap_df <-
  data.frame(Continent=beast_tree_df[, "continent"], row.names=beast_tree_df$seq)
beast_tree_plot %<>%
  gheatmap(gheatmap_df, 
           offset=2, width=0.1, colnames=FALSE) %>%
  scale_x_ggtree(breaks=seq(2005, 2017, by=2))
output$root_dates <- root_dates <-
  slice(beast_tree_plot$data, which.max(height_median)) %>%
  with(c(height_median, height_median-c(height-height_0.95_HPD[[1]]))) %>%
  `*`(365) %>%
  `-`(max(beast_tree_plot$data$date, na.rm=TRUE), .) %>%
  format("%b %Y")
output$root_coords <- root_coords <- slice(beast_tree_plot$data, which.min(x))
beast_tree_plot <- beast_tree_plot + 
  annotate(geom="text", x=root_coords$x-1.5, y=root_coords$y,
           label=paste0(root_dates[1], "\n(", root_dates[3], " - ", root_dates[2], ")"), size=2.5)
output$bangladesh_root_df <- bangladesh_root_df <-
  filter(beast_tree_plot$data, grepl("Bangladesh_2017", label)) %>%
  select("parent") %>%
  unlist() %>%
  slice(beast_tree_plot$data, .) %>%
  slice(which.max(height_median))
output$bangladesh_root_dates <- bangladesh_root_dates <-
  bangladesh_root_df %>%
  with(c(height_median, height_median-c(height-height_0.95_HPD[[1]]))) %>%
  `*`(365) %>%
  `-`(max(beast_tree_plot$data$date, na.rm=TRUE), .) %>%
  format("%b %Y")
# beast_tree_plot <- beast_tree_plot + 
#   annotate(geom="text", x=bangladesh_root_df$x-1.5, y=bangladesh_root_df$y,
#            label=paste0(bangladesh_root_dates[1], "\n(", bangladesh_root_dates[3], " - ", bangladesh_root_dates[2], ")"), size=2.5)
beast_tree_plot <-
  beast_tree_plot + 
  theme(legend.title=element_text(size=12)) +
  scale_color_manual("CHRF", values=c("gray50", "orangered1"), labels=c("Others", "CHRF isolates")) +
  scale_fill_discrete("Continent") +
  scale_alpha_continuous("Probability of \nNode")
output$beast_tree_plot <- beast_tree_plot

ggsave("../../figures/beast_tree_outbreak_e1_plot.pdf", beast_tree_plot,
       width=12, height=20)

save(list=names(output), envir=output, file="parse_beast_output.RData")
