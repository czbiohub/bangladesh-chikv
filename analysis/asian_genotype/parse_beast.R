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


# Skyline -----------------------------------------------------------------

output$skyline_df_relaxed <- skyline_df_relaxed <- 
  read.table("../../data/ncbi_chrf_aln_south_asian_skyline.table", header=TRUE, sep="\t", skip=1)
output$skyline_df_strict <- skyline_df_strict <- 
  read.table("../../data/ncbi_chrf_aln_south_asian_skyline_strictclock.table", header=TRUE, sep="\t", skip=1)
output$skyline_df <- skyline_df <-
  rbind(data.frame(skyline_df_relaxed, clock="relaxed"),
        data.frame(skyline_df_strict, clock="strict"))

output$plotting_func <-
  plotting_func <- function (sky_df, cols=NULL, legend.title=NULL) {
    if (is.null(cols)) cols <- c("mediumorchid3", "steelblue2")
    P <- ggplot(sky_df, aes(x=Time, y=Median, ymin=Lower, ymax=Upper)) + 
      theme_bw() +
      ylab("Pathogen diversity")
    if ("clock" %in% names(sky_df)) {
      if (is.null(legend.title)) legend.title <- "Molecular clock"
      P <- P +
        geom_ribbon(aes(fill=clock), alpha=.4) +
        geom_line(aes(color=clock)) +
        scale_colour_manual(legend.title, values=cols) +
        scale_fill_manual(legend.title, values=cols)
    } else {
      P <- P + geom_ribbon(alpha=.4, fill=cols[1]) + geom_line(colour=cols[1])
    }
    P <- P +
      scale_x_continuous(breaks=seq(2005, 2017, by=2))
    P
  }
output$skyline_plot <- 
  skyline_plot <- 
  plotting_func(skyline_df)
output$skyline_plot_relaxed <- 
  skyline_plot_relaxed <- 
  plotting_func(skyline_df_relaxed)
output$skyline_plot_strict <- 
  skyline_plot_strict <- 
  plotting_func(skyline_df_strict)

ggsave("../../figures/skyline_plot.pdf", skyline_plot, width=6, height=4)
ggsave("../../figures/skyline_plot_relaxed.pdf", skyline_plot_relaxed, width=4.5, height=3)
ggsave("../../figures/skyline_plot_strict.pdf", skyline_plot_strict, width=4.5, height=3)


# Trees -------------------------------------------------------------------

output$beast_tree_relaxed <- beast_tree_relaxed <- 
  ggtree::read.beast("../../data/ncbi_chrf_aln_south_asian_skyline.mcc.tree")
output$beast_tree_strict <- beast_tree_strict <- 
  ggtree::read.beast("../../data/ncbi_chrf_aln_south_asian_skyline_strictclock.mcc.tree")


beast_tree_df <-
  dataenv$mltree_southasian_df
beast_tree_df$seq %<>%
  strsplit("_") %>%
  sapply(function (x) {
    if (x[1]=="CHRF") return (paste(x, collapse="_"))
    x[1]
  }) %>%
  sapply(grep, beast_tree_relaxed@phylo$tip.label, value=TRUE)
output$beast_tree_df <- beast_tree_df

beast_tree_plot_relaxed <-
  ggtree(beast_tree_relaxed, right=TRUE, mrsd=as.character(max(dataenv$mltree_southasian_df$date))) +
  theme_tree2()
beast_tree_plot_relaxed <- beast_tree_plot_relaxed +
  xlim(min(beast_tree_plot_relaxed$data$x)-3, 2045) +
  geom_nodepoint(aes(alpha=posterior))
beast_tree_plot_relaxed <- beast_tree_plot_relaxed %<+% beast_tree_df +
  geom_tiplab(aes(color=chrf), align=TRUE)
beast_tree_plot_relaxed %<>%
  gheatmap(data.frame(Continent=beast_tree_df[, "continent"], row.names=beast_tree_df$seq), 
           offset=10, width=0.1, colnames=FALSE) %>%
  scale_x_ggtree(breaks=seq(2005, 2017, by=2))
root_dates <-
  slice(beast_tree_plot_relaxed$data, which.max(height_median)) %>%
  with(c(height_median, height_median-c(height-height_0.95_HPD[[1]]))) %>%
  `*`(365) %>%
  `-`(max(beast_tree_plot_relaxed$data$date, na.rm=TRUE), .) %>%
  format("%b %Y")
root_coords <- slice(beast_tree_plot_relaxed$data, which.min(x))
beast_tree_plot_relaxed <- beast_tree_plot_relaxed + 
  annotate(geom="text", x=root_coords$x-1.5, y=root_coords$y,
           label=paste0(root_dates[1], "\n(", root_dates[3], " - ", root_dates[2], ")"), size=2.5)
beast_tree_plot_relaxed <-
  beast_tree_plot_relaxed + 
  theme(legend.title=element_text(size=12)) +
  scale_color_manual("CHRF", values=c("gray50", "orangered1"), labels=c("Others", "CHRF isolates")) +
  scale_fill_discrete("Continent") +
  scale_alpha_continuous("Probability of \nNode")
output$beast_tree_plot_relaxed <- beast_tree_plot_relaxed

ggsave("../../figures/beast_tree_plot_relaxed.pdf", beast_tree_plot_relaxed,
       width=12, height=20)

save(list=names(output), envir=output, file="parse_beast_output.RData")
