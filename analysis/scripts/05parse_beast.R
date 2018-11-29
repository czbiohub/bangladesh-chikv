library(ape)
library(ggplot2)
library(ggtree)
library(lubridate)
library(magrittr)
library(countrycode)
library(dplyr)
library(seqinr)
library(coda)
dataenv <- new.env()
load("../parse_raxml_output.RData", envir=dataenv)

output <- new.env()


# summarize output --------------------------------------------------------

lapply(c("log", "trees"), function (x) {
  paste0("/Applications/BEAST\\ 2.5.0/bin/logcombiner -log ../../beast/ncbi_chrf_aln_south_asian_subset_", c("partitioned_", ""), "StrictClock_Skyline_1.", x, 
         " -log ../../beast/ncbi_chrf_aln_south_asian_subset_", c("partitioned_", ""), "StrictClock_Skyline_2.", x, 
         " -burnin 50 -o ../../beast/ncbi_chrf_aln_south_asian_subset_", c("partitioned_", ""), "StrictClock_Skyline_combined.", x) %>%
    sapply(system)
})
paste0("/Applications/BEAST\\ 2.5.0/bin/treeannotator -burnin 50 -heights 'median' ../../beast/ncbi_chrf_aln_south_asian_subset_", c("partitioned_", ""), "StrictClock_Skyline_combined.trees ../../beast/ncbi_chrf_aln_south_asian_subset_", c("partitioned_", ""), "StrictClock_Skyline.mcc.tree") %>%
  sapply(system)

# Skyline -----------------------------------------------------------------

output$skyline_df <- skyline_df <- read.table("../../beast/ncbi_chrf_aln_south_asian_subset_StrictClock_Skyline_skyline.table", header=TRUE, sep="\t", skip=1)
output$skyline_df_partitioned <- skyline_df_partitioned <- read.table("../../beast/ncbi_chrf_aln_south_asian_subset_partitioned_StrictClock_Skyline_skyline.table", header=TRUE, sep="\t", skip=1)

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

ggsave("../../figures/skyline_plot.pdf", skyline_plot, width=6, height=4)

output$skyline_plot_partitioned <- 
  skyline_plot_partitioned <- 
  plotting_func(skyline_df_partitioned)

ggsave("../../figures/skyline_partitioned_plot.pdf", skyline_plot_partitioned, width=6, height=4)


# Estimated rates of evolution --------------------------------------------

output$logfile_partitioned <-
  logfile_partitioned <- 
  read.table("../../beast/ncbi_chrf_aln_south_asian_subset_partitioned_StrictClock_Skyline_combined.log", header=TRUE)
output$logfile <-
  logfile <- 
  read.table("../../beast/ncbi_chrf_aln_south_asian_subset_StrictClock_Skyline_combined.log", header=TRUE)

output$clockrate_plot <- clockrate_plot <-
  select(logfile_partitioned, ends_with("clock.rate")) %>%
  reshape2::melt() %>%
  rbind(., reshape2::melt(select(logfile, clock.rate))) %>%
  ggplot(aes(x=variable, y=value, fill=variable)) %>%
  `+`(list(geom_boxplot(alpha=.4),
           theme_bw(),
           xlab("Rate of substitution per year"),
           ylab("Density"),
           scale_x_discrete(labels=c("non-structural\npolyprotein\n(7424 bases)", "structural\npolyprotein\n(3746 bases)", "non-coding\nregion\n(563 bases)", "whole\ngenome\n(11733 bases)")),
           theme(legend.position="none")
           ))
ggsave("../../figures/clockrate_plot.pdf", clockrate_plot, width=5, height=6)

output$clockrate_table <- clockrate_table <-
  split(clockrate_plot$data$value, clockrate_plot$data$variable) %>% 
  lapply(EpiGenR::hpd) %>%
  lapply(formatC, format="e") %>%
  as.data.frame() %>%
  t()

# Trees -------------------------------------------------------------------

output$beast_tree <- beast_tree <- 
  ggtree::read.beast("../../beast/ncbi_chrf_aln_south_asian_subset_partitioned_StrictClock_Skyline.mcc.tree")

output$beast_tree_df <-
  beast_tree_df <-
  dataenv$mltree_southasian_df %>%
  mutate(seq=gsub("UnitedStatesOfAmerica", "USA", .$seq)%>%gsub("MF773567_Indonesia", "MF773567_Borneo", .))

beast_tree_plot <-
  ggtree(beast_tree, right=TRUE, mrsd=as.character(max(dataenv$mltree_southasian_df$date))) +
  theme_tree2()
beast_tree_plot <- beast_tree_plot +
  xlim(min(beast_tree_plot$data$x)-3, 2045) +
  geom_nodepoint(aes(alpha=posterior))
beast_tree_plot <- beast_tree_plot %<+% beast_tree_df +
  geom_tiplab(aes(color=chrf), align=TRUE)
output$gheatmap_df <- gheatmap_df <-
  data.frame(Continent=beast_tree_df[, "continent"], row.names=beast_tree_df$seq)
beast_tree_plot %<>%
  gheatmap(gheatmap_df, 
           offset=10, width=0.1, colnames=FALSE) %>%
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
beast_tree_plot <- beast_tree_plot + 
  annotate(geom="text", x=bangladesh_root_df$x-1.5, y=bangladesh_root_df$y,
           label=paste0(bangladesh_root_dates[1], "\n(", bangladesh_root_dates[3], " - ", bangladesh_root_dates[2], ")"), size=2.5)
beast_tree_plot <-
  beast_tree_plot + 
  theme(legend.title=element_text(size=12)) +
  scale_color_manual("CHRF", values=c("gray50", "orangered1"), labels=c("Others", "CHRF isolates")) +
  scale_fill_discrete("Continent") +
  scale_alpha_continuous("Probability of \nNode")
output$beast_tree_plot <- beast_tree_plot

ggsave("../../figures/beast_tree_plot.pdf", beast_tree_plot,
       width=12, height=20)


# circular plot -----------------------------------------------------------

beast_tree_radial_plot <- ggtree(beast_tree, right=TRUE, mrsd=as.character(max(dataenv$mltree_southasian_df$date)), layout="circular")
beast_tree_radial_plot <- beast_tree_radial_plot +
  geom_nodepoint(aes(alpha=posterior))
beast_tree_radial_plot <- beast_tree_radial_plot %<+% beast_tree_df +
  geom_tippoint(aes(color=chrf), shape=17)
beast_tree_radial_plot %<>%
  gheatmap(gheatmap_df,
           offset=0.1, width=0.1, colnames=FALSE) %>%
  scale_x_ggtree(breaks=seq(2005, 2017, by=2))
# output$radial_root_coords <- radial_root_coords <- slice(beast_tree_plot$data, which.min(x))
# beast_tree_radial_plot <- beast_tree_radial_plot +
#   annotate(geom="text", x=radial_root_coords$x-1.5, y=radial_root_coords$y,
#            label=paste0(root_dates[1], "\n(", root_dates[3], " - ", root_dates[2], ")"), size=2.5)
beast_tree_radial_plot <-
  beast_tree_radial_plot +
  annotate(geom="text", x=seq(2007, 2017, by=2), y=102, label=seq(2007, 2017, by=2)) + 
  scale_x_continuous(breaks=seq(2007, 2017, by=2)) +
  theme(legend.title=element_text(size=12), 
        panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3)) +
  scale_color_manual("CHRF", values=c("gray50", "orangered1"), labels=c("Others", "CHRF isolates")) +
  scale_fill_discrete("Continent") +
  scale_alpha_continuous("Probability of \nNode")
output$beast_tree_radial_plot <- beast_tree_radial_plot


ggsave("../../figures/beast_tree_radial_plot.pdf", beast_tree_radial_plot,
       width=12, height=12)




write.tree(beast_tree@phylo, file="../beast_tree.nwk")


save(list=names(output), envir=output, file="../parse_beast_output.RData")
