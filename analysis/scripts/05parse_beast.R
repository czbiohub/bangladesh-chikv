library(ape)
library(ggplot2)
library(ggtree)
library(lubridate)
library(magrittr)
library(countrycode)
library(dplyr)
library(seqinr)
library(coda)
library(grid)
library(gridExtra)
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

output$root_dates <- root_dates <-
  slice(beast_tree@data, which.max(height_median)) %>%
  with(c(height_median, height_median-c(height-height_0.95_HPD[[1]]))) %>%
  `*`(365) %>%
  `-`(max(beast_tree_plot$data$date, na.rm=TRUE), .) %>%
  format("%b %Y")

output$root_coords <- root_coords <- slice(beast_tree@data, which.min(x))
output$bangladesh_root_df <- bangladesh_root_df <-
  filter(beast_tree@data, grepl("Bangladesh_2017", label)) %>%
  select("parent") %>%
  unlist() %>%
  slice(beast_tree@data, .) %>%
  slice(which.max(height_median))
output$bangladesh_root_dates <- bangladesh_root_dates <-
  bangladesh_root_df %>%
  with(c(height_median, height_median-c(height-height_0.95_HPD[[1]]))) %>%
  `*`(365) %>%
  `-`(max(beast_tree@data$date, na.rm=TRUE), .) %>%
  format("%b %Y")

region_colours <- c("America"="darkorchid4", 
                    "Africa"="khaki3",
                    "Europe"="deepskyblue2",
                    "Western Pacific"="orange", 
                    "South-East Asia"="springgreen3",
                    "Bangladesh"="red3")
output$beast_tree_df <-
  beast_tree_df <-
  dataenv$mltree_southasian_df %>%
  mutate(seq=gsub("UnitedStatesOfAmerica", "USA", .$seq)%>%gsub("MF773567_Indonesia", "MF773567_Borneo", .)) %>%
  rename(Region=continent, id=seq) %>%
  mutate(acc=strsplit(id, "_")%>%sapply(head, 1)) %>%
  mutate(countries=gsub("United States Of America", "USA", countries) %>%
           gsub("Papua New Guinea", "Papua N G", .)) %>%
  mutate(Region=apply(., 1, function (x) ifelse(x["countries"]=="Bangladesh", "Bangladesh", x["Region"]))) %>%
  mutate(Region=gsub("Americas", "America", Region)) %>%
  mutate(Region=gsub("Oceania", "Western Pacific", Region)) %>%
  mutate(Region=gsub("Asia", "South-East Asia", Region)) %>%
  mutate(Region=factor(Region, levels=names(region_colours)))
output$gheatmap_df <- gheatmap_df <-
  data.frame(Region=factor(beast_tree_df[, "Region"], levels=names(region_colours)),
             row.names=beast_tree_df$id)


# Tree plotting -----------------------------------------------------------

beast_tree_plot_base <-
  ggtree(beast_tree, right=TRUE, mrsd=as.character(max(dataenv$mltree_southasian_df$date))) +
  theme_tree2() +
  geom_tiplab(aes(label=""), align=TRUE)
beast_tree_plot_base %<>%
  gheatmap(gheatmap_df, width=0.06, colnames=FALSE) %>%
  scale_x_ggtree(breaks=seq(2005, 2017, by=2), labels=c(2005, "", 2009, "", 2013, "", 2017))
beast_tree_plot <- beast_tree_plot_base + 
  geom_nodepoint(aes(color=posterior), size=6, shape="bullet", alpha=.75)
beast_tree_plot <- beast_tree_plot +
  annotate(geom="text", x=c(root_coords$x-1.5, bangladesh_root_df$x-1.5),
           y=c(root_coords$y, bangladesh_root_df$y),
           label=c(paste0(root_dates[1], "\n(", root_dates[3], "\n-\n", root_dates[2], ")"),
                   paste0(bangladesh_root_dates[1], "\n(", bangladesh_root_dates[3], "\n-\n", bangladesh_root_dates[2], ")")), 
           size=4.5, lineheight=0.75)
beast_tree_plot <-
  beast_tree_plot +
  coord_cartesian(xlim=c(2002, 2018.8)) +
  scale_fill_manual(name="Region", values=region_colours, breaks=names(region_colours)) +
  scale_color_continuous(name="Probability of\nNode", low="grey99", high="grey20", breaks=seq(0, 1, by=0.25)) +
  theme(axis.text.x=element_text(size=14),
        axis.ticks.length=unit(5, "pt"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14, hjust=0),
        legend.position=c(0.2, 0.15),
        plot.margin=unit(rep(0.1, 4), units="lines"))

get_tree_tiplabs <- function (x, xlim_min, xlim_max) {
  x <- enquo(x)
  beast_tree_plot_base %<+% beast_tree_df +
    geom_tiplab(aes(label=!!x, color=chrf), align=TRUE, offset=1.9) +
    scale_color_manual("CHRF", values=c("gray50", "orangered1"), labels=c("Others", "CHRF isolates")) +
    theme(legend.position="none", plot.margin=unit(c(0.1, 0, 1.45, 0), units="lines")) +
    coord_cartesian(xlim=c(xlim_min, xlim_max))
}
beast_tree_plot_names <- arrangeGrob(get_tree_tiplabs(acc, 2020.19, 2025.05),
                                     get_tree_tiplabs(year, 2019.95, 2021.5),
                                     get_tree_tiplabs(countries, 2020.29, 2029),
                                     nrow=1, widths=c(1.6, 1, 2))
output$beast_tree_plot <- output_beast_tree_plot <- 
  arrangeGrob(beast_tree_plot, beast_tree_plot_names, nrow=1, widths=c(2.5, 1),
              padding=unit(0, "line"))
ggsave("../../figures/beast_tree_plot.pdf", output_beast_tree_plot,
       width=10, height=16)


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
  scale_fill_discrete("Region") +
  scale_alpha_continuous("Probability of \nNode")
output$beast_tree_radial_plot <- beast_tree_radial_plot


ggsave("../../figures/beast_tree_radial_plot.pdf", beast_tree_radial_plot,
       width=12, height=12)




write.tree(beast_tree@phylo, file="../beast_tree.nwk")


save(list=names(output), envir=output, file="../parse_beast_output.RData")
