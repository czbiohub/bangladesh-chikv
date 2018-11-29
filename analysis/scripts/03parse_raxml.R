library(ape)
library(ggplot2)
library(ggtree)
library(lubridate)
library(magrittr)
library(countrycode)
library(dplyr)
library(seqinr)

output <- new.env()

dataenv <- new.env()
load("curate_blast_hits.RData", dataenv)
output$outgroup_names <- outgroup_names <- dataenv$blast_hits %>% 
  dplyr::slice(., sapply(dataenv$genbank_results_filtered_sampled$id, grep, .$sseqid)) %>%
  dplyr::filter(., qcov*pident < 85) %>% 
  `[[`("sseqid") %>% 
  gsub(".1", "", ., fixed=TRUE)


# Plot tree with all NCBI sequences ---------------------------------------------
mltree_raw <- read.nexus("../../data/ncbi_chrf_aln.fasta.raxml.bstree")
mltree_raw$node.label <- 
  slice(read.beast("../../data/ncbi_chrf_aln.fasta.raxml.bstree")@data, order(node))$bootstrap
output$mltree_raw <- mltree_raw %<>% root(., sapply(outgroup_names, grep, .$tip.label, value=TRUE))
strsplit(mltree_raw$tip.label, "_") %>%
  sapply(function (x) {
    paste(x[c(-1, -length(x))], format(as.Date(tail(x, 1)), "%Y"))
  }) %>%
  duplicated %>%
  which %>%
  drop.tip(mltree_raw, .) %>% 
  ggtree() %>%
  `+`(list(theme_tree2(),
           xlim(0, 0.34))) -> mltree_raw_plot
output$mltree_raw_df <-
  mltree_raw_df <- data.frame(seq=mltree_raw_plot$data$label[mltree_raw_plot$data$isTip], stringsAsFactors=FALSE) %>%
  mutate(date=as.Date(strsplit(seq, "_")%>%sapply(tail, 1))) %>%
  mutate(date_=lapply(date, function (x) c(as.Date("2000-01-01"), x)[ifelse(x<as.Date("2000-01-01"), 1, 2)])%>%do.call(what=c)) %>%
  mutate(seq=seq%>%substr(., 1, nchar(.)-6))
mltree_raw_plot$data$label <- mltree_raw_plot$data$label %>% substr(., 1, nchar(.)-6)
mltree_raw_plot <- mltree_raw_plot %<+%
  mltree_raw_df +
  geom_tiplab(aes(color=date_)) +
  theme(legend.position="right", legend.direction = "vertical",
        text=element_text(size=8)) +
  scale_colour_date(name="date", 
                    breaks=as.Date(c("2000-01-01", "2005-01-01", "2010-01-01", "2015-01-01")),
                    labels=c("On or before\n2000", "2005", "2010", "2015"))
output$mltree_raw_plot <- mltree_raw_plot
ggsave("../../figures/mltree_raw_plot.pdf", mltree_raw_plot,
       width=10, height=25)


# Plot tree with just Asian genotype --------------------------------------
mltree_asian <- drop.tip(mltree_raw, sapply(outgroup_names, grep, mltree_raw$tip.label))
mltree_asian$tip.label %<>% gsub("Borneo", "Indonesia", .) %>%
  gsub("USA", "UnitedStatesOfAmerica", .)
output$mltree_asian <- mltree_asian

mltree_asian_plot <- 
  ggtree(mltree_asian) +
  theme_tree2() +
  xlim(0, 0.14)
mltree_asian_df <- data.frame(seq=filter(mltree_asian_plot$data, isTip)$label,
                              stringsAsFactors = FALSE) %>%
  mutate(chrf=factor(grepl("CHRF", seq))) %>%
  mutate(date=as.Date(strsplit(seq, "_")%>%sapply(tail, 1))) %>%
  mutate(year=factor(format(date, "%Y"))) %>%
  mutate(countries=strsplit(seq, "_")%>%sapply(`[`, 2)%>%
           gsub("(and)([A-Z])", " \\1\\2", .)%>%gsub("([a-z])([A-Z])", "\\1 \\2", .)%>%
           gsub("ofthe", " of the", .))
mltree_asian_df %<>% mutate(., continent=factor(countrycode(sourcevar=countries, origin = "country.name", destination = "continent")))
mltree_asian_df[mltree_asian_df$countries=="Micronesia", "continent"] <- "Oceania"  
mltree_asian_df[mltree_asian_df$countries%in%c("Virgin Islands", "Saint Martin"), "continent"] <- "Americas"
rownames(mltree_asian_df) <- mltree_asian_df$seq
output$mltree_asian_df <- mltree_asian_df
mltree_asian_plot <- mltree_asian_plot %<+%
  mltree_asian_df +
  geom_tiplab(aes(colour=chrf), align=TRUE)
mltree_asian_plot %<>%
  gheatmap(data.frame(Continent=mltree_asian_df[, "continent"], row.names=rownames(mltree_asian_df)), 
           offset=.075, width=0.1, colnames=FALSE) %>%
  scale_x_ggtree()
mltree_asian_plot <- mltree_asian_plot + 
  guides(colour = guide_legend(title="CHRF", reverse=TRUE),
         fill = guide_legend(title="Continent"))
output$mltree_asian_plot <- mltree_asian_plot <-
  mltree_asian_plot + 
  theme(legend.title=element_text(size=12)) +
  scale_fill_discrete("Continent") +
  scale_color_manual("CHRF", values=c("gray50", "orangered1"), labels=c("Others", "CHRF isolates"))
ggsave("../../figures/mltree_asian_plot.pdf", mltree_asian_plot,
       width=10, height=30)


# South Asian genotype ----------------------------------------------------

mltree_asian_plot + geom_nodelab(aes(label=node))

output$mltree_southasian <- mltree_southasian <-
  extract.clade(mltree_asian, 244)

mltree_southasian_root2tip <-
  list(df=EpiGenR::root2tip.divergence(mltree_southasian))
mltree_southasian_root2tip$lm <- with(mltree_southasian_root2tip$df, lm(divergence~decimal_date(time)))
mltree_southasian_root2tip$subst.rate <-
  mltree_southasian_root2tip$lm$coeff[[2]]
mltree_southasian_root2tip$R2 <- 
  round(summary(mltree_southasian_root2tip$lm)$adj.r.squared, 2)
output$mltree_southasian_root2tip <- mltree_southasian_root2tip

output$mltree_southasian_root2tip_plot <-
  mltree_southasian_root2tip_plot <-
  ggplot(mltree_southasian_root2tip$df) + theme_bw() + 
    geom_point(aes(x=time, y=divergence)) + 
    stat_smooth(method="lm",aes(x=time, y=divergence)) +
    xlab("Date") + ylab("Divergence") +
    theme(text=element_text(size=12)) +
    annotate(geom="text",
             x=with(mltree_southasian_root2tip$df, min(time)+diff(range(time))*0.3),
             y=max(mltree_southasian_root2tip$df$divergence)*0.9, 
             label=bquote(atop(.("Substitution rate: ")~.(round(mltree_southasian_root2tip$subst.rate*10^4, 2)) %*% 10^-4~.("/site/year"),
                               italic(R)^2 == .(mltree_southasian_root2tip$R2)))%>%
               as.expression%>%as.character,
             parse=TRUE)
ggsave("../../figures/mltree_southasian_root2tip_plot.pdf", mltree_southasian_root2tip_plot,
       width=5.5, height=4)

mltree_southasian_plot <- ggtree(mltree_southasian) + theme_tree2() + xlim(0, 0.05)
mltree_southasian_df <- data.frame(seq=filter(mltree_southasian_plot$data, isTip)$label,
                              stringsAsFactors = FALSE) %>%
  left_join(., mltree_asian_df, by="seq")
rownames(mltree_southasian_df) <- mltree_southasian_df$seq
output$mltree_southasian_df <- mltree_southasian_df
mltree_southasian_plot <- mltree_southasian_plot %<+%
  mltree_southasian_df +
  geom_tiplab(aes(colour=chrf), align=TRUE) +
  scale_color_manual(values=c("gray50", "orangered1"), labels=c("Others", "CHRF isolates")) + 
  guides(colour = guide_legend(reverse=TRUE))
mltree_southasian_plot %<>%
  gheatmap(data.frame(Continent=mltree_southasian_df[, "continent"], row.names=rownames(mltree_southasian_df)), 
           offset=.01, width=0.1, colnames=FALSE) %>%
  scale_x_ggtree()
output$mltree_southasian_plot <- mltree_southasian_plot <-
  mltree_southasian_plot +
  theme(legend.title=element_text(size=12)) +
  scale_fill_discrete("Continent") +
  scale_color_manual("CHRF", values=c("gray50", "orangered1"), labels=c("Others", "CHRF isolates"))
ggsave("../../figures/mltree_southasian_plot.pdf", mltree_southasian_plot,
       width=10, height=30)


# Extract sequences -------------------------------------------------------
sequences <- 
  seqinr::read.alignment("../../data/ncbi_chrf_aln.fasta", format="fasta") %>%
  as.matrix()
sequences <-
  sequences[strsplit(rownames(sequences), "_") %>% 
              sapply(head, 1) %>% 
              lapply(grepl, mltree_southasian$tip.label) %>%
              sapply(any) %>%
              which() %>%
              c(., grep("CHRF", rownames(sequences))) %>%
              unique(), ]
output$sequences <- sequences <- sequences[sequences %>%
  `==`("n") %>%
  rowSums() %>%
  `<`(3000), ]
seqinr::write.fasta(apply(sequences, 1, paste, collapse="")%>%as.list, names=rownames(sequences), 
                    file.out="../../data/ncbi_chrf_aln_south_asian.fasta", as.string=TRUE)
output$partition_text <- partition_text <- readLines("../../data/partition_raxml.txt")
output$partition <- partition <- 
  partition_text %>%
  strsplit(., "= ") %>%
  sapply(tail, 1) %>%
  gsub("-", ":", .) %>%
  paste0("c(", ., ")") %>%
  parse(text=.) %>%
  lapply(eval)
output$partitioned_sequences <- 
  partitioned_sequences <-
  lapply(partition, function (x) sequences[, x])
mapply(function (fn, seqs) write.fasta(sequences=seqs, names=names(seqs), file.out=fn, as.string=TRUE),
       paste0("../../data/ncbi_chrf_aln_south_asian_g", 
              regmatches(partition_text, regexpr("p[0-9]", partition_text)), ".fasta") %>% gsub("gp3", "intergenic", .),
       lapply(partitioned_sequences, apply, 1, paste, collapse="") %>% lapply(as.list))

# save output -------------------------------------------------------------

save(list=names(output), envir=output, file="parse_raxml_output.RData")
