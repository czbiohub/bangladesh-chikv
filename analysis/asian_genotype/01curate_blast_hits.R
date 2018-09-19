library(seqinr)
library(dplyr)
library(reutils)
library(parallel)
library(magrittr)

output.env <- new.env()

qlen <- read.fasta("../../data/CHKV-annotated-Bang2017-trimmed-edited.fasta") %>% sapply(length)
blast_hits <- read.csv("../../data/blast-hit-table.csv", header=FALSE, stringsAsFactors = FALSE,
                       col.names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>%
  mutate(qcov=abs(qend-qstart)/qlen[qseqid]) %>%
  dplyr::filter(!duplicated(sseqid))
ggplot(blast_hits) + theme_bw() + geom_histogram(aes(x=qcov*pident))
output.env$blast_hits <- blast_hits %<>% dplyr::filter(qcov*pident>75)

output.env$qlen <- qlen

genbank_results <- mclapply(blast_hits$sseqid, reutils::efetch, db="nuccore", rettype="gb", retmode="text", mc.cores=4) %>%
  lapply(content)

output.env$genbank_results <- genbank_results

output.env$genbank_results_filtered <-
  genbank_results_filtered <- 
  mclapply(genbank_results, function (result) {
  result <- strsplit(result, "\n")[[1]]
  country <- grep("country", result, value=TRUE)
  if (length(country)==0) return (NULL)
  if (grepl("unknown", country)) return (NULL)
  country <- tail(strsplit(country, '"')[[1]], 1)
  coll_date <- grep("collection_date", result, value=TRUE)
  if (length(coll_date)==0) return (NULL)
  coll_date <- tail(strsplit(coll_date, '"')[[1]], 1)
  coll_date_formatted <- lapply(list(coll_date), function (x) {
    if (is.na(x)) return (x)
    if (grepl("efore", x)) return (x)
    if (grepl("^[0-9][0-9][0-9][0-9]$", x)) return (as.Date(paste0(x, "-07-02")))
    if (grepl("^[a-zA-Z]", x)) return (as.Date(paste0("15-", x), "%d-%b-%Y"))
    if (grepl("^[0-9][0-9]-[a-zA-Z]", x)) return (as.Date(x, "%d-%b-%Y"))
    if (grepl("^[0-9][0-9][0-9][0-9]-[0-9][0-9]$", x)) return (as.Date(paste0(x, "-15")))
    if (class(try(as.Date(x)))=="Date") return (as.Date(x))
    return (NA)
  }) %>%
    do.call(what=c)
  return (data.frame(id=grep("ACCESSION", result, value=TRUE) %>% gsub("ACCESSION   ", "", .),
                     location=country, country=strsplit(country, ":")[[1]][1],
                     collection_date_raw=coll_date,
                     collection_date=coll_date_formatted,
                     stringsAsFactors=FALSE))
}, mc.cores=4) %>%
  do.call(what=rbind) %>%
  mutate(., pident=with(blast_hits[sapply(.$id, grep, blast_hits$sseqid), ], qcov*pident))

# Hall et al (2016): We also confirm that the inclusion of many recent sequences from a single geographical location in an analysis 
# tends to result in a spurious bottleneck effect in the reconstruction and caution against interpreting this as genuine.

set.seed(801820)
output.env$genbank_results_filtered_sampled <-
  genbank_results_filtered_sampled <-
  select(genbank_results_filtered, c(country, collection_date)) %>%
  mutate(collection_date=format(collection_date, "%Y")) %>%
  mutate(x=paste0(country, collection_date)) %>%
  `[[`("x") %>%
  split(genbank_results_filtered, .) %>%
  lapply(function (x) {
    if (nrow(x)<=1) return (x)
    x[sample(1:nrow(x), 2), ]
  }) %>%
  do.call(what=rbind)

ggplot(genbank_results_filtered_sampled) + theme_bw() + geom_histogram(aes(x=pident), color="white") + geom_vline(xintercept=85, color="red")
# use the sequences with less than 85% identity as outgroup

write.table(genbank_results_filtered_sampled, "../../data/blast_filtered_sampled_metadata.txt",
            sep="\t", row.names=FALSE, quote=FALSE)

blast_filtered_sampled_fasta_fn <- "../../data/blast_filtered_sampled.fasta"
efetch(genbank_results_filtered_sampled$id, db="nuccore", rettype="fasta", retmode="text") %>%
  content() %>%
  cat(., file=blast_filtered_sampled_fasta_fn)
genbank_sequences <- read.fasta(blast_filtered_sampled_fasta_fn)
names(genbank_sequences) <- gsub(" ", "", with(genbank_results_filtered_sampled, paste0(id, "_", country, "_", collection_date)))
output.env$genbank_sequences <- genbank_sequences
write.fasta(genbank_sequences, names(genbank_sequences),
            file.out=blast_filtered_sampled_fasta_fn)

writeLines(gsub(" ", "", names(genbank_sequences)), "../../data/list_of_blast_hit_names.txt")



output.env$chrf_seq <- chrf_seq <- read.fasta("../../data/Final-CHIKV-sequences.fasta")

write.fasta(chrf_seq, 
            strsplit(names(chrf_seq), "_") %>% sapply(., function (x) paste0(x[1], x[2], "_Bangladesh_", as.Date(x[3], "%m-%d-%Y"))),
            "../../data/chrf_samples_cleaned.fasta")

fn <- c("../../data/blast_filtered_sampled.fasta", "../../data/chrf_samples_cleaned.fasta", "../../data/inseq.fasta",
        "cut_alignment.R")
system(paste("cat", fn[1], fn[2], ">", fn[3]))
paste("~/anaconda3/bin/aws s3 cp", fn, "s3://lucymli/bangladesh_chikv/") %>%
  sapply(., system)

save(list=names(output.env), envir=output.env, file="curate_blast_hits.RData")

#muscle -in inseq.fasta -out ncbi_chrf_aln.fasta
# lapply(list.files("../../data", ".fasta$", full.names = TRUE), function (x) {
#   seqs <- read.fasta(x)
#   if (!any(names(seqs) %in% names(genbank_sequences))) return (NULL)
#   cat (x, "\n")
#   pos <- sapply(names(seqs), grep, names(genbank_sequences))
#   names(seqs)[sapply(pos, length)>0] <- names(genbank_sequences)[unlist(pos)]
#   write.fasta(seqs, names(seqs), x)
# })
# lapply(list.files("../../data", "[tT]ree$", full.names = TRUE), function (x) {
#   tree <- ape::read.tree(x)
#   if (!any(tree$tip.label %in% names(genbank_sequences))) return (NULL)
#   cat (x, "\n")
#   pos <- sapply(tree$tip.label, grep, names(genbank_sequences))
#   tree$tip.label[sapply(pos, length)>0] <- names(genbank_sequences)[unlist(pos)]
#   ape::write.tree(tree, x)
# })
