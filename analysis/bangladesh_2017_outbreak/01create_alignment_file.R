library(magrittr)
library(seqinr)
library(reutils)
library(parallel)
library(dplyr)

# get_e1_seq <- function (acc_number, name=NULL, Annot=NULL, return_index=FALSE) {
#   tmpf <- tempfile()
#   reutils::efetch(acc_number, db="nucleotide", rettype="gb", outfile=tmpf)
#   gbfile <- readLines(tmpf)
#   full_ranges <- strsplit(grep("CDS", gbfile, value=TRUE), "..", fixed=TRUE) %>%
#     lapply(function (x) as.numeric(gsub("[^0-9]", "", x)))
#   full_ranges_e1 <- full_ranges[[which.min(sapply(full_ranges, diff))]]
#   SEQ <- gsub("[^a-z]", "", paste(gbfile[-1:-grep("ORIGIN", gbfile)], collapse="")) %>%
#     substr(., full_ranges_e1[1], full_ranges_e1[2]) %>%
#     strsplit(., "") %>%
#     `[[`(1) %>%
#     seqinr::as.SeqFastadna(., name, Annot)
#   if (return_index) return (list(seq=SEQ, index=full_ranges_e1))
#   return (SEQ)
# }

seqs <- read.fasta("../../data/ncbi_chrf_aln.fasta")
outbreak_seqs <- seqs[c(grep("Bangladesh_2017", names(seqs)), grep("Italy_2017", names(seqs)))]
tmpf <- tempfile()


# non_chrf_seqnames <- grep("CHRF", names(outbreak_seqs), invert=TRUE, value=TRUE)
# search_acc <- sapply(strsplit(non_chrf_seqnames, "_"), `[`, 1)
# subseq <- get_e1_seq(search_acc[1], name=non_chrf_seqnames[1], Annot=attr(outbreak_seqs[[non_chrf_seqnames[1]]], "Annot"), return_index=TRUE)
# outbreak_seqs <- lapply(outbreak_seqs, function (x) {
#   seqinr::as.SeqFastadna(x[subseq$index[1]:subseq$index[2]], name=names(x), Annot=attr(x, "Annot"))
# })


e1seqs_acc <- paste0("MG", 697262:697282)

e1_blast_hits <- read.table("../../data/outbreak_blast_hits.txt", comment.char ="#", sep="\t", 
                            col.names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"),
                            stringsAsFactors = FALSE) %>%
  dplyr::filter(., !(sseqid %in% e1seqs_acc)) %>%
  mutate(qcov=abs(qend-qstart)/length(e1seqs[[1]])) %>%
  filter(qcov*pident > 99.6)

efetch(c(e1seqs_acc, e1_blast_hits$sseqid), db="nuccore", rettype="fasta", retmode="text", outfile=tmpf)
e1seqs <- read.fasta(tmpf)
e1info <- do.call(rbind, mclapply(names(e1seqs), function (x) {
  outname <- tempfile()
  efetch(x, db="nuccore", rettype="gb", retmode="text", outfile=outname)
  date_output <- system(paste("grep -n collection_date ", outname), intern = TRUE) %>%
    stringi::stri_extract_all_regex(., '(?<=").*?(?=")') %>%
    unlist()
  if (nchar(date_output)==4) date_output <- as.Date(paste0(date_output, "-07-02"))
  else date_output <- as.Date(date_output, "%d-%b-%Y")
  if (grepl("MG049915", x)) date_output <- as.Date(c("2017-09-15"))
  country_output <- system(paste("grep -n country ", outname), intern=TRUE) %>%
    stringi::stri_extract_all_regex(., '(?<=").*?(?=")') %>%
    unlist()
  if (grepl(":", country_output)) country_output <- strsplit(country_output, ":", fixed=TRUE)[[1]][[1]]
  data.frame(date=date_output, country=country_output)
}, mc.cores=8))
names(e1seqs) <- paste(gsub(".1", "", names(e1seqs), fixed=TRUE), e1info$country, e1info$date, sep="_") %>%
  gsub(" ", "", .)

combined_outbreak_seqs <- c(e1seqs, outbreak_seqs)
names(combined_outbreak_seqs) <- c(names(e1seqs), names(outbreak_seqs))
combined_outbreak_seqs <- combined_outbreak_seqs[!duplicated(names(combined_outbreak_seqs))]


# which(sapply(combined_outbreak_seqs, length)>5000)


write.fasta(combined_outbreak_seqs, names(combined_outbreak_seqs), "../../data/inseq_outbreak.fasta")
combined_outbreak_seqs[sapply(combined_outbreak_seqs, length)>5000] %>%
  write.fasta(., names(.), "../../data/inseq_outbreak_long.fasta")
combined_outbreak_seqs[sapply(combined_outbreak_seqs, length)<5000] %>%
  write.fasta(., names(.), "../../data/inseq_outbreak_short.fasta")

system("/usr/local/bin/muscle3.8.31_i86darwin64 -in ../../data/inseq_outbreak_long.fasta -out ../../data/outbreak_aln_long.fasta")
system("/usr/local/bin/muscle3.8.31_i86darwin64 -in ../../data/inseq_outbreak_short.fasta -out ../../data/outbreak_aln_short.fasta")
system("/usr/local/bin/muscle3.8.31_i86darwin64 -profile -in1 ../../data/outbreak_aln_long.fasta -in2 ../../data/outbreak_aln_short.fasta -out ../../data/outbreak_aln.fasta")

outbreak_aln <- read.fasta("../../data/outbreak_aln.fasta")

aln_ranges <- do.call(rbind, lapply(outbreak_aln, function (x) range(which(x!="-"))))
e1_range <- aln_ranges[which.min(apply(aln_ranges, 1, diff)), ]
outbreak_e1_aln <- lapply(outbreak_aln, function (x) seqinr::as.SeqFastadna(x[c(e1_range[1]:e1_range[2])], name=names(x), Annot=attr(x, "Annot")))
write.fasta(outbreak_e1_aln, names(outbreak_e1_aln), "../../data/outbreak_e1_aln.fasta")

outbreak_e1_aln

