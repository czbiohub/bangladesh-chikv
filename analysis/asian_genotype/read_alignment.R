library(seqinr)
library(magrittr)
library(dplyr)
library(lubridate)
library(reutils)
library(parallel)
library(ggplot2)

sequences <- read.fasta("../../data/CHKV_rename-999-filtered-length.fasta")
acc_num <- strsplit(names(sequences), "|", fixed=TRUE) %>%
  sapply(`[`, 4)
search_results <- mclapply(acc_num, efetch, db="nuccore", rettype="gb", retmode="text", mc.cores=4) %>%
  sapply(content) %>%
  strsplit("\n")
collection_dates <- search_results %>%
  mapply(grep, "collection_date", ., value=TRUE, SIMPLIFY=FALSE) %>%
  lapply(., strsplit, '"') %>%
  lapply(unlist) %>%
  lapply(`[`, 2) %>%
  sapply(function (x) {
    if (length(x)==0) return (NA)
    x
  })
countries <- search_results %>%
  mapply(grep, "country", ., value=TRUE, SIMPLIFY=FALSE) %>%
  lapply(., strsplit, '"') %>%
  lapply(unlist) %>%
  lapply(`[`, 2) %>%
  sapply(function (x) {
    if (length(x)==0) return (NA)
    x
  })
nodates_pos <- which(sapply(collection_dates, is.na))
pos <- search_results[nodates_pos] %>%
  mapply(grepl, "Genome-scale phylogenetic analyses of chikungunya virus reveal", ., SIMPLIFY=FALSE) %>%
  sapply(any) %>%
  unlist %>% 
  which %>%
  `[`(nodates_pos, .)
collection_dates[pos] <- c("before 1981", "1995", "before 1986")
nodates_pos <- which(sapply(collection_dates, is.na))
collection_dates[nodates_pos[match(c("JQ067624.1", "EF027136.1", "DQ443544.2"), acc_num[nodates_pos])]] <- c("2010", "2006", "2006")
collection_dates[grep("KT449801", acc_num)] <- "2006"
collection_dates[grep("NC_004162", acc_num)] <- "1953" #NC_004162 == AF369024, Tazania, 1953
countries[grep("NC_004162", acc_num)] <- "Tazania"
acc_num[grep("NC_004162", acc_num)] <- "AF369024"
collection_dates[grep("AB455494", acc_num)] <- "2006"
countries[grep("AB455494", acc_num)] <- "Sri Lanka"
collection_dates[grep("JQ067624", acc_num)] <- "2010"
countries[grep("JQ067624", acc_num)] <- "China"

collection_dates_formatted <- lapply(collection_dates, function (x) {
  if (is.na(x)) return (x)
  if (grepl("efore", x)) return (x)
  if (nchar(x)==4) return (as.Date(paste0(x, "-07-01")))
  if (grepl("^[A-Z]", x)) return (as.Date(paste0("15-", x), "%d-%b-%Y"))
  if (grepl("^[0-9][0-9]-[a-zA-Z]", x)) return (as.Date(x, "%d-%b-%Y"))
  if (grepl("^[0-9][0-9][0-9][0-9]-[0-9][0-9]$", x)) return (as.Date(paste0(x, "-15")))
  if (is.Date(try(as.Date(x)))) return (as.Date(x))
  return (NA)
}) %>%
  do.call(what=c)

remove_ind <- unique(c(which(sapply(countries, is.na)), which(is.na(collection_dates_formatted)), grep("efore", collection_dates_formatted)))

ncbi_newnames <-  paste(acc_num[-remove_ind],
                        strsplit(countries[-remove_ind], ":")%>%sapply(`[`, 1), 
                        collection_dates_formatted[-remove_ind]%>%as.character(), sep="_")
ncbi_sequences <- sequences[-remove_ind]
names(ncbi_sequences) <- ncbi_newnames
write.fasta(ncbi_sequences, 
            ncbi_newnames,
            "../../data/ncbi_seq_cleaned.fasta")

# CHRF samples ------------------------------------------------------------

chrf_samples <- read.fasta("../../data/CHKV-annotated-Bang2017-trimmed-edited.fasta") %>%
  lapply(., function (x) {
    as.SeqFastadna(x[x!="-"], name=attr(x, "name"), Annot=attr(x, "Annot"))
  })
chrf_samples_metadata <- read.csv("../../data/project-rapid_response_007_sample-table.csv", stringsAsFactors=FALSE) %>%
  slice(., sapply(names(chrf_samples), grep, sample_name))
chrf_samples_dates <- strsplit(chrf_samples_metadata$notes, "\n") %>% 
  sapply(function (x) grep("Admission date", x, value=TRUE)) %>% 
  gsub("- Admission date: ", "", .) %>%
  as.Date(., "%d-%b-%y")
chrf_samples_newnames <- paste0(names(chrf_samples), "_Bangladesh_", chrf_samples_dates)
write.fasta(chrf_samples, 
            chrf_samples_newnames,
            "../../data/chrf_samples_cleaned.fasta")


# cleaned_data_meta -------------------------------------------------------

cleaned_data_metadata <- data.frame(id=ncbi_newnames,
           country=strsplit(ncbi_newnames, "_") %>% lapply(tail, 2) %>% sapply(head, 1),
           collection_date=as.Date(strsplit(ncbi_newnames, "_") %>% sapply(tail, 1))) %>%
  rbind(., data.frame(id=chrf_samples_newnames, country="Bangladesh", collection_date=chrf_samples_dates))

write.table(cleaned_data_metadata, "../../data/cleaned_metadata.txt",
            sep="\t", quote=FALSE, row.names=FALSE)


freq_table <- cleaned_data_metadata %>%
  mutate(collection_date=format(collection_date, "%Y")) %>%
  select(-1) %>%
  table() 
freq_table %>%
  reshape2::melt() %>%
  ggplot() %>%
  `+`(geom_tile(aes(x=collection_date, y=country, fill=value)))

country_year <- mapply(substr, ncbi_newnames, regexpr("_", ncbi_newnames)+1, nchar(ncbi_newnames)-6)
ncbi_subsampled <- split(ncbi_sequences, country_year) %>%
  lapply(`[[`, 1)
names(ncbi_subsampled) <- split(names(ncbi_sequences), country_year) %>% sapply(`[`, 1) %>% unname()
write.fasta(ncbi_subsampled, names(ncbi_subsampled), "../../data/ncbi_seq_cleaned_sampled.fasta")




