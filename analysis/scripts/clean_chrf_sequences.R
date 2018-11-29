clean_chrf_sequences <- function (seqs) {
  # Change sequence name format if not already in CHRFXXXX_YYYY-MM-DD format
  if (!grepl("Bangladesh", names(seqs)[1])) {
    split_str <- strsplit(names(seqs), "_")
    newnames <- sapply(split_str, function (x) {
      paste0(x[1], x[2], "_Bangladesh_", as.Date(x[3], "%m-%d-%Y"))
    })
    names(seqs) <- newnames
  }
  return (seqs)
}