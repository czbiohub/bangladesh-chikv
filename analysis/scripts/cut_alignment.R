library(seqinr)#, lib.loc="/mnt/data/")
fastafile <- "../data/alignment/ncbi_chrf_aln.fasta"
aln <- as.matrix(read.alignment(fastafile, format="fasta"))
firstnames <- sapply(strsplit(rownames(aln), "_"), `[`, 1)
firstname <- grep("CHRF", firstnames, invert=TRUE, value=TRUE)[1]
seqchar <- as.character(aln[grep(firstname, rownames(aln)), ])
tmp <- tempfile()
system(paste0('curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=', firstname, '&rettype=gb&retmode=text" > ', tmp))
tmplines <- readLines(tmp)
ranges <- lapply(strsplit(gsub("[ ]+CDS[ ]+", "", grep("CDS", tmplines, value=TRUE)), "..", fixed=TRUE), as.numeric)
CHIKVgp1_pos <- which(seqchar!="-")[ranges[[1]]]
CHIKVgp2_pos <- which(seqchar!="-")[ranges[[2]]]
aln_gp1 <- aln[, CHIKVgp1_pos]
aln_gp2 <- aln[, CHIKVgp2_pos]
# seqnames <- gsub(">", "", grep(">", readLines(fastafile), value=TRUE))
# lapply(1:2, function (i) {
#   aln_gpx <- aln_gp1
#   if (i==2) aln_gpx <- aln_gp2
#   write.fasta(mapply(as.SeqFastadna, apply(aln_gpx, 1, paste, collapse=""), 
#                      name=seqnames, Annot=seqnames, SIMPLIFY=FALSE),
#               seqnames, gsub("aln", paste0("aln_gp", i), fastafile))
# })
#save.image("/mnt/data/data/cut_alignment.RData")

cat("DNA, p1 = ", paste(range(CHIKVgp1_pos), collapse="-"), "\n",
    "DNA, p2 = ", paste(range(CHIKVgp2_pos), collapse="-"), "\n",
    "DNA, p3 = ", paste(c(1, min(CHIKVgp1_pos)-1), collapse="-"), ", ", 
    paste(c(max(CHIKVgp1_pos)+1, min(CHIKVgp2_pos)-1), collapse="-"), ", ",
    paste(c(max(CHIKVgp2_pos)+1, ncol(aln)), collapse="-"),
    sep="", file="../data/modeltest/partition_modeltest.txt")
