# The dada2 processing is based on the offical DADA2 pipeline tutorial 1.16 
#accessible under the following linkhttps://benjjneb.github.io/dada2/tutorial.html


library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")


path1 <- "path to raw reads"
list.files(path1)
path2 <- "path to were processed files should be stored"
fnFs <- sort(list.files(path1, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(path1, pattern="_R2_001.fastq"))
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
fnFs <- file.path(path1, fnFs)
fnRs <- file.path(path1, fnRs)
filt_path <- file.path(path2, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(266,127),
                     maxN=0, maxEE=c(3,6), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, verbose=TRUE) # On Windows set multithread=FALSE
head(out)

# Learn error rates
errF <- learnErrors(filtFs, nbases = 1e+08, randomize = TRUE, multithread=TRUE, verbose = 1)

# Learn error rates
errR <- learnErrors(filtRs, nbases = 1e+08, randomize = TRUE, multithread=TRUE, verbose = 1)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, verbose=1)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE, verbose=1)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 12,
                      verbose=TRUE)

head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

tax <- assignTaxonomy(seqtab.nochim, "path to training fasta file", multithread=TRUE, verbose = TRUE)
#we used the SILVA training set v. 138