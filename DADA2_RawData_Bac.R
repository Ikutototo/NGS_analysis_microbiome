# Setting -----------------------------------
setwd("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome")
dir <- getwd()

# .fastq.gz fileを格納するdirectory
# fastq.gz fileのあるフォルダーを指定する → その都度変更する
path <- paste(dir, "RawData_Bacteria", sep = "/")
list.files(path = path, full.names = TRUE)

lib <- c("dada2", "ggplot2", "dplyr", "Biostrings", "phyloseq", "ShortRead", "crayon")
for (i in lib) {
    print(i)
    print(packageVersion(i))
    library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}

# PlotQualityProfile[BEFORE] ------------------------

fnFs <- sort(list.files(path, pattern = "_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


PlotQC_BF_RData_For <- plotQualityProfile(fnFs[1:9])
PlotQC_BF_RData_Rev <- plotQualityProfile(fnRs[1:9])

# .FASTQ → ShortRead::readFastq() → import .fastq
readfnRs_lengths_summary <- lapply(fnRs, function(file) {
    reads <- ShortRead::readFastq(file)
    read_lengths <- width(ShortRead::sread(reads))
    summary(read_lengths)
})
names(readfnRs_lengths_summary) <- basename(fnRs)
print(head(readfnRs_lengths_summary, 10))


readfnFs_lengths_summary <- lapply(fnFs, function(file) {
    reads <- ShortRead::readFastq(file)
    read_lengths <- width(ShortRead::sread(reads))
    summary(read_lengths)
})
names(readfnFs_lengths_summary) <- basename(fnFs)
print(head(readfnFs_lengths_summary, 10))


save(PlotQC_BF_RData_For, PlotQC_BF_RData_Rev,
    readfnRs_lengths_summary, readfnFs_lengths_summary,
    file = "./RData/PlotQC/PlotQC_BEFORE_Data.RData"
)

# PlotQualityProfile[AFTER] -----------------
# parameterによって結果が異なるため、注意
# オブジェクトとして保存すること( <- )
Results_FilterTrim <- filterAndTrim(
    fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs,
    truncLen = c(219, 219), maxN = 0, maxEE = c(2, 2), truncQ = 2,
    rm.phix = TRUE, compress = TRUE, multithread = TRUE
)

PlotQC_AF_RData_For <- plotQualityProfile(filtFs[1:9])
cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))

PlotQC_AF_RData_Rev <- plotQualityProfile(filtRs[1:9])
cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))

readfiltRs_lengths_summary <- lapply(filtRs, function(file) {
    reads <- ShortRead::readFastq(file)
    read_lengths <- width(ShortRead::sread(reads))
    summary(read_lengths)
})
names(readfiltRs_lengths_summary) <- basename(filtRs)
print(head(readfiltRs_lengths_summary, 10))


readfiltFs_lengths_summary <- lapply(filtFs, function(file) {
    reads <- ShortRead::readFastq(file)
    read_lengths <- width(ShortRead::sread(reads))
    summary(read_lengths)
})
names(readfiltFs_lengths_summary) <- basename(filtFs)
print(head(readfiltFs_lengths_summary, 10))

cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))

save(PlotQC_AF_RData_For, PlotQC_AF_RData_Rev,
    readfiltFs_lengths_summary, readfiltRs_lengths_summary,
    file = "./RData/PlotQC/PlotQC_Filt_Data.RData"
)


# Learn the Error Rates ---------------------
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
dadaFs[[1]]
dadaRs[[1]]


save(errF, errR, dadaFs, dadaRs,
    file = "./RData/PlotQC/dada2_ErrorRateData.RData"
)

cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))



# Merge paired reads ------------------------


mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[, nchar(colnames(seqtab)) %in% 250:257]

dim(seqtab2)
table(nchar(getSequences(seqtab2)))

seqtab.nochim <- removeBimeraDenovo(seqtab2, method = "consensus", multithread = TRUE, verbose = TRUE)

dim(seqtab.nochim)
dim(seqtab.nochim)[2] / dim(seqtab)[2]
sum(seqtab.nochim) / sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(Results_FilterTrim, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track


save(mergers, seqtab, seqtab2, seqtab.nochim, getN, track,
    file = "./RData/PlotQC/MergeData.RData"
)

# Assign Taxonomy ---------------------------

# SILVA database (v.138)
system.time(taxa <- assignTaxonomy(seqtab.nochim,
    refFasta = "~/Documents/RStudio/Novogene/250503/taxa_reference/silva_nr99_v138.2_toGenus_trainset.fa.gz",
    multithread = TRUE
))

system.time(taxa_species <- assignTaxonomy(seqtab.nochim,
    refFasta = "~/Documents/RStudio/Novogene/250503/taxa_reference/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
    multithread = TRUE
))

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print, 20)


save(taxa, taxa_species,
    file = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/Input/phyloseq_Bacteria/TaxaData.RData"
)

cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))


# Inspects SilvaDatabase --------------------
library(DECIPHER)

# SilvaDatabaseのPATHを指定する
Bacteria_dna <- readDNAStringSet("~/Documents/RStudio/Novogene/250503/taxa_reference/silva_nr99_v138.2_toGenus_trainset.fa.gz")
str(Bacteria_dna)
names(Bacteria_dna)[1:3]
Bacteria_dna@ranges@NAMES
