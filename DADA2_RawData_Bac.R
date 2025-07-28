# Setting -----------------------------------
setwd("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome")

ls()
rm(list = ls())
dev.off()

dir <- "/Users/ikutosaito/Documents/RStudio/Novogene/250503/Novorgene"

# .fastq.gz fileを格納するdirectory
# fastq.gz fileのあるフォルダーを指定する → その都度変更する
path <- paste(dir, "RawData_Bacteria", sep = "/")
list.files(path = path, full.names = TRUE)


# lib <- c("dada2", "ggplot2", "dplyr", "Biostrings", "phyloseq", "ShortRead", "crayon")
# for (i in lib) {
#     print(i)
#     print(packageVersion(i))
#     library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
# }

# PlotQualityProfile[BEFORE] ------------------------

library(phyloseq)
library(dada2)
library(dplyr)
library(ShortRead)
library(Biostrings)

fnFs <- sort(list.files(path, pattern = "_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
# filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
filtFs <- file.path(path, "250728_filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "250728_filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


PlotQC_BF_RData_For <- plotQualityProfile(fnFs[1:9])
PlotQC_BF_RData_Rev <- plotQualityProfile(fnRs[1:9])


# Reads_Lengths_Summary[BEFORE] ------------------
# .FASTQ → ShortRead::readFastq() → import .fastq

# Forward → Meanが最も小さいSampleを確認する
readfnFs_lengths_summary <- lapply(fnFs, function(file) {
    reads <- ShortRead::readFastq(file)
    read_lengths <- width(ShortRead::sread(reads))
    summary(read_lengths)
})
names(readfnFs_lengths_summary) <- basename(fnFs)
print(head(readfnFs_lengths_summary, 10))

# Reverse → Meanが最も小さいSampleを確認する
readfnRs_lengths_summary <- lapply(fnRs, function(file) {
    reads <- ShortRead::readFastq(file)
    read_lengths <- width(ShortRead::sread(reads))
    summary(read_lengths)
})
names(readfnRs_lengths_summary) <- basename(fnRs)
print(head(readfnRs_lengths_summary, 10))

# RData Objects
save(PlotQC_BF_RData_For, PlotQC_BF_RData_Rev,
    readfnRs_lengths_summary, readfnFs_lengths_summary,
    file = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/PlotQC/PlotQC_BEFORE_Data.RData"
)


# 重要：下流のすべての解析に影響があるので、慎重に！！
# parameterによって下流の結果が異なるため、注意
# FastQC, SeqenceLengthDistribution → 参考に決める(truncLen())
# Bac5(Forward): Mean=223, Bac9(Reverse): Mean=221
# オブジェクトとして保存すること

# 250728_Filtering
# F/Rにおいて、品質の低下は見られないため、truncLen()は指定しない
Results_FilterTrim <- filterAndTrim(
    fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs,
    maxN = 0, maxEE = c(2, 2), truncQ = 2,
    rm.phix = TRUE, compress = TRUE, multithread = TRUE
)

# Results_FilterTrim <- filterAndTrim(
#     fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs,
#     truncLen = c(219, 219), maxN = 0, maxEE = c(2, 2), truncQ = 2,
#     rm.phix = TRUE, compress = TRUE, multithread = TRUE
# )


# PlotQualityProfile[AFTER] -----------------

PlotQC_AF_RData_For <- plotQualityProfile(filtFs[1:9])
cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))

PlotQC_AF_RData_Rev <- plotQualityProfile(filtRs[1:9])
cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))


# Reads_Lengths_Summary[AFTER] ------------------
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
    file = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/PlotQC/PlotQC_Filt_Data.RData"
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


# save(errF, errR, dadaFs, dadaRs,
#     file = "./RData/PlotQC/dada2_ErrorRateData.RData"
# )

save(errF, errR, dadaFs, dadaRs,
    file = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/PlotQC/250728_dada2_ErrorRateData.RData"
)

cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))


# Merge paired reads ------------------------

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab <- seqtab[, nchar(colnames(seqtab)) %in% 250:257]
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

dim(seqtab.nochim)
dim(seqtab.nochim)[2] / dim(seqtab)[2]
sum(seqtab.nochim) / sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(Results_FilterTrim, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "Nonchim")
rownames(track) <- sample.names
track


# save(mergers, seqtab, seqtab2, seqtab.nochim, getN, track,
#     file = "./RData/PlotQC/MergeData.RData"
# )

save(mergers, seqtab, seqtab.nochim, getN, track,
    file = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Input/250728_MergeData.RData"
)


# Assign Taxonomy ---------------------------

ls()
rm(list = ls())
dev.off()

# dada2::assignTaxonomy() → k-mer頻度とブートストラップによる類似性ベースの分類
load(file = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Input/250728_MergeData.RData")

# SILVA database (v.138)
# taxaオブジェクトのブートストラップ値表示
system.time(taxa <- assignTaxonomy(seqtab.nochim,
    refFasta = "~/Documents/RStudio/Novogene/250503/taxa_reference/silva_nr99_v138.2_toGenus_trainset.fa.gz",
    multithread = TRUE, outputBootstraps = TRUE, verbose = TRUE, minBoot = 50
))

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print, 20)

save(taxa,
    file = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Input/250728_TaxaData.RData"
)

cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))

write.csv(taxa$tax, file = "/Users/ikutosaito/Documents/RStudio/Novogene/250503/export_csv/250728_Taxa.csv")

## Kingom to Species -------------------------
# taxa <- assignTaxonomy(seqtab.nochim,
#                        refFasta = "~/Documents/RStudio/Novogene/Data/NGS_Analysis/sh_general_release_19.02.2025/sh_general_release_dynamic_19.02.2025.fasta",
#                        multithread=TRUE)

# assignTaxonomy()されたオブジェクトに、Speciesを割り当てる
system.time(taxa_species_plus <- dada2::addSpecies(
    taxtab = taxa, refFasta = "~/Documents/RStudio/Novogene/250503/taxa_reference/silva_v138.2_assignSpecies.fa.gz",
    verbose = TRUE, allowMultiple = TRUE,
))


## Species assignment ------------------------

# As of version 1.5.2(dada2), assignSpecies has been reimplemented to be much faster,
# and it now easily scales to handle even very large sets of sequences.

# By default the assignSpecies method only returns species assignments if there is no ambiguity,
# i.e. all exact matches were to the same species.
# However, given that we are generally working with fragments of the 16S gene,
# it is common that exact matches are made to multiple sequences that are identical over the sequenced region.
# This is often still useful information, so to have all sequence hits returned the returnMultiple=TRUE argument can be passed to the assignSpecies function:


# dada2::assignSpecies()：クエリ配列と参照配列との100%完全一致で、種の割り当てを行う

# taxa <- dada2::assignSpecies(
#     seqs = seqtab.nochim,
#     refFasta = "~/Documents/RStudio/Novogene/Data/NGS_Analysis/taxa_reference/silva_v138.2_assignSpecies.fa.gz", # RDPトレーニングセットやSILVA形式など
#     outputBootstraps = TRUE, # BootstrappingScore
#     multithread = TRUE,
#     allowMultiple = TRUE
# )

# If using this workflow on your own data:
# In many environments, few sequences will be assigned to species level.
# That is OK! Reference databases are incomplete,
# and species assignment is at the limit of what is possible from 16S amplicon data.
# Remember, even if the sequenced organism is in the reference database,
# if another species shares the same 16S gene sequence it is impossible to unambiguously assign that sequence to species-level.


# Inspects AmpliconDatabase --------------------
library(DECIPHER)
# General FASTA release (download) [https://unite.ut.ee/repository.php]

## SILVA_SSUfungi_nr99_v138_2_toGenus_trainset.fasta -------------

# SilvaDatabaseのPATHを指定する
Bacteria_dna <- readDNAStringSet("~/Documents/RStudio/Novogene/250503/taxa_reference/silva_nr99_v138.2_toGenus_trainset.fa.gz")
str(Bacteria_dna)


# 分類情報の抽出
# Bacteria_dna@ranges@NAMES
names(Bacteria_dna)[1:3]
as.character(Bacteria_dna)[1:3]


# ;でsplit
split_all <- strsplit(names(Bacteria_dna), "\\;")
split_all[1:3]


max_levels <- max(sapply(split_all, length))

# 欠損値をNAで埋める
split_padded <- lapply(split_all, function(x) {
    c(x, rep(NA, max_levels - length(x)))
})

df <- as.data.frame(do.call(rbind, split_padded), stringsAsFactors = FALSE)

colnames(df) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
head(df)

# 配列データ列を追加
df$Sequence <- as.character(Bacteria_dna)
df$Length <- nchar(df$Sequence)

str(df)


write.csv(
    x = df,
    file = "~/Documents/RStudio/Novogene/250503/export_csv/SilvaDatabase/silva_nr99_v138.2_toGenus_trainset.fa.gz.csv", row.names = TRUE
)

## silva_nr99_v138.2_toSpecies_trainset.fa.gz -------------

# SilvaDatabaseのPATHを指定する
Bacteria_dna <- readDNAStringSet("~/Documents/RStudio/Novogene/250503/taxa_reference/silva_nr99_v138.2_toSpecies_trainset.fa.gz")
str(Bacteria_dna)


# 分類情報の抽出
# Bacteria_dna@ranges@NAMES
names(Bacteria_dna)[1:3]
as.character(Bacteria_dna)[1:3]


# ;でsplit
split_all <- strsplit(names(Bacteria_dna), "\\;")
split_all[1:3]


max_levels <- max(sapply(split_all, length))

# 欠損値をNAで埋める
split_padded <- lapply(split_all, function(x) {
    c(x, rep(NA, max_levels - length(x)))
})

df <- as.data.frame(do.call(rbind, split_padded), stringsAsFactors = FALSE)

colnames(df) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
head(df)

# 配列データ列を追加
df$Sequence <- as.character(Bacteria_dna)
df$Length <- nchar(df$Sequence)

str(df)


write.csv(
    x = df,
    file = "~/Documents/RStudio/Novogene/250503/export_csv/SilvaDatabase/silva_nr99_v138.2_toSpecies_trainset.fa.gz.csv", row.names = TRUE
)


## Phyloseq Object Sequence ------------------
# ASV + Taxa + Sequence → csv

ls()
rm(list = ls())
# dev.off()

library(phyloseq)
library(Biostrings)
load("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Input/PhyloseqData_Bac.RData")
str(PhyseqData)


tax_table <- as.data.frame(tax_table(PhyseqData))
head(tax_table)

# DNAStringSet Objects from Phyloseq Objects
seqs <- refseq(PhyseqData)
head(names(seqs))

# Save as Fasta
writeXStringSet(seqs, filepath = "Seqence_PhyseqData_Fng.fasta")

# 配列データを列に追加
tax_table$Sequence <- as.character(seqs[rownames(tax_table)])
tax_table$Length <- nchar(as.character(seqs[rownames(tax_table)]))

otu_table <- as.data.frame(t(otu_table(PhyseqData)))
otu_tax_table <- cbind(otu_table, tax_table)

write.csv(
    x = otu_tax_table,
    file = "~/Documents/RStudio/Novogene/250503/export_csv/Bacteria_ASV_Sequence.csv", row.names = TRUE
)
