# Setting -----------------------------------
setwd("~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/芝草 2/NGS_consignment/Novogene/Data/250331")
setwd("~/Documents/RStudio/Novogene/Data/250427/RawData_Fungi")
dir <- getwd()

# .fastq.gz fileを格納するdirectory
# fastq.gz fileのあるフォルダーを指定する → その都度変更する
path <- paste(dir,"clean",sep = "/" ) 
list.files(path = path, full.names = TRUE)


# lib <- c("dada2", "ggplot2", "dplyr", "Biostrings", "phyloseq", "ShortRead", "crayon")
# for (i in lib) {
#     print(i)
#     print(packageVersion(i))
#     library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
# }


# PlotQualityProfile[BEFORE] ------------------------

library(dada2)
library(Biostrings)
library(ShortRead)
library(crayon)

fnFs <- sort(list.files(path, pattern="_R1_clean.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_clean.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtFs <- file.path(path, "filtered_byfastp", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered_byfastp", paste0(sample.names, "_R_filt.fastq.gz"))
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
     file = "~/Documents/RStudio/Novogene/Data/250427/RData/Fungi_fastp/Fng_PlotQC_BEFORE_Data.RData")

cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))

# PlotQualityProfile[AFTER] -----------------


# parameterによって結果が異なるため、注意
# オブジェクトとして保存すること( <- )
Results_FilterTrim <- filterAndTrim(fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs,
                                    maxN=0, maxEE=c(2,2), truncQ=2,
                                    rm.phix=TRUE,compress=TRUE, multithread=TRUE)

PlotQC_AF_RData_For <- plotQualityProfile(filtFs[1:9])

PlotQC_AF_RData_Rev <- plotQualityProfile(filtRs[1:9])


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



save(PlotQC_AF_RData_For, PlotQC_AF_RData_Rev,
     readfiltFs_lengths_summary, readfiltRs_lengths_summary,
     file = "~/Documents/RStudio/Novogene/Data/250427/RData/Fungi_fastp/Fng_PlotQC_AFTER_Data.RData")

cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))


# Learn the Error Rates ---------------------
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
dadaRs[[1]]


save(errF, errR, dadaFs, dadaRs,
     file = "~/Documents/RStudio/Novogene/Data/250427/RData/Fungi_fastp/Fng_dada2_ErrorRateData.RData")

cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))



# Merge paired reads ------------------------

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

library(ggplot2)
ggplot(data.frame(nchar(getSequences(seqtab))), aes(x = nchar(getSequences(seqtab)))) +
    geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
    labs(title = "Distribution of Read Lengths", x = "Read Length (bp)", y = "Frequency") +
    theme_minimal()



## その時々でparameterをいじること
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 300:400] ## 今回は広めに設定(ITS2=380bpを考慮すること)

dim(seqtab2)
table(nchar(getSequences(seqtab2)))

ggplot(data.frame(nchar(getSequences(seqtab2))), aes(x = nchar(getSequences(seqtab2)))) +
    geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
    labs(title = "Distribution of Read Lengths", x = "Read Length (bp)", y = "Frequency") +
    theme_minimal()

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)
dim(seqtab.nochim)[2]/dim(seqtab2)[2]
sum(seqtab.nochim)/sum(seqtab2)


table(nchar(getSequences(seqtab.nochim)))
rowSums(seqtab.nochim)



getN <- function(x) sum(getUniques(x))
track <- cbind(Results_FilterTrim, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track


save(mergers, seqtab, seqtab2, seqtab.nochim, track,
     file = "~/Documents/RStudio/Novogene/Data/250427/RData/Fungi_fastp/Fng_MergeData.RData")

# Assign Taxonomy ---------------------------

# By applying QIIME2's classify-sklearn algorithm (Bokulich et al., 2018; Bolyen et al., 2019), 
# a pre-trained Naive Bayes classifier is used for species annotation of each ASV. 
# Annotation database of the project is Unite v9.0. [Nocogene Reports]

rm(list = ls())
load("~/Documents/RStudio/Novogene/Data/250427/RData/Fungi_fastp/Fng_MergeData.RData")


# Unite database (v.10) 

system.time(taxa <- assignTaxonomy(seqtab.nochim,
                                   refFasta = "~/Documents/RStudio/Novogene/Data/NGS_Analysis/sh_general_release_19.02.2025/sh_general_release_dynamic_19.02.2025.fasta",
                                   multithread=TRUE))

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print, 20)


# system.time(taxa <- assignTaxonomy(seqtab.nochim,
#                                    refFasta = "~/Documents/RStudio/Novogene/Data/NGS_Analysis/sh_general_release_s_19.02.2025.tgz",
#                                    multithread=TRUE))



save(taxa,
     file = "~/Documents/RStudio/Novogene/Data/250427/RData/Fungi_fastp/Fng_Taxa.RData")





cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))


# Blastn ------------------------------------

# 250427_ChatGPT
#!/bin/bash

# BLASTデータベースと入力ファイルのパスを指定
db="/path/to/your/blast/db"
input_file="your_sequences.fasta"
output_file="blast_results.txt"

# BLAST検索を実行
blastn -query $input_file -db $db -out $output_file -outfmt 6 -evalue 1e-5 -max_target_seqs 10

echo "BLAST search completed. Results saved to $output_file"


