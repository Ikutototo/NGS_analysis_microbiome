# Setting -----------------------------------
setwd("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome")

ls()
rm(list = ls())
dev.off()


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
# Annotation database of the project is Unite v9.0. [Novogene Reports]

rm(list = ls())
load("~/Documents/RStudio/Novogene/Data/250427/RData/Fungi_fastp/Fng_MergeData.RData")



## Taxonomic assignment (Genus) --------------
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


## Species assignment ------------------------

# As of version 1.5.2(dada2), assignSpecies has been reimplemented to be much faster, 
# and it now easily scales to handle even very large sets of sequences.


# By default the assignSpecies method only returns species assignments if there is no ambiguity, 
# i.e. all exact matches were to the same species. 
# However, given that we are generally working with fragments of the 16S gene, 
# it is common that exact matches are made to multiple sequences that are identical over the sequenced region. 
# This is often still useful information, so to have all sequence hits returned the returnMultiple=TRUE argument can be passed to the assignSpecies function:

taxa <- assignSpecies(seqtab.nochim,
                      refFasta = "~/Documents/RStudio/Novogene/Data/NGS_Analysis/???.fa.gz",     # RDPトレーニングセットやSILVA形式など
                      outputBootstraps = TRUE, # BootstrappingScore
                      multithread = TRUE,
                      allowMultiple = TRUE
)

# If using this workflow on your own data:
# In many environments, few sequences will be assigned to species level. 
# That is OK! Reference databases are incomplete, 
# and species assignment is at the limit of what is possible from 16S amplicon data.
#  Remember, even if the sequenced organism is in the reference database, 
#  if another species shares the same 16S gene sequence it is impossible to unambiguously assign that sequence to species-level.




## Kingom to Species -------------------------
# taxa <- assignTaxonomy(seqtab.nochim,
#                        refFasta = "~/Documents/RStudio/Novogene/Data/NGS_Analysis/sh_general_release_19.02.2025/sh_general_release_dynamic_19.02.2025.fasta",
#                        multithread=TRUE)

# assignTaxonomy()されたオブジェクトに、Speciesを割り当てる
taxa_species_plus <- addSpecies(taxtab = taxa, refFasta = "~/Documents/RStudio/Novogene/Data/NGS_Analysis/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
                                verbose = TRUE, allowMultiple = TRUE)



# Blastn ------------------------------------

# 250427_ChatGPT
#!/bin/bash

# BLASTデータベースと入力ファイルのパスを指定
db="/path/to/your/blast/db"
input_file="your_sequences.fasta"
output_file="blast_results.txt"

# BLAST検索を実行
# blastn -query $input_file -db $db -out $output_file -outfmt 6 -evalue 1e-5 -max_target_seqs 10

echo "BLAST search completed. Results saved to $output_file"


# Inspects AmpliconDatabase --------------------
library(DECIPHER)
# General FASTA release (download) [https://unite.ut.ee/repository.php]
# 250708 URL 無効になっている

# UNITEによるAnnotationが適している場合
# ITS領域を用いた真菌群集解析
# 種レベルでの真菌同定が重要
# 真菌の多様性研究に特化した解析

## sh_general_release_19.02.2025 -------------

# This release consists of a single FASTA file: the RepS/RefS of all SHs, adopting the dynamically use of clustering thresholds whenever available. The format of the FASTA header is:
#     
# Claroideoglomus_sp|AM076567|SH1229972.10FU|reps|k__Fungi;p__Glomeromycota;c__Glomeromycetes;o__Entrophosporales;f__Entrophosporaceae;g__Claroideoglomus;s__Claroideoglomus_sp
# This is the file we recommend for local BLAST searches against the SHs.


# No of RepSが少ないバージョン
# SilvaDatabaseのPATHを指定する
Fungi_dna <- readDNAStringSet("~/Documents/RStudio/Novogene/250503/taxa_reference/sh_general_release_dynamic_19.02.2025.fasta")
str(Fungi_dna)

# 分類情報の抽出
names(Fungi_dna)[1:3]
Fungi_dna@ranges@NAMES

as.character(Fungi_dna)[1:3]

split_all <- strsplit(names(Fungi_dna), "\\|")
taxa <- lapply(split_all, function(x) strsplit(tail(x, 1), ";")[[1]])

# メタ情報（最初の4つ）と分類情報（k__〜s__）を連結
full_info <- mapply(function(meta, tax) {
    length(tax) <- 7  # s__まで不足ならNAで補完
    c(meta[1:4], tax) # names(Fungi_dna)のMetaDataが1~4列存在
}, split_all, taxa)

df <- as.data.frame(t(full_info), stringsAsFactors = FALSE)

colnames(df) <- c("SpeciesName", "Accession", "SH_ID", "Source",
                  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
head(df)

# 配列データ列を追加
df$Sequence <- as.character(Fungi_dna)
df$Length <- nchar(df$Sequence)

str(df)


write.csv(x = df,
          file = "~/Documents/RStudio/Novogene/250503/export_csv/SilvaDatabase/sh_general_release_dynamic_19.02.2025.fasta.csv", row.names = TRUE)


## sh_general_release_s_19.02.2025 -------------

# No of RepSが多いバージョン
# SilvaDatabaseのPATHを指定する (dev = developments??)
Fungi_dna <- readDNAStringSet("~/Documents/RStudio/Novogene/250503/taxa_reference/sh_general_release_s_19.02.2025/sh_general_release_dynamic_s_19.02.2025.fasta")
str(Fungi_dna)

# 分類情報の抽出
names(Fungi_dna)[1:3]
# Fungi_dna@ranges@NAMES

as.character(Fungi_dna)[1:3]


split_all <- strsplit(names(Fungi_dna), "\\|")
taxa <- lapply(split_all, function(x) strsplit(tail(x, 1), ";")[[1]])

# メタ情報（最初の4つ）と分類情報（k__〜s__）を連結
full_info <- mapply(function(meta, tax) {
    length(tax) <- 7  # s__まで不足ならNAで補完
    c(meta[1:4], tax) # names(Fungi_dna)のMetaDataが1~4列存在
}, split_all, taxa)

df <- as.data.frame(t(full_info), stringsAsFactors = FALSE)

colnames(df) <- c("SpeciesName", "Accession", "SH_ID", "Source",
                  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
head(df)

# 配列データ列を追加
df$Sequence <- as.character(Fungi_dna)
df$Length <- nchar(df$Sequence)

str(df)


write.csv(x = df,
          file = "~/Documents/RStudio/Novogene/250503/export_csv/SilvaDatabase/sh_general_release_dynamic_s_19.02.2025.fasta.csv", row.names = TRUE)



## SILVA_SSUfungi_nr99_v138_2_toGenus --------

# Silva taxonomic training data formatted for DADA2 (Silva version 138.2)
# https://zenodo.org/records/14169026
# 

# SilvaDatabaseのPATHを指定する
Fungi_dna <- readDNAStringSet("~/Documents/RStudio/Novogene/250503/taxa_reference/SILVA_SSUfungi_nr99_v138_2_toGenus_trainset.fasta")
str(Fungi_dna)

# 分類情報の抽出
names(Fungi_dna)[1:3]
Fungi_dna@ranges@NAMES

# 配列情報
as.character(Fungi_dna)[1:3]

# ;でsplit
split_all <- strsplit(names(Fungi_dna), "\\;")
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
df$Sequence <- as.character(Fungi_dna)
df$Length <- nchar(df$Sequence)

str(df)


write.csv(x = df,
          file = "~/Documents/RStudio/Novogene/250503/export_csv/SilvaDatabase/SILVA_SSUfungi_nr99_v138_2_toGenus_trainset.fasta.csv", row.names = TRUE)


## SILVA_SSUfungi_nr99_v138_2_toSpecies_trainset.fasta --------

# Silva taxonomic training data formatted for DADA2 (Silva version 138.2)
# https://zenodo.org/records/14169026
# 

# SilvaDatabaseのPATHを指定する
Fungi_dna <- readDNAStringSet("~/Documents/RStudio/Novogene/250503/taxa_reference/SILVA_SSUfungi_nr99_v138_2_toSpecies_trainset.fasta")
str(Fungi_dna)

# 分類情報の抽出
names(Fungi_dna)[1:3]
# Fungi_dna@ranges@NAMES

# 配列情報
as.character(Fungi_dna)[1:3]

# ;でsplit
split_all <- strsplit(names(Fungi_dna), "\\;")
split_all[1:3]


max_levels <- max(sapply(split_all, length))
# 欠損値をNAで埋める
split_padded <- lapply(split_all, function(x) {
    c(x, rep(NA, max_levels - length(x)))
})

df <- as.data.frame(do.call(rbind, split_padded), stringsAsFactors = FALSE)

colnames(df) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
head(df)

# 配列データ列を追加
df$Sequence <- as.character(Fungi_dna)
df$Length <- nchar(df$Sequence)

str(df)


write.csv(x = df,
          file = "~/Documents/RStudio/Novogene/250503/export_csv/SilvaDatabase/SILVA_SSUfungi_nr99_v138_2_toSpecies_trainset.fasta.csv", row.names = TRUE)
# 7908	Fungi	Ascomycota	Leotiomycetes	Helotiales	Sclerotiniaceae	Sclerotinia	homoeocarpa	
# CCTGGTTGATTCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAGCAACTATACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATCGTTTATTTGATAGTACCTTACTACATGGATAACCGTGGTAATTCTAGAGCTAATACATGCTAAAAACCTCGACTTCGGAAGGGGTGTGTTTATTAGATAAAAAACCAATGCCCTTCGGGGCTCCCTGGTGATTCATAATAACCTAACGAATCGCATGGCCTTGTGCCGGCGATGGTTCATTCAAATTTCTGCCCTATCAACTTTCGATGGTAGGATAGTGGCCTACCATGGTTTCAACGGGTAACGGGGAATTAGGGTTCTATTCCGGAGAAGGAGCCTGAGAAACGGCTACTACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCCGACACGGGGAGGTAGTGACAATAAATACTGATACAGGGCTCTTTTGAGTCTTGTAATTGGAATGAGTACAATTTAAATCCCTTAACGAGGAACAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGAACCTTGGGCCTGGCTGACCGGTCCGCCTCACCGCGTGCACTGGTTCGGCCGGGCCTTTCCTTCTGGGGAGCCGCATGCCCTTCACTGGGTGTGTCGGGGAACCAGGACTTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCCTATGCTCGAATACATTAGCATGGAATAATAGAATAGGACGTGTGGTTCTATTTTGTTGGTTTCTAGGACCGCCGTAATGATTAATAGGGATAGTCGGGGGCATCAGTATTCAATTGTCAGAGGTGAAATTCTTGGATTTATTGAAGACTAACTACTGCGAAAGCATTTGCCAAGGATGTTTTCATTAATCAGTGAACGAAAGTTAGGGGATCGAAGACGATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTAGGGATCGGGCGATGTTATCTTTTTGACTCGCTCGGCACCTCACGAGAAATCAAAGTTTTTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTGAAACTTAAAGAAATTGACGGAAAGGCACCACCAGGCGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGGAAACTCACCAGGTTAACCACGGTTGTTACGACCTCTGGGCCTGGAAAAGAAGGGGGGTGGCCCACCTCTCTCTAGTGCTTGTCTTTGTCTGTGTGGGGAAGTCCCCTATTTTGGGCACAGACGCTCCGTAGCGGGAGCGTGACAGGTGCAACACCAGCTGGAACAGAAGACGCCTCCGTTACATGTAACGAAGCCAATTCTGTGGCGAGCCTGGGTCACGCCAGGCCGTCGCAACGCGCGCAAAGCGGTGGGTTCACTGAATGCAGTGGGCTTAAGGTACGTGCTAATCCCGCGAGAAATCGCGCCGCGTGAACAAGGTCCAAAAGCCAAAGTCACGCGGGCCTATCATCTGATAGGCGGTATTTGCGGGGAGTGCCCCAGCACCCTCTCTCGATGGAGGGATGATGCGGGGGGGCTCCTCGACATGCCAGACACAATAAGGATTGACAGATTGAGAGCTCTTTCTTGATTTTGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGTGATTTGTCTGCTTAATTGCGATAACGAACGAGACCTTAACCTGCTAAATAGCCAGGCTAGCTTTGGCTGGTCGCCGGCTTCTTAGAGGGACTATCGGCTCAAGCCGATGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGACAGAGCCAACGAGTTCATCTCCTTGACCGAAAGGTCTGGGTAATCTTGTTAAACTCTGTCGTGCTGGGGATAGAGCATTGCAATTATTGCTCTTCAACGAGGAATGCCTAGTAAGCGCAAGTCATCAGCTTGCGTTGATTACGTCCCTGCCTTTTGTACACACCGCCCGTCGCTACTACCGATTGAATGGCTAAGTGAGGCTTTCGGACTGGCCTAGGGAGGGTGGCAACACCCACCCAGGGCCGGAAAGTTGTCCAAACTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAG
# 2204

## Filtering_Taxa ----------------------------
# 特定の分類群を指定して、検索
id_Ascomycota <- grep("p__Ascomycota", names(Fungi_dna))
Ascomycota_dna <- Fungi_dna[id_Ascomycota]
length(Ascomycota_dna)
names(Ascomycota_dna)[1:5]


write.csv(x = Ascomycota_dna@ranges@NAMES,
          file = "~/Documents/RStudio/Novogene/250503/export_csv/SilvaDatabase/id_Ascomycota.csv", row.names = TRUE)


## FastqData Alignments ---------------------------------
library(Biostrings)

# Query配列をDNAStringとして定義
query_seq <- DNAString("GAAATGCGATACGTAATGTGAATTGCAAATTCAGTGAATCATCGAGTCTTTGAACGCACATTGCGCCCCCTGGTATTCCGGGGGGCATGCCTGTCCGAGCGTCATTGCTGCCCTCAAGCCCGGCTTGTGTGTTGGGCCCCGTCCTCCGATTCCGGGGGACGGGCCCGAAAGGCAGCGGCGGCACCGCGTCCGGTCCTCGAGCGTATGGGGCTTTGTCACCCGCTCTGTAGGCCCGGCCGGCGCTTGCCGATCAACCCAAATTTTTATCCAGGTTGACCTCGGATCAGGTAGGGATACCCGCTGAACTTAA")

# Fungi_dna との局所アラインメントをとる
alignments <- vmatchPattern(query_seq, Fungi_dna, max.mismatch=10)

hit_names <- names(Fungi_dna)[elementNROWS(alignments) > 0]

taxonomy_hits <- sub(".*\\|refs\\|", "", hit_names)
taxonomy_hits <- sub("\\|.*", "", taxonomy_hits)
taxonomy_hits

## Phyloseq Object Sequence ------------------
# ASV + Taxa + Sequence → csv

ls()
rm(list = ls())
# dev.off()

library(phyloseq)
library(Biostrings)
load("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Fungi/Input/PhyseqData_Fng.RData")
str(PhyseqData_Fng)


tax_table <- as.data.frame(tax_table(PhyseqData_Fng))
head(tax_table)

# DNAStringSet Objects from Phyloseq Objects
seqs <- refseq(PhyseqData_Fng)
head(names(seqs))

# Save as Fasta
writeXStringSet(seqs, filepath = "Seqence_PhyseqData_Fng.fasta")

# 配列データを列に追加
tax_table$Sequence <- as.character(seqs[rownames(tax_table)])
tax_table$Length <- nchar(as.character(seqs[rownames(tax_table)]))

otu_table <- as.data.frame(t(otu_table(PhyseqData_Fng)))
otu_tax_table <- cbind(otu_table, tax_table)

write.csv(x = otu_tax_table,
          file = "~/Documents/RStudio/Novogene/250503/export_csv/Fungi_ASV_Sequence.csv", row.names = TRUE)




