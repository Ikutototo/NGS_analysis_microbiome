# Setting -----------------------------------
setwd("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome")

ls()
rm(list = ls())
dev.off()


# load("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Input/MergeData.RData")
# load("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Input/TaxaData.RData")
# load("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Input/phylogenetic.RData")



# ConstructMetaData --------------------
## MetaData ----------------------------------
## MetaData
MetaData <- read.csv("~/Documents/RStudio/Novogene/250503/import_csv/SampleMetaData_Bacteria.csv",
    header = TRUE, sep = ",", fileEncoding = "UTF-8"
)
MetaData$dps <- factor(MetaData$dps, levels = c(0, 3, 7)) # dps列の因子型変換
MetaData$Fungicide.use <- factor(MetaData$Fungicide.use, levels = c("yes", "no")) # Fungicide.use列の因子型変換


## rownames(MetaData) = rownames(seqtab.nochim)で合わせる → phyloseq()で必要
load("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Input/250728_TaxaData.RData")
load("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Input/250728_MergeData.RData")
load("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Input/fitGTR_m4.RData")
rownames(seqtab.nochim)
rownames(MetaData) <- rownames(seqtab.nochim)
sample_names(otu_table(seqtab.nochim, taxa_are_rows = FALSE))
sample_names(sample_data(MetaData))


# cat("Class: ", class(MetaData), "\n")
# str(MetaData)
# cat("Type: ", typeof(MetaData), "\n")
#
# cat("Class: ", class(seqtab.nochim), "\n")
# str(seqtab.nochim)
# cat("Type: ", typeof(seqtab.nochim), "\n")


# Construct the phylogenetic tree -----------
search()
sessionInfo()

library("dada2")
library("msa")
library("phangorn")
library("crayon")

## msa ---------------------------------------
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs

# M4/M1 Macbookpro/Air system.time(second)
# user          system  elapsed
# 26200/43446   29/36   26228/43489
system.time(mult <- msa(seqs, method = "ClustalW", type = "dna", order = "input"))
phang.align <- as.phyDat(mult, type = "DNA", names = getSequence(seqtab.nochim))

# M4/M1 Macbookpro/Air system.time(second)
# user      system  elapsed
# 427/587   12/25   441/614
system.time(dm <- dist.ml(phang.align))


# M4/M1 Macbookpro/Air system.time(second)
# user      system      elapsed
# 593/943   0.867/1.5   593/944
system.time(treeNJ <- NJ(dm)) # Note, tip order != sequence order
system.time(fit <- pml(treeNJ, data = phang.align))


## fitGTR ------------------------------------
## negative edges length changed to 0!
system.time(fitGTR <- update(fit, k = 4, inv = 0.2))

# M4/M1 Macbookpro/Air system.time(second)
# 133768    3126    136868(38 hour)
system.time(fitGTR <- optim.pml(fitGTR,
    model = "GTR", optInv = TRUE, optGamma = TRUE,
    rearrangement = "stochastic", control = pml.control(trace = 0)
))

cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))

detach("package:phangorn", unload = TRUE)

# save(fitGTR,
#      file = "~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/芝草/NGS_consignment/Novogene/Data/250331/RData/fitGTR_done.RData")
# save(mult,phang.align,dm,treeNJ,
#      file = "~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/芝草/NGS_consignment/Novogene/Data/250331/RData/phylogenetic.RData")


# ConstructPhyloseqObjects --------------------

# library(phyloseq)
# PhyseqData <- phyloseq(
#     otu_table(seqtab.nochim, taxa_are_rows = FALSE),
#     sample_data(MetaData),
#     phyloseq::tax_table(taxa), # MicrobiotaProcessとPhyloseqでコンフリクトが起きるので注意
#     phy_tree(fitGTR$tree)
# )

# 250728 (Phy_treeは以前のを使用)
library(phyloseq)
PhyseqData <- phyloseq(
    otu_table(seqtab.nochim, taxa_are_rows = FALSE),
    sample_data(MetaData),
    phyloseq::tax_table(taxa$tax), # MicrobiotaProcessとPhyloseqでコンフリクトが起きるので注意
    phy_tree(fitGTR$tree)
)


dna <- Biostrings::DNAStringSet(taxa_names(PhyseqData))
names(dna) <- taxa_names(PhyseqData)
PhyseqData <- merge_phyloseq(PhyseqData, dna)

# 配列データの行名をASV??に変更
taxa_names(PhyseqData) <- paste0("ASV", seq(ntaxa(PhyseqData)))


# cat("Class: ", class(PhyseqData), "\n")
# str(PhyseqData)
# cat("Type: ", typeof(PhyseqData), "\n")
#
# cat("Class: ", class(dna), "\n")
# str(dna)
# cat("Type: ", typeof(dna), "\n")

## Save_PhyloseqData -------------------------

save(MetaData, PhyseqData, track,
    file = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Output/250728_PhyloseqData_Bacteria.RData"
)


# phyloseqObjectsが構築できているのかcheck
# plot_richness(PhyseqData, x = "dps", measures = c("Shannon", "Simpson"))
