# Setting -----------------------------------
# setwd("~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/芝草/NGS_consignment/Novogene/Data/250331")
setwd("~/Documents/RStudio/Novogene/Data/250331")
rm(list = ls())


# ConstructMetaData --------------------

library("phyloseq")

## Rにおいて、%,()は認識されないので注意
MetaData <- read.csv("~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/Novorgene_Version管理/250428/NGS_Analysis/SampleMetaData_Fungi.csv",
                     header = TRUE, sep = ",", fileEncoding = "UTF-8", stringsAsFactors = TRUE)
# levels(MetaData$dps)
# colnames(MetaData)

## dps列の因子型変換
MetaData$dps <- factor(MetaData$dps, levels = c(0, 3, 7)) 
MetaData$Fungicide.use <- factor(MetaData$Fungicide.use, levels = c("yes", "no")) # Fungicide.use列の因子型変換


## rownamesをMetaData seqtab.nochimで合わせる → phyloseq()で必要
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
# search()
# sessionInfo()


# seqtab.nochimの読み込み
# load("~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/芝草/NGS_consignment/Novogene/Data/250331/RData/PlotQC/msa.RData")
# load("~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/芝草/NGS_consignment/Novogene/Data/250331/RData/phylogenetic.RData")
# load("~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/芝草/NGS_consignment/Novogene/Data/250331/RData/fitGTR_m4.RData")

library("dada2")
library("msa")
library("phangorn")
library("crayon")

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs 

# M4/M1 Macbookpro/Air system.time(second)
## user          system  elapsed 
## 26200/43446   29/36   26228/43489 

# M4 MacPro Fungi(250427) system.time(second)
## user   system  elapsed 
## 2465   3       2468 
system.time(mult <- msa(seqs, method="ClustalW", type="dna", order="input"))


phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab.nochim))

# M4/M1 Macbookpro/Air system.time(second)
# user      system  elapsed 
# 427/587   12/25   441/614
system.time(dm <- dist.ml(phang.align))


# M4/M1 Macbookpro/Air system.time(second)
# user      system      elapsed 
# 593/943   0.867/1.5   593/944
system.time(treeNJ <- NJ(dm)) # Note, tip order != sequence order

system.time(fit <- pml(treeNJ, data = phang.align))

## negative edges length changed to 0!
system.time(fitGTR <- update(fit, k=4, inv=0.2))


# M4/M1 Macbookpro/Air system.time(second)
# 133768    3126    136868(38 hour)
system.time(fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                                rearrangement = "stochastic", control = pml.control(trace = 0)))

cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))


detach("package:phangorn", unload=TRUE)


cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))
# save(fitGTR,
#      file = "~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/芝草/NGS_consignment/Novogene/Data/250331/RData/fitGTR_done.RData")

save(mult,phang.align,dm,treeNJ,fit, fitGTR,
     file = "~/Documents/RStudio/Novogene/Data/250427/RData/Fungi_fastp/Fng_mult.RData")



# ConstructPhyseqObjects --------------------

library(phyloseq)

load("~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/Novorgene_Version管理/250428/RData/Fungi_fastp/Fng_MergeData.RData")
load("~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/Novorgene_Version管理/250428/RData/Fungi_fastp/Fng_Taxa.RData")
load("~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/Novorgene_Version管理/250428/RData/Fungi_fastp/Fng_fitGTR.RData")


PhyseqData_Fng <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                           sample_data(MetaData),
                           phyloseq::tax_table(taxa), # MicrobiotaProcessとPhyloseqでコンフリクトが起きるので注意
                           phy_tree(fitGTR$tree))


dna <- Biostrings::DNAStringSet(taxa_names(PhyseqData_Fng))
names(dna) <- taxa_names(PhyseqData_Fng)
PhyseqData_Fng <- merge_phyloseq(PhyseqData_Fng, dna)

## 配列データの行名をASVsに変更
taxa_names(PhyseqData_Fng) <- paste0("ASV", seq(ntaxa(PhyseqData_Fng)))


# cat("Class: ", class(PhyseqData_Fng), "\n")
# str(PhyseqData_Fng)
# cat("Type: ", typeof(PhyseqData_Fng), "\n")
# 
# cat("Class: ", class(dna), "\n")
# str(dna)
# cat("Type: ", typeof(dna), "\n")


save(MetaData, PhyseqData_Fng, track,
     file = "~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/Novorgene_Version管理/250428/RData/Fungi_fastp/PhyseqData_Fng.RData")


# phyloseqObjectsが構築できているのかcheck
# plot_richness(PhyseqData, x = "dps", measures = c("Shannon", "Simpson"))

# cat(crayon::bgGreen("Processing of plotqualityprofile is complete."))
