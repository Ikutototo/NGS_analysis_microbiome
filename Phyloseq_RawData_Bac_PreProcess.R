# Setting -----------------------------------

## Gitで管理しているなら、Pullして、最新版をローカルにアップデートしてから始めること!!

ls()
rm(list = ls())

setwd("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome")


library(phyloseq)
library(ggplot2)
library(cowplot)

load("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Output/PhyloseqData_Bacteria.RData")


# Prevalence filtering ----------------------

load("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Plot/Prevelence_inTaxa.RData")

## Sample(Bac1~9)において、各ASVsが出現した(>0)Sample数(存在頻度)を出力
## taxa_are_rows → ASVsが行が列かどうかを確認
prev0 = apply(X = otu_table(PhyseqData),
              MARGIN = ifelse(taxa_are_rows(PhyseqData), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

## DataFrameの作成
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(PhyseqData),
                    tax_table(PhyseqData))

prevdf = subset(prevdf, Kingdom == "Bacteria")
# write.csv(x = prevdf, file = "~/Documents/RStudio/Novogene/250503/export_csv/prevdf.csv")

## 全ASVsのベクトル(ID)を取り出す → TotalReads(全SampleのASVs) > 10でフィルタリング
## keep_taxa <- taxa_names(PhyseqData)[taxa_sums(PhyseqData) >= 10]
## ps_filt <- prune_taxa(keep_taxa, PhyseqData)


## 各分類(Phylum)レベルで、tableを出力、各分類(Phylum)レベルで分類されたベクトルに対して、[ASVs数>5]でFiltering
## ASVsが全Sampleの中で、5回以上出現していないといけない 
## 例: ASV1 Bac1~4のみで>0の場合、除去される
## → 比較したい処理区内でのSample数を考慮すること


## 全sample数のうち、何sampleで同一のASVsが出現していれば、良いかの基準を設定
## 0.1の割合だと、9sampleの場合、0.9となり、1sampleで出現していればいいことになるため、1で設定
prevalenceThreshold = 1
prevalenceThreshold

## Execute prevalence filter, using `prune_taxa()` function
# prevdf_phylum_filt = phyloseq::prune_taxa((prev0 > prevalenceThreshold), PhyseqData)
# prevdf_phylum_filt
# 
# 
# prevdf_phylum_subset = subset_taxa(prevdf_phylum_filt, Phylum %in% names(keepPhyla))
# prevdf_phylum_subset

## Phylum ------------------------------------

keepPhyla = table(prevdf$Phylum)[(table(na.omit(prevdf$Phylum)) > 5)] # 全Sample中、5つ以上でPrelvalenceと指定
prevdf_phylum = subset(prevdf, Phylum %in% names(keepPhyla)) # names()は大事
table(prevdf_phylum$Phylum)

write.csv(x = prevdf_phylum,
          file = "~/Documents/RStudio/Novogene/250503/export_csv/prevdf_Phylum.csv", row.names = TRUE)


### Prevalence Plots → Classでラベリングしています
(prevdf_phylum_prevelence_plot <- 
    ggplot(prevdf_phylum, aes(TotalAbundance, Prevalence, color = Class)) +
    geom_hline(yintercept = prevalenceThreshold, alpha = 0.7, linetype = 2, colour = "black") +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_x_log10() +
    scale_y_continuous(breaks = seq(0, 9, by = 1)) +
    xlab("Total Abundance in Phylum") + ylab("Sample Prevalence") +
    facet_wrap(~Phylum) +
    theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
          axis.text = element_text(size = 10, face = "bold", color = "black"),
          strip.text = element_text(size = 8, face = "bold.italic", color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = "none"))

ggsave(filename = "prevdf_phylum_prevelence_plot.png", plot = last_plot(), 
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png/")

Phylum_legend_only <- get_legend(
    ggplot(prevdf_phylum, aes(TotalAbundance, Prevalence, color = Class)) +
        geom_hline(yintercept = prevalenceThreshold, alpha = 0.7, linetype = 2, colour = "black") +
        geom_point(size = 1.5, alpha = 0.8) +
        scale_x_log10() +
        xlab("Total Abundance in Phylum") + ylab("Prevalence") +
        facet_wrap(~Phylum) +
        theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
              axis.text = element_text(size = 12, face = "bold", color = "black"),
              strip.text = element_text(size = 10, face = "bold.italic", color = "black"),
              panel.grid.minor = element_blank()))

ggdraw(Phylum_legend_only)
ggsave(filename = "prevdf_phylum_prevelence_plot_legend_only.png", plot = last_plot(), 
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png/")

dev.off()

## Class -------------------------------------

keepClass = table(prevdf$Class)[(table(prevdf$Class) > 5)]
prevdf_Class = subset(prevdf, Class %in% names(keepClass))
table(prevdf_Class$Class)

write.csv(x = prevdf_Class,
          file = "~/Documents/RStudio/Novogene/250503/export_csv/prevdf_Class.csv", row.names = TRUE)



class_counts_Class <- sort(table(na.omit(prevdf_Class$Class)), decreasing = TRUE)
length(class_counts_Class) 
Class_chunks <- split(class_counts_Class, ceiling(seq_along(class_counts_Class) / 9)) # facet(~??)で分割したい数を指定


plot_list_Class <- list() 
for (i in seq_along(Class_chunks)) {
    subset_df <- subset(prevdf_Class, Class %in% names(Class_chunks[[i]]))
    
    if (nrow(subset_df) == 0) next
    
    p <- ggplot(subset_df, aes(TotalAbundance, Prevalence, color = Family)) +
        geom_hline(yintercept = prevalenceThreshold, alpha = 0.7, linetype = 2, colour = "black") +
        geom_point(size = 1.5, alpha = 0.8) +
        scale_x_log10() +
        xlab("Total Abundance in Class") + ylab("Prevalence") +
        facet_wrap(~Class) +
        theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
              axis.text = element_text(size = 12, face = "bold", color = "black"),
              strip.text = element_text(size = 14, face = "italic", color = "black"),
              panel.grid.minor = element_blank(),
              legend.position = "none")
    
    plot_list_Class[[i]] <- p
}

### Legendのみのオブジェクト
plot_list_Class_legend <- list() 
for (i in seq_along(Class_chunks)) {
    subset_df <- subset(prevdf_Class, Class %in% names(Class_chunks[[i]]))
    
    if (nrow(subset_df) == 0) next
    
    p <- cowplot::get_legend(
        ggplot(subset_df, aes(TotalAbundance, Prevalence, color = Family)) +
        geom_hline(yintercept = prevalenceThreshold, alpha = 0.7, linetype = 2, colour = "black") +
        geom_point(size = 1.5, alpha = 0.8) +
        scale_x_log10() +
        xlab("Total Abundance in Class") + ylab("Prevalence") +
        facet_wrap(~Class) +
        theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
              axis.text = element_text(size = 12, face = "bold", color = "black"),
              strip.text = element_text(size = 14, face = "italic", color = "black"),
              panel.grid.minor = element_blank()))
    
    plot_list_Class_legend[[i]] <- p
}

grid::grid.newpage()
grid::grid.draw(plot_list_Class_legend[[1]]) # ggdraw()にlistを渡さないこと、[[]]で渡す



## Order -------------------------------------

keepOrder = table(prevdf$Order)[(table(prevdf$Order) > 5)]
prevdf_Order = subset(prevdf, Order %in% names(keepOrder))
table(prevdf_Order$Order)

write.csv(x = prevdf_Order,
          file = "~/Documents/RStudio/Novogene/250503/export_csv/prevdf_Order.csv", row.names = TRUE)



class_counts_Order <- sort(table(na.omit(prevdf_Order$Order)), decreasing = TRUE)
length(class_counts_Order) 
Order_chunks <- split(class_counts_Order, ceiling(seq_along(class_counts_Order) / 9)) # facet(~??)で分割したい数を指定


plot_list_Order <- list() 
for (i in seq_along(Order_chunks)) {
    subset_df <- subset(prevdf_Order, Order %in% names(Order_chunks[[i]]))
    
    if (nrow(subset_df) == 0) next
    
    p <- ggplot(subset_df, aes(TotalAbundance, Prevalence, color = Family)) +
        geom_hline(yintercept = prevalenceThreshold, alpha = 0.7, linetype = 2, colour = "black") +
        geom_point(size = 1.5, alpha = 0.8) +
        scale_x_log10() +
        xlab("Total Abundance in Order") + ylab("Prevalence") +
        facet_wrap(~Order) +
        theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
              axis.text = element_text(size = 12, face = "bold", color = "black"),
              strip.text = element_text(size = 14, face = "italic", color = "black"),
              panel.grid.minor = element_blank(),
              legend.position = "none")
    
    plot_list_Order[[i]] <- p
}

### Legendのみのオブジェクト
plot_list_Order_legend <- list() 
for (i in seq_along(Order_chunks)) {
    subset_df <- subset(prevdf_Order, Order %in% names(Order_chunks[[i]]))
    
    if (nrow(subset_df) == 0) next
    
    p <- cowplot::get_legend(
        ggplot(subset_df, aes(TotalAbundance, Prevalence, color = Family)) +
            geom_hline(yintercept = prevalenceThreshold, alpha = 0.7, linetype = 2, colour = "black") +
            geom_point(size = 1.5, alpha = 0.8) +
            scale_x_log10() +
            xlab("Total Abundance in Order") + ylab("Prevalence") +
            facet_wrap(~Order) +
            theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
                  axis.text = element_text(size = 12, face = "bold", color = "black"),
                  strip.text = element_text(size = 14, face = "italic", color = "black"),
                  panel.grid.minor = element_blank()))
    
    plot_list_Order_legend[[i]] <- p
}


## Family ------------------------------------

keepFamily = table(prevdf$Family)[(table(prevdf$Family) > 5)]
prevdf_Family = subset(prevdf, Family %in% names(keepFamily))
table(prevdf_Family$Family)

class_counts_Family <- sort(table(na.omit(prevdf_Family$Family)), decreasing = TRUE)
length(class_counts_Family) 
Family_chunks <- split(class_counts_Family, ceiling(seq_along(class_counts_Family) / 9)) # facet(~??)で分割したい数を指定


plot_list_Family <- list() 
for (i in seq_along(Family_chunks)) {
    subset_df <- subset(prevdf_Family, Family %in% names(Family_chunks[[i]]))
    
    if (nrow(subset_df) == 0) next
    
    p <- ggplot(subset_df, aes(TotalAbundance, Prevalence, color = Family)) +
        geom_hline(yintercept = prevalenceThreshold, alpha = 0.7, linetype = 2, colour = "black") +
        geom_point(size = 1.5, alpha = 0.8) +
        scale_x_log10() +
        xlab("Total Abundance in Family") + ylab("Prevalence") +
        facet_wrap(~Family) +
        theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
              axis.text = element_text(size = 12, face = "bold", color = "black"),
              strip.text = element_text(size = 14, face = "italic", color = "black"),
              panel.grid.minor = element_blank(),
              legend.position = "none")
    
    plot_list_Family[[i]] <- p
}

### Legendのみのオブジェクト
plot_list_Family_legend <- list() 
for (i in seq_along(Family_chunks)) {
    subset_df <- subset(prevdf_Family, Family %in% names(Family_chunks[[i]]))
    
    if (nrow(subset_df) == 0) next
    
    p <- cowplot::get_legend(
        ggplot(subset_df, aes(TotalAbundance, Prevalence, color = Family)) +
            geom_hline(yintercept = prevalenceThreshold, alpha = 0.7, linetype = 2, colour = "black") +
            geom_point(size = 1.5, alpha = 0.8) +
            scale_x_log10() +
            xlab("Total Abundance in Family") + ylab("Prevalence") +
            facet_wrap(~Family) +
            theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
                  axis.text = element_text(size = 12, face = "bold", color = "black"),
                  strip.text = element_text(size = 14, face = "italic", color = "black"),
                  panel.grid.minor = element_blank()))
    
    plot_list_Family_legend[[i]] <- p
}


## Genus -------------------------------------

keepGenus = table(prevdf$Genus)[(table(prevdf$Genus) > 5)] 
prevdf_Genus = subset(prevdf, Genus %in% names(keepGenus))
table(prevdf_Genus$Genus)

class_counts_Genus <- sort(table(na.omit(prevdf_Genus$Genus)), decreasing = TRUE)
length(class_counts_Genus) 
Genus_chunks <- split(class_counts_Genus, ceiling(seq_along(class_counts_Genus) / 9)) # facet(~??)で分割したい数を指定

plot_list_Genus <- list() 
for (i in seq_along(Genus_chunks)) {
    subset_df <- subset(prevdf_Genus, Genus %in% names(Genus_chunks[[i]]))
    
    if (nrow(subset_df) == 0) next
    
    p <- ggplot(subset_df, aes(TotalAbundance, Prevalence, color = Genus)) +
        geom_hline(yintercept = prevalenceThreshold, alpha = 0.7, linetype = 2, colour = "black") +
        geom_point(size = 1.5, alpha = 0.8) +
        scale_x_log10() +
        xlab("Total Abundance in Genus") + ylab("Prevalence") +
        facet_wrap(~Genus) +
        theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
              axis.text = element_text(size = 12, face = "bold", color = "black"),
              strip.text = element_text(size = 14, face = "italic", color = "black"),
              panel.grid.minor = element_blank(),
              legend.position = "none")
    
    plot_list_Genus[[i]] <- p
}



## Save_PrevalenceFilteringPlots ------------------


save(prevdf_phylum_prevelence_plot,Phylum_legend_only,
     plot_list_Class, plot_list_Class_legend,
     plot_list_Order, plot_list_Order_legend,
     plot_list_Family, plot_list_Family_legend, 
     plot_list_Genus,
     file = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Plot/Prevelence_inTaxa.RData")


dev.off()
rm(list = ls())




# Abundance value Counts & refseq ------------

## 有意な差が認められた分類群の存在量を個別で可視化
## ASVsの配列情報を取得
## 全Sampleの内、Countsされた数でFiltering


## Genus == Bryobacter → Filtering → Sequenceを取得 → Data.frame化
write.csv(data.frame(Sequence = refseq(PhyseqData)[taxa_names(subset_taxa(PhyseqData, Genus == "Bryobacter"))],
                     ASV_ID = taxa_names(subset_taxa(PhyseqData, Genus == "Bryobacter")),
                     stringsAsFactors = FALSE),
          file = "~/Documents/RStudio/Novogene/250503/export_csv/Bryobacter_sequences.csv", row.names = FALSE)





## ASV205： -----------------------------------

colnames(psmelt(PhyseqData))
class(psmelt(PhyseqData))
length(psmelt(PhyseqData) |> dplyr::filter(Genus == "Bryobacter"))

PhyseqData |> 
    psmelt() |>
    dplyr::filter(Genus == "Bryobacter") |> 
    write.csv(file = "~/Documents/RStudio/Novogene/250503/export_csv/psmelt_Bryobacter.csv")


prune_taxa("ASV205", PhyseqData) |> 
    psmelt() |> 
    group_by(dps) |> 
    summarise(
        mean_abund = mean(Abundance),
        sd_abund = sd(Abundance)) |> 
ggplot(aes(x = dps, y = mean_abund)) +
    geom_bar(stat = "identity", width = 0.8, fill = "steelblue") +
    geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), width = 0.1) +
    ylab("ASV 205 Bryobacter Abundance") +
    theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
          axis.text = element_text(size = 16, face = "bold", color = "black"),
          strip.text = element_text(size = 14, face = "bold.italic", color = "black"),
          panel.grid.minor = element_blank())

ggsave(filename = "ASV205_Bryobacter_dps_Abundance_Plots.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")

dev.off()



## ASV75： ------------------------------------

prune_taxa("ASV75", PhyseqData) |> 
    psmelt() |> 
    group_by(dps) |> 
    summarise(
        mean_abund = mean(Abundance),
        sd_abund = sd(Abundance)) |> 
    ggplot(aes(x = dps, y = mean_abund)) +
    geom_bar(stat = "identity", width = 0.8, fill = "steelblue") +
    geom_errorbar(aes(ymin = pmax(mean_abund - sd_abund, 0),
                      ymax = mean_abund + sd_abund), width = 0.1) +
    ylab("ASV 75 Thiobacillus Abundance") +
    theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
          axis.text = element_text(size = 16, face = "bold", color = "black"),
          strip.text = element_text(size = 14, face = "bold.italic", color = "black"),
          panel.grid.minor = element_blank())

ggsave(filename = "ASV75_Thiobacillus_dps_Abundance_Plots.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")



### Others --------------------------------


# prune_samples(sample_sums(physeqData_???)>=20, physeqData_???)
# ExportCSVData <- cbind(t(PhyseqData@otu_table@.Data), PhyseqData@tax_table@.Data)

# Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
# filter_taxa(PhyseqData, function(x) sum(x > 3) > (0.2*length(x)), TRUE)


# Standardize abundances to the median sequencing depth
# sample間で異なるリード数(シーケンス深度)を一致させるため
# total = median(sample_sums(physeqData))
# standf = function(x, t=total) round(t * (x / sum(x)))
# physeqData_sd = transform_sample_counts(physeqData, standf)


# Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation
# physeqData_cv = filter_taxa(PhyseqData_???, function(x) sd(x)/mean(x) > 3.0, TRUE)



# FilteringMethods --------------------------
## RelativeAbundunce -------------------------
### Color_Setting -----------------------------

library(phyloseq)
library(ggplot2)
library(cowplot)
library(dplyr)
library(scales)

colors <- c(
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
    "#e6ab02", "#a6761d", "#666666", "#8dd3c7", "#ffffb3",
    "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", 
    "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f", "#7fc97f", 
    "#fdc086", "#ffff99", "#386cb0", "#f0027f", "#bf5b17", 
    "#6a3d9a", "#cab2d6", "#ff7f00", "#b2df8a", "#a6cee3", 
    "#fb9a99", "#1f78b4", "#33a02c", "#b15928")

color_v2 <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B",
              "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#F0E442", "#9E14D2",
              "#7D8B33", "#B5F56D", "#3F51B5", "#D32F2F", "#0288D1", "#7B1FA2", 
              "#388E3C", "#FBC02D", "#0288D1", "#8BC34A", "#FF5722", "#8E24AA", 
              "#795548", "#9C27B0", "#3F51B5", "#4CAF50", "#FF9800", "#E91E63", 
              "#CDDC39")

colors_v3 <- c(colors, "#9999FF", "#FF99CC", "#66FFFF")


### Phylum Level -------------------------------

PhyseqData_Phylum <- PhyseqData  |> 
    subset_taxa(Kingdom == "Bacteria") |> 
    tax_glom(taxrank = "Phylum") %>%                        # agglomerate at phylum level
    transform_sample_counts(function(x) {x/sum(x)} )  |>    # Transform to relative abundance
    psmelt()  |>                                            # Melt to long format
    filter(Abundance > 0.01)  |>                            # Filter out low(>1%) abundance taxa
    arrange(desc(Phylum))

unique(psmelt(PhyseqData)$Phylum)
length(unique(psmelt(PhyseqData)$Phylum))

Descending_Phylum <- PhyseqData_Phylum |> 
    group_by(Phylum)  |> 
    summarise(total_abundance = sum(Abundance, na.rm = TRUE))  |> 
    arrange(total_abundance)  |> 
    pull(Phylum)

PhyseqData_Phylum$Phylum <- factor(PhyseqData_Phylum$Phylum, levels = Descending_Phylum)

ggplot(PhyseqData_Phylum, aes(x = dps, y = Abundance, fill = Phylum)) + 
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1.45)) +
    scale_y_continuous(labels = percent) +
    xlab("Days post Fungicide") +
    ylab("Relative Abundance (Phylum > 1%) \n") +
    theme(axis.title = element_text(size = 14, face = "bold", color = "black"),
          axis.text = element_text(size = 12, face = "bold", color = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
          legend.background = element_rect(fill = "gray90"),
          legend.key = element_rect(fill = "white", color = NA))

ggsave(filename = "Relative_abundance_Phylum.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")


ggplot(PhyseqData_Phylum, aes(x = Sample.Name, y = Abundance, fill = Phylum)) + 
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(reverse = F, keywidth = 1.45, keyheight = 1.45)) +
    scale_y_continuous(labels = percent) +
    xlab("Samples") +
    ylab("Relative Abundance (Phylum > 1%) \n") +
    theme(axis.title = element_text(size = 16, face = "bold",color = "black"),
          axis.text = element_text(size = 10, color = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
          legend.background = element_rect(fill = "gray90"),
          legend.key = element_rect(fill = "white", color = NA),
          legend.key.size = unit(3, "cm"))

ggsave(filename = "Relative_abundance_Phylum_SampleName.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")


### Class Level -------------------------------

PhyseqData_Class <- PhyseqData  |> 
    subset_taxa(Kingdom == "Bacteria") |> 
    tax_glom(taxrank = "Class")  |>                         # agglomerate at phylum level
    transform_sample_counts(function(x) {x/sum(x)} )  |>    # Transform to relative abundance
    psmelt()  |>                                            # Melt to long format
    filter(Abundance > 0.01)  |>                            # Filter out low(>1%) abundance taxa
    arrange(desc(Class))

unique(psmelt(PhyseqData)$Class)
length(unique(psmelt(PhyseqData)$Class))

Descending_Class <- PhyseqData_Class  |> 
    group_by(Class)  |> 
    summarise(total_abundance = sum(Abundance, na.rm = TRUE))  |> 
    arrange(total_abundance)  |> 
    pull(Class)

PhyseqData_Class$Class <- factor(PhyseqData_Class$Class, levels = Descending_Class)

ggplot(PhyseqData_Class, aes(x = dps, y = Abundance, fill = Class)) + 
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1.45)) +
    scale_y_continuous(labels = percent) +
    xlab("Days post Fungicide") +
    ylab("Relative Abundance (Class > 1%) \n") +
    theme(axis.title = element_text(size = 14, face = "bold", color = "black"),
          axis.text = element_text(size = 12, face = "bold", color = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
          legend.background = element_rect(fill = "gray90"),
          legend.key = element_rect(fill = "white", color = NA))

ggsave(filename = "Relative_abundance_Class.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")


ggplot(PhyseqData_Class, aes(x = Sample.Name, y = Abundance, fill = Class)) + 
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(reverse = F, keywidth = 1.45, keyheight = 1.45)) +
    scale_y_continuous(labels = percent) +
    xlab("Samples") +
    ylab("Relative Abundance (Class > 1%) \n") +
    theme(axis.title = element_text(size = 16, face = "bold",color = "black"),
          axis.text = element_text(size = 10, color = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
          legend.background = element_rect(fill = "gray90"),
          legend.key = element_rect(fill = "white", color = NA),
          legend.key.size = unit(3, "cm"))

ggsave(filename = "Relative_abundance_Class_SampleName.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")



### Order Level -------------------------------

PhyseqData_Order <- PhyseqData  |> 
    subset_taxa(Kingdom == "Bacteria") |> 
    tax_glom(taxrank = "Order")  |>                         # agglomerate at phylum level
    transform_sample_counts(function(x) {x/sum(x)} )  |>    # Transform to relative abundance
    psmelt()  |>                                            # Melt to long format
    filter(Abundance > 0.01)  |>                            # Filter out low(>1%) abundance taxa
    arrange(desc(Order))

unique(psmelt(PhyseqData)$Order)
length(unique(psmelt(PhyseqData)$Order))

Descending_Order <- PhyseqData_Order  |> 
    group_by(Order)  |> 
    summarise(total_abundance = sum(Abundance, na.rm = TRUE))  |> 
    arrange(total_abundance)  |> 
    pull(Order)

PhyseqData_Order$Order <- factor(PhyseqData_Order$Order, levels = Descending_Order)

ggplot(PhyseqData_Order, aes(x = dps, y = Abundance, fill = Order)) + 
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1.45)) +
    scale_y_continuous(labels = percent) +
    xlab("Days post Fungicide") +
    ylab("Relative Abundance (Order > 1%) \n") +
    theme(axis.title = element_text(size = 14, face = "bold", color = "black"),
          axis.text = element_text(size = 12, face = "bold", color = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
          legend.background = element_rect(fill = "gray90"),
          legend.key = element_rect(fill = "white", color = NA))

ggsave(filename = "Relative_abundance_Order.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")


ggplot(PhyseqData_Order, aes(x = Sample.Name, y = Abundance, fill = Order)) + 
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(reverse = F, keywidth = 1.45, keyheight = 1.45)) +
    scale_y_continuous(labels = percent) +
    xlab("Samples") +
    ylab("Relative Abundance (Order > 1%) \n") +
    theme(axis.title = element_text(size = 16, face = "bold",color = "black"),
          axis.text = element_text(size = 10, color = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
          legend.background = element_rect(fill = "gray90"),
          legend.key = element_rect(fill = "white", color = NA),
          legend.key.size = unit(3, "cm"))

ggsave(filename = "Relative_abundance_Order_SampleName.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")


### Family Level  -----------------------------

PhyseqData_Family <- PhyseqData  |> 
    subset_taxa(Kingdom == "Bacteria") |> 
    tax_glom(taxrank = "Family")  |>                         # agglomerate at phylum level
    transform_sample_counts(function(x) {x/sum(x)} )  |>    # Transform to relative abundance
    psmelt()  |>                                            # Melt to long format
    filter(Abundance > 0.01)  |>                            # Filter out low(>1%) abundance taxa
    arrange(desc(Family))
    
unique(psmelt(PhyseqData)$Family)
length(unique(psmelt(PhyseqData)$Family))

Descending_Family <- PhyseqData_Family  |> 
    group_by(Family)  |> 
    summarise(total_abundance = sum(Abundance, na.rm = TRUE))  |> 
    arrange(total_abundance)  |> 
    pull(Family)

PhyseqData_Family$Family <- factor(PhyseqData_Family$Family, levels = Descending_Family)

ggplot(PhyseqData_Family, aes(x = dps, y = Abundance, fill = Family)) + 
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1.45)) +
    scale_y_continuous(labels = percent) +
    xlab("Days post Fungicide") +
    ylab("Relative Abundance (Family > 1%) \n") +
    theme(axis.title = element_text(size = 14, face = "bold", color = "black"),
          axis.text = element_text(size = 12, face = "bold", color = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
          legend.background = element_rect(fill = "gray90"),
          legend.key = element_rect(fill = "white", color = NA))

ggsave(filename = "Relative_abundance_Family.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")
       
       
ggplot(PhyseqData_Family, aes(x = Sample.Name, y = Abundance, fill = Family)) + 
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(reverse = F, keywidth = 1.45, keyheight = 1.45)) +
    scale_y_continuous(labels = percent) +
    xlab("Samples") +
    ylab("Relative Abundance (Family > 1%) \n") +
    theme(axis.title = element_text(size = 16, face = "bold",color = "black"),
          axis.text = element_text(size = 10, color = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
          legend.background = element_rect(fill = "gray90"),
          legend.key = element_rect(fill = "white", color = NA),
          legend.key.size = unit(3, "cm"))

ggsave(filename = "Relative_abundance_Family_SampleName.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")


### Genus Level  ------------------------------
PhyseqData_Genus <- PhyseqData  |> 
    subset_taxa(Kingdom == "Bacteria") |> 
    tax_glom(taxrank = "Genus") |>                        
    transform_sample_counts(function(x) {x/sum(x)} )  |>   
    psmelt()  |>                                           
    filter(Abundance > 0.01)  |>                           
    arrange(desc(Genus))


unique(psmelt(PhyseqData)$Genus)
length(unique(psmelt(PhyseqData)$Genus))


Descending_Genus <- PhyseqData_Genus  |> 
    group_by(Genus)  |> 
    summarise(total_abundance = sum(Abundance, na.rm = TRUE))  |> 
    arrange(total_abundance)  |> 
    pull(Genus)

PhyseqData_Genus$Genus <- factor(PhyseqData_Genus$Genus, levels = Descending_Genus)

ggplot(PhyseqData_Genus, aes(x = dps, y = Abundance, fill = Genus)) + 
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colors_v3) +
    guides(fill = guide_legend(reverse = F, keywidth = 0.7, keyheight = 1.45)) +
    scale_y_continuous(labels = percent) +
    xlab("Days post Fungicide") +
    ylab("Relative Abundance (Genus > 1%) \n") +
    theme(axis.title = element_text(size = 14, face = "bold", color = "black"),
          axis.text = element_text(size = 10, face = "bold", color = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
          legend.background = element_rect(fill = "gray90"),
          legend.key = element_rect(fill = "white", color = NA))

ggsave(filename = "Relative_abundance_Genus.png", plot = last_plot(),
       width = 2000, height = 1800, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")


ggplot(PhyseqData_Genus, aes(x = Sample.Name, y = Abundance, fill = Genus)) + 
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colors_v3) +
    guides(fill = guide_legend(reverse = F, keywidth = 1.45, keyheight = 1.45)) +
    scale_y_continuous(labels = percent) +
    xlab("Samples") +
    ylab("Relative Abundance (Genus > 1%) \n") +
    theme(axis.title = element_text(size = 16, face = "bold",color = "black"),
          axis.text = element_text(size = 10, color = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
          legend.background = element_rect(fill = "gray90"),
          legend.key = element_rect(fill = "white", color = NA),
          legend.key.size = unit(3, "cm"))

ggsave(filename = "Relative_abundance_Genus_SampleName.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")




## Top 50 Filtering --------------------------

# Selects Top 50 
top50_taxa <- names(sort(taxa_sums(PhyseqData), decreasing = TRUE)[1:50])
prune_taxa(top_taxa, PhyseqData)



## Top50_Plot_Relative_Abundunce -------------

### 前処理無し → 存在量の可視化
plot_bar(PhyseqData, x = "Sample.Name", fill = "Phylum") + scale_fill_igv()

colnames(psmelt(PhyseqData)) 

### Visualize patterns (scatterplot)
psmelt(PhyseqData) |> 
    select(OTU, Sample, Abundance, dps) |> 
    filter(OTU == top50_taxa[1]) |> # 全Sampleで、最も存在量が多いASVs
    mutate(dps = factor(dps)) |> 
    (\(df) ggplot(df, aes(y = Abundance, x = dps)) + # 無名関数
         geom_boxplot(outlier.shape = NA, width = 0.2) +
         geom_jitter(height = 0, width = 0.2))()


# Richness Index ----------------------------
library(phyloseq)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

## Exports_OTU_Table -------------------------
## ASVs Tableの生成と出力
otu_mat <- t(as(otu_table(PhyseqData), Class = "matrix"))
tax_mat <- as(tax_table(PhyseqData), Class = "matrix")
otu_table <- cbind(otu_mat, tax_mat)
write.csv(x = otu_table, file = "~/Documents/RStudio/Novogene/250503/export_csv/otu_table.csv",
          row.names = TRUE)

## シングルトン有無の確認
sum(as(otu_table(PhyseqData), Class = "matrix") == 1)
tail(phyloseq::taxa_sums(PhyseqData))
plot_richness(PhyseqData, nrow = 3)

ggsave(filename = "RichnessIndex.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")

## Estimate_Richness -------------------------
### 多様性指数は任意で 
### シングルトンが無い場合、警告文が出るが、無視
alpha_df <- estimate_richness(PhyseqData)
colnames(alpha_df)

rownames(data.frame(sample_data(PhyseqData)))
rownames(alpha_df)

### 行数が同じであることを確認した上で、cbind()
alpha_df <- cbind(data.frame(sample_data(PhyseqData)), alpha_df) 


### NonChimera --------------------------------
library(tibble)

rownames(track)
colnames(track)
class(track)

track <- rownames_to_column(as.data.frame(track), var = "Sample.Name")

track |> 
    ggplot(aes(x = Sample.Name, y = nonchim,
               fill = Sample.Name)) +
    scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3",
                                 "#e7298a", "#66a61e", "#e6ab02",
                                 "#a6761d", "#666666", "#8dd3c7")) + 
    geom_bar(stat = "identity", width = 0.7) +
    ylab("Nonchim") +
    xlab("Sample.Name") +
    theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
          axis.text = element_text(size = 12, face = "bold", color = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
          legend.background = element_rect(fill = "gray90"),
          legend.key = element_rect(fill = "white", color = NA))

ggsave(filename = "nonchim.png", plot = last_plot(), 
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")
### Observed ASVs -----------------------------




### Observed ASVs
alpha_df |> 
    ggplot(aes(x = Sample.Name, y = Observed,
               fill = Sample.Name)) +
    scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3",
                                 "#e7298a", "#66a61e", "#e6ab02",
                                 "#a6761d", "#666666", "#8dd3c7")) + 
    geom_bar(stat = "identity", width = 0.7) +
    ylab("Observed ASVs") +
    xlab("Sample.Name") +
    theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
          axis.text = element_text(size = 12, face = "bold", color = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
          legend.background = element_rect(fill = "gray90"),
          legend.key = element_rect(fill = "white", color = NA))


ggsave(filename = "ObservedASVs.png", plot = last_plot(), 
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")



### 一旦断念 (250511)
# alpha_df  |> 
#     group_by(dps)  |> 
#     summarise(
#         Observed_mean = mean(Observed, na.rm = TRUE),
#         Observed_se = sd(Observed, na.rm = TRUE) / sqrt(n()),
#         Shannon_mean = mean(Shannon, na.rm = TRUE),
#         Shannon_se = sd(Shannon, na.rm = TRUE) / sqrt(n()),
#         Chao1_mean = mean(Chao1, na.rm = TRUE) ,
#         Chao1_se = sd(Chao1, na.rm = TRUE) / sqrt(n()),
#         ACE_mean = mean(ACE, na.rm = TRUE) ,
#         ACE_se = sd(ACE, na.rm = TRUE) / sqrt(n()),
#         se.ACE_mean = mean(se.ACE, na.rm = TRUE) ,
#         se.ACE_se = sd(se.ACE, na.rm = TRUE) / sqrt(n()),
#         Simpson_mean = mean(Simpson, na.rm = TRUE) ,
#         Simpson_se = sd(Simpson, na.rm = TRUE) / sqrt(n()),
#         InvSimpson_mean = mean(InvSimpson, na.rm = TRUE) ,
#         InvSimpson_se = sd(InvSimpson, na.rm = TRUE) / sqrt(n()),
#         Fisher_mean = mean(Fisher, na.rm = TRUE) ,
#         Fisher_se = sd(Fisher, na.rm = TRUE) / sqrt(n())) |> 
#     pivot_longer(
#         cols = c(Observed_mean, Shannon_mean, Chao1_mean, ACE_mean, se.ACE_mean),
#         names_to = "index",
#         values_to = "value"
#     )
#     ggplot(aes(x = dps, y = Shannon_mean)) +
#     geom_bar(stat = "identity", fill = "steelblue", width = 0.3) +
#     geom_errorbar(aes(ymin = Shannon_mean - Shannon_se,
#                       ymax = Shannon_mean + Shannon_se),
#                   width = 0.1) +
#     ylab("Shannon Diversity Index") +
#     xlab("dps") +
#     theme_minimal()

## Compare_RichnessIndex ---------------------

### Checking Normality&EquallyDispersed -------
library(car)

#### QQ Plots
qqnorm(alpha_df$Shannon[alpha_df$dps == 0])
qqline(alpha_df$Shannon[alpha_df$dps == 0])

#### 正規性の確認
by(alpha_df$Shannon, alpha_df$dps, shapiro.test)

#### 等分散性の確認
bartlett.test(Shannon ~ dps, data = alpha_df)
leveneTest(Shannon ~ dps, data = alpha_df)

#### → Sample数が3つと小さいため、NonParametricで実施する

### RichnessIndex StatisticalProcessing -------

compare_means(Shannon ~ dps, data = alpha_df,
              method = "wilcox.test", label = "p.format")
compare_means(Shannon ~ Fungicide.use, data = alpha_df,
              method = "wilcox.test", label = "p.format")

compare_means(Shannon ~ dps, data = alpha_df,
              method = "t.test", label = "p.format")
compare_means(Shannon ~ Fungicide.use, data = alpha_df,
              method = "t.test", label = "p.format")

compare_means(Shannon ~ dps, data = alpha_df,
              method = "anova", label = "p.format")
compare_means(Shannon ~ Fungicide.use, data = alpha_df,
              method = "anova", label = "p.format")

compare_means(Shannon ~ dps, data = alpha_df,
              method = "kruskal.test", label = "p.format")
compare_means(Shannon ~ Fungicide.use, data = alpha_df,
              method = "kruskal.test", label = "p.format")
# kruskal.test(Shannon ~ dps, data = alpha_df)




## Plots_Richness ----------------------------
### Compare_dps -------------------------------

# alpha_df |> 
#     ggboxplot(x = "dps", y = "Shannon", color = "dps",
#               palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#               add = "jitter") + 
#     geom_pwc(method = "wilcox_test", label = "{p.format}{p.signif}") +
#     theme(axis.title = element_text(size = 14, face = "bold", color = "black"),
#           axis.text = element_text(size = 12, face = "bold", color = "black"),
#           panel.grid.minor = element_blank(),
#           legend.text = element_text(size = 10, color = "black"),
#           legend.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5))


## 多様性指数の統計処理 
stat.test_dps <- compare_means(Shannon ~ dps, data = alpha_df, 
                               method = "wilcox.test", label = "p.format")
## y.position(p値表示の高さ)
stat.test_dps$y.position <- seq(
    from = max(alpha_df$Shannon, na.rm = TRUE) * 1.05,
    by = max(alpha_df$Shannon, na.rm = TRUE) * 0.05,
    length.out = nrow(stat.test_dps))

## ggplots
ggplot(alpha_df, aes(x = dps, y = Shannon)) +
    geom_boxplot() +
    theme(
        legend.position = "top", 
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray80"), 
        panel.grid.minor = element_line(color = "gray90")) + 
    geom_jitter(aes(color = Sample.Name), width = 0.08, size = 2.5, alpha = 0.6) +
    stat_pvalue_manual(stat.test_dps, label = "WilcoxTest {p.format} {p.signif}") +
    scale_color_manual(values = c("#E64B35",  
                                  "#4DBBD5",  
                                  "#00A087",  
                                  "#3C5488",  
                                  "#F39B7F",  
                                  "#8491B4",  
                                  "#91D1C2",  
                                  "#DC0000",  
                                  "#7E6148",  
                                  "#B09C85",  
                                  "#FFDC91",
                                  "#e7298a")) + 
    
    ggsave(filename = "Richness_Shannon_dps.png", plot = last_plot(),
           width = 2800, height = 2520, dpi = 300, units = "px",
           path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")


### Compare_Fungicide.use ---------------------

alpha_df |> 
    ggboxplot(x = "Fungicide.use", y = "Shannon", color = "Fungicide.use",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              add = "jitter") + 
    geom_pwc(method = "wilcox_test", label = "{p.format}{p.signif}") +
    theme(axis.title = element_text(size = 14, face = "bold", color = "black"),
          axis.text = element_text(size = 12, face = "bold", color = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5))


stat.test <- compare_means(Shannon ~ Fungicide.use, data = alpha_df,
                           method = "wilcox.test", label = "p.format")

## y.position(p値表示の高さ)
stat.test$y.position <- max(alpha_df$Shannon, na.rm = TRUE) * 1.05

## ggplots
ggplot(alpha_df, aes(x = Fungicide.use, y = Shannon)) +
    geom_boxplot() +
    geom_jitter(aes(color = Sample.Name), width = 0.1, size = 3, alpha = 0.6) +
    stat_pvalue_manual(stat.test, label = "【Wilcox.Test】 {p.format}{p.signif}",
                       size = 5, bracket.size = 0.4) +
    scale_color_manual(values = c("#E64B35",  
                                  "#4DBBD5",  
                                  "#00A087",  
                                  "#3C5488",  
                                  "#F39B7F",  
                                  "#8491B4",  
                                  "#91D1C2",  
                                  "#DC0000",  
                                  "#7E6148",  
                                  "#B09C85",  
                                  "#FFDC91")) +
    theme(
        legend.position = "top", 
        axis.title = element_text(size = 18, face = "bold", color = "black"),
        axis.text = element_text(size = 14, face = "bold", color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray80"), 
        panel.grid.minor = element_line(color = "gray90"))


ggsave(filename = "Richness_Shannon_Fungicide.use.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")

# Rarefaction Curve -------------------------

rare <- lapply(rarecurve(as(otu_table(PhyseqData), "matrix"),
                         step = 100, cex = 0.5, label = FALSE),
               function(x){
                   b <- as.data.frame(x)
                   b <- data.frame(ASV = b[,1], raw_read = rownames(b))
                   b$raw_read <- as.numeric(gsub("N", "",  b$raw_read))
                   return(b)
})
names(rare) <- rownames(sample_data(PhyseqData))
rare_df <- purrr::map_dfr(rare, function(x) return(data.frame(x)), .id = "sample")

ggplot(rare_df, aes(x = raw_read, y = ASV, color = sample)) +
    geom_line(linewidth = 1) + scale_color_igv() +
    xlab("Reads") + ylab("The number of ASV")


# Nonmetric Multidimensional Scaling --------

set.seed(1234)
ordinate(PhyseqData, "NMDS", "bray") |> 
plot_ordination(physeq = PhyseqData, color = "dps", shape = "dps") +
    stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill=dps)) +
    geom_point(size = 2) +
    scale_color_startrek() +
    scale_fill_startrek() +
    xlab("Axis 1") +
    ylab("Axis 2") + 
    ggtitle("NMDS")

# MicrobiotaProcessMethods ------------------
## alpha diversity analysis ------------------

library("MicrobiotaProcess")
load("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Output/PhyloseqData_Bacteria.RData")

mpdata <- 
    PhyseqData |> as.MPSE()

# 以下は同じ意味
## x <- x |> log()
## x %<>% log()

mpdata %<>% mp_rrarefy()

mpdata %<>% 
    mp_cal_rarecurve(
        .abundance = RareAbundance,
        chunks = 400
    )
mpdata |> print(width=180)

# p1 <- mpdata |>
#     mp_plot_rarecurve(
#         .rare = RareAbundanceRarecurve,
#         .alpha = Observe,
#     )
# p1



## calculate alpha index and visualization --------
library(ggplot2)
library(MicrobiotaProcess)
library(magrittr)

mpdata %<>% 
    mp_cal_alpha(.abundance=RareAbundance) # α多様性指数の算出
mpdata



f1 <- mpdata |> 
    mp_plot_alpha(
        .group=dps, 
        .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
    ) +
    scale_fill_manual(values=c("#00A087FF", "#3C5488FF","#cab2d6"), guide="none") +
    scale_color_manual(values=c("#00A087FF", "#3C5488FF","#cab2d6"), guide="none")

f2 <- mpdata |>
    mp_plot_alpha(
        .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
    )

f1 / f2


# 関数の内部表示
# MicrobiotaProcess::mp_plot_alpha
# showMethods(mp_plot_alpha)
# getMethod("mp_plot_alpha", signature = "MPSE")
# ggsignif::geom_signif()で、検定手法を定めている


## The visualization of taxonomy abundance --------

mpdata %<>%
    mp_cal_abundance( # for each samples
        .abundance = RareAbundance
    ) |>
    mp_cal_abundance( # for each groups 
        .abundance = RareAbundance,
        .group=dps
    )
mpdata

## MicrobiotaProcess_ModelData ---------------
# ランダムサンプリングされたRareAbundanceを用いていることに注意！
p1 <- mpdata |> 
    mp_plot_abundance(
        .abundance=RareAbundance,
        .group=dps, 
        taxa.class = Phylum, 
        topn = 20,
        relative = TRUE
    )
# visualize the abundance (rarefied) of top 20 phyla for each sample.
p2 <- mpdata |>
    mp_plot_abundance(
        .abundance=RareAbundance,
        .group=dps,
        taxa.class = Phylum,
        topn = 20,
        relative = FALSE
    )
p1 / p2


## The relative abundance and abundance of phyla of all samples --------

h1 <- mpdata |>
    mp_plot_abundance(
        .abundance = RareAbundance,
        .group = dps,
        taxa.class = Phylum,
        relative = TRUE,
        topn = 20,
        geom = 'heatmap',
        features.dist = 'euclidean',
        features.hclust = 'average',
        sample.dist = 'bray',
        sample.hclust = 'average'
    ) 
h2 <- mpdata |>
    mp_plot_abundance(
        .abundance = RareAbundance,
        .group = dps,
        taxa.class = Phylum,
        relative = FALSE,
        topn = 20,
        geom = 'heatmap',
        features.dist = 'euclidean',
        features.hclust = 'average',
        sample.dist = 'bray',
        sample.hclust = 'average'
    )
# the character (scale or theme) of figure can be adjusted by set_scale_theme
# refer to the mp_plot_dist
aplot::plot_list(gglist=list(h1, h2), tag_levels="A")


# visualize the relative abundance of top 20 phyla for each .group (time)
p3 <- mpdata |>
    mp_plot_abundance(
        .abundance=RareAbundance, 
        .group=dps,
        taxa.class = Phylum,
        topn = 20,
        plot.group = TRUE
    )

# visualize the abundance of top 20 phyla for each .group (time)
p4 <- mpdata |>
    mp_plot_abundance(
        .abundance=RareAbundance,
        .group= dps,
        taxa.class = Phylum,
        topn = 20,
        relative = FALSE,
        plot.group = TRUE
    )
p3 / p4

dev.off()

## Beta diversity analysis -------------------

# standardization
# mp_decostand wraps the decostand of vegan, which provides 
# many standardization methods for community ecology.
# default is hellinger, then the abundance processed will
# be stored to the assays slot. 
mpdata %<>% 
    mp_decostand(.abundance=Abundance)
mpdata

# calculate the distance between the samples.
# the distance will be generated a nested tibble and added to the
# colData slot.
mpdata %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
mpdata

detach("package:phyloseq", unload=TRUE)


p5 <- mpdata |> mp_plot_dist(.distmethod = bray)
p5



# when .group is provided, the dot heatmap plot with group information will be return.
p6 <- mpdata |> MicrobiotaProcess::mp_plot_dist(.distmethod = bray, .group = dps)
# The scale or theme of dot heatmap plot can be adjusted using set_scale_theme function.
p6 |> set_scale_theme(
    x = scale_fill_manual(
        values=c("orange","deepskyblue","#1b9e77"), 
        guide = guide_legend(
            keywidth = 1, 
            keyheight = 0.5, 
            title.theme = element_text(size=8),
            label.theme = element_text(size=6))), 
    aes_var = dps) |> # specific the name of variable 
    set_scale_theme(
        x = scale_color_gradient(
            guide = guide_legend(keywidth = 0.5, keyheight = 0.5)),
        aes_var = bray) |> 
    set_scale_theme(
        x = scale_size_continuous(
            range = c(0.1, 3),
            guide = guide_legend(keywidth = 0.5, keyheight = 0.5)),
        aes_var = bray)

dev.off()

# when .group is provided and group.test is TRUE, the comparison of different groups will be returned
p7 <- mpdata |> mp_plot_dist(.distmethod = bray, .group = dps, group.test=TRUE, textsize=2)
p7



## The PCoA analysis -------------------------
# fungicide useでも比較してみること

library(ggplot2)

mpdata %<>% 
    mp_cal_pcoa(.abundance = hellinger, distmethod="bray")

# The dimensions of ordination analysis will be added the colData slot (default).
mpdata

# We also can perform adonis or anosim to check whether it is significant to the dissimilarities of groups.
mpdata %<>%
    mp_adonis(.abundance = hellinger, .formula = ~dps, distmethod = "bray", permutations = 9999, action = "add")
mpdata %>% mp_extract_internal_attr(name = adonis)


p8 <- mpdata %>%
    mp_plot_ord(
        .ord = pcoa, 
        .group = dps, 
        .color = dps, 
        .size = 1.2,
        .alpha = 1,
        ellipse = TRUE,
        show.legend = FALSE # don't display the legend of stat_ellipse
    ) +
    scale_fill_manual(values=c("#00A087FF", "#3C5488FF", "orange")) +
    scale_color_manual(values=c("#00A087FF", "#3C5488FF", "orange")) 

# The size of point also can be mapped to other variables such as Observe, or Shannon 
# Then the alpha diversity and beta diversity will be displayed simultaneously.
p9 <- mpdata %>% 
    mp_plot_ord(
        .ord = pcoa, 
        .group = dps, 
        .color = dps, 
        .size = Shannon,
        ellipse = TRUE,
        show.legend = FALSE # don't display the legend of stat_ellipse 
    ) +
    scale_fill_manual(
        values = c("#00A087FF", "#3C5488FF","orange"), 
        guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
    ) +
    scale_color_manual(
        values=c("#00A087FF", "#3C5488FF", "orange"),
        guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
    ) +
    scale_size_continuous(
        range=c(0.5, 3),
        guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
    )


p8 + p9


## Hierarchical cluster analysis -------------

mpdata %<>%
    mp_cal_clust(
        .abundance = hellinger, 
        distmethod = "bray",
        hclustmethod = "average", # (UPGAE)
        action = "add" # action is used to control which result will be returned
    )

mpdata



# if action = 'add', the result of hierarchical cluster will be added to the MPSE object
# mp_extract_internal_attr can extract it. It is a treedata object, so it can be visualized
# by ggtree.
sample.clust <- mpdata |>
    mp_extract_internal_attr(name='SampleClust')
sample.clust



library(ggtree)
p10 <- ggtree(sample.clust) + 
    geom_tippoint(aes(color=dps)) +
    geom_tiplab(as_ylab = TRUE) +
    ggplot2::scale_x_continuous(expand=c(0, 0.01))
p10


dev.off()

library(ggtreeExtra)
library(ggplot2)

phyla.tb <- mpdata |> 
    mp_extract_abundance(taxa.class = Phylum, topn = 30)

# The abundance of each samples is nested, it can be flatted using the unnest of tidyr.
phyla.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) |> dplyr::rename(Phyla="label")
phyla.tb


p11 <- p10 + 
    geom_fruit(
        data = phyla.tb,
        geom = geom_col,
        mapping = aes(x = RelRareAbundanceBySample, 
                      y = Sample.Name,
                      fill = label
        ),
        orientation = "y",
        #offset = 0.4,
        pwidth = 3, 
        axis.params = list(axis = "x", 
                           title = "The relative abundance of phyla (%)",
                           title.size = 4,
                           text.size = 2, 
                           vjust = 1),
        grid.params = list()
    )
p11



## Biomarker discovery -----------------------
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(MicrobiotaProcess)
library(tidytree)
library(ggstar)
library(forcats)

mpdata |> print(width=150)


### Transform Abundance value -----------------
load("~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/芝草/NGS_consignment/Novogene/Data/250331/RData/SaveObjects/mpdata.RData")

# 相対存在量Countに変換
mpdata %<>%
    mp_cal_abundance(.abundance = Abundance, action = "add",
                     relative = TRUE, force = TRUE)

# 相対存在量が0.1% > でFiltering
# min.propでFilteringを適用される基準を設ける
# Sample数 = 9 + min.prop = 0.12 → 1.08以上のSampleで0.1%以上の存在量が確認されたASVsは、除去されない
# → つまり、2 sampleでmin.abun = 0.1を満たすASVsしか残らない
# min.prop * sample数 = 1未満の場合、少なくとも 1 sampleでmin.prop = 0.1が満たされていれば良い
# つまり、以下のコードと同じ動作をする
# → transform_sample_counts(PhyseqData, function(x) 100*(x / sum(x)))
# → filter_taxa(PhyseqData, function(x) max(x) > 0.1, TRUE)

# 平均ではない点に注意


# 以下のFilteringは、Phyloseq Packageよりも柔軟性が高くて便利
## transform_sample_counts(PhyseqData, function(x) 100*(x / sum(x)))
## filter_taxa(PhyseqData, function(x) max(x) > 0.1, TRUE)
mpdata %<>%
    mp_filter_taxa(.abundance = RelAbundanceBySample,
                   min.prop = 0.1, # 1 sampleでも0.1%>で存在していれば、除去されない
                   min.abun = 0.1 # 最低相対存在量 0.1% > Filtering  
    )

# csvで保存
write.csv(mpdata@assays@data@listData$RelAbundanceBySample,
          file = "~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/芝草/NGS_consignment/Novogene/Data/NGS_Analysis/physeqData_csv/mpdata@assays@data@listData$RelAbundanceBySample.csv")



### Differentially abundant taxa --------------

# 処理区でsampleが最低でも5つ必要!!!! → 断念!!!
mpdata %<>%
    mp_diff_analysis(
       .abundance = Abundance, # 存在量を指定 (相対存在量 or ランダムサンプリング存在量 or 絶対存在量 or 正規化存在量)
       .group = dps, # 比較対象
       first.test.alpha = 0.1, # 有意性の一次判定に用いるp値のしきい値
       second.test.alpha = 0.05, action = "add", relative = TRUE # 多重検定補正後のしきい値
    )


taxa.tree %<>%
    mp_diff_analysis(
        .abundance = RelRareAbundanceBySample,
        .group = dps,
        first.test.alpha = 0.01
    )
# The result is stored to the taxatree or otutree slot, you can use mp_extract_tree to extract the specific slot.
taxa.tree <- mpdata |> 
    mp_extract_tree(type="taxatree")
taxa.tree

# And the result tibble of different analysis can also be extracted 
# with tidytree (>=0.3.5)
taxa.tree |>
    select(label, nodeClass, LDAupper,LDAmean,
           LDAlower,Sign_time, pvalue, fdr) |> 
    dplyr::filter(!is.na(fdr))



### Visualized by ggtree and ggtreeExt --------
# Since taxa.tree is treedata object, it can be visualized by ggtree and ggtreeExtra
p1 <- ggtree(
    taxa.tree,
    layout="radial",
    size = 0.3
) +
    geom_point(
        data = td_filter(!isTip),
        fill="white",
        size=1,
        shape=21
    )


# display the high light of phylum clade.
p2 <- p1 +
    geom_hilight(
        data = td_filter(nodeClass == "Phylum"),
        mapping = aes(node = node, fill = label)
    )


# display the relative abundance of features(OTU)
p3 <- p2 +
    ggnewscale::new_scale("fill") +
    geom_fruit(
        data = td_unnest(RareAbundanceBySample),
        geom = geom_star,
        mapping = aes(
            x = fct_reorder(Sample.Name, dps, .fun=min),
            size = RelRareAbundanceBySample,
            fill = dps,
            subset = RelRareAbundanceBySample > 0
        ),
        starshape = 13,
        starstroke = 0.25,
        offset = 0.04,
        pwidth = 0.8,
        grid.params = list(linetype=2)
    ) +
    scale_size_continuous(
        name="Relative Abundance (%)",
        range = c(.5, 3)
    ) +
    scale_fill_manual(values=c("#1B9E77", "#D95F02"))


# display the tip labels of taxa tree
p4 <- p3 + geom_tiplab(size=2, offset=7.2)


# display the LDA of significant OTU.
p5 <- p4 +
    ggnewscale::new_scale("fill") +
    geom_fruit(
        geom = geom_col,
        mapping = aes(
            x = LDAmean,
            fill = Sign_time,
            subset = !is.na(LDAmean)
        ),
        orientation = "y",
        offset = 0.3,
        pwidth = 0.5,
        axis.params = list(axis = "x",
                           title = "Log10(LDA)",
                           title.height = 0.01,
                           title.size = 2,
                           text.size = 1.8,
                           vjust = 1),
        grid.params = list(linetype = 2)
    )


# display the significant (FDR) taxonomy after kruskal.test (default)
p6 <- p5 +
    ggnewscale::new_scale("size") +
    geom_point(
        data=td_filter(!is.na(Sign_time)),
        mapping = aes(size = -log10(fdr),
                      fill = Sign_time,
        ),
        shape = 21,
    ) +
    scale_size_continuous(range=c(1, 3)) +
    scale_fill_manual(values=c("#1B9E77", "#D95F02"))

p6 + theme(
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 9),
    
# save(p1,p2,p3,p4,p5,p6,
#      file = "~/Library/CloudStorage/GoogleDrive-saito2022@patholab-meiji.jp/My Drive/芝草/NGS_consignment/Novogene/Data/250331/RData")



### Visualized by MicrobiotaProcessMet --------
# Because the released `ggnewscale` modified the internal new aesthetics name,
# The following code is to obtain the new aesthetics name according to version of
# `ggnewscale`


flag <- packageVersion("ggnewscale") >= "0.5.0"
new.fill <- ifelse(flag, "fill_ggnewscale_1", "fill_new_new")
new.fill2 <- ifelse(flag , "fill_ggnewscale_2", "fill_new")


p <- mpdata |>
    mp_plot_diff_res(
        group.abun = TRUE,
        pwidth.abun=0.1
    ) 

# if version of `ggnewscale` is >= 0.5.0, you can also use p$ggnewscale to view the renamed scales.
p <- p  +
    scale_fill_manual(values=c("deepskyblue", "orange")) +
    scale_fill_manual(
        aesthetics = new.fill2, # The fill aes was renamed to `new.fill` for the abundance dotplot layer
        values = c("deepskyblue", "orange")
    ) +
    scale_fill_manual(
        aesthetics = new.fill, # The fill aes for hight light layer of tree was renamed to `new.fill2`
        values = c("#E41A1C", "#377EB8", "#4DAF4A",
                   "#984EA3", "#FF7F00", "#FFFF33",
                   "#A65628", "#F781BF", "#999999"
                   )
        )

p

### mp_plot_diff_cladogram --------------------

f <- mpdata |>
    mp_plot_diff_cladogram(
        label.size = 2.5,
        hilight.alpha = .3,
        bg.tree.size = .5,
        bg.point.size = 2,
        bg.point.stroke = .25
    ) +
    scale_fill_diff_cladogram( # set the color of different group.
        values = c('deepskyblue', 'orange')
    ) +
    scale_size_continuous(range = c(1, 4))

f



### mp_plot_diff_boxplot ----------------------

f.box <- mpdata |> 
    mp_plot_diff_boxplot(.group = dps) |> 
    set_diff_boxplot_color(
        values = c("deepskyblue", "orange"),
        guide = guide_legend(title=NULL)
        )

f.bar <- mpdata |> 
    mp_plot_diff_boxplot(
        taxa.class = c(Genus, OTU), # select the taxonomy level to display
        group.abun = TRUE, # display the mean abundance of each group
        removeUnknown = TRUE, # whether mask the unknown taxa.
    ) |> 
    set_diff_boxplot_color(
        values = c("deepskyblue", "orange"),
        guide = guide_legend(title=NULL)
        )

aplot::plot_list(f.box, f.bar)


### mp_plot_diff_manhattan --------------------

f.mahattan <- mpdata |> 
    mp_plot_diff_manhattan(
        .group = Sign_time,
        .y = fdr,
        .size = 2.4,
        taxa.class = c('OTU', 'Genus'),
        anno.taxa.class = Phylum
        )

f.mahattan





# DESeq2 ------------------------------------
## PhyseqObjects to DESeq2 -----------------------

rm(list = ls())
library(phyloseq)
library(DESeq2)
library(dplyr)


load("~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/RData/phyloseq_Bacteria/Output/PhyloseqData_Bacteria.RData")


### Family Level Filtering ------------------------------

rank_names(PhyseqData)

PhyseqData_Family <- PhyseqData  |> 
    subset_taxa(Kingdom == "Bacteria") |> # Kingdomを"Bacteria"でSubset
    filter_taxa(function(x) mean(x) > 100, TRUE) # Read数が約100000であり、0.1%のReads数 


#### DESeq2 Fungicide.use ---------

dps_dds_Family = phyloseq_to_deseq2(PhyseqData_Family, ~ `Fungicide.use`) 
## DESeq()のParameterは、最適な引数を設定すること
dps_dds_Family = DESeq(dps_dds_Family, test="Wald", fitType="parametric") 

res_Family <- results(dps_dds_Family, cooksCutoff = FALSE)
res_Family <- cbind(as(res_Family, "data.frame"),
                    as(tax_table(PhyseqData)[rownames(res_Family), ], "matrix"))

write.csv(cbind(as(res_Family, "data.frame"), as(tax_table(PhyseqData)[rownames(res_Family), ], "matrix")),
          file = "~/Documents/RStudio/Novogene/250503/export_csv/res_Family.csv")


alpha = 0.01

#### 0.01よりもpadjが小さいASVsをFiltering
sigtab_Family = res_Family[which(res_Family$padj < alpha), ] 

dim(sigtab_Family)


#### Results_ggplot ----------------------------
library(ggplot2)
library(scales)

##### sigtab_Plots ------------------------------

#### Phylum order
x = tapply(sigtab_Family$log2FoldChange, sigtab_Family$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_Family$Phylum = factor(as.character(sigtab_Family$Phylum), levels=names(x))

#### Genus order
x = tapply(sigtab_Family$log2FoldChange, sigtab_Family$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_Family$Genus = factor(as.character(sigtab_Family$Genus), levels=names(x))



#### padjをlog10にする
sigtab_Family$log10value <- -log10(sigtab_Family$padj)

#### log2FoldChange < 0 → Fungicide.use(Bac4~9)において、存在量が増加したことを示す
sigtab_Family$Sign <- ifelse(sigtab_Family$log2FoldChange < 0, "Enriched", "Depleted")

#### 行名をASV列として、追加
sigtab_Family$ASV <- rownames(sigtab_Family)

sigtab_Family <- sigtab_Family |> 
    select(ASV, log10value, log2FoldChange, Sign, everything())

#### sigtab_Familyのcsv保存
write.csv(sigtab_Family, file = "~/Documents/RStudio/Novogene/250503/export_csv/sigtab_Family.csv",
          row.names = FALSE)

unique(sigtab_Family$Family)
unique(sigtab_Family$Phylum)

# unique(sigtab_Family$Family)に基づいて、colorパレットを設定
colors_21 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
               "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
               "#999999", "#66C2A5", "#FC8D62", "#8DA0CB",
               "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", 
               "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3","#E7298A")



ggplot(sigtab_Family, aes(x = Family, y = abs(log2FoldChange), color = Phylum, shape = Sign)) +
    geom_point(mapping = aes(size = 4, alpha = 0.6)) +
    geom_text(aes(label = ASV), vjust = -1, size = 2) + 
    scale_color_manual(values = colors_21) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    scale_shape_manual(values = c("Enriched" = 16, "Depleted" = 17)) +  
    ylab("log2FoldChange") +
    xlab("Family") +
    theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
          axis.text = element_text(size = 12, face = "bold", color = "black"),
          strip.text = element_text(size = 10, face = "bold.italic", color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = c(x=0.65, y=0.9)) +
    labs(shape = "Enriched Circle or Depleted Triangle in FungicideUse") +
    guides(
        color = "none",
        size = "none",
        alpha = "none",
        shape = guide_legend(override.aes = list(size = 5),
                             title.theme = element_text(face = "bold", size = 20)
    ))


ggsave(filename = "DESeq2_Family_log2FoldChange_ColorPhylum_sigtab_Plots.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")


##### res_Plots ---------------------------------

res_Family$log10value <- -log10(res_Family$padj)
res_Family$Sign <- ifelse(res_Family$log2FoldChange < 0, "Enriched", "Depleted")
res_Family$ASV <- rownames(res_Family)

colnames(res_Family)
unique(res_Family$Family)
unique(res_Family$Genus)

### Family数に応じて、Colorを設定
colors_40 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
               "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
               "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
               "#666666", "#8C564B", "#C49C94", "#D62728", "#9467BD", "#2CA02C", "#FFBB78", "#98DF8A",
               "#AEC7E8", "#FF9896", "#C5B0D5", "#F7B6D2", "#C7E9C0", "#9ECAE1", "#FDD0A2", "#DADAEB")

colors_29 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
               "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
               "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
               "#666666", "#8C564B", "#C49C94", "#D62728", "#9467BD")

### log10 Value Plots in Family Colors
ggplot(res_Family, aes(x = Family, y = abs(log10value), color = Family, shape = Sign)) +
    geom_point(mapping = aes(size = 4, alpha = 0.6)) +
    geom_text(aes(label = ASV), vjust = -1, size = 2.5) + 
    geom_hline(mapping = aes(yintercept = -log10(0.05)), linetype = 2) +
    scale_color_manual(values = colors_40) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    scale_shape_manual(values = c("Enriched" = 16, "Depleted" = 17)) +  
    ylab("-log10 Value") +
    xlab("Family") +
    theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
          axis.text = element_text(size = 10,color = "black"),
          strip.text = element_text(size = 10, face = "bold.italic", color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = c(x=0.8, y=0.9)) +
    labs(shape = "Enriched Circle or Depleted Triangle") +
    guides(
        color = "none",
        size = "none",
        alpha = "none",
        shape = guide_legend(override.aes = list(size = 4),
                             title.theme = element_text(face = "bold", size = 16))
    )


ggsave(filename = "DESeq2_Family_log10Value_ColorFamily_res_Family_Plots.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")


#### Legend Plots
library(ggpubr)
library(cowplot)

ggdraw(
    get_legend(
        ggplot(res_Family, aes(x = Family, y = abs(log10value), color = Family)) +
            scale_color_manual(values = colors_40) +
            theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
            ylab("-log10 Value") +
            xlab("Family") +
            theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
                  axis.text = element_text(size = 10,color = "black"),
                  strip.text = element_text(size = 10, face = "bold.italic", color = "black"),
                  panel.grid.minor = element_blank(),
                  legend.position = "top") +
            labs(shape = "Enriched Circle or Depleted Triangle"), return_all = TRUE))

cowplot::ggdraw(
    cowplot::get_legend(
        ggplot(res_Family, aes(x = Family, y = abs(log10value), color = Family)) +
            geom_point(aes(size = 4, alpha = 0.6)) +
            scale_shape_manual(values = c("Enriched" = 16, "Depleted" = 17)) +
            labs(shape = "Enriched Circle or Depleted Triangle") +
            theme_minimal() +
            theme(legend.position = "top")
        )
    )


### log10 Value Plots in Genus Colors
ggplot(res_Family, aes(x = Family, y = abs(log10value), color = Genus, shape = Sign)) +
    geom_point(mapping = aes(size = 4, alpha = 0.6)) +
    geom_text(aes(label = ASV), vjust = -1, size = 2.5) + 
    geom_hline(mapping = aes(yintercept = -log10(0.05)), linetype = 2) +
    scale_color_manual(values = colors_40) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    scale_shape_manual(values = c("Enriched" = 16, "Depleted" = 17)) +  
    ylab("-log10 Value") +
    xlab("Family") +
    theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
          axis.text = element_text(size = 10,color = "black"),
          strip.text = element_text(size = 10, face = "bold.italic", color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = c(x=0.7, y=0.9)) +
    labs(shape = "Enriched Circle or Depleted Triangle") +
    guides(
        color = "none",
        size = "none",
        alpha = "none",
        shape = guide_legend(override.aes = list(size = 4),
                             title.theme = element_text(face = "bold", size = 16))
    )

ggsave(filename = "DESeq2_Family_log10Value_ColorGenus_res_Family_Plots.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")



### log2FoldChange in Genus Colors
ggplot(res_Family, aes(x = Family, y = abs(log2FoldChange), color = Genus, shape = Sign)) +
    geom_point(mapping = aes(size = 4, alpha = 0.6)) +
    geom_text(aes(label = ASV), vjust = -1, size = 2.8) + 
    geom_hline(mapping = aes(yintercept = -log10(0.05)), linetype = 2) + ## 有意水準の設定
    scale_color_manual(values = colors_40) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    scale_shape_manual(values = c("Enriched" = 16, "Depleted" = 17)) +  
    ylab("log2FoldChange Value") +
    xlab("Family") +
    theme(axis.title = element_text(size = 20, face = "bold", color = "black"),
          axis.text = element_text(size = 11,color = "black"),
          strip.text = element_text(size = 11, face = "bold.italic", color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = c(x=0.7, y=0.9)) +
    labs(shape = "Enriched Circle or Depleted Triangle") +
    guides(
        color = "none",
        size = "none",
        alpha = "none",
        shape = guide_legend(override.aes = list(size = 4),
                             title.theme = element_text(face = "bold", size = 16))
    )

ggsave(filename = "DESeq2_Family_log2Value_ColorGenus_res_Family_Plots.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")



### log2FoldChange in Family Colors
ggplot(res_Family, aes(x = Family, y = abs(log10value), color = Family, shape = Sign)) +
    geom_point(mapping = aes(size = 4, alpha = 0.6)) +
    geom_text(aes(label = ASV), vjust = -1, size = 2.5) + 
    geom_hline(mapping = aes(yintercept = -log10(0.05)), linetype = 2) +
    scale_color_manual(values = colors_40) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    scale_shape_manual(values = c("Positive" = 16, "Negative" = 17)) +  
    ylab("-log10 Value") +
    xlab("Family") +
    theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
          axis.text = element_text(size = 10,color = "black"),
          strip.text = element_text(size = 10, face = "bold.italic", color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = "none")

ggsave(filename = "DESeq2_Family_res_Family_Plots.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")





### Genus Level Filtering ---------------------

PhyseqData_Genus <- PhyseqData  |> 
    subset_taxa(Kingdom == "Bacteria") |> 
    filter_taxa(function(x) mean(x) > 100, TRUE) # Read数が約100000であり、0.1%のReads数 

#### DESeq2 Fungicide.use ---------

dps_dds_Genus = phyloseq_to_deseq2(PhyseqData_Genus, ~ `Fungicide.use`) 
## DESeq()のParameterは、最適な引数を設定すること
dps_dds_Genus = DESeq(dps_dds_Genus, test="Wald", fitType="parametric") 


res_Genus <- results(dps_dds_Genus, cooksCutoff = FALSE)
res_Genus <- cbind(as(res_Genus, "data.frame"),
                   as(tax_table(PhyseqData)[rownames(res_Genus), ], "matrix"))

write.csv(cbind(as(res_Genus, "data.frame"), as(tax_table(PhyseqData)[rownames(res_Genus), ], "matrix")),
          file = "~/Documents/RStudio/Novogene/250503/export_csv/res_Genus.csv")


alpha = 0.01

#### 0.01よりもpadjが小さいASVsをFiltering
sigtab_Genus = res_Genus[which(res_Genus$padj < alpha), ] 

#### Taxa Tableもcbind
sigtab_Genus = cbind(as(sigtab_Genus, "data.frame"), as(tax_table(PhyseqData)[rownames(sigtab_Genus), ], "matrix"))
head(sigtab_Genus)
dim(sigtab_Genus)
write.csv(sigtab_Genus, file = "~/Documents/RStudio/Novogene/250503/export_csv/sigtab_Genus.csv")


#### Results_ggplot ----------------------------
library(ggplot2)
library(scales)

##### sigtab_Plots ------------------------------
#### Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

#### Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

sigtab$log10value <- -log10(sigtab$padj)
sigtab$Sign <- ifelse(sigtab$log2FoldChange < 0, "Negative", "Positive")
sigtab$ASV <- rownames(sigtab)

colors_21 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
               "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
               "#999999", "#66C2A5", "#FC8D62", "#8DA0CB",
               "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", 
               "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3","#E7298A")

unique(sigtab$Family)

ggplot(sigtab, aes(x = Family, y = abs(sigtab$log2FoldChange), color = Phylum, shape = Sign)) +
    geom_point() +
    geom_text(aes(label = ASV), vjust = -1, size = 3) + 
    scale_color_manual(values = colors_21) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    scale_shape_manual(values = c("Positive" = 16, "Negative" = 17)) +  
    ylab("log2FoldChange") +
    xlab("Family") +
    theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
          axis.text = element_text(size = 12, face = "bold", color = "black"),
          strip.text = element_text(size = 10, face = "bold.italic", color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = "none")

ggsave(filename = "Relative_abundance_Phylum_SampleName.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")


##### res_Plots ---------------------------------

res$log10value <- -log10(res$padj)
res$Sign <- ifelse(res$log2FoldChange < 0, "Depleted", "Enriched")
res$ASV <- rownames(res)

colnames(res)
unique(res$Family)
unique(res$Genus)

### Family数に応じて、Colorを設定
colors_40 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
               "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
               "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
               "#666666", "#8C564B", "#C49C94", "#D62728", "#9467BD", "#2CA02C", "#FFBB78", "#98DF8A",
               "#AEC7E8", "#FF9896", "#C5B0D5", "#F7B6D2", "#C7E9C0", "#9ECAE1", "#FDD0A2", "#DADAEB")

colors_29 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
               "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
               "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
               "#666666", "#8C564B", "#C49C94", "#D62728", "#9467BD")

### log10 Value Plots in Family Colors
ggplot(res, aes(x = Family, y = abs(log10value), color = Family, shape = Sign)) +
    geom_point(mapping = aes(size = 4, alpha = 0.6)) +
    geom_text(aes(label = ASV), vjust = -1, size = 2.5) + 
    geom_hline(mapping = aes(yintercept = -log10(0.05)), linetype = 2) +
    scale_color_manual(values = colors_40) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    scale_shape_manual(values = c("Enriched" = 16, "Depleted" = 17)) +  
    ylab("-log10 Value") +
    xlab("Family") +
    theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
          axis.text = element_text(size = 10,color = "black"),
          strip.text = element_text(size = 10, face = "bold.italic", color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = c(x=0.8, y=0.9)) +
    labs(shape = "Enriched Circle or Depleted Triangle") +
    guides(
        color = "none",
        size = "none",
        alpha = "none",
        shape = guide_legend(override.aes = list(size = 4),
                             title.theme = element_text(face = "bold", size = 16))
    )


ggsave(filename = "DESeq2_Family_log10Value_ColorFamily_res_Plots.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")


#### Legend Plots
library(ggpubr)
library(cowplot)

ggdraw(
    get_legend(
        ggplot(res, aes(x = Family, y = abs(log10value), color = Family)) +
            scale_color_manual(values = colors_40) +
            theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
            ylab("-log10 Value") +
            xlab("Family") +
            theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
                  axis.text = element_text(size = 10,color = "black"),
                  strip.text = element_text(size = 10, face = "bold.italic", color = "black"),
                  panel.grid.minor = element_blank(),
                  legend.position = "top") +
            labs(shape = "Enriched Circle or Depleted Triangle"), return_all = TRUE))

cowplot::ggdraw(
    cowplot::get_legend(
        ggplot(res, aes(x = Family, y = abs(log10value), color = Family)) +
            geom_point(aes(size = 4, alpha = 0.6)) +
            scale_shape_manual(values = c("Enriched" = 16, "Depleted" = 17)) +
            labs(shape = "Enriched Circle or Depleted Triangle") +
            theme_minimal() +
            theme(legend.position = "top")
    )
)


### log10 Value Plots in Genus Colors
ggplot(res, aes(x = Family, y = abs(log10value), color = Genus, shape = Sign)) +
    geom_point(mapping = aes(size = 4, alpha = 0.6)) +
    geom_text(aes(label = ASV), vjust = -1, size = 2.5) + 
    geom_hline(mapping = aes(yintercept = -log10(0.05)), linetype = 2) +
    scale_color_manual(values = colors_40) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    scale_shape_manual(values = c("Enriched" = 16, "Depleted" = 17)) +  
    ylab("-log10 Value") +
    xlab("Family") +
    theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
          axis.text = element_text(size = 10,color = "black"),
          strip.text = element_text(size = 10, face = "bold.italic", color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = c(x=0.7, y=0.9)) +
    labs(shape = "Enriched Circle or Depleted Triangle") +
    guides(
        color = "none",
        size = "none",
        alpha = "none",
        shape = guide_legend(override.aes = list(size = 4),
                             title.theme = element_text(face = "bold", size = 16))
    )

ggsave(filename = "DESeq2_Family_log10Value_ColorGenus_res_Plots.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")



### log2FoldChange in Genus Colors
ggplot(res, aes(x = Family, y = abs(log2FoldChange), color = Genus, shape = Sign)) +
    geom_point(mapping = aes(size = 4, alpha = 0.6)) +
    geom_text(aes(label = ASV), vjust = -1, size = 2.8) + 
    geom_hline(mapping = aes(yintercept = -log10(0.05)), linetype = 2) + ## 有意水準の設定
    scale_color_manual(values = colors_40) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    scale_shape_manual(values = c("Enriched" = 16, "Depleted" = 17)) +  
    ylab("log2FoldChange Value") +
    xlab("Family") +
    theme(axis.title = element_text(size = 20, face = "bold", color = "black"),
          axis.text = element_text(size = 11,color = "black"),
          strip.text = element_text(size = 11, face = "bold.italic", color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = c(x=0.7, y=0.9)) +
    labs(shape = "Enriched Circle or Depleted Triangle") +
    guides(
        color = "none",
        size = "none",
        alpha = "none",
        shape = guide_legend(override.aes = list(size = 4),
                             title.theme = element_text(face = "bold", size = 16))
    )

ggsave(filename = "DESeq2_Family_log2Value_ColorGenus_res_Plots.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")



### log2FoldChange in Family Colors
ggplot(res, aes(x = Family, y = abs(log10value), color = Family, shape = Sign)) +
    geom_point(mapping = aes(size = 4, alpha = 0.6)) +
    geom_text(aes(label = ASV), vjust = -1, size = 2.5) + 
    geom_hline(mapping = aes(yintercept = -log10(0.05)), linetype = 2) +
    scale_color_manual(values = colors_40) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    scale_shape_manual(values = c("Positive" = 16, "Negative" = 17)) +  
    ylab("-log10 Value") +
    xlab("Family") +
    theme(axis.title = element_text(size = 18, face = "bold", color = "black"),
          axis.text = element_text(size = 10,color = "black"),
          strip.text = element_text(size = 10, face = "bold.italic", color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = "none")

ggsave(filename = "DESeq2_Family_res_Plots.png", plot = last_plot(),
       width = 2800, height = 2520, dpi = 300, units = "px",
       path = "~/Documents/RStudio/Novogene/250503/NGS_analysis_microbiome/png")


## cladogram ---------------------------------
library(phyloseq)
library(DESeq2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(tidytree)
library(dplyr)
library(ggstar)
library(ggplot2)
library(treeio)


tree_family <- phy_tree(PhyseqData_Family) |> as.treedata()

res_df <- as.data.frame(res)

tax_df <- as.data.frame(tax_table(PhyseqData_Family))
tax_df$ASV <- rownames(tax_df)

sigtab <- left_join(res_df, tax_df, by = "ASV")  |> 
    filter(padj < 0.05) |> 
    select(ASV,everything())


p1 <- ggtree(tree_family, layout = "fan", open.angle = 20, size = 0.4)

p1.1 <- p1 + 
    geom_fruit(
        data = sigtab,
        geom = geom_star,
        mapping = aes(y=ASV, fill=Phylum, size=log2FoldChange),
        position="identity",
        starstroke=0.2
    ) + 
    geom_hilight(
        data = td_filter(sigtab == "Phylum"),
        mapping = aes(node = node, fill = label)
    )
    scale_size_continuous(
        range=c(1, 3), # the range of size.
        guide=guide_legend(
            keywidth=0.5, 
            keyheight=0.5,
            override.aes=list(starshape=15),
            order=2
        )
    ) +
    scale_fill_manual(
        values=c("#8DD3C7",  # turquoise
                 "#FFFFB3",  # pale yellow
                 "#BEBADA",  # lavender
                 "#FB8072",  # salmon
                 "#80B1D3",  # sky blue
                 "#FDB462",  # orange
                 "#B3DE69",  # lime green
                 "#FCCDE5",  # pink
                 "#D9D9D9",  # gray
                 "#BC80BD",  # purple
                 "#CCEBC5"),   # mint
        guide="none" 
    ) + 
    scale_starshape_manual(
        values=c(1, 15),
        guide=guide_legend(
            keywidth=0.5,
            keyheight=0.5,
            order=1
        )
    )

p1.1


# Other -------------------------------------

rank_names(PhyseqData)
ps1 = prune_taxa((prev0 > prevalenceThreshold), PhyseqData)
ps2 = subset_taxa(ps1, Phylum %in% names(keepPhyla))

# ps2は、フィルタリングされてきたphyseqObjects
ID <- rank_names(PhyseqData)
for (i in ID) {
    print(length(get_taxa_unique(ps2, taxonomic.rank = i)))
}

PhyseqData_RA <- transform_sample_counts(PhyseqData, function(x) 100* x / sum(x))
PhyseqData_log10 <- transform_sample_counts(PhyseqData, log)


