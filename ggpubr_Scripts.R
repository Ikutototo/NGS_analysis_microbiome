# tidyplots ---------------------------------

library(tidyplots)

study |> 
    tidyplot(x = treatment, y = score, color = treatment) |> 
    add_boxplot() |> 
    add_test_pvalue(ref.group = 1)



# ggpubr ------------------------------------
# https://rpkgs.datanovia.com/ggpubr/index.html

library(ggpubr)

ggboxplot(ToothGrowth, x = "dose", y = "len", color = "dose",
          palette =c("#00AFBB", "#E7B800", "#FC4E07"),
          add = "jitter", shape = "dose") + 
    geom_pwc(method = "t_test", label = "{p.format}{p.signif}")





## Sample1 ----------------------------------

sample <- read.csv("~/Documents/RStudio/Novogene/250503/import_csv/sample_test.csv")

unique(sample$NAME)

library(rstatix)
library(dplyr)

res.aov <- sample  |> 
    dplyr::filter(stringr::str_detect(NAME,"CK")) |>  
    anova_test(dv = X.Mock, between = NAME)

tukey.res <- sample  |> 
    dplyr::filter(stringr::str_detect(NAME,"CK")) |>  
    tukey_hsd(X.Mock ~ NAME)

tukey.res <- tukey.res  |> 
    filter(p.adj < 0.05)


ggboxplot(sample, x = "NAME", y = "X.Mock", add = "jitter") +
    stat_pvalue_manual(tukey.res, label = "p.adj.signif", 
                       y.position = tukey.res$y.position)



custom_p_format <- function(p) {
    rstatix::p_format(p, accuracy = 0.0001, digits = 3, leading.zero = FALSE)
}

(plots <- sample |> 
    dplyr::filter(stringr::str_detect(NAME,"CK")) |> 
    ggboxplot(x = "NAME", y = "X.Mock", color = "NAME", add = "jitter") + 
    geom_pwc(method = "tukey_hsd", label = "p.signif", step.increase = 0.1))


head(sample, 50)
    

## Sample2 ----------------------------------

library(ggpubr)
alpha_df <- read.csv("~/Documents/RStudio/Novogene/250503/export_csv/alpha_df.csv")
colnames(alpha_df)


### Compare_RichnessIndex ---------------------
compare_means(Shannon ~ dps, data = alpha_df,
              method = "wilcox.test", label = "p.format")

compare_means(Shannon ~ dps, data = alpha_df,
              method = "t.test", label = "p.format")

compare_means(Shannon ~ dps, data = alpha_df,
              method = "anova", label = "p.format")

compare_means(Shannon ~ dps, data = alpha_df,
              method = "kruskal.test", label = "p.format")



### Plots_Richness ----------------------------
#### Compare_dps -------------------------------

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


#### Compare_Fungicide.use ---------------------

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




