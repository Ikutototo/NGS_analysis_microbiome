
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




## Sample#1 ----------------------------------

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
    
    
