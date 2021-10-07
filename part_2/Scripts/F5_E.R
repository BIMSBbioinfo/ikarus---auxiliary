#!/opt/R/4.0/bin/Rscript
args = commandArgs(trailingOnly = TRUE)

### INFO: Intersection of Ikarus's signature with Chitars's fused genes
### DATE: 13.04.2021
### AUTHOR: Artem Baranovskii

require(tidyverse)
require(patchwork)
require(broom)

### setup
set.seed(42)
p.dat <- file.path(args[1], "Data/")
p.gs <- file.path(args[1], "Data/signatures/")
p.out <- file.path(args[1], "Figures/")


#### Read in chitars database info on gene fusions
l.fus <- purrr::map2(c("TCGA_ChiTaRS_combined_fusion_information_on_hg19.txt"),
                     c(F),
                     ~ read_tsv(file.path(p.dat, "chitars", .x), col_names = .y)) %>%
  `names<-`(c("Chitars/TCGA"))

#### Load expression matched random gene sets
t.bsd <- purrr::map(list.files(file.path(p.dat, "bootstrapped"), pattern = "bootstrapped", full.names = T),
                    ~ read_tsv(.x,
                               col_names = T,
                               col_types = cols(
           gene_name = col_character(),
           avg.lgexp = col_double(),
           variance = col_double(),
           gset = col_character()
         )) %>%
           mutate(gset = ifelse(str_detect(gset, "Random"), paste0(gset, "_", str_split_fixed(.x, "\\.", Inf)[, 3]), gset))
         ) %>%
  purrr::reduce(., rbind)

#### Test intersections with ChiTaRS
l.tmp <- purrr::map(split(t.bsd, t.bsd$gset), 
           ~ c(intersect(.x$gene_name, str_split_fixed(unique(l.fus$`Chitars/TCGA`$X5), "_", 2)[, 1]) %>% length(), 
               intersect(.x$gene_name, str_split_fixed(unique(l.fus$`Chitars/TCGA`$X9), "_", 2)[, 1]) %>% length()))
t.bsd.int <- tibble(gset = names(l.tmp), 
                    `5p_int` = purrr::map(l.tmp, ~ .x[1]) %>% purrr::reduce(., c), 
                    `3p_int` = purrr::map(l.tmp, ~ .x[2]) %>% purrr::reduce(., c))

#### Morph data to plot histograms
t.bsd.int.tall <- t.bsd.int %>% mutate(gset = str_split_fixed(gset, "_", 2)[, 1]) %>% pivot_longer(., contains("int"), values_to = "ints", names_to = "grp")
t.bsd.int.tall.sum <- t.bsd.int.tall %>% filter(gset == "Tumor")



#### Render in pdf
pdf(file = file.path(p.out, "F5_E.pdf"), onefile = T, width = 7, height = 4.2)
(
  ggplot(data = t.bsd.int.tall %>% filter(gset == "Random", grp == "5p_int"), aes(x = ints)) + 
  geom_density(aes(y = ..count..)) +
  geom_histogram(binwidth = 1, alpha = 0.3, aes(fill = gset)) + 
  geom_vline(xintercept = t.bsd.int.tall.sum[t.bsd.int.tall.sum$grp == "5p_int", ]$ints, 
             color = "#FF0000") +
  #scale_size_manual(values = c(3, 2)) +
  scale_fill_manual(values = c("grey30", "orangered")) + 
  #lims(y = c(0, 6.5)) +
  labs(x = "Intersection size", y = "Frequency", fill = "Gene list", subtitle = "5´ fusion") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 5/4) 
  
  |
    
  ggplot(data = t.bsd.int.tall %>% filter(gset == "Random", grp == "3p_int"), aes(x = ints)) + 
  geom_density(aes(y = ..count..)) +
  geom_histogram(binwidth = 1, alpha = 0.3, aes(fill = gset)) + 
  geom_vline(xintercept = t.bsd.int.tall.sum[t.bsd.int.tall.sum$grp == "3p_int", ]$ints, 
             color = "#FF0000") +
  #scale_size_manual(values = c(3, 2)) +
  scale_fill_manual(values = c("grey30", "orangered")) + 
  #lims(y = c(0, 6.5)) +
  labs(x = "Intersection size", y = "Frequency", fill = "Gene list", subtitle = "3´ fusion") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 5/4)
  ) + 
  plot_layout(guides = "collect")

dev.off()