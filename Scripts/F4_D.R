#!/opt/R/4.0/bin/Rscript
args = commandArgs(trailingOnly = TRUE)

### INFO: Plot Ikarus's intersections with msigdb signatures
### DATE: 13.04.2020
### AUTHOR: Artem Baranovskii

require(tidyverse)
require(msigdbr)

### setup
set.seed(42)
p.gs <- file.path(args[1], "Data/signatures/")
p.out <- file.path(args[1], "Figures/")


#### Read Ikarus signatures into a gene list
l.gs <- map(read_lines(list.files(p.gs, pattern = "\\.gmt", full.names = T)),
            ~ str_split_fixed(.x, "\t", Inf)) %>%
  map(., ~ .x[which(str_length(.x) > 0)])
names(l.gs) <- map(l.gs, ~ .x[[1]])
l.gs <- map(l.gs, ~ .x[3:length(.x)])


### Hallmarks
hs_db <- msigdbr::msigdbr(species = "Homo sapiens")
hs_db.hmrk <- hs_db %>% filter(grepl("hallmark", gs_name, ignore.case = T))
l.hmrk <- purrr::map(split(hs_db.hmrk, hs_db.hmrk$gs_name), ~ .x %>% pull(gene_symbol))
### ALL msigdb
l.msdb <- purrr::map(split(hs_db, hs_db$gs_name), ~ .x %>% pull(gene_symbol))

### Prepare table for plotting
tibble(gset = names(l.msdb), 
       gset.size = purrr::map(l.msdb, ~ length(.x)) %>% simplify(.type = "numeric"),
       intn = purrr::map(l.msdb, ~ intersect(l.gs$Tumor, .x) %>% length()) %>% simplify(.type = "numeric") %>% unname(), 
       unin = purrr::map(l.msdb, ~ union(l.gs$Tumor, .x) %>% length()) %>% simplify(.type = "numeric") %>% unname(), 
       intn_genes = purrr::map(l.msdb, ~ intersect(l.gs$Tumor, .x)) %>% unname()) -> t.msdb
t.msdb <- t.msdb %>% mutate(intn_genes = as.character(intn_genes))

### render in pdf
pdf(file = file.path(p.out, "F4_D.pdf"), onefile = T, width = 7, height = 4.5)
print(
  t.msdb %>%
  dplyr::filter(intn > 0) %>%
  mutate(var.ann = ifelse(intn > 26, gset, ""),
         var.col = ifelse(intn > 26, "red", "")) %>%
ggplot(data = ., aes(x = gset.size, y = intn / unin)) +
  geom_point(aes(color = var.col)) +
  scale_color_manual(values = c("grey30", "red")) +
  ggrepel::geom_text_repel(aes(label = stringi::stri_wrap(var.ann, 7, 0.0)), force = 97, seed = 1488, size = 2.5, max.overlaps = 27) +
  #lims(y = c(0, 45)) +
  labs(x = "MsigDB gene set size", y = "Jaccard similarity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 4/5,
        legend.position = "none")
)
dev.off()