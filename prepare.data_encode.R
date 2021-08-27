#!/opt/R/4.0/bin/Rscript
### INFO: Process ENCODE expression data in Ikarus pipeline fashion
### DATE: 13.04.2020
### AUTHOR: Artem Baranovskii

library(tidyverse)
library(magrittr)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)


w.path <- "/local/Projects/AAkalin_Ikarus/Data/ENCODE/"

##
t.ann <- read_tsv(file.path(w.path, "metadata.tsv"))
### --- Keep only GRCh38
t.ann %>% filter(`File assembly` == "GRCh38", `Output type` == "gene quantifications") -> t.ann
### --- In the downloaded data there were tables processed differently -> remove them from the downstream processing
#### #load full data
t.ann %$% map(`File accession`, 
              ~ fread(file.path(w.path, paste0(.x, ".tsv")), data.table = F)) -> l.dat
#### #identify faulty tables
l.dat %>% map(., ~ any(grepl("gene_id", colnames(.x)))) %>% as.logical() %>% which() -> v.keep
### --- Generate expression table (omitting faulty tables)
map2(l.dat[v.keep], 
     t.ann$`File accession`[v.keep], 
     ~ as_tibble(.x) %>% 
       dplyr::select("gene_id", "TPM") %>% 
       filter(grepl("ENS", gene_id)) %>% 
       dplyr::rename(!!sym(.y) := TPM)) %>% 
  reduce(., inner_join, by = "gene_id") -> t.exp
### --- log TPMs 
t.exp %<>% mutate(across(where(is.numeric), ~ log2(. + 1)))
### --- append gene names
#### ---- prepare HGNC tab
library(org.Hs.eg.db)
uniKeys <- keys(org.Hs.eg.db, keytype = "ENTREZID")
cols <- c("SYMBOL", "ENSEMBL")
t.gns <- AnnotationDbi::select(org.Hs.eg.db, keys = uniKeys, columns = cols, keytype = "ENTREZID") %>% as_tibble()
#### ----
t.exp %>% 
  mutate(gene_id = str_split_fixed(gene_id, "\\.", 2)[, 1]) %>% 
  inner_join(t.gns[, c("SYMBOL", "ENSEMBL")], ., by = c("ENSEMBL" = "gene_id")) -> t.exp
### --- Write tmp annotated table
write_tsv(t.exp, file.path(w.path, "ikarus", paste0(Sys.Date(), "tmp_exp.lg2TPM.tsv")))



## Correct signature sizes if needed
### --- 
regen.sig <- function(t.exp, s.name, p.keep = 0.9) {
  set.seed(42)
  #
  v.gns <- t.exp$SYMBOL
  #
  t.sig <- read_csv("/local/abarano/ads/ikarus-aucell/out/Tumor_gene_list.csv", col_names = F)
  n.sig <- read_csv("/local/abarano/ads/ikarus-aucell/out/Normal_gene_list.csv", col_names = F)
  #
  g.mis <- list()
  g.mis[["t"]] <- t.sig %>% dplyr::filter(!X1 %in% v.gns) %>% pull(X1)
  g.mis[["n"]] <- n.sig %>% dplyr::filter(!X1 %in% v.gns) %>% pull(X1)
  sigs.upd <- list()
  sigs.upd[["t"]] <- c(t.sig$X1[!t.sig$X1 %in% g.mis[["t"]]], 
                       sample(g.mis[["t"]], round((length(t.sig$X1) - length(g.mis[["t"]])) / p.keep, 0) - (length(t.sig$X1) - length(g.mis[["t"]]))))
  sigs.upd[["n"]] <- c(n.sig$X1[!n.sig$X1 %in% g.mis[["n"]]], 
                       sample(g.mis[["n"]], round((length(n.sig$X1) - length(g.mis[["n"]])) / p.keep, 0) - (length(n.sig$X1) - length(g.mis[["n"]]))))
  #
  dummy <- purrr::map2(sigs.upd, 
                       c("Tumor", "Normal"), 
                       ~ write_csv(tibble(X1 = .x), paste0("/local/abarano/ads/ikarus-aucell/out/", .y, "_gene_list_", s.name, ".csv"), col_names = F))
}
regen.sig(t.exp, "ENCODE", p.keep = 0.83)



# OUTPUT 
# table, gene names. annotation
m.exp <- t(as.matrix(t.exp %>% dplyr::select(-SYMBOL, -ENSEMBL)))
colnames(m.exp) <- 0:(dim(m.exp)[2] - 1)
rownames(m.exp) <- 0:(dim(m.exp)[1] - 1)

t.gns <- tibble(X1 = 0:(dim(m.exp)[2] - 1), X2 = t.exp$SYMBOL)
t.ann.out <- t.ann %>% 
  dplyr::filter(`File accession` %in% colnames(t.exp)) %>% 
  dplyr::select(raw = `Biosample term name`, tier_0 = `Biosample type`) %>% 
  mutate(tier_1 = "", tier_2 = "", tier_3 = "")

## write
write.table(m.exp, sep = ",", file = file.path(w.path, "ikarus", "counts.csv"))
write.table(t.gns, sep = ",", file = file.path(w.path, "ikarus", "genes_symbol.csv"), col.names = F, row.names = F)
write.table(t.ann.out, sep = ",", file = file.path(w.path, "ikarus", "cell_labels.csv"), row.names = F)

