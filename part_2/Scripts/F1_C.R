#!/opt/R/4.0/bin/Rscript
args = commandArgs(trailingOnly = TRUE)

### INFO: Heatmap for the figure F1C
### DATE: 10.03.2021
### AUTHOR: Artem Baranovskii

require(tidyverse)
require(data.table)
require(ComplexHeatmap)
require(circlize)
require(tictoc)

### setup
set.seed(42)
dset <- "tirosh17_headneck"
p.wd <- file.path(args[1], "Data/sc_cancer/", dset)        # <---------- Change paths here to locate other scRNAseq datasets
p.gs <- file.path(args[1], "Data/signatures/")             # <---------- Change paths here to read other gene lists
p.out <- file.path(args[1], "Figures/")


#### Read signatures into a gene list
l.gs <- map(read_lines(list.files(p.gs, pattern = "\\.gmt", full.names = T)),
            ~ str_split_fixed(.x, "\t", Inf)) %>% 
  map(., ~ .x[which(str_length(.x) > 0)])
names(l.gs) <- map(l.gs, ~ .x[[1]])
l.gs <- map(l.gs, ~ .x[3:length(.x)])


#### Expression data and annotation
gnms <- read_csv(file.path(p.wd, "genes_symbol.csv"), col_names = F)
etab <- data.table::fread(file.path(p.wd, "counts_proc.csv"), skip = 1)
etab <- etab[, -1] # remove dummy variable
colnames(etab) <- gnms[[ifelse(dim(gnms)[2] == 1, 1, 2)]]
annot <- data.table::fread(file.path(p.wd, "cell_labels.csv"),
                           data.table = F)
annot <- annot %>% mutate(dummy = 1:dim(annot)[1]) # dummy var for cell subsetting
annot <- annot %>% mutate(raw = str_remove(raw, "-")) %>% 
  mutate(raw = paste0(raw, "_", tier_0))
# subset cells with reasonable representability in a dataset OR keep all of them, if the dataset of good quality
c.sel <- names(table(annot$raw))


#### print heatmap(s) into multi-page pdf
## open a connection to a pdf file where the plots will be saved
tic() # track time
pdf(file = paste0(p.out, "/", "F1_C.pdf"),
    onefile = T,
    width = 11,
    height = 14)
## iterate over gene signatures in the list
for (i in names(l.gs)[2]) {                                # <---------- keep only tumor signature
  # subset jan's signature genes and compute average expresion per cell type
  ss <- intersect(l.gs[[i]], colnames(etab))
  etab.ss <- etab[, ..ss]
  t.avge <- purrr::map(c.sel, 
                       ~ etab.ss[annot[annot$raw == .x, ]$dummy, ] %>% as.matrix() %>% matrixStats::colMeans2()) %>% 
    purrr::reduce(., rbind)
  t.avge <- apply(t.avge, 2, scale)
  t.avge <- t.avge[, colSums(is.na(t.avge)) == 0]
  t.avge <- t(t.avge) %>% `colnames<-`(c.sel) #%>% `rownames<-`(ss)
  ## Prepare HM annotation data, some minor corrections
  hm <- Heatmap(t.avge)
  an.col <- annot[annot$raw %in% colnames(t.avge), c("raw", "tier_0", "tier_1", "tier_2")] %>% 
    distinct(raw, .keep_all = T) %>% 
    mutate(tier_1 = ifelse(raw %in% c("granulocytes", "Langerhans cells"), "Immune", tier_1)) %>% 
    `[`(match(colnames(t.avge), .$raw), ) %>% 
    mutate(across(.fns = ~ ifelse(. == "", "NA", .)))
  ## some workaround 
  colnames(t.avge) <- str_split_fixed(c.sel, "_", 2)[, 1]
  ## build annotation
  ha <- HeatmapAnnotation(
    `Malignancy` = an.col$tier_0,
    `Cell type` = an.col$tier_1, 
    `Cell lineage` = an.col$tier_2, 
    col = list(`Malignancy` = c("Tumor" = "orangered", "Normal" = "#a1cae2"),
               `Cell type` = c("NA" = "#ffffff", 
                               "Endothelial" = "#f14668", "Epithelial" = "#ffd880", "Fibroblast" = "#6b011f", "Immune" = "steelblue3", 
                               "Other" = "darkolivegreen3", "Pericyte" = "aquamarine3"),
               `Cell lineage` = c("NA" = "#ffffff", 
                                  "Lymphoid" = "#a2d0c1", "Myeloid" = "#f8dc81")),
    gp = gpar(col = "black"), 
    border = F, 
    annotation_name_side = "left"
  )
  
  ## color function for z-scores
  col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  ## print is necesssary for iteration, will not work otherwise
  print(
    Heatmap(t.avge, name = "z-score", column_title = paste0(i, "; n genes = ", dim(t.avge)[1]),
            #rect_gp = gpar(col = "white", lwd = 0.3), 
            col = col_fun, top_annotation = ha,
            column_names_gp = gpar(fontsize = 12),
            column_names_max_height = max_text_width(
              colnames(t.avge), 
              gp = gpar(fontsize = 12)
            )
    )
  )
}
## close sonncetion to the .pdf
dev.off()
toc() # track time


