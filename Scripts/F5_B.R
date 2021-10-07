#!/opt/R/4.0/bin/Rscript
args = commandArgs(trailingOnly = TRUE)

### INFO: Performance of Ikarus's tumor signature on the human protein atlas's (THPA) pathology data (per cancer)
### DATE: 12.09.2021
### AUTHOR: Artem Baranovskii

require(tidyverse)

### setup
set.seed(42)
p.thp <- file.path(args[1], "Data/THPA/")
p.gs <- file.path(args[1], "Data/signatures/")
p.out <- file.path(args[1], "Figures/")


#### Read Ikarus signatures into a gene list
l.gs <- map(read_lines(list.files(p.gs, pattern = "\\.gmt", full.names = T)),
            ~ str_split_fixed(.x, "\t", Inf)) %>%
  map(., ~ .x[which(str_length(.x) > 0)])
names(l.gs) <- map(l.gs, ~ .x[[1]])
l.gs <- map(l.gs, ~ .x[3:length(.x)])

#### Process THPA pathology data
t.pat <- data.table::fread(file.path(p.thp, "pathology.tsv"), data.table = F)

rbind(t.pat %>% filter(`Gene name` %in% l.gs$Tumor) %>% mutate(t.gl = "Cancer list"),
      t.pat %>% mutate(t.gl = "All genes")) %>% 
  mutate(t.gl = factor(t.gl, levels = c("Cancer list", "All genes"))) -> t.pat

### Prognostic in at least one cancer?
t.pat.tmp <- as_tibble(t.pat) %>% 
  filter(!is.na(`prognostic - favorable`) | !is.na(`prognostic - unfavorable`)) %>% 
  filter(`prognostic - favorable` < 0.01 | `prognostic - unfavorable` < 0.01) %>% 
  dplyr::select(1:3, "prognostic - favorable", "prognostic - unfavorable", t.gl) %>% 
  pivot_longer(., contains("prognostic"), names_to = "grp", values_to = "km_pval") %>% 
  filter(!is.na(km_pval)) %>% 
  mutate(grp = str_split_fixed(grp, " - ", 2)[, 2])
### Summed
t.pat.tmp %>% group_by(t.gl, grp) %>% dplyr::count() %>% ungroup() %>% group_by(t.gl) %>% mutate(perc = (n / sum(n)) * 100) -> t.pat.sum

### Per cancer
t.pat.tmp %>% group_by(t.gl, Cancer, grp) %>% dplyr::count() %>% ungroup() %>% group_by(t.gl, Cancer) %>% mutate(perc = (n / sum(n)) * 100) %>% ungroup() -> t.pat.pcr
## Append NAs so the barplots look fancy
t.pat.pcr <- rbind(t.pat.pcr, 
                   tribble(~t.gl, ~Cancer, ~grp, ~n, ~perc, 
                           "Cancer list", "colorectal cancer", "unfavorable", NA, NA, 
                           "Cancer list", "glioma", "favorable", NA, NA,
                           "Cancer list", "head and neck cancer", "favorable", NA, NA,
                           "Cancer list", "lung cancer", "favorable", NA, NA,
                           "Cancer list", "melanoma", "favorable", NA, NA,
                           "Cancer list", "ovarian cancer", "unfavorable", NA, NA))
ordr <- t.pat.pcr %>% filter(t.gl == "Cancer list", grp == "unfavorable") %>% arrange(desc(n)) %>% pull(Cancer) %>% unique()


pdf(file = file.path(p.out, "F5_B.pdf"), onefile = T, width = 5, height = 5)
  print(
    ggplot(data = t.pat.pcr %>% filter(t.gl == "Cancer list") %>% mutate(Cancer = factor(Cancer, levels = ordr)),
        aes(x = Cancer, y = n, fill = grp)) + 
    geom_bar(stat = "identity", position = "dodge", color = "grey30", width = 0.7) + 
    scale_fill_manual(values = c("#aaaaaa", "#FF0000")) + 
    scale_y_continuous(breaks = seq(0, 35, by = 5)) +
    labs(y = "Number of prognostic genes", x = "Cancer type", fill = "Prognosis") +
    theme_bw(base_size = 14) +
    theme(aspect.ratio = 4/5, 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = c(.75, .75))
  )
dev.off()