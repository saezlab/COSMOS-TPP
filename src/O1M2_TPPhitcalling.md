---
title: "O1M2: TPP hit calling and scoring"
output:
  html_document:
    keep_md: yes
    toc: yes
    theme: united
    code_folding: hide  
---

# General settings



```r
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  include = TRUE,
  cache = TRUE,
  cache.lazy = FALSE,
  eval = TRUE
  # fig.width = 6*(1+sqrt(5))/2, fig.height = 6 #golden ratio
)

knitr::opts_knit$set(root.dir = "C:/Users/burtsche/Documents/COSMOS-TPP_paper")
```

## Load packages


```r
library(tidyverse)
library("reshape2")
library(OmnipathR)
#library(xlsx)
library(viridis)
library(ggrepel)
options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(visNetwork)
library(knitr)
library(UpSetR)
library(TPP2D)
library(cowplot)
library(ReactomePA)

mutate <- dplyr::mutate
select <- dplyr::select
group_by <- dplyr::group_by
```

## Functions


```r
GSE_RPA <- function(geneList, universe) {
  geneList <- mapIds(org.Hs.eg.db, geneList, "ENTREZID", "SYMBOL")
  universe <- mapIds(org.Hs.eg.db, universe,'ENTREZID','SYMBOL')
  pathway_enrichment <- enrichPathway(
    gene = geneList,
    organism = "human",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    universe = universe,
    minGSSize = 5,
    maxGSSize = 500,
    readable = TRUE
  )
}
```



```r
values_effect <- c(
  "stabilized" = "darkcyan",
  "destabilized" = "aquamarine3",
  "increased abundance" = "maroon",
  "decreased abundance" = "salmon"
)
```

## Load data


```r
TPP2D_UWB <- read.table("data/2DTPP_UWB_Ola_2021-03-24.txt", row.names = 1, sep = "\t")
phosphoproteomics <- read_tsv("data/LIMMA_results_PhosphoFollowUp_Phospho.txt")


load("data/KSN.RData")
load("data/210802_limma_UWB1.289_initialdataset_correctedphospho.RData")
load("data/220510_viper_footprints.RData")
```

Aggregate and filter TPP UWB data to considered condition and quality criteria.


```r
proteins_merged <- TPP2D_UWB %>%
  filter(condition == "Olaparib 24h") %>%
  # convert temperature to numeric
  mutate(
    qupm = as.numeric(qupm),
    temperature = as.numeric(gsub("C", "", temperature)),
    gene_name = paste(representative, gene_name, sep = ".")
  ) %>%
  # filter proteins which were found with qupm > 1 at least once
  group_by(gene_name) %>%
  filter(any(qupm > 1)) %>%
  ungroup() %>%
  # extract relevant FCs
  dplyr::select(gene_name, matches("rel_fc_protein_.+uM$"), temperature, hit_call) %>%
  melt(id.vars = c("gene_name", "temperature", "hit_call")) %>%
  mutate(concentration = as.numeric(sub("^.*_(.+)uM$", "\\1", variable))) %>%
  # convert FC to log2fc
  mutate(log2_value = log2(value))
```

# TPP hit calling

The permutation takes one night!!! Don't run this chunk per default (eval=F).


```r
TPPstatistical_input <- list()
for (file in c(
  "TPPraw_uwb_ola24_48.csv", "TPPraw_uwb_ola24_49.csv", "TPPraw_uwb_ola24_50.csv",
  "TPPraw_uwb_ola24_51.csv", "TPPraw_uwb_ola24_52.csv", "TPPraw_uwb_ola24_53.csv"
)) {
  tmp_file <- read_csv(paste0("data/raw/", file))
  tmp_file <- tmp_file %>%
    mutate(
      gene_name = `Gene name`, # paste(`Gene name`, Accession, sep = "."),
      protein_id = Accession
    ) %>% # c("1": as.character(nrow(tmp_file))),
    # protein_id = as.character(protein_id))  %>%
    dplyr::select(protein_id, gene_name, QUPM, matches("sumionarea"), matches("fold change")) %>%
    setNames(c(
      "protein_id", "gene_name", "qupm", "signal_sum_126", "signal_sum_127L", "signal_sum_127H", "signal_sum_128L", "signal_sum_128H",
      "signal_sum_129L", "signal_sum_129H", "signal_sum_130L", "signal_sum_130H", "signal_sum_131L", "rel_fc_126", "rel_fc_127L", "rel_fc_127H",
      "rel_fc_128L", "rel_fc_128H", "rel_fc_129L", "rel_fc_129H", "rel_fc_130L", "rel_fc_130H", "rel_fc_131L"
    )) %>%
    mutate_at(.funs = ~ . / rel_fc_128H, .vars = c("rel_fc_126", "rel_fc_127L", "rel_fc_127H", "rel_fc_128L", "rel_fc_128H"))

  TPPstatistical_input <- append(TPPstatistical_input, list(tmp_file))
}


data("config_tab")
data("raw_dat_list")

names(TPPstatistical_input) <- c("exp1", "exp2", "exp3", "exp4", "exp5", "exp6")

configuration <- config_tab %>%
  mutate(
    Compound = "Olaparib",
    Temperature = c(42.1, 44.1, 46.2, 48.1, 50.4, 51.9, 54, 56.1, 58.2, 60.1, 62.4, 63.9),
    RefCol = c("128H", "131L", "128H", "131L", "128H", "131L", "128H", "131L", "128H", "131L", "128H", "131L"),
    `126` = c("10", "-", "10", "-", "10", "-", "10", "-", "10", "-", "10", "-"),
    `127L` = c("4", "-", "4", "-", "4", "-", "4", "-", "4", "-", "4", "-"),
    `127H` = c("1", "-", "1", "-", "1", "-", "1", "-", "1", "-", "1", "-"),
    `128L` = c("0.4", "-", "0.4", "-", "0.4", "-", "0.4", "-", "0.4", "-", "0.4", "-"),
    `128H` = c("0", "-", "0", "-", "0", "-", "0", "-", "0", "-", "0", "-"),
    `129L` = c("-", "10", "-", "10", "-", "10", "-", "10", "-", "10", "-", "10"),
    `129H` = c("-", "4", "-", "4", "-", "4", "-", "4", "-", "4", "-", "4"),
    `130L` = c("-", "1", "-", "1", "-", "1", "-", "1", "-", "1", "-", "1"),
    `130H` = c("-", "0.4", "-", "0.4", "-", "0.4", "-", "0.4", "-", "0.4", "-", "0.4"),
    `131L` = c("-", "0", "-", "0", "-", "0", "-", "0", "-", "0", "-", "0")
  )

import_df <- import2dDataset(
  configTable = configuration,
  data = TPPstatistical_input,
  idVar = "protein_id",
  intensityStr = "signal_sum_",
  fcStr = "rel_fc_",
  nonZeroCols = "qupm",
  geneNameVar = "gene_name",
  addCol = NULL,
  qualColName = "qupm",
  naStrs = c("NA", "n/d", "NaN"),
  concFactor = 1e6,
  medianNormalizeFC = TRUE,
  filterContaminants = TRUE
)

recomp_sig_df <- recomputeSignalFromRatios(import_df)
preproc_df <- resolveAmbiguousProteinNames(recomp_sig_df)
model_params_df <- getModelParamsDf(
  df = preproc_df
)

fstat_df <- computeFStatFromParams(model_params_df)

set.seed(12, kind = "L'Ecuyer-CMRG")
null_model <- bootstrapNullAlternativeModel(
  df = preproc_df,
  params_df = model_params_df,
  B = 20
)

fdr_tab <- getFDR(
  df_out = fstat_df,
  df_null = null_model
)

hits <- findHits(
  fdr_df = fdr_tab,
  alpha = 0.1
)
TPPFstatistic_hits <- hits %>%
  dplyr::select(clustername, nObs, F_statistic, FDR, slopeH1, detected_effectH1, rssH0, rssH1)

plot2dTppFit(recomp_sig_df, "CHEK2",
  model_type = "H1"
)

save(TPPFstatistic_hits, file = "results/211007_TPPfstathits_UWB24h")
```

Load the results to avoid recalculation


```r
load("data/211007_TPPfstathits_UWB24h.RData")

TPPFstatistic_hits_adap <- TPPFstatistic_hits %>%
  mutate(
    mode = ifelse(detected_effectH1 == "stability" & slopeH1 > 0, "stabilized", NA),
    mode = ifelse(detected_effectH1 == "stability" & slopeH1 < 0, "destabilized", mode),
    mode = ifelse(detected_effectH1 != "stability" & slopeH1 < 0, "decreased abundance", mode),
    mode = ifelse(detected_effectH1 != "stability" & slopeH1 > 0, "increased abundance", mode),
    score = sign(slopeH1) * sqrt(rssH0 - rssH1),
    p_adj = FDR,
    pval = NA,
    hit_call = "TPPFstatistic",
    gene_name = clustername,
    representative = mapIds(org.Hs.eg.db, clustername, "UNIPROT", "SYMBOL")
  ) %>%
  dplyr::select(representative, gene_name, hit_call, mode, score, pval, p_adj)
```

```r
ggplot(TPPFstatistic_hits_adap, aes(x= mode, fill=mode)) +
  geom_bar() +
  scale_fill_manual(values = values_effect) +
  cowplot::theme_cowplot() *
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](O1M2_TPPhitcalling_files/figure-html/number_TPPhits-1.png)<!-- -->


# TPP hit evaluation Evaluation

## Pathways

Check pathway enrichment for de/stabilized hits


```r
pathway_check_input <- TPPFstatistic_hits_adap %>%
  filter(p_adj < 0.05)


pathways_universe <- proteins_merged %>%
  separate(gene_name, into = c("representative", "gene_name"), sep = "\\.")%>% 
  distinct(gene_name)

pathway_enrichment <- pathway_check_input %>%
  mutate(enrichment = list(as.data.frame(GSE_RPA(
    geneList = gene_name,
    universe = pathways_universe$gene_name
  ))))

pathway_enrichment_unnest <- pathway_enrichment %>%
  distinct(gene_name, enrichment) %>%
  distinct(enrichment) %>% 
  unnest()

pathway_enrichment_significant <- pathway_enrichment_unnest %>%
  filter(qvalue < 0.05)
```

Plot top 10 pathways for stability score pathway enrichment.


```r
# extract top10 pathways for plot
top_TPP_pathways_TPPFstatistic <- pathway_enrichment_significant %>%
  as.data.frame() %>%
  top_n(10, wt = -log10(qvalue))

p_pathway_UWB1.289_TPPFstatistic <- ggplot(top_TPP_pathways_TPPFstatistic, aes(y = -log10(qvalue), x = reorder(Description, -log10(qvalue)))) +
  geom_bar(stat = "identity", position = position_dodge(), fill = "indianred1") +
  geom_hline(yintercept = -log10(0.05)) +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 10)) +
  labs(x = "", title = "Significant pathways\nTPP F-statistic hits")
p_pathway_UWB1.289_TPPFstatistic
```

![](O1M2_TPPhitcalling_files/figure-html/toppathways-1.png)<!-- -->



## Phosphorylation

Assess changes on phosphoproteomic level and their impact on TPP proteins


```r
psites <- phosphoproteomics %>%
  mutate(
    gene_name = mapIds(org.Hs.eg.db, representative, "SYMBOL", "UNIPROT"),
    condition_long = paste(condition, comparison, sep = "_"),
    UP = representative,
    GeneSymbol = gene_name,
    aa = str_extract(psite, "[[:upper:]]"),
    pos = as.numeric(str_extract(psite, "\\d{1,4}")),
    ID_cosmos = paste(GeneSymbol, "_", aa, pos, sep = ""),
    ID_cosmos = gsub(",.+$", "", ID_cosmos))%>%
  filter(!is.na(logFC) & !is.na(ID_cosmos)&comparison == "Ola") %>%
  dplyr::mutate(
    condition_long = stringr::str_replace(condition_long, "_Ola", "uM"),
    condition_long = stringr::str_replace(condition_long, "BRCA_hom", "UWB1.289_BRCA1"),
    condition_long = stringr::str_replace(condition_long, "BRCA_nul", "UWB1.289"),
    condition_long = stringr::str_replace(condition_long, "04", "0.4")
  ) %>% 
  filter(condition_long == "UWB1.289_24h_4uM") %>% 
  distinct(ID_cosmos, condition_long, .keep_all = T) %>% 
  dplyr::select(ID_cosmos, gene_name, logFC, t, condition_long, adj.P.Val)

TPP_hits_phosphositelevel <- inner_join(TPPFstatistic_hits_adap, psites, by = "gene_name") %>%
  #filter(hit_call != "footprinting") %>%3
  mutate(significance = ifelse(adj.P.Val < 0.05, "significant", "not significant"))
```


```r
volcano_affected_psites <- ggplot(TPP_hits_phosphositelevel,
                                  aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_hline(yintercept = 0, colour = "gray28")+
  geom_vline(xintercept = 0, colour = "gray28")+
  geom_point(aes(colour = significance))+
  geom_text_repel(label = TPP_hits_phosphositelevel$gene_name, size = 3, colour = "gray2") +
  scale_colour_manual(values = c("grey", "darkcyan")) +
  theme_bw(base_size = 14)
volcano_affected_psites
```

![](O1M2_TPPhitcalling_files/figure-html/TPP phosphorylation-1.png)<!-- -->

# TPP activities

To estimate the activity of TPP hits we searched for consistent edges between upstream kinases identified by viper and downstream TPP proteins.


```r
# load viper results
viper_phospho_top <- viperRes_phospho_corrected %>% 
  group_by(condition_long) %>% 
  top_n(30, wt = absNES)

viper_trans_top <- viperRes_trans %>% 
  group_by(condition_long) %>% 
  top_n(30, wt = absNES)
```



```r
# filter prior knowledge links
annoated_ptms_measured <- KSN %>%
  separate(substrate_genesymbol, into = c("substrate_genesymbol", "site"), sep = "_") %>% 
  # subset prior knowledge to subset with measured kinase and target
  filter(enzyme_genesymbol %in% viper_phospho_top$regulator &
    substrate_genesymbol %in% TPPFstatistic_hits_adap$gene_name) %>%
  left_join(viper_phospho_top, by = c("enzyme_genesymbol" = "regulator")) %>%
  # predict for each psite an activation state based on the upstream kinase activity
  mutate(Kinase_activity = sign(NES)) %>%
  dplyr::select(enzyme_genesymbol, sign, substrate_genesymbol, site, Kinase_activity) %>%
  mutate(substrate_activity_predicted =sign*Kinase_activity)
```


## Combine TPP activity information for COSMOS


```r
TPP_cosmos_input <- annoated_ptms_measured %>% 
  # only keep TPP proteins for which there is consistent activation or inactivation information 
  # necessary because we do not no which upstream kinase is potentially regulating the TPP protein
  # the others can be modelled active or inactive supported by (de)stabilization information
  distinct(TPP_protein = substrate_genesymbol, substrate_activity_predicted) %>% 
  group_by(TPP_protein) %>%
  summarise(substrate_activity_predicted_summ = sign(mean(substrate_activity_predicted)))

# aggregate final COSMOS input df
TPP_cosmos_input_final <- TPP_cosmos_input %>%
  mutate(substrate_activity_predicted_summ = ifelse(substrate_activity_predicted_summ == 0, NaN, substrate_activity_predicted_summ)) %>%
  bind_rows(data.frame(TPP_protein = TPPFstatistic_hits_adap$gene_name, substrate_activity_predicted_summ = NaN)) %>%
  distinct(TPP_protein, .keep_all = T)
```

# Remove inconsistent edges

To support the optimization to model consistent links between kinases and TPP proteins as assumed for the activity estimation, we removed incosistent prior knowledge edges. 

For example: Kinase (active) --activates_via_S10--> TPP protein (S10 dephosphorylated) would be removed, because the dephosphorylation of this psite can not be explained using this link (but this information is not always reflected in TPP activity).

To do this correction we use the abundance corrected version of the phosphoproteomics data because we are interested in exclusively phosphorylation driven changes.


```r
# Format phosphoproteomics data after intensity correction
psites_corrected <- limma_UWB1.289_initialdataset_correctedphospho %>%
  filter(condition == "BRCA_nul_24h_4" & comparison == "Ola") %>%
  mutate(
    gene_name = mapIds(org.Hs.eg.db, representative, "SYMBOL", "UNIPROT"),
    condition_long = paste(condition, comparison, sep = "_"),
    UP = representative,
    GeneSymbol = gene_name,
    aa = str_extract(psite, "[[:upper:]]"),
    pos = as.numeric(str_extract(psite, "\\d{1,4}"))
  ) %>%
  mutate(ID_cosmos = paste(GeneSymbol, "_", aa, pos, sep = "")) %>%
  dplyr::select(ID_cosmos, logFC, t, condition_long) %>%
  filter(!is.na(logFC) & !is.na(ID_cosmos)) %>%
  mutate(ID_cosmos = gsub(",.+$", "", ID_cosmos))


# format KSN

KSN <- KSN %>% 
  separate(col = substrate_genesymbol, c("substrate_protein", "psite"), remove = FALSE)
```

Format footprinting and TPP results and prepare for end object which is exported for COSMOS.


```r
Kinases_input_COSMOS <- sign(viper_phospho_top$NES)
names(Kinases_input_COSMOS) <- paste0(
  "X",
  mapIds(org.Hs.eg.db, viper_phospho_top$regulator, "ENTREZID", "SYMBOL")
)
Kinases_input_COSMOS <- Kinases_input_COSMOS[names(Kinases_input_COSMOS) != "XNA"]

TPPFstatistic_hits_adap <- TPPFstatistic_hits_adap %>% 
  mutate(signed_score = ifelse(mode == "destabilized" | mode == "decreased abundance", -1, 1))
```

## Combine prior knowledge, phosphorylation and TPP

Perform the kinase sign filtering to remove links btw kianses and TPP proteins


```r
# extract kinase targets from TPP proteins
kinase_targets <- TPPFstatistic_hits_adap %>% 
  filter(gene_name %in% KSN$substrate_protein) %>% 
  mutate(cosmosID = paste0("X", mapIds(org.Hs.eg.db,gene_name, "ENTREZID", "SYMBOL")))

# merge with limma results
KSN_withFCs <- KSN %>%
  left_join(psites, by = c("substrate_genesymbol" = "ID_cosmos")) %>%
  dplyr::select(-condition_long)
KSN_TPPtargets <- KSN_withFCs %>%
  filter(substrate_protein %in% kinase_targets$gene_name) %>% 
  mutate(source_cosmos =  paste0("X", mapIds(org.Hs.eg.db, enzyme_genesymbol, "ENTREZID", "SYMBOL")))

# use edges with viper kinase
inputs <- Kinases_input_COSMOS[names(Kinases_input_COSMOS) %in% KSN_TPPtargets$source_cosmos]

# merge
KSN_TPPtargets$kinase_sign <- sapply(KSN_TPPtargets$source_cosmos, function(x, inputs) {
  if (x %in% names(inputs)) {
    return(as.numeric(inputs[x]))
  } else {
    return(NA)
  }
}, inputs = inputs)

# remove edges without kinase measurement
KSN_TPPtargets_filt <- filter(KSN_TPPtargets, !is.na(logFC) & !is.na(kinase_sign))
# get inconsistent edges
KSN_TPPtargets_filt <- KSN_TPPtargets_filt %>% 
  filter(sign(logFC * sign) != sign(kinase_sign)) %>% 
  select(substrate_genesymbol, enzyme_genesymbol, sign, logFC, kinase_sign)

KSN_TPPtargets_filt
```

```
##   substrate_genesymbol enzyme_genesymbol sign       logFC kinase_sign
## 1           CDCA5_T159              CDK1    1 -0.06879140           1
## 2            CDCA5_S75              CDK1    1 -0.11694422           1
## 3             RRN3_S44              CDK2    1 -0.02706761           1
```


# Save


```r
save(TPP_cosmos_input_final, viper_trans_top, viper_phospho_top, KSN_TPPtargets_filt,
    file = "220512_cosmosinput_merged.RData")
```
