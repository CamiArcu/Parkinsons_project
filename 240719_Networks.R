# Netwroks

library(tidyverse)
# library("babelgene")

setwd("/data_husihuilke/shared/SMB/carcu/Parkinsons_project")

Mouse_sc_whole_brain <- read_tsv("Data/Networks/GRNdb_mouse_sc_whole_Brain-regulons.txt")
Human_gtex_basal <- read_tsv("Data/Networks/GRAND_Brain_Basal_Ganglia.csv")
Human_gtex_cerebellum <- read_tsv("Data/Networks/GRAND_Brain_Cerebellum.csv")
Human_gtex_other <- read_tsv("Data/Networks/GRAND_Brain_Other.csv")
Human_gtex_brain <- read_tsv("Data/Networks/GRNdb_Human_Brain_GTEx-regulons.txt")

####

geneID.SNup.GO.0050808.upvsdown <- readRDS("Outputs/geneID.SNup.GO.0050808.upvsdown.rds")
geneID.SNdown.GO.0050808.upvsdown <- readRDS("Outputs/geneID.SNdown.GO.0050808.upvsdown.rds")
genes.SNup <- readRDS("Outputs/genes.SNup.rds")
genes.SNdown <- readRDS("Outputs/genes.SNdown.rds")
genes.SNup.down <- readRDS("Outputs/genes.SNup.down.rds")


Mouse_sc_whole_brain <- mutate(Mouse_sc_whole_brain,
                               is.TF.SNup.GO.0050808 = TF %in% geneID.SNup.GO.0050808.upvsdown,
                               is.TF.SNdown.GO.0050808 = TF %in% geneID.SNdown.GO.0050808.upvsdown,
                               is.gene.SNup.GO.0050808 = gene %in% geneID.SNup.GO.0050808.upvsdown,
                               is.gene.SNdown.GO.0050808 = gene %in% geneID.SNdown.GO.0050808.upvsdown,
                               gene_class = if_else(is.gene.SNup.GO.0050808, "SNup", "none"),
                               gene_class = if_else(is.gene.SNdown.GO.0050808, "SNdown", gene_class)
                               )

# Only regulons in which our disered genes participates
Mouse_sc_whole_brain_restricted_TF <- filter(Mouse_sc_whole_brain, 
                                          (gene %in% c(geneID.SNup.GO.0050808.upvsdown, geneID.SNdown.GO.0050808.upvsdown) |
                                             TF %in% c(geneID.SNup.GO.0050808.upvsdown, geneID.SNdown.GO.0050808.upvsdown)))$TF
Mouse_sc_whole_brain_restricted <- filter(Mouse_sc_whole_brain, TF %in% Mouse_sc_whole_brain_restricted_TF)
# This still returns a very broad network

# Only TF (and its targets) or targets that are our desired genes
Mouse_sc_whole_brain_hyperestricted <- filter(Mouse_sc_whole_brain, 
                                             (gene %in% c(geneID.SNup.GO.0050808.upvsdown, geneID.SNdown.GO.0050808.upvsdown) |
                                                TF %in% c(geneID.SNup.GO.0050808.upvsdown, geneID.SNdown.GO.0050808.upvsdown)))
write_tsv(Mouse_sc_whole_brain_hyperestricted, "Outputs/GRNdb_mouse_sc_whole_Brain-regulons_hyperrestricted.txt")

# Only TF (and its targets) or targets that are our desired genes (which are all up or down regulñated genes associated with a common term in the up and down dataset)
Mouse_sc_whole_brain <- mutate(Mouse_sc_whole_brain,
                               is.TF.SNup = TF %in% genes.SNup,
                               is.TF.SNdown = TF %in% genes.SNdown,
                               is.gene.SNup = gene %in% genes.SNup,
                               is.gene.SNdown = gene %in% genes.SNdown,
                               gene_class_allupdown = if_else(is.gene.SNup, "SNup", "none"),
                               gene_class_allupdown = if_else(is.gene.SNdown, "SNdown", gene_class_allupdown)
                               )

Mouse_sc_whole_brain_hyperestricted_allupdown <- filter(Mouse_sc_whole_brain, 
                                              (gene %in% genes.SNup.down |
                                                 TF %in% genes.SNup.down) )

write_tsv(Mouse_sc_whole_brain_hyperestricted_allupdown, "Outputs/GRNdb_mouse_sc_whole_Brain-regulons_hyperrestricted_allupdown.txt")


# saveRDS(LC_SN_2mo_2mo_upboth_genes, "Outputs/LC_SN_2mo_2mo_upboth_genes.rds")
# saveRDS(LC_SN_2mo_2mo_downboth_genes, "Outputs/LC_SN_2mo_2mo_downboth_genes.rds")
# 
# saveRDS(LC_SN_1mo_1mo_upboth_genes, "Outputs/LC_SN_1mo_1mo_upboth_genes.rds")
# saveRDS(LC_SN_1mo_1mo_downboth_genes, "Outputs/LC_SN_1mo_1mo_downboth_genes.rds")

# Physical STRING interactors from mouse
aSyn_STRING <- read_tsv("Data/Networks/STRING/Snca_string_interactions_1st_layer_physical_texmining_database_experiment.tsv")
aSyn_STRING_asyn <- filter(aSyn_STRING, `#node1` == "Snca") #65 interactors
write_tsv(aSyn_STRING_asyn, "Outputs/aSyn_STRING_asyn.tsv")

# Physical interaction in HEK proteomics experiment.
nucleolar_aSyn_interactors <- readxl::read_xlsx("Data/Networks/PRMT5_Baf_aSyn/febs17037-sup-0001-supinfo_aS_Interacting_Proteins.xlsx", sheet = 1)
nucleolar_aSyn_interactors$`Uniprot Name  Gene symbol` <- str_split_fixed(nucleolar_aSyn_interactors$`Uniprot Name  Gene symbol`, " ", 2)[,2] %>% gsub(pattern = " ", replacement = "")
nucleolar_aSyn_interactors <- c(nucleolar_aSyn_interactors$`Gene symbol...4`, nucleolar_aSyn_interactors$`Gene symbol...8`, nucleolar_aSyn_interactors$`Gene symbol...12`, nucleolar_aSyn_interactors$`Uniprot Name  Gene symbol`) %>% na.omit()
non_nucleolar_aSyn_interactors <- readxl::read_xlsx("Data/Networks/PRMT5_Baf_aSyn/febs17037-sup-0001-supinfo_aS_Interacting_Proteins.xlsx", sheet = 2)
non_nucleolar_aSyn_interactors <- c(non_nucleolar_aSyn_interactors$`Gene symbol...4`, non_nucleolar_aSyn_interactors$`Gene symbol...8`, non_nucleolar_aSyn_interactors$`Gene symbol...12`) %>% na.omit()
aSyn_interactors <- data.frame(Gene = c(nucleolar_aSyn_interactors, non_nucleolar_aSyn_interactors),
                               Nucleolar = rep(c(T,F), 
                                               times = c(length(nucleolar_aSyn_interactors), length(non_nucleolar_aSyn_interactors))
                                               ))

orthologs <- read_tsv("Data/Orthologs/mart_export_human2mouse.txt") %>% dplyr::select(`Gene name`, `Mouse gene name`) %>%  unique()
orthologs_mouse <- read_tsv("Data/Orthologs/Jax_HMD_HumanPhenotype.txt", col_names = F)
names(orthologs_mouse) <- c("SYMBOL", "ENTREZ", "MGI_SYMBOL", "MGI", "MP", "none")
orthologs_mouse <- dplyr::select(orthologs_mouse, SYMBOL, MGI_SYMBOL) %>% unique()

aSyn_interactors_equiv <- left_join(aSyn_interactors, orthologs, by = c("Gene" = "Gene name"))
aSyn_interactors_equiv <- left_join(aSyn_interactors_equiv, orthologs_mouse, by = c("Gene" = "SYMBOL"))
aSyn_interactors_equiv <- mutate(aSyn_interactors_equiv,
                                 Ortholog = if_else(is.na(aSyn_interactors_equiv$MGI_SYMBOL), `Mouse gene name`, MGI_SYMBOL))

aSyn_interactors_equiv$Ortholog %>%  is.na() %>% which()
aSyn_interactors_equiv <- na.omit(aSyn_interactors_equiv) %>% 
  dplyr::select(Ortholog, Nucleolar) %>% unique()
write_tsv(aSyn_interactors_equiv, "Outputs/aSyn_interactors_PRMT5_paper.tsv")

####
LC_1mo_updown_genes <- readRDS("Outputs/LC_1mo_updown_genes.rds")
LC_1mo_up_genes <- readRDS("Outputs/LC_1mo_up_genes.rds")
LC_1mo_down_genes <- readRDS("Outputs/LC_1mo_down_genes.rds")

Mouse_sc_whole_brain <- read_tsv("Data/Networks/GRNdb_mouse_sc_whole_Brain-regulons.txt")

Mouse_sc_whole_brain <- mutate(Mouse_sc_whole_brain,
                               is.TF.LCup = TF %in% LC_1mo_up_genes,
                               is.TF.LCdown = TF %in% LC_1mo_down_genes,
                               is.gene.LCup = gene %in% LC_1mo_up_genes,
                               is.gene.LCdown = gene %in% LC_1mo_down_genes,
                               gene_class_allupdown = if_else(is.gene.LCup, "LCup", "none"),
                               gene_class_allupdown = if_else(is.gene.LCdown, "LCdown", gene_class_allupdown)
)

Mouse_sc_whole_brain_hyperestricted_allupdown_LC <- filter(Mouse_sc_whole_brain, 
                                                        ((gene %in% LC_1mo_updown_genes) |
                                                           (TF %in% LC_1mo_updown_genes)) )

write_tsv(Mouse_sc_whole_brain_hyperestricted_allupdown_LC, "Outputs/GRNdb_mouse_sc_whole_Brain-regulons_hyperrestricted_allupdown_LC.txt")


