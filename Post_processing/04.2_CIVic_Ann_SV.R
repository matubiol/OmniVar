#!/usr/bin/env Rscript

# Load/install packages
pkgs <- c("tidyverse", "optparse", "writexl")

for(p in pkgs){
  if(!require(p, character.only = TRUE, quietly = TRUE)){
    install.packages(p, verbose = FALSE)
    library(p, character.only = TRUE, verbose = FALSE)
  }
}

# ============================================================================ #
#                         Setup command line arguments                         #
# ============================================================================ #
if (exists("snakemake")) {
  opt <- list(
    tsv = snakemake@input[["tsv"]],
    civic = snakemake@input[["civic"]],
    xlsx = snakemake@output[["xlsx"]]
  )
} else {
  option_list <- list(
    make_option(c("-i", "--tsv_file"), type="character"),
    make_option(c("-c", "--civic"), type="character"),
    make_option(c("-o", "--output"), type="character")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
}

vars <- read.delim(opt$tsv)
civic <- read.delim(opt$civic) %>%
  separate(
    col = molecular_profile,
    sep = " ",
    into = c("Var", "Var_type"),
  ) %>%
  mutate(
    Var_type_short = case_when(
      grepl("delins", Var_type, ignore.case = T) ~ "COMPLEX",
      Var_type == "Deleterious" ~ "Harmful",
      grepl("del", Var_type, ignore.case = T) ~ "DEL",
      grepl("ins", Var_type, ignore.case = T) ~ "INS",
      grepl("Amplification", Var_type, ignore.case = T) | grepl("dup", Var_type, ignore.case = T) ~ "DUP",
      TRUE ~ Var_type
    ),
    Gene = case_when(
      Var_type == "Fusion" ~ Var,
      Var_type_short == "DUP" ~ paste(Var, "DUP"),
      Var_type_short == "DEL" ~ paste(Var, "DEL"),
      Var_type_short == "INS" ~ paste(Var, "INS"),
      TRUE ~ Var
    )
  ) %>%
  select(Gene, Var_type, disease)

# Filter by list of genes of interest and annotate with diseases
genes_fus = c(
  "ABL1", "ABL2", "BCR", "CBFA2T3", "CBFB", "CREBBP", "CSF1R",
  "EPOR", "ETV6", "FGFR1", "FUS", "IL3", "JAK2", "KMT2A", "MEF2D",
  "MNX1", "NPM1", "NUP214", "NUP98", "PDGFRA", "PDGFRB", "RARA",
  "RBM15", "RUNX1", "TAL1", "TCF3", "TP53", "FIP1L1", "MECOM"
)

vars_ann <- vars %>%
  separate(col = "Gene", into = c("Gene1", "Gene2"), sep = "::", remove = F) %>%
  filter(any(Gene1 %in% genes_fus, Gene2 %in% genes_fus)) %>%
  mutate(
    Gene = if_else(
      SVType %in% c("DUP", "DEL", "INS"),
      paste(Gene, SVType),
      Gene
      )
    ) %>%
  left_join(civic, by = join_by(Gene),
    relationship = "many-to-many") %>%
  mutate(
    Gene = if_else(grepl(" ", Gene), sub(" .*", "", Gene), Gene),
    disease = if_else(disease == "", NA_character_, disease)
  ) %>%
  group_by(Gene, SVType, Support, Coordinates) %>%
  summarise(Disease = paste(disease, collapse = ", "),
            .groups = "drop") %>%
  mutate(Disease = sub(", NA", "", Disease),
         Disease = if_else(Disease == "NA", NA_character_, Disease)) %>%
  select(Gene, SVType, Support, Coordinates, Disease) %>%
  arrange(desc(Support), Disease)

write_xlsx(vars_ann, opt$xlsx)