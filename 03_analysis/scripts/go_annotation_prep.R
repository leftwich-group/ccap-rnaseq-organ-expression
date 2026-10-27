#!/usr/bin/env Rscript
# go_annotation_prep.R  (optional, run once)
# Parses eggNOG-mapper output into 02_preprocessing/04_go_annotation/geneID2GO_*.rds,
# which R/load_data.R then reads. Only needed if those .rds files are missing or
# the annotation inputs changed. Not part of run_pipeline.sh.
# Run from the repository root:  Rscript 03_analysis/scripts/go_annotation_prep.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/setup.R")) stop(
  "Run from the repository root.", call. = FALSE)
suppressPackageStartupMessages(library(here))
source(here::here("03_analysis/R/setup.R"))
# ---------------------------------------------------------------------------

# ---- eggNOG-mapper annotation (primary) -----------------------------------
# Expected input: eggNOG output TSV for both short and long genome proteins
# Parsed to gene2go list for topGO

outdir_go <- here("01_data/GO_terms/eggNOG-mapper/parsed_gene2go")
dir.create(outdir_go, recursive = TRUE, showWarnings = FALSE)

parse_eggnog <- function(eggnog_tsv, gff_file, genome = "short") {

  # Read GFF to build protein to gene map
  gff <- read_tsv(gff_file, comment = "#",
                  col_names = c("seqid","source","type","start","end",
                                "score","strand","phase","attributes"),
                  show_col_types = FALSE) %>%
    filter(type == "CDS") %>%
    mutate(
      protein = str_extract(attributes, "(?<=protein_id=)[^;]+"),
      gene    = str_extract(attributes, "(?<=gene=)[^;]+"),
      locus   = str_extract(attributes, "(?<=locus_tag=)[^;]+"),
      gene_id = coalesce(gene, locus)
    ) %>%
    filter(!is.na(protein), !is.na(gene_id)) %>%
    # strip CCAP1982_ prefix for long genome
    mutate(gene_id = if (genome == "long")
             sub("^CCAP1982_", "", gene_id) else gene_id) %>%
    distinct(protein, gene_id)

  prot2gene <- setNames(gff$gene_id, gff$protein)

  # Read eggNOG TSV
  egg <- read_tsv(eggnog_tsv, comment = "##", show_col_types = FALSE) %>%
    dplyr::select(query = 1, GOs = contains("GOs")) %>%
    filter(!is.na(GOs), GOs != "-")

  # Expand GO terms
  egg_long <- egg %>%
    mutate(GO_list = str_split(GOs, ",")) %>%
    unnest(GO_list) %>%
    mutate(
      GO_list  = trimws(GO_list),
      gene_id  = prot2gene[query]
    ) %>%
    filter(!is.na(gene_id), str_starts(GO_list, "GO:")) %>%
    distinct(gene_id, go_id = GO_list)

  write_tsv(egg_long,
            file.path(outdir_go,
                      paste0("medfly_", genome, "_gene2go_compact.tsv")))

  split(egg_long$go_id, egg_long$gene_id)
}

geneID2GO_short <- parse_eggnog(
 here("01_data/GO_terms/eggNOG-mapper/medfly_short.emapper.annotations"),
  here("01_data/Genomic_files/GCF_000347755.3_Ccap_2.1_genomic.gff.gz"),
  genome = "short"
)
geneID2GO_long  <- parse_eggnog(
  here("01_data/GO_terms/eggNOG-mapper/medfly_long.emapper.annotations"),
  here("01_data/Genomic_files/GCA_905071925.1_EGII-3.2.1_genomic.gff.gz"),
  genome = "long"
)

saveRDS(geneID2GO_short, here("02_preprocessing/04_go_annotation/geneID2GO_short.rds"))
saveRDS(geneID2GO_long,  here("02_preprocessing/04_go_annotation/geneID2GO_long.rds"))
message("GO annotation setup complete.")
