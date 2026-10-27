#!/usr/bin/env Rscript
# Builds exon and intron SAF annotations for featureCounts from a GTF or GFF3.
# Exons: union of exons across all transcripts of a gene.
# Introns: gene span minus the exons of every gene on either strand (the
# libraries are unstranded), keeping segments of at least min_intron_bp.
# Command line:  Rscript 01_build_exon_intron_features.R <annotation> <out_dir> [min_intron_bp]
# Interactive (RStudio, Positron): edit the three settings below and source the file.
# GFF3 input takes gene IDs from the gene attribute when gene_id is absent.
# Output: <out_dir>/exons.saf, introns.saf, intron_lengths.tsv
# Set FORCE=1 to rebuild existing outputs (Sys.setenv(FORCE = 1) when interactive).
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  gtf_path      <- args[1]
  out_dir       <- args[2]
  min_intron_bp <- if (length(args) >= 3) as.integer(args[3]) else 100L
} else if (interactive()) {
  gtf_path      <- here::here("01_data/Genomic_files/GCF_000347755.3_Ccap_2.1_genomic.gff.gz")
  out_dir       <- here::here("02_preprocessing/06_intron_exon_counts/results")
  min_intron_bp <- 100L
} else {
  stop("Usage: Rscript 01_build_exon_intron_features.R <annotation> <out_dir> [min_intron_bp]",
       call. = FALSE)
}
rebuild <- nzchar(Sys.getenv("FORCE"))

first_value <- function(gr, nm) {
  if (!nm %in% names(mcols(gr))) return(rep(NA_character_, length(gr)))
  vapply(as.list(mcols(gr)[[nm]]),
         function(v) if (length(v)) as.character(v[1]) else NA_character_, "")
}

build_features <- function(gtf_path, out_dir, min_intron_bp, rebuild) {
  if (!file.exists(gtf_path)) stop("Annotation not found: ", gtf_path, call. = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_files <- file.path(out_dir, c("exons.saf", "introns.saf", "intron_lengths.tsv"))
  if (all(file.exists(out_files) & file.size(out_files) > 0) && !rebuild) {
    message("[OK] outputs already exist in ", out_dir, " - skipping (FORCE=1 to rebuild)")
    return(invisible(out_files))
  }

  gtf <- import(gtf_path)
  ex  <- gtf[gtf$type == "exon"]
  if (length(ex) == 0) stop("No exon features in ", gtf_path, call. = FALSE)
  ex$gene_id <- fcoalesce(first_value(ex, "gene_id"), first_value(ex, "gene"))
  if (anyNA(ex$gene_id)) {
    message("Dropping ", sum(is.na(ex$gene_id)), " exons without a gene ID")
    ex <- ex[!is.na(ex$gene_id)]
  }
  message("Exons: ", length(ex), "  Genes: ", length(unique(ex$gene_id)))

  mcols(ex) <- mcols(ex)["gene_id"]
  strand(ex) <- "*"
  by_gene   <- split(ex, ex$gene_id)
  all_exons <- reduce(ex)

  flatten <- function(grl) {
    gr <- unlist(grl, use.names = FALSE)
    names(gr) <- rep(names(grl), elementNROWS(grl))
    gr
  }
  exon_union <- flatten(reduce(by_gene))
  gene_span  <- flatten(range(by_gene))

  # Exons of any gene overlapping each gene span, removed from that span.
  ov     <- findOverlaps(gene_span, all_exons)
  masked <- split(all_exons[subjectHits(ov)],
                  factor(queryHits(ov), levels = seq_along(gene_span)))
  intron_list <- psetdiff(gene_span, masked)
  names(intron_list) <- names(gene_span)
  introns <- flatten(intron_list)
  introns <- introns[width(introns) >= min_intron_bp]
  if (length(introns) == 0) stop("No intron segments of at least ", min_intron_bp, " bp", call. = FALSE)

  to_saf <- function(gr) data.table(GeneID = names(gr), Chr = as.character(seqnames(gr)),
                                    Start = start(gr), End = end(gr), Strand = ".")
  fwrite(to_saf(exon_union), out_files[1], sep = "\t")
  fwrite(to_saf(introns),    out_files[2], sep = "\t")

  lens <- merge(
    data.table(gene_id = names(exon_union), exon_bp = width(exon_union))[, .(exon_bp = sum(exon_bp)), by = gene_id],
    data.table(gene_id = names(introns), intron_bp = width(introns))[, .(intron_bp = sum(intron_bp), n_intron_segments = .N), by = gene_id],
    by = "gene_id", all.x = TRUE)
  lens[is.na(intron_bp), `:=`(intron_bp = 0L, n_intron_segments = 0L)]
  fwrite(lens, out_files[3], sep = "\t")

  message("[OK] exon rows: ", length(exon_union), "  intron rows: ", length(introns),
          "  genes with introns: ", sum(lens$intron_bp > 0))
  invisible(out_files)
}

build_features(gtf_path, out_dir, min_intron_bp, rebuild)
