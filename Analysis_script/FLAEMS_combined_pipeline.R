setwd("/stornext/Genomics/data/CLL_venetoclax/workspace/FLAMES/")
source("/stornext/Genomics/data/CLL_venetoclax/workspace/FLAMES/renv/activate.R")

library("FLAMES")

sces <- FLAMES::sc_long_multisample_pipeline(annotation = "/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/gencode.v33.annotation.gtf",
  fastqs = "/stornext/Genomics/data/CLL_venetoclax/RaCHseq_paper/re_basecall/rachseq/fastq_polyA/filtered_fastq",
  outdir = "/stornext/Genomics/data/CLL_venetoclax/RaCHseq_paper/re_basecall/rachseq/fastq_polyA/filtered_fastq/flames_out",
  genome_fa = "/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/GRCh38.primary_assembly.genome.fa",
  minimap2_dir = locate_minimap2_dir(),
  match_barcode = FALSE,
  config_file = "/stornext/Genomics/data/CLL_venetoclax/RaCHseq_paper/re_basecall/rachseq/config_sclr_capture.json")

if (!is.null(sces)) {
  saveRDS(sces, "/stornext/Genomics/data/CLL_venetoclax/RaCHseq_paper/re_basecall/rachseq/fastq_polyA/filtered_fastq/flames_out/sces.rds")
}
