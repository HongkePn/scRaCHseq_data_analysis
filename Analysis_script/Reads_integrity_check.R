#setwd("/stornext/Genomics/data/CLL_venetoclax/workspace/PB_LN_comparison/")
#source("/stornext/Genomics/data/CLL_venetoclax/workspace/PB_LN_comparison/renv/activate.R")

library(dplyr)
library(ggplot2)
library(Biostrings)

# QC for FLT-seq data

fq <- list.files("/vast/scratch/users/peng.h/rachseq/demultiplex/")
fq <- fq[grep("fq.gz", fq)]

#seq.ls <- list()
QC_plot.ls <- list()

#collapse reads with same bc and umi
getID <- function(df) {
  df <- df[order(df$length, decreasing = T), ]
  df$num <- nrow(df)
  df <- df[1, ]
  return(df)
}

for (i in fq) {
  print(paste(i, '=================================='))
  x <- Biostrings::readDNAStringSet(paste0("/vast/scratch/users/peng.h/rachseq/demultiplex/", i), format = "fastq")
  
  info <- data.frame(length = width(x), 
                     ID = names(x))
  #subset the beginning of each read
  #rm reads that are less than 15 and longer than 20000
  print(paste(i, "reads shorter than 15bp:", sum(info$length < 15)))
  print(paste(i, "reads longer than 20000bp:", sum(info$length > 20000)))
  keep_df <- info$ID[info$length > 15 & info$length < 20000] 
  
  
  print(paste(i, 'find tso ========================='))
  #find TSO in reads
  x.sub <- x[keep_df]
  rm(x)
  x.start <- subseq(x.sub, start = 1, end = 15)
  ##set mis-match as 3
  y_length_df <- data.frame(length = width(x.sub),
                            category = "demultiplexed")
  
  x.start <- x.start[vcountPattern("CTTATATGGG", x.start, max.mismatch = 3) > 0]
  x.sub <- x.sub[names(x.start)] ## subset read with TSO
  rm(x.start)
  
  ##set mis-match as 3
  tmp <- data.frame(length = width(x.sub),
                    category = "with_TSO")
  y_length_df <- rbind(tmp, y_length_df)
  
  print(paste(i, 'find polyA ========================='))
  #find polyA tails
  x.polyA <- x.sub[vcountPattern("AAAAAAAAAAGTACTCTGCGTTGATACCACTGCTT", x.sub, max.mismatch = 3) > 0]
  
  rm(x.sub)
  ## save reads with polyA
  print("save x.polyA ===========================")
  writeXStringSet(x.polyA, paste0("/vast/scratch/users/peng.h/rachseq/fastq_polyA/", i), format = "fastq", compress = T)
  
  tmp <- data.frame(length = width(x.polyA),
                    category = "with_polyA")
  y_length_df <- rbind(tmp, y_length_df)
  y_length_df$category <- factor(y_length_df$category, levels = c("demultiplexed", "with_TSO", "with_polyA"))
  QC_plot.ls[[i]] <- ggplot(y_length_df, aes(x = length, fill = category)) + geom_histogram(bins = 70, alpha = 0.7, position="identity") + 
    scale_x_continuous(limits = c(0, 5000), breaks = seq(0, 5000, 500)) +
    ggtitle(paste("Read Length Distribution of", strsplit(i, "\\.")[[1]][1], "FLT-seq data")) +
    xlab("length (bp)") + 
    ylab("nCount") + 
    theme_bw() + 
    scale_fill_manual(values = c("grey80", "grey70", "grey40")) +
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 15))
  
  saveRDS(QC_plot.ls, "/vast/scratch/users/peng.h/rachseq/fastq_polyA/QC_plot.rds")
  
}