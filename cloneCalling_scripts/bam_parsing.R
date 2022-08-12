#!/usr/bin/env Rscript

library(data.table)
library(parallel)

args = commandArgs(trailingOnly=TRUE)




arg_mapper <- function(celltag.version) {
  celltag.df <- data.frame(version = c("8N-v1", "8N-v2", "8N-v3","multi-v1"),
                           grep.pattern = c("GGT[ACTG]{8}GAATTC",
                                            "GTGATG[ACTG]{8}GAATTC",
                                            "TGTACG[ACTG]{8}GAATTC",
                                            "[ACTG]{3}GT[ACTG]{3}CT[ACTG]{3}AG[ACTG]{3}TG[ACTG]{3}CA[ACTG]{3}"),
                           ct.len = c(17,17,17,28),
                           sec.pattern = c("CAGCTTCCGAG", "CAGCTTCCGAG", "CAGCTTCCGAG","CCGGTAA"),
                           stringsAsFactors = F)

  rownames(celltag.df) <- celltag.df$version

  if (!(any(celltag.version %in% celltag.df$version))) {
    stop('Supported CellTag versions include: 8N-v1, 8N-v2, 8N-v3 and multi-v1')
  }
  return(c(celltag.df[celltag.version, "grep.pattern"], celltag.df[celltag.version, "ct.len"], celltag.df[celltag.version, "sec.pattern"]))
}

custom.bam.process <- function(arg_list) {
  # Install Rsamtools
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    BiocManager::install("Rsamtools")
  }
  library(Rsamtools)

  bam.file <- arg_list$bam
  pattern <- arg_list$p.calling
  bar.length <- arg_list$bar.length
  technique <- arg_list$technique
  assay <- arg_list$assay
  sample.curr <- arg_list$sample
  cb.file <- arg_list$cb
  sec.grep <- arg_list$sec.grep

  # Get the bam file
  bamFile <- BamFile(bam.file)

  # Get the size of the bam file
  bam.size <- file.size(bam.file)
  total <- bam.size/(1000000 * 82.99)
  print(paste0("Reading ", bam.file, " ..."))

  # Initialize the progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  # Initialize the number of lines to read at once
  yieldSize(bamFile) <- 1000000
  open(bamFile)
  if (tolower(assay) == 'rna') {
    if (tolower(technique) == "10x") {
      parameters <- ScanBamParam(what = scanBamWhat(), tag = c("CB", "GN", "UB", "CR"))
    } else {
      if (tolower(technique) == "dropseq") {
        parameters <- ScanBamParam(what = scanBamWhat(), tag = c("XC", "GN", "XM", "GE"))
      } else {
        stop("We don't support your current single-cell sequencing technology. Please contact us to add.")
      }
    }
  }

  else if (tolower(assay) == 'atac') {
    if (tolower(technique) == '10x') {
      parameters <- ScanBamParam(what = scanBamWhat(), tag = c("CB", "CR"))
    }
    else {
      stop("We don't support your current single-cell technology with scATAC")
    }
  }

  else {
    stop("Only single cell RNA and ATAC assays are supported for now")
  }
  bam.parsed.df <- data.table()
  count <- 0
  while(TRUE) {
    curr.read <- scanBam(bamFile, param = parameters)[[1]]
    #    print(count)
    if (length(curr.read$qname) <= 0) {
      break
    } else {
      # Read in all information
      curr.seqs <- as.character(curr.read$seq)

      # Check if the sequences contain the celltag motif
      reg.rslt <- regexpr(pattern, curr.seqs, ignore.case = TRUE, perl = TRUE)
      contain.idx <- which(reg.rslt > 0)
      if (length(contain.idx) > 0) {
        if (tolower(technique) == "10x") {
          curr.cell.bc <- curr.read$tag$CB
          if (tolower(assay) == 'rna') {
            curr.umi <- curr.read$tag$UB
          }
          else {
            curr.umi <- "NNNNNNNN"
          }
          curr.chr <- curr.read$rname
        } else {
          if ((tolower(technique) == "dropseq") & (tolower(assay) == 'rna')){
            curr.cell.bc <- curr.read$tag$XC
            curr.umi <- curr.read$tag$XM
          }
        }
        curr.cell.tag <- rep(NA, length(curr.read$qname))
        if (!(is.null(curr.cell.bc) | is.null(curr.umi))) {

          # Initialize the current data table
          curr.df <- data.table(Cell.BC = curr.cell.bc, UMI = curr.umi, Cell.Tag = curr.cell.tag, Chr = curr.chr, Read = curr.seqs)
          curr.f.seq <- curr.seqs[contain.idx]
          start.loc <- reg.rslt[contain.idx]
          end.loc <- start.loc + bar.length - 1

          curr.full.tag <- substr(curr.f.seq, start = start.loc, stop = end.loc)
          # only.tag <- substr(curr.full.tag, start = (nchar(short.nt.before.tag) + 1), stop = (nchar(short.nt.before.tag) + 8))

          curr.df$Cell.Tag[contain.idx] <- curr.full.tag

          # Add to the current data frame
          if (nrow(bam.parsed.df) <= 0) {
            bam.parsed.df <- curr.df[contain.idx,]
          } else {
            bam.parsed.df <- rbind(bam.parsed.df, curr.df[contain.idx, ])
          }
        }
      }
    }
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
  close(bamFile)
  close(pb)

  dir.create(file.path("temp",sample.curr), recursive = TRUE)

  write.table(bam.parsed.df,
              file.path("temp",sample.curr,paste0(sample.curr,"_bam_parse_raw.txt")),
              quote = FALSE, sep = "\t", row.names = FALSE)
  print(paste('Finished processing: ', sample.curr))  

  #remove non-CB
  cb <- read.csv(cb.file, header = FALSE)
  bam.parsed.df <- bam.parsed.df[bam.parsed.df$Cell.BC %in% cb$V1,]

  #secondary grepping 
  new_table_1 <- bam.parsed.df[is.na(bam.parsed.df$Chr),]
  new_table_2 <- bam.parsed.df[!is.na(bam.parsed.df$Chr),]
  new_table_2 <- new_table_2[grep(sec.grep, new_table_2$Read),]

  new_table_1 <- rbind(new_table_1, new_table_2)
  new_table_1 <- new_table_1[,c("Cell.BC","UMI","Cell.Tag")]


  #remove NAs
  new_table_1 <- new_table_1[rowSums(is.na(new_table_1))==0,]

  dir.create(file.path("celltag_reads",sample.curr), recursive = TRUE)

  write.table(new_table_1,
              file.path("celltag_reads",sample.curr,paste0(sample.curr,"_bam_parse.txt")),
              quote = FALSE, sep = "\t", row.names = FALSE)
  print(paste('Finished processing: ', sample.curr))
}

parse_args <- function(x) {
  list_curr <- list()
  
  if(length(x) == 5) {
    temp <- arg_mapper(x[[5]])
    list_curr$p.calling <- temp[1]
    list_curr$bar.length <- as.integer(temp[2])
    list_curr$sec.grep <- temp[3]
  } else if(length(x) == 7) {
      list_curr$p.calling <- x[[5]]
      list_curr$bar.length <- as.integer(x[[6]])
      list_curr$sec.grep <- x[[7]]
  } else {
  stop("config file not valid")
  }
  
  list_curr$sample <- x[[1]]
  list_curr$bam <- x[[2]]
  list_curr$cb <- x[[3]]
  list_curr$assay <- x[[4]]
  list_curr$technique <- '10x'
  
  return(list_curr)
  
}




#read in config file and parse arguments
config_file <- read.csv(args[1], header=TRUE)
all_args <- apply(config_file, 1, parse_args)

mclapply(all_args, custom.bam.process, mc.preschedule=FALSE)

