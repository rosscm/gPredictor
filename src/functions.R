# Extract sequence context of specified CRISPR guide set
#
#' @param input (char) guide library input file (output from prepIn.R)
#' @param genome_build (char) genome build of guide library
#' @param seq_length (int) length of sequence to extract
#' @param strand (char) run context extraction negative or positive DNA strand

getSeqContext <- function(input, genome_build, seq_length, strand) {

  # Need to format CHROMOSOME column differently when using hg38 coordinates
  if (genome_build == "hg38") { input$CHROMOSOME <- gsub("chr", "", input$CHROMOSOME) }

  # Filter data by +/- strand
  input_strand <- filter(input, STRAND == strand)

  if (seq_length == 30) { # 30 bp gRNA + sequence context

    pam_start = 25
    pam_end = 27

    if (strand == "+") {
      if (genome_build == "hg19") {
        start <- input_strand$START-4
        end <- input_strand$STOP+6
      }
      if (genome_build == "hg38") {
        start <- input_strand$START-3
        end <- input_strand$STOP+6
      }
    }
    if (strand == "-") {
      if (genome_build == "hg19") {
        start <- input_strand$START-6
        end <- input_strand$STOP+4
      }
      if (genome_build == "hg38") {
        start <- input_strand$START-5
        end <- input_strand$STOP+4
      }
    }
  }

  if (seq_length == 79) { # 79 bp gRNA + sequence context

    pam_start = 44
    pam_end = 46

    if (strand == "+") {
      if (genome_build == "hg19") {
        start <- input_strand$START-23
        end <- input_strand$STOP+36
      }
      if (genome_build == "hg38") {
        start <- input_strand$START-22
        end <- input_strand$STOP+36
      }
    }
    if (strand == "-") {
      if (genome_build == "hg19") {
        start <- input_strand$START-36
        end <- input_strand$STOP+23
      }
      if (genome_build == "hg38") {
        start <- input_strand$START-35
        end <- input_strand$STOP+23
      }
    }
  }

  # Extracting set of sequences from BSgenome object
  seqs <- BSgenome::getSeq(
    Hsapiens,
    input_strand$CHROMOSOME,
    strand=input_strand$STRAND,
    start=start,
    end=end
  )

  # Convert to dataframe format
  seqs <- as.data.frame(seqs) # gRNA with sequence context
  pam <- substr(seqs$x, pam_start, pam_end) # PAM sequence
  seq_df <- data.frame(
    gene=input_strand$GENE,
    strand=input_strand$STRAND,
    guide=as.character(input_strand$SEQUENCE),
    seqs,
    pam=pam)
  return(seq_df)
}
