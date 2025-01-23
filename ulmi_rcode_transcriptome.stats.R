# Load fasta files
fasta_file1 <- "/home/uni01/UFFF/chano/ULMI/ULMI.GEA/ULMI.LOCAL/ulmi_preliminary_assembly_74593.fasta"
fasta_file2 <- "/home/uni01/UFFF/chano/ULMI/ULMI.GEA/ULMI.LOCAL/ulmi_final_assembly_51904.fasta"

library(seqinr)

# Function to calculate GC content
calculate_gc_content <- function(seq) {
  gc_count <- sum(seq == "G" | seq == "C")
  total_count <- length(seq)
  return((gc_count / total_count) * 100)
}

# Function to calculate N50 and L50
calculate_n50_l50 <- function(contig_lengths) {
  sorted_lengths <- sort(contig_lengths, decreasing = TRUE)
  cumulative_lengths <- cumsum(sorted_lengths)
  half_genome_length <- sum(contig_lengths) / 2
  
  n50 <- sorted_lengths[which(cumulative_lengths >= half_genome_length)[1]]
  l50 <- which(cumulative_lengths >= half_genome_length)[1]
  
  return(list(N50 = n50, L50 = l50))
}

# Function to analyze a FASTA file
analyze_fasta <- function(fasta_file) {
  if (!file.exists(fasta_file)) {
    stop("The file does not exist: ", fasta_file)
  }
  
  # Read the FASTA file
  sequences <- read.fasta(fasta_file)
  
  # Number of contigs
  num_contigs <- length(sequences)
  
  # Total length of all contigs
  total_length <- sum(sapply(sequences, length))
  
  # ID of the largest contig
  contig_lengths <- sapply(sequences, length)
  largest_contig_id <- names(sequences)[which.max(contig_lengths)]
  
  # Number of contigs with more than 500 bp
  contigs_over_500 <- sum(contig_lengths > 500)
  
  # Number of contigs with more than 1000 bp
  contigs_over_1000 <- sum(contig_lengths > 1000)
  
  # GC content
  gc_content <- sapply(sequences, calculate_gc_content)
  avg_gc_content <- mean(gc_content)
  
  # N50 and L50
  n50_l50 <- calculate_n50_l50(contig_lengths)
  
  return(list(
    num_contigs = num_contigs,
    total_length = total_length,
    largest_contig_id = largest_contig_id,
    contigs_over_500 = contigs_over_500,
    contigs_over_1000 = contigs_over_1000,
    gc_content = avg_gc_content,
    N50 = n50_l50$N50,
    L50 = n50_l50$L50
  ))
}

if (file.exists(fasta_file1) & file.exists(fasta_file2)) {
  result1 <- analyze_fasta(fasta_file1)
  result2 <- analyze_fasta(fasta_file2)
  
  # Create a data frame to summarize the results
  summary_table <- data.frame(
    Parameter = c("Number of Contigs", "Total Length", "Largest Contig ID", "Contigs > 500 bp", "Contigs > 1000 bp", "GC Content (%)", "N50", "L50"),
    `Preliminary Assembly` = c(result1$num_contigs, result1$total_length, result1$largest_contig_id, result1$contigs_over_500, result1$contigs_over_1000, result1$gc_content, result1$N50, result1$L50),
    `Final Assembly` = c(result2$num_contigs, result2$total_length, result2$largest_contig_id, result2$contigs_over_500, result2$contigs_over_1000, result2$gc_content, result2$N50, result2$L50)
  )
  
  # Write the summary table to a text file
  write.table(summary_table, file = "assembly_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
  print("Summary table has been written to 'assembly_summary.txt'")
} else {
  stop("One or both of the specified FASTA files do not exist.")
}



