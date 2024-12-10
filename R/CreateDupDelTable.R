#' Create Dup Del Table
#'
#' @param window
#' @param accession
#' @param vcf_path
#' @param output_path
#' @param isPlot True or False, default is True
#' @return
#' @export
CreateDupDelTable<-function(window, accession, vcf_path, output_path, isPlot=TRUE, isTable=TRUE) {
  chr_data <- CalculateCentromereSize()
  vcf_table <- CreateVcfTable(accession, vcf_path, chr_data)

  deletion_vcf_table=MajorMinorCalc(Table = vcf_table,minDP = 1,maxDP = 100000000,minAF = 0)
  duplication_vcf_table = MajorMinorCalc(Table = vcf_table,minDP = 20,maxDP = 100000000,minAF = 0.1)
  duplication_vcf_table$position2=duplication_vcf_table$position+chr_data$chr_total[duplication_vcf_table$chr]

  if(isTable == TRUE) {
    duplication_df <- CreateDuplicationTable(window, duplication_vcf_table, chr_data)
    deletion_df <- CreateDeletionTable(deletion_vcf_table, chr_data)

    merged <- merge(deletion_df, duplication_df, by = "chr")
    merged$accession = accession

    merged <- merged %>% select(accession, chr, p_value, significance, deletion)
    write.csv(merged, paste0(output_path, accession, "_", "DupDel.csv"))
  }
  if(isPlot == TRUE) {
    jpeg(paste0(output_path, accession, "_", window, ".jpg"), width = 800, height = 300, quality = 100)
    PlotGenome(orderedTable = duplication_vcf_table,Window = window,Ylim = 3,PValue = "TRUE")
    dev.off()
  }

}

