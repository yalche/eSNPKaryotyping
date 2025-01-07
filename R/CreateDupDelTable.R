#' Create Dup Del Table
#'
#' @param window
#' @param accession
#' @param vcf_path
#' @param output_path
#' @param isPlot True or False, default is True
#' @return
#' @export
CreateDupDelTable<-function(windows_list, accessions_list, vcf_path, output_path, isPlot=TRUE, isTable=TRUE) {
  chr_data <- CalculateCentromereSize()
  final_table <- data.frame(accession = c(), chr = c(), window = c(), p_value = c(), significance = c(), deletion = c())

  for (accession in accessions_list) {
    file_path = paste0(vcf_path, accession, "/VcfTable.csv")
    if (file.exists(file_path)) {
      vcf_table <- read.csv(file_path)
    }
    else {
      vcf_table <- CreateVcfTable(accession, vcf_path, chr_data)
    }

    deletion_vcf_table=MajorMinorCalc(Table = vcf_table,minDP = 1,maxDP = 100000000,minAF = 0)
    duplication_vcf_table = MajorMinorCalc(Table = vcf_table,minDP = 20,maxDP = 100000000,minAF = 0.2)
    duplication_vcf_table$position2=duplication_vcf_table$position+chr_data$chr_total[duplication_vcf_table$chr]


    for (window in windows_list) {
      if(isTable == TRUE) {
        duplication_df <- CreateDuplicationTable(window, duplication_vcf_table, chr_data)
        deletion_df <- CreateDeletionTable(deletion_vcf_table, chr_data)
        if (nrow(duplication_df) != 0) {
          merged <- merge(deletion_df, duplication_df, by = "chr")
          merged$accession = accession
          merged$window = window

          merged <- merged %>% select(accession, chr, window, p_value, significance, deletion)
          final_table <- rbind(final_table, merged)
        }
      }

      if(isPlot == TRUE) {
        jpeg(paste0(output_path, accession, "_", window, ".jpg"), width = 800, height = 300, quality = 100)
        PlotGenome(orderedTable = duplication_vcf_table,Window = window,Ylim = 3,PValue = "TRUE")
        dev.off()
      }
    }
  }

  write.csv(final_table, paste0(output_path, "DupDel.csv"))

}

