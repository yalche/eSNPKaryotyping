
#' CreateVcfTable
#'
#' This function returns the square of a number.
#'
#' @param x A numeric value.
#' @return The square of the input number.
#' @export
CreateVcfTable<-function(accession, vcf_path, chr_data) {
  table=EditVCF(Directory = paste0(vcf_path, accession, '/'), Organism = "Human")
  table$chr=as.numeric(table$chr)
  table=table[order(table$chr,table$position),]
  table=table[table$chr>0,]

  return(table)
}


