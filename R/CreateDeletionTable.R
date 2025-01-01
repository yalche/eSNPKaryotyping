#' Create Deletion Table
#'
#' This function returns the deletion table
#'
#' @param vcf_table
#' @param chr_data
#' @return
#' @export
CreateDeletionTable<-function(vcf_table, chr_data) {

  tbl = vcf_table
  tbl$snp = tbl$AF2>0

  total_homo=sum(tbl$snp==FALSE)
  total_hetro=sum(tbl$snp==TRUE)
  pval=NULL
  rat=NULL

  mx = 24
  cent = c()
  rat = c()

  for (i in 1:mx){
    if (i<(mx-1)){
      chr=i
      tbl2=tbl[tbl$chr==chr,]
      d=tbl2[tbl2$position<chr_data$centromere_pos[i],]
      homo=sum(d$snp==FALSE)
      hetro=sum(d$snp==TRUE)
      if (hetro == 0) {
        hetro <- 1
      }
      rat=c(rat,homo/hetro)

      d=tbl2[tbl2$position>chr_data$centromere_pos[i],]
      homo=sum(d$snp==FALSE)
      hetro=sum(d$snp==TRUE)
      if (hetro == 0) {
        hetro <- 1
      }
      rat=c(rat,homo/hetro)

      cent = c(cent, paste0(i, " p arm"))
      cent = c(cent, paste0(i, " q arm"))

    }
  }

  homo=sum(tbl2$snp==FALSE)
  pval=NULL

  for (i in 1:((mx-2)*2)){
    np=1
    if(is.na(rat[i])==FALSE){np=t.test(rat,mu=rat[i])$p.value}
    pval = cbind(pval,np)}

  pval=p.adjust(pval,method  ="fdr")
  pval=log(pval,10)*(-1)
  rat[is.na(rat)==TRUE]=total_homo/total_hetro

  is_deleted=pval>3 & rat>(total_homo/total_hetro*5)
  df = data.frame(chr=cent, deletion=is_deleted)

  return(df)
}
