
chr_index <- function(index, chr_total, centromere_pos_total) {
  # Check if inputs are valid
  for (i in 1:length(chr_total)) {
    if (chr_total[i] > index) {
      if (index < centromere_pos_total[i - 1]) {
        return(paste0(i-1, " p arm"))
      }
      return (paste0(i-1, " q arm")) # Return the index (adjusted by 1)
    }
  }
  if (index < centromere_pos_total[24]) {
    return(paste0(24, " p arm"))
  }
  return (paste0(24, " q arm")) 
}

ttest=function(x, vcf_table){
  ttt=t.test(x,vcf_table$MajorMinor,alternative = "greater")$p.value
  return(ttt)
}

#' Create Duplication Table
#'
#' This function returns the duplication table
#'
#' @param window 
#' @param vcf_table 
#' @param chr_data
#' @return 
#' @export
CreateDuplicationTable<-function(window, vcf_table, chr_data) {

  df_all = data.frame(pos = c(), ratio = c(), p_value = c())
  for (i in 1:24) {
    temp <- vcf_table[vcf_table$chr==i,]
    pos <- rollmedian(temp$position2,window)
    ratio <- rollmedian(temp$MajorMinor,window)
    p_value <- rollapply(temp$MajorMinor, width = window, FUN = function(x) ttest(x, vcf_table))
    temp_df = data.frame(pos, ratio, p_value)
    df_all = rbind(df_all, temp_df)
  }
  
  df_all$p_value=p.adjust(df_all$p_value,"fdr")
  df_all$p_value=-1*log(df_all$p_value,10)
  
  df_all <- df_all %>%
    rowwise() %>%
    mutate(chr = chr_index(pos, chr_data$chr_total, chr_data$centromere_pos_total)) %>%
    ungroup()
  
  aggregated_df <- df_all %>%
    group_by(chr) %>% 
    summarize(p_value = mean(p_value), genes = n()) %>%   mutate(significance = ifelse(p_value > 2, "High", "Low"))

  return(aggregated_df)
}

