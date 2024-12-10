#' Calculate Centromere Size
#'
#' This function Calculate Centromere Size
#'
#' @return
#' @export
CalculateCentromereSize<-function() {
  centromere_pos=c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6,12.5)
  centromere_pos=centromere_pos*1000000
  chr_size_list=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
  chr_total=0
  centromere_pos_total = c()

  for (i in 1:(length(chr_size_list)-1)){
    chr_total=c(chr_total,sum(chr_size_list[1:i]))
    centromere_pos_total <- c(centromere_pos_total, (chr_total[i] + centromere_pos[i]))
  }
  return(list(chr_total = chr_total, centromere_pos_total = centromere_pos_total, sum_chr = sum(chr_size_list)))

}

