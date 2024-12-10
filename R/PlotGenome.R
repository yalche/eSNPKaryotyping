#' PlotGenome
#'
#' Plot Allelic ratio along the genome for duplication detection.
#' @param orderedTable The variable containing the output of the MajorMinorCalc function
#' @param Window an odd number determining the window of the moving average plot, usually 151
#' @param Ylim the plot y axes maximal limit, usually 3
#' @param PValue Determine if to add the P-Value bar to the graph, accepts "TRUE" or "FALSE" values
#' @param Organism "Human" or "Mouse"
#' @export
#' @return None

PlotGenome<-function(orderedTable,Window,Ylim,PValue){
  #par(mfrow=c(1,1))
  par(mar=c(5, 4, 4, 2))
  window=Window
  #moving median
  

  chr_data=CalculateCentromereSize()
  chr_total <- chr_data$chr_total
  
  genome_size <- chr_data$chr_sum
  
  orderedTable$position=orderedTable$position+chr_total[orderedTable$chr]
  
  col=rollmedian(orderedTable$chr,window)%%2+1
  col[col==1]="dodgerblue4"
  col[col==2]="dodgerblue1"
  
  plot(rollmedian(orderedTable$position,window),rollmedian(orderedTable$MajorMinor,window),col="dimgrey",pch=15,cex=0.4,ylim=c(1,Ylim), ylab="Allelic Ratio",typ="l",xlab="",xaxt = "n",xlim=c(1,genome_size))
  points(rollmedian(orderedTable$position,window),rollmedian(orderedTable$MajorMinor,window),col=col,pch=15,cex=0.4,ylim=c(1,Ylim), ,xlab="Chromosomal Position")
  
  
  if(PValue==TRUE){
    g=orderedTable$MajorMinor
    
    ttest=function(x){
      ttt=t.test(x,g,alternative = "greater")$p.value
      return(ttt)
    }
    
    tt=rollapply(g,width=window,FUN=ttest)
    
    
    
    tt=p.adjust(tt,"fdr")
    tt=-1*log(tt,10)
    cl=colorpanel(50000, "white","grey", "red")
    pmax=ceiling(max(tt))
    if(pmax<2){pmax=2}
    l=length(rollmedian(orderedTable$MajorMinor,window))
    xx=rollmedian(orderedTable$position,window)
    
    for(i in 1:l){
      lines(c(xx[i],xx[i]),c(2.5,2.7),col=cl[round(tt[i]*(50000/pmax))])
      if (tt[i]>2){lines(c(xx[i],xx[i]),c(2.85,2.75),col="gray34")}
    }
    
    left=50000000
    for(i in 1:100){
      rect(left,3.1,left+5000000,3.3,xpd=T,col=cl[i*500],border =NA )
      left=left+5000000
    }
    left=50000000
    lines(c(left,550000000),c(3.1,3.1),xpd=T)
    lines(c(left,550000000),c(3.3,3.3),xpd=T)
    lines(c(left,left),c(3.1,3.3),xpd=T)
    lines(c(550000000,550000000),c(3.1,3.3),xpd=T)
    
    
    mtext(side = 3,at=0,"0",line=0.8,cex=0.8)
    mtext(side = 3,at=600000000,pmax,line=0.8,cex=0.8)
    mtext(side = 3,at=600000000/2,"-Log10(P-Value)",line=2.5,cex=0.8)
    
  }
  

  for(i in 1 :24){
    if (i>1){lines(c(chr_total[i],chr_total[i]),c(1,Ylim),col="gray48")}
    lines(c(chr_total[i]+centromere_pos[i],chr_total[i]+centromere_pos[i]),c(1,Ylim),col="gray55",lty=4)
    text(chr_total[i]+chr_size[i]/2,Ylim-0.1,i,cex=0.8)
  }


  sum=0
  for (i in 1:(length(chr_size)-1)){
    #lines(c(sum,sum+chr_size[i]),c(1,1),col=i,lwd=10)
    
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=10,col="black")  
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=8,col="gray50")
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=7,col="gray53")
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=6,col="gray59")
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=4,col="gray75")
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=2,col="gray85")
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=1,col="gray90")
    lines(c(sum+centromere_pos[i],sum+centromere_pos[i]),c(1.01,0.99),col="grey13",lwd=2)
    
    sum=sum+chr_size[i]
  }
  
  lines(c(sum+20000000,genome_size-20000000),c(1,1),lwd=10,col="black")  
  lines(c(sum+20000000,genome_size-20000000),c(1,1),lwd=8,col="gray50")
  lines(c(sum+20000000,genome_size-20000000),c(1,1),lwd=7,col="gray53")
  lines(c(sum+20000000,genome_size-20000000),c(1,1),lwd=6,col="gray59")
  lines(c(sum+20000000,genome_size-20000000),c(1,1),lwd=4,col="gray75")
  lines(c(sum+20000000,genome_size-20000000),c(1,1),lwd=2,col="gray85")
  lines(c(sum+20000000,genome_size-20000000),c(1,1),lwd=1,col="gray90")
  lines(c(sum+centromere_pos[24],sum+centromere_pos[24]),c(1.01,0.99),col="grey13",lwd=2)
  
  
}