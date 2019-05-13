plot_matrix <- function(A, type="other"){
  library(reshape2)
  library(ggplot2)
  P<-nrow(A)
  K<-ncol(A)
  longData<-cbind(which(!is.na(A),arr.ind = TRUE),as.vector(A))
  longData<-as.data.frame(longData)
  colnames(longData) <- c("row", "col", "value")
  if( type=="Lambda" ){
    yLab = "d"; xLab = "k"
    tit = expression(Lambda~entries)
  } else if( type=="Theta" ){
    yLab = "s"; xLab = "k"
    tit = expression(Theta~entries)
  } else{
    yLab = "p"; xLab = "k"
    tit = "Matrix entries"
  }
  ggplot(longData, aes(x = col, y = row)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0) +
    labs(title=tit) + theme_bw() +
    theme_minimal() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                            axis.text.y=element_text(size=9),
                            plot.title=element_text(size=11)) +
    scale_y_reverse() + ylab(yLab) + xlab(xLab)
}

plot_Lambda_true <- function(Lambda, doses=1:nrow(Lambda)/nrow(Lambda), inds=1:3){
  library(gridExtra)
  library(ggplot2)
  Lambda = as.data.frame(Lambda)
  xLim = range(doses)
  yLim = range(Lambda[,inds])
  if(length(inds)>3){stop("inds must be of length <= 3")}
  if(length(inds)==1){
    p1 = qplot(doses, Lambda[,inds[1]], geom="line") + xlim(xLim) + ylim(yLim) + 
      theme_minimal() + theme_bw() + xlab("dose") + ylab("loading") + 
      ggtitle(paste("Loading",inds[1]))
    grid.arrange(p1, nrow = 1)
  }else if(length(inds)==2){
    p2 = qplot(doses, Lambda[,inds[2]], geom="line") + xlim(xLim) + ylim(yLim) +
      theme_minimal() + theme_bw() + xlab("dose") + ylab("loading") + 
      ggtitle(paste("Loading",inds[2]))
    grid.arrange(p1, p2, nrow = 1)
  }else{
    p3 = qplot(doses, Lambda[,inds[3]], geom="line") + xlim(xLim) + ylim(yLim) +
      theme_minimal() + theme_bw() + xlab("dose") + ylab("loading") + 
      ggtitle(paste("Loading",inds[3]))
    grid.arrange(p1, p2, p3, nrow = 1)
  }
}

plot_Lambda <- function(Lambda_low, Lambda, Lambda_upp, Lambda_true, k, dose=1:nrow(Lambda)/nrow(Lambda)){
  library(reshape2)
  library(ggplot2)
  df = as.data.frame(cbind(dose, Lambda_low[,k], Lambda[,k], Lambda_upp[,k], Lambda_true[,k]))
  colnames(df) = c("dose", "ll", "est", "ul", "truth")
  
  ggplot(data = df, aes(dose, est)) +
    geom_ribbon(aes(x=dose, ymax=ll, ymin=ul), fill="grey", alpha=.5) +
    geom_line(aes(y = ul), lty = 3, colour = 'black') +
    geom_line(aes(y = ll), lty = 3, colour = 'black')+
    geom_line(lty = 2) + geom_line(aes(y = truth), colour = 'red') +
    ylab("loading") + ggtitle(paste("Column",k,"of loadings matrix for Y")) + 
    theme_minimal() + theme_bw()
}

plot_data <- function(Y, true_curve, avg_dose_resp=rep(0,nrow(Y)), 
                      doses=1:nrow(Y)/nrow(Y), inds=1:3){
  library(gridExtra)
  library(ggplot2)
  Y = as.data.frame(Y)
  for(j in 1:ncol(Y)){
    Y[,j] = Y[,j]+avg_dose_resp
    true_curve[,j] = true_curve[,j]+avg_dose_resp
  }
  xLim = range(doses)
  yLim = range(c(Y[,inds],true_curve[,inds]), na.rm=T)
  p1 = qplot(doses, Y[,inds[1]]) + xlim(xLim) + ylim(yLim) + 
    theme_minimal() + theme_bw() + xlab("dose") + ylab("response") + 
    ggtitle("Chemical 1") + geom_line(aes(x=doses,y=true_curve[,inds[1]]))
  p2 = qplot(doses, Y[,inds[2]]) + xlim(xLim) + ylim(yLim) +
    theme_minimal() + theme_bw() + xlab("dose") + ylab("response") + 
    ggtitle("Chemical 2") + geom_line(aes(x=doses,y=true_curve[,inds[2]]))
  p3 = qplot(doses, Y[,inds[3]]) + xlim(xLim) + ylim(yLim) +
    theme_minimal() + theme_bw() + xlab("dose") + ylab("response") + 
    ggtitle("Chemical 3") + geom_line(aes(x=doses,y=true_curve[,inds[3]]))
  
  grid.arrange(p1, p2, p3, nrow = 1)
}
