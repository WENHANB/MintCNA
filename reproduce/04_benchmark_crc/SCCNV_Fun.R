SCCNV<-function(dat1,
                minmap=0.3,
                maxN = 0.1,
                ploidy = 2,
                windowsize = 10,
                binsize =1000){
  
  cell_number = ncol(dat1)-7
  genomesize = nrow(dat1)*binsize
  ### Normalize by mappability: dat2 ###
  dat2<-dat1[which(dat1[,2]>=minmap & dat1[,7]<=maxN),]
  dat2[,8:ncol(dat2)]<-dat2[,8:ncol(dat2)]/dat2[,4]
  
  
  # colnames(dat2)[1:10]
  # dat2[1:10,8]
  # dat1[1:10,8]/dat1[1:10,4]
  
  ### Normalize by GC content: dat3 ###
  print("Normalize by GC content")
  gc_percentile <- quantile(dat2[,6],seq(0,100)/100)
  # gc_percentile[1]
  # min(dat2[,6])
  
  # the quantile of each RC
  gc_pp<-rep(-1,nrow(dat2))
  # i<-2
  for (i in 1:nrow(dat2)){
    
    for (j in 1:(length(gc_percentile)-1)){
      if (dat2[i,6]>=gc_percentile[j] & dat2[i,6]<=gc_percentile[j+1]){
        gc_pp[i]<-(j-1)
      }
    }
  }
  # which(gc_pp==(-1))
  
  
  # sort(unique(gc_pp))
  # ii<-1
  # RC<-dat2[,13]
  # sort(unique(gc_pp))
  gcNorm<-function(RC,gc_pp=gc_pp){
    tmp.ave.seq<-c()
    # ii<-32
    for (ii in 0:100){
      tmp.ind<-which(gc_pp==ii)
      if (length(tmp.ind)==0){tmp.ave.seq<-c(tmp.ave.seq,1)}
      if (length(tmp.ind)>0){tmp.ave.seq<-c(tmp.ave.seq,mean(RC[tmp.ind]))}
    }
    # if there are no corresponding quantiles, use the average of adjacent quantile means
    com.seq<-which(tmp.ave.seq==1)
    candi.seq<-which(tmp.ave.seq!=1)
    # k<-com.seq[1]
    for (kk in com.seq){
      tmp.ind2<-which.min(abs(kk-candi.seq))
      tmp.ave.seq[kk]<-mean(tmp.ave.seq[c(tmp.ind2,tmp.ind2+1)])
    }
    
    
    overall.mean<-mean(RC)
    
    tmp.nr_n<-c()
    # jj=15
    for (jj in 1:length(RC)){
      tmp.nr_n<-c(tmp.nr_n,overall.mean*RC[jj]/tmp.ave.seq[gc_pp[jj]+1])
    }
    # which(is.na(tmp.nr_n))
    
    return(tmp.nr_n)
  }
  
  dat3<-dat2
  # i<-8
  for (i in 8:ncol(dat2)){
    # print(i)
    dat3[,i]<-gcNorm(RC=dat2[,i],gc_pp=gc_pp)
  }
  ### Cross batch normalize ###
  print("Cross batch normalize")
  # convert to raw CN estimate: dat4 #
  # tmp.dat4<-dat3[-which(dat3[,1] %in% c("X","Y")),8:ncol(dat3)]
  tmp.dat4<-dat3[,8:ncol(dat3)]
  cell_nr<-apply(tmp.dat4,2,mean)
  
  
  dat4    <-dat3
  for (i in 1:ncol(tmp.dat4)){
    # i<-1
    dat4[,i+7]<-ploidy*tmp.dat4[,i]/cell_nr[i]
  }
  
  # across sample normalization
  
  dat5<-dat4
  temp.dat5<-dat4[,-c(1:7)]
  for (i in 1:nrow(temp.dat5)){
    for (j in 1:ncol(temp.dat5)){
      temp2<-(sum(temp.dat5[i,])-temp.dat5[i,j])/(ncol(temp.dat5)-1)
      if (temp2==0){dat5[i,j+7]=0}else {
        dat5[i,j+7] = (temp.dat5[i,j]/temp2)*ploidy
      }
    }
  }
  
  
  ## For each sample, normalize to the median cnv equals to 2 ##
  nfactor<-apply(dat5[,-c(1:7)],2,median)
  
  dat51<-dat5
  for (i in 8:ncol(dat5)){
    dat51[,i]<-(ploidy/nfactor[i-7])*dat5[,i]
  }
  dat5<-dat51
  # dat51[1:10,8]
  
  print("Futher smoothing")
  ### Futher smoothing: dat6 ###
  step=windowsize/2
  dat6<-dat5
  
  for (i in 1:nrow(dat5)){
    
    for (j in 8:ncol(dat5)){
      p1=i-step
      p2=i+step+1
      if (p1<=0) {p1=1}
      if (p2>nrow(dat5)){p2=nrow(dat5)}
      temp<-c()
      for (k in p1:p2){
        if (dat5[k,1] == dat5[i,1]){
          temp <- c(temp,dat5[k,j])
        }
      }
      dat6[i,j] <- mean(temp)
      
    }
  }
  
  ### CNV calling: dat7 ###
  print("CNV calling")
  # Estimate sigma #
  nr<-dat6[,-c(1:7)]
  sigma<-c()
  for (i in 1:ncol(nr)){
    tmp.sigma<-abs(quantile(nr[,i],0.309)-ploidy)+abs(quantile(nr[,i],0.691)-ploidy)
    sigma <- c(sigma,tmp.sigma)
  }
  
  norm_pdf<-function(x,mm,sg){
    
    prob<-dnorm(x=(x-mm)/sg,mean = 0,sd=1)
    return(prob)
  }
  
  # Define criteria #
  cri=1 - 1 / genomesize * windowsize * binsize
  # cri=0.95
  # cri=0.9
  dat7<-dat6
  # post.prob.list: store the posterior probabilities for each cell in each list element
  post.prob.list<-list()
  for (j in 8:ncol(dat6)){
    
    tmp.prob.data<-data.frame()
    for (i in 1:nrow(dat6)){
      pdf=norm_pdf(x=dat6[i,j],mm=seq(0,20),sg=sigma[j-7])
      ln=pdf/sum(pdf)
      # a<-which(ln>cri)
      # if (length(a)==0){dat7[i,j]=NA}
      # if (length(a)==1){dat7[i,j]=a-1}
      tmp.prob.data<-rbind(tmp.prob.data,ln)
    }
    colnames(tmp.prob.data)<-paste0("CN",seq(0,20))
    post.prob.list[[j-7]]<-tmp.prob.data
    
  }
  
  # CNV res
  return(post.prob.list)
}

