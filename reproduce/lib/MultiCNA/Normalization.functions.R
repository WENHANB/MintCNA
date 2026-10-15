CorrectGC<-function(RC,GCContent,step){
  stepseq<-seq(min(GCContent),max(GCContent),by=step)
  #stepseq<-seq(0,100,by=step)
  MasterMedian<-median(RC,na.rm=T)
  MedianGC<-rep(0,length(stepseq)-1)
  RCNormMedian<-RC
  for (i in 1:(length(stepseq)-1)){
    if (i==1){
      ind<-which(GCContent>=stepseq[i] & GCContent<=stepseq[i+1])
    }
    if (i!=1){
      ind<-which(GCContent>stepseq[i] & GCContent<=stepseq[i+1])
    }
    if (length(ind)>0){
      m<-median(RC[ind],na.rm=T)
      if (m>0){
        MedianGC[i]<-m
        RCNormMedian[ind]<-RC[ind]*MasterMedian/m
      }
    }
  }
  RCNormList<-list()
  RCNormList$Median<-MedianGC
  RCNormList$StepGC<-stepseq[1:(length(stepseq)-1)]
  RCNormList$RCNorm<-RCNormMedian
  RCNormList
}

################ Correzione dei RC dalla Mappability ####################
CorrectMAP<-function(RC,MAPContent,step){
  stepseq<-seq(min(MAPContent),max(MAPContent),by=step)
  #stepseq<-seq(0,100,by=step)
  MasterMedian<-median(RC,na.rm=T)
  MedianMAP<-rep(0,length(stepseq)-1)
  RCNormMedian<-RC
  for (i in 1:(length(stepseq)-1)){
    if (i==1){
      ind<-which(MAPContent>=stepseq[i] & MAPContent<=stepseq[i+1])
    }
    if (i!=1){
      ind<-which(MAPContent>stepseq[i] & MAPContent<=stepseq[i+1])
    }
    
    if (length(ind)>0){
      m<-median(RC[ind],na.rm=T)
      if (m>0){
        MedianMAP[i]<-m
        RCNormMedian[ind]<-RC[ind]*MasterMedian/m
      }
    }
  }
  RCNormList<-list()
  RCNormList$Median<-MedianMAP
  RCNormList$StepMAP<-stepseq[1:(length(stepseq)-1)]
  RCNormList$RCNorm<-RCNormMedian
  RCNormList
}

# -------------------------------------------------------------------------------------------- #
# below are normalization code for normalization for GC content and mappability
GC.MAP.normalization<-function(Y, #raw matrix
                               ref_qc#info including GC, mappability positions, etc.
){
  
  sampname <- colnames(Y)
  start=ref_qc$start
  end=ref_qc$end
  exon_size=end-start
  gc=ref_qc$gc
  map=ref_qc$mapp  
  
  RCTL <- Y
  RC_norm=matrix(data=NA,nrow = dim(RCTL)[1],ncol = dim(RCTL)[2])
  rownames(RC_norm)=rownames(RCTL)
  colnames(RC_norm)=colnames(RCTL)
  for(sub in 1:dim(RCTL)[2]){
    ### Mappability normalization ###
    step <- 0.01
    RCMAPNormListIn <- CorrectMAP(RCTL[,sub],MAPContent=map,step)
    RCMAPNormIn <- RCMAPNormListIn$RCNorm
    
    ### GC-content Normalization ###
    step <- 0.05
    RCGCNormListIn <- CorrectGC(RCMAPNormIn,GCContent=gc,step)
    RCGCNormIn <- RCGCNormListIn$RCNorm
    
    RC_norm[,sub]=RCGCNormIn
  }
  
  lrr=matrix(data=NA,nrow = dim(RC_norm)[1],ncol = dim(RC_norm)[2])
  rownames(lrr)=rownames(RC_norm)
  colnames(lrr)=colnames(RC_norm)
  for(i in 1:dim(RC_norm)[2]){
    mean=median(RC_norm[,i])
    rc=RC_norm[,i]
    lrr[,i]=log2((rc+0.01)/mean)
  }
  log2Rdata=lrr
  return(list(RC_norm=RC_norm,log2Rdata=log2Rdata))
  
}

############# Library size normalization #############
libsize_normalize <- function(X, method = c("cpm", "prop"),
                              log1p = FALSE, scale = 1e6,
                              eps = 1e-8) {
  method <- match.arg(method)
  X <- as.matrix(X)
  
  lib <- colSums(X, na.rm = TRUE)
  lib_safe <- pmax(lib, eps)
  
  if (method == "prop") scale <- 1
  
  Xn <- sweep(X, 2, lib_safe, "/") * scale
  
  if (log1p) Xn <- log1p(Xn)
  return(list(X_norm = Xn, libsize = lib))
}