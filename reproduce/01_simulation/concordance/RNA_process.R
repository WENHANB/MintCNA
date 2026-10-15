library(data.table)

args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)  

for (p in c(1:3)) {
  
file.path(CNA_DIR,"rna_counts_gene_by_cell_g1_p",p,".tsv")
label <- fread(file.path(CNA_DIR,"cell_types.tsv"))
Y <- fread(file.path(CNA_DIR,paste0("rna_counts_gene_by_cell_g1_p",p,".tsv")))
ref_info_rna <- fread(file.path(CNA_DIR,"genes_hg38.tsv"))

## CPM 
library(edgeR)
rna_read <- as.matrix(Y[,-1])
dge <- DGEList(counts = rna_read)
cpm_matrix <- cpm(dge)
rna_read <- cpm_matrix

# CopyKAT - Preprocess
rawmat <- rna_read

## Step 1
print("step1: read and filter data ...")
print(paste(nrow(rawmat), " genes, ", ncol(rawmat), " cells in raw data", sep=""))
genes.raw <- apply(rawmat, 2, function(x)(sum(x>0)))
if(sum(genes.raw> 200)==0) stop("none cells have more than min.gene.per.cell")
if(sum(genes.raw<200)>1){
  rawmat <- rawmat[, -which(genes.raw< 200)]
  print(paste("filtered out ", sum(genes.raw<=200), " cells with less than min.gene.per.cell; remaining ", ncol(rawmat), " cells", sep=""))
}
##
der<- apply(rawmat,1,function(x)(sum(x>0)))/ncol(rawmat)
if(sum(der>0.05)>=1){
  rawmat <- rawmat[which(der > 0.05), ]; print(paste(nrow(rawmat)," genes preseved after LOW.DR filtering", sep=""))
}
ref <- ref_info_rna[which(der > 0.05),]
WNS1 <- "data quality is ok"
if(nrow(rawmat) < 7000){
  WNS1 <- "low data quality"
  UP.DR<- 0.05
  print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
}

## Step 2 
print("step2: log-Freeman-Tukey transformation ...")
norm.mat<- log(sqrt(rawmat)+sqrt(rawmat+1))
#norm.mat<- apply(norm.mat,2,function(x)(x <- x-mean(x)))

## Step 3
print("step 3: smoothing data with dlm ...")
dlm.sm <- function(c){
  model <- dlm::dlmModPoly(order=1, dV=0.16, dW=0.001)
  x <- dlm::dlmSmooth(norm.mat[, c], model)$s
  x<- x[2:length(x)]
  x <- x-mean(x)
}

test.mc <-parallel::mclapply(1:ncol(norm.mat), dlm.sm, mc.cores = 1)
norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat), byrow = FALSE)
colnames(norm.mat.smooth) <- colnames(norm.mat)
rownames(norm.mat.smooth) <- ref$gene_id

nc_samples <- label[which(label$normal == "normal"), root]
nc_median <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% nc_samples)],1,median)
norm.mat.smooth <- norm.mat.smooth-nc_median
l2r <- norm.mat.smooth[,setdiff(colnames(norm.mat.smooth),nc_samples)]

## Median shift normalization
log2R_matrix <- apply(l2r, 2, function(col) {col - median(col, na.rm = TRUE)})

write.csv(log2R_matrix, file = file.path(CNA_DIR,paste0("l2r_all_cells_CopyKAT_p",p,".csv")))
write.csv(ref, file = file.path(CNA_DIR,paste0("ref_all_cells_CopyKAT_p",p,".csv")))
}

