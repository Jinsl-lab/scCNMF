library(Matrix)
library(Rcpp)
library(lsa)
library(umap)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ComplexHeatmap)
library(circlize)



#' Define multiplication between matrices
#' @param A Matrix
#' @param B Matrix
#' @export
eigenMapMatMult <- function(A, B) {
  C <- A %*% B
  return(C)
}

#' Define multiplication between matrices
#'
#' @param A Matrix
#' @param B Matrix
#' @export
eigenMapMatcrossprod <- function(A, B) {
  C <- t(A) %*% B
  return(C)
}

#' Define multiplication between matrices
#'
#' @param A Matrix
#' @param B Matrix
#' @export
eigenMapMattcrossprod <- function(A, B) {
  C <- A %*% t(B)
  return(C)
}

#' Define multiplication between matrices
#'
#' @param A Matrix
#' @export
crossprodcpp <- function(A){
  C <- t(A) %*% A
  return(C)
}

#' Read raw data in 10X format
#' X1: scRNA-seq data
#' X2: scATAC-seq data
#' PeakName: Name of the ATAC peaks
#' Barcode: Cell Name
#' GeneLoc: Location of genes on chromosomes
#' GeneName: Gene name
#' PeakLoc: Location of peaks on chromosomes
#' chr: Chromosomes name
#'
#' @param foldername Data storage address, required to contain features.tsv, barcodes.tsv and matrix.mtx
#' @export
readRawData<-function(foldername){
    filename=paste0(foldername,"matrix.mtx")
    mtx=read.table(filename,sep = "",header = T,comment.char = '%')
    filename=paste0(foldername,"features.tsv")
    feature=read.table(filename,sep='\t')
    if_Peaks=match(feature[,3],"Peaks")
    if_Peaks[is.na(if_Peaks)]=0
    if_Peaks=as.logical(if_Peaks)
    PeakName=feature[if_Peaks,2]
    GeneName=feature[!if_Peaks,2]
    FeatureName=feature[,2]
    chr=unique(feature[,4])
    chr=chr[grep("^chr",chr)]
    chr_idx=match(feature[,4],chr)
    chr_idx[is.na(chr_idx)]=0
    GeneLoc=matrix(chr_idx[!if_Peaks],nrow = length(chr_idx[!if_Peaks]))
    GeneLoc=cbind(GeneLoc,feature[!if_Peaks,5])
    PeakLoc=matrix(chr_idx[if_Peaks],nrow = length(chr_idx[if_Peaks]))
    PeakLoc=cbind(PeakLoc,feature[if_Peaks,5])
    filename=paste0(foldername,"barcodes.tsv")
    Barcode=read.table(filename,sep="\t")
    # RNA
    if_RNA=match(FeatureName,GeneName)
    isRNA=if_RNA;
    isRNA[is.na(isRNA)]=0
    isRNA=as.logical(isRNA)
    exp=match(mtx[,1],which(isRNA==TRUE))
    exp[is.na(exp)]=0
    exp=as.logical(exp)
    RNAData=mtx[exp,]
    RNAData[,1]=if_RNA[RNAData[,1]]
    X1=sparseMatrix(RNAData[,1],RNAData[,2],x=log2(1+RNAData[,3]),dims = c(length(GeneName),length(Barcode[,1])))
    # ATAC
    if_ATAC=match(FeatureName,PeakName)
    isATAC=if_ATAC;
    isATAC[is.na(isATAC)]=0
    isATAC=as.logical(isATAC)
    exp=match(mtx[,1],which(isATAC==TRUE))
    exp[is.na(exp)]=0
    exp=as.logical(exp)
    ATACData=mtx[exp,]
    ATACData[,1]=if_ATAC[ATACData[,1]]
    X2=sparseMatrix(ATACData[,1],ATACData[,2],x=log10(1+ATACData[,3]),dims = c(length(PeakName),length(Barcode[,1])))
    # Data preprocessing
    tag=GeneLoc[,1]>0
    GeneName=GeneName[tag];
    GeneLoc=GeneLoc[tag,]
    X1=X1[tag,]
    tag=PeakLoc[,1]>0
    PeakName=PeakName[tag]
    PeakLoc=PeakLoc[tag,]
    X2=X2[tag,]
    tag=rowSums(as.matrix(X1)>0 )>0;
    GeneName=GeneName[tag];
    GeneLoc=GeneLoc[tag,]
    X1=X1[tag,]
    tag=rowSums(as.matrix(X2)>0)>10;
    PeakName=PeakName[tag]
    PeakLoc=PeakLoc[tag,]
    X2=X2[tag,]
    rownames(X1)<-GeneName
    colnames(X1)<-Barcode$V1
    rownames(X2)<-PeakName
    colnames(X2)<-Barcode$V1
    return(list(X1=X1,
                X2=X2,
                PeakName=PeakName,
                Barcode=Barcode,
                GeneLoc=GeneLoc,
                GeneName=GeneName,
                PeakLoc=PeakLoc,
                chr=chr))
}

#' Normalize the input matrix
#'
#' @param X Matrix
#' @export
Normlize<- function (X){
  X2=as.matrix(X);
  X2[X2==0]=NaN;
  b=apply(X2,2,median,na.rm=TRUE);
  return(as(t(t(as.matrix(X))/b*mean(b)),"sparseMatrix"));
}

#' Z-score initialization
#'
#' @param x Matrix
#' @export
z_score_normalize <- function(x) {
  apply(x, 2, function(col) (col - mean(col)) / sd(col))
}

#' Fold_RE_TG_multiAdjust
#'
#' @param X1 scRNA-seq data matrix
#' @param X2 scATAC-seq data matrix
#' @param GeneLoc Location of genes on chromosomes
#' @param PeakLoc Location of peaks on chromosomes
#'
#' @export
Fold_RE_TG_multiAdjust<- function (X1,X2,GeneLoc,PeakLoc){
  X1=as.matrix(X1);
  X2=as.matrix(X2);
  P1=Fold_RE_TG_MultiAdjustCore(X1, X2, GeneLoc, PeakLoc);
  return(P1);
}

#' Perform TF-IDF conversion operation
#'
#' @param X Data matrix
#' @export
TFIDF_trans <-function(X){
  X=as.matrix((X>0));
  tf1=t(t(X)/log(1+colSums(X)));
  idf=log(1+dim(X)[2]/(1+rowSums(X>0)));
  X1=tf1*idf;
  X1[is.na(X1)] <- 0;
  return(as(X1,"sparseMatrix"))
}

#' Generate an initialized jaccard similarity matrix
#'
#' @param matrixmtx Data matrix
#' @export
jaccard_matrix<-function(matrixmtx){
a=read.table(matrixmtx,sep = "",header = T,comment.char = '%')
  mean_value <- mean(a[[3]])
  binary_matrix <- X1
  binary_matrix@x[binary_matrix@x < mean_value] <- 0
  binary_matrix@x[binary_matrix@x > mean_value] <- 1
  binary_matrix_dense <- as.matrix(binary_matrix)
  jaccard_matrix <- as.matrix(proxy::dist(t(binary_matrix_dense), method = "Jaccard"))
  jaccard_similarity <- 1 - jaccard_matrix
  return(jaccard_similarity)
}

#' Integration operations are performed on the two matrices X1,X2 obtained after initialization, and the rcpp program is called in the middle.
#' Function 'ccNMF_cpp_function' calls c++ program 'ccNMF.cpp'
#'
#' @param X1 scRNA-seq data matrix
#' @param X2 scATAC-seq data matrix
#' @param s Probability of Bernoulli distribution
#' @param K The rank of the latent factor
#' @param alpha Model parameter
#' @param beta Model parameter
#' @param gamma Model parameter
#' @param maxIter Maximum number of iterations
#' @param stop_rule Iteration stopping conditions, taking a fixed number of iterations when equal to 1 and adaptive iterations when equal to 2
#' @param GeneName Name of genes
#' @param PeakName Name of the ATAC peaks
#' @param GeneLoc Location of genes on chromosomes
#' @param PeakLoc Location of peaks on chromosomes
#' @param feature_cut_perc Feature cut percent
#' @param corr_cut_k Top peak-gene pairs
#' @param core Core
#' @export
ccNMFmodel<- function(X1, X2,
                   s=0.25,K=20, alpha = 1, beta = 1, gamma = 10000, maxIter=500, stop_rule=1,
                   GeneName,PeakName,GeneLoc,PeakLoc,
                   feature_cut_perc=0.01,corr_cut_k=100000,core=8) {
  X11=Normlize(X1);
  X21=Normlize(X2);
  numCut=feature_cut_perc*dim(X11)[2];
  a=Matrix::rowSums(X11>0);
  a1=Matrix::rowSums(X21>0);
  X111=X11[a>numCut,];
  X211=X21[a1>numCut,];
  GeneLoc=GeneLoc[a>numCut,];
  PeakLoc=PeakLoc[a1>numCut,];
  GeneName=GeneName[a>numCut];
  PeakName=PeakName[a1>numCut];
  R=Fold_RE_TG_multiAdjust(X111,X211,GeneLoc,PeakLoc)
  R=as(R,"sparseMatrix");
  mk=sort(R[R>0],decreasing = TRUE);
  corr_cut_k=min(length(mk),corr_cut_k)
  mk=mk[1:corr_cut_k];
  corr_cut=mk[length(mk)];
  a=Matrix::colSums(R>=corr_cut);
  X211=X211[a>0,];
  a2=Matrix::rowSums(R>=corr_cut);
  X111=X111[a2>0,];
  R=R[a2>0,a>0];
  GeneLoc=GeneLoc[a2>0,];
  PeakLoc=PeakLoc[a>0,]
  GeneName=GeneName[a2>0];
  PeakName=PeakName[a>0];
  c=which(((as.matrix(R)>=corr_cut)==1),arr.ind = TRUE);
  Reg=X211[c[,2],]+X111[c[,1],];
  Reg_gene_location=GeneLoc[c[,1],];
  Reg_peak_location=PeakLoc[c[,2],];
  Reg_dis=abs(Reg_gene_location[,2]-Reg_peak_location[,2]);
  d0=2*10^5;
  Reg_w=exp(-Reg_dis/d0);
  Reg_adj=Reg*Reg_w;
  X212=TFIDF_trans(X211);
  Reg_info=cbind(c[,1],c[,2]);
  Reg_info=cbind(Reg_info,Reg_dis);
  X1<-X111
  X2<-X212
  W10=matrix(runif(dim(X1)[1]*K,min = 0,max=1),dim(X1)[1],K);
  W20=matrix(runif(dim(X2)[1]*K,min = 0,max=1),dim(X2)[1],K);
  W30=matrix(runif(dim(Reg_adj)[1]*K,min = 0,max=1),dim(Reg_adj)[1],K);
  H0=matrix(runif(dim(X1)[2]*K,min=0,max=1),K,dim(X1)[2]);
  # θ0=jaccard_matrix("filtered_feature_bc_matrix/matrix.mtx")
  θ0=matrix(runif(dim(X1)[2]*dim(X1)[2],min=0,max=1),dim(X1)[2],dim(X1)[2])
  O0=matrix(as.double(rbinom(dim(X1)[2]*dim(X1)[2],1,s)),dim(X1)[2],dim(X1)[2])
  object=ccNMF_cpp_function(as.matrix(X1),as.matrix(X2),as.matrix(Reg_adj),K,maxIter,stop_rule,alpha,beta,gamma,W10,W20,W30,H0,θ0,O0,c[,1],c[,2],Reg_w,core);
  object[["Reg_gene_name"]]=GeneName[c[,1]];
  object[["Reg_peak_name"]]=PeakName[c[,2]];
  object[["Gene_name"]]=GeneName;
  object[["peak_name"]]=PeakName;
  object[["Reg_adj"]]=Reg_adj;
  result <- list(W1 = object$W1, W2 = object$W2, W3= object$W3, H = object$H, θ = object$θ, Reg_gene_name = object$Reg_gene_name, Reg_peak_name = object$Reg_peak_name, Gene_name = object$Gene_name, peak_name = object$peak_name, R = object$Reg_adj,
               options = list(s = s, alpha = alpha, beta = beta, gamma = gamma, maxIter = maxIter))
  return(result);
}

#' Calculation of the Jaccard factor
#'
#' @param X Matrix
#'
#' @export
Jaccard_result <- function(X) {
  M <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  for (i in 1:nrow(X)) {
    for (j in 1:ncol(X)) {
      M[i, j] <- X[i, j] / (X[i, i] + X[j, j] - X[i, j])
    }
  }
  return(M)
}

#' Realizing cell clustering
#'
#' @param H Matrix for integrating results
#' @param k The number of nearest neighbors chosen by the KNN
#'
#' @export
Simbased_clustering<- function(H,k=20){
  sim=cosine(H);
  Sim=colSums(sim);
  twom=sum(Sim);
  sim_norm=sim-(as.matrix(Sim) %*% (t(as.matrix(Sim/twom))))
  f=matrix(0,dim(sim_norm)[1],dim(sim_norm)[2]);
  for (i in (1:dim(sim_norm)[2])) {
    f[,i]=order(sim_norm[,i],decreasing = TRUE);
  }
  KNN=matrix(0,dim(sim_norm)[2],dim(sim_norm)[2]);
  for (i in (1:dim(sim_norm)[2])) {
    KNN[f[2:(k+1),i],i]=1;
  }
  KNN2=t(KNN)%*%KNN;
  A=Jaccard_result(KNN2);
  A[is.na(A)] <- 0;
  Sim=colSums(as.matrix(A));
  twom=sum(Sim);
  B=A-((as.matrix(Sim)%*%t(as.matrix(Sim)))/twom);
  a=neighbor_default(B);
  return(a);
}

#' Construct a new network matrix by transforming the input matrix using sparse matrices
#'
#' @param A Matrix
#' @param B Index vector
#'
#' @export
metanet <- function(A,B){
  X <- sparseMatrix(i=1:length(B),j=B,x=1)
  t(X) %*% A %*% X
}

#' Grouping operations
#'
#' @param X Vector
#'
#' @export
tidyresult <- function(X){
  Y <- rep(0,length(X))
  for (i in 1:length(X))
  {
    if (Y[i]==0) Y[X==X[i]] <- max(Y)+1
  }
  Y
}

#' Dividing nodes in the network into communities makes connections within communities tight and connections between communities sparse
#' The algorithm keeps merging nodes with similar connectivity patterns in an iterative manner until it reaches the optimal community segmentation result
#'
#' @param X Matrix
#' @export
neighbor_default<-function(X){
  eps <- 2e-16
  if (sum(X-t(X))>0) X <- (X+t(X))/2  #force symmetric matrix
  M <- X
  n <- nrow(X)
  S <- t(1:n)
  dtot <- 0
  S2 <- 1:n
  Sb <- NULL
  n.outer <- 0
  while(!identical(Sb,S2)){
    n.outer <- n.outer+1
    y <- unique(S2)
    y <- y[order(y,decreasing=F)]
    Sb <- S2
    print(paste(c("Merging",length(y),"communities"),collapse=" "))
    yb <- NULL
    G <- sparseMatrix(i=1:length(y),j=y,x=1)
    dstep <- 1
    nsteps <- 0
    while((!identical(yb,y)) && (dstep/dtot>2*eps)){
      yb <- y
      dstep <- 0
      nsteps <- nsteps+1
      ord.i <- 1:nrow(M)
      for (i in ord.i)  {
        u <- unique(c(y[i],y[M[,i]>0]))
        u <- u[order(u,decreasing=F)]
        dH <- t(M[,i]) %*% G[,u]
        yi <- which(u==y[i])
        dH[yi] <- dH[yi] - M[i,i]
        k <- max.col(dH)
        if (dH[k]>dH[yi]){
          dtot <- dtot+dH[k]-dH[yi]
          dstep <- dstep+dH[k]-dH[yi]
          G[i,y[i]] <- 0
          G[i,u[k]] <- 1
          y[i] <- u[k]
        }
      }
    }
    y <- tidyresult(y)
    for (i in 1:length(y)){
      S[S==i] <- y[i]
      S2[S2==i] <- y[i]
    }
    if (identical(Sb,S2)){
      return(S)
    }
    M <- metanet(X,S2)
  }
}

#' Implementing dimensionality reduction clustering and visualizing low-dimensional plots
#'
#' @param H Matrix for integrating results
#' @param cell_num_smooth Calculation parameters
#' @param celltype Celltype
#' @param method Dimensionality reduction method
#'
#' @export
clustering<-function(H,cell_num_smooth=sqrt(dim(H)[2]),celltype,method="umap"){
  HH=t(t(H)/colSums(H));
  S=matrix(ncol = 3);
  S=na.omit(S);
  S_tmp=Simbased_clustering(HH,cell_num_smooth);
  if(method =="umap"){
    print("Run UMAP")
    umap_result <- umap(t(HH))
    ans <- data.frame(x = umap_result$layout[, 1], y = umap_result$layout[, 2], col = S_tmp[1,])
  }
  else{
    print("Run t-SNE")
    T=Rtsne::Rtsne(t(HH));
    ans <- data.frame(x=T$Y[,1],y=T$Y[,2],col=S_tmp[1,]);
  }
  # s<-celltype
  s<-as.character(ans$col)
  colourCount = length(unique(s))
  getPalette = colorRampPalette(c("#f8766d","#e68613","#cd9600","#aba300","#7cae00","#0cb702","#00be67","#00c19a",
                                  "#00bfc4","#00b8e7","#00a9ff","#8494ff","#c77cff","#ed68ed","#ff61cc","#ff68a1"))
  plot<-ggplot(ans,aes(x=x,y=y,colour=factor(s)))+
  geom_point(shape=16,size=0.2,aes(col = s))+
  #scale_colour_viridis(option="turbo",discrete=TRUE)+
  scale_colour_manual(values = getPalette(colourCount))+
  labs(x = "UMAP_1", y = "UMAP_2") +  # ggtitle("Liger_human_CellType")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_test()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),legend.title = element_blank())+
  theme(legend.position = "right" ,legend.box = "horizontal")
  return(list(HH=HH,
              S=S_tmp,
              ans=ans,
              plot=plot))
}


#' Drawing format
#'
#' @export
theme_opts <- function() {
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme_classic() +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(color = "black")) +
    theme(axis.line.y = element_line(color = "black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(legend.key = element_blank()) + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
}

#' Custom color palettes
#'
#' @param n Num
#' @export
scPalette <- function(n) {
    colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC',
                    '#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A','#E3BE00','#FB9A99',
                    '#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B',
                    '#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
    if (n <= length(colorSpace)) {
        colors <- colorSpace[1:n]
    } else {
        colors <- grDevices::colorRampPalette(colorSpace)(n)
    }
    return(colors)
}

#' Drawing feature maps to visualize cells in 2D space with superimposed features
#'
#' @param clusting_result UMAP/t-SNE results
#' @param feature_using Eigenvector
#' @param feature_scores Matrix containing feature scores
#' @param method Dimensionality reduction methods, e.g. tsne, umap
#' @param colormap RColorbrewer color palette to use
#' @param color_direction Sets the order of the colors in the scale. If 1 (default), the colors will be in the order of the output of RColorBrewer::brewer.pal(). If -1, the colors are in reverse order
#' @param nCol Drawing columns
#' @param show.axes Whether to display axes
#' @param cell_size Size of the point
#' @param show_legend Whether to display a single legend
#' @param show_legend_combined Whether to display only one legend
#' @export
feature_score_visual <- function(clusting_result, feature_using = NULL, feature_scores, method = "umap",
                                      colormap = "GnBu", color_direction = 1,nCol = NULL,
                                      show.axes = T,  cell_size = 0.3,
                                      show_legend = T, show_legend_combined = F) {
  data_using <- as.matrix(feature_scores[ ,feature_using])
  if (is.null(nCol)) {
    if (length(feature_using) > 9) {
      nCol <- 4
    } else {
      nCol <- min(length(feature_using), 3)
    }
  }
  if (method == "tsne") {
    xlabel = "t-SNE_1"
    ylabel = "t-SNE_2"
  } else if (method == "umap") {
    xlabel = "UMAP_1"
    ylabel = "UMAP_2"
  }
  df <- data.frame(x = clusting_result$ans$x, y = clusting_result$ans$y)
  numFeature = length(feature_using)
  gg <- vector("list", numFeature)
  for (i in seq_len(numFeature)) {
    feature.name <- feature_using[i]
    df$feature.data <- data_using[ ,i]
    g <- ggplot(data = df, aes(x, y)) +
      geom_point(aes(colour = feature.data), size = cell_size) +
      scale_color_distiller(palette = colormap, direction = color_direction, guide = guide_colorbar(title = NULL, ticks = T, label = T, barwidth = 0.5), na.value = "lightgrey") +
      labs(title = feature.name) + theme_opts() +
      theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in")) + labs(x = xlabel, y = ylabel)
    if (!show_legend) {
      g <- g + theme(legend.position = "none")
    }
    if (show_legend_combined & i == numFeature) {
      g <- g + theme(legend.position = "right", legend.key.height = grid::unit(0.15, "in"), legend.key.width = grid::unit(0.5, "in"), legend.title = NULL)
    }
    if (!show.axes) {
      g <- g + theme_void()
    }
    gg[[i]] <- g
  }
  gg.combined <- cowplot::plot_grid(plotlist = gg, ncol = nCol)
  gg.combined
}

#' Sort features (genes/genomes) and show top markers in each factor
#'
#' @param results UMAP/t-SNE results
#' @param assay Define the assay to be displayed, e.g., assay = “RNA”
#' @param factor_show The set of factors shown
#' @param ncol Drawing columns
#' @param feature_show Eigenvector
#' @param top Show top-ranked features
#' @param features_diff Differential characteristics
#' @param ylabel Labels
#' @param n Default display of the first n genes
#' @export
feature_rank <- function(results, assay, factor_show = NULL, ncol = NULL, feature_show = NULL, top = 0.5, features_diff = NULL, ylabel = "Weight",n=10) {
  if (assay == "RNA") {
    W <- results$W1
    rownames(W)<-rownames(results$W1)
    colnames(W)<-colnames(results$W1)
  } else if (assay == "ATAC") {
    W <- results$W2
    rownames(W)<-rownames(results$W2)
    colnames(W)<-colnames(results$W2)
  }
  features <- rownames(W)
  if (!is.null(factor_show)) {
    W <- W[, factor_show]
  }
  K = ncol(W)
  W <- sweep(W,1,rowSums(W),FUN = `/`)
  W[is.na(W)] <- 0
  Wg <- vector("list", K)
  for (i in 1:K) {
    W_order <- sort(W[,i],decreasing=F, index.return = T)
    features_ordered <- features[W_order$ix]
    if (!is.null(features_diff)) {
      features_diffi <- as.character(features_diff$features[features_diff$factors == i])
    }else {
      features_diffi <- as.character(features)
    }
    if (!is.null(feature_show)) {
      features_diffi <- intersect(features_diffi, feature_show)
    }
    idx <- match(features_diffi, features_ordered)
    data_show <- matrix(0, nrow(W), 1); data_show[idx] <- 1
    if (!is.null(top) & top < 1) {
      idx_bottom <- seq_len(floor((1-top)*nrow(W))); data_show[idx_bottom] <- 0
    }
    Wg[[i]] <- cbind(Weight =  as.numeric(W_order$x), factor = colnames(W)[i], Ranking = seq_len(nrow(W)), Show = as.numeric(data_show), Genes = features_ordered)
  }
  data <- Wg[[1]]
  for (i in 2:K) {
    data <- rbind(data, Wg[[i]])
  }
  df <- as.data.frame(data, stringsAsFactors=FALSE)
  colnames(df) <- c("Weight", "factor", "Ranking", "Show","Genes")
  df$Weight <- as.numeric(as.character(df$Weight))
  df$Ranking <- as.numeric(as.character(df$Ranking))
  df$Show <- as.numeric(as.character(df$Show))
  data_topFeature <- df %>%
    filter(!is.na(Genes) & Genes != "") %>%
    group_by(factor) %>%
    slice_max(Weight, n = n, with_ties = FALSE) %>%
    ungroup()
  print(table(data_topFeature$factor))
  factor_order <- paste0("factor", factor_show)
  df$factor <- factor(df$factor, levels = factor_order)
  data_topFeature$factor <- factor(data_topFeature$factor, levels = factor_order)
  gg <- ggplot(df, aes(Ranking, Weight)) +
    geom_line(colour = "grey80", linewidth = 1) +
    facet_wrap(~ factor, ncol = ncol, scales = "free") +
    theme_opts() +
    theme(text = element_text(size = 8),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14)) +
    theme(strip.background = element_rect(fill="grey80")) +
    ylab(ylabel) +
    geom_point(size = 3, shape = 1, data = data_topFeature) +
    ggrepel::geom_text_repel(aes(label = Genes),
                             data = data_topFeature,
                             segment.color = "grey50",
                             segment.alpha = 1,
                             direction = "y",
                             nudge_x = -300,
                             hjust = 1,
                             size = 4,
                             segment.size = 0.5,
                             max.overlaps = n)
  return(gg)
}

#' Identify markers in each cell population
#'
#' @param element Eigenexpression matrix
#' @param clusting_results UMAP/t-SNE results
#' @param assay Define the assay to be displayed, e.g., assay = “RNA”
#' @param features Eigenvector
#' @param test_using Test to be used (“bimod” or “wilcox”)
#' @param thresh_pc Threshold
#' @param thresh_fc Threshold
#' @param thresh_p Threshold
#' @export
cluster_markers <- function(element,clusting_results, assay, features = NULL, test_using = "bimod",
                                   thresh_pc = 0.25, thresh_fc = 0.25, thresh_p = 0.01) {
  if (assay == "RNA") {
    X <- element$X1
  } else {
      X <- element$X2
    }
  if (is.null(features)) {
    features_using <- row.names(X)
  } else {
    features_using <- intersect(features, row.names(X))
  }
  data_using <- X[features_using,]
  data_using <- as.matrix(data_using)
  groups <- factor(clusting_results$ans$col)
  labels <- groups
  level_using <- levels(labels)
  numCluster <- length(level_using)
  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )
  mean.fxn <- function(x) {
    return(log(x = mean(x = expm1(x = x)) + 1))
  }
  genes.de <- vector("list", length = numCluster)
  for (i in 1:numCluster) {
    features <- features_using
    cell_using1 <- which(labels %in% level_using[i])
    cell_using2 <- base::setdiff(1:length(labels), cell_using1)
    thresh_min <- 0
    pct.1 <- round(
      x = rowSums(x = data_using[features, cell_using1, drop = FALSE] > thresh_min) /
        length(x = cell_using1),
      digits = 3
    )
    pct.2 <- round(
      x = rowSums(x = data_using[features, cell_using2, drop = FALSE] > thresh_min) /
        length(x = cell_using2),
      digits = 3
    )
    data.alpha <- cbind(pct.1, pct.2)
    colnames(x = data.alpha) <- c("pct.1", "pct.2")
    alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
    names(x = alpha.min) <- rownames(x = data.alpha)
    features <- names(x = which(x = alpha.min > thresh_pc))
    if (length(x = features) == 0) {
      stop("No features pass thresh_pc threshold")
    }
    data.1 <- apply(X = data_using[features, cell_using1, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
    data.2 <- apply(X = data_using[features, cell_using2, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
    FC <- abs(data.1 - data.2)
    features_diff <- names(x = which(x = FC > thresh_fc))
    features <- intersect(x = features, y = features_diff)
    if (length(x = features) == 0) {
      stop("No features pass thresh_fc threshold")
    }
    data1 <- data_using[features, cell_using1, drop = FALSE]
    data2 <- data_using[features, cell_using2, drop = FALSE]
    if (test_using == "wilcox") {
      pvalues <- unlist(
        x = my.sapply(
          X = 1:nrow(x = data1),
          FUN = function(x) {
            return(wilcox.test(data1[x, ], data2[x, ], alternative = "greater")$p.value)
          }
        )
      )
    } else if (test_using == 'bimod') {
      pvalues <- unlist(
        x = my.sapply(
          X = 1:nrow(x = data1),
          FUN = function(x) {
            return(DifferentialLRT(
              x = as.numeric(x = data1[x,]),
              y = as.numeric(x = data2[x,])
            ))
          }
        )
      )
    }
    pval.adj = stats::p.adjust(
      p = pvalues,
      method = "bonferroni",
      n = nrow(X)
    )
    genes.de[[i]] <- data.frame(clusters = level_using[i], features = as.character(rownames(data1)), pvalues = pvalues, pval_adj = pval.adj, logFC = FC[features], data.alpha[features,])
  }
  markers.all <- data.frame()
  for (i in 1:numCluster) {
    gde <- genes.de[[i]]
    gde <- gde[order(gde$pvalues, -gde$logFC), ]
    gde <- subset(gde, subset = pvalues < thresh_p)
    if (nrow(gde) > 0) {
      markers.all <- rbind(markers.all, gde)
    }
  }
  markers.all$features <- as.character(markers.all$features)
  print(mean(pval.adj))
  return(markers.all)
}


#' Likelihood ratio test
#'
#' @param x Data
#' @param y Data
#' @param xmin parameter
#' @export
DifferentialLRT <- function(x, y, xmin = 0) {
  lrtX <- bimoddata(x = x)
  lrtY <- bimoddata(x = y)
  lrtZ <- bimoddata(x = c(x, y))
  lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
  return(pchisq(q = lrt_diff, df = 3, lower.tail = F))
}

#' Calculate the log-likelihood
#'
#' @param x Data
#' @param xmin Parameter
#'
#' @export
bimoddata <- function(x, xmin = 0) {
  x1 <- x[x <= xmin]
  x2 <- x[x > xmin]
  x11 <- length(x = x2) / length(x = x)
  x11[x11 > 1 - 1e-5] <- 1 - 1e-5
  x11[x11 < 1e-5] <- 1e-5
  likA <- length(x = x1) * log(x = 1 - x11)
  if (length(x = x2) < 2) {
    mysd <- 1
  } else {
    mysd <- sd(x = x2)
  }
  likB <- length(x = x2) *
    log(x = x11) +
    sum(dnorm(x = x2, mean = mean(x = x2), sd = mysd, log = TRUE))
  return(likA + likB)
}

#' Generation of heat maps of differentially characterized expression of different cell populations
#'
#' @param element Eigenexpression matrix
#' @param clusting_results UMAP/t-SNE results
#' @param assay Define the assay to be displayed, e.g., assay = “RNA”
#' @param feature_using Eigenvector
#' @param class Basis for clustering
#' @param color Color
#' @export
feature_heatmap <- function(element,clusting_results, assay, feature_using, class="celltype", color = NULL) {
  if (assay == "RNA") {
    data <- element$X1
  } else {
      data <- element$X2
    }
  if (class == "cluster") {
    groups=factor(clusting_results$ans$col)
  } else {
    groups=clusting_results$celltype
  }
  feature_using <- feature_using[feature_using %in% rownames(data)]
  data_using <- data[feature_using,]
    data_using = Matrix::t(scale(Matrix::t(data_using), center = T))
  data_using <- z_score_normalize(as.matrix(data_using))
  cell_order <- order(groups)
  data_using <- data_using[,cell_order]
  numCluster <- length(unique(groups))
  if (is.null(color)) {
    color <- scPalette(numCluster)
  }
  colorGate = structure(color, names = as.character(levels(groups)))
  col_annotation = HeatmapAnnotation(group = sort(groups),col = list(group = colorGate),
                                     annotation_name_side = "left",simple_anno_size = unit(0.2, "cm"))
  Heatmap(data_using,name = "zscore",
          col = colorRamp2(c(-2,-1,0,1,2), c("#A279C7", "#6F88C2", "#6FA2D5", "#79C5D2", "#F2E273"),space = "LAB"),
          cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = FALSE,
          show_row_dend = FALSE,
          show_row_names = FALSE, row_names_side = "left", row_names_rot = 0,row_names_gp = gpar(fontsize = 6),
          width = unit(10, "cm"),
          bottom_annotation = col_annotation,
          heatmap_legend_param = list(title = NULL, legend_width = unit(0.0001, "cm"),labels_gp = gpar(font = 6))
  )
}




