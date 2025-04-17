# scCNMF
Use "scCNMF" to integrate the paired scRNA-seq data and scATAC-seq data and cluster cells.

## example

### Download
devtools::install_github("Jinsl-lab/scCNMF/scCNMF")
### Read data
element=readRawData("filtered_feature_bc_matrix/")
### Run ccNMF
results <- scCNMFmodel(X1 = element$X1, X2 = element$X2, K=20, s = 0.25, alpha = 1, beta = 1, gamma = 10000, maxIter = 500,
                 GeneName = element$GeneName, PeakName = element$PeakName, GeneLoc = element$GeneLoc, PeakLoc = element$PeakLoc)
### Dimensionality reduction and clustering
ans_umap=ClusteringVisual(as.matrix(results$H),method="umap",celltype=celltype)  
ans_tsne=ClusteringVisual(as.matrix(results$H),method="tsne",celltype=celltype)
### Adding Comments
colnames(results$H)<-element$Barcode$V1  
rownames(results$H)<- paste0("factor", 1:20)  
rownames(results$W1)<-results$Gene_name  
colnames(results$W1)<- paste0("factor", 1:20)  
rownames(results$W2)<-results$peak_name  
colnames(results$W2) <- paste0("factor", 1:20)  
### Plotting
FactorVisual(ans, feature_scores = t(results$H), feature_using = paste0("factor",c(2,5,10,15,19,20)),  
                     method = "umap", nCol = 3, cell_size = 0.1, show_legend = T, show_legend_combined = F)  
feature_genes = c('CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'IL7R', 'CD19', 'CD79A', 
                  'CD79B', 'MS4A1', 'BANK1', 'CR2', 'CD14', 'CD68', 'FCGR3A', 'ITGAM', 
                  'ITGAX', 'C1QA', 'NCAM1', 'KLRC1', 'KLRB1', 'GNLY', 'NKG7', 'CD34', 
                  'CD38', 'CD45', 'CD235A', 'GATA1', 'HBB', 'HBA1', 'HBA2', 'CD11b', 
                  'LYZ', 'CD163', 'CD206', 'CD209', 'CD11c', 'CD1c', 'CD141', 'CLEC9A', 
                  'MS4A4A', 'FCER1A', 'MS4A6A', 'MS4A7', 'CD56', 'CD161', 'CD94', 'NKG7', 
                  'GNLY', 'PRF1')   
FeatureRank(results,assay = 'RNA', factor_show = c(2,5,10,15,19,20) , ncol = 2, feature_show = feature_genes, top = 0.25, ylabel = "Gene score")  
options(future.globals.maxSize = 10*1024 * 1024^2)  
future::plan("multisession", workers = 4)  
markers.RNA.cluster <- ClusterMarkers(element = element,clusting_results = ans, assay = 'RNA')  
n.top = 10  
markers.RNA.clusterTop <- markers.RNA.cluster %>% group_by(clusters) %>% top_n(n.top, logFC) %>% slice(1:n.top)  
FeatureHeatmap(element = element,clusting_results = ans, assay = "RNA", feature_using = markers.RNA.clusterTop$features, class = "cluster")  
