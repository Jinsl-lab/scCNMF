# scCNMF
Use "scCNMF" to integrate the paired scRNA-seq data and scATAC-seq data and cluster cells.

## Example

### Download
```
devtools::install_github("Jinsl-lab/scCNMF/scCNMF")  
library(scCNMF)
```
### Read data
The scCNMF inputs have the same format as the 10x multiome data and can be downloaded from [10x Genomics](https://www.10xgenomics.com/datasets?configure%5BhitsPerPage%5D=50&configure%5BmaxValuesPerFacet%5D=1000&refinementList%5Bplatform%5D%5B0%5D=Chromium%20Single%20Cell&refinementList%5Bproduct.name%5D%5B0%5D=Epi%20Multiome%20ATAC%20%2B%20Gene%20Expression). Downloads the “Filtered feature barcode matrix MEX (DIR)” file for reading. For other sequencing formats, such as from the GEO database, the data should be converted to a 10x-like format as input. Here, we provide preprocessed PBMC data as a demonstration example, which can be downloaded from [Jinsl-lab/scCNMF/data](https://github.com/Jinsl-lab/scCNMF/tree/main/data).
#### Read raw data and preprocess it
```
element=readRawData("filtered_feature_bc_matrix/")
```
#### For the preprocessed example data we provide, it is possible to read it directly
```
element=readRDS("pbmcdemo.rds")
```

### Run scCNMF
The preprocessed data is input into the scCNMF framework for integration and the integration results are obtained.
```
results <- scCNMFmodel(X1 = element$X1, X2 = element$X2, K=20, s = 0.25, alpha = 1, beta = 1, gamma = 10000, maxIter = 500,
                 GeneName = element$GeneName, PeakName = element$PeakName, GeneLoc = element$GeneLoc, PeakLoc = element$PeakLoc)
```
### Dimensionality reduction and clustering
The integration results were subjected to UMAP and t-SNE downscaling, respectively, and when the real cell labels were missing, the cell clusters obtained by scCNMF could be used as substitutes.
```
celltype<-element$celltype
ans_umap=ClusteringVisual(as.matrix(results$H),method="umap",celltype=celltype)   
ans_tsne=ClusteringVisual(as.matrix(results$H),method="tsne",celltype=celltype)
```
### Adding Comments
```
colnames(results$H)<-element$Barcode  
rownames(results$H)<- paste0("factor", 1:20)  
rownames(results$W1)<-results$Gene_name  
colnames(results$W1)<- paste0("factor", 1:20)  
rownames(results$W2)<-results$peak_name  
colnames(results$W2) <- paste0("factor", 1:20)
```
### Plotting
```
FactorVisual(ans_umap, feature_scores = t(results$H), feature_using = paste0("factor",c(2,5,10,15,19,20)),  
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
markers.RNA.cluster <- ClusterMarkers(element = element,clusting_results = ans_umap, assay = 'RNA')  
n.top = 10  
markers.RNA.clusterTop <- markers.RNA.cluster %>% group_by(clusters) %>% top_n(n.top, logFC) %>% slice(1:n.top)  
FeatureHeatmap(element = element,clusting_results = ans_umap, assay = "RNA", feature_using = markers.RNA.clusterTop$features, class = "cluster")
```
