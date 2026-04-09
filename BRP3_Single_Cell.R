
#single cell expression plot generation by seurat object


install.packages("remotes")
remotes::install_version("ggplot2", version = "3.4.4", repos = "http://cran.us.r-project.org")

library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(tibble)
library(RColorBrewer)
library(scCustomize)
library(viridis)
library(Nebulosa)
library(Matrix)

# Co-expression of Brinp3 and Pcsk2 in mouse hypothalamus
# Data were collected from https://datasets.cellxgene.cziscience.com/d3be7423-d664-4913-89a9-a506cae4c28f.h5ad
seurat_object <- readRDS("mouse_hypothalamus.rds")
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- ScaleData(seurat_object, features = c("ENSMUSG00000035131", "ENSMUSG00000027419"))
#For Brinp3 expressoin
df <- FetchData(seurat_object, vars = c("UMAP_1", "UMAP_2", "ENSMUSG00000035131"))
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = ENSMUSG00000035131)) + geom_point(size = 1) + scale_color_viridis_c(option = "magma", name = " ") + theme_minimal()+ labs(title = "Brinp3 expression") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))
#For Pcsk2 expression
df <- FetchData(seurat_object, vars = c("UMAP_1", "UMAP_2", "ENSMUSG00000027419"))
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = ENSMUSG00000027419)) + geom_point(size = 1) + scale_color_viridis_c(option = "magma", name = " ") + theme_minimal()+ labs(title = "Pcsk2 expression") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))
#For co-expression plot
df <- FetchData(seurat_object, vars = c("UMAP_1", "UMAP_2", "ENSMUSG00000035131", "ENSMUSG00000027419"))
df$coexpression_score <- (df$ENSMUSG00000035131 + df$ENSMUSG00000027419) / 2
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = coexpression_score)) + geom_point(size = 1) + scale_color_viridis_c(option = "magma", name = " ") + theme_minimal() + labs(title = "Co-expression of Brinp3 and Pcsk2") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))

##Dot plot mouse hypothalamus
seurat_object <- readRDS("mouse_hypothalamus.rds") #Mouse hypothalamus
DotPlot_scCustom(seurat_object, features = c("ENSMUSG00000035131", "ENSMUSG00000027419"), group.by = "region", colors_use = rev(viridis_plasma_dark_high))

##DOT PLOT for Brinp3+ cells in mouse VMH
seurat_object <- readRDS("mouse_hypothalamus.rds")#Mouse hypothalamus
seurat_object <- NormalizeData(seurat_object, 
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)

vmh_cells <- subset(seurat_object, subset = Region_summarized == 
                      "Ventromedial hypothalamic nucleus" 
                    & cell_type == "neuron")

vmh_cells$Brinp3_status <- ifelse(vmh_cells@assays$RNA@data["ENSMUSG00000035131", ] > 0, 
                                  "Brinp3+", "Brinp3-")
table(vmh_cells$Brinp3_status)
vmh_cells$only_Brinp3_positive <- ifelse(vmh_cells$Brinp3_status == "Brinp3+" ,"Brinp3+","Other")
Brinp3_positive <- subset(vmh_cells, subset = only_Brinp3_positive == "Brinp3+")

seurat_object <- ScaleData(seurat_object, 
                           features = c("ENSMUSG00000035131", "ENSMUSG00000027419", "ENSMUSG00000048482",
                                        "ENSMUSG00000049796", "ENSMUSG00000019890", 
                                        "ENSMUSG00000027400", "ENSMUSG00000021647", 
                                        "ENSMUSG00000057722", "ENSMUSG00000005534" ,
                                        "ENSMUSG00000024027", "ENSMUSG00000055254", "ENSMUSG00000019768", 
                                        "ENSMUSG00000030500", "ENSMUSG00000026751"))
table(vmh_cells$Brinp3_status)
DotPlot_scCustom(
  vmh_cells,
  features = c("ENSMUSG00000035131", "ENSMUSG00000027419", "ENSMUSG00000048482",
               "ENSMUSG00000049796", "ENSMUSG00000019890", 
               "ENSMUSG00000027400", "ENSMUSG00000021647", 
               "ENSMUSG00000057722", "ENSMUSG00000005534" ,
               "ENSMUSG00000024027", "ENSMUSG00000055254", "ENSMUSG00000019768", 
               "ENSMUSG00000030500", "ENSMUSG00000026751"),
  group.by = "Brinp3_status", dot.scale = 6,
  colors_use = rev(viridis_plasma_dark_high)) + 
  coord_flip()


# Co-expression of BRINP3 and PCSK2 in human hypothalamus
# https://datasets.cellxgene.cziscience.com/480e89e7-84ad-4fa8-adc3-f7c562a77a78.h5ad
seurat_object <- readRDS("human_hypothalamus.rds")
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- ScaleData(seurat_object, features = c("ENSG00000162670", "ENSG00000125851"))
#For BRINP3 expression
df <- FetchData(seurat_object, vars = c("UMAP_1", "UMAP_2", "ENSG00000162670"))
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = ENSG00000162670)) + geom_point(size = 1) + scale_color_viridis_c(option = "magma", name = " ") + theme_minimal()+ labs(title = "BRINP3 expression") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))
#For PCSK2 expression
df <- FetchData(seurat_object, vars = c("UMAP_1", "UMAP_2", "ENSG00000125851"))
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = ENSG00000125851)) + geom_point(size = 1) + scale_color_viridis_c(option = "magma", name = " ") + theme_minimal()+ labs(title = "PCSK2 expression") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))
#For co-expression plot
df <- FetchData(seurat_object, vars = c("UMAP_1", "UMAP_2", "ENSG00000162670", "ENSG00000125851"))
df$coexpression_score <- (df$ENSG00000162670 + df$ENSG00000125851) / 2
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = coexpression_score)) + geom_point(size = 1) +scale_color_viridis(option = "magma", name = " ") + theme_minimal() + labs(title = "Co-expression of BRINP3 and PCSK2") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))


##Dot plot human hypothalamus
seurat_object <- readRDS("human_hypothalamus.rds")#Human hypothalamus
levels(seurat_object@meta.data$region)
seurat_object@meta.data$region <- as.character(seurat_object@meta.data$region) # Combine "Thalamaus" and "Thalamus" → "Thalamus"
seurat_object@meta.data$region[seurat_object@meta.data$region == "Thalamaus"] <- "Thalamus"
seurat_object@meta.data$region <- factor(seurat_object@meta.data$region) # Convert back to factor
levels(seurat_object@meta.data$region)
DotPlot_scCustom(seurat_object, features = c("ENSG00000162670", "ENSG00000125851"), 
                 group.by = "region", colors_use = rev(viridis_plasma_dark_high))

##DOT PLOT for co-expression in BRINP3+ cells in human VMH
seurat_object <- readRDS("human_hypothalamus.rds") #Human hypothalamus
seurat_object <- NormalizeData(seurat_object, 
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)
hum_vmh_cells <- subset(seurat_object, subset = region ==
                          "VMH" & cell_type == "neuron")
hum_vmh_cells$BRINP3_status <- ifelse(hum_vmh_cells@assays$RNA@data["ENSG00000162670", ] > 0,
                                      "BRINP3+", "BRINP3-")
hum_vmh_cells$only_BRINP3_positive <- ifelse(hum_vmh_cells$BRINP3_status == "BRINP3+" ,"BRINP3+","Other")
BRINP3_positive <- subset(hum_vmh_cells, subset = only_BRINP3_positive == "BRINP3+")
seurat_object <- ScaleData(seurat_object,
                           features = c("ENSG00000162670", "ENSG00000125851", "ENSG00000176697",
                                        "ENSG00000147571", "ENSG00000133636",
                                        "ENSG00000101327", "ENSG00000164326",
                                        "ENSG00000116678", "ENSG00000171105",
                                        "ENSG00000112164", "ENSG00000148053", "ENSG00000091831",
                                        "ENSG00000091664", "ENSG00000136931"))
table(hum_vmh_cells$BRINP3_status)
DotPlot_scCustom( hum_vmh_cells,
                  features = c("ENSG00000162670", "ENSG00000125851", "ENSG00000176697",
                               "ENSG00000147571", "ENSG00000133636",
                               "ENSG00000101327", "ENSG00000164326",
                               "ENSG00000116678", "ENSG00000171105",
                               "ENSG00000112164", "ENSG00000148053", "ENSG00000091831",
                               "ENSG00000091664", "ENSG00000136931"),
                  group.by = "BRINP3_status", dot.scale = 6,
                  colors_use = rev(viridis_plasma_dark_high)) + 
  coord_flip()


# Co-expression of BRINP3 and PCSK2 in human colon, mouse colon and mouse ileum
# Data were collected from https://singlecell.broadinstitute.org/single_cell/study/SCP1038/the-human-and-mouse-enteric-nervous-system-at-single-cell-resolution#study-download

#human colon
#Expression matrix generation
expression_matrix <- readMM("gene_sorted-hli.matrix.mtx")
genes <- read.delim("hli.genes.tsv", header = FALSE)
barcodes <- read.delim("hli.barcodes.tsv", header = FALSE)
rownames(expression_matrix) <- genes$V1
colnames(expression_matrix) <- barcodes$V1
head(expression_matrix[,1:5])
#head(expression_matrix[1:5,])

#tsne file
data.tsne <- read.table("hli.tsne2.txt", header = TRUE, sep = "\t")
head(data.tsne$NAME)  
data.tsne <- data.tsne[-1,]
colnames(data.tsne)
rownames(data.tsne) <- data.tsne$NAME
#head(rownames(data.tsne))
data.tsne <- data.tsne[,-1]
data.tsne <- data.tsne[,-3]
data.tsne[1:5,]
data.tsne_mat <- as(data.tsne, "matrix")
colnames(data.tsne_mat)
colnames(data.tsne_mat) <- paste0("tSNE_", 1:ncol(data.tsne_mat))
colnames(data.tsne_mat)
str(data.tsne_mat)
data.tsne_mat <- apply(data.tsne_mat, 2, as.numeric)
data.tsne_mat <- as.matrix(data.tsne_mat)
#str(data.tsne_mat)

#metadata preparation
metadata <- read.table("all.meta.txt", header = TRUE, sep = "\t")
head(metadata)
head(metadata$NAME)
metadata <- metadata[-1,]
colnames(metadata)
rownames(metadata) <- metadata$NAME
sum(duplicated(metadata$NAME))
metadata$NAME[duplicated(metadata$NAME)]
duplicated_rows <- metadata[duplicated(metadata$NAME) | duplicated(metadata$NAME, fromLast = TRUE), ]
head(duplicated_rows)
View(duplicated_rows)
specific_duplicates <- metadata[metadata$NAME == "H10_S46.AAAGATGAGGTACTCT", ]
head(specific_duplicates)
rownames(metadata) <- make.unique(as.character(metadata$NAME))
metadata[1:5,]
metadata <- metadata[,-1]
metadata[1:5,]

#Creating seurat object
seurat_object <- CreateSeuratObject(counts = expression_matrix,meta.data = metadata)
rownames(data.tsne_mat) <- rownames(seurat_object@meta.data)
length(rownames(seurat_object@meta.data))
length(data.tsne_mat[,1])
rownames(data.tsne_mat) <- rownames(seurat_object@meta.data)
length(rownames(seurat_object@meta.data))
length(data.tsne_mat[,1])
seurat_object[["tsne"]] <- CreateDimReducObject(embeddings = data.tsne_mat,key = "tSNE_",assay = DefaultAssay(seurat_object))
DimPlot(seurat_object, reduction = "tsne")
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- ScaleData(seurat_object, features = c("BRINP3", "PCSK2"))

df <- FetchData(seurat_object, vars = c("tSNE_1", "tSNE_2", "BRINP3", "PCSK2")) #Co-expression plot
df$coexpression_score <- (df$BRINP3 + df$PCSK2) / 2
ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = coexpression_score)) + geom_point(size = 1) + 
  scale_color_viridis(option = "magma", name = " ") + theme_minimal() + labs(title = "Co-expression of BRINP3 and PCSK2") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))

df <- FetchData(seurat_object, vars = c("tSNE_1", "tSNE_2", "BRINP3"))  #BRINP3 expression plot
ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = BRINP3)) + geom_point(size = 1) + scale_color_viridis_c(option = "magma", name = " ") + 
  theme_minimal()+ labs(title = "BRINP3 expression") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))

df <- FetchData(seurat_object, vars = c("tSNE_1", "tSNE_2", "PCSK2"))  #PCSK2 expression plot
ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = PCSK2)) + geom_point(size = 1) + scale_color_viridis_c(option = "magma", name = " ") + 
  theme_minimal()+ labs(title = "PCSK2 expression") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))

DotPlot_scCustom(seurat_object, features = c("BRINP3", "PCSK2"), group.by = "Annotation", colors_use = rev(viridis_plasma_dark_high))

----------------------
  #Mouse Colon
  #Expression matrix generation
  expression_matrix <- readMM("gene_sorted-mli.matrix.mtx")
genes <- read.delim("mli.genes.tsv", header = FALSE)
barcodes <- read.delim("mli.barcodes.tsv", header = FALSE)
rownames(expression_matrix) <- genes$V1
colnames(expression_matrix) <- barcodes$V1
head(expression_matrix[,1:5])

#tsne file
data.tsne <- read.table("mli.tsne2.txt", header = TRUE, sep = "\t")
head(data.tsne$NAME)  
data.tsne <- data.tsne[-1,]
colnames(data.tsne)
rownames(data.tsne) <- data.tsne$NAME
data.tsne <- data.tsne[,-1]
data.tsne <- data.tsne[,-3]
data.tsne[1:5,]
data.tsne_mat <- as(data.tsne, "matrix")
colnames(data.tsne_mat)
colnames(data.tsne_mat) <- paste0("tSNE_", 1:ncol(data.tsne_mat))
colnames(data.tsne_mat)
str(data.tsne_mat)
data.tsne_mat <- apply(data.tsne_mat, 2, as.numeric)
data.tsne_mat <- as.matrix(data.tsne_mat)

#metadata preparation
metadata <- read.table("all.meta.txt", header = TRUE, sep = "\t")
head(metadata)
head(metadata$NAME)
metadata <- metadata[-1,]
colnames(metadata)
rownames(metadata) <- metadata$NAME
sum(duplicated(metadata$NAME))
metadata$NAME[duplicated(metadata$NAME)]
duplicated_rows <- metadata[duplicated(metadata$NAME) | duplicated(metadata$NAME, fromLast = TRUE), ]
head(duplicated_rows)
View(duplicated_rows)
specific_duplicates <- metadata[metadata$NAME == "H10_S46.AAAGATGAGGTACTCT", ]
head(specific_duplicates)
rownames(metadata) <- make.unique(as.character(metadata$NAME))
metadata[1:5,]
metadata <- metadata[,-1]
metadata[1:5,]

#Creating seurat object
seurat_object <- CreateSeuratObject(counts = expression_matrix,meta.data = metadata)
rownames(data.tsne_mat) <- rownames(seurat_object@meta.data)
length(rownames(seurat_object@meta.data))
length(data.tsne_mat[,1])
rownames(data.tsne_mat) <- rownames(seurat_object@meta.data)
length(rownames(seurat_object@meta.data))
length(data.tsne_mat[,1])
seurat_object[["tsne"]] <- CreateDimReducObject(embeddings = data.tsne_mat,key = "tSNE_",assay = DefaultAssay(seurat_object))
DimPlot(seurat_object, reduction = "tsne")
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- ScaleData(seurat_object, features = c("Brinp3", "Pcsk2"))

df <- FetchData(seurat_object, vars = c("tSNE_1", "tSNE_2", "Brinp3", "Pcsk2")) #Co-expression plot
df$coexpression_score <- (df$Brinp3 + df$Pcsk2) / 2
ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = coexpression_score)) + geom_point(size = 1) +scale_color_viridis(option = "magma", name = " ") + 
  theme_minimal() + labs(title = "Co-expression of Brinp3 and Pcsk2") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))

df <- FetchData(seurat_object, vars = c("tSNE_1", "tSNE_2", "Brinp3")) #Brinp3 expression plot
ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = Brinp3)) + geom_point(size = 1) + scale_color_viridis_c(option = "magma", name = " ") + 
  theme_minimal()+ labs(title = "Brinp3 expression") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))

df <- FetchData(seurat_object, vars = c("tSNE_1", "tSNE_2", "Pcsk2")) #Pcsk2 expression plot
ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = Pcsk2)) + geom_point(size = 1) + scale_color_viridis_c(option = "magma", name = " ") + 
  theme_minimal()+ labs(title = "Pcsk2 expression") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))

#mouse ileum
#Expression matrix generation
expression_matrix <- readMM("gene_sorted-msi.matrix.mtx")
genes <- read.delim("msi.genes.tsv", header = FALSE)
barcodes <- read.delim("msi.barcodes.tsv", header = FALSE)
rownames(expression_matrix) <- genes$V1
colnames(expression_matrix) <- barcodes$V1
head(expression_matrix[,1:5])

#tsne file
data.tsne <- read.table("msi.tsne2.txt", header = TRUE, sep = "\t")
head(data.tsne$NAME)  
data.tsne <- data.tsne[-1,]
colnames(data.tsne)
rownames(data.tsne) <- data.tsne$NAME
data.tsne <- data.tsne[,-1]
data.tsne <- data.tsne[,-3]
data.tsne[1:5,]
data.tsne_mat <- as(data.tsne, "matrix")
colnames(data.tsne_mat)
colnames(data.tsne_mat) <- paste0("tSNE_", 1:ncol(data.tsne_mat))
colnames(data.tsne_mat)
str(data.tsne_mat)
data.tsne_mat <- apply(data.tsne_mat, 2, as.numeric)
data.tsne_mat <- as.matrix(data.tsne_mat)

#metadata preparation
metadata <- read.table("all.meta.txt", header = TRUE, sep = "\t")
head(metadata)
head(metadata$NAME)
metadata <- metadata[-1,]
colnames(metadata)
rownames(metadata) <- metadata$NAME
sum(duplicated(metadata$NAME))
metadata$NAME[duplicated(metadata$NAME)]
duplicated_rows <- metadata[duplicated(metadata$NAME) | duplicated(metadata$NAME, fromLast = TRUE), ]
head(duplicated_rows)
View(duplicated_rows)
specific_duplicates <- metadata[metadata$NAME == "H10_S46.AAAGATGAGGTACTCT", ]
head(specific_duplicates)
rownames(metadata) <- make.unique(as.character(metadata$NAME))
metadata[1:5,]
metadata <- metadata[,-1]
metadata[1:5,]

#Creating seurat object
seurat_object <- CreateSeuratObject(counts = expression_matrix,meta.data = metadata)
rownames(data.tsne_mat) <- rownames(seurat_object@meta.data)
length(rownames(seurat_object@meta.data))
length(data.tsne_mat[,1])
rownames(data.tsne_mat) <- rownames(seurat_object@meta.data)
length(rownames(seurat_object@meta.data))
length(data.tsne_mat[,1])
seurat_object[["tsne"]] <- CreateDimReducObject(embeddings = data.tsne_mat,key = "tSNE_",assay = DefaultAssay(seurat_object))
DimPlot(seurat_object, reduction = "tsne")
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- ScaleData(seurat_object, features = c("Brinp3", "Pcsk2"))

df <- FetchData(seurat_object, vars = c("tSNE_1", "tSNE_2", "Brinp3", "Pcsk2")) #Co-expression plot
df$coexpression_score <- (df$Brinp3 + df$Pcsk2) / 2
ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = coexpression_score)) + geom_point(size = 1) + scale_color_viridis(option = "magma", name = " ") + 
  theme_minimal() + labs(title = "Co-expression of Brinp3 and Pcsk2") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))

df <- FetchData(seurat_object, vars = c("tSNE_1", "tSNE_2", "Brinp3")) #Brinp3 expression plot
ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = Brinp3)) + geom_point(size = 1) + scale_color_viridis_c(option = "magma", name = " ") + 
  theme_minimal()+ labs(title = "Brinp3 expression") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))

df <- FetchData(seurat_object, vars = c("tSNE_1", "tSNE_2", "Pcsk2")) #Pcsk2 expression plot
ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = Pcsk2)) + geom_point(size = 1) + scale_color_viridis_c(option = "magma", name = " ") + 
  theme_minimal()+ labs(title = "Pcsk2 expression") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = .5))

DotPlot_scCustom(seurat_object, features = c("Brinp3", "Pcsk2"), group.by = "Annotation", colors_use = rev(viridis_plasma_dark_high))



