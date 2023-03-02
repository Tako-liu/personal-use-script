fpkm <- read.delim("E:/Tomato/tpm/transcript_FPKM_matrix.csv", row.names = 1, sep = ',',check.names = FALSE)
sample_info<-read.csv("D:/mrna-seq/meiling-tomato/sample_info.csv",header = T,row.names = 1)
#rownames(sample_info) <- sample_info[,1]
#sample_info <- sample_info[,-1]
#BiocManager::install("PCAtools")
library(PCAtools)
sample_cor<-cor(fpkm)#计算样本相关性
sample_cor2<-round(sample_cor,digits = 2)#保留两位小数
library(pheatmap)
pheatmap(sample_cor2,display_numbers = TRUE,cluster_rows = FALSE,cluster_cols = TRUE,fontsize = 8)
##计算距离矩阵
sample_dist<-dist(t(fpkm))#先转置在计算距离矩阵


##聚类
sample_hc<-hclust(sample_dist)
plot(sample_hc)
library(PCAtools)
fpkm2 <- scale(fpkm)
set.seed(111)
##主成分分析
pca<-pca(fpkm2,metadata=sample_info,removeVar = 0.1)#样本名字顺序要完全一致！
pca_loadings<-pca$loadings#提取loadings信息
pca_loadings[1:4,1:4]#打印出1-4行/1-4列
pca_rotated<-pca$rotated#提取rotated信息
pca_rotated[1:4,1:4]

screeplot(pca)#画图:主成分对样本差异的解释度
biplot(pca,#画图:PC1/PC2图
       x="PC1",
       y="PC3",#PC1/PC2可换
       shape = "group",)#形状=分组信息
