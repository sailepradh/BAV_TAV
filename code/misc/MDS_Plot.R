library("edgeR")

setwd ("/Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/PCA_heatMaps")

SP_BAV_PD = read.table( "BAV_all_SP.txt",sep ="\t", header= FALSE)
SP_TAV_PD = read.table( "TAV_all_SP.txt",sep ="\t", header= FALSE)

probes_BAV = as.character((SP_BAV_PD$V1))
probes_TAV = as.character((SP_TAV_PD$V1))

Common_genes = intersect(probes_BAV,probes_TAV)

SP_BAV_BAVTAV <-SP_BAV_PD[which(SP_BAV_PD$V1 %in% Common_genes),]
SP_TAV_BAVTAV <-SP_TAV_PD[which(SP_TAV_PD$V1 %in% Common_genes),]


counts_BAV <- SP_BAV_BAVTAV[,c(12,13,14)]
row.names(counts_BAV) <- as.character(SP_BAV_BAVTAV$V1)

counts_TAV <- SP_TAV_BAVTAV[,c(12,13,14)]
row.names(counts_TAV) <- as.character(SP_TAV_BAVTAV$V1)

counts_SP_BAVTAV <- cbind (counts_TAV,counts_BAV)
colnames (counts_SP_BAVTAV) <-  c("TAV2431","TAV2515","TAV2709",
                                  "BAV2375", "BAV2424","BAV2714")
row.names (counts_SP_BAVTAV)

mobDataGroups <- c("TAV", "TAV", "TAV",
                   "BAV","BAV","BAV")

d <- DGEList(counts=counts_SP_BAVTAV,
             group=factor(mobDataGroups))

dim(d)
d.full <- d
head(d$counts)
head(cpm(d))
d <- calcNormFactors(d)
d
plotMDS(d, col=as.numeric(d$samples$group))

############################
### PCA plot of the above MDS plot
library(RColorBrewer)
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/factoextra")
library(factoextra)

row_sub = apply(counts_SP_BAVTAV[1:6], 1, function(row) all(row !=0 ))
counts_SP_BAVTAV_refined  <- counts_SP_BAVTAV[row_sub,]
dim (counts_SP_BAVTAV)
dim (counts_SP_BAVTAV_refined)

log_TF <- log2(counts_SP_BAVTAV_refined [,1:6])
head (log_TF)
tr <- t(log_TF)
ir.pca <- prcomp(tr )
ir.genes <- row.names(tr)
cal= c("TAV", "TAV", "TAV",
       "BAV","BAV","BAV")

fviz_eig(ir.pca)
print(ir.pca)
plot(ir.pca, type = "l")
summary(ir.pca)
plot(ir.pca$x[,1:2],xlab = "PC1", ylab = "PC2", col =as.factor(cal), pch =19)
lab= c("TAV", "TAV", "TAV",
       "BAV","BAV","BAV")
p = fviz_pca_ind(ir.pca, pointsize = 10,
                 labelsize =10, habillage=lab,
                 mean.point = FALSE, legend.title = "cell type")+
  labs(title ="Biological replicates BAV and TAV") + 
  xlim(-200,200) + ylim (-200, 200)

cols <-  c("#009E73","#0072B2")

pdf("PCA_plot_replicates.pdf", width = 15, height = 15)
p+ scale_color_manual(values = cols) +
  theme(legend.text=element_text(size=30),
        legend.title = element_text(size=30),
        legend.position = 'bottom',
        legend.key.size = unit(5, 'lines'),
        axis.title.x = element_text(size=30,face="bold", margin = margin(t = 20, r = 0, b =0, l = 0)),
        axis.text.x = element_text(size=20),
        axis.title.y = element_text(size=30,face="bold",margin = margin(t = 0, r = 20, b =0, l = 0)),
        axis.text.y = element_text(size=20),
        plot.title = element_text(colour="grey20",size=30,hjust=0.5),
        plot.margin = unit(c(2, 2, 2, 2), "cm"))
dev.off()
#ggtitle("PD Interactome profile of Unexpressed genes \n in all three cell types and replicates \n")
#dev.off()

##########################################################################################
## Dendogram and Pearson coefficient calculation
library(gplots)
data=as.matrix(tr)
head (tr)

dim (tr)
dd <- dist (scale (tr), method = "euclidean")
hc <- hclust (dd, method = "ward.D2")
plot(hc)
plot(hc, hang = -1, cex = 0.6)
hcd <- as.dendrogram(hc)
plot(hcd, type = "rectangle", ylab = "Height")

nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
plot(hcd,  xlab = "Height",
     nodePar = nodePar, horiz = TRUE)


### Visulization using default theme and theme_dendro()
#install.packages("ggdendro")

library(ggplot2)
library(ggdendro)
ggdendrogram(hc)
ggdendrogram(hc, rotate = TRUE, theme_dendro = TRUE)
ggplot(hc)

###############################################
dd <- dist (scale (tr), method = "euclidean")
clus <- hcluster(dd, method = 'corr')

##############################################################################################

library(mvtnorm)
#nstall.packages("psych")
library(psych)
pairs.panels(log_TF[,1:6], 
             method = "pearson",
             cor.coef = TRUE
             # correlation method
             #hist.col = "#00AFBB",
             #density = TRUE,  # show density plots
             #ellipses = TRUE # show correlation ellipses
)


#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
library(ggplot2)
library(ggpubr)
library(gridExtra)

p <- list()
for (i in c(1:6)){
  for (j in c(1:6)) {
    names <- c ("TAV2431","TAV2515",  "TAV2709","BAV2375","BAV2424","BAV2714")
    x_name <- names[i]
    y_name <- names[j]
    print (x_name)
    print (y_name)
    
    p1 <- ggscatter(log_TF, x = x_name, y = y_name,
                    add = "reg.line" , conf.int =  TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    xlab = x_name , ylab = y_name)
    print (p1)
  }
}

grid.arrange(grobs = p[1])

multiplot(plotlist = p,cols = 6)
grid.arrange(p,ncol=6, nrow = 6)


install.packages("ggcorrplot")
library(ggcorrplot)
ggcorrplot(log_TF)



corr <-cor(log_TF[,1:6])
gcorrplot(corr)