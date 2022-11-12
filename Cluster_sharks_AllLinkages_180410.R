#######################################
#Clustering code
#SR only
#Chapter 4
#180407
########################################

#in this code I
#(a) load in database and convert the betadiv ouput to a dist matrix
#(b) run all of the linkages and convert to dendrograms
#(C)

#(a)
#install.packages("dplyr")
#install.packages("yaml")
library(dplyr)
library(data.table)
#install.packages("hydroGOF")
library(dendextend)
#see clValid for the comparison plot between clustering linkages
library(tidyr)
library(reshape2)

##############################
#GET JOB ID FROM THE .SH FILE
##############################
# s <- Sys.getenv("params")
# print(s)
# s=1
#filename <- paste("/home/lnkdavi/sr_180404/cellcombo_sp_jobid/sharks_betadiv_jobid", s, ".csv", sep="")
# filename <- paste("C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/cluster/sharks_betadiv_jobid", s, ".csv", sep="")
# betadiv <- read.csv(file.path(filename))
# betadiv2 <- betadiv[,c("cellid1", "cellid2", "holt.PD")]

setwd("C:/C.drive/Chapter 4/Routputfiles_florent_SR_180403/cluster/")
temp2 = list.files(pattern="*.csv")

for (i in 1:length(temp2)) {
  assign(temp2[i], read.csv(temp2[i]))
}

temp3 <- list()
for (i in 1:length(temp2)) {
  temp3[[i]] <- read.csv(temp2[i])
  type <- strsplit(temp2[i], "\\.") %>% unlist
  #temp3[[i]]<- [,c("cellid1", "cellid2", "holt.PD")] 
}

temp4 <- do.call(what = rbind, args = temp3)
temp5 <- temp4[,c("cellid1", "cellid2", "holt.SR")]
head(temp5)

#just doing this to make this dataset symmetrical
cellid1 <- c(1:100)
cellid2 <- c(1:100)
combo <- expand.grid (cellid1, cellid2)
combo$beta <- rnorm(1:10000, mean = 2, sd =2)
head(combo)
names(combo) <- c("cellid1", "cellid2", "holt.SR")
betadiv<- combo
#turn into a matrix and dist object
matrix<-tidyr::spread(betadiv, cellid2, holt.SR)
rownames(matrix) <- matrix$cellid1
dim(matrix)
matrix2 <- matrix[,-1]
class(matrix2)
matrix3 <- as.matrix(matrix2)
head(matrix3)

matrix3 <- as.dist(matrix2)
class(matrix3)

#(b) create a list of dendrograms from the different linkage methods
###################################################################################################
#code from here:https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html
###################################################################################################
hclust_methods <- c("ward.D", "single", "complete", "average", "mcquitty", 
                    "median", "centroid", "ward.D2")
###############################################################################
#loop over each of the linkage methods and save the output of each dendrogram
###############################################################################
dendlistSR <- dendlist()
clustSRlist <- list()
clustmethods <- function(i,hclust_methods, matrix3, dendlistSR){
  clustSR <- hclust(matrix3, method = hclust_methods[i])
  clustSRlist[i]<<- list(clustSR)
  dendlistSR[i] <<- dendlist(dendlistSR, as.dendrogram(clustSRlist[[i]]))
  }

clustsrdend <- lapply(1:length(hclust_methods), clustmethods,
              hclust_methods,
              matrix3, dendlistSR)
# names(clustsrdend) <- hclust_methods
names(dendlistSR) <- hclust_methods

#(c) calculate the cophenetic correlation coefficient for each of the linkage methods
##############################################################################
#Cophenetic score and correlation between original data and the linkage data
##############################################################################
#this is to create a scatterplot of the Betadiv scores and the linkage scores
#however these plots are not a priority, so I have # out the code to keep the data
#and just kept the global data that will let me calculate the correlation
# coph<- list()
# 
# cophmethods <- function(i,clustSRlist,
#                         hcust_methods){
#   coph[[i]] <<- cophenetic(clustSRlist[[i]])
#   #give a table of cophenetic correlation for each linkage method
#   cordata <- setNames(melt(as.matrix(coph[[i]])), c('rows', 'vars', 'values'))
#   cordata$hclust_method <- as.character(hclust_methods[i])
#   return(cordata)
# }
# 
# cophendatabase <- lapply(1:length(hclust_methods), cophmethods,
#                          clustSRlist,
#                          hcust_methods)
#dismatrixcorr=as.data.frame(do.call(rbind, cophendatabase))
#head(dismatrixcorr)
#fwrite(dismatrixcorr, file = "/home/lnkdavi/sr_180404/copheneticcorr.csv")

#######################################################################
#This is to compare overall correlation scores between linkage methods
#this is Fig S9 Cluster algorithm performance analysis
#######################################################################
coph<- list()

corvaluefunc <- function(i,clustSRlist,matrix3, coph, hcust_methods){
  coph[[i]] <- cophenetic(clustSRlist[[i]]) 
  corvalue<- as.data.frame(cor(matrix3, coph[[i]]))
  corvalue$hclust_method <- as.character(hclust_methods[i])
  return(corvalue)
}

cophcorrvalue <- lapply(1:length(hclust_methods), corvaluefunc, clustSRlist,
                        matrix3, coph, hcust_methods)
class(cophcorrvalue)
cophcorrvaluedb=as.data.frame(do.call(rbind, cophcorrvalue))
head(cophcorrvaluedb)
fwrite(cophcorrvaluedb, file = "/home/lnkdavi/sr_180404/cophenetic.csv")

#(d) Cluster optimizing, height versus number of clusters 
################################################################################
#calculate the height at each number of clusters to optimize number of clusters
#this is the evaluation plot in Kreft and Jetz
#do this just for all linkage methods
#######################################
clust_heights <- list()

optclusfunc <- function(i, dendlistSR, hclust_methods) {
  clust_heights <- as.data.frame(heights_per_k.dendrogram(dendlistSR[[i]]))
  names(clust_heights) <- "dendoheight"
  clust_heights$numclusters <- rownames(clust_heights)
  clust_heights$method <- as.character(hclust_methods[i])
  #plot(dendoheight~ numclusters, UPGMA_heights)
  clust_heights$numclusters <- as.numeric(clust_heights$numclusters)
  clust_heights$ClusterOpt <- round(hydroGOF::rmse(clust_heights$dendoheight, clust_heights$numclusters),0)
  return(clust_heights)
}

optclusouput<- lapply(1:length(hclust_methods), optclusfunc, dendlistSR, 
                      hclust_methods)
outclusall=as.data.frame(do.call(rbind, optclusouput))
head(outclusall) #ClusterOpt
fwrite(cophcorrvaluedb, file = "/home/lnkdavi/sr_180404/cophenetic.csv")

#(e) assign clusters to cell ids
###################################################
#Assign cluster number to each grid cell
###################################################
i=1
head(outclusall)

assignclust <- function(optclusouput) {
  optimimall <- as.numeric(round(hydroGOF::rmse(outclusall$dendoheight, outclusall$numclusters),0))
  optimimall <- as.numeric(round(hydroGOF::rmse(outclusall$dendoheight, outclusall$numclusters),0))
  matrix$clustergroups<-cutree(dendlistSR$average, k=50) # k=upgmaoptimum)
  head(matrix)
}

fwrite(matrix[,"clustergroupsUPGMA"], file = "/home/lnkdavi/sr_180404/copheneticcorr.csv")

#######################################
#do this just for UPGMA
#######################################
dend_average <- dendlistSR$average
dend_average
#heights_per_k.dendrogram(dend_average) #this gives the heights for each cut in the tree
UPGMA_heights <- as.data.frame(heights_per_k.dendrogram(dend_average))
head(UPGMA_heights)
names(UPGMA_heights) <- "dendoheight"
UPGMA_heights$numclusters <- rownames(UPGMA_heights)
plot(dendoheight~ numclusters, UPGMA_heights)

#root mean square error criterion?
#rmse gives the standard deviation of the model prediction error a smaller
#value indicates better model performace


UPGMA_heights$numclusters <- as.numeric(UPGMA_heights$numclusters)
UPGMA_heights$ClusterOpt <- round(hydroGOF::rmse(UPGMA_heights$dendoheight, UPGMA_heights$numclusters),0)
print(UPGMA_heights)
fwrite(UPGMA_heights, file = "/home/lnkdavi/sr_180404/copheneticcorr.csv")



#########################################
#how similar are clusters?
#############################################
#calculate the average beta diversity of other cells outside of each of the clusters
head(outclusall)
average <- filter(outclusall, method == "average")
head(average)




#########################################
#########################################
#########################################
##ROUGH CODE BELOW
#this is a good exercise
library(vegan) # code from here: http://cc.oulu.fi/~jarioksa/opetus/metodi/sessio3.pdf
data(dune) #community matrix
dune
d <- vegdist(dune) #betadiv database 
par(mfrow=c(1,3))
caver <- hclust(d, method="aver")
plot(caver, hang=-1)
caver
summary(caver)
caver$height
caver$order
caver$labels
caver$merge #don't know what this is
plot(caver, hang=-1)
#my_dend <- as.dendrogram(hclust(dist(d)))
my_dend <- as.dendrogram(caver)
plot(my_dend)
#use this to append a column of group onto the comm database
dune$groups<-cutree(caver, k=c(1:10))
dim(dune)
dune[,20:31]
?cutree
heights_per_k.dendrogram(my_dend) #this gives the heights for each cut in the tree
diag <- as.data.frame(heights_per_k.dendrogram(my_dend))
head(diag)
names(diag) <- "dendoheight"
diag$numclusters <- rownames(diag)
plot(dendoheight~ numclusters, diag)
#how do I calculate RMSE?
#root mean square error criterion?
#rmse gives the standard deviation of the model prediction error a smaller
#value indicates better model performace


str(diag)
diag$numclusters <- as.numeric(diag$numclusters)
hydroGOF::rmse(diag$dendoheight, diag$numclusters)



dune %>%
  mutate(cluster = groups) %>%
  head
library(dendextend)
dend_h <- heights_per_k.dendrogram(my_dend)
par(mfrow = c(1,2))
plot(my_dend)
plot(my_dend, ylim = c(dend_h["3"], dend_h["1"]))

a <- fviz_nbclust(dune, FUN = hcut, method = "wss", print.summary = TRUE)

my_dend %>% nleaves
my_dend %>% nnodes #number of nodes includes leaves
my_dend %>% head # A combination of "str" with "head"
sapply(caver, '[')$height #height of each of the nodes
sapply(caver, '[')$order
get_nodes_attr(my_dend, "height")
sort(unique(cophenetic(caver)))

#see this for representation of dendrograms
#https://cran.r-project.org/web/packages/dendextend/vignettes/introduction.html#how-to-explore-a-dendrograms-parameters
par(mfrow = c(1,2))
my_dend %>% set("branches_k_color", k = 3) %>% plot(main = "Nice defaults")
my_dend %>% set("branches_k_color", value = 3:1, k = 3) %>% 
  plot(main = "Controlling branches' colors\n(via clustering)")

#We can collapse branches under a tolerance level using the 
#collapse_branch function:

#sapply(hclust(dist(d)), '[')$height #height of each of the nodes
#sapply(hclust(dist(d)), '[')$order
#library(dendextend)
get_nodes_attr(my_dend, "height")
sort(unique(cophenetic(caver)))


t<-vegemite(dune, caver)
plot(csin, hang = -1)
rect.hclust(csin, 3)
cl <- cutree(csin, 3) #cuts the tree based on the number of clusters you want
cl
table(cl) #summarizes the number of species per cluster
table(cl, cutree(csin, 3))
table(cl, cutree(caver, 3)) #confusion matrix of num spp per different clustering
#cophenetic distance
plot(d, cophenetic(csin), asp=1)
abline(0, 1)
plot(d, cophenetic(ccom), asp=1)
abline(0, 1)
plot(d, cophenetic(caver), asp=1)
abline(0, 1)
#













head(betadiv3)
betadiv3$holt.PD <- rnorm(937, mean = 0 , sd =1)
cellid1 <- c(1,2,3,4,5)
cellid2 <- c(1,2,3,4,5)
dataframe <- expand.grid(cellid1, cellid2)
names(dataframe) <- c("cellid1", "cellid2")
dataframe2<- filter(dataframe, as.numeric(cellid1)> as.numeric(cellid2))
dataframe2$value <- rnorm(10, mean = 0 , sd =1)

matrix<-spread(dataframe2, cellid2, value)
rownames(matrix) <- matrix$cellid1
matrix2 <- matrix[,-1]
head(matrix2)
matrix3 <- as.dist(matrix2)
test <- hclust(matrix3, method = "complete", members = NULL)
summary(test)
test$height
test$order
test$labels
test$merge
plot(test)
?cophenetic
library(picante)
?mpd
mpd(samp, dis, abundance.weighted=FALSE)
distmatrix <- cophenetic(test)
cor(matrix3, distmatrix) #correlation between the original data and the linkage dist matrix

com = hclust(matrix3, method = "complete", members = NULL)
plot(com)

cluster <- function (matrix3)
  
data(UScitiesD)

require(graphics)

### Example 1: Violent crime rates by US state

hc <- hclust(dist(USArrests), "ave")
plot(hc)
plot(hc, hang = -1)

## Do the same with centroid clustering and *squared* Euclidean distance,
## cut the tree into ten clusters and reconstruct the upper part of the
## tree from the cluster centers.
hc <- hclust(dist(USArrests)^2, "cen")
memb <- cutree(hc, k = 10)
class(memb)
plot(memb)
cent <- NULL
for(k in 1:10){
  cent <- rbind(cent, colMeans(USArrests[memb == k, , drop = FALSE]))
}
hc1 <- hclust(dist(cent)^2, method = "cen", members = table(memb))
opar <- par(mfrow = c(1, 2))
plot(hc,  labels = FALSE, hang = -1, main = "Original Tree")
plot(hc1, labels = FALSE, hang = -1, main = "Re-start from 10 clusters")
par(opar)


{
  com = hclust(matrix3, method = "complete", members = NULL)
  upgma = hclust(matrix3, method = "average", members = NULL)
  wardd = hclust(matrix3, method = "ward.D", members = NULL)
  wardd2 = hclust(matrix3, method = "ward.D2", members = NULL)
  wpgma = hclust(matrix3, method = "mcquitty", members = NULL)
  wpgmc = hclust(matrix3, method = "median", members = NULL)
  upgmc = hclust(matrix3, method = "centroid", members = NULL)
  cluster=c(complete = com, upgma=upgma,wardd=wardd,wardd2=wardd2,
            wpgma = wpgma, wpgmc = wpgmc, upgmc = upgmc, 
            cellid1 = cellid1)
  return(cluster)
}


#fwrite(resfinal, file = paste0("C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/test", args, ".csv"))
#CEDAR
fwrite(resfinal, file = paste0("/home/lnkdavi/sr_180404/sr_betadiv_jobid", s, ".csv"))
