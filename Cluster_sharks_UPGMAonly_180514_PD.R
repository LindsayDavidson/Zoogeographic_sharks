#######################################
#Clustering code
#taxonomic and pd clusters for sharks
#UPGMA and then all linkage methods
#Chapter 4
#180407
########################################

#in this code I
#(a) convert the betadiv ouput to a dist matrix
#(b) 
#(C)

#(a)
library(dplyr)
library(data.table)
library(dendextend)
#see clValid for the comparison plot between clustering linkages
library(tidyr)
library(reshape2)
library(factoextra)


# ########################################################
# #read in beta diversity files and join into one file
# ########################################################
# setwd("/home/lnkdavi/sharks_180328/sharks_betadiv_output/")
# setwd("C:/C.drive/Chapter 4/Routputfiles_florent_SR_180403/cluster/")
# temp2 = list.files(pattern="*.csv")
# for (i in 1:length(temp2)) {
#   assign(temp2[i], read.csv(temp2[i]))
# }
# temp3 <- list()
# for (i in 1:length(temp2)) {
#   temp3[[i]] <- read.csv(temp2[i])[,c("cellid1", "cellid2", "holt.SR")]
# }
# temp4 <- do.call(what = rbind, args = temp3)
# temp3[[i]] <- read.csv(temp2[i])[,c("cellid1", "cellid2", "holt.SR")]
# #just doing this to make this dataset symmetrical
# cellid1 <- c(1:100)
# cellid2 <- c(1:100)
# combo <- expand.grid (cellid1, cellid2)
# combo$beta <- rnorm(1:10000, mean = 2, sd =2)
# head(combo)
# names(combo) <- c("cellid1", "cellid2", "holt.SR")
# betadiv<- combo

#########################################
#Betdiv database that is joined together and cells with less than 5 species filtered out
#create two new database of PD and SR
########################################
betadiv2 <- read.csv("C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/cluster/sharks_beta_PDSR.csv")
betadiv2 <- read.csv("/home/lnkdavi/scratch/sharks_beta_PDSR_greaterthan5.csv")
tail(betadiv2)
max(betadiv2$cellid2)
betadiv <- betadiv2[,c("cellid1", "cellid2", "holt.PD")]
head(betadiv)

###########################################################
#PD clusters
##########################################################
max(betadiv$cellid2)
#betadiv[order]
tail(betadiv)
#tail(matrix$cellid1)
matrix<-tidyr::spread(betadiv, cellid2, holt.PD)
rownames(matrix) <- matrix$cellid1
dim(matrix)
matrix[9440:9450,9440:9450]
matrix2 <- matrix[,-1]
matrixtest <- as.matrix(matrix2)
dim(matrixtest)
matrixtest[1:10, 1:10]
matrixtest[9440:9450, 9440:9450]
matrix3 <- as.dist(matrix2)
matrix3[9440:9450]

###########################################
#PD cluster
#code from here:https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html
# ###############################################
hclust_methods <- c( "average")
clustSR <- hclust(matrix3, method = hclust_methods)
class(clustSR)
dend <- as.dendrogram(clustSR)
nodeheight <- as.data.frame(heights_per_k.dendrogram(dend))
head(nodeheight)
dim(nodeheight)
tail(nodeheight)
head(nodeheight)
names(nodeheight) <- "dendoheight"
nodeheight$numclusters <- rownames(nodeheight)

#############################
#do this at a different time
############################
#calculate the between sum of squares
#install.packages("GMD")
css.obj <- GMD::css.hclust(dist.obj=matrix3,
                           hclust.obj=clustSR,
                           k=50)
head(css.obj)
css.obj$percent <- css.obj$totbss/css.obj$tss*100
elbow.obj <-elbow(css.obj, inc.thres = 0.1, ev.thres=0.89)
print(elbow.obj)
?elbow
fwrite(test, file = "/home/lnkdavi/scratch/wssrays_plotthis.csv")
fwrite(css.obj, file="C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/cluster/clusterwsssharks_plotthis.csv")


css.obj100 <- GMD::css.hclust(dist.obj=matrix3,
                           hclust.obj=clustSR,
                           k=100)
css.obj100$percent <- css.obj100$totbss/css.obj100$tss*100
head(css.obj100)
elbow.obj100 <-elbow(css.obj100, inc.thres = 0.1, ev.thres=0.99)
print(elbow.obj100)

fwrite(test, file = "/home/lnkdavi/scratch/wssrays_plotthis.csv")
fwrite(css.obj100, file="C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/cluster/clusterwsssharks100_plotthis.csv")
#use this one, 100 to define the number of clusters

#elbow.obj <-elbow.batch(css.obj, ev.thres=0.90, inc.thres=0.05)
k <- elbow.obj$k
cutree.obj <- cutree(clustSR, k=c(30,41,42,43,45,46,47,48,k))
head(betadiv2)
head(cutree.obj)
groups2<- as.data.frame(cutree.obj)
head(groups2)
groups2$cellid <- rownames(groups2)
head(groups2)

#write.csv(groups2, "/home/lnkdavi/scratch/betadivsharks_cluster_tomap.csv")
write.csv(groups2, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/cluster/betadivsharks_cluster_tomap.csv")

dev.new(width=12, height=6)
# par(mfcol=c(1,2),mar=c(4,5,3,3),omi=c(0.75,0,0,0))
# plot(betadiv3$cellid1,betadiv3$holt.PD,pch=as.character(css.obj$k),col=css.obj$k,cex=0.75,
#      main="Clusters of simulated data")
plot.elbow(css.obj,elbow.obj,if.plot.new=FALSE)
#plot 5 by 5
ggplot(css.obj, aes(x=k, y =  ev))+ 
  scale_y_continuous( name = "Explained variance")+
  geom_point(colour="black", size=1) + geom_line()+
  scale_x_continuous(name = "Number of clusters\nsharks PD")+
  geom_vline(xintercept=49, colour="red", lwd=0.5) +geom_hline(yintercept=0.89, colour="red", lwd=0.5)+
  # annotate("text", label ="49")  + 
  # annotate("rect", xmin = mean(CL$nclsum), xmax = Inf, ymin = -Inf, ymax = mean(CL$mmtsum) , fill= "grey40") + 
  #scale_colour_manual(values = c("#78B7C5","#3B9AB2", "#E1AF00","#F21A00", "#EBCC2A" ))+
  #scale_alpha_manual(values = c(0.95, 0.95, 0.95, 0.95, 0.95))+
  theme (plot.background = element_rect(fill = "NA", colour = "NA"), 
         axis.line = element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), 
         axis.ticks.length = unit(0,"lines"),
         panel.background = element_rect(fill = "white", colour = "grey50"), 
         axis.text.x= element_text(size = 20, vjust=1,colour = "grey20"),
         axis.text.y= element_text(size = 20, colour= c("grey20")),
         axis.title.x= element_text(size = 20, colour = "grey20"),
         axis.title.y= element_text(size = 20, colour = "grey20"),
         legend.key = element_rect(fill = "grey60", colour = "grey60"))

t +     theme(legend.position="none", 
              panel.border = element_rect(fill="NA",colour = "NA", size=1, linetype="solid"),
              legend.text = element_blank())

      

####################
#run a broken stick model, do this at a different time
######################

# #run a broken stick model on this, get the rnumber of clusters. 
# str(nodeheight)
# nodeheight$numclusters <- as.numeric(nodeheight$numclusters)
# tail(nodeheight)
# dim(nodeheight)
# nodeheight2<- filter(nodeheight, numclusters != 9450)
# nodeheight2<- filter(nodeheight, numclusters <200)
# tail(nodeheight2)
# out.lm <- lm(dendoheight~numclusters, nodeheight2)
# osegm <- segmented::segmented(out.lm, seg.Z=~numclusters)
# 
# qplot(numclusters, dendoheight, group = numclusters > 12, geom = c('point', 'smooth'), 
#       method = 'lm', se= FALSE, lwd=I(1), data = nodeheight2) +
#     xlab("Number of clusters" ) +ylab( "Merging height") +
#   geom_vline(xintercept=c(12), col=I("red"), lwd=I(1)) +
#   theme (plot.background = element_rect(fill = "NA", colour = "NA"), 
#          axis.line.x = element_line(colour = "grey60"),
#          panel.grid.major = element_blank(), 
#          panel.grid.minor = element_blank(), 
#          #axis.ticks.length = unit(0,"lines"),
#          panel.background = element_rect(fill = "white", colour = "grey60"), 
#          axis.text.x= element_text(size = 20, vjust=1,colour = "grey20"),
#          axis.text.y= element_text(size = 20, colour= c("grey20")),
#          # axis.text.y= element_blank(),
#          axis.title.x= element_text(size = 20, colour = "grey20"),
#          axis.title.y= element_text(size = 20, colour = "grey20"), 
#          legend.title=element_blank(), 
#          legend.text = element_text(), 
#          legend.key = element_blank(),
#          legend.position="none", 
#          axis.ticks.length =unit(0.15, "cm"),
#          axis.ticks.x = element_line(colour="grey60"),
#          legend.background = element_blank(), 
#          plot.margin=unit(c(1,1,1,1),"cm"))

#################################################################
#run a second diagnotic that is the MSE per cluster and calcualte the difference in the MSE values, that should be the number of clusters
#calculate the difference inbetween the MSE values
##################################################################
#root mean square error criterion?
#rmse gives the standard deviation of the model prediction error a smaller
#value indicates better model performace

nodeheight$numclusters <- as.numeric(nodeheight$numclusters)
nodeheight$clusterrmse <- round(hydroGOF::rmse(nodeheight$dendoheight, nodeheight$numclusters),0)
head(nodeheight)
fwrite(nodeheight, file = "/home/lnkdavi/scratch/sharkstaxonomicclusterheight_optcluster_toplot.csv")

diffmse <- data.frame()
diffmse <- as.data.frame(as.numeric(nodeheight$numclusters))
names(diffmse) <- "numofclusters"
for (i in 2:nrow(nodeheight)) {
  nodeheight[i,"diff"] <- round(as.numeric(nodeheight[i+1,"dendoheight"] - nodeheight[i,"dendoheight"]),2)*-1
}
head(nodeheight)

plot(diff~ numclusters, nodeheight, xlim = c(0,200)) #elbow plot of dissimilarity and number of clusters
#loess_fit <- loess(diff~ numclusters, nodeheight)
lines(nodeheight$numclusters, predict(loess_fit), col = "blue")
nodeheight2 <- nodeheight[-1,]
test <- nodeheight2[order(nodeheight2$numclusters),]
head(test)

#export 6*8
# qplot(numclusters, diff, geom = c('point', 'smooth'),
qplot(numclusters, diff, geom = c('point'),
     # method = 'loess', se= FALSE,
     size=I(3), data = test) +
  xlab("Number of clusters" ) +ylab( "Height difference\nbetween clusters, pB") +
  scale_x_continuous(limits = c(0,100))+
  geom_vline(xintercept=c(37), col=I("red"), lwd=I(1)) +
  theme (plot.background = element_rect(fill = "NA", colour = "NA"),
         axis.line.x = element_line(colour = "grey60"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         #axis.ticks.length = unit(0,"lines"),
         panel.background = element_rect(fill = "white", colour = "grey60"),
         axis.text.x= element_text(size = 20, vjust=1,colour = "grey20"),
         axis.text.y= element_text(size = 20, colour= c("grey20")),
         # axis.text.y= element_blank(),
         axis.title.x= element_text(size = 20, colour = "grey20"),
         axis.title.y= element_text(size = 20, colour = "grey20"),
         legend.title=element_blank(),
         legend.text = element_text(),
         legend.key = element_blank(),
         legend.position="none",
         axis.ticks.length =unit(0.15, "cm"),
         axis.ticks.x = element_line(colour="grey60"),
         legend.background = element_blank(),
           plot.margin=unit(c(1,1,1,1),"cm"))

#################################################################
#################################################################
#these are different if they are PD or SR, this is PD
#################################################################
#Cut the tree based on the number of clusters identified
groups <- cutree(clustSR, k=c(2,3,4,5,6,7,8,9,10,11,12,13, 18, 29,30,31, 32, 33, 34, 35, 36, 37))
head(groups)
class(groups)
groups2<- as.data.frame(groups)
head(groups2)
groups2$cellid <- rownames(groups2)
head(groups2)
colnames(groups2) <- paste0("X", colnames(groups2))
head(groups2)

######count the number of cells that fall into the 22 cluster category. Get rid of the clusters with less than 10 cells
str(groups2$X37)
#groups3 <- groups2[!duplicated(groups2),]
groups2 <- groups2 %>% group_by(X37) %>% dplyr::mutate(countX37=length(X37)) %>% ungroup()
head(groups2[,15:24])
erase<- groups2 %>% group_by(X37) %>% dplyr::summarise(countX37=length(X37)) #put in results
head(erase)
write.csv(groups2,"C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/sharkPDclusters_tomap.csv" )


#####################################################################
#plot at a certain height that matches the number of clusters,for SOM
#you have to cut then merge back together
#this works - keep
#####################################################################
dend_h <- heights_per_k.dendrogram(dend) # (this can take some time)
length(dend_h)
tail(dend_h)
k=37
dend_h["37"] #what is the height at the 37th cluster
dend_h2 <- as.data.frame(dend_h)
names(dend_h2) <- "dendoheight"
tail(dend_h2)
dend_h2$numclusters <- rownames(dend_h2)

plot(dend, las=1, cex=15) #0.1661958
abline(h=dend_h[k],col="red", lwd=2)

#############################################
#cut the tree then assign a colour, use this for GIS and for making the dendrogram for
#the figure
#Create a function to generate a continuous color palette
###########################################################
#read in the file that has the collapsed branches. 
#and then assign a value based on the dendoheight
# 
# 
# #This adds a column of color values
# # based on the y values
# str(dend_h2)
# dend_PD37 <- dplyr::filter(dend_h2, as.numeric(numclusters) <38)
# dim(dend_PD37)
# dend_PD37$col <- rbPal(37)[as.numeric(cut(dend_PD37$dendoheight,breaks = 37))]
# head(dend_PD37)
# length(unique(dend_PD37$Col))


# dendplot <- dend %>% 
#   set("branches_k_color", k=35) %>% set("branches_lwd", 1.2) %>%
#   set("labels_colors") %>% set("labels_cex", c(.9,1.2)) 
# #%>% 
#  # rotate(as.character(5:1)) %>% #rotate to match labels new order %>%
#   #set("leaves_pch", 19) %>% set("leaves_col", col)
# plot(dendplot, las =1, cex=18) #som figure 2
# 
# dendlad <- ladderize(dendplot)
# 
# plot(dendlad, ylim = c(dend_h["20"], dend_h["1"]), xlim=c(6600,6900))
# plot(dendlad, ylim = c(dend_h["20"], dend_h["1"]), xlim=c(6500,8000))
# 
# dendplot <- dend %>% 
#   set("branches_k_color", k=20) %>% set("branches_lwd", 1.2)
# #dendlad <- ladderize(dendplot)
# plot(dendplot) #plot the full dendrogram
# rect.hclust(clustSR, k=20, border="red")
# #points(3481.843119, 0.4262861, pch=21, bg=2, cex=5)
# #points(8000, 0.4262861, pch=21, bg=2, cex=5)
# str(dendlad)


#####################################################################
#plot at a certain height that matches the number of clusters, main text figure
#you have to cut then merge back together
#this works - keep
#####################################################################
# dend_h <- heights_per_k.dendrogram(dend) # (this can take some time)
# testnode<- get_nodes_attr(dend, "height")

# dend_h["37"] #0.2438443 , 0.1661958 for sharks PD
# dend_h2 <- as.data.frame(dend_h)
# names(dend_h2) <- "dendoheight"
# head(dend_h2)
# dend_h2$numclusters <- rownames(dend_h2)
# 
# plot(dend, las=1, cex=15, type="triangle", edge.root=TRUE) #0.1661958
# abline(h=0.1661958,col="red", lwd=2)

#BASED ON ABOVE NUMBER OF CLUSTERS (56) BUT WITH SOME REMOVED, TRIM THE TREE AND THEN DROP.TIP FOR SOME OF THE CLUSTER LEAVES

#this is the clusterID and cellnumber within each cluster
dend_h["37"]

d <- cut(as.dendrogram(clustSR), h=0.1661958)
plot(d$upper) #this is the plot of the 51 clusters,
dupper$label[dupper$order]
#want to make it look betters and less crowded
#i got rid of the branches that had less than 10 cells in them
num <- c(4,5,6,8,9,10,11,12,14,15,17,18,19,20,21,22,23,26,27,28,32,34,35,37)
branchnum <-paste("Branch", num)
head(branchnum)
dprune <- prune(d$upper,branchnum)
dupper <- as.hclust(dprune)
dupper$label[dupper$order] #labels
dupper.node <- as.matrix(cophenetic(dupper))
dupper.node[1:3,1:3]
dupper.node[upper.tri(dupper.node, diag=TRUE) ] <- NA
plot(hclust(as.dist(dupper.node)), hang = -1)
abline(h= 0.1661958)
duppernode <-hclust(as.dist(dupper.node))
class(duppernode)
duppernode2 <- as.dendrogram(duppernode)
dendprune <- data.frame(get_leaves_attr(duppernode2, "label"))
duppernodes <- data.frame(get_nodes_attr(duppernode2, "height"))
write.csv(duppernodes,"C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/sharks.pd.nodes.csv" )
data.frame(get_leaves_attr(d2, "height"))
get_branches_heights(d2, sort = F)
# testnode<- get_nodes_attr(d2, "height", include_leaves = FALSE, include_branches = TRUE)
# data.frame(get_leaves_attr(d2, "height"))
# attributes(d$upper[[1]])$label
#plot(d2, nodePar = list(lab.cex = 0.6, lab.col = "black", pch = NA))
#abline(h= 0.1718593)

#dendprune <- data.frame(get_leaves_attr(d2, "label")) #export this and use to colour in GIS
#data.frame(get_nodes_attr(d2, "height"))
#data.frame(get_leaves_attr(d2, "height"))
#get_branches_heights(d2, sort = F)

plot(highlight_branches(duppernode2))
abline(h= 0.4320194, col="blue")

firstcol <- "#045A8DFF"
brightgreen <- "#A3FF73FF"
col<- viridisLite::viridis(12,  option = "plasma")
col2 <-viridisLite::viridis(20, option = "viridis")
col3<- c(brightgreen, col2[12], col[2:12] )
col3
#export as letter
duppernode2 %>%set("branches_k_color", value=col3, k = 13) %>%
  set("branches_lwd", 8) %>% plot(main = "Shark PD")

# plot(d2, las=1, cex=15, type="rectangle",edge.root=FALSE) #0.1661958
# dendprune <- data.frame(get_leaves_attr(dprune, "label"))
# dendprune$nodeheight <- data.frame(get_leaves_attr(dprune, "height", include_branches = TRUE))
# dendprune$nodeheight <- data.frame(get_nodes_attr(dprune, "midpoint"))

#what if I tr


trimmedPDnodes<- get_nodes_attr(d$upper , "height")
plot(d$upper, xlim=c(0,105))
abline(h=0.1661958, col="red")
plot(d$upper, xlim=c(700,975)) plot(d$upper, xlim=c(600,875)) 


####use this code to make a plot of the dendrogram
plot(d$upper) #this is the plot of the 51 clusters, 
#want to make it look betters and less crowded
#i got rid of the branches that had less than 10 cells in them
dprune <- prune(d$upper,c("Branch 4", "Branch 5", "Branch 9", "Branch 10", 
                          "Branch 15", "Branch 22", "Branch 26", 
                          "Branch 37"))
plot(dprune)
abline(h=0.1661904, col="red")
class(dprune)
get_nodes_attr(dprune, "label")
plot(dprune)
abline(h=0.07300050, col="red")

dendprune <- data.frame(get_nodes_attr(dprune, "label"))
dendprune$nodeheight <- data.frame(get_nodes_attr(dprune, "height", include_branches = TRUE))
dendprune$nodeheight <- data.frame(get_nodes_attr(dprune, "midpoint"))
#dendprune <- get_nodes_attr(dprune, "height", include_branches = FALSE)

# (this can take some time)
dendprune2 <- as.data.frame(dendprune)
names(dendprune2) <- "dendoheight"
head(dendprune2)
dendprune2$numclusters <- rownames(dendprune2)
plot(dprune, las=1, cex=15) #0.1661958

col=rainbow
rbPal <- colorRampPalette(c('red','blue'))
dendprune2$col <- rainbow(37)[as.numeric(cut(dendprune2$dendoheight,breaks = 37))]
head(dendprune2)
#export this to gis, join based on cluster number and assign colour based onteh col 
#column
write.csv(dendprune2, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/sharks_PDclustersandColours_forgis.csv")


plot(dprune, xlim=c(0,105))
plot(dprune, xlim=c(8000,16075)) 
dprune





dend_h["37"] #0.2438443 , 0.1661958 for sharks PD
dend_h2 <- as.data.frame(dend_h)
names(dend_h2) <- "dendoheight"
head(dend_h2)
dend_h2$numclusters <- rownames(dend_h2)

plot(dprune, xlim=c(0,105))
plot(dprune, xlim=c(700,975)) #this is the plot of the 51 clusters, want to make it look betters and less crowded
plot(dprune, xlim=c(600,875)) 
dprune
# colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))
# col=(colfunc(51))
# 
# dendplot <- dend %>% 
#   set("branches_k_color", k=51) %>% set("branches_lwd", 1.2) %>%
#   set("labels_colors") %>% set("labels_cex", c(.9,1.2)) 
# #%>% 
# # rotate(as.character(5:1)) %>% #rotate to match labels new order %>%
# #set("leaves_pch", 19) %>% set("leaves_col", col)
# plot(dendplot, las =1, cex=18) #som figure 2
# 
# 
# plot(dendlad, ylim = c(dend_h["51"], dend_h["1"]), xlim=c(6600,6900))
# plot(dendlad, ylim = c(dend_h["51"], dend_h["1"]), xlim=c(6500,8000))


##################################################################3
#end of taxonomic analysis
####################################################################



##############################################
#code for calculating the mean phylobeta outside of each cluster

#read in the database that has the cells and number of clusters and phylobeta values?






######################################################
#CODE below for running different linkage methods
#####################################################
dendlistSR <- dendlist()
clustSRlist <- list()
hclust_methods <- c( "average")

###############################################################################
#loop over each of the linkage methods and save the output of each dendrogram
###############################################################################
clustmethods <- function(i,hclust_methods, matrix3, dendlistSR){
  clustSR <- hclust(matrix3, method = hclust_methods[i])
  clustSRlist[i] <<- list(clustSR)
  dendlistSR[i] <<- dendlist(dendlistSR, as.dendrogram(clustSRlist[[i]]))
  }

clustsrdend <- lapply(1:length(hclust_methods), clustmethods,
              hclust_methods,
              matrix3, dendlistSR)
#names(clustsrdend) <- hclust_methods
names(dendlistSR) <- hclust_methods

#######################################################################
#This is to compare overall correlation scores between linkage methods
#this is Fig S9 Cluster algorithm performance analysis
#I don't need this for now as I am just calculating UPGMA only
#######################################################################
# corvaluefunc <- function(i,matrix3, coph, hcust_methods){
#   corvalue<- as.data.frame(cor(matrix3, coph[[i]]))
#   corvalue$hclust_method <- as.character(hclust_methods[i])
#   return(corvalue)
# }
# 
# cophcorrvalue <- lapply(1:length(hclust_methods), corvaluefunc,
#                         matrix3, coph, hcust_methods)
# class(cophcorrvalue)
# cophcorrvaluedb=as.data.frame(do.call(rbind, cophcorrvalue))
# head(cophcorrvaluedb)
# fwrite(cophcorrvaluedb, file = "/home/lnkdavi/sr_180404/cophenetic.csv")

##############################################################################
#Cophenetic score and correlation between original data and the linkage data
##############################################################################
#this is to create a plot of the Betadiv scores and the linkage scores
coph<- list()

cophmethods <- function(i,clustSRlist, matrix3, hcust_methods){
   coph[[i]] <<- cophenetic(clustSRlist[[i]]) 
   #give a table of cophenetic correlation for each linkage method
  cordata <- setNames(melt(as.matrix(coph[[i]])), c('rows', 'vars', 'linkage_values'))
  cordata$hclust_method <- as.character(hclust_methods[i])
  return(cordata)
}

cophendatabase <- lapply(1:length(hclust_methods), cophmethods,
                         clustSRlist, matrix3, hcust_methods)
class(cophendatabase)
dismatrixcorr=as.data.frame(do.call(rbind, cophendatabase))
head(dismatrixcorr)
dim(dismatrixcorr)
originalBetaDiv <- setNames(melt(as.matrix(matrix3)), c('rows', 'vars', 'original_values'))
dismatrixcorr <- inner_join(originalBetaDiv, dismatrixcorr, by=c("rows", "vars" ))
head(dismatrixcorr)
plot(dismatrixcorr$original_values, dismatrixcorr$linkage_values)
summary(lm(dismatrixcorr$original_values~dismatrixcorr$linkage_values))
#need to add the original value to this database so that we can correlate them
fwrite(dismatrixcorr, file = "/home/lnkdavi/sr_180404/copheneticcorr.csv")

################################################################################
#calculate the height at each number of clusters to optimize number of clusters
################################################################################
###########################################################
#create the elbow plot with height and number of clusters
###########################################################
dendlistSR

#######################################
#do this just for UPGMA
#######################################
dend_average <- dendlistSR$average
dend_average
#heights_per_k.dendrogram(dend_average) #this gives the heights for each cut in the tree

#this is a way to calculate wss for different number of clusters
# dend <- as.dendrogram(clustSR)
# wss <- function(d) {
#   sum(scale(d, scale = FALSE)^2)
# }
# 
# wrap <- function(i, hc, x) {
#   cl <- cutree(hc, i)
#   spl <- split(x, cl)
#   wss <- sum(sapply(spl, wss))
#   wss
# }
# ## Takes a little while as we evaluate all implied clustering up to 150 groups
# res <- sapply(seq.int(1, nrow(betadiv)), wrap, hc = clustSR, x = betadiv)
# head(res)
# plot(seq_along(res), res, type = "b", pch = 19)

UPGMA_heights <- data.frame()
length(dend_average)

UPGMA_heights <- as.data.frame(heights_per_k.dendrogram(dend_average))
tail(UPGMA_heights)
dim(UPGMA_heights)
names(UPGMA_heights) <- "dendoheight"
UPGMA_heights$numclusters <- rownames(UPGMA_heights)
plot(dendoheight~ numclusters, UPGMA_heights, xlim = c(0,2000)) #elbow plot of dissimilarity and number of clusters


#ASSIGN CLUSTERS
clustercells <- as.data.frame(rownames(matrix2))
clustercells$clustergroupsUPGMA<-cutree(dendlistSR$average, k=c(20, 138))
tail(clustercells)
clustercells$cellid <- rownames(clustercells)
write.csv(clustercells, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharksclusterandcellsid_toplot.csv")

 ##########################################
# #root mean square error criterion - did not use this
# #rmse gives the standard deviation of the model prediction error a smaller
# #value indicates better model performace
# UPGMA_heights$numclusters <- as.numeric(UPGMA_heights$numclusters)
# tail(UPGMA_heights)
# dim(UPGMA_heights)
# UPGMA_heights$ClusterOpt <- round(hydroGOF::rmse(UPGMA_heights[c(1:1049), "dendoheight"], UPGMA_heights[c(1:1049), "numclusters"],0))
# print(UPGMA_heights)
# fwrite(UPGMA_heights, file = "/home/lnkdavi/sr_180404/clusterheight_optcluster.csv")
# plot(UPGMA_heights$numclusters, UPGMA_heights$dendoheight, xlim=c(0,653))
# #ASSIGN CLUSTERS
# upgmaoptimum <- as.numeric(round(hydroGOF::rmse(UPGMA_heights$dendoheight, UPGMA_heights$numclusters),0))
# head(originalBetaDiv)
# clustercells <- as.data.frame(rownames(matrix2))
# clustercells$clustergroupsUPGMA<-cutree(dendlistSR$average, k=upgmaoptimum)
# tail(clustercells)
# fwrite(clustercells, file = "/home/lnkdavi/sr_180404/clusterandcellsid.csv")

#########################################
#how similar are clusters?
#############################################


#######################################
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
head(outclusall)


###################################################
#Assign cluster number to each grid cell
###################################################
dune$groups<-cutree(caver, k=c(1:10))
dim(dune)




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
