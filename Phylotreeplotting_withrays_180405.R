library(ape)
library(Biostrings)
library(ggplot2)
library(ggtree)
#install.packages("geiger")
library(geiger)
library(reshape2)
library(here)

# install.apackages("ggtree")
# source("https://bioconductor.org/biocLite.R")
# biocLite("ggtree")

sharkphy1 <- read.nexus("data/100.Shark.Tree.nex")
class(sharkphy1) 
sharkphy2 <- read.tree("C:/C.drive/Chapter 4/PhylogenyChris/sharknew")
sharkphy_one<-sample(sharkphy2,size=1)[[1]]
str(sharkphy_one)
class(sharkphy_one)
#write.tree(sharkphy_one, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharknew_onetree")   # newick format
#sharkphy <- read.tree("C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharknew_onetree")
ch_trees<- sharkphy2

ch_taxa<- read.csv('C:/C.drive/Chapter 4/PhylogenyChris/SharkTaxonomyForDan.csv')
head(ch_taxa)
family<- unique(ch_taxa$Family)
order<- unique(ch_taxa$Order)
order<- data.frame(unique(ch_taxa$Order))
head(order)
order<- factor(order$unique.ch_taxa.Order.)
class(order)
str(order)
#rays <- c("Chimaeriformes", "Rajiformes", "Myliobatiformes",
#          "Torpediniformes","Rhinopristiformes")
#sharkorder <- order[!order$unique.ch_taxa.Order. %in% rays,]
sharkorder <- family

class(sharkorder)
str(sharkorder)
sharkch_taxa <- ch_taxa

#######################################################
#SHARK ORDER TREE
#######################################################
ord.tree<- list()
for(i in 1:length(ch_trees)){
  m<- match(ch_trees[[i]]$tip.label,ch_taxa$EDGE_Scientific_name)
  ch_trees[[i]]$tip.label<- paste(ch_trees[[i]]$tip.label,ch_taxa$Order[m],sep='.')
  m2<- match(order,gsub(".*\\.","",ch_trees[[i]]$tip.label))
  ord.tree[[i]]<- drop.tip(ch_trees[[i]],setdiff(ch_trees[[i]]$tip.label,ch_trees[[i]]$tip.label[m2]))
  ord.tree[[i]]$tip.label<- gsub(".*\\.","",ord.tree[[i]]$tip.label)
}

ord.tree_one<-sample(ord.tree,size=1)[[1]]
plot.new()
class(ord.tree_one)
laddord.tree_one<-ladderize(ord.tree_one)

#OLD CODE erase if I don't need it
# ordsharktree<- list()
# 
# head(sharkch_taxa)
# rm(i)
# for(i in 1:length(ch_trees)){
#   m<- match(ch_trees[[i]]$tip.label,sharkch_taxa$EDGE_Scientific_name)
#   ch_trees[[i]]$tip.label<- paste(ch_trees[[i]]$tip.label,sharkch_taxa$Family[m],sep='.')
#   m2<- match(sharkorder,gsub(".*\\.","",ch_trees[[i]]$tip.label))
#   ordsharktree[[i]]<- drop.tip(ch_trees[[i]],setdiff(ch_trees[[i]]$tip.label,ch_trees[[i]]$tip.label[m2]))
#   ordsharktree[[i]]$tip.label<- gsub(".*\\.","",ordsharktree[[i]]$tip.label)
# }

# ordshark.tree_one<-sample(ord.tree_one,size=1)[[1]]
# laddord.tree_one_sharks<-ladderize(sample(ord.tree_one, size=1))

#draw in blocks at certain time periods
#plot 5 by 8 landscape
laddord.tree_one_sharks$node.ages
root(laddord.tree_one_sharks)
max(node.depth.edgelength(laddord.tree_one_sharks))

dev.off()
plot.phylo(laddord.tree_one)
edgelabels(round(laddord.tree_one$edge.length,3),cex=0.7)

# create extra margin room on the right for an axis
dev.off()
par(mar=c(1,1,1,1)+2 )
plot.new()
plot.phylo(laddord.tree_one)
plot.phylo(laddord.tree_one, type = "phylogram" ,
           use.edge.length= TRUE, x.lim = 391.914,
           show.node.label=TRUE, show.tip.label=TRUE, plot = FALSE)
axisPhylo(1, las=1, backward = TRUE)
#axisPhylo(1, las=1, backward = FALSE)
#devonian 360-408
rect(391.914-391.914,0, 391.914-360, 60, col= ("#002060"), border = NA) #cabonifoerous #col= ("#002060")
#Carboniferous 286-360
rect(391.914-360,0, 391.914-290, 60, col= ("#00206099"), border = NA) #permian
#Permian 248-286
rect(391.914-290,0, 391.914-251, 60, col= ("#00206088"), border = NA) #triassic
#Triassic 213-248
rect(391.914-251,0, 391.914-201, 60, col= ("#1f4e79"), border = NA)
#Jurassic 144-213
rect(391.914-201,0, 391.914-145, 60, col= ("#2e75b6"), border = NA)
#Cretaceous 65-144
rect(391.914-145,0, 391.914-66, 60, col= ("#9dc3e6"), border = NA)
#Paleogene 24-65
rect(391.914-66,0, 391.914-3, 60, col=("#deebf7"), border = NA)
#Neogene 0-24
rect(391.914-3,0, 391.914-0, 60, col= c("#deebf744"), border = NA)
par(new=TRUE)

plot.phylo(laddord.tree_one, type = "phylogram" ,
           use.edge.length= TRUE, x.lim = 391.914,edge.width =2,
           show.node.label=TRUE, show.tip.label=TRUE, cex.lab = 18)
axisPhylo(1, las=1, backward = TRUE, cex.lab=22)
#export as a 5*8

edgelabels(round(laddord.tree_one$edge.length,3),cex=0.7)

ggtree(laddord.tree_one, layout= "rectangular", 
       branch.length = laddord.tree_one$edge.length)

##############################################
##Get the node age for each family across 100 trees
##############################################
#get the node numbers of the tips across 100 trees
#then assign the edge lengths to the order name
edge.lengths<- list()
tips <- ordsharktree[[100]]$tip.label
for(i in 1:length(ordsharktree)){
  nodes<-sapply(tips,function(x,y) which(y==x),
                y=ordsharktree[[i]]$tip.label)
  ## then get the edge lengths for those nodes
  edge.lengths[[i]]<-setNames(ordsharktree[[i]]$edge.length[sapply(
    nodes,function(x,y) which(y==x),y=ordsharktree[[i]]$edge[,2])],
    names(nodes))
}



library(dplyr)

#then create a database and calculate the mean and std deviation for 
#each family
edgelengthsfinal=as.data.frame(do.call(rbind, edge.lengths))
edgelengthsfinal$nodeest <- "nodest"
#turn from wide to long
edgelengthmelt <- melt(edgelengthsfinal, id.var = "nodeest")
head(edgelengthmelt)
mnedge <- edgelengthmelt %>% group_by(variable) %>% dplyr::summarise(mean = mean(value), stddev = sd(value)) %>%
  ungroup()
head(mnedge)
write.csv(mnedge, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/OFamily_sahrkrays_AgeTable1.csv")

node.height(laddord.tree_one_sharks)
max(node.height(laddord.tree_one_sharks)) #height is the number of edges on the longest path from node to leaf
max(node.depth(laddord.tree_one_sharks)) #depth is the number of edges to root node, 1 is given to tips
node.depth.edgelength(laddord.tree_one_sharks)
max(node.depth(laddord.tree_one_sharks))

#######################################################
#Species TREE
#######################################################

class(ch_trees[[1]])
ch_trees1 <- ch_trees[[1]]
ch_trees1$node.ages
root(ch_trees1)
max(node.depth.edgelength(ch_trees1))

dev.off()
# create extra margin room on the right for an axis
par(mar=c(0,0,0,0) + 0.1)
plot.new()
plot.phylo(ch_trees1, type = "phylogram" ,
           use.edge.length= TRUE, 
           show.node.label=TRUE, show.tip.label=TRUE, plot = FALSE)
axisPhylo(1, las=1, backward = TRUE)
#axisPhylo(1, las=1, backward = FALSE)
rect(365.2529-365.2529,0, 365.2529-360, 1200, col= ("#002060"), border = NA) #cabonifoerous #col= ("#002060")
rect(365.2529-360,0, 365.2529-290, 1200, col= ("#00206099"), border = NA) #permian
rect(365.2529-290,0, 365.2529-251, 1200, col= ("#00206088"), border = NA) #triassic
rect(365.2529-251,0, 365.2529-201, 1200, col= ("#1f4e79"), border = NA)
rect(365.2529-201,0, 365.2529-145, 1200, col= ("#2e75b6"), border = NA)
rect(365.2529-145,0, 365.2529-66, 1200, col= ("#9dc3e6"), border = NA)
rect(365.2529-66,0, 365.2529-3, 1200, col=("#deebf7"), border = NA)
rect(365.2529-3,0, 365.2529-0, 1200, col= c("#deebf744"), border = NA)
par(new=TRUE)

plot.phylo(ch_trees1, type = "phylogram" ,
           use.edge.length= TRUE, 
           show.node.label=TRUE, show.tip.label=TRUE, cex.lab = 18)
axisPhylo(1, las=1, backward = TRUE, cex.lab=22)

edgelabels(round(laddord.tree_one_sharks$edge.length,3),cex=0.7)

ggtree(laddord.tree_one_sharks, layout= "rectangular", 
       branch.length = laddord.tree_one_sharks$edge.length)

##############################################
##Get the node age for each family across 100 trees
##############################################
#get the node numbers of the tips across 100 trees
#then assign the edge lengths to the order name
edge.lengths<- list()
tips <- ordsharktree[[100]]$tip.label
for(i in 1:length(ordsharktree)){
  nodes<-sapply(tips,function(x,y) which(y==x),
                y=ordsharktree[[i]]$tip.label)
  ## then get the edge lengths for those nodes
  edge.lengths[[i]]<-setNames(ordsharktree[[i]]$edge.length[sapply(
    nodes,function(x,y) which(y==x),y=ordsharktree[[i]]$edge[,2])],
    names(nodes))
}



library(dplyr)

#then create a database and calculate the mean and std deviation for 
#each family
edgelengthsfinal=as.data.frame(do.call(rbind, edge.lengths))
edgelengthsfinal$nodeest <- "nodest"
#turn from wide to long
edgelengthmelt <- melt(edgelengthsfinal, id.var = "nodeest")
head(edgelengthmelt)
mnedge <- edgelengthmelt %>% group_by(variable) %>% dplyr::summarise(mean = mean(value), stddev = sd(value)) %>%
  ungroup()
head(mnedge)
write.csv(mnedge, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/OFamily_sahrkrays_AgeTable1.csv")

node.height(laddord.tree_one_sharks)
max(node.height(laddord.tree_one_sharks)) #height is the number of edges on the longest path from node to leaf
max(node.depth(laddord.tree_one_sharks)) #depth is the number of edges to root node, 1 is given to tips
node.depth.edgelength(laddord.tree_one_sharks)
max(node.depth(laddord.tree_one_sharks))

{
#######################################################
#FAMILY TREE
#######################################################
fam.tree<- list()
for(i in 1:length(ch_trees)){
  m<- match(ch_trees[[i]]$tip.label,ch_taxa$EDGE_Scientific_name)
  ch_trees[[i]]$tip.label<- paste(ch_trees[[i]]$tip.label,ch_taxa$Family[m],sep='.')
  m2<- match(family,gsub(".*\\.","",ch_trees[[i]]$tip.label))
  fam.tree[[i]]<- drop.tip(ch_trees[[i]],setdiff(ch_trees[[i]]$tip.label,ch_trees[[i]]$tip.label[m2]))
  fam.tree[[i]]$tip.label<- gsub(".*\\.","",fam.tree[[i]]$tip.label)
}

fam.tree_one<-sample(fam.tree,size=1)[[1]]
plot.new()
class(fam.tree_one)
laddfam.tree_one<-ladderize(fam.tree_one)
plot.phylo(laddfam.tree_one, type = "phylogram" ,
           use.edge.length= TRUE, add.scale.bar(length=500),
           show.node.label=TRUE, show.tip.label=TRUE)
axisPhylo(1, las=1, backward = TRUE)
}

{
#######################################################
#ORDER TREE
#######################################################
ord.tree<- list()
  i=1
for(i in 1:length(ch_trees)){
  m<- match(ch_trees[[i]]$tip.label,ch_taxa$EDGE_Scientific_name)
  ch_trees[[i]]$tip.label<- paste(ch_trees[[i]]$tip.label,ch_taxa$Order[m],sep='.')
  m2<- match(order,gsub(".*\\.","",ch_trees[[i]]$tip.label))
  ord.tree[[i]]<- drop.tip(ch_trees[[i]],setdiff(ch_trees[[i]]$tip.label,ch_trees[[i]]$tip.label[m2]))
  ord.tree[[i]]$tip.label<- gsub(".*\\.","",ord.tree[[i]]$tip.label)
}

ord.tree_one<-sample(ord.tree,size=1)[[1]]
plot.new()
class(ord.tree_one)
laddord.tree_one<-ladderize(ord.tree_one)

plot.phylo(laddord.tree_one, type = "phylogram" ,
           use.edge.length= TRUE, 
           show.node.label=TRUE, show.tip.label=TRUE)
axisPhylo(1, las=1, backward = TRUE)

#draw in blocks at certain time periods
#plot 5 by 8 landscape
laddord.tree_one$node.ages
root(laddord.tree_one)
max(node.depth.edgelength(laddord.tree_one))

plot.new()
plot.phylo(laddord.tree_one, type = "phylogram" ,
           use.edge.length= TRUE, 
           show.node.label=TRUE, show.tip.label=TRUE, plot = FALSE)
axisPhylo(1, las=1, backward = TRUE)
rect(393.6343-400,0, 393.6343-360, 30, col= ("#002060"), border = NA) #cabonifoerous
rect(393.6343-360,0, 393.6343-290, 30, col= ("#00206099"), border = NA) #permian
rect(393.6343-290,0, 393.6343-251, 30, col= ("#00206088"), border = NA) #triassic
rect(393.6343-251,0, 393.6343-201, 20, col= ("#1f4e79"), border = NA)
rect(393.6343-201,0, 393.6343-145, 20, col= ("#2e75b6"), border = NA)
rect(393.6343-145,0, 393.6343-66, 20, col= ("#9dc3e6"), border = NA)
rect(393.6343-66,0, 393.6343-3, 20, col=("#deebf7"), border = NA)
rect(393.6343-3,0, 393.6343-0, 20, col= c("#deebf744"), border = NA)
par(new=TRUE)
plot.phylo(laddord.tree_one, type = "phylogram" ,
           use.edge.length= TRUE, edge.width=2,
           show.node.label=TRUE, show.tip.label=TRUE, bg=NA)
axisPhylo(1, las=1, backward = TRUE)

#then create a database and calculate the mean and std deviation for 
#each family
##############################################
##Get the node age for each family across 100 trees
##############################################
#get the node numbers of the tips across 100 trees
#then assign the edge lengths to the order name
edge.lengths<- list()
tips <- ord.tree[[100]]$tip.label
for(i in 1:length(ord.tree)){
  nodes<-sapply(tips,function(x,y) which(y==x),
                y=ord.tree[[i]]$tip.label)
  ## then get the edge lengths for those nodes
  edge.lengths[[i]]<-setNames(ord.tree[[i]]$edge.length[sapply(
    nodes,function(x,y) which(y==x),y=ord.tree[[i]]$edge[,2])],
    names(nodes))
}
edgelengthsfinal=as.data.frame(do.call(rbind, edge.lengths))
edgelengthsfinal$nodeest <- "nodest"
#turn from wide to long
edgelengthmelt <- melt(edgelengthsfinal, id.var = "nodeest")
head(edgelengthmelt)
mnedge <- edgelengthmelt %>% group_by(variable) %>% dplyr::summarise(mean = mean(value), stddev = sd(value)) %>%
  ungroup()
head(mnedge)
  write.csv(mnedge, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/OrderAgeTable1_sharksandrays.csv")

}

#######################################################
#SHARK ORDER TREE
#######################################################
ordsharktree<- list()

for(i in 1:length(ch_trees)){
  m<- match(ch_trees[[i]]$tip.label,sharkch_taxa$EDGE_Scientific_name)
  ch_trees[[i]]$tip.label<- paste(ch_trees[[i]]$tip.label,sharkch_taxa$Order[m],sep='.')
  m2<- match(sharkorder,gsub(".*\\.","",ch_trees[[i]]$tip.label))
  ordsharktree[[i]]<- drop.tip(ch_trees[[i]],setdiff(ch_trees[[i]]$tip.label,ch_trees[[i]]$tip.label[m2]))
  ordsharktree[[i]]$tip.label<- gsub(".*\\.","",ordsharktree[[i]]$tip.label)
}

ordshark.tree_one<-sample(ordsharktree,size=1)[[1]]
laddord.tree_one_sharks<-ladderize(ordshark.tree_one)
plot.phylo(laddord.tree_one_sharks, type = "phylogram" ,
           use.edge.length= TRUE, 
           show.node.label=TRUE, show.tip.label=TRUE, cex.lab = 18)
axisPhylo(1, las=1, backward = TRUE, cex.lab=22)
edgelabels(round(laddord.tree_one_sharks$edge.length,3),cex=0.7)

ggtree(laddord.tree_one_sharks, layout= "rectangular", 
       branch.length = laddord.tree_one_sharks$edge.length)

##############################################
##Get the node age for each family across 100 trees
##############################################
#get the node numbers of the tips across 100 trees
#then assign the edge lengths to the order name
edge.lengths<- list()
tips <- ordsharktree[[100]]$tip.label
for(i in 1:length(ordsharktree)){
  nodes<-sapply(tips,function(x,y) which(y==x),
                y=ordsharktree[[i]]$tip.label)
  ## then get the edge lengths for those nodes
  edge.lengths[[i]]<-setNames(ordsharktree[[i]]$edge.length[sapply(
    nodes,function(x,y) which(y==x),y=ordsharktree[[i]]$edge[,2])],
    names(nodes))
}

#then create a database and calculate the mean and std deviation for 
#each family
edgelengthsfinal=as.data.frame(do.call(rbind, edge.lengths))
edgelengthsfinal$nodeest <- "nodest"
#turn from wide to long
edgelengthmelt <- melt(edgelengthsfinal, id.var = "nodeest")
head(edgelengthmelt)
mnedge <- edgelengthmelt %>% group_by(variable) %>% dplyr::summarise(mean = mean(value), stddev = sd(value)) %>%
  ungroup()
head(mnedge)
write.csv(mnedge, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/OrderAgeTable1.csv")

node.height(laddord.tree_one_sharks)
max(node.height(laddord.tree_one_sharks)) #height is the number of edges on the longest path from node to leaf
max(node.depth(laddord.tree_one_sharks)) #depth is the number of edges to root node, 1 is given to tips
node.depth.edgelength(laddord.tree_one_sharks)
max(node.depth(laddord.tree_one_sharks))



# ####################################
# #previous code for plotting
# tree <- read.tree("C:/C.drive/Chapter 4/PhylogenyChris/sharknew_onetree")
# tree
# str(tree)
# class(tree)
# max(tree$edge.length)
# ggtree(tree)
# ladtree<-ladderize(tree)
# node.height(tree)
# max(node.height(tree))
# max(node.depth(tree))
# node.depth.edgelength(tree)
# max(node.depth(tree))
# 
# # (7) Designating clades using the node.leaves command in geiger. If you want to name the clade that incorporates Taxon 26 and 27
# #http://www.eve.ucdavis.edu/~wainwrightlab/Roi/Site/Teaching_files/Intro2Phylo_S2.R
# cladeA <- node.leaves(mytree, mrca(mytree)["Taxon_26", "Taxon_27"])
# # Alternatively you can just give the node.leaves command the node number of the basal node of the 
# clade<- node.leaves(mytree, 51) #find the nodes of the families
# # (8) Identifying the distance of a node from the tips of a phylogeny (node height) in an ultrametric tree
# branching.times(tree)
# max(branching.times(tree))
# # (9) Determining if you have polytomies/multichotomies and randomly resolving polytomies/multichotomies with zero length branches
# is.binary.tree(tree) # Returns true or false depending on whether the tree contains polytomies
# 
# plot.phylo(ladtree, type = "phylogram",use.edge.length= TRUE, add.scale.bar(length=500), show.node.label=TRUE, show.tip.label=FALSE)
# axisPhylo(1, las=1, backward = TRUE)
# 
# plot.phylo(tree, type = "phylogram", use.edge.length= TRUE, add.scale.bar(length=500))
# plot.phylo(ladtree, add.scale.bar(length=1000))
# 
# 
# node.depth.edgelength(tree)
# ?node.depth
# get.fields(tree)
# #phylogram, where the x-axis shows the genetic change / evolutionary distance
# #cladogram - is
# #cladogram
# ggtree(tree, branch.length="none") + geom_treescale()
# ggtree(tree, aes(color=group))
# ggtree(tree, mrsd=-400) + theme_tree2() 
# ggtitle("Divergence time")
