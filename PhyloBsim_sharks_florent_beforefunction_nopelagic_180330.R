library(picante) #READ.NEXUS FUNCTION
library(dplyr)
library(reshape2)
library(betapart)

###################################################
#LOAD UP THE TREE FROM CHRIS
###################################################
#sharkphy1 <- read.nexus("C:/C.drive/Chapter 4/PhylogenyChris/100.Shark.Tree.nex")
#class(sharkphy1) #has 100 trees to account for the uncertainty
#write.tree(sharkphy1, "C:/C.drive/Chapter 4/PhylogenyChris/sharknew")   # newick format
#this was tree number 1
sharkphy2 <- read.tree("C:/C.drive/Chapter 4/PhylogenyChris/sharknew")
sharkphy_one<-sample(sharkphy2,size=1)[[1]]
write.tree(sharkphy_one, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharknew_onetree")   # newick format
sharkphy <- read.tree("C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharknew_onetree")

# #tree two
# sample(2:100, 1) #sub2 was 3
# #sharkphy2 <- read.tree("C:/C.drive/Chapter 4/PhylogenyChris/sharknew")
# class(sharkphy_sub2)
# sharkphy_sub2<-sharkphy2[[3]] #use the sub2 value here
# write.tree(sharkphy_sub2, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharknew_secondsubset")   # newick format
# sharkphysecond <- read.tree("C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharknew_secondsubset")
# 
# #tree three
# sample(2:100, 1) #this value was 79
# sharkphy_sub3<-sharkphy2[[79]]
# write.tree(sharkphy_sub3, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharknew_thirdsubset")   # newick format
# sharkphythird <- read.tree("C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharknew_thirdsubset")
# 
# #tree four
# sample(2:100, 1) #this value was 34
# sharkphy_sub4<-sharkphy2[[34]]
# sharkphy_sub4$edge.length
# test<- sharkphy_sub4$edge.length - sharkphy_sub3$edge.length #look at differences in edge lengths
# write.tree(sharkphy_sub4, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharknew_fourthsubset")   # newick format
# sharkphyfourth <- read.tree("C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharknew_fourthsubset")

################################################
#CEDAR code
#################################################
# sharkphy1 <- read.nexus("/home/lnkdavi/100.Shark.Tree.nex")
# #class(sharkphy1) #has 100 trees to account for the uncertainty
# write.tree(sharkphy1, "/home/lnkdavi/sharknew")   # newick format
# sharkphy2 <- read.tree("/home/lnkdavi/sharknew")
# #subset the database for one tree
# # sharkphy2
# # class(sharkphy2)
# sharkphy<-sample(sharkphy2,size=1)[[1]]

#######################################################################################################
#NOW THAT I HAVE DEALT WITH THE TREE MOVE ON TO CALCULATIONS,
#THIS CODE IS SIMILAR TO SPECIESBSIM_DATA.R
#######################################################################################################
treehex_unfiltered <- read.csv("C:/C.drive/Chapter 4/GISfiles_ver3/GISexportfiles/SJ_Hex4DegChon_sharks_forR_171005.csv") #sj file of chon and grid cells, also added the Tree names as IUCN and tree names differ because of spelling and taxonomic changes
#treehex_unfiltered <- read.csv("/home/lnkdavi/SJ_HexChond4Deg_withoutIUCNforR_RAYS_171005.csv")
print(length(unique(treehex_unfiltered$TREEname))) #474 species
head(treehex_unfiltered)

treehex <- filter(treehex_unfiltered,TREEname != " " )
#remove the species that don't match Chris's phylogeny list (SppList_170926-CGM.csv)
treehex <- filter(treehex_unfiltered,TREEvalue == 1 )

print(length(unique(treehex$TREEname))) #474 species, no rays or chimaeras
print(length(unique(treehex$Unique_ID))) #number of grid cells 9452 grid cells

###################################
#get rid of pelagic species
##################################
pelagic<- read.csv("C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks_nopelagic/pelagic.shark.list.csv")
head(pelagic)
pel2 <- pelagic[,c("genus_name", "species_na")]
pel2$name <- paste0(pel2$genus_name,"_", pel2$species_na)
head(pel2)

treehex2 <- treehex[!treehex$TREEname %in% pel2$name,]
print(length(unique(treehex2$TREEname))) #457 species, no rays or chimaeras
print(length(unique(treehex2$Unique_ID))) #number of grid cells 3934 grid cells

treehex <- treehex2

widedcast2 <- dcast(treehex,  Unique_ID~TREEname) #convert from long to wide format (ie. pivot table)
widedcast2[1:2,1:2]
dim(widedcast2)
rownames(widedcast2) <- widedcast2$Unique_ID
widedcast2[is.na(widedcast2)] <-0
head(widedcast2)
wide <- widedcast2[,c(-1)]
head(wide)
wide[1:2,1:2]

# ###############################
# #cleaning and matching datasets
# ###############################
combined <- match.phylo.comm(sharkphy, wide)
head(combined)
sharkphy <- combined$phy
widematch <- combined$comm
all.equal(rownames(combined$comm), rownames(combined$phy))
sharkphy #451
dim(widematch) #451
widematch[1:2,1:2]
test <- widematch

write.tree(sharkphy, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks_nopelagic/sharks_onetree")   # newick format
write.csv(widematch, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks_nopelagic/sharkwide_match.csv")   # newick format

###################################
#PD per cell
###################################
pdpercell <- pd(widematch, sharkphy)
head(pdpercell)
pdpercell$cellid <- rownames(pdpercell)
pdpercell$cellid <- as.character(pdpercell$cellid)
dim(pdpercell) #9452
str(pdpercell)
head(pdpercell)
#pdpercell <- filter(pdpercell, SR!=0) #cellid 7150 has no species removed
write.csv(pdpercell, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks_nopelagic/pdpercell_sharks.csv" )
# pdpercell_read <- read.csv("C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/pdpercell_sharks.csv")

################################################
#Create a vector of every cell combination
################################################
pdpercell <- pdpercell_read
cell1 <- pdpercell$cellid
cell2<- pdpercell$cellid
cellcombo <- expand.grid (cell1, cell2)
dim(cellcombo) #89,340,304  combinations of the cells. 
colnames(cellcombo) <-c("cell1", "cell2")
#filter out the top diagonal of the pairwise metric, you don't need this as it is redundant
head(cellcombo)
cellcombo$cell1 <- as.numeric(levels(cellcombo$cell1))[cellcombo$cell1]
cellcombo$cell2 <- as.numeric(levels(cellcombo$cell2))[cellcombo$cell2]
cellcombo2 <- filter(cellcombo, cell1>cell2)
dim(cellcombo2)
head(cellcombo2)
str(cellcombo2)
cellcombo <- cellcombo2

write.csv(cellcombo2, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks_nopelagic/cellcombotest_sharks.csv" )
#use this if yu don't want to run the code to create the database
#cellcombo_read <- read.csv("C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks/cellcombotest_sharks.csv")
# cellcombo <- cellcombo_read[,c("cell1", "cell2")]
# head(cellcombo)
# str(cellcombo)

#create a vector of cellid values
#should be 9452 cells long
#9452/20
#470+470+470+470+470+470+470+470+470+470+470+470+470+470+470+470+470+470+470+470+52

tail(cellcombo) #44,665,425 entries
dim(cellcombo)
cellcombo <- cellcombo[,c("cell1", "cell2")]
uniquecell <- as.data.frame(unique(cellcombo$cell2))
dim(uniquecell) #3933
tail(cellcombo)

# lastvalue <- 23794
# uniquecell2<- rbind(uniquecell, lastvalue)
tail(uniquecell2)
uniquecell$group <- rep(1:20,c(200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,180,153)) 
head(uniquecell)
names(uniquecell) <- c("cellid1", "group")
#match this to cellcombo, then split based on the group. 
head(cellcombo)
head(uniquecell)
cellcombo_grouped<- inner_join(cellcombo, uniquecell, by = c("cell1" = "cellid1"))
tail(cellcombo_grouped)
str(cellcombo_grouped)
jobid <- split(cellcombo_grouped, cellcombo_grouped$group) 
lapply(names(jobid), function(x){write.csv(jobid[[x]], file = paste0("C:/C.drive/Chapter 4/Routputfiles_florent_180222/sharks_nopelagic/cellcombo_sharks_jobid/sharkscellcombo_", x, ".csv"), row.names=FALSE)})

# pdsum = double counting
# Pdcombined = outline of the venn, counts the total, no double counting
# Pdshared = pdsum - pdcombined = shared only
# Pdunique2 = pdcombined - pd1
# Pdunique1 = pdcombined - pd2

str(pdpercell)
str(test)
str(cellcombo)
cellcombo$cell1 <- as.factor(cellcombo$cell1)
cellcombo$cell2 <- as.factor(cellcombo$cell2)


betafunc <- function (pdpercell, i, j, test,sharkphy)
{
  pd.sum <- (pdpercell[i,c("PD", "SR")]+pdpercell[j,c("PD", "SR")])#double counts
  pd.combined <- as.matrix(pd((as.matrix(test[i,]+ test[j,])), sharkphy))[1,]#total PD between two cells, this doesn't count shared species twice
  #pd.combined is JUST the shared species
  # pd.combined <- as.matrix(pd(t(as.matrix(test[i,]+ test[j,])), sharkphy))[1,]#total PD between two cells, this doesn't count shared species twice
  min.pd=c(PD=min(pd.combined["PD"]-pdpercell[c(i,j),"PD"]),
           SR=min(pd.combined["SR"]-pdpercell[c(i,j),"SR"]))
  pd.shared = (pd.sum-pd.combined)
  uniquebc <- pd.combined - pd.shared
  simpson = min.pd/(min.pd + pd.shared) #b/(b+a)
  sorenson = (uniquebc/(2*pd.shared +uniquebc))#b+c/(2a+b+c) 
  # sorenson = (pd.sum/(2*pd.shared +pd.sum))#b+c/(2a+b+c) 
  holt = (1-(pd.shared / (min.pd + pd.shared))) 
  beta=c(site1=i,site2=j,cellid1 = as.numeric(pdpercell[i,"cellid"]), cellid2 = as.numeric(pdpercell[j,"cellid"]), 
         pd.sum=pd.sum, pd.combined=pd.combined,
         min.pd= min.pd, simpson.pd=simpson["PD"],
         simpson.sr=simpson["SR"],sorenson.pd=sorenson["PD"], sorenson.sr=sorenson["SR"], holt = holt)
  # beta=c(site1=i,site2=j,cellid1 = as.pdpercell[i,3], cellid2 = pdpercell[j,3], pd.sum=pd.sum, pd.combined=pd.combined, 
  #        min.pd= min.pd, simpson.pd=simpson["PD"],sorenson=sorenson, holt = holt)
  
  return(beta)
}

#here we are creating a function to populate cellcombo (all possible combinations)
populate=function(k,cellcombo,pdpercell, test,sharkphy){
  betafunc(pdpercell, i=cellcombo[k,"cell1"],  j=cellcombo[k,"cell2"], test,sharkphy)  
}

#cellcombo <- cellcombo[1:10,]

res=lapply(1:nrow(cellcombo),populate,cellcombo,pdpercell,test,sharkphy)
class(res)
resfinal=as.data.frame(do.call(rbind, res))

#check that this gives the same results as the betapart package
#BETA FUNCTIONS FOR THE BETAPART PACKAGE
#see betapart: an R package for the study of beta diversity
#Baselga and Orme MEE

#correlate the two 
sorenson_manual <- as.data.frame(resfinal[,c("cellid1", "cellid2", "sorenson.pd.PD")])
head(sorenson_manual)
?phylo.beta.pair
sorenson <- phylo.beta.pair(test[1:10,], sharkphy, index.family="sorensen")
betasor <- sorenson$phylo.beta.sor
betasor2 <- melt(as.matrix(betasor), varnames = c("cellid1", "cellid2"))
names(betasor2) <- c("cellid1", "cellid2", "sorenson_betapart")
class(betasor2)
sorenson_manual$cellid1<- as.integer(sorenson_manual$cellid1)
sorenson_manual$cellid2<- as.integer(sorenson_manual$cellid2)
sorenson_manual$sorenson.pd.PD<- as.numeric(sorenson_manual$sorenson.pd.PD)
betamatch <- inner_join(betasor2, sorenson_manual, by=c("cellid1", "cellid2"))
head(betamatch)

linearmodel<-lm(sorenson_betapart~sorenson.pd.PD, data=betamatch)
summary(linearmodel)
plot(sorenson_betapart~sorenson.pd.PD, data=betamatch)

