library(picante) #READ.NEXUS FUNCTION
library(dplyr)
#install.packages("dplyr")
library(reshape2)
library(betapart)
library(here)
#beta part package

#LOAD UP THE TREE FROM CHRIS, THIS HAS 100 trees
#sharkphy1 <- read.nexus("C:/C.drive/Chapter 4/PhylogenyChris/100.Shark.Tree.nex")
#class(sharkphy1) #has 100 trees to account for the uncertainty
#write.tree(sharkphy1, "C:/C.drive/Chapter 4/PhylogenyChris/sharknew")   # newick format
sharkphy2 <- read.tree("data/sharknew")
sharkphy_one<-sample(sharkphy2,size=1)[[1]]
write.tree(sharkphy_one, "data/sharknew_onetree")   # newick format
sharkphy <- read.tree("data/sharknew_onetree")

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
treehex_unfiltered <- read.csv("data/SJ_Hex4DegChon_rays_forR_171005.csv") #sj file of chon and grid cells, also added the Tree names as IUCN and tree names differ because of spelling and taxonomic changes
#treehex_unfiltered <- read.csv("/home/lnkdavi/SJ_HexChond4Deg_withoutIUCNforR_RAYS_171005.csv")
print(length(unique(treehex_unfiltered$TREEname)))
head(treehex_unfiltered)

treehex_unfiltered$TREEname

treehex <- dplyr::filter(treehex_unfiltered,TREEname != " " )
#remove the species that don't match Chris's phylogeny list (SppList_170926-CGM.csv)
treehex <- filter(treehex_unfiltered,TREEvalue == 1 )
unique(treehex$TREEvalue)

print(length(unique(treehex$TREEname))) #542 species??
print(length(unique(treehex$Unique_ID))) #number of grid cells 6082 grid cells

####get rid of pelagic species
head(treehex)
treehex2<- filter(treehex, TREEname != "Manta_birostris" )
print(length(unique(treehex2$TREEname)))
treehex3<- filter(treehex2, TREEname != "Manta_alfredi"  )
treehex <- treehex3
treehex3<- filter(treehex, TREEname != "Aetobatus_narinari"  )
treehex <- treehex3
print(length(unique(treehex$TREEname))) #539 species??
print(length(unique(treehex$Unique_ID))) #number of grid cells 3430 grid cells

widedcast <- dcast(treehex,  Unique_ID~TREEname) #convert from long to wide format (ie. pivot table)
head(widedcast)
rownames(widedcast) <- widedcast$Unique_ID
widedcast[is.na(widedcast)] <-0
widedcast[1:2,1:2]
wide <- widedcast[,c(-1)]
#apply(wide,1,sum)

# ###############################
# #cleaning and matching datasets
# ###############################
combined <- match.phylo.comm(sharkphy, wide)
#head(combined)
sharkphy <- combined$phy
wide <-combined$comm
head(wide)
length(wide)
write.tree(sharkphy, "data/rays_onetree")   # newick format
write.csv(wide, "data/rayspecieshexmatrix.csv" )

#widematch <- combined$comm
#all.equal(rownames(combined$comm), rownames(combined$phy))
#test <- widematch[c(1:10),] #RUN THE code ON THIS TEST DF
#test<- widematch
#head(test)
#unique(rownames(test)) #50 rownames
test <- wide

# #make a database of cellcombos
# df <- as.data.frame()
# cellid <- as.data.frame(rownames(test))
# #Create a vector of every cell combination
# cell1 <- cellid
# rownames(cell1)<-"cell1"
# cell2<- cellid
# cellcombo <- expand.grid (cell1, cell1)
# dim(cellcombo) #100 combinations of the cells. 
# colnames(cellcombo) <-c("cell1", "cell2")
# dim(cellcombo)
# #cellcombotest <- cellcombo[c(1:10),]
# write.csv(cellcomboall, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/cellcombotest_ray.csv" )
# #write.csv(cellcombo, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/cellcombo.csv" )

#PD per cell
pdpercell <- pd(test, sharkphy)
pdpercell$cellid <- rownames(pdpercell)
dim(pdpercell) #10
str(pdpercell)
head(pdpercell)
write.csv(pdpercell,"data/rays_pdpercell.csv")

#Create a vector of every cell combination
cell1 <- pdpercell$cellid
cell2<- pdpercell$cellid
cellcombo <- expand.grid (cell1, cell2)
dim(cellcombo) #100 combinations of the cells.
colnames(cellcombo) <-c("cell1", "cell2")
dim(cellcombo)
#cellcombotest <- cellcombo[c(1:10),]
write.csv(cellcombo, "data/cellcombotest_ray.csv" )
#write.csv(cellcombo, "C:/C.drive/Chapter 4/Routputfiles_florent_180222/cellcombo.csv" )


#I believe this is all for creating a bunch of files to iterate over. Not sure I need this over. I used this when I sent jobs to the server. 
#create a vector of cellid values
#should be 6,082 cells long
#608+608+608+608+608+608+608+608+608+610
#create 10 different vectors of the different cell ids, divided up by the numbers above
#ie so the first 608, then the next 608 etc.
head(cellcombo) #36990724 entries
cellcombo$cell1 <- as.numeric(levels(cellcombo$cell1))[cellcombo$cell1]
cellcombo$cell2 <- as.numeric(levels(cellcombo$cell2))[cellcombo$cell2]
#filter out the top diagonal of the pairwise metric, you don't need this as it is redundant
cellcombo <- filter(cellcombo, cell1>cell2)

uniquecell2 <- as.data.frame(unique(cellcombo$cell1))
dim(uniquecell2) #3429
head(uniquecell2)
head(cellcombo)
tail(uniquecell2)
addone<-3778
uniquecell <-rbind(uniquecell2, addone)
head(uniquecell)
tail(uniquecell)

308*10+310 #waht is this???
#trying to get to 3199??
#original
#uniquecell$group <- 
#  rep(1:12,c(308,308,308,308,308,308,308,308,308,308,308,42))  #why 12 groups? what is this?
uniquecell$group <- 
  rep(1:11,c(308,308,308,308,308,308,308,308,308,308,119))  #why 12 groups? what is this?
dim(uniquecell)
head(uniquecell)
names(uniquecell) <- c("cellid1", "group")

#match this to cellcombo, then split based on the group. 
head(cellcombo)
head(uniquecell)
cellcombo_grouped<- inner_join(cellcombo, uniquecell, by = c("cell1" = "cellid1"))
head(cellcombo_grouped)
jobid <- split(cellcombo_grouped, cellcombo_grouped$group) 

lapply(names(jobid), function(x){write.csv(jobid[[x]], file = paste0("output/rayscellcombo_", x, ".csv"), row.names=FALSE)})

#and pdcombined, which is the shared
#i added the cell id values here

# pdsum = double counting
# Pdcombined = outline of the venn, counts the total, no double counting
# Pdshared = pdsum - pdcombined = shared only
# Pdunique2 = pdcombined - pd1
# Pdunique1 = pdcombined - pd2



betafunc <- function (i, j, test,sharkphy)
{
  pd1 <- pd(as.matrix(test[rownames(test) ==i,]),sharkphy)
  pd2 <- pd(as.matrix(test[rownames(test) ==j,]),sharkphy)
  pd.sum<- pd1+pd2
  #pd.sum <- (pdpercell[i,c("PD", "SR")]+pdpercell[j,c("PD", "SR")])#double counts
  pd.combined <- as.matrix(pd((as.matrix(test[rownames(test) ==i,]+ test[rownames(test) ==j,])), sharkphy))[1,]#total PD between two cells, this doesn't count shared species twice
  #pd.combined is JUST the shared species
  # pd.combined <- as.matrix(pd(t(as.matrix(test[i,]+ test[j,])), sharkphy))[1,]#total PD between two cells, this doesn't count shared species twice
  min.pd=c(PD=min(pd.combined["PD"]-c(pd1[,"PD"],pd2[,"PD"])),
           SR=min(pd.combined["SR"]-c(pd1[,"SR"],pd2[,"SR"])))
  # min.pd=c(PD=min(pd.combined["PD"]-pdpercell[c(i,j),"PD"]),
  #          SR=min(pd.combined["SR"]-pdpercell[c(i,j),"SR"]))
  pd.shared = (pd.sum-pd.combined)
  uniquebc <- pd.combined - pd.shared
  simpson = min.pd/(min.pd + pd.shared) #b/(b+a)
  sorenson = (uniquebc/(2*pd.shared +uniquebc))#b+c/(2a+b+c) 
 # sorenson = (pd.sum/(2*pd.shared +pd.sum))#b+c/(2a+b+c) 
  holt = (1-(pd.shared / (min.pd + pd.shared))) 
  beta=c(pd1 = pd1, pd2=pd2,site1=i,site2=j,cellid1 =as.numeric(rownames(test)[row.names(test) %in% i]), cellid2 = as.numeric(rownames(test)[row.names(test) %in% j]), 
         pd.sum=pd.sum, pd.combined=pd.combined,
         min.pd= min.pd, simpson.pd=simpson["PD"],
         simpson.sr=simpson["SR"],sorenson.pd=sorenson["PD"], sorenson.sr=sorenson["SR"], holt = holt)
  # beta=c(site1=i,site2=j,cellid1 = as.pdpercell[i,3], cellid2 = pdpercell[j,3], pd.sum=pd.sum, pd.combined=pd.combined, 
  #        min.pd= min.pd, simpson.pd=simpson["PD"],sorenson=sorenson, holt = holt)
  
  return(beta)
}


#here we are creating a function to populate cellcombo (all possible combinations)
populate=function(k,cellcombo, test,sharkphy){
  betafunc(i=cellcombo[k,"cell1"], j=cellcombo[k,"cell2"], test,sharkphy)  
}

#then define a function that applys pdsumfunc over n number of combination
dim(cellcombo)
head(cellcombo)
#is this were we define the number of iterations per job?
res=lapply(1:nrow(cellcombo),populate,cellcombo,test,sharkphy)
res=lapply(1:100,populate,cellcombo,test,sharkphy)
class(res)
resfinal=as.data.frame(do.call(rbind, res))

#check that this gives the same results as the betapart package
#BETA FUNCTIONS FOR THE BETAPART PACKAGE
#see betapart: an R package for the study of beta diversity
#Baselga and Orme MEE

#correlate the two 
sorenson_manual <- as.data.frame(resfinal[,c("cellid1", "cellid2", "sorenson.pd.PD")])
head(sorenson_manual)
sorenson <- phylo.beta.pair(test, sharkphy, index.family="sorensen")
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

