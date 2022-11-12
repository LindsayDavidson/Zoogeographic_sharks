library(picante) #READ.NEXUS FUNCTION
library(dplyr)

# #GET JOB ID FROM THE .SH FILE
# args <- commandArgs(trailingOnly = TRUE)
# job.id <- as.numeric(args[1])
# print(args) #just to make sure it was not the step that failed  
# 
#cellcombo
args =2
head(cellcombo)
dim(cellcombo)
allcombo <- cellcombo
dim(allcombo)

cellcombo2 <- read.csv("C:/C.drive/Chapter 4/Routputfiles_florent_180222/rays/cellcombo_rays_jobid/rayscellcombo_1.csv")
head(cellcombo2)
dim(cellcombo2) #3,697,856
cellcombo <- cellcombo2
#cellcombo2 <- read.csv(file.path(paste0("C:/C.drive/Chapter 4/Routputfiles_florent_180222/rays_test/rayscellcombo_",args, ".csv")))
# #cellcombo <- cellcombo2[1:50,]
# head(cellcombo2) #7 unique ell ids.
# dim(cellcombo2)
# 

#LOAD DATABASES
#species matrix
test3 <- read.csv("C:/C.drive/Chapter 4/Routputfiles_florent_180222/rays/rayspecieshexmatrix.csv")
head(test3)
dim(test3) #6082
#test4 <- filter(test3, test3$Unique_ID  %in%  cellcombo$cell1 |
#                  test3$Unique_ID %in%  cellcombo$cell2)
rownames(test3) <- test3$Unique_ID
test3[is.na(test3)] <-0
test5 <- test3[,c(-1,-2)]
head(test5)
dim(test5)
str(test5)
# dim(test4)
# test <- test4[1:10, ]
# dim(test)
# head(test)
# str(test)
# 
# dim(test2) #51
# rownames(test2) <- test2$Unique_ID
# test2[is.na(test2)] <-0
# test <- test2[,c(-1,-2)]
# head(test)

#pdpercell
pdpercell3 <- read.csv("C:/C.drive/Chapter 4/Routputfiles_florent_180222/rays/rays_pdpercell.csv")
head(pdpercell3)
pdpercell2 <- pdpercell3[, -1]
rownames(pdpercell2) <- pdpercell2$X
head(pdpercell2)
dim(pdpercell2)
str(pdpercell2)
pdpercell2$cellid <- as.factor(pdpercell2$cellid)
#pdpercell_test <- filter(pdpercell2, pdpercell2$cellid  %in%  cellcombo$cell1 |
#                  pdpercell2$cellid %in%  cellcombo$cell2)
pdpercell <- pdpercell2
#pdpercell <- pdpercell_test
head(pdpercell)
pdpercell$SR <- as.numeric(pdpercell$SR)
pdpercell$cellid <- as.character(pdpercell$cellid)
str(pdpercell)

#tree
sharkphy <- read.tree("C:/C.drive/Chapter 4/Routputfiles_florent_180222/rays/rays_onetree")

# #Create a vector of every cell combination
# dim(pdpercell)
# cell1 <- pdpercell$cellid
# cell2<- pdpercell$cellid
# cellcombo <- expand.grid (cell1, cell2)
# dim(cellcombo) #36,990,724 #
# colnames(cellcombo) <-c("cell1", "cell2")
# dim(cellcombo)
# str(cellcombo)

# # #I think I need to match all of the datasets, as is, this is the complete pdpercell
# # #complete tree, and a subset of the cellcombo
# # #need to match test and cellcombo first then this??
# head(pdpercell)
# pdpercell <- filter(pdpercell1, pdpercell1$cellid  %in%  cellcombo$cell1 |
#                        pdpercell1$cellid  %in%  cellcombo$cell2)
# 
# dim(test)
# combined <- match.phylo.comm(sharkphy, test)
# head(combined)
# sharkphy <- combined$phy
# test <- combined$comm
#all.equal(rownames(combined$comm), rownames(combined$phy))

# #they are matched, but if they weren't they you can run metadata <- metadata[rownames(comm),]


# #cedar
# test <- read.csv("/home/lnkdavi/scratch/betaraysFlorent/specieshexmatrix.csv")
# #pdpercell
# pdpercell <- read.csv("/home/lnkdavi/scratch/betaraysFlorent/pdpercell.csv")
# #cellcombo
# cellcombo <- read.csv("/home/lnkdavi/scratch/betaraysFlorent/cellcombo.csv")
# sharkphy <- read.tree("/home/lnkdavi/sharknew_onetree")
# i=1
# j=2

betafunc <- function (pdpercell, i, j, test,sharkphy)
{
  pd.sum <- (pdpercell[i,1:2]+pdpercell[j,1:2])#double counts
  pd.combined <- as.matrix(pd((as.matrix(test[i,]+ test[j,])), sharkphy))[1,]#total PD between two cells, this doesn't count shared species twice
  # pd.combined <- as.matrix(pd(t(as.matrix(test[i,]+ test[j,])), sharkphy))[1,]#total PD between two cells, this doesn't count shared species twice
  min.pd=c(PD=min(pd.combined["PD"]-pdpercell[c(i,j),"PD"]),
           SR=min(pd.combined["SR"]-pdpercell[c(i,j),"SR"]))
  pd.shared = (pd.sum-pd.combined)
  simpson = min.pd/(min.pd + pd.shared) #b/(b+a)
  sorenson = (pd.sum/(2*pd.shared +pd.sum))#b+c/(2a+b+c) 
  holt = (1-(pd.shared / (min.pd + pd.shared))) 
  beta=c(site1=i,site2=j,cellid1 = as.numeric(pdpercell[i,3]), cellid2 = as.numeric(pdpercell[j,3]), 
         pd.sum=pd.sum, pd.combined=pd.combined,
         min.pd= min.pd, 
         simpson.pd=simpson["PD"],sorenson.pd=sorenson["PD"], holt = holt)
  # beta=c(site1=i,site2=j,cellid1 = as.pdpercell[i,3], cellid2 = pdpercell[j,3], pd.sum=pd.sum, pd.combined=pd.combined, 
  #        min.pd= min.pd, simpson.pd=simpson["PD"],sorenson=sorenson, holt = holt)
  
  return(beta)
}

#here we are creating a function to populate cellcombo (all possible combinations)
populate=function(k,cellcombo,pdpercell, test,sharkphy){
  betafunc(pdpercell, i=cellcombo[k,1],  j=cellcombo[k,2], test,sharkphy)  
}
  
#then define a function that applys pdsumfunc over n number of combination
#is this were we define the number of iterations per job?

head(cellcombo)
dim(cellcombo)
str(cellcombo)
cellcombo <- cellcombo[,-3]
cellcombo$cell1 <- as.factor(cellcombo$cell1)
cellcombo$cell2 <- as.factor(cellcombo$cell2)

head(pdpercell)
length(unique(pdpercell$cellid))
dim(pdpercell)
str(pdpercell)
sharkphy
dim(test)
head(test)
str(test)

#res=lapply(1:dim(cellcombo_jobid),populate,cellcombo,pdpercell,test,sharkphy)
res=lapply(1:nrow(cellcombo),populate,cellcombo,pdpercell,test,sharkphy)
class(res)
resfinal=as.data.frame(do.call(rbind, res))
head(resfinal)
class(resfinal)




