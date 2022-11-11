library(picante) #READ.NEXUS FUNCTION
#library(dplyr)
#library(reshape2)
#library(betapart)

#LOAD UP THE TREE FROM CHRIS, THIS HAS 100 SPECIES
# newick format
sharkphy2 <- read.tree("/home/lnkdavi/scratch/betaraysFlorent/")
sharkphy<-sample(sharkphy2,size=1)[[1]]

#species matrix
test <- read.csv("/home/lnkdavi/scratch/betaraysFlorent/specieshexmatrix.csv")
#pdpercell
pdpercell <- read.csv("/home/lnkdavi/scratch/betaraysFlorent/pdpercell.csv")
#cellcombo
cellcombo <- read.csv("/home/lnkdavi/scratch/betaraysFlorent/cellcombo.csv")

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
         min.pd= min.pd, simpson.pd=simpson["PD"],
         simpson.sr=simpson["SR"],sorenson.pd=sorenson["PD"], sorenson.sr=sorenson["SR"], holt = holt)
  # beta=c(site1=i,site2=j,cellid1 = as.pdpercell[i,3], cellid2 = pdpercell[j,3], pd.sum=pd.sum, pd.combined=pd.combined, 
  #        min.pd= min.pd, simpson.pd=simpson["PD"],sorenson=sorenson, holt = holt)
  
  return(beta)
}

#here we are creating a function to populate cellcombo (all possible combinations)
populate=function(k,cellcombo,pdpercell, test,sharkphy){
  betafunc(pdpercell, i=cellcombo[k,1], cellcombo[k,2], test,sharkphy)  
}
  
#then define a function that applys pdsumfunc over n number of combination
#is this were we define the number of iterations per job?
dim(cellcombo)
dim(pdpercell)

res=lapply(1:30,populate,cellcombo,pdpercell,test,sharkphy)
class(res)
resfinal=as.data.frame(do.call(rbind, res))
head(resfinal)
class(resfinal)




