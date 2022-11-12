#before you run this on cedar you need to load in R and load the packages each time
#to do this type: module load r, then R

library(picante) #READ.NEXUS FUNCTION
library(dplyr)
#library(betapart) #betapart will compute species beta diversity
#phylo-beta is also in betapart but Flo has written up code that is faster 
#hclust - use this to bluster beta-div matrix
#several clusting methods should be used and test (which clusterization procedure give you the best tree
library(reshape2)
#install.packages("plyr")
#install.packages("geometry")
#library(plyr)
#install.packages("betapart")
library(betapart)

# #LOAD UP THE TREE FROM CHRIS, THIS HAS 100 SPECIES
# sharkphy1 <- read.nexus("C:/C.drive/Chapter 4/PhylogenyChris/100.Shark.Tree.nex")
# #class(sharkphy1) #has 100 trees to account for the uncertainty
# write.tree(sharkphy1, "C:/C.drive/Chapter 4/PhylogenyChris/sharknew")   # newick format
# sharkphy2 <- read.tree("C:/C.drive/Chapter 4/PhylogenyChris/sharknew")
# #subset the database for one tree
# # sharkphy2
# # class(sharkphy2)
# sharkphy<-sample(sharkphy2,size=1)[[1]]

sharkphy1 <- read.nexus("/home/lnkdavi/scratch/100.Shark.Tree.nex")
#class(sharkphy1) #has 100 trees to account for the uncertainty
write.tree(sharkphy1, "/home/lnkdavi/scratch/sharknew")   # newick format
sharkphy2 <- read.tree("/home/lnkdavi/scratch/sharknew")
#subset the database for one tree
# sharkphy2
# class(sharkphy2)
sharkphy<-sample(sharkphy2,size=1)[[1]]


#######################################################################################################
#NOW THAT I HAVE DEALT WITH THE TREE MOVE ON TO CALCULATIONS,THIS CODE IS SIMILAR TO SPECIESBSIM_DATA.R
#######################################################################################################
#treehex_unfiltered <- read.csv("C:/C.drive/Chapter 4/GISfiles_ver2/SJ_HexChond4Deg_withoutIUCNforR_171005.csv") #sj file of chon and grid cells, also added the Tree names as IUCN and tree names differ because of spelling and taxonomic changes
treehex_unfiltered <- read.csv("/home/lnkdavi/scratch/betasharksandrays/SJ_HexChond4Deg_withoutIUCNforR_171005.csv")

treehex <- filter(treehex_unfiltered,TREEname != " ")
head(treehex)
print(length(unique(treehex$TREEname)))
print(length(unique(treehex$Unique_ID)))

widedcast <- dcast(treehex,  Unique_ID~TREEname) #convert from long to wide format (ie. pivot table)
rownames(widedcast) <- widedcast$Unique_ID
widedcast[is.na(widedcast)] <-0
wide <- widedcast[,c(-1,-2)]
print(str(wide))
#print(class(wide))

# ###############################
# #cleaning and matching datasets
# ###############################
combined <- match.phylo.comm(sharkphy, wide)
#head(combined)
sharkphy <- combined$phy
wide <- combined$comm

all.equal(rownames(combined$comm), rownames(combined$phy))
# #they are matched, but if they weren't they you can run metadata <- metadata[rownames(comm),]

#getAnywhere(phylo.betapart.core) get the code behind this function

dim(wide) #5385
#wide <- wide[c(1:10),]
test <- wide

# #BETA FUNCTIONS FOR THE BETAPART PACKAGE 
# #see betapart: an R package for the study of beta diversity
# #Baselga and Orme MEE
# min <- (phylo.betapart.core(wide, sharkphy))$min.not.shared
# head(min)
# #write.table(min, "/home/lnkdavi/min.csv", sep = ",", 
# #            row.names = TRUE, col.names = TRUE)
# shared <- (phylo.betapart.core(wide, sharkphy))$shared
# head(shared)
# #write.table(min, "/home/lnkdavi/shared.csv", sep = ",", 
# #            row.names = TRUE, col.names = TRUE )
# bsimPDbetapart <- 1-(shared / (min + shared))
# write.table(bsimPD2, "/home/lnkdavi/bsim_rays.csv", sep = ",", 
#             row.names = TRUE, col.names = TRUE)
# as.matrix(bsimPDbetapart)
# write.table(as.matrix(bsimPDbetapart), "C:/C.drive/Chapter 4/RoutputFiles/bsim_raystest.csv", sep = ",", 
#             row.names = TRUE, col.names = TRUE)

#################################################################
##USING THE PICANTE PACKAGE AND LOOPS TO GET THE SAME RESULTS
#################################################################
sharedloop = matrix(nrow = nrow(test),ncol =nrow(test),0)
total = matrix(nrow = nrow(test),ncol =nrow(test),0)

for (i in (1:nrow(test))){ #
  for (j in (1:nrow(test))) {
    if (i>j) {
      pd.1 <- pd(as.matrix(test[i,]), sharkphy) [,1] #drop the SR column, ignore this for now, as if I keep it I can keep the hexagon name in the rownames
      pd.2 <- pd(as.matrix(test[j,]), sharkphy)[,1] #drop the SR column, so these columns have the same species in them,
      pd.total <- pd(as.matrix(test[i,])+(as.matrix(test[j,])), sharkphy)[,1] #toal PD between two cells, this doesn't count shared species twice
      pd.shared <- (pd.1+pd.2) - pd.total
      sharedloop[i,j] <-pd.shared
      total [i,j] <- pd.total
      #write.table(total,"C:/C.drive/Chapter 4/rayphylo_total.csv", sep = ",", row.names = TRUE)
      #write.table(sharedloop,"C:/C.drive/Chapter 4/rayphylo_shared.csv", sep = ",", row.names = TRUE)}
      #write.table(vec,"/home/lnkdavi/rayphylo_shared.csv", sep = ",", row.names = TRUE)}
    }
  }
}
class(total)
rownames(sharedloop) <- rownames(test)
colnames(sharedloop) <- rownames(test)
rownames(total) <- rownames(test)
colnames(total) <- rownames(test)
pdcombined <-as.matrix(sharedloop)
write.table(sharedloop,"/home/lnkdavi/scratch/betasharksandrays/sharksandraysphylo_sharedPD.csv", sep = ",", row.names = TRUE)
#write.table(total,"C:/C.drive/Chapter 4/RoutputFiles/rayphylo_totalPD.csv", sep = ",", row.names = TRUE)
#write.table(sharedloop,"C:/C.drive/Chapter 4/RoutputFiles/rayphylo_sharedPD.csv", sep = ",", row.names = TRUE)
write.table(total,"/home/lnkdavi/scratch/betasharksandrays/sharkphylo_totalPD.csv", sep = ",", row.names = TRUE)

# #tail(pdcombined)

#################################################################################
#CALCUALTE B - UNIQUE PD IN FIRST CELL
##################################################################################
#calculate PD per cell, make a list
PDpercell <- data.frame()
for (i in (1:nrow(test))){
  #pdpercellrun <- pd(comm[i,], phy) #[,1] #drop the SR column, ignore this for now, as if I keep it I can keep the hexagon name in the rownames
  pdpercellrun <- pd(test[i,], sharkphy) #[,1] #if the loop doesn't run try the above line instead
  PDpercell <- rbind(pdpercellrun, PDpercell)
  #pd[i,j] <- pd.1
} 
PDpercell$Unique_ID <- row.names(PDpercell)
head(PDpercell)

#merge this database with the total PD between two cells that does not double count shared species
#first melt the matrix
meltedtotal <- melt(total)
names(meltedtotal) <- c("cell.row", "cell.col", "totalshared")
PDpercell$Unique_ID<- as.numeric(PDpercell$Unique_ID)

#join the database, total.shared pd with the pd per cell
meltedtotalcol <- full_join(meltedtotal, PDpercell, by=c("cell.col"="Unique_ID"))
head(meltedtotalcol)
str(meltedtotalcol$PD)
meltedtotalcol$uniquePDpercellcol <- meltedtotalcol$totalshared - meltedtotalcol$PD
head(meltedtotalcol)
meltedtotalcol$uniquePDpercellcol[meltedtotalcol$uniquePDpercellcol < 0] <- 0 #get rid of the negative values

#join the database, total.shared pd with the pd per cell
meltedtotalrow <- full_join(meltedtotal, PDpercell, by=c("cell.row"="Unique_ID"))
head(meltedtotalrow)
str(meltedtotalrow$PD)
meltedtotalrow$uniquePDpercellrow <- meltedtotalrow$totalshared - meltedtotalrow$PD
head(meltedtotalrow)
meltedtotalrow$uniquePDpercellrow[meltedtotalrow$uniquePDpercellrow < 0] <- 0 #get rid of the negative values

#find the min between the meltedtotalrow and the meltedtotalcol
meltedtotalboth <- full_join(meltedtotalrow, meltedtotalcol, by=c("cell.row"="cell.row","cell.col"="cell.col"))
head(meltedtotalboth)
meltedtotalboth$minvalue <- with(meltedtotalboth, pmin(uniquePDpercellrow, uniquePDpercellcol))

#transform back to a matrix
mincol <- dcast(meltedtotalboth[, c("cell.row", "cell.col","minvalue")], cell.row~cell.col)
rownames(mincol)<- mincol$cell.row
mincol <- as.matrix(mincol[,-1])
#write.table(mincol, "C:/C.drive/Chapter 4/RoutputFiles/rayphylo_min.csv", sep = ",", 
#            row.names = TRUE, col.names = TRUE)
write.table(mincol, "/home/lnkdavi/scratch/betasharksandrays/sharksphylo_min.csv", sep = ",", 
            row.names = TRUE, col.names = TRUE)
head(mincol)

###########################################################
# #NOW CALCULATE BETA DIVERSITY
# ###########################################################
# #shared PD, minimum value of unique PD between cells, same results but from different packages, 
bsimPDpicante <- 1-(pdcombined / (mincol + pdcombined))
bsimPDpicante2<- as.dist(bsimPDpicante)
#bsimPDbetapart <- 1-(shared / (min + shared))
#write.table(bsimPDpicante, "C:/C.drive/Chapter 4/RoutputFiles/rayphylo_bsimPDpicante.csv", sep = ",", 
#            row.names = TRUE, col.names = TRUE)

#write.table(bsimPDpicante, "/home/lnkdavi/scratch/betasharksandrays/rayphylo_bsimPDpicante.csv", sep = ",", 
#            row.names = TRUE, col.names = TRUE)

x<- as.numeric(round(quantile(1:nrow(bsimPDpicante)),0))
x[2]

write.table(bsimPDpicante[x[1]:x[2],], "/home/lnkdavi/scratch/betasharksandrays/sharkandrayphylo_bsimPDpicante1.csv", sep = ",", 
            row.names = TRUE, col.names = TRUE)
write.table(bsimPDpicante[x[2]:x[3],], "/home/lnkdavi/scratch/betasharksandrays/sharkandrayphylo_bsimPDpicante2.csv", sep = ",", 
            row.names = TRUE, col.names = TRUE)
write.table(bsimPDpicante[x[3]:x[4],], "/home/lnkdavi/scratch/betasharksandrays/sharkandrayphylo_bsimPDpicante3.csv", sep = ",", 
            row.names = TRUE, col.names = TRUE)
write.table(bsimPDpicante[x[4]:x[5],], "/home/lnkdavi/scratch/betasharksandrays/sharkandrayphylo_bsimPDpicante4.csv", sep = ",", 
            row.names = TRUE, col.names = TRUE)

write.table(bsimPDpicante, "/home/lnkdavi/scratch/betasharksandrays/sharkandrayphylo_bsimPDpicante.csv", sep = ",", 
            row.names = TRUE, col.names = TRUE)

# #########################################################################
# #NOW CLUSTER THE CELLS 
# #########################################################################
# #dendro1 = hclust(bsimdistPD, method = "average")# method average is UPGMA
# #class(dendro1)
# #summary(dendro1)
# #plot(dendro1)
# 
