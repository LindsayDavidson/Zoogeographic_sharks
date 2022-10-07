 #shark species per cluster

sharkcomm <- read.csv("C:/C.drive/Chapter 4/Species_within_clusters/SJ_Hex4degChon_171005_sharkstrimmed.csv")
head(sharkcomm)

sharkcluster <- read.csv("C:/C.drive/Chapter 4/Species_within_clusters/sharkPD_clusters.csv")
head(sharkcluster)

tax <- read.csv("C:/C.drive/Chapter 4/PhylogenyChris/SharkTaxonomyForDan.csv")
head(tax)


sharkcommtax <- inner_join(sharkcomm, tax[,c(1,2,3,4,5,7)], by = c("binomial"="IUCN_Scientific_name"))
#check for duplicates
head(sharkcommtax)
unique(sharkcommtax$Superorder)
sharkcommtax3 <- filter(sharkcommtax, Superorder !="Batoidea")
sharkcommtax2 <- filter(sharkcommtax3, Superorder !="Holocephali")

#now join in the cluster ID
sharkcommtax3 <- inner_join(sharkcluster, sharkcommtax2, by = c("Unique_ID"="Unique_ID") )

length(unique(sharkcommtax3$binomial))
length(unique(sharkcommtax3$X12))

head(sharkcommtax3)

sharkcommdup <- sharkcommtax3[!duplicated(sharkcommtax3[c("Unique_ID", "binomial")]),]

sharkcommdup <- sharkcommdup %>% group_by(X12) %>% dplyr::mutate(cellspercluster = length(Unique_ID))
head(sharkcommdup)

sharkcommdup <- sharkcommdup %>% group_by(X12,binomial ) %>% dplyr::mutate(speciescount = length(binomial))
head(sharkcommdup)

sharkcommdup <- sharkcommdup %>% group_by(X12,Family ) %>% dplyr::mutate(familycount = length(Family))
head(sharkcommdup)

sharkcommdup <- sharkcommdup %>% group_by(X12,Order ) %>% dplyr::mutate(ordercount = length(Order))
head(sharkcommdup)


#top 5 species per cluster
sharkcounts <- sharkcommdup[!duplicated(sharkcommdup[c( "binomial", "X12")]),]
#how find the top ten species per cluster, in terms of column speciescount
maxspecies <- sharkcounts %>% arrange_(~desc(speciescount)) %>%
  group_by(X12) %>%
  top_n(n = 5, wt = speciescount)
head(maxspecies)

#top 5 families per cluster
sharkcountsfam <- sharkcommdup[!duplicated(sharkcommdup[c( "Family", "X12")]),]
#how find the top ten species per cluster, in terms of column speciescount
maxfamily <- sharkcountsfam %>% arrange_(~desc(familycount)) %>%
  group_by(X12) %>%
  top_n(n = 5, wt = familycount)
head(maxfamily)

#top 5 orders per cluster
sharkcountsord <- sharkcommdup[!duplicated(sharkcommdup[c( "Order", "X12")]),]
#how find the top ten species per cluster, in terms of column speciescount
maxorder <- sharkcountsord %>% arrange_(~desc(ordercount)) %>%
  group_by(X12) %>%
  top_n(n = 5, wt = ordercount)
head(maxorder)

#join together
final <- cbind(maxfamily[,c("X12","Family", "familycount")], 
               maxspecies[,c("X12","binomial", "speciescount")])#, 
               #maxorder[,c("sharkPDclusterID", "Order", "ordercount")])
head(final)

family <- maxfamily[,c("X12","Family", "familycount")]
species <- maxspecies[,c("X12","binomial", "speciescount")]
finalorder <- maxorder[,c("X12", "Order", "ordercount")]
head(finalorder)

write.csv(family, "C:/C.drive/Chapter 4/Species_within_clusters/sharksPD_Sumoffamilies.csv")
write.csv(species, "C:/C.drive/Chapter 4/Species_within_clusters/sharksPD_Sumofspecies.csv")
write.csv(finalorder, "C:/C.drive/Chapter 4/Species_within_clusters/sharksPD_sumoforders.csv")
