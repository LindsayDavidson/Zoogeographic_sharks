 #ray species per cluster

raycomm <- read.csv("C:/C.drive/Chapter 4/Species_within_clusters/SJ_Hex4degChon_171005_trimmed_nopelagic.csv")
head(raycomm)

tax <- read.csv("C:/C.drive/Chapter 4/PhylogenyChris/SharkTaxonomyForDan.csv")
head(tax)

raycommtax <- inner_join(raycomm, tax[,c(1,2,3,4,5,7)], by = c("binomial"="IUCN_Scientific_name"))
#check for duplicates
head(raycommtax)
raycommtax2 <- filter(raycommtax, Superorder =="Batoidea")

length(unique(raycommtax2$binomial))

raycommdup <- raycommtax2[!duplicated(raycommtax2[c("Unique_ID", "binomial")]),]

raycommdup <- raycommdup %>% group_by(rayPDclusterID) %>% dplyr::mutate(cellspercluster = length(Unique_ID))
head(raycommdup)

raycommdup <- raycommdup %>% group_by(rayPDclusterID,binomial ) %>% dplyr::mutate(speciescount = length(binomial))
head(raycommdup)

raycommdup <- raycommdup %>% group_by(rayPDclusterID,Family ) %>% dplyr::mutate(familycount = length(Family))
head(raycommdup)

raycommdup <- raycommdup %>% group_by(rayPDclusterID,Order ) %>% dplyr::mutate(ordercount = length(Order))
head(raycommdup)


#top 5 species per cluster
raycounts <- raycommdup[!duplicated(raycommdup[c( "binomial", "rayPDclusterID")]),]
#how find the top ten species per cluster, in terms of column speciescount
maxspecies <- raycounts %>% arrange_(~desc(speciescount)) %>%
  group_by(rayPDclusterID) %>%
  top_n(n = 5, wt = speciescount)
head(maxspecies)

#top 5 families per cluster
raycountsfam <- raycommdup[!duplicated(raycommdup[c( "Family", "rayPDclusterID")]),]
#how find the top ten species per cluster, in terms of column speciescount
maxfamily <- raycountsfam %>% arrange_(~desc(familycount)) %>%
  group_by(rayPDclusterID) %>%
  top_n(n = 5, wt = familycount)
head(maxfamily)

#top 5 orders per cluster
raycountsord <- raycommdup[!duplicated(raycommdup[c( "Order", "rayPDclusterID")]),]
#how find the top ten species per cluster, in terms of column speciescount
maxorder <- raycountsord %>% arrange_(~desc(ordercount)) %>%
  group_by(rayPDclusterID) %>%
  top_n(n = 5, wt = ordercount)
head(maxorder)

#join together
final <- cbind(maxfamily[,c("rayPDclusterID","Family", "familycount")], 
               maxspecies[,c("rayPDclusterID","binomial", "speciescount")])#, 
               #maxorder[,c("rayPDclusterID", "Order", "ordercount")])
head(final)
finalorder <- maxorder[,c("rayPDclusterID", "Order", "ordercount")]
head(finalorder)

write.csv(finalorder, "C:/C.drive/Chapter 4/Species_within_clusters/RaysPD_SumofOrders.csv")
write.csv(final, "C:/C.drive/Chapter 4/Species_within_clusters/RaysPD_sumofspeciesfamilies.csv")
