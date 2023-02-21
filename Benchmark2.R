#####
#Libraries loading
library(tidyverse)
library(factoextra)
##### 
par(mfrow=c(2,2))

#IntaRNA

#separate antisense from intergenic
inta_AS <- inta_s %>% filter(Relative.position.to.CDSs == "AS")
inta_OT <- inta_s %>% filter(Relative.position.to.CDSs != "AS")

#AS
hist(inta_AS$E, breaks = 25,
     main = "A",
     xlab = "Hybridization energy") #before normalization

hist(inta_AS$E / inta_AS$`sRNA length`, breaks = 100,
     main = "A",
     xlab = "Hybridization energy") #after normalization

#Other
hist(inta_OT$E, breaks = 25,
     main = "Intergenic & RBS sRNAs: hybridization energy (pre-normalized) values - IntaRNA",
     xlab = "Hybridization energy") #before normalization

hist(inta_OT$E / inta_OT$`sRNA length`, breaks = 100,
     main = "B",
     xlab = "Hybridization energy") #after normalization


#####
#sTar

#separate antisense from intergenic
sTar_AS <- sTar_s %>% filter(Relative.position.to.CDSs == "AS")
sTar_OT <- sTar_s %>% filter(Relative.position.to.CDSs != "AS")

#AS
hist(sTar_AS$E, breaks = 25,
     main = "Antisense sRNAs: hybridization energy (pre-normalized) values - IntaRNAsTar",
     xlab = "Hybridization energy") #before normalization

hist(sTar_AS$E / sTar_AS$`sRNA length`, breaks = 100,
     main = "A",
     xlab = "Hybridization energy") #after normalization

#Other
hist(sTar_OT$E, breaks = 25,
     main = "Intergenic & RBS sRNAs: hybridization energy (pre-normalized) values - IntaRNAsTar",
     xlab = "Hybridization energy") #before normalization


hist(sTar_OT$E / sTar_OT$`sRNA length`, breaks = 100,
     main = "B",
     xlab = "Hybridization energy") #after normalization

#####
#sRNAfTarget

#separate antisense from intergenic
srnar_AS <- srnar_s %>% filter(Relative.position.to.CDSs == "AS")
srnar_OT <- srnar_s %>% filter(Relative.position.to.CDSs != "AS")

#AS
hist(srnar_AS$E, breaks = 40,
     main = "Hybridization energy before normalization - sRNAfTarget",
     xlab = "Hybridization energy")#before normalization

hist(srnar_AS$E / srnar$`sRNA length`, breaks = 40,
     main = "A",
     xlab = "Hybridization energy")#after normalization

#Other sRNAs
hist(srnar_OT$E, breaks = 100,
     main = "Hybridization energy before normalization - sRNAfTarget",
     xlab = "Hybridization energy") #before normalization

hist(srnar_OT$E / srnar_OT$`sRNA length`, breaks = 100,
     main = "B",
     xlab = "Hybridization energy") #after normalization: as expected there's a significant variation

#####
#Clustering analysis, using unsupervised k-means on 
#top40 data, by taking a sample of 5000 records
set.seed(400)
inta_kmeans <- inta_s %>% select("id2","E", "sRNA length")
inta_kmeans <- inta_kmeans[sample(nrow(inta_kmeans), 5000),] #taking 5000 samples
inta_kmeans <- inta_kmeans[,-1]
inta_kmeans_scale <- scale(inta_kmeans) #scale data

inta_kmeans_data <- dist(inta_kmeans_scale)#distance matrix

fit <- kmeans(inta_kmeans_scale, centers = 3, nstart = 100)
fviz_cluster(fit, data = inta_kmeans_scale, geom = "point",
             main = "B") #first try on pre-assumption of three clusters

#calculate how many clusters
fviz_nbclust(inta_kmeans_scale, kmeans, method = "wss") + 
  labs(subtitle = "Elbow method") +
  geom_vline(xintercept = 4, linetype = 2) +
  geom_vline(xintercept = 3, linetype = 1, color = "red") +
  ggtitle("A. Optimal number of clusters") #we selected 5 centroids 

fit2 <- kmeans(inta_kmeans_scale, centers = 5, nstart = 100) #fitting model
fviz_cluster(fit2, data = inta_kmeans_scale, geom = "point", main = "C")
