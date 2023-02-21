#####
#Libraries loading
library(readr)
library(tidyverse)
library(Biostrings)
library(stringr)
library(msa)
#####
#Data upload + gene list creation
MTB_Resistance_Mediating <- read_delim("MTB_Resistance_Mediating.csv", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)

lista.geni <- unique(MTB_Resistance_Mediating$`Gene ID`)[-107] #gene list

#####
#Interaction list - IntaRNA
inta_filt_res <- inta[inta$id1%in%lista.geni,]
top1pergene_inta <- inta_filt_res %>% group_by(id1) %>% mutate(ranks=rank(E)) %>% 
  arrange(id1,ranks) %>% group_by(id2) %>% filter(ranks==1)

write.csv2(top1pergene_inta, file = "Tabella geni e candidati.csv") #creazione tabella

#Distributions pre-normalization
hist(top1pergene_inta$E, breaks = 50, main = "Pre-normalization per length", 
     xlab = "H.E.")

#Normalization + subsequent distribution
hist(top1pergene_inta$E / top1pergene_inta$`sRNA length`, breaks = 50,
     main = "Post-normalization per length", xlab = "H.E. normalized")

#####
#Interaction list - IntaRNAsTar
star_filt_res <- sTar[sTar$id1%in%lista.geni,]
top1pergene_star <- star_filt_res %>% group_by(id1) %>% mutate(ranks=rank(E)) %>% 
  arrange(id1,ranks) %>% group_by(id1) %>% filter(ranks==1)

write.csv2(top1pergene_star, file = "DRS sTar.csv")

#Distributions pre-normalization
hist(top1pergene_star$E, breaks = 50, main = "Pre-normalization per length", 
     xlab = "H.E.")

#Normalization + subsequent distribution
hist(top1pergene_star$E / top1pergene_star$`sRNA length`, breaks = 50,
     main = "Post-normalization per length", xlab = "H.E. normalized")

#####
#Multiple sequence alignment
seqs <- readRNAStringSet("cand_drs_pred.fasta", format = "fasta")
alignment <- msa(seqs)
print(alignment, show = "complete")

msaPrettyPrint(alignment, output="pdf", showNames="none",
               showLogo="none", askForOverwrite=FALSE, verbose=FALSE)

#####
#extract DRS interactions predicted from IntaRNA top40
#we need these .csv to produce coverage and region plots per each candidate

cand1405 <- top1pergene_inta %>% filter(id2 == "candidate_1405")
write.csv2(cand1405, file = "cand1405.csv")
listageniper1405 <- cand1405$id1
cand2023 <- top1pergene_inta %>% filter(id2 == "candidate_2023")
write.csv2(cand2023, file = "cand2023.csv")
listageniper2023 <- cand2023$id1
cand2223 <- top1pergene_inta %>% filter(id2 == "candidate_2223")
write.csv2(cand2223, file = "cand2223.csv")
cand63 <- top1pergene_inta %>% filter(id2 == "candidate_63")
write.csv2(cand63, file = "cand63.csv")
cand899 <- top1pergene_inta %>% filter(id2 == "candidate_899")
write.csv2(cand899, file = "cand899.csv")

#annex table, at the end of thesis
tabella_annex <- inta_s %>% filter(id2 == "candidate_63" | id2 =="candidate_899" | id2 == "candidate_1405" | id2=="candidate_2023" | id2 == "candidate_2223") 
write.csv2(tabella_annex, file = "tabella_annex.csv")



