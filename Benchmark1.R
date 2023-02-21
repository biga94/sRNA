#####
#Libraries loading
library(tidyverse)
library(readr)
library(EnvStats)
######## 
#Data import

#IntaRNA
inta <- read_delim("Predictions/res_inta2.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE)

inta$id1 <- gsub(".*-","", inta$id1)
inta$id2 <- str_remove(inta$id2,"\\:.*")

#IntaRNAsTar
sTar <- read_delim("Predictions/res_tar.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE)

sTar$id1 <- gsub(".*-","", sTar$id1)
sTar$id2 <- str_remove(sTar$id2,"\\:.*")

#sRNARFTarget
srnar <- read.delim('Predictions/Prediction_probabilities.csv',header=T,
                    stringsAsFactors = F,sep='\t') %>% 
  rename(id1=mRNA_ID,id2=sRNA_ID,E=Prediction_Probability) %>% 
  mutate(E=1-E,alg='srnar')

srnar$id1 <- gsub(".*-","",srnar$id1)
srnar$id2 <- str_remove(srnar$id2,"\\:.*")

#####

#Shortlist

ShortListed <- read.csv2("ShortListed.csv", header = T, sep = ";") #importo shortlist
ShortListed <- ShortListed %>% rename(Candidate = Candidate.ID) %>% rename("sRNA length" = Length.1) #cambio nome variabili

#sottraggo informazioni dalla shortlist riguardo la lunghezza, da inserire nei dataset predittivi
cand_length <- ShortListed %>% select(Candidate, `sRNA length`) %>% rename(id2=Candidate)
inta <- inta %>% inner_join(cand_length, by="id2") #intarna
sTar <- sTar %>% inner_join(cand_length, by="id2") #sTar
srnar <- srnar %>% inner_join(cand_length, by="id2")

#tipo di sRNA 
cand_type <- ShortListed %>% select(Candidate, Relative.position.to.CDSs) %>% rename(id2=Candidate)
inta <- inta %>% inner_join(cand_type, by="id2") #intarna
sTar <- sTar %>% inner_join(cand_type, by="id2") #sTar
srnar <- srnar %>% inner_join(cand_type, by="id2") #sRNARFTarget

#####
#Filtering dati: provo a selezionare top40 e lavorare sui pvalue

#Selezionato i top40 per il filtering. Da cambiare il cutoff
inta_s <- inta %>% group_by(id2) %>% filter(rank(E) < 41)
srnar_s <- srnar %>% group_by(id2) %>% filter(rank(E) < 41) %>% relocate(id2, .after = "id1")
sTar_s <- sTar %>% group_by(id2) %>% filter(rank(E) < 41)

inta_filt <- inta %>% filter(p.value < 0.05) #Selezionato le top interazioni per valore di pvalue. cut-off p < 0.05

#####
#Comparison degli algoritmi

#creazione dataset con tutti e tre gli algoritmi
ISS <- rbind(inta_s,srnar_s,sTar_s)

#quanti sRNA sono presenti in 
shared.ISS <- ISS %>% select(-alg, -E) %>%  
  group_by(id1, id2) %>% filter(n()==1) %>% 
  distinct() %>% group_by(id2) %>% summarise(n=n())

mean(shared.ISS$n)
boxplot(shared.ISS$n, main = "Condivisi tra i 3 algoritmi")
hist(shared.ISS$n)

#Provo a farlo su coppie di algoritmi

ISr <- rbind(inta_s, srnar_s)
ISt <- rbind(inta_s, sTar_s)
SrSt <- rbind(srnar_s, sTar_s)

#intarna + srnarftarget
shared.ISr <- ISr %>% select(-alg, -E) %>% 
  group_by(id1,id2) %>% filter(n()==2) %>% 
  distinct() %>% group_by(id2) %>% summarise(n=n())

boxplot(shared.ISr$n, main = "Condivisi tra Intarna e sRNArfTarget")
hist(shared.ISr$n)

### intarna + sTar

shared.ISt <- ISt %>% select(-alg, -E) %>% 
  group_by(id1, id2) %>% filter(n()==2) %>% 
  distinct() %>% group_by(id2) %>% summarise(n=n())

mean(shared.ISt$n)
boxplot(shared.ISt$n, main="Condivisi tra Intarna e sTar")
hist(shared.ISt$n)

### srnarftarget + sTar

shared.SrSt <- SrSt %>% select(-alg, -E) %>% 
  group_by(id1, id2) %>% filter(n()==2) %>% 
  distinct() %>% group_by(id2) %>% summarise(n=n())

mean(shared.SrSt$n)
boxplot(shared.SrSt$n, main="Condivisi tra sRNArfTarget e sTar")
hist(shared.SrSt$n)

#Ranking
#Rank effettuato da raw data, ogni sRNA è stato rankato per ogni target

inta_rank <- inta %>% group_by(id2) %>% mutate(ranks = rank(E)) %>% arrange(id2)
srnar_rank <- srnar %>% group_by(id2) %>% mutate(ranks=rank(E)) %>% arrange(id2)
sTar_rank <- sTar %>% group_by(id2) %>% mutate(ranks=rank(E)) %>% arrange(id2,ranks)

#Inta + srnar
joined_alg <- inta_rank %>% rbind(srnar_rank) %>% group_by(id1, id2) %>% filter(n() == 2) #shared

joined_alg %>% group_by(id1, id2) %>% summarise(new_rank=mean(ranks)) %>%
  group_by(id2) %>% mutate(rank=rank(new_rank)) %>% 
  filter(rank <41) -> joined_ISr

#Inta + sTar
joined_alg2 <- inta_rank %>% rbind(sTar_rank) %>% group_by(id1, id2) %>% filter(n() == 2)

joined_alg2 %>% group_by(id1, id2) %>% summarise(new_rank=mean(ranks)) %>% 
  group_by(id2) %>% mutate(rank=rank(new_rank)) %>% 
  filter(rank < 41) %>% relocate(id2, .after="id1") -> joined_ISt

#srnarftarget + star
joined_alg3 <- srnar_rank %>% rbind(sTar_rank) %>% group_by(id1, id2) %>% filter(n()==2)

joined_alg3 %>% group_by(id1,id2) %>% summarise(new_rank=mean(ranks)) %>% 
  group_by(id2) %>% mutate(rank=rank(new_rank)) %>% 
  filter(rank < 41) -> joined_SrSt

#####
#modelli GEV
inta_AS <- inta %>% filter(Relative.position.to.CDSs == "AS" | Relative.position.to.CDSs == "AS to 5'/3' UTRs")
fit01 <- egevd((inta_AS$E / inta_AS$`sRNA length`) * -1, ci = T, conf.level = 0.95)
inta_notAS <- inta %>% filter(Relative.position.to.CDSs != "AS" & Relative.position.to.CDSs!= "AS to 5'/3' UTRs")
fit02 <- egevd((inta_notAS$E / inta_notAS$`sRNA length`) * -1, ci = T, conf.level = 0.95)

star_as <- sTar %>% filter(Relative.position.to.CDSs == "AS" | Relative.position.to.CDSs == "AS to 5'/3' UTRs")
fit03 <- egevd((star_as$E / star_as$`sRNA length`) * -1, ci = T, conf.level = 0.95)
star_notas <- sTar %>% filter(Relative.position.to.CDSs != "AS" & Relative.position.to.CDSs!= "AS to 5'/3' UTRs")
fit04 <- egevd((star_notas$E / star_notas$`sRNA length`) * -1, ci = T, conf.level = 0.95)
