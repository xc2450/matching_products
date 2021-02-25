

rm(list=ls())

library(readxl)
library(dplyr) 
library(bit64)
library(data.table)
library(descr)
library(readr)
library(checkmate)
library(lme4)
library(ggplot2)
library(gganimate)
library(gifski)
library(png)
library(transformr)
library(zoo)
library(foreach)
library(reshape2)
library(tidyr)
library(car)
library(geobr)
library(scales)
library(maptools)
library(RColorBrewer)
library(fuzzyjoin)
library(stringdist)
library(Rfast)


base <- read_excel("amostra.xlsx",skip=2)

base1 <- base %>% select(c(CODIGOLICITACAO,num_seq_item,descr_item,nm_un_fornec))

base2 <- base %>% select(c(CODIGOLICITACAO,num_seq_item,descr_prod,unid_comercial,cod_gtin,cod_ncm))

base1 <- base1 %>% mutate(cod_lic_item=paste(CODIGOLICITACAO,num_seq_item,sep="_")) %>%
  distinct(cod_lic_item,.keep_all = T) %>% select(-cod_lic_item)

base2 <- base2 %>% group_by(CODIGOLICITACAO) %>% mutate(min_seq= min(num_seq_item)) %>%
  filter(num_seq_item == min_seq) %>% select(-num_seq_item)

lics1 <- base1 %>% group_by(CODIGOLICITACAO) %>% summarise(n1=n())
lics2 <- base2 %>% group_by(CODIGOLICITACAO) %>% summarise(n2=n())

lics_c <- full_join(lics1,lics2) %>% mutate(dif_n = n1-n2) %>% filter(dif_n == 0)
rm(lics1,lics2)

base1 <- base1 %>% right_join(lics_c) %>% select(-c(n1,n2,dif_n)) %>% mutate(index1=descr_item,
                                                                             index1_1=nm_un_fornec)
base2 <- base2 %>% right_join(lics_c) %>% select(-c(n1,n2,dif_n)) %>% mutate(index2=descr_prod,
                                                                             index2_1=unid_comercial) 


base2$index2 <- gsub("Á","A",base2$index2)
base2$index2 <- gsub("À","A",base2$index2)
base2$index2 <- gsub("Ã","A",base2$index2)
base2$index2 <- gsub("Â","A",base2$index2)
base2$index2 <- gsub("É","E",base2$index2)
base2$index2 <- gsub("È","E",base2$index2)
base2$index2 <- gsub("Ê","E",base2$index2)
base2$index2 <- gsub("Í","I",base2$index2)
base2$index2 <- gsub("Ì","I",base2$index2)
base2$index2 <- gsub("Î","I",base2$index2)
base2$index2 <- gsub("Ó","O",base2$index2)
base2$index2 <- gsub("Ò","O",base2$index2)
base2$index2 <- gsub("Õ","O",base2$index2)
base2$index2 <- gsub("Ô","O",base2$index2)
base2$index2 <- gsub("Ú","U",base2$index2)
base2$index2 <- gsub("Ù","U",base2$index2)
base2$index2 <- gsub("Û","U",base2$index2)
base2$index2 <- gsub("Ç","C",base2$index2)



base1$index1 <- gsub(" ","",base1$index1)
base2$index2 <- gsub(" ","",base2$index2)

base1$index1_1 <- gsub(" ","",base1$index1_1)
base2$index2_1 <- gsub(" ","",base2$index2_1)

methods <- c("osa","lv","dl","lcs","qgram","cosine","jaccard","jw")

base_f <- data.frame()

for(y in 1:length(methods)){
base_f1 <- data.frame()

for(x in lics_c$CODIGOLICITACAO){
  
  base_n1 <- base1 %>% filter(CODIGOLICITACAO==x) 
  base_n2 <- base2 %>% filter(CODIGOLICITACAO==x)
  
  base_n <- full_join(base_n1,base_n2)
  
  base_n$distancia1 <- stringdist(base_n$index1,base_n$index2, method=methods[y])
  base_n$distancia2 <- stringdist(base_n$index1_1,base_n$index2_1, method=methods[y])
  
  base_n <- base_n %>% arrange(distancia1)

base_n <- base_n %>% ungroup() %>% mutate(rank1 = rank(distancia1)) %>%
  group_by(rank1) %>% mutate(rank2 = rank(distancia2)) %>% ungroup() %>% mutate(rank=rank1*100+rank2)

for(z in 1:nrow(base_n1)){
base_n <- base_n %>% filter(index1 != base_n$index1[which(base_n$rank == Rfast::nth(base_n$rank,z, descending = F))] | rank == Rfast::nth(base_n$rank,z, descending = F)) %>%
  filter(index2 != base_n$index2[which(base_n$rank == Rfast::nth(base_n$rank,z, descending = F))] | rank == Rfast::nth(base_n$rank,z, descending = F))
}
  
  base_n <- base_n %>% group_by(num_seq_item) %>% slice_min(distancia1)   
  base_n <- base_n %>% group_by(num_seq_item) %>% slice_min(distancia2)   
  
  base_n <- base_n %>%
    group_by(CODIGOLICITACAO,index2) %>%
    mutate(dup=ifelse(duplicated(index2),1,0),
           dup=max(dup))
  
  base_n <- base_n %>%
    group_by(CODIGOLICITACAO,num_seq_item) %>%
    mutate(dup2=ifelse(duplicated(num_seq_item),1,0),
           dup2=max(dup2))

  base_n <- base_n %>% mutate(dup_un=case_when(
    dup == 1 ~ 1,
    dup2 == 1 ~ 1,
    TRUE ~ 0
  ))
    
  
  base_f1 <- bind_rows(base_f1,base_n)
  
  base_f1 <- base_f1 %>% mutate(method=methods[y])
  
  rm(base_n1,base_n2,base_n)
  

}
base_f <- bind_rows(base_f,base_f1)
rm(base_f1)

}

probs <- base_f %>% group_by(method) %>% summarise(probs=sum(dup,dup2),
                                                   prob1=sum(dup),
                                                   prob2=sum(dup2),
                                                   good=n()-sum(dup_un))

base_f_mf <- base_f %>% filter(method==probs$method[which(probs$probs==min(probs$probs))])
