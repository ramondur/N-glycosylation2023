###PREPARATION DATASET HUMAN PROTEOME FINAL 4 of July 2022
####################
####LOAD PACKAGES###
####################
library(dplyr)
library(ggplot2)
library(rlist)
library(tidyr)
library(ggpubr)
library(stringr)
library(seqinr)

###############################
####PREPARE MASTER DATA SET####
###############################
##MASTER DATASET##
#Merge information from output TANGO and correct protein IDs
tango_APRs <- read.delim("Data/Human/TANGO_output/PSX_tangoresidue.out", header=TRUE)
colnames(tango_APRs) <- c("Protein","Position","Residue","tango_Score","length")
tango_APRs$length <- NULL
list_names <- strsplit(as.character(tango_APRs$Protein), "\\|")
res <- matrix(nrow = length(tango_APRs$Protein), ncol = 1)
for (i in 1:length(tango_APRs$Protein)){
  res[i,] <- list_names[[i]][2]
}
tango_APRs$Protein <- res

##FILTERING
#Remove redundant proteins with CD HIT 90
#cd-hit using the web server with default options
names90 <- read.delim("Data/Human/UniProt/human_90.fasta", header=FALSE)
list_names <- strsplit(as.character(names90$V1), "\\|")
res <- matrix(nrow = length(names90$V1), ncol = 1)
for (i in 1:length(names90$V1)){
  res[i,] <- list_names[[i]][2]
}
names90$V1 <- res
names90 <- unique(names90)
tango_APRs <- subset(tango_APRs, Protein %in% names90$V1) #from 20344 to 19408
#Remove proteins with less than 25 aa and more that 10,000 aa and define APRs score >10
tango_APRs$length <- 1
proteinsum <- aggregate(length ~ Protein, data=tango_APRs, FUN="sum") #19408 proteins
colnames(proteinsum) <- c("Protein","protein_length")
tango_APRs <- merge(tango_APRs,proteinsum, by="Protein")
length(which(proteinsum$protein_length < 25)) # 29 < 25 aa
shorties <- tango_APRs[tango_APRs$Protein %in% proteinsum$Protein[proteinsum$protein_length<25],]
master_df <- tango_APRs[!(tango_APRs$Protein %in% shorties$Protein),]
##SAVE HUMAN PROTEOME
human_proteome <- aggregate (Residue ~ Protein, data=master_df, FUN="paste",collapse = "")
colnames(human_proteome) <- c("Protein","Sequence")
save(human_proteome, file = "Data/RData/Human/human_proteome.RData")
#Create files for SignalP and TMHMM
write.fasta(as.list(human_proteome$Sequence[1:5000]), human_proteome$Protein[1:5000], "Data/Human/UniProt/1-5000.fa", open = "w", nbchar = 60, as.string = FALSE)
write.fasta(as.list(human_proteome$Sequence[5001:10000]), human_proteome$Protein[5001:10000], "Data/Human/UniProt/5001-10000.fa", open = "w", nbchar = 60, as.string = FALSE)
write.fasta(as.list(human_proteome$Sequence[10001:15000]), human_proteome$Protein[10001:15000], "Data/Human/UniProt/10001-15000.fa", open = "w", nbchar = 60, as.string = FALSE)
write.fasta(as.list(human_proteome$Sequence[15001:nrow(human_proteome)]), human_proteome$Protein[15001:nrow(human_proteome)], "Data/Human/UniProt/15001-end.fa", open = "w", nbchar = 60, as.string = FALSE)

#Transmembrane domains from TMHMM (deepTMHMM using biolib webserver)
trans1.5000 <- read.delim("Data/TMHMM/Human/trans1-5000.txt", header=FALSE)
trans5001.10000 <- read.delim("Data/TMHMM/Human/trans5001-10000.txt", header=FALSE)
trans10001.15000 <- read.delim("Data/TMHMM/Human/trans10001-15000.txt", header=FALSE)
transend <- read.delim("Data/TMHMM/Human/15001-end.txt", header=FALSE)
TMHMM <- rbind(trans1.5000,trans5001.10000,trans10001.15000,transend)
TMHMM <- subset(TMHMM, V5 != "PredHel=0" & V1 %in% human_proteome$Protein, select = c(1,6))
colnames(TMHMM) <- c("Protein","Domains")
TMHMM$Domains <- as.character(TMHMM$Domains)
list_names <- strsplit(TMHMM$Domains, "=")
res <- matrix(nrow = length(TMHMM$Domains), ncol = 1)
for (i in 1:length(TMHMM$Domains)){
  res[i,] <- list_names[[i]][2]
}
TMHMM$Domains <- res
s <- strsplit(TMHMM$Domains, split = "i|o")
TMHMM <- data.frame(Protein = rep(TMHMM$Protein, sapply(s, length)), Domains = unlist(s))
TMHMM$Domains <- as.character(TMHMM$Domains)
TMHMM <-subset(TMHMM, Domains != '')
list_names <- strsplit(as.character(TMHMM$Domains), "-")
res1 <- matrix(nrow = length(TMHMM$Domains), ncol = 1)
res2 <- matrix(nrow = length(TMHMM$Domains), ncol = 1)
for (i in 1:length(TMHMM$Domains)){
  res1[i,] <- list_names[[i]][1]
  res2[i,] <- list_names[[i]][2]
}
TMHMM$start <- res1
TMHMM$end <- res2
master_df$TMHMM <- "No"
master_df$TMHMM[which(master_df$Protein %in% TMHMM$Protein)] <- "Yes" #5120 proteins with at least one trasmembrane domain
master_df$TMHMM_domain <- "No"
for (i in 1:length(TMHMM$Domains)){
  prot <- as.character(TMHMM$Protein[i])
  transmembrane <-subset(master_df, Protein == prot)
  start <- TMHMM$start[i]
  end <- TMHMM$end[i]
  transmembrane$TMHMM_domain[start:end] <- "Yes"
  master_df$TMHMM_domain[which(master_df$Protein == as.character(TMHMM$Protein[i]))] <- transmembrane$TMHMM_domain
  print(i)
}
#Signal peptide (SignalP 6.0 using biolib webserver)
sigp1.5000 <- read.delim("Data/SignalP/Human/SigP_1-5000.txt", header=FALSE)
sigp5001.10000 <- read.delim("Data/SignalP/Human/SigP_5001-10000.txt", header=FALSE)
sigp10001.15000 <- read.delim("Data/SignalP/Human/SigP_10001-15000.txt", header=FALSE)
sigp15001.20000 <- read.delim("Data/SignalP/Human/Sigp_15001-20000.txt", header=FALSE)
sigpend <- read.delim("Data/SignalP/Human/SigP_20001-end.txt", header=FALSE)
SignalP <- rbind(sigp1.5000,sigp5001.10000,sigp10001.15000,sigp15001.20000,sigpend)
SignalP <- subset(SignalP, V1 %in% human_proteome$Protein)
error <- human_proteome[which(!human_proteome$Protein %in% SignalP$V1),]
error <- read.delim("Data/SignalP/Human/error.txt", header=FALSE)
SignalP <- rbind(SignalP,error)
SignalP <- unique(subset(SignalP, V2 == "SP(Sec/SPI)" & V1 %in% human_proteome$Protein, select = c(1,5)))
colnames(SignalP) <- c("Protein","Cleaveage_site")
site <- strsplit(as.character(SignalP$Cleaveage_site), " ")
res <- matrix(nrow = length(SignalP$Cleaveage_site), ncol = 1)
for (i in 1:length(SignalP$Cleaveage_site)){
  res[i,] <- site[[i]][3]
}
SignalP$Cleaveage_site <- res
site <- strsplit(as.character(SignalP$Cleaveage_site), "-")
res <- matrix(nrow = length(SignalP$Cleaveage_site), ncol = 1)
for (i in 1:length(SignalP$Cleaveage_site)){
  res[i,] <- site[[i]][1]
}
SignalP$Cleaveage_site <- res
master_df$SignalP <- "No"
for (i in 1:length(SignalP$Protein)){
  uniprot <- as.character(SignalP$Protein[i])
  protein <-subset(master_df, Protein == uniprot)
  site <- SignalP$Cleaveage_site[i]
  protein$SignalP[1:site] <- "Yes"
  master_df$SignalP[which(master_df$Protein == uniprot)] <- protein$SignalP
  print(i)
}
save(master_df, file = "Data/RData/TANGO_PTM_HUMAN_15/FILTERED/all_tango_res_filtered.RData")
##ANNOTATE CATEGORIES (APR, GKs)
#Define TANGO APRs ids
master_df$APRdef_tango <- "No"
master_df$APRdef_tango[which(master_df$tango_Score > 10)] <- "APR" 
tango_APRs <- subset(master_df, master_df$APRdef_tango=="APR")
j <- 1
tango_APRs$APRcount_tango <- 0
tango_APRs$APRcount_tango[1] <- 1
for (i in 2:length(tango_APRs$Protein)){
  if (tango_APRs$Protein[i] != tango_APRs$Protein[i-1] | tango_APRs$Position[i] != (tango_APRs$Position[i-1]+1)) {
    j <- j + 1
    tango_APRs$APRcount_tango[i] <- j 
  }
  else{
    tango_APRs$APRcount_tango[i] <- j 
  }
  print(i)
}
aprsum <- aggregate(length ~ APRcount_tango * Protein, data=tango_APRs, FUN="sum")
shortAPRs <- subset (aprsum, aprsum$length < 5) 
tango_APRs <- tango_APRs[-which(tango_APRs$APRcount_tango %in% shortAPRs$APRcount_tango),]
tango_APRs2 <- subset(tango_APRs, select = c(1,2,12))
master_df2 <- merge(master_df,tango_APRs2, by=c("Protein","Position"),all = TRUE)
master_df2$APRcount_tango[which(is.na(master_df2$APRcount_tango))] <- 0
master_df2$APRdef_tango[which(master_df2$APRcount_tango == 0)] <- "No"
onlyapr_tango <- tango_APRs

#Annotate TANGO categories and add avgTANGO per APR and maxTANGO protein
load("/Users/u0125634/Documents/PhD/Projects/PTM_project/Data/RData/TANGO_PTM_HUMAN_15/FILTERED/all_tango_res_annotated.RData")
annotateDF <- function(prepare_annotation, a,b,c,d) {
  prepare_annotation$APRdef2_tango <- gsub(a, b, prepare_annotation$APRdef2_tango)
  prepare_annotation$APRdef2_tango <- gsub(c, d, prepare_annotation$APRdef2_tango)
  return(prepare_annotation)
}
master_df3$APRdef2_tango <-as.numeric(as.factor(master_df3$APRdef_tango))
master_df3$APRdef2_tango[which(master_df3$APRdef2_tango == 2)] <- 0
prepare_annotation <- aggregate (APRdef2_tango ~ Protein, data=master_df3, FUN="paste",collapse = "")
prepare_annotation2 <- annotateDF(prepare_annotation,"01","21","10","12") #GK_1
prepare_annotation2 <- annotateDF(prepare_annotation2,"02","32","20","23") #GK_2
prepare_annotation2 <- annotateDF(prepare_annotation2,"03","43","30","34") #GK_3
prepare_annotation2 <- paste(prepare_annotation2$APRdef2_tango, collapse="")
prepare_annotation2 <- unlist(strsplit(prepare_annotation2,""))
master_df3$APRdef2_tango <- prepare_annotation2
master_df3$APRdef_tango[which(master_df3$APRdef2_tango == "0")] <- "Distal region"
master_df3$APRdef_tango[which(master_df3$APRdef2_tango == "2")] <- "GK_1"
master_df3$APRdef_tango[which(master_df3$APRdef2_tango == "3")] <- "GK_2"
master_df3$APRdef_tango[which(master_df3$APRdef2_tango == "4")] <- "GK_3"
master_df3$Protein <- as.factor(master_df3$Protein)
prepare_annotation3 <- annotateDF(prepare_annotation,"01","A1","10","1D") #GK_1
prepare_annotation3 <- annotateDF(prepare_annotation3,"0A","BA","D0","DE") #GK_2
prepare_annotation3 <- annotateDF(prepare_annotation3,"0B","CB","E0","EF") #GK_3
prepare_annotation3 <- paste(prepare_annotation3$APRdef2_tango, collapse="")
prepare_annotation3 <- unlist(strsplit(prepare_annotation3,""))
master_df3$Side <- prepare_annotation3
master_df3$Side[which(master_df3$Side %in% c("0","1"))] <- "None"
master_df3$Side[which(master_df3$Side %in% c("A","B","C"))] <- "N_ter"
master_df3$Side[which(master_df3$Side %in% c("D","E","F"))] <- "C_ter"
save(master_df3, file = "Data/RData/TANGO_PTM_HUMAN_15/FILTERED/all_tango_res_annotated2.RData")
#Summarize TANGO information
load("/Users/u0125634/Documents/PhD/Projects/PTM_project/Data/RData/TANGO_PTM_HUMAN_15/FILTERED/all_tango_res_annotated2.RData")
onlyapr_tango <- subset(master_df3, master_df3$APRdef_tango=="APR")
tango_APRs <- aggregate (Residue ~ APRcount_tango, data=onlyapr_tango, FUN="paste",collapse = "")
colnames(tango_APRs) <- c("APRcount_tango","APR")
TotalTANGO <- aggregate (tango_Score ~ APRcount_tango, data=onlyapr_tango, FUN="sum")
colnames(TotalTANGO) <- c("APRcount_tango","TotalTANGO")
lengthAPR <- aggregate (length ~ APRcount_tango, data=onlyapr_tango, FUN="sum")
colnames(lengthAPR) <- c("APRcount_tango","lengthAPR")
startAPR <- aggregate (Position ~ APRcount_tango, data=onlyapr_tango, FUN="min")
colnames(startAPR) <- c("APRcount_tango","startAPR")
endAPR <- aggregate (Position ~ APRcount_tango, data=onlyapr_tango, FUN="max")
colnames(endAPR) <- c("APRcount_tango","endAPR")
proteins_apr <- unique(subset(onlyapr_tango, select = c(1,12)))
signp_apr <- aggregate (SignalP ~ APRcount_tango, data=onlyapr_tango, FUN="unique")
tmhmm_apr <- aggregate (TMHMM_domain ~ APRcount_tango, data=onlyapr_tango, FUN="unique")
tango_APRs <- Reduce(function(x, y) merge(x, y, all=TRUE), list(proteins_apr,tango_APRs,TotalTANGO,lengthAPR,startAPR,endAPR,signp_apr,tmhmm_apr))
tango_APRs$avgTANGO <- tango_APRs$TotalTANGO / tango_APRs$lengthAPR
tango_APRs2 <- subset(tango_APRs, lengthAPR < 16 & SignalP == "No" & TMHMM_domain == "No") #From 114134 to 85989
onlyapr_tango <- onlyapr_tango[which(onlyapr_tango$APRcount_tango %in% tango_APRs2$APRcount_tango),]
len_APR <- aggregate (length ~ APRcount_tango, data=onlyapr_tango, FUN="sum")
onlyapr_tango$length <- NULL
onlyapr_tango <- merge(onlyapr_tango,len_APR)
master_df3$APRcount_tango[which(!master_df3$APRcount_tango %in% onlyapr_tango$APRcount_tango)] <- 0
master_df3$APRdef_tango[which(master_df3$APRcount_tango == 0)] <- "Distal region"
master_df3$APRdef2_tango[which(master_df3$APRcount_tango == 0)] <- 0
master_df3$Side[which(master_df3$APRcount_tango == 0)] <- "None"
total_scoreAPR <- aggregate (tango_Score ~ APRcount_tango, data=onlyapr_tango, sum)
colnames(total_scoreAPR) <- c("APRcount_tango","TotalScore")
onlyapr_tango <- merge(onlyapr_tango,total_scoreAPR,by="APRcount_tango")
onlyapr_tango$avgScore <- onlyapr_tango$TotalScore / onlyapr_tango$length
maxAPR_Protein <- aggregate (avgScore ~ Protein, data=onlyapr_tango, max)
colnames(maxAPR_Protein) <- c("Protein","maxProtscore")
onlyapr_tango <- merge(onlyapr_tango,maxAPR_Protein,by="Protein")
temp <- unique(subset(onlyapr_tango, select = c(1,2,20)))
temp2 <- unique(subset(onlyapr_tango, select = c(1,21)))
master_df3 <- merge(master_df3,temp,by = c("Protein","APRcount_tango"), all.x = TRUE)
master_df3 <- merge(master_df3,temp2,by = c("Protein"), all.x = TRUE)
master_df3 <- master_df3[order(master_df3$Protein, master_df3$Position),]
rownames(master_df3) <- NULL
master_df3$maxProtscore[which(is.na(master_df3$maxProtscore))] <- 0
master_df3$avgScore[which(is.na(master_df3$avgScore))] <- 0
save(master_df3, file = "Data/RData/TANGO_PTM_HUMAN_15/FILTERED/master_df.RData")
save(onlyapr_tango, file = "Data/RData/TANGO_PTM_HUMAN_15/FILTERED/onlyapr_tango.RData")
save(tango_APRs, file = "Data/RData/TANGO_PTM_HUMAN_15/FILTERED/tango_APRs_ALL.RData")
save(tango_APRs2, file = "Data/RData/TANGO_PTM_HUMAN_15/FILTERED/tango_APRs.RData")

