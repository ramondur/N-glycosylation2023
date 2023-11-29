#PTM dataset obtained from dbPTM. Accessed 24 of June of 2022. Files are stored in ./Data/dbPTM
#Parse files into a single dataset
Acetylation <- read.delim("Data/dbPTM/Acetylation", header=FALSE)
names <- t(as.data.frame(strsplit(as.character(Acetylation$V1),'_')))[,2]
residues <- strsplit(as.character(Acetylation$V6),'')
res <- NULL
res <- matrix(nrow = length(Acetylation$V6), ncol = 1)
for (i in 1:length(Acetylation$V6)){
  res[i,] <- residues[[i]][11]
}
Acetylation <- data.frame(Acetylation$V2,Acetylation$V3,Acetylation$V4)
colnames(Acetylation) <- c("Protein","Position","Type")
Acetylation$Residues <- as.factor(res) 
Acetylation$names <- names

NlinkedGlycosylation <- read.delim("Data/dbPTM/N-linked Glycosylation", header=FALSE)
names <- t(as.data.frame(strsplit(as.character(NlinkedGlycosylation$V1),'_')))[,2]
residues <- strsplit(as.character(NlinkedGlycosylation$V6),'')
res <- NULL
for (i in 1:length(NlinkedGlycosylation$V6)){
  res <- append(res,residues[[i]][11])
}
NlinkedGlycosylation <- data.frame(NlinkedGlycosylation$V2,NlinkedGlycosylation$V3,NlinkedGlycosylation$V4)
colnames(NlinkedGlycosylation) <- c("Protein","Position","Type")
NlinkedGlycosylation$Residues <- as.factor(res) 
NlinkedGlycosylation$names <- names

OlinkedGlycosylation <- read.delim("Data/dbPTM/O-linked Glycosylation", header=FALSE)
names <- t(as.data.frame(strsplit(as.character(OlinkedGlycosylation$V1),'_')))[,2]
residues <- strsplit(as.character(OlinkedGlycosylation$V6),'')
res <- NULL
for (i in 1:length(OlinkedGlycosylation$V6)){
  res <- append(res,residues[[i]][11])
}
OlinkedGlycosylation <- data.frame(OlinkedGlycosylation$V2,OlinkedGlycosylation$V3,OlinkedGlycosylation$V4)
colnames(OlinkedGlycosylation) <- c("Protein","Position","Type")
OlinkedGlycosylation$Residues <- as.factor(res) 
OlinkedGlycosylation$names <- names

Palmitoylation <- read.delim("Data/dbPTM/S-palmitoylation", header=FALSE)
names <- t(as.data.frame(strsplit(as.character(Palmitoylation$V1),'_')))[,2]
residues <- strsplit(as.character(Palmitoylation$V6),'')
res <- NULL
for (i in 1:length(Palmitoylation$V6)){
  res <- append(res,residues[[i]][11])
}
Palmitoylation <- data.frame(Palmitoylation$V2,Palmitoylation$V3,Palmitoylation$V4)
colnames(Palmitoylation) <- c("Protein","Position","Type")
Palmitoylation$Residues <- as.factor(res) 
Palmitoylation$names <- names

Phosphorylation <- read.delim("Data/dbPTM/Phosphorylation", header=FALSE)
list_names <- strsplit(as.character(Phosphorylation$V1), "_")
res <- matrix(nrow = length(Phosphorylation$V1), ncol = 1)
for (i in 1:length(Phosphorylation$V1)){
  res[i,] <- list_names[[i]][2]
}
names <- res
residues <- strsplit(as.character(Phosphorylation$V6),'')
res <- matrix(nrow = length(Phosphorylation$V6), ncol = 1)
for (i in 1:length(Phosphorylation$V6)){
  res[i,] <- residues[[i]][11]
}
Phosphorylation <- data.frame(Phosphorylation$V2,Phosphorylation$V3,Phosphorylation$V4)
colnames(Phosphorylation) <- c("Protein","Position","Type")
Phosphorylation$Residues <- as.factor(res) 
Phosphorylation$names <- names

Succinylation <- read.delim("Data/dbPTM/Succinylation", header=FALSE)
names <- t(as.data.frame(strsplit(as.character(Succinylation$V1),'_')))[,2]
residues <- strsplit(as.character(Succinylation$V6),'')
res <- NULL
for (i in 1:length(Succinylation$V6)){
  res <- append(res,residues[[i]][11])
}
Succinylation <- data.frame(Succinylation$V2,Succinylation$V3,Succinylation$V4)
colnames(Succinylation) <- c("Protein","Position","Type")
Succinylation$Residues <- as.factor(res) 
Succinylation$names <- names

Sulfoxidation <- read.delim("Data/dbPTM/Sulfoxidation", header=FALSE)
names <- t(as.data.frame(strsplit(as.character(Sulfoxidation$V1),'_')))[,2]
residues <- strsplit(as.character(Sulfoxidation$V6),'')
res <- NULL
for (i in 1:length(Sulfoxidation$V6)){
  res <- append(res,residues[[i]][11])
}
Sulfoxidation <- data.frame(Sulfoxidation$V2,Sulfoxidation$V3,Sulfoxidation$V4)
colnames(Sulfoxidation) <- c("Protein","Position","Type")
Sulfoxidation$Residues <- as.factor(res) 
Sulfoxidation$names <- names

Ubiquitination <- read.delim("Data/dbPTM/Ubiquitination", header=FALSE)
names <- t(as.data.frame(strsplit(as.character(Ubiquitination$V1),'_')))[,2]
residues <- strsplit(as.character(Ubiquitination$V6),'')
res <- matrix(nrow = length(Ubiquitination$V6), ncol = 1)
for (i in 1:length(Ubiquitination$V6)){
  res[i,] <- residues[[i]][11]
}
Ubiquitination <- data.frame(Ubiquitination$V2,Ubiquitination$V3,Ubiquitination$V4)
colnames(Ubiquitination) <- c("Protein","Position","Type")
Ubiquitination$Residues <- as.factor(res) 
Ubiquitination$names <- names

Methylation <- read.delim("Data/dbPTM/Methylation", header=FALSE)
names <- t(as.data.frame(strsplit(as.character(Methylation$V1),'_')))[,2]
residues <- strsplit(as.character(Methylation$V6),'')
res <- NULL
for (i in 1:length(Methylation$V6)){
  res <- append(res,residues[[i]][11])
}
Methylation <- data.frame(Methylation$V2,Methylation$V3,Methylation$V4)
colnames(Methylation) <- c("Protein","Position","Type")
Methylation$Residues <- as.factor(res) 
Methylation$names <- names

Malonylation <- read.delim("Data/dbPTM/Malonylation", header=FALSE)
names <- t(as.data.frame(strsplit(as.character(Malonylation$V1),'_')))[,2]
residues <- strsplit(as.character(Malonylation$V6),'')
res <- NULL
for (i in 1:length(Malonylation$V6)){
  res <- append(res,residues[[i]][11])
}
Malonylation <- data.frame(Malonylation$V2,Malonylation$V3,Malonylation$V4)
colnames(Malonylation) <- c("Protein","Position","Type")
Malonylation$Residues <- as.factor(res) 
Malonylation$names <- names

Glutathionylation <- read.delim("Data/dbPTM/Glutathionylation", header=FALSE)
names <- t(as.data.frame(strsplit(as.character(Glutathionylation$V1),'_')))[,2]
residues <- strsplit(as.character(Glutathionylation$V6),'')
res <- NULL
for (i in 1:length(Glutathionylation$V6)){
  res <- append(res,residues[[i]][11])
}
Glutathionylation <- data.frame(Glutathionylation$V2,Glutathionylation$V3,Glutathionylation$V4)
colnames(Glutathionylation) <- c("Protein","Position","Type")
Glutathionylation$Residues <- as.factor(res) 
Glutathionylation$names <- names

Sumoylation <- read.delim("Data/dbPTM/Sumoylation", header=FALSE)
names <- t(as.data.frame(strsplit(as.character(Sumoylation$V1),'_')))[,2]
residues <- strsplit(as.character(Sumoylation$V6),'')
res <- NULL
for (i in 1:length(Sumoylation$V6)){
  res <- append(res,residues[[i]][11])
}
Sumoylation <- data.frame(Sumoylation$V2,Sumoylation$V3,Sumoylation$V4)
colnames(Sumoylation) <- c("Protein","Position","Type")
Sumoylation$Residues <- as.factor(res) 
Sumoylation$names <- names

Snitrosylation <- read.delim("Data/dbPTM/S-nitrosylation", header=FALSE)
names <- t(as.data.frame(strsplit(as.character(Snitrosylation$V1),'_')))[,2]
residues <- strsplit(as.character(Snitrosylation$V6),'')
res <- NULL
for (i in 1:length(Snitrosylation$V6)){
  res <- append(res,residues[[i]][11])
}
Snitrosylation <- data.frame(Snitrosylation$V2,Snitrosylation$V3,Snitrosylation$V4)
colnames(Snitrosylation) <- c("Protein","Position","Type")
Snitrosylation$Residues <- as.factor(res) 
Snitrosylation$names <- names

Neddylation <- read.delim("Data/dbPTM/Neddylation", header=FALSE)
names <- t(as.data.frame(strsplit(as.character(Neddylation$V1),'_')))[,2]
residues <- strsplit(as.character(Neddylation$V6),'')
res <- NULL
for (i in 1:length(Neddylation$V6)){
  res <- append(res,residues[[i]][11])
}
Neddylation <- data.frame(Neddylation$V2,Neddylation$V3,Neddylation$V4)
colnames(Neddylation) <- c("Protein","Position","Type")
Neddylation$Residues <- as.factor(res) 
Neddylation$names <- names

PTMs_df <- rbind(Acetylation,NlinkedGlycosylation,OlinkedGlycosylation,Palmitoylation,Phosphorylation,Succinylation,Ubiquitination,
                 Methylation,Malonylation,Glutathionylation,Sumoylation,Snitrosylation,Neddylation,Sulfoxidation)
PTMs_df <- unique(PTMs_df)
PTMs_df$Type2 <- paste0(PTMs_df$Type,"_",PTMs_df$Residue)
PTMs_df_human <- subset(PTMs_df, names == "HUMAN", select = c(1,2,4,3,6))
PTMs_df_human_1000 <- names(which(table(PTMs_df_human$Type2) > 1000))
PTMs_df_human <- subset(PTMs_df_human, Type2 %in% PTMs_df_human_1000)
table(PTMs_df_human$Type2)

save(PTMs_df, file = "Data/dbPTM/all_PTM_data.RData")
save(PTMs_df_human, file = "Data/dbPTM/human_PTM_data.RData")