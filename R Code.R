library(dplyr)
library(tidyverse)
library(ggplot2)

# Make sure you do analysis in terms of frequencies because we send different number of cells for different samples, 
# so templates number is just relative
setwd("//172.25.60.69/CCTI_Labs/CCTI_SYKES/EXP FF7- XKLT-001 Exp/Bulk TCRb sequincing Adaptive Data XKLT001/All MLRs, Official Analysis")
# or, if working on Mac
setwd("/Volumes/CCTI_Labs/CCTI_SYKES/EXP FF7- XKLT-001 Exp/Bulk TCRb sequincing Adaptive Data XKLT001/All MLRs, Official Analysis")

CombinedRearrangements <- read_tsv("Raw MLR Data/CombinedRearrangements.tsv")
CombinedRearrangements = as.data.frame(CombinedRearrangements)
rownames(CombinedRearrangements)=CombinedRearrangements[,1]

#remove MART1 & other known lab contamination

contaminants.sequences <- c(
  "CASSAETQYF", #20D11
  "CASSIEGPTGELFF", #D222D2
  "CASSSFWGSDTGELFF", #clone 5
  "CASSLAGGANSPLHF", # #BRI164 R164 
  "CASSLSFGTEAFF",#DMF5 MART1:26-35
  "CASSSTGLPYGYTF",#HAI.7
  "CASSTRLAGDEQFF", #1.C8
  "CASSLSASGGATDTQYF",#A3.10
  "CASWTVSYNEQFF", #1F3 
  "CASSLGPGQRETQYF", #A5.5
  "CASSLERDGYTF" #A1.9
)

CombinedRearrangements = CombinedRearrangements[setdiff(rownames(CombinedRearrangements), contaminants.sequences),]

# resolve ambiguous, ie. remove putative clones between CD4 & CD8 that are present in one condition 
# in a ratio less than 10 compared to the other
# need to open alloreactivity_pintoprofile R script in WD and save functions

CombinedRearrangements$CD4 = CombinedRearrangements[,4]/7 + CombinedRearrangements[,6]/7 + 
  CombinedRearrangements[,8]/7 + CombinedRearrangements[,10]/7 + CombinedRearrangements[,12]/7 + 
  CombinedRearrangements[,14]/7 + CombinedRearrangements[,16]/7 
CombinedRearrangements$CD8 = CombinedRearrangements[,5]/7 + CombinedRearrangements[,7]/7 + 
  CombinedRearrangements[,9]/7 + CombinedRearrangements[,11]/7 + CombinedRearrangements[,13]/7 + 
  CombinedRearrangements[,15]/7 + CombinedRearrangements[,17]/7  
CombinedRearrangements=resolveambiguous(CombinedRearrangements, 19, 20)
CD4=which(CombinedRearrangements[,19]>0)
CD8=which(CombinedRearrangements[,20]>0)

#removing ambig from all other CD4/CD8 columns
CombinedRearrangements[setdiff(1:nrow(CombinedRearrangements),CD4),c(4,6,8,10,12,14,16)]=0
CombinedRearrangements[setdiff(1:nrow(CombinedRearrangements),CD8),c(5,7,9,11,13,15,17)]=0

# Update Present In

CombinedRearrangements$`Present In` <- rowSums(CombinedRearrangements[, 4:18] > 0)

# save file 

write.csv(CombinedRearrangements, file="CombinedRearrangementsAllMLRs.csv")

# Identify DRTCC's

# Start with PreTx Direct LN MLR 
cd4.HVG=CombinedRearrangements[,c(18, 8)]
cd8.HVG=CombinedRearrangements[,c(18,9)]
cd4.HVG = cd4.HVG[cd4.HVG[,2] > 0.00001, ]
cd8.HVG = cd8.HVG[cd8.HVG[,2] > 0.00001, ]
cd4.HVG = cd4.HVG[cd4.HVG[,2] >= cd4.HVG[,1]*2, ]
cd8.HVG = cd8.HVG[cd8.HVG[,2] >= cd8.HVG[,1]*2, ]
cd4.HVG.PreTxDirectLN = cd4.HVG
cd8.HVG.PreTxDirectLN = cd8.HVG

# PreTx Direct PBMC MLR 
cd4.HVG=CombinedRearrangements[,c(18, 10)]
cd8.HVG=CombinedRearrangements[,c(18,11)]
cd4.HVG = cd4.HVG[cd4.HVG[,2] > 0.00001, ]
cd8.HVG = cd8.HVG[cd8.HVG[,2] > 0.00001, ]
cd4.HVG = cd4.HVG[cd4.HVG[,2] >= cd4.HVG[,1]*2, ]
cd8.HVG = cd8.HVG[cd8.HVG[,2] >= cd8.HVG[,1]*2, ]
cd4.HVG.PreTxDirectPBMC = cd4.HVG
cd8.HVG.PreTxDirectPBMC = cd8.HVG

# PreTx Indirect

cd4.HVG=CombinedRearrangements[,c(18, 12)]
cd8.HVG=CombinedRearrangements[,c(18,13)]
cd4.HVG = cd4.HVG[cd4.HVG[,2] > 0.00001, ]
cd8.HVG = cd8.HVG[cd8.HVG[,2] > 0.00001, ]
cd4.HVG = cd4.HVG[cd4.HVG[,2] >= cd4.HVG[,1]*2, ]
cd8.HVG = cd8.HVG[cd8.HVG[,2] >= cd8.HVG[,1]*2, ]
cd4.HVG.PreTxIndirectPBMC = cd4.HVG
cd8.HVG.PreTxIndirectPBMC = cd8.HVG

# POD28 MLR

cd4.HVG=CombinedRearrangements[,c(4, 16)]
cd8.HVG=CombinedRearrangements[,c(5,17)]
cd4.HVG = cd4.HVG[cd4.HVG[,2] > 0.00001, ]
cd8.HVG = cd8.HVG[cd8.HVG[,2] > 0.00001, ]
cd4.HVG = cd4.HVG[cd4.HVG[,2] >= cd4.HVG[,1]*2, ]
cd8.HVG = cd8.HVG[cd8.HVG[,2] >= cd8.HVG[,1]*2, ]
cd4.HVG.POD28 = cd4.HVG
cd8.HVG.POD28 = cd8.HVG

# POD49 MLR

cd4.HVG=CombinedRearrangements[,c(14, 6)]
cd8.HVG=CombinedRearrangements[,c(15,7)]
cd4.HVG = cd4.HVG[cd4.HVG[,2] > 0.00001, ]
cd8.HVG = cd8.HVG[cd8.HVG[,2] > 0.00001, ]
cd4.HVG = cd4.HVG[cd4.HVG[,2] >= cd4.HVG[,1]*2, ]
cd8.HVG = cd8.HVG[cd8.HVG[,2] >= cd8.HVG[,1]*2, ]
cd4.HVG.POD49 = cd4.HVG
cd8.HVG.POD49 = cd8.HVG

# Create string of unique CD4 DRTCC sequences 

cd4.HVG = list(cd4.HVG.PreTxDirectLN, cd4.HVG.PreTxDirectPBMC, cd4.HVG.PreTxIndirectPBMC, cd4.HVG.POD28, cd4.HVG.POD49)
row_names_list <- lapply(cd4.HVG, rownames)
cd4.DRTCCs <- unlist(row_names_list)
cd4.DRTCCs = unique(cd4.DRTCCs)

# Create string of unique CD8 DRTCC sequences 

cd8.HVG = list(cd8.HVG.PreTxDirectLN, cd8.HVG.PreTxDirectPBMC, cd8.HVG.PreTxIndirectPBMC, cd8.HVG.POD28, cd8.HVG.POD49)
row_names_list <- lapply(cd8.HVG, rownames)
cd8.DRTCCs <- unlist(row_names_list)
cd8.DRTCCs = unique(cd8.DRTCCs)

# Both CD4 & CD8

DRTCCs = c(cd4.DRTCCs, cd8.DRTCCs)
DRTCCs = unique(DRTCCs)

# Now create master table with all of the DRTCCs and their frequencies in the MLRs

MasterDRTCC = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% DRTCCs, c(1,3,8,9,10,11,12,13,16,17,6,7)]

# Update Present In

MasterDRTCC$`Present In` <- rowSums(MasterDRTCC[, 3:12] > 0)

write.csv(MasterDRTCC, file="MasterDRTCC.csv")

##################### Now track DRTCCs in Unstim populations

CombinedRearrangements <- read_tsv("Unstim Raw Data/CombinedRearrangements.tsv")
CombinedRearrangements = as.data.frame(CombinedRearrangements)
rownames(CombinedRearrangements)=CombinedRearrangements[,1]

#remove MART1 & other known lab contamination

contaminants.sequences <- c(
  "CASSAETQYF", #20D11
  "CASSIEGPTGELFF", #D222D2
  "CASSSFWGSDTGELFF", #clone 5
  "CASSLAGGANSPLHF", # #BRI164 R164 
  "CASSLSFGTEAFF",#DMF5 MART1:26-35
  "CASSSTGLPYGYTF",#HAI.7
  "CASSTRLAGDEQFF", #1.C8
  "CASSLSASGGATDTQYF",#A3.10
  "CASWTVSYNEQFF", #1F3 
  "CASSLGPGQRETQYF", #A5.5
  "CASSLERDGYTF" #A1.9
)

CombinedRearrangements = CombinedRearrangements[setdiff(rownames(CombinedRearrangements), contaminants.sequences),]

# resolve ambiguous, ie. remove putative clones between CD4 & CD8 that are present in one condition 
# in a ratio less than 10 compared to the other
# need to open alloreactivity_pintoprofile R script in WD and save functions

CombinedRearrangements$CD4 = CombinedRearrangements[,4]/3 + CombinedRearrangements[,7]/3 + 
  CombinedRearrangements[,9]/3 
CombinedRearrangements$CD8 = CombinedRearrangements[,5]/3 + CombinedRearrangements[,8]/3 + 
  CombinedRearrangements[,10]/3
CombinedRearrangements=resolveambiguous(CombinedRearrangements, 12, 13)
CD4=which(CombinedRearrangements[,12]>0)
CD8=which(CombinedRearrangements[,13]>0)

#removing ambig from all other CD4/CD8 columns
CombinedRearrangements[setdiff(1:nrow(CombinedRearrangements),CD4),c(4,7,9)]=0
CombinedRearrangements[setdiff(1:nrow(CombinedRearrangements),CD8),c(5,8,10)]=0

# Update Present In

CombinedRearrangements$`Present In` <- rowSums(CombinedRearrangements[, 4:11] > 0)

# save file 

write.csv(CombinedRearrangements, file="CombinedRearrangementsAllUnstim.csv")

# Create Master CD4 DRTCC table

MasterCD4DRTCC = MasterDRTCC[,c(2,3,4,6,8,10,12)]
# Remove rows where at least one of the specified columns is 0
MasterCD4DRTCC <- MasterCD4DRTCC[rowSums(MasterCD4DRTCC[,c(3:7)] == 0) < 5, ]
write.csv(MasterCD4DRTCC, file="MasterCD4DRTCC.csv")

# Create Master CD8 DRTCC table

MasterCD8DRTCC = MasterDRTCC[,c(2,3,5,7,9,11,13)]
# Remove rows where at least one of the specified columns is 0
MasterCD8DRTCC <- MasterCD8DRTCC[rowSums(MasterCD8DRTCC[,c(3:7)] == 0) < 5, ]
write.csv(MasterCD8DRTCC, file="MasterCD8DRTCC.csv")

# Create Master NonDRTCC table

MasterNonDRTCC = CombinedRearrangements[!(CombinedRearrangements$`Amino Acid` %in% MasterDRTCC$`Amino Acid`), c(1,11)]
MasterNonDRTCC = MasterNonDRTCC[MasterNonDRTCC$XKLT001_Unstim_PBMC_TCRB>0,]
write.csv(MasterNonDRTCC, file="MasterNonDRTCC.csv")
sum(MasterNonDRTCC$XKLT001_Unstim_PBMC_TCRB)*19783

# Now track DRTCCs 
setwd("//172.25.60.69/CCTI_Labs/CCTI_SYKES/EXP FF7- XKLT-001 Exp/Bulk TCRb sequincing Adaptive Data XKLT001/All MLRs, Official Analysis/Results")
# Start with PreTx Unstim PBMCs (Unsorted)

# Unique Clones in PreTx Unstim PBMC
nrow(CombinedRearrangements[CombinedRearrangements$XKLT001_Unstim_PBMC_TCRB>0,])
# Total Freq of Unique Clones in PreTx Unstim PBMC -- then multiply by total productive templates 
# (found on Adaptive website) in order to get the number of templates
sum(CombinedRearrangements[CombinedRearrangements$XKLT001_Unstim_PBMC_TCRB>0,11])*19783

CD4DRTCCsinPreTxUnstimPBMC = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterCD4DRTCC$`Amino Acid`, c(1,11)]
CD4DRTCCsinPreTxUnstimPBMC = CD4DRTCCsinPreTxUnstimPBMC[CD4DRTCCsinPreTxUnstimPBMC$XKLT001_Unstim_PBMC_TCRB>0,]
sum(CD4DRTCCsinPreTxUnstimPBMC$XKLT001_Unstim_PBMC_TCRB)*19783
write.csv(CD4DRTCCsinPreTxUnstimPBMC, file="CD4DRTCCsinPreTxUnstimPBMC.csv")

CD8DRTCCsinPreTxUnstimPBMC = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterCD8DRTCC$`Amino Acid`, c(1,11)]
CD8DRTCCsinPreTxUnstimPBMC = CD8DRTCCsinPreTxUnstimPBMC[CD8DRTCCsinPreTxUnstimPBMC$XKLT001_Unstim_PBMC_TCRB>0,]
sum(CD8DRTCCsinPreTxUnstimPBMC$XKLT001_Unstim_PBMC_TCRB)*19783
write.csv(CD8DRTCCsinPreTxUnstimPBMC, file="CD8DRTCCsinPreTxUnstimPBMC.csv")

TotalDRTCCsinPreTxUnstimPBMC = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterDRTCC$`Amino Acid`, c(1,11)]
TotalDRTCCsinPreTxUnstimPBMC = TotalDRTCCsinPreTxUnstimPBMC[TotalDRTCCsinPreTxUnstimPBMC$XKLT001_Unstim_PBMC_TCRB>0,]
sum(TotalDRTCCsinPreTxUnstimPBMC$XKLT001_Unstim_PBMC_TCRB)*19783
sum(TotalDRTCCsinPreTxUnstimPBMC$XKLT001_Unstim_PBMC_TCRB)*100
write.csv(TotalDRTCCsinPreTxUnstimPBMC, file="TotalDRTCCsinPreTxUnstimPBMC.csv")

NonDRTCCsinPreTxUnstimPBMC = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterNonDRTCC$`Amino Acid`, c(1,11)]
NonDRTCCsinPreTxUnstimPBMC = NonDRTCCsinPreTxUnstimPBMC[NonDRTCCsinPreTxUnstimPBMC$XKLT001_Unstim_PBMC_TCRB>0,]
sum(NonDRTCCsinPreTxUnstimPBMC$XKLT001_Unstim_PBMC_TCRB)*19783
sum(NonDRTCCsinPreTxUnstimPBMC$XKLT001_Unstim_PBMC_TCRB)*100
write.csv(NonDRTCCsinPreTxUnstimPBMC, file="NonDRTCCsinPreTxUnstimPBMC.csv")

# Now with POD28 Unstim
# CD4
# Unique Clones in POD28 PBMC CD4
nrow(CombinedRearrangements[CombinedRearrangements$`XKLT001-POD28-Unstim-CD4_TCRB`>0,])
# Total Freq of Unique Clones in PreTx Unstim PBMC -- then multiply by total productive templates 
# (found on Adaptive website) in order to get the number of templates
sum(CombinedRearrangements[CombinedRearrangements$`XKLT001-POD28-Unstim-CD4_TCRB`>0,4])*29551

CD4DRTCCsinPOD28CD4 = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterCD4DRTCC$`Amino Acid`, c(1,4)]
CD4DRTCCsinPOD28CD4 = CD4DRTCCsinPOD28CD4[CD4DRTCCsinPOD28CD4$`XKLT001-POD28-Unstim-CD4_TCRB`>0,]
sum(CD4DRTCCsinPOD28CD4$`XKLT001-POD28-Unstim-CD4_TCRB`)*29551
sum(CD4DRTCCsinPOD28CD4$`XKLT001-POD28-Unstim-CD4_TCRB`)*100
write.csv(CD4DRTCCsinPOD28CD4, file="CD4DRTCCsinPOD28CD4.csv")

NonDRTCCsinPOD28CD4 = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterNonDRTCC$`Amino Acid`, c(1,4)]
NonDRTCCsinPOD28CD4 = NonDRTCCsinPOD28CD4[NonDRTCCsinPOD28CD4$`XKLT001-POD28-Unstim-CD4_TCRB`>0,]
sum(NonDRTCCsinPOD28CD4$`XKLT001-POD28-Unstim-CD4_TCRB`)*29551
sum(NonDRTCCsinPOD28CD4$`XKLT001-POD28-Unstim-CD4_TCRB`)*100
write.csv(NonDRTCCsinPOD28CD4, file="NonDRTCCsinPOD28CD4.csv")

# Relative Clonal Detection Rate
print((nrow(CD4DRTCCsinPOD28CD4)/48)/(nrow(NonDRTCCsinPOD28CD4)/18183))
# Relative Template Detection Rate -- or Normalized Cumulative Frequency Fold Change
print(((sum(CD4DRTCCsinPOD28CD4$`XKLT001-POD28-Unstim-CD4_TCRB`)*29551)/70)/((sum(NonDRTCCsinPOD28CD4$`XKLT001-POD28-Unstim-CD4_TCRB`)*29551)/19658))

# CD8

# Unique Clones in POD28 PBMC CD8
nrow(CombinedRearrangements[CombinedRearrangements$`XKLT001-POD28-Unstim-CD8_TCRB`>0,])
# Total Freq of Unique Clones in PreTx Unstim PBMC -- then multiply by total productive templates 
# (found on Adaptive website) in order to get the number of templates
sum(CombinedRearrangements[CombinedRearrangements$`XKLT001-POD28-Unstim-CD4_TCRB`>0,4])*41186

CD8DRTCCsinPOD28CD8 = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterCD8DRTCC$`Amino Acid`, c(1,5)]
CD8DRTCCsinPOD28CD8 = CD8DRTCCsinPOD28CD8[CD8DRTCCsinPOD28CD8$`XKLT001-POD28-Unstim-CD8_TCRB`>0,]
sum(CD8DRTCCsinPOD28CD8$`XKLT001-POD28-Unstim-CD8_TCRB`)*41186
sum(CD8DRTCCsinPOD28CD8$`XKLT001-POD28-Unstim-CD8_TCRB`)*100
write.csv(CD8DRTCCsinPOD28CD8, file="CD8DRTCCsinPOD28CD8.csv")

NonDRTCCsinPOD28CD8 = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterNonDRTCC$`Amino Acid`, c(1,5)]
NonDRTCCsinPOD28CD8 = NonDRTCCsinPOD28CD8[NonDRTCCsinPOD28CD8$`XKLT001-POD28-Unstim-CD8_TCRB`>0,]
sum(NonDRTCCsinPOD28CD8$`XKLT001-POD28-Unstim-CD8_TCRB`)*41186
sum(NonDRTCCsinPOD28CD8$`XKLT001-POD28-Unstim-CD8_TCRB`)*100
write.csv(NonDRTCCsinPOD28CD8, file="NonDRTCCsinPOD28CD8.csv")

# Relative Clonal Detection Rate
print((nrow(CD8DRTCCsinPOD28CD8)/23)/(nrow(NonDRTCCsinPOD28CD8)/18183))
# Relative Template Detection Rate -- or Normalized Cumulative Frequency Fold Change
print(((sum(CD8DRTCCsinPOD28CD8$`XKLT001-POD28-Unstim-CD8_TCRB`)*41186)/55)/((sum(NonDRTCCsinPOD28CD8$`XKLT001-POD28-Unstim-CD8_TCRB`)*41186)/19658))

# Now, POD33 PBMC
# CD4

# Unique Clones in POD33 PBMC CD4
nrow(CombinedRearrangements[CombinedRearrangements$XKLT001_POD33_CD4_Unstim_TCRB>0,])
# Total Freq of Unique Clones in PreTx Unstim PBMC -- then multiply by total productive templates 
# (found on Adaptive website) in order to get the number of templates
sum(CombinedRearrangements[CombinedRearrangements$XKLT001_POD33_CD4_Unstim_TCRB>0,7])*8248

CD4DRTCCsinPOD33CD4 = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterCD4DRTCC$`Amino Acid`, c(1,7)]
CD4DRTCCsinPOD33CD4 = CD4DRTCCsinPOD33CD4[CD4DRTCCsinPOD33CD4$XKLT001_POD33_CD4_Unstim_TCRB>0,]
nrow(CD4DRTCCsinPOD33CD4)
sum(CD4DRTCCsinPOD33CD4$XKLT001_POD33_CD4_Unstim_TCRB)*8248
sum(CD4DRTCCsinPOD33CD4$XKLT001_POD33_CD4_Unstim_TCRB)*100
write.csv(CD4DRTCCsinPOD33CD4, file="CD4DRTCCsinPOD33CD4.csv")

NonDRTCCsinPOD33CD4 = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterNonDRTCC$`Amino Acid`, c(1,7)]
NonDRTCCsinPOD33CD4 = NonDRTCCsinPOD33CD4[NonDRTCCsinPOD33CD4$XKLT001_POD33_CD4_Unstim_TCRB>0,]
nrow(NonDRTCCsinPOD33CD4)
sum(NonDRTCCsinPOD33CD4$XKLT001_POD33_CD4_Unstim_TCRB)*8248
sum(NonDRTCCsinPOD33CD4$XKLT001_POD33_CD4_Unstim_TCRB)*100
write.csv(NonDRTCCsinPOD33CD4, file="NonDRTCCsinPOD33CD4.csv")

# Relative Clonal Detection Rate
print((nrow(CD4DRTCCsinPOD33CD4)/48)/(nrow(NonDRTCCsinPOD33CD4)/18183))
# Relative Template Detection Rate -- or Normalized Cumulative Frequency Fold Change
print(((sum(CD4DRTCCsinPOD33CD4$XKLT001_POD33_CD4_Unstim_TCRB)*8248)/70)/((sum(NonDRTCCsinPOD33CD4$XKLT001_POD33_CD4_Unstim_TCRB)*8248)/19658))

# CD8

# Unique Clones in POD33 PBMC CD8
nrow(CombinedRearrangements[CombinedRearrangements$XKLT001_POD33_CD8_Unstim_TCRB>0,])
# Total Freq of Unique Clones in PreTx Unstim PBMC -- then multiply by total productive templates 
# (found on Adaptive website) in order to get the number of templates
sum(CombinedRearrangements[CombinedRearrangements$XKLT001_POD33_CD8_Unstim_TCRB>0,8])*12942

CD8DRTCCsinPOD33CD8 = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterCD8DRTCC$`Amino Acid`, c(1,8)]
CD8DRTCCsinPOD33CD8 = CD8DRTCCsinPOD33CD8[CD8DRTCCsinPOD33CD8$XKLT001_POD33_CD8_Unstim_TCRB>0,]
nrow(CD8DRTCCsinPOD33CD8)
sum(CD8DRTCCsinPOD33CD8$XKLT001_POD33_CD8_Unstim_TCRB)*12942
sum(CD8DRTCCsinPOD33CD8$XKLT001_POD33_CD8_Unstim_TCRB)*100
write.csv(CD8DRTCCsinPOD33CD8, file="CD8DRTCCsinPOD33CD8.csv")

NonDRTCCsinPOD33CD8 = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterNonDRTCC$`Amino Acid`, c(1,8)]
NonDRTCCsinPOD33CD8 = NonDRTCCsinPOD33CD8[NonDRTCCsinPOD33CD8$XKLT001_POD33_CD8_Unstim_TCRB>0,]
nrow(NonDRTCCsinPOD33CD8)
sum(NonDRTCCsinPOD33CD8$XKLT001_POD33_CD8_Unstim_TCRB)*12942
sum(NonDRTCCsinPOD33CD8$XKLT001_POD33_CD8_Unstim_TCRB)*100
write.csv(NonDRTCCsinPOD33CD8, file="NonDRTCCsinPOD33CD8.csv")

# Relative Clonal Detection Rate
print((nrow(CD8DRTCCsinPOD33CD8)/23)/(nrow(NonDRTCCsinPOD33CD8)/18183))
# Relative Template Detection Rate -- or Normalized Cumulative Frequency Fold Change
print(((sum(CD8DRTCCsinPOD33CD8$XKLT001_POD33_CD8_Unstim_TCRB)*12942)/55)/((sum(NonDRTCCsinPOD33CD8$XKLT001_POD33_CD8_Unstim_TCRB)*12942)/19658))

# POD49

# CD4

# Unique Clones in POD49 PBMC CD4
nrow(CombinedRearrangements[CombinedRearrangements$XKLT001_POD49_CD4_Unstim_TCRB>0,])
# Total Freq of Unique Clones in PreTx Unstim PBMC -- then multiply by total productive templates 
# (found on Adaptive website) in order to get the number of templates
sum(CombinedRearrangements[CombinedRearrangements$XKLT001_POD49_CD4_Unstim_TCRB>0,9])*24722

CD4DRTCCsinPOD49CD4 = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterCD4DRTCC$`Amino Acid`, c(1,9)]
CD4DRTCCsinPOD49CD4 = CD4DRTCCsinPOD49CD4[CD4DRTCCsinPOD49CD4$XKLT001_POD49_CD4_Unstim_TCRB>0,]
nrow(CD4DRTCCsinPOD49CD4)
sum(CD4DRTCCsinPOD49CD4$XKLT001_POD49_CD4_Unstim_TCRB)*24722
sum(CD4DRTCCsinPOD49CD4$XKLT001_POD49_CD4_Unstim_TCRB)*100
write.csv(CD4DRTCCsinPOD49CD4, file="CD4DRTCCsinPOD49CD4.csv")

NonDRTCCsinPOD49CD4 = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterNonDRTCC$`Amino Acid`, c(1,9)]
NonDRTCCsinPOD49CD4 = NonDRTCCsinPOD49CD4[NonDRTCCsinPOD49CD4$XKLT001_POD49_CD4_Unstim_TCRB>0,]
nrow(NonDRTCCsinPOD49CD4)
sum(NonDRTCCsinPOD49CD4$XKLT001_POD49_CD4_Unstim_TCRB)*24722
sum(NonDRTCCsinPOD49CD4$XKLT001_POD49_CD4_Unstim_TCRB)*100
write.csv(NonDRTCCsinPOD49CD4, file="NonDRTCCsinPOD49CD4.csv")

# Relative Clonal Detection Rate
print((nrow(CD4DRTCCsinPOD49CD4)/48)/(nrow(NonDRTCCsinPOD49CD4)/18183))
# Relative Template Detection Rate -- or Normalized Cumulative Frequency Fold Change
print(((sum(CD4DRTCCsinPOD49CD4$XKLT001_POD49_CD4_Unstim_TCRB)*24722)/70)/((sum(NonDRTCCsinPOD49CD4$XKLT001_POD49_CD4_Unstim_TCRB)*24722)/19658))

# CD8

# Unique Clones in POD49 PBMC CD8
nrow(CombinedRearrangements[CombinedRearrangements$XKLT001_POD49_CD8_Unstim_TCRB>0,])
# Total Freq of Unique Clones in PreTx Unstim PBMC -- then multiply by total productive templates 
# (found on Adaptive website) in order to get the number of templates
sum(CombinedRearrangements[CombinedRearrangements$XKLT001_POD49_CD8_Unstim_TCRB>0,10])*60311

CD8DRTCCsinPOD49CD8 = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterCD8DRTCC$`Amino Acid`, c(1,10)]
CD8DRTCCsinPOD49CD8 = CD8DRTCCsinPOD49CD8[CD8DRTCCsinPOD49CD8$XKLT001_POD49_CD8_Unstim_TCRB>0,]
nrow(CD8DRTCCsinPOD49CD8)
sum(CD8DRTCCsinPOD49CD8$XKLT001_POD49_CD8_Unstim_TCRB)*60311
sum(CD8DRTCCsinPOD49CD8$XKLT001_POD49_CD8_Unstim_TCRB)*100
write.csv(CD8DRTCCsinPOD49CD8, file="CD8DRTCCsinPOD49CD8.csv")

NonDRTCCsinPOD49CD8 = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterNonDRTCC$`Amino Acid`, c(1,10)]
NonDRTCCsinPOD49CD8 = NonDRTCCsinPOD49CD8[NonDRTCCsinPOD49CD8$XKLT001_POD49_CD8_Unstim_TCRB>0,]
nrow(NonDRTCCsinPOD49CD8)
sum(NonDRTCCsinPOD49CD8$XKLT001_POD49_CD8_Unstim_TCRB)*60311
sum(NonDRTCCsinPOD49CD8$XKLT001_POD49_CD8_Unstim_TCRB)*100
write.csv(NonDRTCCsinPOD49CD8, file="NonDRTCCsinPOD49CD8.csv")

# Relative Clonal Detection Rate
print((nrow(CD8DRTCCsinPOD49CD8)/23)/(nrow(NonDRTCCsinPOD49CD8)/18183))
# Relative Template Detection Rate -- or Normalized Cumulative Frequency Fold Change
print(((sum(CD8DRTCCsinPOD49CD8$XKLT001_POD49_CD8_Unstim_TCRB)*60311)/55)/((sum(NonDRTCCsinPOD49CD8$XKLT001_POD49_CD8_Unstim_TCRB)*60311)/19658))

# Now, with POD49 Kidney Bx

# Unique Clones in POD49 Kidney Bx
nrow(CombinedRearrangements[CombinedRearrangements$`XKLT001-POD49-Kidney-biopsy_TCRB`>0,])
# Total Freq of Unique Clones in PreTx Unstim PBMC -- then multiply by total productive templates 
# (found on Adaptive website) in order to get the number of templates
sum(CombinedRearrangements[CombinedRearrangements$`XKLT001-POD49-Kidney-biopsy_TCRB`>0,6])*66

CD4DRTCCsinPOD49KidneyBx = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterCD4DRTCC$`Amino Acid`, c(1,6)]
CD4DRTCCsinPOD49KidneyBx = CD4DRTCCsinPOD49KidneyBx[CD4DRTCCsinPOD49KidneyBx$`XKLT001-POD49-Kidney-biopsy_TCRB`>0,]
sum(CD4DRTCCsinPOD49KidneyBx$`XKLT001-POD49-Kidney-biopsy_TCRB`)*66
write.csv(CD4DRTCCsinPOD49KidneyBx, file="CD4DRTCCsinPOD49KidneyBx.csv")

CD8DRTCCsinPOD49KidneyBx = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterCD8DRTCC$`Amino Acid`, c(1,6)]
CD8DRTCCsinPOD49KidneyBx = CD8DRTCCsinPOD49KidneyBx[CD8DRTCCsinPOD49KidneyBx$`XKLT001-POD49-Kidney-biopsy_TCRB`>0,]
sum(CD8DRTCCsinPOD49KidneyBx$`XKLT001-POD49-Kidney-biopsy_TCRB`)*66
write.csv(CD8DRTCCsinPOD49KidneyBx, file="CD8DRTCCsinPOD49KidneyBx.csv")

NonDRTCCsinPOD49KidneyBx = CombinedRearrangements[CombinedRearrangements$`Amino Acid` %in% MasterNonDRTCC$`Amino Acid`, c(1,6)]
NonDRTCCsinPOD49KidneyBx = NonDRTCCsinPOD49KidneyBx[NonDRTCCsinPOD49KidneyBx$`XKLT001-POD49-Kidney-biopsy_TCRB`>0,]
sum(NonDRTCCsinPOD49KidneyBx$`XKLT001-POD49-Kidney-biopsy_TCRB`)*66
sum(NonDRTCCsinPOD49KidneyBx$`XKLT001-POD49-Kidney-biopsy_TCRB`)*100
write.csv(NonDRTCCsinPOD49KidneyBx, file="NonDRTCCsinPOD49KidneyBx.csv")





























###############################

CD8DRTCCMo8 = CombinedRearrangementsPt1Freq[, c(1,5,7,13)]
rownames(CD8DRTCCMo8) = CD8DRTCCMo8[,1]
CD8DRTCCMo8 = CD8DRTCCMo8[CD8DRTCCMo8[,4] > 0,]
CD8DRTCCMo8 = CD8DRTCCMo8[CD8DRTCCFrequency[,1],]
CD8DRTCCMo8 = na.omit(CD8DRTCCMo8)
write.csv(CD8DRTCCMo8, file="CD8DRTCCMo8.csv")
sum(CD8DRTCCMo8[,3])

sum(CD8DRTCCFrequency[,9])

# Using templates

CombinedRearrangementsTemplates <- read_csv("CombinedRearrangementsTemplates.csv")
CombinedRearrangementsTemplates = as.data.frame(CombinedRearrangementsTemplates)
rownames(CombinedRearrangementsTemplates) = CombinedRearrangementsTemplates[,1]
CD4DRTCCTemplates = CombinedRearrangementsTemplates[rownames(cd4.HVG),]
CD8DRTCCTemplates = CombinedRearrangementsTemplates[rownames(cd8.HVG),]

write.csv(CD4DRTCCTemplates, file="CD4DRTCCTemplates.csv")
write.csv(CD8DRTCCTemplates, file="CD8DRTCCTemplates.csv")

sum(CD8DRTCCTemplates[,10])


#################### REAL CODE STARTS HERE

write.csv(POD33CD4, file="POD33CD4.csv")
write.csv(POD33CD8, file="POD33CD8.csv")
write.csv(POD49CD4, file="POD49CD4.csv")
write.csv(POD49CD8, file="POD49CD8.csv")

DRTCCinPOD33CD4 = POD33CD4[POD33CD4$`Amino Acid` %in% All_DRTCCs$CDR3.aa,]
sum(DRTCCinPOD33CD4[,4])
nonDRTCCinPOD33CD4 = POD33CD4[POD33CD4$`Amino Acid` %in% nonDRTCC$CDR3.aa,]
sum(nonDRTCCinPOD33CD4[,4])
DRTCCinPOD33CD4_Parent = All_DRTCCs[All_DRTCCs$CDR3.aa %in% DRTCCinPOD33CD4$`Amino Acid`,]
sum(DRTCCinPOD33CD4_Parent$Clones)
nonDRTCCinPOD33CD4_Parent = All_Unstim[All_Unstim$CDR3.aa %in% nonDRTCCinPOD33CD4$`Amino Acid`,]
sum(nonDRTCCinPOD33CD4_Parent$Clones)

write.csv(DRTCCinPOD33CD4, file="DRTCCinPOD33CD4.csv")
write.csv(nonDRTCCinPOD33CD4, file="nonDRTCCinPOD33CD4.csv")
write.csv(DRTCCinPOD33CD4_Parent, file="DRTCCinPOD33CD4_Parent.csv")
write.csv(nonDRTCCinPOD33CD4_Parent, file="nonDRTCCinPOD33CD4_Parent.csv")

DRTCCinPOD33CD8 = POD33CD8[POD33CD8$`Amino Acid` %in% All_DRTCCs$CDR3.aa,]
sum(DRTCCinPOD33CD8[,4])
nonDRTCCinPOD33CD8 = POD33CD8[POD33CD8$`Amino Acid` %in% nonDRTCC$CDR3.aa,]
sum(nonDRTCCinPOD33CD8[,4])
DRTCCinPOD33CD8_Parent = All_DRTCCs[All_DRTCCs$CDR3.aa %in% DRTCCinPOD33CD8$`Amino Acid`,]
sum(DRTCCinPOD33CD8_Parent$Clones)
nonDRTCCinPOD33CD8_Parent = All_Unstim[All_Unstim$CDR3.aa %in% nonDRTCCinPOD33CD8$`Amino Acid`,]
sum(nonDRTCCinPOD33CD8_Parent$Clones)

write.csv(DRTCCinPOD33CD8, file="DRTCCinPOD33CD8.csv")
write.csv(nonDRTCCinPOD33CD8, file="nonDRTCCinPOD33CD8.csv")
write.csv(DRTCCinPOD33CD8_Parent, file="DRTCCinPOD33CD8_Parent.csv")
write.csv(nonDRTCCinPOD33CD8_Parent, file="nonDRTCCinPOD33CD8_Parent.csv")

DRTCCinPOD49CD4 = POD49CD4[POD49CD4$`Amino Acid` %in% All_DRTCCs$CDR3.aa,]
sum(DRTCCinPOD49CD4[,4])
nonDRTCCinPOD49CD4 = POD49CD4[POD49CD4$`Amino Acid` %in% nonDRTCC$CDR3.aa,]
sum(nonDRTCCinPOD49CD4[,4])
DRTCCinPOD49CD4_Parent = All_DRTCCs[All_DRTCCs$CDR3.aa %in% DRTCCinPOD49CD4$`Amino Acid`,]
sum(DRTCCinPOD49CD4_Parent$Clones)
nonDRTCCinPOD49CD4_Parent = All_Unstim[All_Unstim$CDR3.aa %in% nonDRTCCinPOD49CD4$`Amino Acid`,]
sum(nonDRTCCinPOD49CD4_Parent$Clones)

write.csv(DRTCCinPOD49CD4, file="DRTCCinPOD49CD4.csv")
write.csv(nonDRTCCinPOD49CD4, file="nonDRTCCinPOD49CD4.csv")
write.csv(DRTCCinPOD49CD4_Parent, file="DRTCCinPOD49CD4_Parent.csv")
write.csv(nonDRTCCinPOD49CD4_Parent, file="nonDRTCCinPOD49CD4_Parent.csv")

DRTCCinPOD49CD8 = POD49CD8[POD49CD8$`Amino Acid` %in% All_DRTCCs$CDR3.aa,]
sum(DRTCCinPOD49CD8[,4])
nonDRTCCinPOD49CD8 = POD49CD8[POD49CD8$`Amino Acid` %in% nonDRTCC$CDR3.aa,]
sum(nonDRTCCinPOD49CD8[,4])
DRTCCinPOD49CD8_Parent = All_DRTCCs[All_DRTCCs$CDR3.aa %in% DRTCCinPOD49CD8$`Amino Acid`,]
sum(DRTCCinPOD49CD8_Parent$Clones)
nonDRTCCinPOD49CD8_Parent = All_Unstim[All_Unstim$CDR3.aa %in% nonDRTCCinPOD49CD8$`Amino Acid`,]
sum(nonDRTCCinPOD49CD8_Parent$Clones)

write.csv(DRTCCinPOD49CD8, file="DRTCCinPOD49CD8.csv")
write.csv(nonDRTCCinPOD49CD8, file="nonDRTCCinPOD49CD8.csv")
write.csv(DRTCCinPOD49CD8_Parent, file="DRTCCinPOD49CD8_Parent.csv")
write.csv(nonDRTCCinPOD49CD8_Parent, file="nonDRTCCinPOD49CD8_Parent.csv")


# Identifying non-DRTCCs

All_DRTCCs <- read_excel("//172.25.60.69/CCTI_Labs/CCTI_SYKES/EXP FF7- XKLT-001 Exp/Bulk TCRb sequincing Adaptive Data XKLT001/POD33 & 49/All DRTCCs.xlsx")
All_Unstim <- read_excel("//172.25.60.69/CCTI_Labs/CCTI_SYKES/EXP FF7- XKLT-001 Exp/Bulk TCRb sequincing Adaptive Data XKLT001/POD33 & 49/All Unstim.xlsx")
nonDRTCC = All_Unstim[!(All_Unstim[,5] %in% All_DRTCCs[,5]),]
CD4nonDRTCCFrequency = CombinedRearrangements[CD4nonDRTCC,]
write.csv(CD4nonDRTCCFrequency, file="CD4nonDRTCCFrequencyPt1.csv")

UnstimCD8 = CombinedRearrangements[,c(1,7)]
CD8nonDRTCC = rownames(UnstimCD8)[!(rownames(UnstimCD8) %in% rownames(cd8.HVG))]
CD8nonDRTCCFrequency = CombinedRearrangements[CD8nonDRTCC,]
write.csv(CD8nonDRTCCFrequency, file="CD8nonDRTCCFrequency.csv")
nonDRTCCMo8 = CD8nonDRTCCFrequency[,9]
nonDRTCCMo8 = as.data.frame(nonDRTCCMo8)
nonDRTCCMo8 = nonDRTCCMo8[nonDRTCCMo8[,1] > 0,]
nonDRTCCMo8 = as.data.frame(nonDRTCCMo8)
sum(nonDRTCCMo8[,1])

# non DRTCC using templates

CombinedRearrangementsTemplates <- read_csv("CombinedRearrangementsTemplates.csv")
CombinedRearrangementsTemplates = as.data.frame(CombinedRearrangementsTemplates)
rownames(CombinedRearrangementsTemplates) = CombinedRearrangementsTemplates[,1]
CD8nonDRTCCTemplates = CombinedRearrangementsTemplates[CD8nonDRTCC,]
write.csv(CD8nonDRTCCTemplates, file="CD8nonDRTCCTemplates.csv")
sum(CD8nonDRTCCTemplates[,10])

CD8nonDRTCCMo8 = CombinedRearrangementsTemplates[, c(1,10,12)]
CD8nonDRTCCMo8 = CD8nonDRTCCMo8[CD8nonDRTCCMo8[,2] > 0,]
CD8nonDRTCCMo8 = CD8nonDRTCCMo8[rownames(CD8nonDRTCCTemplates),]
CD8nonDRTCCMo8 = na.omit(CD8nonDRTCCMo8)
write.csv(CD8nonDRTCCMo8, file="CD8nonDRTCCMo8.csv")
sum(CD8nonDRTCCMo8[,3])


CD4nonDRTCCTemplates = CombinedRearrangementsTemplates[rownames(CD4nonDRTCCFrequency),]
write.csv(CD4nonDRTCCTemplates, file="CD4nonDRTCCTemplates.csv")
# Alternatively, if we don't want frequency threshold
CD4nonDRTCCTemplates = anti_join(CD4_Unstim, CD4DRTCCTemplates, by = "Amino Acid")
write.csv(CD4nonDRTCCTemplates, file="CD4nonDRTCCTemplates.csv")

# Count number of unique sequences 

CombinedRearrangements_Templates <- read_tsv("//172.25.60.69/CCTI_Labs/CCTI_PANO/PANORAMA Study/PANORAMA Study/TCR sequencing data/Pt 101-001/POD216/CombinedRearrangements_12-22-2023_3-53-27_PM (templates).tsv")
CombinedRearrangements_Templates = as.data.frame(CombinedRearrangements_Templates)
rownames(CombinedRearrangements_Templates)=CombinedRearrangements_Templates[,1]
CombinedRearrangements_Templates = CombinedRearrangements_Templates[setdiff(rownames(CombinedRearrangements_Templates), contaminants.sequences),]
CombinedRearrangements_Templates$CD4 = CombinedRearrangements_Templates[,7]/3 + CombinedRearrangements_Templates[,10]/3 + CombinedRearrangements_Templates[,12]/3
sum(CombinedRearrangements_Templates[,14])
CombinedRearrangements_Templates$CD8 = CombinedRearrangements_Templates[,8]/4 + CombinedRearrangements_Templates[,9]/4 + CombinedRearrangements_Templates[,11]/4 + CombinedRearrangements_Templates[,13]/4
sum(CombinedRearrangements_Templates[,15])
CombinedRearrangements_Templates=resolveambiguous(CombinedRearrangements_Templates, 14, 15)
CD4=which(CombinedRearrangements_Templates[,14]>0)
CD8=which(CombinedRearrangements_Templates[,15]>0)
CombinedRearrangements_Templates[setdiff(1:nrow(CombinedRearrangements_Templates),CD4),c(7,10,12)]=0
CombinedRearrangements_Templates[setdiff(1:nrow(CombinedRearrangements_Templates),CD8),c(8,9,11,13)]=0

length(which(CombinedRearrangements_Templates[,13]>0))
sum(CombinedRearrangements_Templates[,13])

### Statistical analysis of detection rates

# Using unique clones
n1 <- 71165 # DRTCC total unique clones
p1 <- 15695/71165 # Proportion of DRTCC identified in Mo8 CD8

n2 <- 71165 # nonDRTCC total unique clones
p2 <- 26003/71165 # proportion of non DRTCC identified in Mo8 CD8

# Pooled proportion
p_pooled <- (n1 * p1 + n2 * p2) / (n1 + n2)

# Calculate Z-statistic
z_stat <- (p1 - p2) / sqrt(p_pooled * (1 - p_pooled) * (1/n1 + 1/n2))

# Two-tailed test, find critical Z-value for alpha = 0.05
critical_value <- qnorm(1 - 0.05/2)

# Compare Z-statistic to critical value
if (abs(z_stat) > critical_value) {
  cat("Reject the null hypothesis\n")
} else {
  cat("Fail to reject the null hypothesis\n")
}

# Alternatively, calculate the p-value
p_value <- 2 * (1 - pnorm(abs(z_stat)))

# Compare p-value to alpha
if (p_value < 0.05) {
  cat("Reject the null hypothesis\n")
} else {
  cat("Fail to reject the null hypothesis\n")
}

### Using number of templates

n1 <- 179968  # DRTCC total number of templates
p1 <- 251/179968  # Proportion of DRTCC identified in Mo8 CD8

n2 <- 179968   # nonDRTCC total unique clones
p2 <- 161147/179968  # proportion of non DRTCC identified in Mo8 CD8

# Pooled proportion
p_pooled <- (n1 * p1 + n2 * p2) / (n1 + n2)

# Calculate Z-statistic
z_stat <- (p1 - p2) / sqrt(p_pooled * (1 - p_pooled) * (1/n1 + 1/n2))

# Two-tailed test, find critical Z-value for alpha = 0.05
critical_value <- qnorm(1 - 0.05/2)

# Compare Z-statistic to critical value
if (abs(z_stat) > critical_value) {
  cat("Reject the null hypothesis\n")
} else {
  cat("Fail to reject the null hypothesis\n")
}

# Alternatively, calculate the p-value
p_value <- 2 * (1 - pnorm(abs(z_stat)))

# Compare p-value to alpha
if (p_value < 0.05) {
  cat("Reject the null hypothesis\n")
} else {
  cat("Fail to reject the null hypothesis\n")
}


# CFSE Low that don't meet DRTCC Criteria

CD8DRTCCTemplates <- read_csv("CD8DRTCCTemplates.csv")
CombinedRearrangementsTemplates <- read_csv("CombinedRearrangementsTemplates.csv")
CFSELow = CombinedRearrangementsTemplates[CombinedRearrangementsTemplates[,12]>0,]
CFSELow_nonDRTCC = anti_join(CFSELow, CD8DRTCCTemplates, by = "Amino Acid")

write.csv(CFSELow, file="CFSELow.csv")
write.csv(CFSELow_nonDRTCC, file="CFSELow_nonDRTCC.csv")

names(PANORAMA_Pt1)
#alloreactive
cd4.HVG=CombinedRearrangements[,c(12,10)]
cd8.HVG=CombinedRearrangements[,c(13,11)]
allo.HVG=listAlloreactive(cd4.HVG,cd8.HVG, fold=2, freq=0.0001)
allo.HVG = na.omit(allo.HVG)
length(allo.HVG[[1]])
length(allo.HVG[[2]])
sum(PANORAMA_Pt1[listAlloreactive(cd4.HVG,cd8.HVG)[[1]],4])
sum(PANORAMA_Pt1[listAlloreactive(cd4.HVG,cd8.HVG)[[2]],5])
View(PANORAMA_Pt1)

#unstim
unstimCD4=PANORAMA_Pt1[PANORAMA_Pt1[,6]>0,1]
unstimCD8=PANORAMA_Pt1[PANORAMA_Pt1[,7]>0,1]

length(intersect(unstimCD4,allo.HVG[[1]]))
length(intersect(unstimCD8,allo.HVG[[2]]))

unstimCD4=setdiff(unstimCD4,allo.HVG[[1]])
unstimCD8=setdiff(unstimCD8,allo.HVG[[2]])
length(unstimCD4)
length(unstimCD8)
sum(PANORAMA_Pt1[unstimCD4,6])
sum(PANORAMA_Pt1[unstimCD8,7])
names(PANORAMA_Pt1)

#abundanceplots
data=PANORAMA_Pt1
abundancePlot(normalize(data[,c(4,6)]))
abundancePlot(normalize(data[,c(5,7)]))

#anothermethodbutKevinhasneverused
get_abundance(data = data, output = 'plot', method = 'Obradovic')

#clonality
cloneCal(PANORAMA_Pt1[,4])
cloneCal(PANORAMA_Pt1[,6])
cloneCal(PANORAMA_Pt1[,5])
cloneCal(PANORAMA_Pt1[,7])

#R20
getR20(PANORAMA_Pt1[,4])
getR20(PANORAMA_Pt1[,6])
getR20(PANORAMA_Pt1[,5])
getR20(PANORAMA_Pt1[,7])

#nonHVGmappable
rCD4mappable=rownames(PANORAMA_Pt1[(PANORAMA_Pt1[,12]+PANORAMA_Pt1[,14])>0,])
rCD8mappable=rownames(PANORAMA_Pt1[(PANORAMA_Pt1[,13]+PANORAMA_Pt1[,15])>0,])
CD4nonHVG=setdiff(rCD4mappable,allo.HVG[[1]])
CD8nonHVG=setdiff(rCD8mappable,allo.HVG[[2]])
length(rCD4mappable)
length(rCD8mappable)
length(CD4nonHVG)
length(CD8nonHVG)

write.table(CD4.allo)

write.table(CD4.allo,file ="P:/CCTI_USERS/Kevin Breen/UPitt DCreg Project/Kevin R code/1001/Alloclones.txt",quote=F,row.names=F,col.names=F, sep="\t")

#post tx blood definitions and overlap with pre-tx
#6mth blood
CD4_6mth=PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,1]
CD8_6mth=PANORAMA_Pt1[PANORAMA_Pt1[,11]>0,1]
length(CD4_6mth)
length(CD8_6mth)
sum(PANORAMA_Pt1[CD4_6mth,10])
sum(PANORAMA_Pt1[CD8_6mth,11])

length(intersect(CD4_6mth,unstimCD4))
length(intersect(CD8_6mth,unstimCD8))
length(intersect(CD4_6mth,allo.HVG[[1]]))
length(intersect(CD8_6mth,allo.HVG[[2]]))
length(intersect(CD4_6mth,CD4nonHVG))
length(intersect(CD8_6mth,CD8nonHVG))
sum(PANORAMA_Pt1[unstimCD4,10])
sum(PANORAMA_Pt1[unstimCD8,11])
sum(PANORAMA_Pt1[allo.HVG[[1]],10])
sum(PANORAMA_Pt1[allo.HVG[[2]],11])
sum(PANORAMA_Pt1[CD4nonHVG,10])
sum(PANORAMA_Pt1[CD8nonHVG,11])

cd4.allo.freq= sum(PANORAMA_Pt1[CD4.allo,10])/ sum(PANORAMA_Pt1[,10])

#6m posttx blood HvG clone

cd4= PANORAMA_Pt1[,c(14,12,10)]

cd8= PANORAMA_Pt1[,c(15,13,11)]

run(cd4,cd8, fold=2,freq=0.0001, filename='P:/CCTI_USERS/Kevin Breen/UPitt DCreg Project/Kevin R code/1001/6mHvG.txt')



rownames(PANORAMA_Pt1)=PANORAMA_Pt1[,1]

cd4= PANORAMA_Pt1[,c(14,12)]

cd8= PANORAMA_Pt1[,c(15,13)]

allo=listAlloreactive(cd4, cd8, fold=2,freq=0.0001)

intersect(allo[[1]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))

intersect(allo[[2]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))


#blood 12mth
CD4_12mth=PANORAMA_Pt1[PANORAMA_Pt1[,4]>0,1]
CD8_12mth=PANORAMA_Pt1[PANORAMA_Pt1[,5]>0,1]
length(CD4_12mth)
length(CD8_12mth)
sum(PANORAMA_Pt1[CD4_12mth,4])
sum(PANORAMA_Pt1[CD8_12mth,5])

length(intersect(CD4_12mth,unstimCD4))
length(intersect(CD8_12mth,unstimCD8))
length(intersect(CD4_12mth,allo.HVG[[1]]))
length(intersect(CD8_12mth,allo.HVG[[2]]))
length(intersect(CD4_12mth,CD4nonHVG))
length(intersect(CD8_12mth,CD8nonHVG))
sum(PANORAMA_Pt1[unstimCD4,4])
sum(PANORAMA_Pt1[unstimCD8,5])
sum(PANORAMA_Pt1[allo.HVG[[1]],4])
sum(PANORAMA_Pt1[allo.HVG[[2]],5])
sum(PANORAMA_Pt1[CD4nonHVG,4])
sum(PANORAMA_Pt1[CD8nonHVG,5])

#12m posttx blood HvG clone

cd4= PANORAMA_Pt1[,c(14,12,4)]

cd8= PANORAMA_Pt1[,c(15,13,5)]

run(cd4,cd8, fold=2,freq=0.0001, filename='P:/CCTI_USERS/Kevin Breen/UPitt DCreg Project/Kevin R code/1001/12mHvG.txt')


rownames(PANORAMA_Pt1)=PANORAMA_Pt1[,1]

cd4= PANORAMA_Pt1[,c(14,12)]

cd8= PANORAMA_Pt1[,c(15,13)]

allo=listAlloreactive(cd4, cd8, fold=2,freq=0.0001)

intersect(allo[[1]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))

intersect(allo[[2]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))

#blood 18mth
CD4_18mth=PANORAMA_Pt1[PANORAMA_Pt1[,7]>0,1]
CD8_18mth=PANORAMA_Pt1[PANORAMA_Pt1[,8]>0,1]
length(CD4_18mth)
length(CD8_18mth)
sum(PANORAMA_Pt1[CD4_18mth,7])
sum(PANORAMA_Pt1[CD8_18mth,8])

length(intersect(CD4_18mth,unstimCD4))
length(intersect(CD8_18mth,unstimCD8))
length(intersect(CD4_18mth,allo.HVG[[1]]))
length(intersect(CD8_18mth,allo.HVG[[2]]))
length(intersect(CD4_18mth,CD4nonHVG))
length(intersect(CD8_18mth,CD8nonHVG))
sum(PANORAMA_Pt1[unstimCD4,7])
sum(PANORAMA_Pt1[unstimCD8,8])
sum(PANORAMA_Pt1[allo.HVG[[1]],7])
sum(PANORAMA_Pt1[allo.HVG[[2]],8])
sum(PANORAMA_Pt1[CD4nonHVG,7])
sum(PANORAMA_Pt1[CD8nonHVG,8])


names(PANORAMA_Pt1)

#18m posttx blood HvG clone

cd4= PANORAMA_Pt1[,c(14,12,7)]

cd8= PANORAMA_Pt1[,c(15,13,8)]

run(cd4,cd8, fold=2,freq=0.0001, filename='P:/CCTI_USERS/Kevin Breen/UPitt DCreg Project/Kevin R code/1001/18mHvG.txt')


rownames(PANORAMA_Pt1)=PANORAMA_Pt1[,1]

cd4= PANORAMA_Pt1[,c(14,12)]

cd8= PANORAMA_Pt1[,c(15,13)]

allo=listAlloreactive(cd4, cd8, fold=2,freq=0.0001)

intersect(allo[[1]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))

intersect(allo[[2]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))


#liver bx 12 mth
liverBx_12mth=PANORAMA_Pt1[PANORAMA_Pt1[,6]>0,1]
length(liverBx_12mth)
sum(PANORAMA_Pt1[liverBx_12mth,6])
length(intersect(liverBx_12mth,unstimCD4))
length(intersect(liverBx_12mth,unstimCD8))
length(intersect(liverBx_12mth,CD4.allo))
length(intersect(liverBx_12mth,CD8.allo))
length(intersect(liverBx_12mth,CD4_6mth))
length(intersect(liverBx_12mth,CD8_6mth))
length(intersect(liverBx_12mth,CD4_12mth))
length(intersect(liverBx_12mth,CD8_12mth))
length(intersect(liverBx_12mth,CD4_18mth))
length(intersect(liverBx_12mth,CD8_18mth))
sum(PANORAMA_Pt1[unstimCD4,6])
sum(PANORAMA_Pt1[unstimCD8,6])
sum(PANORAMA_Pt1[CD4.allo,6])
sum(PANORAMA_Pt1[CD8.allo,6])
sum(PANORAMA_Pt1[CD4_6mth,6])
sum(PANORAMA_Pt1[CD8_6mth,6])
sum(PANORAMA_Pt1[CD4_12mth,6])
sum(PANORAMA_Pt1[CD8_12mth,6])
sum(PANORAMA_Pt1[CD4_18mth,6])
sum(PANORAMA_Pt1[CD8_18mth,6])

#liver bx 18 mth
liverBx_18mth=PANORAMA_Pt1[PANORAMA_Pt1[,9]>0,1]
length(liverBx_18mth)
sum(PANORAMA_Pt1[liverBx_18mth,9])
length(intersect(liverBx_18mth,unstimCD4))
length(intersect(liverBx_18mth,unstimCD8))
length(intersect(liverBx_18mth,CD4.allo))
length(intersect(liverBx_18mth,CD8.allo))
length(intersect(liverBx_18mth,CD4_6mth))
length(intersect(liverBx_18mth,CD8_6mth))
length(intersect(liverBx_18mth,CD4_12mth))
length(intersect(liverBx_18mth,CD8_12mth))
length(intersect(liverBx_18mth,CD4_18mth))
length(intersect(liverBx_18mth,CD8_18mth))
sum(PANORAMA_Pt1[unstimCD4,9])
sum(PANORAMA_Pt1[unstimCD8,9])
sum(PANORAMA_Pt1[CD4.allo,9])
sum(PANORAMA_Pt1[CD8.allo,9])
sum(PANORAMA_Pt1[CD4_6mth,9])
sum(PANORAMA_Pt1[CD8_6mth,9])
sum(PANORAMA_Pt1[CD4_12mth,9])
sum(PANORAMA_Pt1[CD8_12mth,9])
sum(PANORAMA_Pt1[CD4_18mth,9])
sum(PANORAMA_Pt1[CD8_18mth,9])

#combine columns for JSD analysis
PANORAMA_Pt1$unstim = PANORAMA_Pt1[,14]+PANORAMA_Pt1[,15]
PANORAMA_Pt1$blood6 = PANORAMA_Pt1[,10]+PANORAMA_Pt1[,11]
PANORAMA_Pt1$blood12 = PANORAMA_Pt1[,4]+PANORAMA_Pt1[,5]
PANORAMA_Pt1$blood18 = PANORAMA_Pt1[,7]+PANORAMA_Pt1[,8]

names(PANORAMA_Pt1)

jsdReport(PANORAMA_Pt1[,c(18,6,9,19,20,21)],topN=1000)
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:
  
  ```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:
  
  ```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

View(Pt_101_001_Pre_Tx_MLR_and_Unstim)


library(dplyr)
library(tidyverse)
library(ggplot2)

PANORAMA_Pt1= as.data.frame(Pt_101_001_Pre_Tx_MLR_and_Unstim)

rownames(PANORAMA_Pt1)=PANORAMA_Pt1[,1]

#remove MART1 contamination
PANORAMA_Pt1=PANORAMA_Pt1[setdiff(rownames(PANORAMA_Pt1),"ACGTTGGCGTCTGCTGTACCCTCTCAGACATCTGTGTACTTCTGTGCCAGCAGCCTAAGTTTCGGCACTGAAGCTTTCTTTGGACAA"),]

names(PANORAMA_Pt1)
PANORAMA_Pt1$CD4 = PANORAMA_Pt1[,4]+PANORAMA_Pt1[,6]
PANORAMA_Pt1$CD8 = PANORAMA_Pt1[,5]+PANORAMA_Pt1[,7]
PANORAMA_Pt1=resolveambiguous(PANORAMA_Pt1,8,9)
CD4=which(PANORAMA_Pt1[,8]>0)
CD8=which(PANORAMA_Pt1[,9]>0)

#removing ambig from all other CD4/CD8 columns
PANORAMA_Pt1[setdiff(1:nrow(PANORAMA_Pt1),CD4),c(4,6)]=0
PANORAMA_Pt1[setdiff(1:nrow(PANORAMA_Pt1),CD8),c(5,7)]=0

names(PANORAMA_Pt1)
#alloreactive
cd4.HVG=PANORAMA_Pt1[,c(6,4)]
cd8.HVG=PANORAMA_Pt1[,c(7,5)]
allo.HVG=listAlloreactive(cd4.HVG,cd8.HVG, fold=2, freq=0.00002)
length(allo.HVG[[1]])
length(allo.HVG[[2]])
sum(PANORAMA_Pt1[listAlloreactive(cd4.HVG,cd8.HVG)[[1]],4])
sum(PANORAMA_Pt1[listAlloreactive(cd4.HVG,cd8.HVG)[[2]],5])
View(PANORAMA_Pt1)

#unstim
unstimCD4=PANORAMA_Pt1[PANORAMA_Pt1[,6]>0,1]
unstimCD8=PANORAMA_Pt1[PANORAMA_Pt1[,7]>0,1]

length(intersect(unstimCD4,allo.HVG[[1]]))
length(intersect(unstimCD8,allo.HVG[[2]]))

unstimCD4=setdiff(unstimCD4,allo.HVG[[1]])
unstimCD8=setdiff(unstimCD8,allo.HVG[[2]])
length(unstimCD4)
length(unstimCD8)
sum(PANORAMA_Pt1[unstimCD4,6])
sum(PANORAMA_Pt1[unstimCD8,7])
names(PANORAMA_Pt1)

#abundanceplots
data=PANORAMA_Pt1
abundancePlot(normalize(data[,c(4,6)]))
abundancePlot(normalize(data[,c(5,7)]))

#anothermethodbutKevinhasneverused
get_abundance(data = data, output = 'plot', method = 'Obradovic')

#clonality
cloneCal(PANORAMA_Pt1[,4])
cloneCal(PANORAMA_Pt1[,6])
cloneCal(PANORAMA_Pt1[,5])
cloneCal(PANORAMA_Pt1[,7])

#R20
getR20(PANORAMA_Pt1[,4])
getR20(PANORAMA_Pt1[,6])
getR20(PANORAMA_Pt1[,5])
getR20(PANORAMA_Pt1[,7])

#nonHVGmappable
rCD4mappable=rownames(PANORAMA_Pt1[(PANORAMA_Pt1[,12]+PANORAMA_Pt1[,14])>0,])
rCD8mappable=rownames(PANORAMA_Pt1[(PANORAMA_Pt1[,13]+PANORAMA_Pt1[,15])>0,])
CD4nonHVG=setdiff(rCD4mappable,allo.HVG[[1]])
CD8nonHVG=setdiff(rCD8mappable,allo.HVG[[2]])
length(rCD4mappable)
length(rCD8mappable)
length(CD4nonHVG)
length(CD8nonHVG)

write.table(CD4.allo)

write.table(CD4.allo,file ="P:/CCTI_USERS/Kevin Breen/UPitt DCreg Project/Kevin R code/1001/Alloclones.txt",quote=F,row.names=F,col.names=F, sep="\t")

#post tx blood definitions and overlap with pre-tx
#6mth blood
CD4_6mth=PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,1]
CD8_6mth=PANORAMA_Pt1[PANORAMA_Pt1[,11]>0,1]
length(CD4_6mth)
length(CD8_6mth)
sum(PANORAMA_Pt1[CD4_6mth,10])
sum(PANORAMA_Pt1[CD8_6mth,11])

length(intersect(CD4_6mth,unstimCD4))
length(intersect(CD8_6mth,unstimCD8))
length(intersect(CD4_6mth,allo.HVG[[1]]))
length(intersect(CD8_6mth,allo.HVG[[2]]))
length(intersect(CD4_6mth,CD4nonHVG))
length(intersect(CD8_6mth,CD8nonHVG))
sum(PANORAMA_Pt1[unstimCD4,10])
sum(PANORAMA_Pt1[unstimCD8,11])
sum(PANORAMA_Pt1[allo.HVG[[1]],10])
sum(PANORAMA_Pt1[allo.HVG[[2]],11])
sum(PANORAMA_Pt1[CD4nonHVG,10])
sum(PANORAMA_Pt1[CD8nonHVG,11])

cd4.allo.freq= sum(PANORAMA_Pt1[CD4.allo,10])/ sum(PANORAMA_Pt1[,10])

#6m posttx blood HvG clone

cd4= PANORAMA_Pt1[,c(14,12,10)]

cd8= PANORAMA_Pt1[,c(15,13,11)]

run(cd4,cd8, fold=2,freq=0.0001, filename='P:/CCTI_USERS/Kevin Breen/UPitt DCreg Project/Kevin R code/1001/6mHvG.txt')



rownames(PANORAMA_Pt1)=PANORAMA_Pt1[,1]

cd4= PANORAMA_Pt1[,c(14,12)]

cd8= PANORAMA_Pt1[,c(15,13)]

allo=listAlloreactive(cd4, cd8, fold=2,freq=0.0001)

intersect(allo[[1]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))

intersect(allo[[2]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))


#blood 12mth
CD4_12mth=PANORAMA_Pt1[PANORAMA_Pt1[,4]>0,1]
CD8_12mth=PANORAMA_Pt1[PANORAMA_Pt1[,5]>0,1]
length(CD4_12mth)
length(CD8_12mth)
sum(PANORAMA_Pt1[CD4_12mth,4])
sum(PANORAMA_Pt1[CD8_12mth,5])

length(intersect(CD4_12mth,unstimCD4))
length(intersect(CD8_12mth,unstimCD8))
length(intersect(CD4_12mth,allo.HVG[[1]]))
length(intersect(CD8_12mth,allo.HVG[[2]]))
length(intersect(CD4_12mth,CD4nonHVG))
length(intersect(CD8_12mth,CD8nonHVG))
sum(PANORAMA_Pt1[unstimCD4,4])
sum(PANORAMA_Pt1[unstimCD8,5])
sum(PANORAMA_Pt1[allo.HVG[[1]],4])
sum(PANORAMA_Pt1[allo.HVG[[2]],5])
sum(PANORAMA_Pt1[CD4nonHVG,4])
sum(PANORAMA_Pt1[CD8nonHVG,5])

#12m posttx blood HvG clone

cd4= PANORAMA_Pt1[,c(14,12,4)]

cd8= PANORAMA_Pt1[,c(15,13,5)]

run(cd4,cd8, fold=2,freq=0.0001, filename='P:/CCTI_USERS/Kevin Breen/UPitt DCreg Project/Kevin R code/1001/12mHvG.txt')


rownames(PANORAMA_Pt1)=PANORAMA_Pt1[,1]

cd4= PANORAMA_Pt1[,c(14,12)]

cd8= PANORAMA_Pt1[,c(15,13)]

allo=listAlloreactive(cd4, cd8, fold=2,freq=0.0001)

intersect(allo[[1]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))

intersect(allo[[2]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))

#blood 18mth
CD4_18mth=PANORAMA_Pt1[PANORAMA_Pt1[,7]>0,1]
CD8_18mth=PANORAMA_Pt1[PANORAMA_Pt1[,8]>0,1]
length(CD4_18mth)
length(CD8_18mth)
sum(PANORAMA_Pt1[CD4_18mth,7])
sum(PANORAMA_Pt1[CD8_18mth,8])

length(intersect(CD4_18mth,unstimCD4))
length(intersect(CD8_18mth,unstimCD8))
length(intersect(CD4_18mth,allo.HVG[[1]]))
length(intersect(CD8_18mth,allo.HVG[[2]]))
length(intersect(CD4_18mth,CD4nonHVG))
length(intersect(CD8_18mth,CD8nonHVG))
sum(PANORAMA_Pt1[unstimCD4,7])
sum(PANORAMA_Pt1[unstimCD8,8])
sum(PANORAMA_Pt1[allo.HVG[[1]],7])
sum(PANORAMA_Pt1[allo.HVG[[2]],8])
sum(PANORAMA_Pt1[CD4nonHVG,7])
sum(PANORAMA_Pt1[CD8nonHVG,8])


names(PANORAMA_Pt1)

#18m posttx blood HvG clone

cd4= PANORAMA_Pt1[,c(14,12,7)]

cd8= PANORAMA_Pt1[,c(15,13,8)]

run(cd4,cd8, fold=2,freq=0.0001, filename='P:/CCTI_USERS/Kevin Breen/UPitt DCreg Project/Kevin R code/1001/18mHvG.txt')


rownames(PANORAMA_Pt1)=PANORAMA_Pt1[,1]

cd4= PANORAMA_Pt1[,c(14,12)]

cd8= PANORAMA_Pt1[,c(15,13)]

allo=listAlloreactive(cd4, cd8, fold=2,freq=0.0001)

intersect(allo[[1]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))

intersect(allo[[2]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))


#liver bx 12 mth
liverBx_12mth=PANORAMA_Pt1[PANORAMA_Pt1[,6]>0,1]
length(liverBx_12mth)
sum(PANORAMA_Pt1[liverBx_12mth,6])
length(intersect(liverBx_12mth,unstimCD4))
length(intersect(liverBx_12mth,unstimCD8))
length(intersect(liverBx_12mth,CD4.allo))
length(intersect(liverBx_12mth,CD8.allo))
length(intersect(liverBx_12mth,CD4_6mth))
length(intersect(liverBx_12mth,CD8_6mth))
length(intersect(liverBx_12mth,CD4_12mth))
length(intersect(liverBx_12mth,CD8_12mth))
length(intersect(liverBx_12mth,CD4_18mth))
length(intersect(liverBx_12mth,CD8_18mth))
sum(PANORAMA_Pt1[unstimCD4,6])
sum(PANORAMA_Pt1[unstimCD8,6])
sum(PANORAMA_Pt1[CD4.allo,6])
sum(PANORAMA_Pt1[CD8.allo,6])
sum(PANORAMA_Pt1[CD4_6mth,6])
sum(PANORAMA_Pt1[CD8_6mth,6])
sum(PANORAMA_Pt1[CD4_12mth,6])
sum(PANORAMA_Pt1[CD8_12mth,6])
sum(PANORAMA_Pt1[CD4_18mth,6])
sum(PANORAMA_Pt1[CD8_18mth,6])

#liver bx 18 mth
liverBx_18mth=PANORAMA_Pt1[PANORAMA_Pt1[,9]>0,1]
length(liverBx_18mth)
sum(PANORAMA_Pt1[liverBx_18mth,9])
length(intersect(liverBx_18mth,unstimCD4))
length(intersect(liverBx_18mth,unstimCD8))
length(intersect(liverBx_18mth,CD4.allo))
length(intersect(liverBx_18mth,CD8.allo))
length(intersect(liverBx_18mth,CD4_6mth))
length(intersect(liverBx_18mth,CD8_6mth))
length(intersect(liverBx_18mth,CD4_12mth))
length(intersect(liverBx_18mth,CD8_12mth))
length(intersect(liverBx_18mth,CD4_18mth))
length(intersect(liverBx_18mth,CD8_18mth))
sum(PANORAMA_Pt1[unstimCD4,9])
sum(PANORAMA_Pt1[unstimCD8,9])
sum(PANORAMA_Pt1[CD4.allo,9])
sum(PANORAMA_Pt1[CD8.allo,9])
sum(PANORAMA_Pt1[CD4_6mth,9])
sum(PANORAMA_Pt1[CD8_6mth,9])
sum(PANORAMA_Pt1[CD4_12mth,9])
sum(PANORAMA_Pt1[CD8_12mth,9])
sum(PANORAMA_Pt1[CD4_18mth,9])
sum(PANORAMA_Pt1[CD8_18mth,9])

#combine columns for JSD analysis
PANORAMA_Pt1$unstim = PANORAMA_Pt1[,14]+PANORAMA_Pt1[,15]
PANORAMA_Pt1$blood6 = PANORAMA_Pt1[,10]+PANORAMA_Pt1[,11]
PANORAMA_Pt1$blood12 = PANORAMA_Pt1[,4]+PANORAMA_Pt1[,5]
PANORAMA_Pt1$blood18 = PANORAMA_Pt1[,7]+PANORAMA_Pt1[,8]

names(PANORAMA_Pt1)

jsdReport(PANORAMA_Pt1[,c(18,6,9,19,20,21)],topN=1000)
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:
  
  ```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:
  
  ```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

View(Pt_101_001_Pre_Tx_MLR_and_Unstim)


library(dplyr)
library(tidyverse)
library(ggplot2)

PANORAMA_Pt1= as.data.frame(Pt_101_001_Pre_Tx_MLR_and_Unstim)

rownames(PANORAMA_Pt1)=PANORAMA_Pt1[,1]

#remove MART1 contamination
PANORAMA_Pt1=PANORAMA_Pt1[setdiff(rownames(PANORAMA_Pt1),"ACGTTGGCGTCTGCTGTACCCTCTCAGACATCTGTGTACTTCTGTGCCAGCAGCCTAAGTTTCGGCACTGAAGCTTTCTTTGGACAA"),]

names(PANORAMA_Pt1)
PANORAMA_Pt1$CD4 = PANORAMA_Pt1[,4]+PANORAMA_Pt1[,6]
PANORAMA_Pt1$CD8 = PANORAMA_Pt1[,5]+PANORAMA_Pt1[,7]
PANORAMA_Pt1=resolveambiguous(PANORAMA_Pt1,8,9)
CD4=which(PANORAMA_Pt1[,8]>0)
CD8=which(PANORAMA_Pt1[,9]>0)

#removing ambig from all other CD4/CD8 columns
PANORAMA_Pt1[setdiff(1:nrow(PANORAMA_Pt1),CD4),c(4,6)]=0
PANORAMA_Pt1[setdiff(1:nrow(PANORAMA_Pt1),CD8),c(5,7)]=0

names(PANORAMA_Pt1)
#alloreactive
cd4.HVG=PANORAMA_Pt1[,c(6,4)]
cd8.HVG=PANORAMA_Pt1[,c(7,5)]
allo.HVG=listAlloreactive(cd4.HVG,cd8.HVG, fold=6, freq=0.001)
length(allo.HVG[[1]])
length(allo.HVG[[2]])
sum(PANORAMA_Pt1[listAlloreactive(cd4.HVG,cd8.HVG)[[1]],4])
sum(PANORAMA_Pt1[listAlloreactive(cd4.HVG,cd8.HVG)[[2]],5])
View(PANORAMA_Pt1)

#unstim
unstimCD4=PANORAMA_Pt1[PANORAMA_Pt1[,6]>0,1]
unstimCD8=PANORAMA_Pt1[PANORAMA_Pt1[,7]>0,1]

length(intersect(unstimCD4,allo.HVG[[1]]))
length(intersect(unstimCD8,allo.HVG[[2]]))

unstimCD4=setdiff(unstimCD4,allo.HVG[[1]])
unstimCD8=setdiff(unstimCD8,allo.HVG[[2]])
length(unstimCD4)
length(unstimCD8)
sum(PANORAMA_Pt1[unstimCD4,6])
sum(PANORAMA_Pt1[unstimCD8,7])
names(PANORAMA_Pt1)

#abundanceplots
data=PANORAMA_Pt1
abundancePlot(normalize(data[,c(4,6)]))
abundancePlot(normalize(data[,c(5,7)]))

#anothermethodbutKevinhasneverused
get_abundance(data = data, output = 'plot', method = 'Obradovic')

#clonality
cloneCal(PANORAMA_Pt1[,4])
cloneCal(PANORAMA_Pt1[,6])
cloneCal(PANORAMA_Pt1[,5])
cloneCal(PANORAMA_Pt1[,7])

#R20
getR20(PANORAMA_Pt1[,4])
getR20(PANORAMA_Pt1[,6])
getR20(PANORAMA_Pt1[,5])
getR20(PANORAMA_Pt1[,7])

#nonHVGmappable
rCD4mappable=rownames(PANORAMA_Pt1[(PANORAMA_Pt1[,12]+PANORAMA_Pt1[,14])>0,])
rCD8mappable=rownames(PANORAMA_Pt1[(PANORAMA_Pt1[,13]+PANORAMA_Pt1[,15])>0,])
CD4nonHVG=setdiff(rCD4mappable,allo.HVG[[1]])
CD8nonHVG=setdiff(rCD8mappable,allo.HVG[[2]])
length(rCD4mappable)
length(rCD8mappable)
length(CD4nonHVG)
length(CD8nonHVG)

write.table(CD4.allo)

write.table(CD4.allo,file ="P:/CCTI_USERS/Kevin Breen/UPitt DCreg Project/Kevin R code/1001/Alloclones.txt",quote=F,row.names=F,col.names=F, sep="\t")

#post tx blood definitions and overlap with pre-tx
#6mth blood
CD4_6mth=PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,1]
CD8_6mth=PANORAMA_Pt1[PANORAMA_Pt1[,11]>0,1]
length(CD4_6mth)
length(CD8_6mth)
sum(PANORAMA_Pt1[CD4_6mth,10])
sum(PANORAMA_Pt1[CD8_6mth,11])

length(intersect(CD4_6mth,unstimCD4))
length(intersect(CD8_6mth,unstimCD8))
length(intersect(CD4_6mth,allo.HVG[[1]]))
length(intersect(CD8_6mth,allo.HVG[[2]]))
length(intersect(CD4_6mth,CD4nonHVG))
length(intersect(CD8_6mth,CD8nonHVG))
sum(PANORAMA_Pt1[unstimCD4,10])
sum(PANORAMA_Pt1[unstimCD8,11])
sum(PANORAMA_Pt1[allo.HVG[[1]],10])
sum(PANORAMA_Pt1[allo.HVG[[2]],11])
sum(PANORAMA_Pt1[CD4nonHVG,10])
sum(PANORAMA_Pt1[CD8nonHVG,11])

cd4.allo.freq= sum(PANORAMA_Pt1[CD4.allo,10])/ sum(PANORAMA_Pt1[,10])

#6m posttx blood HvG clone

cd4= PANORAMA_Pt1[,c(14,12,10)]

cd8= PANORAMA_Pt1[,c(15,13,11)]

run(cd4,cd8, fold=2,freq=0.0001, filename='P:/CCTI_USERS/Kevin Breen/UPitt DCreg Project/Kevin R code/1001/6mHvG.txt')



rownames(PANORAMA_Pt1)=PANORAMA_Pt1[,1]

cd4= PANORAMA_Pt1[,c(14,12)]

cd8= PANORAMA_Pt1[,c(15,13)]

allo=listAlloreactive(cd4, cd8, fold=2,freq=0.0001)

intersect(allo[[1]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))

intersect(allo[[2]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))


#blood 12mth
CD4_12mth=PANORAMA_Pt1[PANORAMA_Pt1[,4]>0,1]
CD8_12mth=PANORAMA_Pt1[PANORAMA_Pt1[,5]>0,1]
length(CD4_12mth)
length(CD8_12mth)
sum(PANORAMA_Pt1[CD4_12mth,4])
sum(PANORAMA_Pt1[CD8_12mth,5])

length(intersect(CD4_12mth,unstimCD4))
length(intersect(CD8_12mth,unstimCD8))
length(intersect(CD4_12mth,allo.HVG[[1]]))
length(intersect(CD8_12mth,allo.HVG[[2]]))
length(intersect(CD4_12mth,CD4nonHVG))
length(intersect(CD8_12mth,CD8nonHVG))
sum(PANORAMA_Pt1[unstimCD4,4])
sum(PANORAMA_Pt1[unstimCD8,5])
sum(PANORAMA_Pt1[allo.HVG[[1]],4])
sum(PANORAMA_Pt1[allo.HVG[[2]],5])
sum(PANORAMA_Pt1[CD4nonHVG,4])
sum(PANORAMA_Pt1[CD8nonHVG,5])

#12m posttx blood HvG clone

cd4= PANORAMA_Pt1[,c(14,12,4)]

cd8= PANORAMA_Pt1[,c(15,13,5)]

run(cd4,cd8, fold=2,freq=0.0001, filename='P:/CCTI_USERS/Kevin Breen/UPitt DCreg Project/Kevin R code/1001/12mHvG.txt')


rownames(PANORAMA_Pt1)=PANORAMA_Pt1[,1]

cd4= PANORAMA_Pt1[,c(14,12)]

cd8= PANORAMA_Pt1[,c(15,13)]

allo=listAlloreactive(cd4, cd8, fold=2,freq=0.0001)

intersect(allo[[1]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))

intersect(allo[[2]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))

#blood 18mth
CD4_18mth=PANORAMA_Pt1[PANORAMA_Pt1[,7]>0,1]
CD8_18mth=PANORAMA_Pt1[PANORAMA_Pt1[,8]>0,1]
length(CD4_18mth)
length(CD8_18mth)
sum(PANORAMA_Pt1[CD4_18mth,7])
sum(PANORAMA_Pt1[CD8_18mth,8])

length(intersect(CD4_18mth,unstimCD4))
length(intersect(CD8_18mth,unstimCD8))
length(intersect(CD4_18mth,allo.HVG[[1]]))
length(intersect(CD8_18mth,allo.HVG[[2]]))
length(intersect(CD4_18mth,CD4nonHVG))
length(intersect(CD8_18mth,CD8nonHVG))
sum(PANORAMA_Pt1[unstimCD4,7])
sum(PANORAMA_Pt1[unstimCD8,8])
sum(PANORAMA_Pt1[allo.HVG[[1]],7])
sum(PANORAMA_Pt1[allo.HVG[[2]],8])
sum(PANORAMA_Pt1[CD4nonHVG,7])
sum(PANORAMA_Pt1[CD8nonHVG,8])


names(PANORAMA_Pt1)

#18m posttx blood HvG clone

cd4= PANORAMA_Pt1[,c(14,12,7)]

cd8= PANORAMA_Pt1[,c(15,13,8)]

run(cd4,cd8, fold=2,freq=0.0001, filename='P:/CCTI_USERS/Kevin Breen/UPitt DCreg Project/Kevin R code/1001/18mHvG.txt')


rownames(PANORAMA_Pt1)=PANORAMA_Pt1[,1]

cd4= PANORAMA_Pt1[,c(14,12)]

cd8= PANORAMA_Pt1[,c(15,13)]

allo=listAlloreactive(cd4, cd8, fold=2,freq=0.0001)

intersect(allo[[1]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))

intersect(allo[[2]], rownames(PANORAMA_Pt1[PANORAMA_Pt1[,10]>0,]))


#liver bx 12 mth
liverBx_12mth=PANORAMA_Pt1[PANORAMA_Pt1[,6]>0,1]
length(liverBx_12mth)
sum(PANORAMA_Pt1[liverBx_12mth,6])
length(intersect(liverBx_12mth,unstimCD4))
length(intersect(liverBx_12mth,unstimCD8))
length(intersect(liverBx_12mth,CD4.allo))
length(intersect(liverBx_12mth,CD8.allo))
length(intersect(liverBx_12mth,CD4_6mth))
length(intersect(liverBx_12mth,CD8_6mth))
length(intersect(liverBx_12mth,CD4_12mth))
length(intersect(liverBx_12mth,CD8_12mth))
length(intersect(liverBx_12mth,CD4_18mth))
length(intersect(liverBx_12mth,CD8_18mth))
sum(PANORAMA_Pt1[unstimCD4,6])
sum(PANORAMA_Pt1[unstimCD8,6])
sum(PANORAMA_Pt1[CD4.allo,6])
sum(PANORAMA_Pt1[CD8.allo,6])
sum(PANORAMA_Pt1[CD4_6mth,6])
sum(PANORAMA_Pt1[CD8_6mth,6])
sum(PANORAMA_Pt1[CD4_12mth,6])
sum(PANORAMA_Pt1[CD8_12mth,6])
sum(PANORAMA_Pt1[CD4_18mth,6])
sum(PANORAMA_Pt1[CD8_18mth,6])

#liver bx 18 mth
liverBx_18mth=PANORAMA_Pt1[PANORAMA_Pt1[,9]>0,1]
length(liverBx_18mth)
sum(PANORAMA_Pt1[liverBx_18mth,9])
length(intersect(liverBx_18mth,unstimCD4))
length(intersect(liverBx_18mth,unstimCD8))
length(intersect(liverBx_18mth,CD4.allo))
length(intersect(liverBx_18mth,CD8.allo))
length(intersect(liverBx_18mth,CD4_6mth))
length(intersect(liverBx_18mth,CD8_6mth))
length(intersect(liverBx_18mth,CD4_12mth))
length(intersect(liverBx_18mth,CD8_12mth))
length(intersect(liverBx_18mth,CD4_18mth))
length(intersect(liverBx_18mth,CD8_18mth))
sum(PANORAMA_Pt1[unstimCD4,9])
sum(PANORAMA_Pt1[unstimCD8,9])
sum(PANORAMA_Pt1[CD4.allo,9])
sum(PANORAMA_Pt1[CD8.allo,9])
sum(PANORAMA_Pt1[CD4_6mth,9])
sum(PANORAMA_Pt1[CD8_6mth,9])
sum(PANORAMA_Pt1[CD4_12mth,9])
sum(PANORAMA_Pt1[CD8_12mth,9])
sum(PANORAMA_Pt1[CD4_18mth,9])
sum(PANORAMA_Pt1[CD8_18mth,9])

#combine columns for JSD analysis
PANORAMA_Pt1$unstim = PANORAMA_Pt1[,14]+PANORAMA_Pt1[,15]
PANORAMA_Pt1$blood6 = PANORAMA_Pt1[,10]+PANORAMA_Pt1[,11]
PANORAMA_Pt1$blood12 = PANORAMA_Pt1[,4]+PANORAMA_Pt1[,5]
PANORAMA_Pt1$blood18 = PANORAMA_Pt1[,7]+PANORAMA_Pt1[,8]

names(PANORAMA_Pt1)

jsdReport(PANORAMA_Pt1[,c(18,6,9,19,20,21)],topN=1000)
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:
  
  ```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:
  
  ```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

# Z test

# Using unique clones
n1 <-  100961 # DRTCC total unique clones
p1 <- .1125 # Proportion of DRTCC identified in Mo8 CD8

n2 <- 25326 # nonDRTCC total unique clones
p2 <- .0775 # proportion of non DRTCC identified in Mo8 CD8

# Pooled proportion
p_pooled <- (n1 * p1 + n2 * p2) / (n1 + n2)

# Calculate Z-statistic
z_stat <- (p1 - p2) / sqrt(p_pooled * (1 - p_pooled) * (1/n1 + 1/n2))

# Two-tailed test, find critical Z-value for alpha = 0.05
critical_value <- qnorm(1 - 0.05/2)

# Compare Z-statistic to critical value
if (abs(z_stat) > critical_value) {
  cat("Reject the null hypothesis\n")
} else {
  cat("Fail to reject the null hypothesis\n")
}

# Alternatively, calculate the p-value
p_value <- 2 * (1 - pnorm(abs(z_stat)))

# Compare p-value to alpha
if (p_value < 0.05) {
  cat("Reject the null hypothesis\n")
} else {
  cat("Fail to reject the null hypothesis\n")
}

