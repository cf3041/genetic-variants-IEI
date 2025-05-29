## Variant mapping to structural metrics

library(dplyr)
library(stringr)
library(readxl)
library(car)
library(tidyr)
library(ggplot2)

## Load screen data
sgrna_lib <-readxl::read_xlsx("~/Table S1 Variome_sgRNA_Library.xlsx")
screen_df <-readxl::read_xlsx("~/Table S2 Variome_Screen_Results.xlsx",sheet="All Donors - pAKT_pS6")

## Alpha Missense Pathogenicity
## Load Alpha Missense Data
missense_pik3cd <- read.csv("~/AF-O00329-F1-aa-substitutions.csv") 
missense_pik3r1 <- read.csv("~/AF-P27986-F1-aa-substitutions.csv")

amino_acid_mapping <- c(
  'A' = 'Ala', 'C' = 'Cys', 'D' = 'Asp', 'E' = 'Glu', 'F' = 'Phe',
  'G' = 'Gly', 'H' = 'His', 'I' = 'Ile', 'K' = 'Lys', 'L' = 'Leu',
  'M' = 'Met', 'N' = 'Asn', 'P' = 'Pro', 'Q' = 'Gln', 'R' = 'Arg',
  'S' = 'Ser', 'T' = 'Thr', 'V' = 'Val', 'W' = 'Trp', 'Y' = 'Tyr'
)

convert_mutation <- function(mutation) {
  # Extract original and mutated amino acids
  original_amino_acid <- substr(mutation, 1, 1)
  mutated_amino_acid <- substr(mutation, nchar(mutation), nchar(mutation))
  
  # Convert to three-letter codes
  original_three_letter <- amino_acid_mapping[original_amino_acid]
  mutated_three_letter <- amino_acid_mapping[mutated_amino_acid]
  
  # Combine the three-letter codes with the position in between
  position <- substr(mutation, 2, nchar(mutation) - 1)
  formatted_mutation <- paste0(original_three_letter, position, mutated_three_letter)
  
  return(formatted_mutation)
}
missense_pik3cd$protein.variant <- paste("p.",sapply(missense_pik3cd$protein_variant, convert_mutation),sep="")
missense_pik3r1$protein.variant <- paste("p.",sapply(missense_pik3r1$protein_variant, convert_mutation),sep="")


missense_pik3cd$position <- as.numeric(gsub("[^0-9]", "", missense_pik3cd$protein.variant))
missense_pik3r1$position <- as.numeric(gsub("[^0-9]", "", missense_pik3r1$protein.variant))
missense_pik3cd$pathogenicity.score <- missense_pik3cd$am_pathogenicity
missense_pik3r1$pathogenicity.score <- missense_pik3r1$am_pathogenicity


## Correlation with screen results
screen_temp <- screen_df[which(screen_df$Gene=="PIK3CD"),]
screen_temp <- screen_temp[(grep("Missense",screen_temp$`Mutation category`)),] ## subset to missense mutations

merge_df <- c()
for(i in 1:nrow(screen_temp)){
  row_temp <- screen_temp[i,]
  cat <- strsplit(row_temp$`Mutation category`, ";")[[1]]
  mut <- strsplit(row_temp$`Amino acid edits (3-9 window)`, ";")[[1]]
  
  for(n in 1:length(cat)){
    if(cat[n]=="Missense"){ #filter to only missense mutations
      temp = data.frame(mutation=paste("p.",mut[n],sep=""),sgrna = row_temp$sgrna, LFC=row_temp$LFC,FDR=row_temp$FDR,p.twosided=row_temp$p.twosided,ctrl =row_temp$control_median)
      merge_df <- rbind(merge_df,temp)
    }
  }
}

merge_df <- merge(merge_df,missense_pik3cd,by.x="mutation",by.y="protein.variant")
merge_df <- unique(merge_df)
merge_df <- merge_df[which(merge_df$ctrl > 30),]
merge_df$color <- "Not significant"
merge_df[which(merge_df$p.twosided<0.05),"color"] <- "Significant"


merge_df_filt <- data.frame()
sgrna <- unique(merge_df$sgrna)
for(i in 1:length(sgrna)){
  temp <- merge_df[which(merge_df$sgrna == sgrna[i]),]
  temp <- temp[temp$pathogenicity.score == max(temp$pathogenicity.score), ]
  merge_df_filt <- rbind(merge_df_filt,temp)
}

## Visualize screen results (LFC) vs Alpha Missense Predicted pathogenicity
ggplot(merge_df_filt, aes(x = (LFC), y = pathogenicity.score,color=color)) + scale_color_manual(name = "",values = c("gray","firebrick")) +
  geom_point(size = 2, alpha = 1,aes(text=text)) +  # Scatter points 
  #geom_smooth(method = "lm", se = FALSE, color = "red") +  # Linear regression line
  labs(title = "PIK3CD",
       x = "Screen LFC",
       y = "Alpha Missense Pathogenicity Score",color = "") +
  theme_minimal()


## Bin alpha missense pathogenicity 
merge_df_filt$path_bin <- cut(merge_df_filt$pathogenicity.score, breaks = 10) 

## Levene's Test for homogeneity of variance
leveneTest(merge_df_filt$LFC ~ merge_df_filt$path_bin)[1,3] # p 3.932425e-22




## Structural protein comparisons (G2P data) 

## Load data
struc_pik3cd <- read.delim("~/G2P_PIK3CD_O00329_protein_features.tsv")
struc_pik3r1 <- read.delim("~/G2P_PIK3R1_P27986_protein_features.tsv")


screen_temp <- screen_df[which(screen_df$Gene=="PIK3CD"),]
screen_temp <- screen_temp[(grep("Missense",screen_temp$`Mutation category`)),] #subset to missense mutations

merge_df <- c()
for(i in 1:nrow(screen_temp)){
  row_temp <- screen_temp[i,]
  cat <- strsplit(row_temp$`Mutation category`, ";")[[1]]
  mut <- strsplit(row_temp$`Amino acid edits (3-9 window)`, ";")[[1]]
  
  for(n in 1:length(cat)){
    if(cat[n]=="Missense"){ #filter to only missense mutations
      temp = data.frame(mutation=paste("p.",mut[n],sep=""),sgrna = row_temp$sgrna, LFC=row_temp$LFC,FDR=row_temp$FDR,p.twosided=row_temp$p.twosided,ctrl = row_temp$control_median)
      merge_df <- rbind(merge_df,temp)
    }
  }
}
merge_df$aa.position <- as.numeric(gsub("[^0-9]", "", merge_df$mutation))

merge_df <- merge(merge_df,struc_pik3cd[,c("residueId","AlphaFold.confidence..pLDDT.",
                                           "Accessible.surface.area..Å...",
                                           "Druggability.score..fpocket..", "Hydropathy","Phi.angle..degrees..","Psi.angle..degrees..")],
                  by.x="aa.position",by.y = "residueId")

merge_df <- merge_df[which(merge_df$ctrl > 30),]
merge_df$color <- "Not significant"
merge_df[which(merge_df$p.twosided<0.05),"color"] <- "Significant"


ggplot(merge_df, aes(x = LFC, y = AlphaFold.confidence..pLDDT.,color=color)) + scale_color_manual(name = "",values = c("gray","firebrick")) +
  geom_point( size = 2, alpha = 1) +  # Scatter points
  labs(title = "PIK3CD",
       x = "Screen LFC",
       y = "AlphaFold pLDDT") +
  theme_minimal()

ggplot(merge_df, aes(x = LFC, y = Accessible.surface.area..Å...,color=color)) + scale_color_manual(name = "",values = c("gray","firebrick")) +
  geom_point( size = 2, alpha = 1) +  # Scatter points
  labs(title = "PIK3CD",
       x = "Screen LFC",
       y = "Accessible surface area") +
  theme_minimal()


## Bin pLDDT
merge_df$pLDDT_bin <- cut(merge_df$AlphaFold.confidence..pLDDT., breaks = 10)  # Adjust the number of bins as needed
merge_df$access_bin <- cut(merge_df$Accessible.surface.area..Å..., breaks = 10)  # Adjust the number of bins as needed

## Levene's Test for homogeneity of variance
leveneTest(merge_df$LFC ~ merge_df$pLDDT_bin)[1,3]
leveneTest(merge_df$LFC ~ merge_df$access_bin)[1,3]

