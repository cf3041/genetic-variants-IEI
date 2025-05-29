## Variant Level Scoring

library(dplyr)
library(gridExtra)
library(stringr)
library(aod)


## Load data
sgrna_lib <-readxl::read_xlsx("~/Table S1 Variome_sgRNA_Library.xlsx")
screen_df <-readxl::read_xlsx("~/Table S2 Variome_Screen_Results.xlsx",sheet="All Donors - pAKT_pS6")

## Subset to gene of interest, variant scoring only relevant for missense variants
screen_temp <- screen_df[which(screen_df$Gene=="PIK3CD" & screen_df$`Mutation category`%in% c("Missense;Missense;","Missense;Missense;Missense;","Missense;",
                                                                                              "Missense;Missense;Silent;","Missense;Silent;",
                                                                                              "Missense;Silent;Missense;","Silent;Missense;","Silent;Missense;Missense;")),]

## Reformat data
merge_df <- c()
for(i in 1:nrow(screen_temp)){
  row_temp <- screen_temp[i,]
  cat <- strsplit(row_temp$`Mutation category`, ";")[[1]]
  mut <- strsplit(row_temp$`Amino acid edits (3-9 window)`, ";")[[1]]
  
  for(n in 1:length(cat)){
    if(cat[n]=="Missense"){ #filter to only missense mutations
      temp = data.frame(mutation=paste("p.",mut[n],sep=""),sgrna = row_temp$sgrna, LFC=row_temp$LFC,FDR=row_temp$FDR,p.twosided=row_temp$p.twosided,ctrl_count =row_temp$control_median)
      merge_df <- rbind(merge_df,temp)
    }
  }
}

merge_df <- merge_df[which(merge_df$ctrl_count > 30),]

## One hot encoding for base edits by sgRNA
design_matrix <- merge_df[,c("mutation","sgrna","LFC")] %>%
  mutate(present = 1) %>%
  tidyr::spread(mutation, present, fill = 0)

n = ncol(design_matrix) 
qr(as.matrix(design_matrix[,3:n]))$rank ## rank 869 < 929 possible edits suggesting linear dependence


## Fit multiple linear regression model
lm_model <- lm(LFC ~ . + 0, data = design_matrix[,2:n])

lm_results <- data.frame(coef = lm_model$coefficients)
lm_results$edit <- rownames(lm_results)

lm_results$amino <-as.numeric(gsub(".*?([0-9]+).*", "\\1", lm_results$edit))  
lm_results$gene <- "PIK3CD"


## Wald test for significance testing
coef <- coef(lm_model)
coef <- coef[!is.na(coef)]
vcov_matrix <- vcov(lm_model)
vcov_matrix <- vcov_matrix[names(coef),names(coef)]
predictor_names <- names(coef)

p <- c()
names <- c()
for(i in 1:length(predictor_names)){
  term_index <- which(names(coef) == predictor_names[i])  # Get index of the variable
  
  test_result <- wald.test(b = coef, 
                           Sigma = vcov_matrix, 
                           Terms = term_index)
  p <- c(p,test_result[["result"]][["chi2"]][["P"]])
  names <- c(names,predictor_names[i])
}
wald <- data.frame(pval = p, edit = names)
wald$pval_adj <- p.adjust(wald$pval, method = "BH")

lm_results <- merge(lm_results,wald, by = 'edit')





