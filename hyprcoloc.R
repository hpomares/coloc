#Colocalization with HyPrColoc from http://127.0.0.1:30129/library/hyprcoloc/doc/hyprcoloc.html
library(MendelianRandomization)
library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(tidyverse)
library(rio)
library(hyprcoloc)

# load
setwd("/Users/hugopomaresmillan/Desktop/mediation/DIET_MR")

# T2D # Mahajan
t2dadjbmi_exp_data <- extract_instruments(outcomes='ebi-a-GCST007516')# already clumps 58

# DIET. Download files from https://www.thessgac.org/data
carb_out_data <- read_outcome_data(
    snps = t2dadjbmi_exp_data$rsid,
    filename = "Diet_Carbohydrate_GWAS_MA_SSGAC_2020_MolPsych.txt",
    sep = "\t",
    snp_col = "rsID",
    beta_col = "Beta",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "FREQA1_HRC",
    pval_col = "Pval",
    samplesize_col = "N")

sugar_out_data <- read_outcome_data(
    snps = t2dadjbmi_exp_data$rsid,
    filename = "Diet_Sugar_GWAS_MA_SSGAC_2020_MolPsych.txt",
    sep = "\t",
    snp_col = "rsID",
    beta_col = "Beta",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "FREQA1_HRC",
    pval_col = "Pval",
    samplesize_col = "N")

fat_out_data <- read_outcome_data(
    snps = t2dadjbmi_exp_data$rsid,
    filename = "Diet_Fat_GWAS_MA_SSGAC_2020_MolPsych.txt",
    sep = "\t",
    snp_col = "rsID",
    beta_col = "Beta",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "FREQA1_HRC",
    pval_col = "Pval",
    samplesize_col = "N")

prot_out_data <- read_outcome_data(
    snps = t2dadjbmi_exp_data$rsid,
    filename = "Diet_Protein_GWAS_MA_SSGAC_2020_MolPsych.txt",
    sep = "\t",
    snp_col = "rsID",
    beta_col = "Beta",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "FREQA1_HRC",
    pval_col = "Pval",
    samplesize_col = "N")

# harmonise to forward strand
MergedDF_1 <- harmonise_data(
    exposure_dat = t2dadjbmi_exp_data, 
    outcome_dat = c(carb_out_data,sugar_out_data, fat_out_data, prot_out_data))
dim(MergedDF_1) #57

# format for coloc
names(MergedDF_1)
betas <- MergedDF_1[,c(1,6,7,27,40,53)]
rownames(betas) <- betas[,1]
betas[,1] <- NULL
names(betas)[1] <- "T1"
names(betas)[2] <- "T2"
names(betas)[3] <- "T3"
names(betas)[4] <- "T4"
names(betas)[5] <- "T5"
betas <- as.matrix(betas)
#
ses <- MergedDF_1[,c(1,61,16,28,41,54)]
rownames(ses) <- ses[,1]
ses[,1] <- NULL
names(ses)[1] <- "T1"
names(ses)[2] <- "T2"
names(ses)[3] <- "T3"
names(ses)[4] <- "T4"
names(ses)[5] <- "T5"
ses <- as.matrix(ses)
#
binary.traits = c(1,0,0,0,0);
traits <- paste0("T", 1:dim(betas)[2]);
rsid <- rownames(betas);
res <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid, snpscores = TRUE);
res #rs1558902
sink("res_df_score.txt")           
res$snpscores[[2]]
sink() 

res_binary <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid, binary.outcomes = binary.traits);
res_binary

# Assessing sensitivity to  prior configuration parameters 
prior.options = c(1e-4, 1e-10, 1e-15, 1e-20, 1e-25, 1e-100);
for(i in prior.options){
    res_binary <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid, 
                            uniform.priors = TRUE, prior.1 = i, reg.steps = 2);
    print(paste0("prior.1 = ",i));
    print(res_binary);
}

#  prior.1 = 1e-10
prior.1 = 1e-10;
prior.c.options = c(0.05, 0.02, 0.01, 0.005);
for(i in prior.c.options){
    res_binary <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid,
                            uniform.priors = FALSE, prior.1 = prior.1, prior.c = i);
    print(c(paste0("prior.1 = ",prior.1), paste0("prior.c = ",i)));
    print(res_binary);
}

# BB algorithm off
res_off <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid, bb.alg = FALSE);
res_off 
res_off_df<-res_off$results

# Plot
library(locuscomparer)
MergedDF_1_plot<-MergedDF_1
names(MergedDF_1_plot)[names(MergedDF_1_plot) == "SNP"] <- "rsid"
names(MergedDF_1_plot)[names(MergedDF_1_plot) == "pval.exposure"] <- "pval"

MergedDF_1_filtered_16<-filter(MergedDF_1_plot, chr.exposure == 16)
locuscompare(in_fn1 = MergedDF_1_filtered_16[,c(1,63)], 
             in_fn2 = MergedDF_1_filtered_16[,c(1,63)],
             title1 = 'T2D adjusted BMI GWAS',
             title2 = 'Carbohydrate-Sugar (E%) intake GWAS')
#
MergedDF_1_filtered_10<-filter(MergedDF_1_plot, chr.exposure == 10)
locuscompare(in_fn1 = MergedDF_1_filtered_10[,c(1,5)], 
             in_fn2 = MergedDF_1_filtered_10[,c(1,5)],
             title1 = 'T2D adjusted BMI GWAS',
             title2 = 'Carbohydrate-Sugar (E%) intake GWAS')

######################### ######################### ######################### 
######################### Pairwise: T2D-CARB #########################
######################### ######################### ######################### 
MergedDF_t2dadjbmi_carb <- harmonise_data(
    exposure_dat = t2dadjbmi_exp_data, 
    outcome_dat = carb_out_data)
dim(MergedDF_t2dadjbmi_carb)
#coloc
names(MergedDF_t2dadjbmi_carb)
betas_t2d_carb <- MergedDF_t2dadjbmi_carb[,c(1,6,7)]
rownames(betas_t2d_carb) <- betas_t2d_carb[,1]
betas_t2d_carb[,1] <- NULL
names(betas_t2d_carb)[1] <- "T1"
names(betas_t2d_carb)[2] <- "T2"
betas_t2d_carb <- as.matrix(betas_t2d_carb)
#
ses_t2d_carb <- MergedDF_t2dadjbmi_carb[,c(1,22,16)]
rownames(ses_t2d_carb) <- ses_t2d_carb[,1]
ses_t2d_carb[,1] <- NULL
names(ses_t2d_carb)[1] <- "T1"
names(ses_t2d_carb)[2] <- "T2"
ses_t2d_carb <- as.matrix(ses_t2d_carb)
#
traits_t2d_carb <- paste0("T", 1:dim(betas_t2d_carb)[2]);
rsid_t2d_carb <- rownames(betas_t2d_carb);
res_pairwise_t2d_carb <- hyprcoloc(betas_t2d_carb, ses_t2d_carb, trait.names=traits_t2d_carb, snp.id=rsid_t2d_carb,snpscores = TRUE);
res_pairwise_t2d_carb #
#
prior.options = c(1e-4, 1e-10, 1e-15, 1e-20, 1e-25, 1e-100);
for(i in prior.options){
    res_pairwise_t2d_carb <- hyprcoloc(betas_t2d_carb, ses_t2d_carb, trait.names=traits_t2d_carb, snp.id=rsid_t2d_carb, 
                                         uniform.priors = TRUE, prior.1 = i, reg.steps = 2);
    print(paste0("prior.1 = ",i));
    print(res_pairwise_t2d_carb);
}
######################### Pairwise: T2D-SUGAR #########################
######################### ######################### ######################### 
MergedDF_t2dadjbmi_sugar <- harmonise_data(
    exposure_dat = t2dadjbmi_exp_data, 
    outcome_dat = sugar_out_data)
dim(MergedDF_t2dadjbmi_sugar)
#coloc
names(MergedDF_t2dadjbmi_sugar)
betas_t2d_sugar <- MergedDF_t2dadjbmi_sugar[,c(1,6,7)]
rownames(betas_t2d_sugar) <- betas_t2d_sugar[,1]
betas_t2d_sugar[,1] <- NULL
names(betas_t2d_sugar)[1] <- "T1"
names(betas_t2d_sugar)[2] <- "T2"
betas_t2d_sugar <- as.matrix(betas_t2d_sugar)
#
ses_t2d_sugar <- MergedDF_t2dadjbmi_sugar[,c(1,22,16)]
rownames(ses_t2d_sugar) <- ses_t2d_sugar[,1]
ses_t2d_sugar[,1] <- NULL
names(ses_t2d_sugar)[1] <- "T1"
names(ses_t2d_sugar)[2] <- "T2"
ses_t2d_sugar <- as.matrix(ses_t2d_sugar)
#
traits_t2d_sugar <- paste0("T", 1:dim(betas_t2d_sugar)[2]);
rsid_t2d_sugar <- rownames(betas_t2d_sugar);
res_pairwise_t2d_sugar <- hyprcoloc(betas_t2d_sugar, ses_t2d_sugar, trait.names=traits_t2d_sugar, snp.id=rsid_t2d_sugar,snpscores = TRUE);
res_pairwise_t2d_sugar #
#
prior.options = c(1e-4, 1e-10, 1e-15, 1e-20, 1e-25, 1e-100);
for(i in prior.options){
    res_pairwise_t2d_sugar <- hyprcoloc(betas_t2d_sugar, ses_t2d_sugar, trait.names=traits_t2d_sugar, snp.id=rsid_t2d_sugar, 
                                       uniform.priors = TRUE, prior.1 = i, reg.steps = 2);
    print(paste0("prior.1 = ",i));
    print(res_pairwise_t2d_sugar);
}
######################### Pairwise: T2D-prot #########################
######################### ######################### ######################### 
MergedDF_t2dadjbmi_prot <- harmonise_data(
    exposure_dat = t2dadjbmi_exp_data, 
    outcome_dat = prot_out_data)
dim(MergedDF_t2dadjbmi_prot)
#coloc
names(MergedDF_t2dadjbmi_prot)
betas_t2d_prot <- MergedDF_t2dadjbmi_prot[,c(1,6,7)]
rownames(betas_t2d_prot) <- betas_t2d_prot[,1]
betas_t2d_prot[,1] <- NULL
names(betas_t2d_prot)[1] <- "T1"
names(betas_t2d_prot)[2] <- "T2"
betas_t2d_prot <- as.matrix(betas_t2d_prot)
#
ses_t2d_prot <- MergedDF_t2dadjbmi_prot[,c(1,22,16)]
rownames(ses_t2d_prot) <- ses_t2d_prot[,1]
ses_t2d_prot[,1] <- NULL
names(ses_t2d_prot)[1] <- "T1"
names(ses_t2d_prot)[2] <- "T2"
ses_t2d_prot <- as.matrix(ses_t2d_prot)
#
traits_t2d_prot <- paste0("T", 1:dim(betas_t2d_prot)[2]);
rsid_t2d_prot <- rownames(betas_t2d_prot);
res_pairwise_t2d_prot <- hyprcoloc(betas_t2d_prot, ses_t2d_prot, trait.names=traits_t2d_prot, snp.id=rsid_t2d_prot,snpscores = TRUE);
res_pairwise_t2d_prot #
#
prior.options = c(1e-4, 1e-10, 1e-15, 1e-20, 1e-25, 1e-100);
for(i in prior.options){
    res_pairwise_t2d_prot <- hyprcoloc(betas_t2d_prot, ses_t2d_prot, trait.names=traits_t2d_prot, snp.id=rsid_t2d_prot, 
                                        uniform.priors = TRUE, prior.1 = i, reg.steps = 2);
    print(paste0("prior.1 = ",i));
    print(res_pairwise_t2d_prot);
}
######################### Pairwise: T2D-fat #########################
######################### ######################### ######################### 
MergedDF_t2dadjbmi_fat <- harmonise_data(
    exposure_dat = t2dadjbmi_exp_data, 
    outcome_dat = fat_out_data)
dim(MergedDF_t2dadjbmi_fat)
#coloc
names(MergedDF_t2dadjbmi_fat)
betas_t2d_fat <- MergedDF_t2dadjbmi_fat[,c(1,6,7)]
rownames(betas_t2d_fat) <- betas_t2d_fat[,1]
betas_t2d_fat[,1] <- NULL
names(betas_t2d_fat)[1] <- "T1"
names(betas_t2d_fat)[2] <- "T2"
betas_t2d_fat <- as.matrix(betas_t2d_fat)
#
ses_t2d_fat <- MergedDF_t2dadjbmi_fat[,c(1,22,16)]
rownames(ses_t2d_fat) <- ses_t2d_fat[,1]
ses_t2d_fat[,1] <- NULL
names(ses_t2d_fat)[1] <- "T1"
names(ses_t2d_fat)[2] <- "T2"
ses_t2d_fat <- as.matrix(ses_t2d_fat)
#
traits_t2d_fat <- paste0("T", 1:dim(betas_t2d_fat)[2]);
rsid_t2d_fat <- rownames(betas_t2d_fat);
res_pairwise_t2d_fat <- hyprcoloc(betas_t2d_fat, ses_t2d_fat, trait.names=traits_t2d_fat, snp.id=rsid_t2d_fat,snpscores = TRUE);
res_pairwise_t2d_fat #
#
prior.options = c(1e-4, 1e-10, 1e-15, 1e-20, 1e-25, 1e-100);
for(i in prior.options){
    res_pairwise_t2d_fat <- hyprcoloc(betas_t2d_fat, ses_t2d_fat, trait.names=traits_t2d_fat, snp.id=rsid_t2d_fat, 
                                       uniform.priors = TRUE, prior.1 = i, reg.steps = 2);
    print(paste0("prior.1 = ",i));
    print(res_pairwise_t2d_fat);
}
############### END #################################
