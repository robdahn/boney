#Bone and Fat Measures Validation (using UKB data)

#April 2024

library(dplyr)
library(party)
library(flexplot)
library(psych)
library(stringr)
library(readr)
library(caret)

out_dir <- "/Volumes/OneTouch/MRData/202301_MRIbones/derivatives/"


rm("spm_boney", "ukb_bone", "ukb_1000", "bones_ukb", 
   "boney_measures", "cor_boney", "ukb_r1000")


#load SPM-preprocessed Boney measures 

spm_boney <- read.csv(paste0("/Volumes/OneTouch/MRData/202301_MRIbones/derivatives/",
                             "boney20240418_ref3_WM/Boney_xml_ref3_WM.csv"), header=T)
fnames <- basename(spm_boney$filename)
subi <- unlist( str_locate_all(fnames[1],'sub-') )
eid_ukb <- as.data.frame(str_sub(fnames, subi[2]+1, subi[2]+7))
names(eid_ukb) <- c("eid")

spm_boney <- cbind(eid_ukb, spm_boney)
spm_boney$eid <- as.numeric(spm_boney$eid)


#load bone measures from UKB
ukb_bone <- read_tsv('/Users/polona/Downloads/UKB_November_2022/ukb670018_fieldset3_MRI20230822.tsv')
exclusion_UKB <- read.csv('/Users/polona/Downloads/exclusion_UKB.csv', header = F) 
ukb_bone <- ukb_bone[!(ukb_bone$eid %in% exclusion_UKB$V1), ] #withdrawal from the study
names(ukb_bone) <- c("eid", "waist", "sit_height", "VAT", "ASAT", "fat_percent", "BMI", 
                     "BMD_head", "BMD_total", "T_BMD_total", "femur_BMD")
ukb_bone <- ukb_bone[(ukb_bone$eid %in% ukb_2000$eid), ]


demographics <- ukb_2000%>%dplyr::select(eid, age, sex)
demographics <- demographics%>%distinct(eid, .keep_all = T)
ukb_bone <- left_join(demographics, ukb_bone, by="eid")

#do the split-half
set.seed(45)
ukb_1000 <- ukb_2000%>%group_by(sex)%>%slice_sample(n=500)#"randomly" pick 1000 subjects 

ukb1000eid <- ukb_1000$eid

#pick 1000 subjects that are not in the previous sample
ukb_explore <- ukb_bone %>% dplyr::filter(!eid %in% ukb1000eid) %>%
                            dplyr::select(eid,BMD_head)

bones_ukb <- left_join(ukb_explore, spm_boney, by = "eid")#merge the file with the extracted measures

boney_measures <- select_if(bones_ukb, is.numeric)  #get only values that are numeric
boney_measures <- boney_measures[,-c(1,nearZeroVar(boney_measures, freqCut = 100/0))]
#R is adding sex as a grouping var, remove it, remove also variables without variance

#exploratory analysis of best predictors for BMD head in the big sample
set.seed(88)
rf_bones<-cforest(BMD_head~., data=boney_measures)

val <- estimates(rf_bones) #best predictor: sROI_bonecortex03 == -(main_sBMDH)

wkdir <- "/Volumes/OneTouch/MRData/202301_MRIbones/derivatives/scripts/"

# exploratory analysis of best predictors for BMD head in 5 folds
source(paste0(wkdir,"best_measure.R"))

#Add validation measures from Juelich
source(paste0(wkdir,"juelich_vars_boney.R"))

#test the selected measure on the other sample
#ukb_1000$eid <- as.character(ukb_1000$eid)
ukb_1000 <- left_join(ukb_1000, spm_boney, by="eid")


#create a dataset with relevant measures & vars
#main_sBMDH is the opposite measure of the sROI_bonecortex03
cor_boney <- ukb_1000%>%select(eid, sex, age, tismri_TIV, tismri_volr01, tismri_volr02,
                               tismri_volr03, classic_bone_med, main_sBMDH, sROI_bonecortex03, 
                               sROI_bonethickness03, BMD_head, BMD_total, T_BMD_total, femur_BMD,
                               sit_height, sROI_headthickness03, BMI, VAT, ASAT,fat_percent, waist,
                               walk, moderate_PA, vigorous_PA, smoking, r_alcohol) 
names(cor_boney) <- c("eid", "sex","age", "TIV", "rGMV", "rWMV", 
                      "rCSFV", "classic_bone_med", "main_sBMDH", "occipital_intensity",
                      "bone thickness_occ", "BMD head", "BMD total", "TBMD total", "BMD femur", 
                      "sit height", "head thickness", "BMI","VAT", "ASAT", "fat %", "waist",
                      "walk", "moderate_PA", "vigorous_PA", "smoking", "alcohol")


#Validate: correlation analysis 
ukb_r1000 <- cor_boney[,-1]

ukb_r1000$sex <- as.numeric(ukb_r1000$sex)


r_bones <- psych::corr.test(ukb_r1000, method="spearman", adjust="holm")


write.csv(r_bones$r, "boney_correlation.csv")
write.csv(r_bones$p, "boney_cor_pval.csv")
write("Probability values (Entries above the diagonal are adjusted for multiple tests.) ", "boney_cor_pval.txt")

#write.csv(r_bones$r, file.path(out_dir, paste0("correlation_", basename(csvfiles[ti])"))

#write.csv(r_bones$r, paste0("correlation_", basename(csvfiles[ti]))
#write.csv(r_bones$p, paste0("pvals_", basename(csvfiles[ti]))

#cor.plot(ukb_r1000,numbers=TRUE,upper=T,diag=TRUE,cex=.7)

#separate analysis for men and women
cor_bones_m <- cor_boney%>%filter(sex==1)
cor_bones_m <- cor_bones_m[,-(1:2)]
cor_bones_f <- cor_boney%>%filter(sex==2)
cor_bones_f <- cor_bones_f[,-(1:2)]

#Holm adjusted for multiple comparisons
r_bones_m <- psych::corr.test(cor_bones_m, method = "spearman", adjust = "holm")
r_bones_f <- psych::corr.test(cor_bones_f, method = "spearman", adjust = "holm")

write.csv(r_bones_m$r, "men_cor_Boney.csv")
write.csv(r_bones_f$r, "women_cor_Boney.csv")
write.csv(r_bones_m$p, "men_pval_Boney.csv")
write.csv(r_bones_f$p, "women_pval_Boney.csv")