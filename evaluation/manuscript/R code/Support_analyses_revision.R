#Support analyses for the revision of the manuscript

#July 2024

library(readr)
library(dplyr)
library(naniar)
library(janitor)



#import the Juelich data  
UKBJuelich <- read_tsv('/Users/polona/Downloads/UKB_November_2022/ukb670018_fieldset3_MRI20230822.tsv')
exclusion_UKB <- read.csv('/Users/polona/Downloads/exclusion_UKB.csv', header = F)

#clean the names
colnames(UKBJuelich) <- paste0("v_", colnames(UKBJuelich))
UKBJuelich <- UKBJuelich%>%dplyr::rename("eid"="v_eid")
UKBJuelich <- clean_names(UKBJuelich)
UKBJuelich <- UKBJuelich[!(UKBJuelich$eid %in% exclusion_UKB$V1), ]
#UKBJuelich <- left_join(UKBJuelich, age_proper, by="eid", keep = F, multiple = "first")

#get variables of interest
ukb_validation <- UKBJuelich%>%dplyr::select(1, 312, 2:12, 15, 21:24, 43:45, 48:49, 51:67)


#clean the missing responses
ukb_validation <- ukb_validation%>% replace_with_na_all(condition = ~.x == c(-3))
ukb_validation <- ukb_validation%>% replace_with_na_all(condition = ~.x == c(-1))

#recode the alcohol intake (0 - not drink, 5- drink every day )
ukb_validation <- ukb_validation %>%
  mutate(r_alcohol=6-v_1558_2_0) %>%
  select(-v_1558_2_0)

#NA means the person doesn't smoke, set to 0  
ukb_validation$v_20162_2_0[is.na(ukb_validation$v_20162_2_0)] <- 0


#calculate WHR, WHtR
ukb_validation <- ukb_validation%>%dplyr::mutate(
                                          WHR = v_48_2_0/v_49_2_0, # waist to hip ratio
                                          WHtR = v_48_2_0/v_50_2_0, # waist to height ratio
                                          SitHt_Ht = v_51_2_0/v_50_2_0) # sitting height to standing height ratio


#load SPM-preprocessed Boney measures 
spm_boney <- read.csv(paste0("/Volumes/OneTouch/MRData/202301_MRIbones/derivatives/",
                             "boney20240418_ref3_WM/Boney_xml_ref3_WM.csv"), header=T)
fnames <- basename(spm_boney$filename)
subi <- unlist( str_locate_all(fnames[1],'sub-') )
eid_ukb <- as.data.frame(str_sub(fnames, subi[2]+1, subi[2]+7))
names(eid_ukb) <- c("eid")

spm_boney <- cbind(eid_ukb, spm_boney)
spm_boney$eid <- as.numeric(spm_boney$eid)


#merge the UKB data with the Boney measures (old sample)
df_validation <- left_join(ukb_1000, ukb_validation, by="eid")
df_validation <- left_join(df_validation, spm_boney, by="eid")

cor_boney <- df_validation%>%select(eid, v_31_0_0, age, tismri_TIV, tismri_volr01, tismri_volr02,
                                tismri_volr03, classic_bone_med, main_sBMDH, sROI_bonecortex01:sROI_bonecortex09,
                                vROI_bonecortex03,sROI_bonethickness01:sROI_bonethickness09, v_23226_2_0, v_23239_2_0,
                                v_23300_2_0, v_51_2_0, v_50_2_0, sROI_headthickness01:sROI_headthickness09, tis_head,
                                v_21001_2_0, v_22407_2_0, v_22408_2_0, v_23099_2_0, WHR, WHtR,
                                v_874_2_0, v_894_2_0, v_914_2_0, v_20162_2_0, r_alcohol) 

#names(cor_boney) <- c("eid", "sex","age", "TIV", "rGMV", "rWMV", 
#                      "rCSFV", "classic_bone_med", "main_sBMDH", "occipital_intensity",
#                      "bone thickness_occ", "BMD head", "BMD total", "BMD femur", 
#                      "sit height", "height", "TAP", "IAP", "BMI","VAT", "ASAT", "fat %", "WHR",
#                      "WHtR", "walk", "moderate_PA", "vigorous_PA", "smoking", "alcohol")


#Validate: correlation analysis 
ukb_r1000 <- cor_boney[,-1]

r_bones <- psych::corr.test(ukb_r1000, method="spearman", adjust="holm")

write.csv(r_bones$r, "boney_correlation_revision.csv")
write.csv(r_bones$p, "boney_cor_pval_revision.csv")

#separate analysis for men and women
cor_bones_m <- cor_boney%>%filter(v_31_0_0==1)
cor_bones_m <- cor_bones_m[,-(1:2)]
cor_bones_f <- cor_boney%>%filter(v_31_0_0==0)
cor_bones_f <- cor_bones_f[,-(1:2)]

#Holm adjusted for multiple comparisons
r_bones_m <- psych::corr.test(cor_bones_m, method = "spearman", adjust = "holm")
r_bones_f <- psych::corr.test(cor_bones_f, method = "spearman", adjust = "holm")

write.csv(r_bones_m$r, "men_cor_Boney_revision.csv")
write.csv(r_bones_f$r, "women_cor_Boney_revision.csv")
write.csv(r_bones_m$p, "men_pval_Boney_revision.csv")
write.csv(r_bones_f$p, "women_pval_Boney_revision.csv")











#prepare the validation sample dataset for logistic regression
falls <- ukb_validation%>%dplyr::select(eid, v_3005_3_0)#v3005-3.0 UKB "fractures resulting from simple fall at the measurement point 3"
old_eid <- left_join(old_eid, falls, by="eid")

sum(old_eid$v_3005_3_0, na.rm = T) #only 8 people reported a fall and a fracture from it
old_eid$v_3005_3_0[is.na(old_eid$v_3005_3_0)] <- 0

fract <- old_eid%>%dplyr::select(eid, sex, age, main_sBMDH, v_3005_3_0, `BMD head`, `BMD femur`)
fract$sex <- ifelse(fract$sex == 2, 0, 1)#recode women to 0

#run the logistic model
fract_model <- glm(v_3005_3_0 ~ sex + age + main_sBMDH,family=binomial(link='logit'),data=fract)
summary(fract_model)# no parameter is statistically significantly predicting a fracture
#Call:
#glm(formula = v_3005_3_0 ~ sex + age + main_sBMDH, family = binomial(link = "logit"), 
#    data = fract)

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)
#(Intercept) -2.68808    3.38650  -0.794    0.427
#sex         -1.18789    0.84045  -1.413    0.158
#age         -0.01035    0.05656  -0.183    0.855
#main_sBMDH   1.24659    2.10696   0.592    0.554

#(Dispersion parameter for binomial family taken to be 1)

#Null deviance: 93.189  on 999  degrees of freedom
#Residual deviance: 90.551  on 996  degrees of freedom
#AIC: 98.551

#Number of Fisher Scoring iterations: 8
#



#Association with serum markers
#library(flexplot)
#flexplot(main_sBMDH~v_30710_0_0, data= ukb_validation) #CRP
#flexplot(main_sBMDH~v_30680_0_0, data= ukb_validation) #Calcium
#flexplot(main_sBMDH~v_30610_0_0, data= ukb_validation) #Alkaline phosphatase
#flexplot(main_sBMDH~v_30770_0_0, data= ukb_validation) #IGF-1

