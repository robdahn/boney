### OASIS-3 Retest case for Boney ###

# July 2024


library(dplyr)
library(stringr)
library(ggplot2)
library(psych)


retest_Boney <- read.csv(paste0("/Volumes/OneTouch/MRData/202301_MRIbones/OASIS3/",
                                "BIDS/sub-OAS30002/ses-d2340/anat/Boney_xm2csv_report_20240427-003047.csv")) 

fnames <- basename(retest_Boney$filename)
subi <- unlist( str_locate_all(fnames[1],'sub-') )
OASISID <- as.data.frame(str_sub(fnames, subi[2]+1, subi[2]+8))
names(OASISID) <- c("OASISID")

retest_Boney <- cbind(OASISID, retest_Boney)

#reshape data to get the second time points in a wide format
for (i in 1:(length(retest_Boney$OASISID))) {
  if (i < length(retest_Boney$OASISID) && substr(retest_Boney$OASISID[i+1],1,8) == substr(retest_Boney$OASISID[i],1,8) ) {
    retest_Boney$bonecortex03_t2[i] <- retest_Boney$sROI_bonecortex03[i+1] 
    retest_Boney$sBMDH_t2[i] <- retest_Boney$main_sBMDH[i+1]
    retest_Boney$bonethickness03_t2[i] <- retest_Boney$sROI_bonethickness03[i+1]
    retest_Boney$tis_head_t2[i] <- retest_Boney$tis_head[i+1]
    retest_Boney$headthickness03_t2[i] <- retest_Boney$sROI_headthickness03[i+1]
    retest_Boney$tismri_volr01_t2[i] <- retest_Boney$tismri_volr01[i+1]#grey matter
  }
  else {
    retest_Boney$bonecortex03_t2[i] <- 0
    retest_Boney$sBMDH_t2[i] <-  0
    retest_Boney$bonethickness03_t2[i] <- 0
    retest_Boney$tis_head_t2[i] <- 0
    retest_Boney$headthickness03_t2[i] <- 0
    retest_Boney$tismri_volr01_t2[i] <- 0
  }
}

retest_Boney <- retest_Boney%>%dplyr::select(OASISID, sROI_bonecortex03, bonecortex03_t2, 
                                       main_sBMDH, sBMDH_t2, sROI_bonethickness03,
                                       bonethickness03_t2, tis_head, tis_head_t2, 
                                       sROI_headthickness03, headthickness03_t2,
                                       tismri_volr01, tismri_volr01_t2, tis_res_vx_vol1)

retest <- filter(retest_Boney, bonecortex03_t2!= 0)

BMDH.icc <- ICC(retest[c(2:3)]) #Single_random_raters     ICC2 0.57 4.0 157 157 1.9e-17        0.40        0.69
bone_thick.icc <- ICC(retest[c(6:7)]) #Single_random_raters     ICC2 0.81 9.8 157 157 1.9e-39        0.75        0.86
IAP.icc <- ICC(retest[c(8:9)])#Single_random_raters     ICC2 0.50 3 157 157 7.7e-12        0.38        0.61
TAP.icc <- ICC(retest[c(10:11)])#Single_random_raters     ICC2 0.16 1.4 157 157 0.021      0.0058        0.31
GMV.icc <- ICC(retest[c(12:13)])#Single_random_raters     ICC2 0.93 31 157 157 3.8e-74        0.91        0.95


#testRetest(retest_try[,c(2)], retest_try[,c(3)], time = NULL, check.keys = F)


#check cases with the same protocol/resolution
res_12 <- retest_Boney%>%filter(tis_res_vx_vol1 == 1.2)
res_1 <- retest_Boney%>%filter(tis_res_vx_vol1 == 1)
same_protocol <- setdiff(res_12$OASISID, res_1$OASISID)#we get 64 participants with two 1.2 mm res scans

same_res <- retest_Boney[retest_Boney$OASISID %in% same_protocol,]
same_res <- filter(same_res, bonecortex03_t2!= 0)

BMDHsame.icc <- ICC(same_res[c(2:3)])#Single_random_raters     ICC2 0.95 45  63  63 9.3e-36        0.91        0.97
bone_thick_same.icc <- ICC(same_res[c(6:7)])#Single_random_raters     ICC2 0.83 11  63  63 5.9e-18        0.74        0.89
IAP_same.icc <- ICC(same_res[c(8:9)])#Single_random_raters     ICC2 0.66 5.7  63  63 3.2e-11        0.43        0.80
TAP_same.icc <- ICC(same_res[c(10:11)])#Single_random_raters     ICC2 0.53 3.6  63  63 4.6e-07        0.32        0.69
GMV_same.icc <- ICC(same_res[c(12:13)])#Single_random_raters     ICC2 0.97 76  63  63 1.0e-42        0.96        0.98

#check cases with different protocols between both time points
diffpro <- retest_Boney[!retest_Boney$OASISID %in% same_protocol,]

diffpro <- filter(diffpro, bonecortex03_t2!= 0)
diffpro <- diffpro %>% select(OASISID, sROI_bonecortex03, bonecortex03_t2, 
                              main_sBMDH, sBMDH_t2, sROI_bonethickness03,
                              bonethickness03_t2, tis_head, tis_head_t2, 
                              sROI_headthickness03, headthickness03_t2,
                              tismri_volr01, tismri_volr01_t2, tis_res_vx_vol1)

BMDHdiff.icc <- ICC(diffpro[c(2:3)])#
bone_thick_diff.icc <- ICC(diffpro[c(6:7)])#Single_fixed_raters      ICC3 0.81 9.4  93  93 2.6e-23        0.72        0.87
IAP_diff.icc <- ICC(diffpro[c(8:9)])#Single_fixed_raters      ICC3 0.32 1.9  93  93 0.00077       0.127        0.49
TAP_diff.icc <- ICC(diffpro[c(10:11)])
GMV_diff.icc <- ICC(diffpro[c(12:13)])#Single_fixed_raters      ICC3 0.89 17  93  93 5.1e-34        0.84        0.93







