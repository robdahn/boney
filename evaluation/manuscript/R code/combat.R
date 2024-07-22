#Trial with harmonisation

#July 2024

library(sva)
library(tidyverse)
library(psych)


#prepare a smaller dataset
data_combat <- retest_Boney%>%dplyr::select(OASISID, sROI_bonecortex03, sROI_bonethickness03, 
                                            sROI_headthickness03, tis_head, tismri_volr01, res1)
dat1 <- data_combat[,c(1,7)]
dat2 <- data_combat[,-c(1,7)] #remove the IDs and the resolution (quasi batch)
dat2 <- dat2 %>% t %>% as.matrix() 


batch <- dat1$res1

#modcombat <- model.matrix(~as.factor(res1), data=dat1)
harm_dat <- sva::ComBat(dat2, batch = batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

#get the data back to the long format
harm_dat <- harm_dat %>% t %>% as.data.frame()
harm_dat <- cbind(OASISID, harm_dat, batch)
harm_dat <-  harm_dat %>% dplyr::mutate(main_sBMDH = -1*sROI_bonecortex03)

#reshape data to get the second time points in a wide format
for (i in 1:(length(harm_dat$OASISID))) {
  if (i < length(harm_dat$OASISID) && substr(harm_dat$OASISID[i+1],1,8) == substr(harm_dat$OASISID[i],1,8) ) {
    harm_dat$bonecortex03_t2[i] <- harm_dat$sROI_bonecortex03[i+1] 
    harm_dat$sBMDH_t2[i] <- harm_dat$main_sBMDH[i+1]
    harm_dat$bonethickness03_t2[i] <- harm_dat$sROI_bonethickness03[i+1]
    harm_dat$tis_head_t2[i] <- harm_dat$tis_head[i+1]
    harm_dat$headthickness03_t2[i] <- harm_dat$sROI_headthickness03[i+1]
    harm_dat$tismri_volr01_t2[i] <- harm_dat$tismri_volr01[i+1]#grey matter
  }
  else {
    harm_dat$bonecortex03_t2[i] <- 0
    harm_dat$sBMDH_t2[i] <-  0
    harm_dat$bonethickness03_t2[i] <- 0
    harm_dat$tis_head_t2[i] <- 0
    harm_dat$headthickness03_t2[i] <- 0
    harm_dat$tismri_volr01_t2[i] <- 0
  }
}

# remove double subject rows
retest_harm <-  harm_dat[seq(1, nrow(harm_dat), 2),]


retest_harm <- retest_harm%>%dplyr::select(OASISID, sROI_bonecortex03, bonecortex03_t2, 
                                             main_sBMDH, sBMDH_t2, sROI_bonethickness03,
                                             bonethickness03_t2, tis_head, tis_head_t2, 
                                             sROI_headthickness03, headthickness03_t2,
                                             tismri_volr01, tismri_volr01_t2, batch)#batch is the resolution


#ICC on harmonised data
BMDH.icc <- ICC(retest_harm[c(2:3)]) #Single_random_raters     ICC2 0.88 15 157 157 3.8e-52        0.83        0.91
bone_thick.icc <- ICC(retest_harm[c(6:7)]) #Single_random_raters     ICC2 0.80 9 157 157 5.4e-37        0.74        0.85
IAP.icc <- ICC(retest_harm[c(8:9)])#Single_random_raters     ICC2 0.63 4.6 157 157 4.2e-20        0.53        0.72
TAP.icc <- ICC(retest_harm[c(10:11)])#Single_random_raters     ICC2 0.56 3.6 157 157 4.1e-15        0.44        0.66
GMV.icc <- ICC(retest_harm[c(12:13)])#Single_random_raters     ICC2 0.92 24 157 157 3.3e-66        0.89        0.94

#confirm it's the same for different protocol/resolution
res_12 <- retest_Boney%>%filter(tis_res_vx_vol1 == 1.2)
res_1 <- retest_Boney%>%filter(tis_res_vx_vol1 == 1)
same_protocol <- setdiff(res_12$OASISID, res_1$OASISID)#we get 64 participants with two 1.2 mm res scans

dif_res <- retest_harm[!retest_harm$OASISID %in% same_protocol,]

BMDHdiff.icc <- ICC(retest_harm[c(2:3)])#Single_random_raters     ICC2 0.88 15 157 157 3.8e-52        0.83        0.91
bone_thick_diff.icc <- ICC(retest_harm[c(6:7)])#Single_random_raters     ICC2 0.80 9 157 157 5.4e-37        0.74        0.85
IAP_diff.icc <- ICC(retest_harm[c(8:9)])#Single_random_raters     ICC2 0.63 4.6 157 157 4.2e-20        0.53        0.72
TAP_diff.icc <- ICC(retest_harm[c(10:11)])#Single_random_raters     ICC2 0.56 3.6 157 157 4.1e-15        0.44        0.66
GMV_diff.icc <- ICC(retest_harm[c(12:13)])#Single_random_raters     ICC2 0.92 24 157 157 3.3e-66        0.89        0.94
