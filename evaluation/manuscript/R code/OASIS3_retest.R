### OASIS-3 Retest case for Boney ###

# April 2024


library(dplyr)
library(stringr)
library(ggplot2)


OASIS3 <- read.csv("/Volumes/OneTouch/MRData/202301_MRIbones/OASIS3/OASIS3all.csv",
                   header=F)

#######sample selection##############
#Keep MR data 
OASIS3_MR <-  as.data.frame(str_subset(OASIS3$V1,"MR"))
names(OASIS3_MR) <- c("sub")

#Extract day from the file name
OASIS3_MR$day <- as.numeric(str_match(OASIS3_MR$sub, '_d(\\d+)')[, 2])

#Calculate the diff. between the time-points within the subjects
for (i in 1:(length(OASIS3_MR$sub)-1) ) {
  if (substr(OASIS3_MR$sub[i+1],1,8) == substr(OASIS3_MR$sub[i],1,8)) {
    OASIS3_MR$dif[i] <- OASIS3_MR$day[i+1]-OASIS3_MR$day[i]
  } 
  else {
    OASIS3_MR$dif[i] <- 0
  }
} 

#which images have less than 3 months interscan time
OASIS3_MR$retest <- ifelse(OASIS3_MR$dif < 90 & OASIS3_MR$dif > 0, 1, 0)
OASIS3_MR <- OASIS3_MR%>%dplyr::mutate(followup = lag(OASIS3_MR$retest) == 1)

#pick the test-retest images for the total list

OASIS3_retest <- OASIS3_MR%>%dplyr::filter(retest == 1 | followup == T)



#Same for the CT data
OASIS3_CT <-  as.data.frame(str_subset(OASIS3$V1,"CT"))
names(OASIS3_CT) <- c("sub")
OASIS3_CT$day <- as.numeric(str_match(OASIS3_CT$sub, '_d(\\d+)')[, 2])
OASIS3_CT$ID <- substr(OASIS3_CT$sub,1,8)

OASIS3_retest <- left_join(OASIS3_retest, OASIS3_CT, by=c("ID", "day"))

#save the sample
OASIS_retest <- OASIS3_retest%>%dplyr::select(sub.x,sub.y, ID)
write.csv2(OASIS_retest, "OASIS_retest.csv")


#Re-test evaluation 
retest_Boney <- read.csv(paste0("/Volumes/OneTouch/MRData/202301_MRIbones/OASIS3/",
                                "BIDS/sub-OAS30002/ses-d2340/anat/Boney_xm2csv_report_20240427-003047.csv")) 

#/Volumes/OneTouch/MRData/202301_MRIbones/OASIS3/Boney_xm2csv_report_OASIS-20240421-005156.csv (norm_musc)
#Boney_xm2csv_report_20231103-140924.csv old
#Boney_xm2csv_report_20231101-174901.csv #failed fat-suppression correction

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

#check the protocol properties
retest_Boney$res1 <- ifelse(retest_Boney$tis_res_vx_vol1==1, 1, 0)


#check all cases (mixed and matched protocols together)
retest <- filter(retest_Boney, bonecortex03_t2!= 0)

retest <- retest%>%select(OASISID, tis_res_vx_vol1, sROI_bonecortex03, bonecortex03_t2,
                          main_sBMDH, sBMDH_t2, sROI_bonethickness03, bonethickness03_t2,
                          sROI_headthickness03, headthickness03_t2, tis_head, tis_head_t2,
                          tismri_volr01, tismri_volr01_t2, res1)

ggplot(retest, aes(sROI_bonecortex03, bonecortex03_t2))+geom_point()

cor_retest <- retest[,-(1:2)]
psych::corr.test(cor_retest, adjust="holm")


#check cases with the same protocol/resolution
res_12 <- retest_Boney%>%filter(tis_res_vx_vol1 == 1.2)
res_1 <- retest_Boney%>%filter(tis_res_vx_vol1 == 1)
same_protocol <- setdiff(res_12$OASISID, res_1$OASISID)#we get 64 participants with two 1.2 mm res scans

same_res <- retest_Boney[retest_Boney$OASISID %in% same_protocol,]
same_res <- filter(same_res, bonecortex03_t2!= 0)

same_res <- same_res%>%select(OASISID, tis_res_vx_vol1, sROI_bonecortex03, bonecortex03_t2,
                          main_sBMDH, sBMDH_t2, sROI_bonethickness03, bonethickness03_t2,
                          sROI_headthickness03, headthickness03_t2, tis_head, tis_head_t2,
                          tismri_volr01, tismri_volr01_t2)

#remove an outlier due to failed SPM segmentation: OAS31462
same_res <- same_res[-c(63),]
cor_same_res <- same_res[,-(1:2)]
psych::corr.test(cor_same_res, adjust="holm")

#RMSE
sqrt(mean((same_res$sROI_bonecortex03 - same_res$bonecortex03_t2)^2))

ggplot(same_res, aes(sROI_bonecortex03, bonecortex03_t2))+geom_point()+
      geom_smooth(method="lm")+
      xlab("proxy BMD (t1)")+
      ylab("proxy BMD (t2)")+
      theme_bw(base_size = 14)+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(same_res, aes(sROI_bonethickness03, bonethickness03_t2))+geom_point()+
  geom_smooth(method="lm")+
  xlab("Occipital bone thickness (t1)")+
  ylab("Occipital bone thickness (t2)")+
  theme_bw(base_size = 14)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(same_res, aes(sROI_headthickness03, headthickness03_t2))+geom_point()+
  geom_smooth(method="lm")+
  xlab("Head tissue thickness (t1)")+
  ylab("Head tissue thickness (t2)")+
  theme_bw(base_size = 14)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



#check cases with different protocols between both time points
diffpro <- retest_Boney[!retest_Boney$OASISID %in% same_protocol,]

diffpro <- filter(diffpro, bonecortex03_t2!= 0)
diffpro <- diffpro %>% select(OASISID, tis_res_vx_vol1, sROI_bonecortex03, bonecortex03_t2,
                              main_sBMDH, sBMDH_t2, sROI_bonethickness03, bonethickness03_t2,
                              sROI_headthickness03, headthickness03_t2, tis_head, tis_head_t2,
                              tismri_volr01, tismri_volr01_t2, res1)


ggplot(diffpro, aes(sROI_bonethickness03, bonethickness03_t2, colour=as.factor(res1)))+geom_point()
cor_diffpro <- diffpro[,-(1:2)]
psych::corr.test(cor_diffpro, adjust="holm")

# the time-point and protocol are not matched 
# we have to correct this offset by flipping the first with the second measurement
diffpro2 <- diffpro 
for ( coli in seq(from=3, to=length(diffpro2)-2, by=2) )  {
  for ( i in 1:length(diffpro2$OASISID)) {
    if ( diffpro$res1[i]==1 ) { # use diffpro without 2 here!
      tmp <- diffpro2[i,coli]
      diffpro2[i,coli] <- diffpro2[i,coli+1]
      diffpro2[i,coli+1] <- tmp 
      }  
    }
}  


#plot
ggplot(diffpro2, aes(sROI_bonethickness03, bonethickness03_t2, colour=as.factor(res1)))+geom_point()
         
cor_diffpro2 <- diffpro2[,-(1:2)]
psych::corr.test(cor_diffpro2, adjust="holm")


#we should do the same as above for the full sample correlation
retest2 <- retest 
for ( coli in seq(from=3, to=length(retest2)-2, by=2) )  {
  for ( i in 1:length(retest2$OASISID)) {
    if ( retest$res1[i]==1 && !retest$OASISID[i] %in% same_protocol ) { # use retest without 2 here!
      tmp <- retest2[i,coli]
      retest2[i,coli] <- retest2[i,coli+1]
      retest2[i,coli+1] <- tmp 
    }  
  }
} 

cor_retest2 <- retest2[,-(1:2)]
psych::corr.test(cor_retest2, adjust="holm")


#Scanner info
oasis_scanner <- read.csv(paste0("/Volumes/OneTouch/MRData/202301_MRIbones/OASIS3/",
"+tables/OASIS3_data_files/scans/MRI-json-MRI_json_information/resources/csv/files/OASIS3_MR_json.csv"))
oasis_scanner <- oasis_scanner%>%rename(OASISID = subject_id)%>%filter(scan.category=="T1w")



retest_Boney$label <- substr(retest_Boney$P_org, 86, 103)
retest_Boney$label <- str_replace(retest_Boney$label, "sess-", "MR_")
retest_Boney$label <- str_replace(retest_Boney$label, "ses-", "MR_")


scanner <- retest_Boney[retest_Boney$OASISID %in% same_res$OASISID,]
scanner<- left_join(scanner, oasis_scanner, by="label", multiple = "all")


#BMI OASIS retest
#extract the file names for merging data
retest_Boney$OASIS_session_label <-  str_sub(retest_Boney$filename, 58, 75)
retest_Boney$OASIS_session_label <-  str_replace(retest_Boney$OASIS_session_label,"/ses-", "_MR_")

#merge files
small_BMI <- left_join(retest_Boney, OASIS_BMI, by ="OASIS_session_label")

#take every second case
small_BMI <- small_BMI[seq(1, nrow(small_BMI), 2), ]

#correlate head_thickness and BMI
small_BMI <- small_BMI%>%dplyr::select(BPSYS, BPDIAS, HRATE, BMI, c(sROI_headthickness01:sROI_headthickness09),
                                       c(vROI_headthickness01:vROI_headthickness09))

corr.test(small_BMI, method="spearman")

