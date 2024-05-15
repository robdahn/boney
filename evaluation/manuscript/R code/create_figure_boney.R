#Data for the supplement
cor_boney_wm <- bones_ukb_WM %>% dplyr::filter(eid %in% old_id$eid)

cor_boney_wm <- cor_boney_wm%>%dplyr::select(
                             eid, sex, age, sit_height, VAT, ASAT, fat_percent, BMI,
                             BMD_head, BMD_total, T_BMD_total, femur_BMD,
                             tismri_volr01, tismri_volr02, tismri_volr03, tismri_TIV,
                             tis_bone, vROI_bonecortex03, sROI_bonecortex03, main_sBMDH,
                             sROI_bonethickness01:sROI_bonethickness09, tis_head,
                             sROI_headthickness01:sROI_headthickness09, waist)

cor_boney_wm <- left_join(cor_boney_wm, ukb_validation, by="eid")

ukb_r1000_WM <- cor_boney_wm[,-1]
ukb_r1000_WM$sex <- as.numeric(ukb_r1000_WM$sex)


r_bones_WM <- psych::corr.test(ukb_r1000_WM, method="spearman", adjust="holm")


write.csv(r_bones_WM$r, "boney_correlation_WM.csv")
write.csv(r_bones_WM$p, "boney_cor_pval_WM.csv")

cor_bones_m <- ukb_r1000_WM%>%filter(sex==1)
cor_bones_m <- cor_bones_m[,-1]
cor_bones_f <- ukb_r1000_WM%>%filter(sex==2)
cor_bones_f <- cor_bones_f[,-1]

#Holm adjusted for multiple comparisons
r_bones_m <- psych::corr.test(cor_bones_m, method = "spearman", adjust = "holm")
r_bones_f <- psych::corr.test(cor_bones_f, method = "spearman", adjust = "holm")

write.csv(r_bones_m$r, "men_cor_Boney_WM.csv")
write.csv(r_bones_f$r, "women_cor_Boney_WM.csv")
write.csv(r_bones_m$p, "men_pval_Boney_WM.csv")
write.csv(r_bones_f$p, "women_pval_Boney_WM.csv")


#Figures
ukb_r1000_WM$sex <- factor(ukb_r1000_WM$sex, levels=c('1','2'))
levels(ukb_r1000_WM$sex) <- list("men"= 1, "women"= 2)

ggplot(data=ukb_r1000_WM, aes(main_sBMDH, BMD_head, colour=sex))+
  geom_point(aes(shape =sex, colour=sex), size = 4, alpha=0.5)+ 
  scale_colour_manual(values=c("#0185ff", "#fb7602"))+
  scale_shape_manual(values = c(15, 16))+
  #geom_smooth(method="lm")+
  theme_classic(base_size = 18)+
  theme(legend.position = c(0.9, 0.90))+
  xlab("Proxy BMD")+
  ylab("Head BMD")+
  scale_x_reverse()


ggplot(data=ukb_r1000_WM, aes(y = main_sBMDH, x = smoking.x, colour=sex))+
  geom_point(aes(shape =sex, colour=sex), size = 3, alpha=0.5)+ 
  scale_colour_manual(values=c("#0185ff", "#fb7602"))+
  scale_shape_manual(values = c(15, 16))+
  geom_smooth(method="loess", se = F)+
  theme_classic(base_size = 14)+
  theme(legend.position = c(0.90, 0.90))+
  xlab("Proportion of smoking")+
  ylab("Proxy BMD")+
  scale_y_reverse()+
  geom_jitter(width = 0.08)

ggplot(data=ukb_r1000_WM, aes(y = main_sBMDH, x = moderate_PA.x, colour=sex))+
  geom_point(aes(shape =sex, colour=sex), size = 3, alpha=0.5)+ 
  scale_colour_manual(values=c("#0185ff", "#fb7602"))+
  scale_shape_manual(values = c(15, 16))+
  geom_smooth(method="loess", se = F)+
  theme_classic(base_size = 14)+
  theme(legend.position = c(0.90, 0.90))+
  xlab("Duration of moderate PA [min]")+
  ylab("Proxy BMD")+
  scale_y_reverse()+
  geom_jitter(width = 0.08)

ggplot(data=ukb_r1000_WM, aes(y = sROI_headthickness03, x = ASAT, colour=sex))+
  geom_point(aes(shape =sex, colour=sex), size = 3, alpha=0.5)+ 
  scale_colour_manual(values=c("#0185ff", "#fb7602"))+
  scale_shape_manual(values = c(15, 16))+
  geom_smooth(method="lm", se = F)+
  theme_classic(base_size = 14)+
  theme(legend.position = c(0.90, 0.90))+
  xlab("Abdominal subcutaneous adipose tissue")+
  ylab("Head thickness (occ.)")