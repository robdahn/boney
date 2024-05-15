#variables Juelich for Boney validation
library(dplyr)
library(naniar)
library(janitor)

if (!exists("ukb_validation")){
  
  #import the Juelich data  
  UKBJuelich <- read_tsv('/Users/polona/Downloads/UKB_November_2022/ukb670018_fieldset2_MRI.tsv')
  exclusion_UKB <- read.csv('/Users/polona/Downloads/exclusion_UKB.csv', header = F)
  
  #clean the names and remove the double cases
  colnames(UKBJuelich) <- paste0("v_", colnames(UKBJuelich))
  UKBJuelich <- UKBJuelich%>%dplyr::rename("eid"="v_eid")
  UKBJuelich <- clean_names(UKBJuelich)
  UKBJuelich <- UKBJuelich[!(UKBJuelich$eid %in% exclusion_UKB$V1), ]
  UKBJuelich <- UKBJuelich%>%distinct(eid, .keep_all = T)
  
  #select the lifestyle factors
  ukb_validation <- UKBJuelich%>%dplyr::select(eid, v_874_2_0, v_894_2_0, v_914_2_0, v_20162_2_0, v_1558_2_0)
  names(ukb_validation) <- c("eid", "walk", "moderate_PA", "vigorous_PA", "smoking", "alcohol")
  
  #clean the missing responses
  ukb_validation <- ukb_validation%>% replace_with_na_all(condition = ~.x == c(-3))
  ukb_validation <- ukb_validation%>% replace_with_na_all(condition = ~.x == c(-1))
  
  #recode the alcohol intake (0 - not drink, 5- drink every day )
  ukb_validation <- ukb_validation %>%
    mutate(r_alcohol = 6-alcohol) %>%
    select(-alcohol)
  
  #NA means the person doesn't smoke, set to 0  
  ukb_validation$smoking[is.na(ukb_validation$smoking)] <- 0
  
} else {
  
  #bind the validation values to the ukb_1000
  ukb_1000 <- left_join(ukb_1000, ukb_validation, by="eid")
}

