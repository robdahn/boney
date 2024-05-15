#bootstrap of the correlation coefficient between the two versions of the measures
library(ggplot2)
#library(CarletonStats)
#library(boot)
library(bootcorci)



#test normality
shapiro.test(ukb_1000$proxyBMD)

ggplot(ukb_1000, aes(proxyBMD)) +
  geom_histogram()


#######bootstrap
set.seed(123)
b_BMD <- bootCor(ukb_3000$BMD_head, ukb_3000$BMD_total, conf.level = 0.95, B = 10000,
                 plot.hist = TRUE)


###alternative
install.packages("cocor")
library(cocor)
cocor.dep.groups.overlap( 0.75, 0.68, 0.98, 1000, 
                          alternative = "greater", 
                          test = "all", 
                          alpha = 0.05, 
                          conf.level = 0.95, 
                          null.value = 0 
                          )
#ref: 10.1371/journal.pone.0121945

###bootstrap for Spearman correlation coef.

cor_BMD <- corci(ukb_1000$BMD_head, ukb_1000$proxyBMD, method = "spearman")






# not working
#x <- ukb_1000$BMD_head
#y <- ukb_1000$proxyBMD
#dat <- data.frame(x,y)

#set.seed(88)

#b_BMD <- boot::boot(dat, 
#                    statistic = function(dat, i) {
#                      cor(dat[1,2], method='pearson') 
#                    },
#                    R = 1000
#)

#boot.ci(b_BMD, type = c("norm", "basic", "perc", "bca"))

#plot(density(b_BMD$t))
#abline(v = 0, lty = "dashed", col = "grey60")