#cross validation of exploratory best measure
library(caret)
library(party)
library(flexplot)


if (exists("bones_ukb")) {
  
  bones_ukb$eid <- as.numeric(bones_ukb$eid)
  bones_ukb <- select_if(bones_ukb, is.numeric)  #get only numeric vars
  bones_ukb<- bones_ukb[,-c(nearZeroVar(bones_ukb, freqCut = 100/0))]#remove vars with zero variance
  
  
  set.seed(88)
  folds <- createFolds(y=bones_ukb$eid, k = 5)
  
  #apply random forest to all folds
  for (i in 1:length(folds)) {
    rfi <- cforest(BMD_head~., data = bones_ukb [folds[[i]],]) 
    
    # save cforest model
    #rf  <- paste("rf_bones_", i, sep = "")
    #assign(rf, rfi) 
    
    val <- paste("val_", i, sep = "")
    print(assign(val, flexplot::estimates(rfi) )$importance[1:3]) #print the top 3 most important predictors 
    #flexplot calculates the variable importance of each predictor in the RF regression
    #defined as RMSE of predicted versus permuted observations
    
    
  }
} else {
  print("Missing bones_ukb data frame.")
}



