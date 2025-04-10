## ------------------------------------------------------------------------------- ##
#  ---------------------- BeefClass Prediction -----------------------------------  #
#                                                                                   #
# - This code produces results/ figures etc used in the paper:
#   "Objective carcass grading for bovine animals applying carcass weight, length, 
#   sex and breed information as predictors". 
# - Authors (code): Andrew Heggli, Lars Erik Gangsei, Hilde Vinje
# - Date: 26th of January 2021

## 1) Clean workspace, install packages, functions and load data ----------------- ##
rm(list=ls())

## 2) Install packages ----------------------------------------------------------- ##
packages <- c('splines','nlme','mixlm','xtable')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

for(pp in packages){library(pp,character.only=TRUE)}

## 3) Load data and functions ---------------------------------------------------- ##
## a) Load data - file name: 20200916_Beefclass_data.csv 
Beefclass_data <- read.csv2(choose.files("",
                                         caption = 'choose the file "20200916_Beefclass_data.csv"'))

## b) Run functions - file name Heggli_Gangsei_Vinje_RFunctions_BeefClass.R
source(choose.files("",
                    caption = 'choose the file "OBKLAS_Rfunctions.R"'))


## 4) Descriptive data used in Table1. ------------------------------------------- ##

## a) y_c: EUROP conformity as average of fore, middle and back parts.
Beefclass_data$y_c <- (Beefclass_data$EUROPConformBack + 
                         Beefclass_data$EUROPConformMiddle +
                         Beefclass_data$EUROPConformFront)/3

Beefclass_data$y_c[is.na(Beefclass_data$y_c)] <- 
  Beefclass_data$EUROPConformTot[is.na(Beefclass_data$y_c)]

## b) y_f: EUROPE Fat Cover.
Beefclass_data$y_f <- Beefclass_data$EUROPFat

## c) See file "Heggli_Gangsei_Vinje_RFunctions_BeefClass.R" for "summary_data()"
Table1 <- summary_data(Beefclass_data)

## 5) Create response variables and predictor variables in dataframe "Input_data" ##

## a) y_c and y_f(responses): see over. 
##    A Total registration is applied if registrations for parts are missing.
Input_data <- data.frame(y_c = Beefclass_data$y_c,
                         y_f = Beefclass_data$y_f)


## b) x_k (predictor): KFactor. 
Input_data$x_k <- Beefclass_data$Weight/((Beefclass_data$Length/100)^3)

## c) X_B (predictor): Breed defined as 9 different groups.
BreedGrops <- list(Group1 = c('00','01','02','03','04','05','06','07','08',
                              as.character(c(10:12,14:17,91,99))),
                   Group2 = as.character(c(13,27,28,29,31,33,35,39,40,90,98)),
                   Group3 ='09',Group4 = '21',Group5 = '22',Group6 = '23',
                   Group7 = c('24','30'),Group8 = '25',Group9 = '26')

# Create NA matrix with columns for each breed group
Input_data$BreedHR <- matrix(NA,dim(Beefclass_data)[1],length(BreedGrops))
colnames(Input_data$BreedHR) <- names(BreedGrops)

# Fill matrix with percentage (0.0 - 1.0) of heritage from each group
# Rows are individuals, Columns are groups. All rows sum to 1.0
# i.e. an individual could be: 0.125 Group 1, 0.125 Group 3, and 0.75 Group 8

for(ii in 1:length(BreedGrops))
{
  ci <- which(is.element(names(Beefclass_data),
                         paste('Breed',BreedGrops[[ii]],sep=''))) 
  Input_data$BreedHR[,ii] <- rowSums(as.matrix(Beefclass_data[,ci]))/100 
}

# See file "Heggli_Gangsei_Vinje_RFunctions_BeefClass.R" for "design_matrix()" 
Input_data$X_b <- design_matrix(Input_data$BreedHR)


## e) Z_a, Z_w and Z_l (predictors): B-spline bases from on ln(age), weight and 
##    length
Knots <- list(Age = log(c(100,300,730,1460,10*365)),
              Weight = c(100,275,450),
              Length = c(150,190,230))

# See file "Heggli_Gangsei_Vinje_RFunctions_BeefClass.R" for "design_matrix_BS()" 
Input_data$Z_a <- I(design_matrix_BS(log(Beefclass_data$Age),Knots$Age))
Input_data$Z_w <- I(design_matrix_BS(Beefclass_data$Weight,Knots$Weight))
Input_data$Z_l <- I(design_matrix_BS(Beefclass_data$Length,Knots$Length))

## f) CarcassID and Sex, variables used in analysis.
Input_data$CarcassID <- Beefclass_data$CarcassID
Input_data$Sex <- Beefclass_data$Sex


## 6) Split "Input_data" into training and test data ---------------------------- ##
Input_data <- list(Train = Input_data[Beefclass_data$SubSample=='Train',],
                   Test1 = Input_data[Beefclass_data$SubSample=='Test1',],
                   Test2 = Input_data[Beefclass_data$SubSample=='Test2',])

## 7) Model selection, best subset, AIC ------------------------------------------ ##

## a) Define variables that might or might not be part of the model
Pred_names <- c('X_b','Z_w','Z_l','Z_a')

## b) Matrix for storing results
Pred_var <- cbind(rep(0:1,each=8),rep(rep(0:1,each=4),2),
                  rep(rep(0:1,each=2),4),rep(0:1,8),matrix(NA,16,8))

colnames(Pred_var) <- c(Pred_names,'AIC_BC','AIC_CC','AIC_BF','AIC_CF','RMSE_BC','RMSE_CC','RMSE_BF','RMSE_CF')

Pred_var

## c) Loop over input combinations. Calculate AIC and RMSE for models crossed over 
##    EUROP-classes (conformity and fat cover), and the two sexes. 

for (ii in 1:dim(Pred_var)[1]){
  Mod_form_C <- paste('y_c ~ ', paste(c('x_k',Pred_names[Pred_var[ii,1:4]==1]),
                                      collapse = '+'))
  Pred_var[ii,5] <- AIC(lm(Mod_form_C, 
                           data = Input_data$Train[Input_data$Train$Sex==0,]))
  Pred_var[ii,9] <- summary(lm(Mod_form_C, 
                               data = Input_data$Train[Input_data$Train$Sex==0,]))$sigma
  Pred_var[ii,6] <- AIC(lm(Mod_form_C, 
                           data = Input_data$Train[Input_data$Train$Sex==1,]))
  Pred_var[ii,10] <- summary(lm(Mod_form_C, 
                                data = Input_data$Train[Input_data$Train$Sex==1,]))$sigma
  
  Mod_form_F <- paste('y_f ~ ', paste(c('x_k',Pred_names[Pred_var[ii,1:4]==1]),
                                      collapse = '+'))
  Pred_var[ii,7] <- AIC(lm(Mod_form_F, 
                           data = Input_data$Train[Input_data$Train$Sex==0,]))
  Pred_var[ii,11] <- summary(lm(Mod_form_F, 
                                data = Input_data$Train[Input_data$Train$Sex==0,]))$sigma
  Pred_var[ii,8] <- AIC(lm(Mod_form_F, 
                           data = Input_data$Train[Input_data$Train$Sex==1,]))
  Pred_var[ii,12] <- summary(lm(Mod_form_F, 
                                data = Input_data$Train[Input_data$Train$Sex==1,]))$sigma
}

#By running the following commands the output will show the model with the lowest AIC and RMSE:

#Pred_var
#For the full model:    
# apply(Pred_var,2,min)

#   X_b          Z_w          Z_l          Z_a       AIC_BC       AIC_CC 
#0.0000000    0.0000000    0.0000000    0.0000000 3762.2480086 2647.2900962 
#   AIC_BF       AIC_CF      RMSE_BC      RMSE_CC      RMSE_BF      RMSE_CF 
#6576.3523754 5535.4547314  0.5800948    0.5865486    1.1187845    1.5541281 


#For the reduced model: 
# apply(Pred_var[1:8,],2,min)


## 8) Fit prediction models ------------------------------------------------------ ##

## a) Model form for the full model (including all predictors). 

Mod_full <- formula( ~ x_k + X_b + Z_w + Z_l + Z_a)

## b) Model formulas for the different combinations of EUROP classes crossed over
##    sex and model complexity (reduced or full). See calculations in Pred_var for
##    model selection. Note that Z_l (B-spline base on length) is omitted for 
##    combinations "Conformity_males" (full) and "Fat_females_Red"

Mod_forms <- list(Conformity_Males = update(Mod_full,'y_c~.-Z_l'),
                  Conformity_Females = update(Mod_full,'y_c~.'),
                  Fat_Males = update(Mod_full,'y_f~.'),
                  Fat_Females = update(Mod_full,'y_f~.'),
                  Conformity_Males_Red = update(Mod_full,'y_c~ . -X_b'),
                  Conformity_Females_Red = update(Mod_full,'y_c~ . -X_b'),
                  Fat_Males_Red = update(Mod_full,'y_f~ . -X_b - Z_l'),
                  Fat_Females_Red = update(Mod_full,'y_f~ . -X_b'))

## c) Fit the models
Models <- vector('list',length(Mod_forms)) 
names(Models)  <- names(Mod_forms)
for(mm in 1:length(Mod_forms))
{
  t_sex <- ifelse(is.element(mm,seq(1,7,by=2)),0,1)
  Models[[mm]] <- lm(Mod_forms[[mm]],
                     data = Input_data$Train[Input_data$Train$Sex==t_sex,])
}


## d) Get summary statistics, for models, i.e. Estimated error variance (MSE), 
##    R-squared values and AIC values
SumStats <- NULL
for(ii in 1:length(Models))
{
  SumStats <- cbind(SumStats,c(summary(Models[[ii]])$sigma^2,
                               summary(Models[[ii]])$r.squared,
                               AIC(Models[[ii]])))
}

rownames(SumStats) <- c('MSE','R2','AIC')
colnames(SumStats) <- names(Models)

## 9) Test the model, i.e., based on Eq.2 in Heggli et.al. (2020) ---------------- ##

## a) Construct predicted and adjusted values in the two test sets (loop over jj)
for(jj in 1:3)
{
  Input_data[[jj]]$Y_hat <- I(matrix(NA,dim(Input_data[[jj]])[1],length(Models)/2))
  colnames(Input_data[[jj]]$Y_hat) <- gsub('_Males','',names(Models)[seq(1,8,by=2)])
  Input_data[[jj]]$Y_adj <- Input_data[[jj]]$Y_hat
  
  for(ii in 1:4)
  {
    Input_data[[jj]]$Y_hat[Input_data[[jj]]$Sex == 0,ii] <- 
      predict(Models[[(ii-1)*2+1]],newdata = Input_data[[jj]][Input_data[[jj]]$Sex == 0,])
    Input_data[[jj]]$Y_hat[Input_data[[jj]]$Sex == 1,ii] <- 
      predict(Models[[ii*2]],newdata = Input_data[[jj]][Input_data[[jj]]$Sex == 1,])
    if(is.element(ii,c(1,3)))
    {
      Input_data[[jj]]$Y_adj[,ii] <- Input_data[[jj]]$y_c - Input_data[[jj]]$Y_hat[,ii]
    }else{
      Input_data[[jj]]$Y_adj[,ii] <- Input_data[[jj]]$y_f - Input_data[[jj]]$Y_hat[,ii]
    }
  }
}


## RMSE per BreedGroup
Input_data$Train$SQ_Res <- (cbind(Input_data$Train$y_c,
                                  Input_data$Train$y_f,
                                  Input_data$Train$y_c,
                                  Input_data$Train$y_f)-
                              Input_data$Train$Y_hat)^2

RMSE_group <- matrix(NA,dim(Input_data$Train$BreedHR)[2],5)
colnames(RMSE_group) <- c(colnames(Input_data$Train$Y_hat),'n')
RMSE_group[,5] <- colSums(Input_data$Train$BreedHR) 

for(ii in 1:4)
{
  RMSE_group[,ii] <- (colSums(matrix(Input_data$Train$SQ_Res[,ii],
                                     dim(Input_data$Train$BreedHR)[1],
                                     dim(Input_data$Train$BreedHR)[2],
                                     byrow=FALSE)*Input_data$Train$BreedHR)/
                        colSums(Input_data$Train$BreedHR))
}

# Run the following command to show the root mean squared error per breed group: 
# RMSE_group


## b) Make Array with Estimates and Confidence Intervals
Table2 <- array(NA,dim=c(2,16,4))
dimnames(Table2)[[3]] <- colnames(Input_data[[2]]$Y_adj)
dimnames(Table2)[[2]] <- c(paste(rep(c('mu','sigmaSQ_u','sigmaSQ_t','Corr'),each=3),
                                 rep(c('Estimate','2.5%','97.5%'),4)),
                           'EU - % point', 'EU - bias', 'EU - slope', 'RMSEP')
dimnames(Table2)[[1]] <- c('Test1','Test2')

EU_point_func <- function(diff,Penalty){return(Penalty[min(c(abs(diff)+1,5))])}

for(jj in 1:2)
{
  for(ii in 1:4)
  {
    Table2[jj,1:12,ii] <- as.vector(Confint_REMod(lm(Y_adj[,ii]~as.factor(CarcassID),data = 
                                                       Input_data[[jj+1]]),CIlim = 95))
    if(is.element(ii,c(1,3)))
    {
      EU_median <- round(tapply(round(Input_data[[jj+1]]$y_c),
                                Input_data[[jj+1]]$CarcassID,median))
      Penalty <- c(10,6,-9,-27,-48)
    }else{
      EU_median <- round(tapply(round(Input_data[[jj+1]]$y_f),
                                Input_data[[jj+1]]$CarcassID,median))
      Penalty <- c(10,9,0,-13,-30)
    }
    EU_Auto <- round(tapply(Input_data[[jj+1]]$Y_hat[,ii],
                            Input_data[[jj+1]]$CarcassID,median))
    
    Table2[jj,13,ii] <- (100*sum(apply(EU_median-EU_Auto,1,EU_point_func,Penalty = Penalty))/
                           (10*length(EU_median)))
    
    Table2[jj,14,ii] <- mean(EU_median-EU_Auto)
    
    Table2[jj,15,ii] <- coef(lm(EU_median~EU_Auto))[2]
    
    Table2[jj,16,ii] <- sqrt(mean(((Input_data[[jj+1]]$Y_adj[,ii])^2)))
    
  }
}


## 10) Figure and pearson correlation for observed vs predicted values ----------- ##
## a) Create a matrix for storing Pearson correlations
main_vec <- c('Conformity full','Fat full','Conformity reduced','Fat reduced')
Corr_mat <- matrix(NA,2,4)
colnames(Corr_mat) <- main_vec
rownames(Corr_mat) <- c('National','International')

## b) Figure adjustments
dev.off()
axislim <- c(1,15)
par(mfrow=c(2,2))

## c) Make the plots, and calculate correlations

for(ii in c(1,3,2,4))
{
  for(jj in 1:2)
  {
    if(is.element(ii,c(1,3)))
    {
      y_vec <-Input_data[[jj+1]]$y_c
      MGroups <- c('P','O','R','U','E')
    }else{
      y_vec <-Input_data[[jj+1]]$y_f
      MGroups <- as.character(1:5)
    }
    Corr_mat[jj,ii] <- cor(y_vec,Input_data[[jj+1]]$Y_hat[,ii],
                           use = 'pairwise.complete.obs',method = 'pearson') 
    if(jj==1)
    {
      plot(Input_data$Test1$Y_hat[,ii],y_vec,
           ylim = axislim,xlim = axislim, axes = FALSE,main = main_vec[ii],
           xlab = 'Predicted', ylab = 'Observed',pch='.',cex=5,col='gray') 
    }else{
      points(Input_data$Test2$Y_hat[,ii],y_vec,pch=17) 
    }
  }
  curve(1*x,add=TRUE,col='black',lwd =1.5, lty=3)
  axis(1, at = seq(2,15, by = 3), labels = MGroups, 
       cex.axis= 0.8, tick = FALSE)
  axis(1, at = 1:15,labels = rep(c('-','','+'),5), tick = FALSE)
  axis(1, at = seq(3.5,15, by = 3), labels = rep('',4))
  axis(2, at = seq(2,15, by = 3), labels = MGroups, 
       cex.axis= 0.8, tick = FALSE,las = 2)
  axis(2, at = 1:15,labels = rep(c('-','','+'),5),las = 2, tick = FALSE)
  axis(2, at = seq(3.5,15, by = 3), labels = rep('',4))
  box()
  abline(h = seq(3.5,15, by = 3),col = 'grey', lty = 3)
  abline(v = seq(3.5,15, by = 3),col = 'grey', lty = 3)
  text(x=4,y=14,labels = paste(c('N:','I:'),as.character(round(Corr_mat[,ii],2)),
                               sep='',collapse=' and '))
}


## 11) View results used in the article ------------------------------------------ ##

## a) Table 1
# > print(Table1$Train)
#         n_obs  n_ind  y_c             y_f             Weight            
# Males   "2140" "2140" "6.56 +/- 2.2"  "6.9 +/- 1.54"  "323.04 +/- 79.65"
# Females "1482" "1482" "4.59 +/- 2.12" "7.67 +/- 2.63" "263.16 +/- 71.79"
#          Length             Age                 
# Males   "204.38 +/- 14.1"  "541.33 +/- 173.83" 
# Females "207.06 +/- 17.65" "1399.79 +/- 958.66"

# > print(Table1$Test1)
#         n_obs n_ind y_c             y_f             Weight            
# Males   "187" "86"  "7.58 +/- 2.25" "7.41 +/- 1.71" "341.61 +/- 64.22"
# Females "76"  "38"  "5.19 +/- 1.94" "8.14 +/- 1.87" "302.23 +/- 53.63"
#         Length             Age                  
# Males   "202.48 +/- 11.73" "523.21 +/- 96.39"   
# Females "211.39 +/- 10.5"  "2149.32 +/- 1172.64"

# > print(Table1$Test2)
#         n_obs n_ind y_c             y_f             Weight            
# Males   "35"  "7"   "6.63 +/- 1.26" "7.83 +/- 2.12" "316.29 +/- 68.73"
# Females "15"  "3"   "3.33 +/- 0.5"  "8.4 +/- 2.03"  "276.7 +/- 34.47" 
#         Length             Age                 
# Males   "203.14 +/- 12.27" "565 +/- 56.82"     
# Females "219.67 +/- 10.41" "1346.33 +/- 465.07"

## b) Summary statistics from Training data used in text in article
# > print(round(SumStats,digits=2))
#       Conformity_Males Conformity_Females Fat_Males Fat_Females Conformity_Males_Red
# MSE             0.34               0.34      1.25        2.42                 0.59
# R2              0.93               0.92      0.47        0.65                 0.88
# AIC          3762.25            2647.29   6576.35     5535.45              4961.40
#     Conformity_Females_Red Fat_Males_Red Fat_Females_Red
# MSE             0.46          1.70            3.57
# R2              0.90          0.28            0.49
# AIC          3054.39       7218.15         6105.81

## c) Table 2
# > print(round((Table2[1,,]), digits = 4))
#                    Conformity   Fat   Conformity_Red Fat_Red
# mu Estimate            0.1056 -0.1779        -0.0504  0.0057
# mu 2.5%                0.0035 -0.4268        -0.2006 -0.3079
# mu 97.5%               0.2077  0.0710         0.0998  0.3194
# sigmaSQ_u Estimate     0.0985  0.1643         0.0985  0.1643
# sigmaSQ_u 2.5%         0.0789  0.1316         0.0789  0.1316
# sigmaSQ_u 97.5%        0.1264  0.2109         0.1264  0.2109
# sigmaSQ_t Estimate     0.2843  1.8884         0.6692  3.0440
# sigmaSQ_t 2.5%         0.2164  1.4792         0.5208  2.3928
# sigmaSQ_t 97.5%        0.3903  2.4955         0.8919  4.0041
# Corr Estimate          0.7428  0.9200         0.8718  0.9488
# Corr 2.5%              0.6567  0.8892         0.8242  0.9287
# Corr 97.5%             0.8106  0.9427         0.9075  0.9635
# EU - % point          69.2742 53.9516        54.7581 31.3710
# EU - bias              0.1210 -0.2097        -0.0645  0.0403
# EU - slope             1.0130  0.7833         0.9886  0.4300
# RMSEP                  0.6232  1.4410         0.8739  1.7843


# > print(round((Table2[2,,]), digits = 4))
#                    Conformity    Fat  Conformity_Red Fat_Red
# mu Estimate            0.3097  0.4075         0.1654  0.9419
# mu 2.5%                0.0556 -0.0926        -0.1879  0.1524
# mu 97.5%               0.5638  0.9077         0.5187  1.7313
# sigmaSQ_u Estimate     0.1900  0.5000         0.1900  0.5000
# sigmaSQ_u 2.5%         0.1281  0.3370         0.1281  0.3370
# sigmaSQ_u 97.5%        0.3111  0.8186         0.3111  0.8186
# sigmaSQ_t Estimate     0.1200  0.5124         0.2676  1.4258
# sigmaSQ_t 2.5%         0.0472  0.2159         0.1163  0.6467
# sigmaSQ_t 97.5%        0.6997  2.3732         1.1277  5.3186
# Corr Estimate          0.3872  0.5061         0.5848  0.7404
# Corr 2.5%              0.1222  0.2305         0.3132  0.5109
# Corr 97.5%             0.7309  0.8037         0.8447  0.9130
# EU - % point          88.0000 84.0000        84.0000 38.0000
# EU - bias              0.1000  0.4000         0.0000  1.0000
# EU - slope             0.9012  1.4674         0.8841  1.5833
# RMSEP                  0.6246  1.0570         0.6741  1.6311

## d) Correlations between automatic and human classification (test sets):
# > print(round(Corr_mat,digits=2))
#               Conformity full Fat full Conformity reduced Fat reduced
# National                 0.97     0.61               0.94        0.30
# International            0.96     0.91               0.94        0.87

## e) RMSE per breed group based on training data set
# > print(round(RMSE_group, digits = 4))
# Conformity    Fat  Conformity_Red Fat_Red      n
#     0.3097  1.5804         0.4432  1.8306   2135.3125
#     0.4124  1.8789         0.5470  2.8946   148.1875
#     0.3033  1.4881         0.5872  1.8575   72.8125
#     0.3316  1.9405         0.3922  3.8054   218.5625
#     0.4201  1.9087         0.7155  3.0355   459.0625
#     0.3967  1.9397         0.5395  3.8252   122.3125
#     0.3419  1.8043         0.9164  3.2633   335.5625
#     0.3278  1.8306         0.5271  2.7481   94.5000
#     0.3518  2.9045         0.4896  9.4548   35.6875
