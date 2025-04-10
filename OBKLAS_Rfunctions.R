## ------------------------------------------------------------------------------- ##
#  ---------------------- BeefClass Prediction -----------------------------------  #
#                                                                                   #
# - This code contains functions used in:
#   "Objective carcass grading for bovine animals applying carcass weight, length, 
#   sex and breed information as predictors"
# - Authors (code): Andrew Heggli, Lars Erik Gangsei, Hilde Vinje
# - Date: 17th of November 2020


## 1) summary_data --------------------------------------------------------------- ##
# Coded as a function to make the main script file easier readable
# Input: A data frame, Beefclass_data, containing factor elements 'Sex','CarcassID',
#         'SubSample',and numeric elements 'y_c','y_f','Weight','Length','Age'
# Output: A list with one element for each subsample. Containing number of carcasses
#       , number of observations, and mean and standarddeviation for numeric variables
#         
summary_data <- function(Beefclass_data)
{
  Beefclass_data$Sex_v <- Beefclass_data$Sex
  Beefclass_data$Sex_v[is.na(Beefclass_data$Sex_v)] <- 99
  Beefclass_data$Sex_v <- as.factor(Beefclass_data$Sex_v)
  Subsamples <- as.character(unique(Beefclass_data$SubSample))
  res <- vector('list',length(Subsamples))
  names(res) <- Subsamples
  num_var <- c('y_c','y_f','Weight','Length','Age')
  
  mean_std <- function(df,ind)
  {
    xx <- tapply(df[,ind],list(df$CarcassID,df$Sex_v),mean)
    return(apply(cbind(as.character(round(apply(xx,2,mean,na.rm=TRUE),2)),
                       as.character(round(sqrt(apply(xx,2,var,na.rm=TRUE)),2))),
                 1,paste,collapse=' +/- '))
  }
  
  for(ii in 1:length(res))
  {
    Delta_BC <- Beefclass_data[Beefclass_data$SubSample==Subsamples[ii],]
    res[[ii]] <- cbind(as.character(tapply(Delta_BC$CarcassID,Delta_BC$Sex_v,length)),
                       as.character(rowSums(table(Delta_BC$Sex_v,Delta_BC$CarcassID)>0)))
    for(jj in 1:length(num_var))
    {
      res[[ii]] <- cbind(res[[ii]],mean_std(Delta_BC,
                                            which(names(Delta_BC)==num_var[jj])))
    }
    rownames(res[[ii]]) <- c('Males','Females','Castrats','Unknown Sex')
    colnames(res[[ii]]) <- c('n_obs','n_ind',num_var)
    res[[ii]] <- res[[ii]][!is.na(res[[ii]][,1]),]
    
  }
  return(res)
}


## 2) design_matrix -------------------------------------------------------------- ##
# Innput: A numeric matrix XX (n x (p+1)) whose elements all have value 
#         between 0 and 1, and where each row sums to 1.
# Output: A numeric matrix (n x p). The first column of XX is deleted, and the 
#         p columns returned are equal to difference between the p last columns 
#         of the original matrix and corresponding elements of the first column.
#         
design_matrix <- function(XX)
{
  XX <- XX[,colSums(XX)>0]
  pp <- dim(XX)[2]-1
  return(XX[,-1]-XX[,1]%*%t(as.matrix(rep(1,pp))))
}

## 3) design_matrix based on B-splines ------------------------------------------- ##
# Innput: A numeric vector xx of length n for which splines are evaluated.
#         kk: A knot sequence for quadratic B-splines of length p
# Output: A numeric matrix (n x p). 
#         
design_matrix_BS <- function(xx,kk)
{
  XX <- bs(pmax(kk[1],pmin(kk[length(kk)],xx)),Boundary.knots = 
             kk[c(1,length(kk))],knots = kk,degree=2)
  return(design_matrix(XX))
}


## 4) Confint_REMod ------------------------------------------------------------- ##
# Innput: REObj. Object fittet by lm(y~u) function, where y is a numeric response
#                and u is a categorical predictor.
# Output: A 3 x 4 matrix containing estimates and confidence intervalls for 
#         overall mean (mu), sigma_sq, sigma_u_sq, and interclass correlation 
#         (sigma_u_sq/(sigma_u_sq+sigma_sq)) in the random effects modell defined
#         by y_ij = u_i + e_ij, u_i ~ N(0,sigma_u_sq) and e_ij ~ N(0,sigma_sq)
#         See. Chap. 17.3 in Dean, Voss & Draguljic (2017). Design and Analysis
#         of Experiments. Springer International Publishing AG. Cham, Switzerland.
#         2nd Ed.
#         Also: Montgomery, D. (2017).Design and Analysis of Experiments. John 
#         Wiley & Sons, Inc. Hoboken, NJ. 9th Ed. 

Confint_REMod <- function(REObj,CIlim = 95) 
{
  
  alpha2 <- (100-CIlim)/200 # alpha2 is the significance level divided by two
  
  # Get information from REObj
  MSTr <- anova(REObj)[[3]][1] # MSTr is the treatment  mean sum of squares
  MSE <- anova(REObj)[[3]][2] # MSE is the error mean sum of squares
  df_MSTr <- anova(REObj)[[1]][1] # Degrees of freedom for treatments
  df_MSE <- anova(REObj)[[1]][2] # Degrees of freedom for errors
  NN <- df_MSTr + df_MSE + 1 # Total number of observations 
  
  # Number of replicates per observation:
  r_vec <- as.vector(table(as.factor(REObj$model[[2]])))
  
  # Degrees of freedom associated with correlation ("c" p. 623 in Dean et.al.)
  df_corr <- (sum(r_vec)^2-sum(r_vec^2))/(sum(r_vec)*(df_MSTr)) 
  
  # Degrees of freedom associated estimated sigmaT (Eq. 17.3.14 p.627 in Dean et.al.)
  df_SigmaT <- (MSTr-MSE)^2/((MSTr^2/df_MSTr)+(MSE^2/df_MSE)) 
  
  # Creates a 3x4 matrix which shall be returned
  res <- matrix(NA,3,4) 
  colnames(res) <- c('mu','sigmaSQ','sigmaSQ_u','Corr')  
  rownames(res) <- c('Estimate',paste(as.character(100*c(alpha2,1-alpha2)),
                                      '% CI',sep = '')) 
  
  # Point estimates
  res[1,] <- c(coef(REObj)[1],MSE,(MSTr-MSE)/df_corr,NA) 
  res[1,4] <- res[1,3]/(res[1,3]+res[1,2]) #Intraclass correlation coefficient 
  
  # Confidence interval grand mean. Ref Montgomery Eq. 3.58. If unbalanced data only
  # approximate
  res[2:3,1] <- res[1,1]+qt(c(alpha2,1-alpha2),df_MSE)*sqrt(MSTr/NN) 
  
  # Confidence interval sigma^2 (Eq. 17.3.6 p. 623 in Dean et.al.) 
  res[2:3,2] <- MSE*df_MSE/qchisq(c(alpha2,1-alpha2),df_MSE,lower.tail=FALSE) 
  
  # Approximate confidence interval sigma_T^2 (Eq. 17.3.15 p. 627 in Dean et.al.) 
  res[2:3,3] <- df_SigmaT*res[1,3]/qchisq(c(alpha2,1-alpha2),df_SigmaT,lower.tail=FALSE)
  
  # Confidence interval for correlation. Only approximate if unbalanced data. Ref.
  # Eq.3.57 p. 116 in Montgomery.
  LL <- (1/df_corr)*(((MSTr/MSE)*(1/(qf(alpha2,df_MSTr,df_MSE,lower.tail=FALSE))))-1)
  UU <- (1/df_corr)*(((MSTr/MSE)*(1/(qf(1-alpha2,df_MSTr,df_MSE,lower.tail=FALSE))))-1)
  res[2,4] <- LL/(1+LL) 
  res[3,4] <- UU/(1+UU)
  
  return(res)
}
