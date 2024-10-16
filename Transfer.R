rm(list = ls())  

set.seed(133)

load("Wheat_1.RData")

library(glmnet)
library(ggplot2)

Partitions <- 10
nfold<-10
datalabel <- "Wheat_1"
Proxy_index <-2
Gold_index <-1
Pheno<- dat_ls$Pheno
Geno<- dat_ls$Geno
Env_names = unique(Pheno$Env)
Trait_names=colnames(Pheno)[3:ncol(Pheno)]
transferlabel<-paste0(Env_names[Proxy_index],"to",Env_names[Gold_index])
Y_proxydata<-Pheno[Pheno$Env==Env_names[Proxy_index],]
Y_golddata<-Pheno[Pheno$Env==Env_names[Gold_index],]
Y_proxydata<-Y_proxydata[Y_proxydata$GID %in% intersect(Y_proxydata$GID, Y_golddata$GID),]
Y_golddata<-Y_golddata[Y_golddata$GID %in% Y_proxydata$GID,]
Geno<-Geno[intersect(Y_proxydata$GID, Y_golddata$GID),intersect(Y_proxydata$GID, Y_golddata$GID)]


#### Variable and data preparation section####
###SVD of Geno
X=svd(Geno)
U=X$u
d=X$d
Q_var=quantile(d,probs=0.00)
Q_var
Pos_Q=which(d>Q_var)
U_red=U[,Pos_Q]
d_red=d[Pos_Q]
D=diag(sqrt(d_red))
xD1=U_red%*%D
X=cbind(1, xD1)
dim(X)

n=dim(X)[1]
p=dim(xD1)[2]


#I_1 is the identity matrix but with 0 at the a_00, sometimes referred to as D 
I_1<- diag(c(0, rep(1, p)))

#### Functions section ####
mse <- function(observed, predicted) mean((observed - predicted)^2)

Betahat=function(lam, X,Y){
  #using svd to calculate the inverse
  A<-(t(X)%*%X)+lam*I_1
  
  svd_A <- tryCatch(svd(A), error = function(e) NULL)
  if (is.null(svd_A)) {
    message("Error in SVD computation.")
    return(NULL)
  }
  U <- svd_A$u
  Sigma <- svd_A$d
  V <- svd_A$v
  
  Sigma_inv <- diag(1 / Sigma)
  
  # Compute the inverse of A using the formula A^-1 = V * Sigma^-1 * U^T
  A_inv <- V %*% Sigma_inv %*% t(U)
  return(A_inv%*%(t(X)%*%Y))}

Delta=function(lam, X,Y, Betaproxy){
  #Transforming Y into Y*
  Y_star=Y-X%*%Betaproxy
  #using svd to calculate the inverse
  A<-(t(X)%*%X)+lam*I_1
  svd_A <- tryCatch(svd(A), error = function(e) NULL)
  
  if (is.null(svd_A)) {
    message("Error in SVD computation.")
    return(NULL)
  }
  U <- svd_A$u
  Sigma <- svd_A$d
  V <- svd_A$v
  Sigma_inv <- diag(1 / Sigma)
  # Compute pseudoinverse of A, using A^-1 = V * Sigma^-1 * U^T
  A_inv <- V %*% Sigma_inv %*% t(U)
  return(A_inv%*%(t(X)%*%Y_star))}

Prediction_ARR = function(Xtst, Beta){
  Xtst%*%Beta }

Prediction_Trans = function(Xtst, Delta, Betaproxy){
  Xtst%*%(Betaproxy+Delta) }
#Dataframes for saving the results

newResults <- data.frame()
newPredictions <- data.frame()
newCoefs <- data.frame()


#Generating the lambda grid (100)

R2v =  exp(seq(log(0.01),log(0.99),length=100))
psxx = mean(apply(X[,-1],1,function(x)sum(x^2)))
plambv = (1-R2v)/(R2v/psxx)#R2v 

### for loop for each trait

for (t in 1:length(Trait_names)) {
  Y_proxy<-Y_proxydata[,t+2]
  Y_gold<-Y_golddata[,t+2]
  trait<-Trait_names[t]
  

  #### Proxy section (Transfer_RR) #### 
  
  #Notice that since we are not partioning the data, we dont have predictions nor results, only coefficients
  pRidge <- cv.glmnet(
    x = X[,-1],
    y = Y_proxy,
    family = "gaussian",
    alpha = 0,
    nfolds = nfold,
    lambda = plambv,
    type.measure = 'mse',
    standardize = FALSE
  )
  
  betaproxy<-coef(pRidge, s = "lambda.min")
  
  newbetacoef_row<- data.frame(Class="Proxy",Env=Env_names[Proxy_index], coefType="Beta_Glmnet", Trait=trait, Partition= "All", as.matrix(t(betaproxy)))
  colnames(newbetacoef_row)<-c("Class","Env", "coefType", "Trait", "Partition", 0:p)
  newCoefs<-rbind(newCoefs, newbetacoef_row)
  
  
  
  ## We do the 10 partitions of X and Y_gold 
  for (i in 1:Partitions) {
    set.seed(i)
    training_proportion <- 0.8
    #For Gold
    training_indices <- sample(n, n * training_proportion)
    X_tr <-X[training_indices,]
    Y_tr <- Y_gold[training_indices ]
    X_tst<-X[-training_indices, ]
    Y_tst <- Y_gold[-training_indices ]
    
    sxx = mean(apply(X_tr[,-1],1,function(x)sum(x^2)))
    lambv = (1-R2v)/(R2v/sxx)#R2v 
    
    
    
    n_p=length(training_indices)
    Grpv_k = findInterval(cut(sample(1:n_p, n_p), breaks=nfold), 1:n_p)
    #pGrpv_k = findInterval(cut(sample(1:n_p, n_p), breaks=nfold), 1:n_p)
    
    #Scaling step
    scalingfactor<-  mean(Y_tr)/mean(Y_proxy) 
    Y_tr<-Y_tr/scalingfactor
    
    
    #### Gold section (Transfer_RR) #####
    #vector for saving the averages MSEs for each lambda
    tavgMSElambda<-c()
    
    #Matrix for saving innerMSEs
    tinnerMSEs<-matrix(NA, nrow = 100, ncol = nfold)
    #For each lambda we perform 10-fold cross validation and obtain and average MSE
    
    for (k in 1:nfold)  {
      innertraining_indices<- which(Grpv_k!=k)
      X_innertraining <-X_tr[innertraining_indices,]
      Y_innertraining <- Y_tr[innertraining_indices ]
      X_innertesting<- X_tr[-innertraining_indices, ]
      Y_innertesting <- Y_tr[-innertraining_indices ]
      for (l in 1:length(lambv)){
        lambda=lambv[l]
        
        A<-(t(X_innertraining)%*%X_innertraining)+lambda*I_1
        
        svd_A <- tryCatch(svd(A), error = function(e) NULL)
        if (is.null(svd_A)) {
          message("Error in SVD computation. Skipping this iteration")
          next
        }
        U <- svd_A$u
        Sigma <- svd_A$d
        V <- svd_A$v
        
        Sigma_inv <- diag(1 / Sigma)
        
        # Compute the pseudoinverse of A using the formula A^-1 = V * Sigma^-1 * U^T
        A_inv <- V %*% Sigma_inv %*% t(U)
        
        #Transfer part
        
        Y_star <- Y_innertraining-X_innertraining%*%betaproxy
        delta<-A_inv%*%(t(X_innertraining)%*%Y_star)
        
        itpred<-Prediction_Trans(X_innertesting, delta, betaproxy)
        tMSE_inner<-mse(Y_innertesting, itpred)
        tinnerMSEs[l,k]<-tMSE_inner
        cat(transferlabel, trait, i, k, l,  "tMSE_inner:", tMSE_inner, "\n")
      }
    }
    for (l in 1:length(lambv)){
      tavgMSElambda[l]<- mean(tinnerMSEs[l,], na.rm=TRUE)
    }
    
    #Now we have 100 avgMSE (one for each lambda) so we can obtain lambda.min and train using this
    tlambda.min=lambv[which.min(tavgMSElambda)]
    
    
    ##Transfer Training ##
    
    Delta_final=Delta(tlambda.min, X_tr, Y_tr, betaproxy)
    
    newdeltacoef_row<- data.frame(Class="Gold", Env=Env_names[Gold_index], coefType="Delta", Trait=trait, Partition= i, t(as.matrix(Delta_final)))
    colnames(newdeltacoef_row)<-c("Class","Env", "coefType", "Trait", "Partition", 0:p)
    newCoefs<-rbind(newCoefs, newdeltacoef_row)
    
    transferprediction<-Prediction_Trans(X_tst,Delta_final, betaproxy)
    #Scaling back step
    transferprediction<- transferprediction*scalingfactor
    
    #We prepare transfer prediction and add it to allPredictions
    tnewprediction_row <- data.frame(Env=Env_names[Gold_index], Trait=trait, predType="Transfer_RR", Partition=i, t(as.matrix(transferprediction)))
    
    newPredictions <- rbind(newPredictions, tnewprediction_row,tnewprediction_row)
    
    # We evaluate and report on newResults
    
    tNRMSE<-(sqrt(mse(Y_tst, transferprediction)))/mean(Y_tst)
    tCor <- cor(Y_tst, as.matrix(transferprediction))
    
    cat(i, "tNRMSE:", tNRMSE, "\n")
    
    tnewresult_row <- data.frame(
      Class="Gold",
      Env=Env_names[Gold_index],
      Trait=trait,
      Model="Transfer_RR",
      Partition = i, 
      NRMSE = tNRMSE, 
      Cor = tCor, 
      Lambda.min = tlambda.min)
    colnames(tnewresult_row)=c("Class", "Env", "Trait", "Model", "Partition", "NRMSE", "Cor", "Lambda.min")
    newResults <- rbind(newResults, tnewresult_row)
  }
  
  
  
  #### Proxy section (Transfer_ARR) #### 
  
  #vector for saving the averages MSEs for each lambda
  pavgMSElambda<-c()
  #Matrix for saving innerMSEs
  pinnerMSEs<-matrix(NA, nrow = 100, ncol = nfold)
  
  #For each lambda we perform 10-fold cross validation and obtain and average MSE

  pGrpv_k = findInterval(cut(sample(1:n, n), breaks=nfold), 1:n)
  
  for (k_p in 1:nfold)  {
    
    innertraining_indices<- which(pGrpv_k!=k_p)
    X_innertraining_a <-X[innertraining_indices,]
    Y_innertraining_a <- Y_proxy[innertraining_indices ]
    X_innertesting_a<-X[-innertraining_indices, ]
    Y_innertesting_a <- Y_proxy[-innertraining_indices ]
    for (l in 1:length(plambv)){
      lambda=plambv[l]
      betahat <- tryCatch(Betahat(lambda, X_innertraining_a, Y_innertraining_a), error = function(e) NULL)
      if (is.null(betahat)) {
        message("Error in Betahat computation. Skipping this iteration.")
        next
      }
      pinnerpred<-Prediction_ARR(X_innertesting_a, betahat)
      aMSE_inner<-mse(Y_innertesting, pinnerpred)
      pinnerMSEs[l,k_p]<-aMSE_inner
      cat(trait, k_p, l,  "pinnerMSE:", aMSE_inner, "\n")
    }
  }
  for (l in 1:length(plambv)){
    pavgMSElambda[l]<- mean(pinnerMSEs[l,], na.rm = TRUE)}
  
  #Now we have 100 MSE (one for each lambda) so we can obtain lambda.min and retrain using this
  plambda.min=plambv[which.min(pavgMSElambda)]
  #This is outside the kfold cv, hence we have one betaproxy per trait
  betaproxy_a=Betahat(plambda.min, X, Y_proxy)
  #We save betaproxy 
  newbetacoef_row<- data.frame(Class="Proxy", Env=Env_names[Proxy_index], coefType="Beta_ARR", Trait=trait, Partition= "All", t(as.matrix(betaproxy_a)))
  colnames(newbetacoef_row)<-c("Class","Env", "coefType", "Trait", "Partition", 0:p)
  newCoefs <- rbind(newCoefs, newbetacoef_row)

  
  for (i in 1:Partitions) {
    set.seed(i)
    training_proportion <- 0.8
    #For Gold
    training_indices <- sample(n, n * training_proportion)
    X_tr_a <-X[training_indices,]
    Y_tr_a <- Y_gold[training_indices ]
    X_tst_a<-X[-training_indices, ]
    Y_tst_a <- Y_gold[-training_indices ]
    #Generating the lambda grid (100)
    
    R2v =  exp(seq(log(0.01),log(0.99),length=100))
    sxx = mean(apply(X_tr_a[,-1],1,function(x)sum(x^2)))
    lambv = (1-R2v)/(R2v/sxx)#R2v 
    
  
    n_p=length(training_indices)
    Grpv_k = findInterval(cut(sample(1:n_p, n_p), breaks=nfold), 1:n_p)
    
    #Scaling step
    scalingfactor<-  mean(Y_tr_a)/mean(Y_proxy) 
    Y_tr_a<-Y_tr_a/scalingfactor
      
    #### Gold  section (Transfer_ARR) #####
    #vector for saving the averages MSEs for each lambda

    tavgMSElambda_a<-c()
    
    #Matrices for saving innerMSEs

    tinnerMSEs_a<-matrix(NA, nrow = 100, ncol = nfold)
    #For each lambda we perform 10-fold cross validation and obtain and average MSE
    
    for (k in 1:nfold)  {
      innertraining_indices<- which(Grpv_k!=k)
      X_innertraining_a <-X_tr_a[innertraining_indices,]
      Y_innertraining_a <- Y_tr_a[innertraining_indices ]
      X_innertesting_a<- X_tr_a[-innertraining_indices, ]
      Y_innertesting_a <- Y_tr_a[-innertraining_indices ]
      for (l in 1:length(lambv)){
        lambda=lambv[l]

        A<-(t(X_innertraining_a)%*%X_innertraining_a)+lambda*I_1
        
        # Compute the inverse using solve
        A_inv <- tryCatch(solve(A), error = function(e) NULL)
        if (is.null(A_inv)) {
          message("Error in solve computation. Skipping this iteration")
          next
        }
        #Transfer part
        
        Y_star_a <- Y_innertraining_a-X_innertraining_a%*%betaproxy_a
        delta_a<-A_inv%*%(t(X_innertraining_a)%*%Y_star_a)
        
        itpred_a<-Prediction_Trans(X_innertesting_a, delta_a, betaproxy_a)
        tMSE_inner_a<-mse(Y_innertesting_a, itpred_a)
        tinnerMSEs_a[l,k]<-tMSE_inner_a
        cat(transferlabel, trait, i, k, l,  "tMSE_inner:", tMSE_inner_a, "\n")
      }
    }
    for (l in 1:length(lambv)){
      tavgMSElambda_a[l]<- mean(tinnerMSEs_a[l,], na.rm=TRUE)
    }
    
    #Now we have 100 avgMSE (one for each lambda) so we can obtain lambda.min and retrain using this
    tlambda.min_a=lambv[which.min(tavgMSElambda_a)]
    
    
    ##Transfer Retraining ##
    
    Delta_final_a=Delta(tlambda.min_a, X_tr_a, Y_tr_a, betaproxy_a)
    
    newdeltacoef_row<- data.frame(Class="Gold", Env=Env_names[Gold_index], coefType="Delta_ARR", Trait=trait, Partition= i, t(as.matrix(Delta_final_a)))
    colnames(newdeltacoef_row)<-c("Class","Env", "coefType", "Trait", "Partition", 0:p)
    newCoefs <- rbind(newCoefs, newdeltacoef_row)
    
    transferprediction_a<-Prediction_Trans(X_tst_a,Delta_final_a, betaproxy_a)
    #Scaling back step
    transferprediction_a<- transferprediction_a*scalingfactor
    
    #We prepare transfer prediction and add it to newPredictions
    tnewprediction_row <- data.frame(Env=Env_names[Gold_index], Trait=trait, predType="Transfer_ARR", Partition=i, t(as.matrix(transferprediction_a)))
    
    newPredictions <- rbind(newPredictions, tnewprediction_row)
    
    # We evaluate and report on newResults
    
    tNRMSE_a<-(sqrt(mse(Y_tst_a, transferprediction_a)))/mean(Y_tst_a)
    tCor_a <- cor(Y_tst_a, as.matrix(transferprediction_a))
    
    cat(i, "tNRMSE:", tNRMSE_a, "\n")
    
    
    tnewresult_row <- data.frame(
      Class="Gold",
      Env=Env_names[Gold_index],
      Trait=trait,
      Model="Transfer_ARR",
      Partition = i, 
      NRMSE = tNRMSE_a, 
      Cor = tCor_a, 
      Lambda.min = tlambda.min_a)
    colnames(tnewresult_row)=c("Class", "Env", "Trait", "Model", "Partition", "NRMSE", "Cor", "Lambda.min")
    newResults <- rbind(newResults, tnewresult_row)
  }
  
}

#Printing out the results


#We load coefficients, results and predictions to the existing csvs for this particular transfer

write.table(newCoefs, file = paste0(datalabel,"_TraditionalRidge_",transferlabel,"_Coef_",nfold,"cv.csv"), row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
write.table(newResults, file = paste0(datalabel,"_TraditionalRidge_",transferlabel,"_Result_",nfold,"cv.csv"), row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
write.table(newPredictions, file = paste0(datalabel,"_TraditionalRidge_",transferlabel,"Predictions_",nfold,"cv.csv"), row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
#Now we load the results and plot them

allResults<-read.csv(paste0(datalabel,"_TraditionalRidge_",transferlabel,"_Result_",nfold,"cv.csv"))
allResults$Model <- factor(allResults$Model, levels = c("RR", "ARR", "Transfer_RR", "Transfer_ARR"))

#### Plotting section ####


fullcorPlot<-ggplot(allResults[allResults$Class=="Gold",], aes(x = Model, y = Cor)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1, height = 0), colour = "blue", alpha = 0.5) +
  ggtitle(paste0("Dataset: ",datalabel)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12),  
    axis.title.x = element_text(face = "bold", size=12), 
    axis.title.y = element_text(face = "bold", size=13),
    strip.text = element_text(size = 12)
  ) + 
  facet_wrap(~Trait, labeller = label_both, nrow=1)

ggsave(paste0(datalabel,"_fullcompareTransfer_",transferlabel,"_Correlation_Plot.png"), 
       plot = fullcorPlot, device = "png", width = 2.5*length(Trait_names), height = 5, units = "in", dpi = 300)

fullNRMSEPlot<-ggplot(allResults[allResults$Class=="Gold",], aes(x = Model, y = NRMSE)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1, height = 0), colour = "blue", alpha = 0.5) +
  ggtitle(paste0("Dataset: ",datalabel)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12),  
    axis.title.x = element_text(face = "bold", size=12), 
    axis.title.y = element_text(face = "bold", size=13),
    strip.text = element_text(size = 12)
  ) + 
  facet_wrap(~Trait, labeller = label_both, nrow=1)

ggsave(paste0(datalabel,"_fullcompareTransfer_",transferlabel,"_NRMSE_Plot.png"), 
       plot = fullNRMSEPlot, device = "png", width = 2.5*length(Trait_names), height = 5, units = "in", dpi = 300)

file.rename(paste0(datalabel,"_TraditionalRidge_",transferlabel,"_Coef_",nfold,"cv.csv"),
            gsub("TraditionalRidge", "PostTransfer", paste0(datalabel,"_TraditionalRidge_",transferlabel,"_Coef_",nfold,"cv.csv")))

file.rename(paste0(datalabel,"_TraditionalRidge_",transferlabel,"_Result_",nfold,"cv.csv"),
            gsub("TraditionalRidge", "PostTransfer", paste0(datalabel,"_TraditionalRidge_",transferlabel,"_Result_",nfold,"cv.csv")))

file.rename(paste0(datalabel,"_TraditionalRidge_",transferlabel,"Predictions_",nfold,"cv.csv"),
            gsub("TraditionalRidge", "PostTransfer", paste0(datalabel,"_TraditionalRidge_",transferlabel,"Predictions_",nfold,"cv.csv")))
