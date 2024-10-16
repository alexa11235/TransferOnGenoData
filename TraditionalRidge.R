rm(list = ls())  


library(glmnet)
library(ggplot2)

set.seed(2713)
load("Wheat_1.RData")

Markers<-dat_ls$Markers
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

fold_indices<- data.frame()
coef_df<-data.frame()
allResults<-data.frame()
allPredictions<-data.frame()

R2v =  exp(seq(log(0.01),log(0.99),length=100))
psxx = mean(apply(X[, -1],1,function(x)sum(x^2)))
plambv = (1-R2v)/(R2v/psxx)#R2v

### for loop for each trait
for (t in 1:length(Trait_names)) {
  Y_proxy<-Y_proxydata[,t+2]
  Y_gold<-Y_golddata[,t+2]
  trait<-Trait_names[t]
  
  ## We do the 10 partitions of X, Y_gold and Y_proxy
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
    newfoldrow<-data.frame(Trait_index=t, Partition=i, t(Grpv_k))
    fold_indices<-rbind(fold_indices, newfoldrow)
    
    scalingfactor<-  mean(Y_tr)/mean(Y_proxy) 
    Y_tr<-Y_tr/scalingfactor

    #### Gold Section (RR) #####
    
    Ridge <- cv.glmnet(
      x = X_tr[,-1],
      y = Y_tr,
      lambda = lambv,
      family = "gaussian",
      alpha = 0,
      type.measure = 'mse',
      nfolds = nfold, 
      foldid = Grpv_k, 
      standardize = FALSE
    )
    
    Ridgeprediction <- predict(Ridge, newx = X_tst[,-1], s = "lambda.min", type = "response")*scalingfactor
    
    newbetacoef_row<- data.frame(Class="Gold",Env=Env_names[Gold_index] , coefType="Beta_Glmnet", Trait=trait, Partition= i, as.matrix(t(coef(Ridge, s = "lambda.min"))))
    colnames(newbetacoef_row)<-c("Class","Env", "coefType", "Trait", "Partition", 0:p)
    coef_df<-rbind(coef_df, newbetacoef_row)
    
    
    NRMSE <- (sqrt(mse(Y_tst, Ridgeprediction)))/mean(Y_tst)
    Cor <- cor(Y_tst, Ridgeprediction)
    
    #We prepare and save the observed and predicted values to  allPredictions
    newobserved_row <- data.frame(Env=Env_names[Gold_index], Trait=trait, predType="Observed", Partition=i, t(Y_tst)) 
    newprediction_row <- data.frame(Env=Env_names[Gold_index],Trait=trait, predType="RR", Partition=i, t(Ridgeprediction))
    allPredictions <- rbind(allPredictions, newobserved_row, newprediction_row)
    newresult_row <- data.frame(
      Class="Gold",
      Env=Env_names[Gold_index],
      Trait=trait,
      Model="RR",
      Partition = i,
      NRMSE = NRMSE,
      Cor = Cor,
      Lambda.min = Ridge$lambda.min
    )
    colnames(newresult_row) <- c("Class","Env", "Trait", "Model", "Partition", "NRMSE", "Cor", "Lambda.min")
    allResults <- rbind(allResults, newresult_row)
    
    
    #### Gold section (ARR) #####
    #vector for saving the averages MSEs for each lambda
    avgMSElambda<-c()

    #Matrices for saving innerMSEs
    innerMSEs<-matrix(NA, nrow = 100, ncol = nfold)
    #For each lambda we perform 10-fold cross validation and obtain an average MSE
    
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
        betahat <- A_inv%*%(t(X_innertraining)%*%Y_innertraining)
        iapred<-Prediction_ARR(X_innertesting, betahat)
        aMSE_inner<-mse(Y_innertesting, iapred)
        innerMSEs[l,k]<-aMSE_inner
        cat("ARR for", as.character(Env_names[Gold_index]), trait, i, k, l,  "innerMSE:", aMSE_inner, "\n")
  
      }
    }
    for (l in 1:length(lambv)){
      avgMSElambda[l]<- mean(innerMSEs[l,], na.rm=TRUE)
    }
    
    #Now we have 100 avgMSE (one for each lambda) so we can obtain lambda.min and retrain using this
    lambda.min=lambv[which.min(avgMSElambda)]
    
    #### ARR retraining ####
    
    Beta_final=Betahat(lambda.min, X_tr, Y_tr)
    
    anewbetacoef_row<- data.frame(Class="Gold",Env=Env_names[Gold_index], coefType="Beta_ARR", Trait=trait, Partition= i, t(as.matrix(Beta_final)))
    colnames(anewbetacoef_row)<-c("Class","Env", "coefType", "Trait", "Partition", 0:p)
    coef_df<-rbind(coef_df, anewbetacoef_row)
    
    analyticalprediction<-Prediction_ARR(X_tst,Beta_final)*scalingfactor
    
    #We prepare analytical prediction and add it 
    anewprediction_row <- data.frame(Env=Env_names[Gold_index],Trait=trait, predType="ARR", Partition=i, t(analyticalprediction))
    allPredictions <- rbind(allPredictions, anewprediction_row)
    
    # We evaluate and report on allResults  
    aNRMSE<-(sqrt(mse(Y_tst, analyticalprediction)))/mean(Y_tst)
    aCor <- cor(Y_tst, analyticalprediction)
    
    cat(i, "aNRMSE:", aNRMSE, "\n")
    
    anewresult_row <- data.frame(
      Class="Gold",
      Env=Env_names[Gold_index],
      Trait=trait,
      Model="ARR",
      Partition = i, 
      NRMSE = aNRMSE, 
      Cor = aCor, 
      Lambda.min = lambda.min)
    colnames(anewresult_row)=c("Class","Env", "Trait", "Model", "Partition", "NRMSE", "Cor", "Lambda.min")
    allResults <- rbind(allResults, anewresult_row)
    
  }
  
}

#Printing out the results
colnames(allPredictions) <- c("Env", "Trait", "predType", "Partition", seq_along(Y_tst))
rownames(allPredictions) <- NULL
allResults$Model <- factor(allResults$Model, levels = c("RR", "ARR"))
write.csv(allResults, file = paste0(datalabel,"_TraditionalRidge_",transferlabel,"_Result_",nfold,"cv.csv"), row.names = FALSE)
write.csv(allPredictions, file = paste0(datalabel,"_TraditionalRidge_",transferlabel,"Predictions_",nfold,"cv.csv"), row.names = FALSE)
write.csv(coef_df, file = paste0(datalabel,"_TraditionalRidge_",transferlabel,"_Coef_",nfold,"cv.csv"), row.names = FALSE)
write.csv(fold_indices, file = paste0(datalabel,"_TraditionalRidge_",transferlabel,"_Grpv_",nfold,"cv.csv"), row.names = FALSE)

#### Plotting section ####

corPlot<-ggplot(allResults[allResults$Class=="Gold",], aes(x = Model, y = Cor)) +
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

ggsave(paste0(datalabel,"_TraditionalRidge_",transferlabel,"_Correlation_Plot.png"), 
       plot = corPlot, device = "png", width = 2.5*length(Trait_names), height = 5, units = "in", dpi = 300)

NRMSEPlot<-ggplot(allResults[allResults$Class=="Gold",], aes(x = Model, y = NRMSE)) +
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

ggsave(paste0(datalabel,"_TraditionalRidge_",transferlabel,"_NRMSE_Plot.png"), 
       plot = NRMSEPlot, device = "png", width = 2.5*length(Trait_names), height = 5, units = "in", dpi = 300)
