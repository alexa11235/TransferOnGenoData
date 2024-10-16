This repository includes the code and datasets for the paper titled:

 "Harnessing Transfer Learning for Genomic Prediction Optimization"

-----------------------------------------------------------------------------------------

Last edited by Barajas Ramírez Eduardo 
               Montesinos Lopez Osval

For any questions please contact oamontes1@ucol.mx or ebarajas2@ucol.mx
--------------------------------------------------------------------------------------------

We include 2 R files

1. TraditionalRidge.R
2. Transfer.R

and 11 datasets 

EYT_1
EYT_2
EYT_3
Indica
Japonica
Wheat_1
Wheat_2 
Wheat_3 
Wheat_4
Wheat_5 
Wheat_6 

The R files are written for the Wheat_1 dataset but are easily customizable for any dataset.

For both files:

Line 5: load("Wheat_1.RData") 
Can be modified to load a different dataset

Line 10: Partitions <-10 
Can be modified to change the number of outer partitions performed previous to the cross validation.

Line 11: nfold <- 10
Can be changed to customize the number of folds in the k-fold cross validation throughout each code. It ssuld be noted though that the same number of folds are used in both models (RR and ARR in file 1, and Transfer_RR and Transfer_ARR in file 2)

Line 12: datalabel <- "Wheat_1"
Can and should be customized to the label that the user wants to see on the plots and output files. It should be indicative of the dataset.

Line 13 and 14: 
Proxy_index <-2
Gold_index <-1
Indicate the indices of the environments to be used as proxy and gold respectively. The index is the position in which they appear on the dataset. A trivial restriction is that these numbers should be different.
Whenever this are modified, a new transfer instance is being trained

Line 19:
transferlabel<-paste0(Env_names[Proxy_index],"to",Env_names[Gold_index])
Can also technically be modified but one should be careful since as written, depends on line 17.

Lines 15 to 18 and Lines 20-24 can vary slightly depending on the dataset, but the general idea remains the same.
------------------------------------------------------------------------------------------

To reproduce the results showed in the paper (for the "Wheat_1" dataset) you can follow the following steps:

Step 1) Run file 1 (TraditionalRidge.R) on the same location of the dataset to train the RR and ARR modelst. This outputs 2 images and 4 csv files:

-Image 1 ("Wheat_1_TraditionalRidge_YT_14_15toYT_13_14_Correlation_Plot.png") which compares the RR model to the ARR model in terms of Pearson's correlation
-Image 2 (Wheat_1_TraditionalRidge_YT_14_15toYT_13_14_NRMSE_Plot.png) which compares the the RR model to the ARR model in terms of NRMSE. 



-CSV 1 (Wheat_1_TraditionalRidge_YT_14_15toYT_13_14_Coef_10cv.csv) which contains the coeficients obtained for each partition and for each model
-CSV 2 (Wheat_1_TraditionalRidge_YT_14_15toYT_13_14_Grpv_10cv.csv) which contains the fold indices used in the cross validation
-CSV 3 (Wheat_1_TraditionalRidge_YT_14_15toYT_13_14_Result_10cv.csv) which contains the results (Cor, NRMSE and Lambda.min) for each model, per partition and per trait.
-CSV 4 (Wheat_1_TraditionalRidge_YT_14_15toYT_13_14_Result_10cv.csv) which contains the predictions per trait, per partition and per model.


Step 2) Run file 2 (Transfer.R) on the same location of the dataset ("Wheat_1") to train the Transfer_RR and Transfer_ARR models. This outputs two more images -which now compare all 4 models- and completes and renames CSV 1, CSV 3, CSV 4.

-Image 3 (Wheat_1_fullcompareTransfer_YT_14_15toYT_13_14_Correlation_Plot.png) Compares the 4 models (RR, ARR, Transfer_RR and Transfer_ARR) in terms of Pearson's correlation and should coincide statistically with Appendix Figure A4 (b)
-Image 4 (Wheat_1_fullcompareTransfer_YT_14_15toYT_13_14_NRMSE_Plot.png) Compares the 4 models (RR, ARR, Transfer_RR and Transfer_ARR) in terms of NRMSE and should coincide statistically with Appendix Figure A4 (b)

-CSV 1 is now Wheat_1_PostTransfer_YT_14_15toYT_13_14_Coef_10cv and contains coeficients of all 4 models.
-CSV 3 is now Wheat_1_PostTransfer_YT_14_15toYT_13_14_Result_10cv and contains results for all 4 models.
-CSV 4 is now Wheat_1_PostTransfer_YT_14_15toYT_13_14Predictions_10cv and contains predictions for all 4 models.

It should be mentioned that in order to reproduce all the figures included on the paper, multiple transfer instances were ran and then unified in illustrative images. The tables were similarly refactored, but the general essence and improvement of the models is kept. 


