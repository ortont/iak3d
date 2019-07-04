# iak3d
Increment-averaged kriging for 3D prediction of soil properties


Save files to a directory, and edit the following four lines near the top of the egEdgeroi.R file:

wDir <- 'U://scripts/'
iakDir <- 'U://scripts/iakTests'
lmm2Dir <- 'U://scripts/fLMM2'
dataDir <- 'Z://ortont/Data/edgeroiEg'

wDir is your working directory
iakDir is where you save the iak files, and lmm2Dir is where you save the fLMM2.R file (could be the same place)
dataDir is where you want output to be saved

Example from Edgeroi dataset (data from GSIF and ithir packages) predicting clay content
