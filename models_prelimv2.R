#Utilizing free surfer pull script, we now initialize baseline + neuroimaging models for Anxiety and Depression
data_dir = '/Users/rhuang/Documents/Classes/Fall_16/Stat_Research/R/data' #freesurfer_out.csv
script_dir = '/Users/rhuang/Documents/Classes/Fall_16/Stat_Research/R' #directory of this Rscript

##### 
setwd(data_dir)
#reading in data
data <- read.csv("freesurfer_out.csv", header = T)[,-1] #free-surfer
diag = load(file = 'menninger_clincial.rda') #Clinical data should be in data directory
clinical = get(diag); rm(dataset) 

#Formatting
data$MIND_ID <- as.numeric(substring(data$MIND_ID, 6, 9)) #format mind_ID to just the number so we can match with diagnosis data

#Begin experiment / Naming conventions fixed from here on out
# total_ID <- subset(data, MIND_ID %in% c(1:125,clinical$MIND_ID))$MIND_ID
# split_samples <- sample(total_ID, floor(2/6 * length(total_ID)))
# test_ID <- split_samples[1:80]
# valid_ID <- split_samples[81:161]
# write.csv(test_ID, 'test_ID.csv')
# write.csv(valid_ID, 'valid_ID.csv')


test_ID <- read.csv('test_ID.csv')[,2] #load test ids saved randomly splitting
valid_ID <- read.csv('valid_ID.csv')[,2] #load test ids saved randomly splitting

data_mri = subset(data, MIND_ID %in% c(1:125,clinical$MIND_ID) & !(MIND_ID %in% test_ID)) #subsets to the patients that have both clinical and neuroimaging results
data_clinical = subset(clinical, MIND_ID %in% data$MIND_ID & !(MIND_ID %in% test_ID))

#Create disease indicator matrix based on raw label formatting
subset.regex <- function(data, label = list("name", "ID"), pattern = c("\\([1234567890]{3}","NOS", "Disorder"), exceptions = c("SCID2PDNOS")) {
  #Initialize
  name <- label[[1]]
  ID <- label[[2]]
  #Grab disease codes [+1 for PST, +2 for CUR]
  index1 <- grep(pattern[1], name, perl=TRUE, value=FALSE)
  index1 <- index1[which(name[index1+2] == "Current")]
  
  #Grab NOS [+1 for PST, +2 for CUR]
  index2.tmp <- ID[grep(pattern[2], sub('.*(?=.{3}$)', '', ID, perl=T))]
  #Grab Disorders without disease code [No PST or CUR, assuming all current for now]
  #index3 <- grep(pattern[3], ID, perl=TRUE, value=FALSE)
  
  #Remove exceptions
  index2 <- index2.tmp[index2.tmp!=exceptions]
  index2 <- match(index2, ID)
  #print(colnames(data[,c(index1,index2)]))
  
  #Aggregate index w/ MIND_ID
  index <- unique(c(index1 + 2, index2 + 2))
  data[,index]
}
disease_raw <- (subset.regex(data_clinical, label =list(labels,colnames(data_clinical)))) #all current diseases

#DSM-formatting conversion
depression_dsm_names <- c("DYSTHYMICCUR","MDDSINGLECUR","MDDRECURCUR","DEPRESSCUR")
anxiety_dsm_names <- c("PANICWOAGORCUR","PANICWITHAGORCUR","AGORAPHOBIACUR","SOCIALPHOBIACUR",
                       "SPECPHOBIACUR","OCDCUR","GADCUR","ANXGMCCUR","PTSDCUR","ANXNOSCUR")
  
num_control <- sum(data_mri$MIND_ID <= 125)

depression_indicator <- apply(disease_raw[,depression_dsm_names],1, sum, na.rm = TRUE)
depression_indicator[depression_indicator>0] <- 1
depression_indicator <- c(rep(0,times = num_control), depression_indicator)


anxiety_indicator <- apply(disease_raw[,anxiety_dsm_names],1, sum, na.rm = TRUE)
anxiety_indicator[anxiety_indicator>0] <- 1
anxiety_indicator <- c(rep(0,times = num_control), anxiety_indicator)


sum(anxiety_indicator) #190 of 369 with anxiety
sum(depression_indicator)  #208 of 369 with anxiety
sum((anxiety_indicator == 1)  & (depression_indicator == 1)) #132 of 369 have anxiety and depression
sum((anxiety_indicator == 0)  & (depression_indicator == 0)) #132 of 369 have anxiety and depression

disease_dsm <- data.frame(MIND_ID = data_mri$MIND_ID,anxiety = anxiety_indicator, depression = depression_indicator) #dsm9 grouped-disease matrix; depression and anxiety only for now
#Demographic data pulling

#Data pre-processing for model input
na_mri <- which(colSums(data_mri) == 0) ; data_mri <- data_mri[,-na_mri]#remove NA from MRI data
data_mri_scaled <- scale(data_mri[,-1], center = T, scale = T)

#Clean outliers from mri data: Pipeline system
outlier_index <- which(data_mri_scaled > 5, arr.ind = T)
data_mri_scaled[outlier_index]
outlier_MINDID <- data_mri[outlier_index[,1],1] #MIND ID of outliers
outlier_feature <- colnames(data_mri_scaled)[outlier_index[,2]] #feature name of outlier values
outlier_values_scaled <- data_mri_scaled[outlier_index] ; outlier_values_unscaled <- data_mri[,-1][outlier_index]

outlier_df <- data.frame(MIND_ID = outlier_MINDID, feature = outlier_feature,
                         scaled_value = outlier_values_scaled, unscaled_value = outlier_values_unscaled)

length(unique(outlier_MINDID))
###Final renames
data_mri_clean <- subset(data_mri, !(MIND_ID %in% outlier_MINDID))
disease_dsm_clean <- subset(disease_dsm, !(MIND_ID %in% outlier_MINDID))

#####
#initialize libraries
library(glmnet)
#intialize input data
data_mri_tr <- subset(data_mri_clean, !(MIND_ID %in% valid_ID))
data_mri_v <- subset(data_mri_clean, (MIND_ID %in% valid_ID))

data_mri_scaled_tr <- scale(data_mri_tr[,-1], center = T, scale = T)
data_mri_scaled_v <- scale(data_mri_v[,-1], center = T, scale = T)

Yv <- subset(disease_dsm_clean, (MIND_ID %in% valid_ID))[,3]; Ytr <- subset(disease_dsm_clean, !(MIND_ID %in% valid_ID))[,3] #Depression indicator based on DSM
Xtr <- data_mri_scaled_tr ; Xv <- data_mri_scaled_v
summary(Xtr)
##### Data exploration/summary
# summary(data_mri_scaled_v) #only finding is that there are some extreme outliers
# colnames(Xtr_demo)
# ncol(Xtr_full) - length(colnames(Xtr_demo))
# colnames(Xtr_full)
# summary(Xtr_full)
#####
### cycle for doing 100 cross validations
### and take the average of the mean error curves
### initialize vector for final data.frame with Mean Standard Errors
Lambdas <- function(...) {
  cv <- cv.glmnet(...)
  return(data.table(cvm=cv$cvm, lambda=cv$lambda))
}

OptimLambda <- function(k, ...) {
  # Returns optimal lambda for glmnet.
  #
  # Args:
  #   k: # times to loop through cv.glmnet.
  #   ...: Other args passed to cv.glmnet.
  #
  # Returns:
  #   Lambda associated with minimum average CV error over runs.
  #
  # Example:
  #   OptimLambda(k=100, y=y, x=x, alpha=alpha, nfolds=k)
  #
  require(parallel)
  require(data.table)
  require(plyr)
  MSEs <- data.table(rbind.fill(mclapply(seq(k), function(dummy) Lambdas(...))))
  return(MSEs[, list(mean.cvm=mean(cvm)), lambda][order(mean.cvm)][1])
}

min_lam <- OptimLambda(k=10, y=as.factor(Ytr), x=Xtr, family = 'binomial',type.measure="auc",
                         standardize = FALSE, alpha=1)

dep_lam <- cv.glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,type.measure = 'auc', alpha = 1)
depfit <- glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,
                  lambda = dep_lam$lambda.min, alpha = 1)

coef(depfit)
dep_coef <- as.matrix(depfit$beta)
dep_coef[dep_coef!=0,]

sum(predict(depfit, newx = Xv, type = 'class') == Yv)/length(Yv)
dep_pred = predict(depfit, newx = Xv)
sum(Yv)/length(Yv)

table(data.frame(predict = predict(fullfit, newx = Xv, type = 'class'), actual = Yv))

#df100
print(glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE, alpha = 1)) #Check lambda where df = 100
df100_features <- as.matrix(glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,
                                   alpha = 1, lambda = .001157)$beta)
write.csv(df100_features, "features_df100.csv")

# 
# write.csv(data.frame(predict = predict(fullfit, newx = Xv_full, type = 'class'), actual = Yv), "depress_full251016.csv")
# g_full_depression <- glmnet(Xtr_full, as.factor(Ytr),family = "binomial", standardize = FALSE, alpha = 1)


#key_predictor_plot(g_full_depression, alpha_min = .1) + ggtitle("Elastic net (alpha .5) regularization path - Binary Depression Model")

#####Anxiety Binary Models
Yv <- subset(disease_dsm_clean, (MIND_ID %in% valid_ID))[,2]; Ytr <- subset(disease_dsm_clean, !(MIND_ID %in% valid_ID))[,2] #Anxiety indicator based on DSM

anx_lam <- cv.glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,type.measure = 'auc', alpha = 1)
anxfit <- glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,
                  lambda = min_lam$lambda.min, alpha = 1)

coef(anxfit)
anx_coef <- as.matrix(anxfit$beta)
anx_coef[anx_coef!=0,]
sum(predict(anxfit, newx = Xv, type = 'class') == Yv)/length(Yv)
anx_pred = predict(anxfit, newx = Xv)
sum(Yv)/length(Yv)

table(data.frame(predict = predict(anxfit, newx = Xv, type = 'class'), actual = Yv))
#df100
print(glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE, alpha = 1)) #Check lambda where df = 100
df100_features <- as.matrix(glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,
                                   alpha = 1, lambda = .001157)$beta)
write.csv(df100_features, "features_anxiety_df100.csv")

#ROCR metrics
library(ROCR)
#Anxiety
pred <- prediction(anx_pred, Yv)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, col=rainbow(10), main ="Anxiety ROC")
abline(a=0, b= 1)

acc.perf = performance(pred, measure = 'acc')
plot(acc.perf, main = "Anxiety Accuracy")
#Depression
pred <- prediction(dep_pred, Yv)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, col=rainbow(10), main = "Depression ROC")
abline(a=0, b= 1)

acc.perf = performance(pred, measure = 'acc')
plot(acc.perf, main = "Depression Accuracy")
#Plotter function
key_predictor_plot <- function(x, ..., num_predictors=10, alpha_min=0.4){
  stopifnot(require("ggplot2"))
  stopifnot(require("reshape2"))
  stopifnot(require("gtools"))
  
  if(inherits(x, "glmnet")){
    beta <- x$beta # Direct access of beta rather than coef() avoids intercept
  } else if(inherits(x, "ncvreg")){
    ## ncvreg includes an intercept term (unlike glmnet) so we omit it
    beta <- if (length(x$penalty.factor) == NROW(x$beta)){
      coef(x)
    } else {
      coef(x)[-1, , drop = FALSE] ## (intercept)
    }
  }
  lambda <- x$lambda
  
  active = glmnet::nonzeroCoef(beta)
  nactive = length(active)
  beta <- as.matrix(beta[active, , drop=FALSE])
  index = log(lambda)
  xlab=expression(log(lambda))
  
  ## Identify the first num_predictors to enter the model
  ## (or fewer, in the case when we add two variables simultaneously)
  key_predictors <- names(which(beta[,max(which(colSums(beta != 0) <= num_predictors))] != 0))
  if((length(key_predictors) > num_predictors) | (length(key_predictors) < 2)){
    ## If we have lots of predictors (e.g., ridge)
    ## just pull out the predictors with the largest |\hat{\beta}_j|
    ## and carry on
    key_predictors <- rownames(beta)[rank(abs(beta[,ncol(beta)]))][1:num_predictors]
  }
  
  
  ## Now we rearrange our data in a 'long' form preferred by ggplot2
  mb <- melt(beta)
  colnames(mb) <- c("Variable", "Step", "Beta")
  mb$Variable <- as.character(mb$Variable)
  mb <- cbind(mb, key_predictor=mb$Variable %in% key_predictors)
  mb <- cbind(mb, LogLambda=rep(index, each=nrow(beta)))
  mb <- cbind(mb, plot_color=ifelse(mb$key_predictor, mb$Variable, "Other"))
  
  ## And add some extra columns to use as plot parameters.
  g <- ggplot(mb, aes(x=LogLambda, y=Beta, group=Variable))
  if(!all(mb$key_predictor)){
    g <- g + geom_line(aes(col=plot_color, alpha=key_predictor))
    g <- g + scale_alpha_discrete(aes(alpha=key_predictor), guide="none", range=c(alpha_min, 1))
  } else {
    g <- g + geom_line(aes(col=plot_color))
  }
  
  g <- g + xlab(expression(log(lambda))) + ylab(expression(hat(beta[j])))
  g <- g + scale_color_discrete(guide=guide_legend(title="Predictor"), 
                                breaks=c(mixedsort(unique(mb$Variable)), "Other"))
  
  ## Conventionally, these are plotted with decreasing LogLambda
  g <- g + scale_x_reverse()
  
  g 
}