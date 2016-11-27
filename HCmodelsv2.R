#Utilizing free surfer pull script, we now initialize baseline + neuroimaging models for Anxiety and Depression
data_dir = '/Users/rhuang/Documents/Classes/Fall_16/Stat_Research/R/data' #freesurfer_out.csv
script_dir = '/Users/rhuang/Documents/Classes/Fall_16/Stat_Research/R' #directory of this Rscript

##### 
setwd(data_dir)
#Load cleaned data
disease_dsm_clean = read.csv('master_list_clean.csv')
data_mri_clean = read.csv('mri_clean.csv')

# HC_ID <- subset(data_mri_clean, (MIND_ID %in% 1:125))$MIND_ID #114 healthy controls post outlier removal
# dep_ID <- subset(disease_dsm_clean, (depression == 1) & (anxiety == 0))$MIND_ID
# anx_ID <- subset(disease_dsm_clean, (anxiety == 1) & (depression == 0))$MIND_ID
# 
# HC_ID_v <- sample(HC_ID, floor(1/6*115)) 
# dep_ID_v <- sample(dep_ID, floor(1/6*length(dep_ID)))
# anx_ID_v <- sample(anx_ID, floor(1/6*length(anx_ID)))
# 
# # saveRDS(list(HC = HC_ID_v, dep = dep_ID_v, anx= anx_ID_v), file = 'valid_set.rds')
# 
# 
# HC_ID_ts <- sample(setdiff(HC_ID, HC_ID_v), floor(1/6*115))
# dep_ID_ts <- sample(setdiff(dep_ID,dep_ID_v), floor(1/6*length(dep_ID)))
# anx_ID_ts <- sample(setdiff(anx_ID,anx_ID_v), floor(1/6*length(anx_ID)))
# 
# # saveRDS(list(HC = HC_ID_ts, dep = dep_ID_ts, anx= anx_ID_ts), file = 'ts_set.rds')
# 
# 
# HC_ID_tr <- setdiff(HC_ID, c(HC_ID_v,HC_ID_ts))
# dep_ID_tr <- setdiff(dep_ID, c(dep_ID_v,dep_ID_ts))
# anx_ID_tr <- setdiff(anx_ID, c(anx_ID_v,anx_ID_ts))
# saveRDS(list(HC = HC_ID_tr, dep = dep_ID_tr, anx= anx_ID_tr), file = 'tr_set.rds')

##
readRDS("tr_set.rds");readRDS("valid_set.rds"); readRDS("ts_set.rds")


HC_ID_v <- readRDS("valid_set.rds")[[1]]
dep_ID_v <- readRDS("valid_set.rds")[[2]]
anx_ID_v <- readRDS("valid_set.rds")[[3]]

HC_ID_ts <- readRDS("ts_set.rds")[[1]]
dep_ID_ts <- readRDS("ts_set.rds")[[2]]
anx_ID_ts <- readRDS("ts_set.rds")[[3]]

HC_ID_tr <- readRDS("tr_set.rds")[[1]]
dep_ID_tr <- readRDS("tr_set.rds")[[2]]
anx_ID_tr <- readRDS("tr_set.rds")[[3]]

#initialize libraries
library(glmnet)
#intialize input data
valid_ID <- c(HC_ID_v, dep_ID_v)
train_ID <- c(HC_ID_tr, dep_ID_tr)
test_ID <-  c(HC_ID_ts, dep_ID_ts)
  
valid_ID <- c(HC_ID_v, anx_ID_v)
train_ID <- c(HC_ID_tr, anx_ID_tr)
test_ID <-  c(HC_ID_ts, anx_ID_ts)
 
data_mri_tr <- subset(data_mri_clean, (MIND_ID %in% train_ID))
data_mri_v <- subset(data_mri_clean, (MIND_ID %in% valid_ID))
data_mri_ts <- subset(data_mri_clean, (MIND_ID %in% test_ID))

data_mri_scaled_tr <- scale(data_mri_tr[,-1], center = T, scale = T)
data_mri_scaled_v <- scale(data_mri_v[,-1], center = T, scale = T)
data_mri_scaled_ts <- scale(data_mri_ts[,-1], center = T, scale = T)

Yv <- subset(disease_dsm_clean, (MIND_ID %in% valid_ID))[,3]; Ytr <- subset(disease_dsm_clean, (MIND_ID %in% train_ID))[,3] #Depression indicator based on DSM
Yts <- subset(disease_dsm_clean, (MIND_ID %in% test_ID))[,3]
#anx
Yv <- subset(disease_dsm_clean, (MIND_ID %in% valid_ID))[,2]; Ytr <- subset(disease_dsm_clean, (MIND_ID %in% train_ID))[,2] #Depression indicator based on DSM
Yts <- subset(disease_dsm_clean, (MIND_ID %in% test_ID))[,2]

Xtr <- data_mri_scaled_tr ; Xv <- data_mri_scaled_v; Xts <- data_mri_scaled_ts
summary(Xtr)
##### Data exploration/summary

#####
str(Xtr)
##Feature Selection
randomLasso <- function(x = Xtr, y = as.factor(Ytr), times = 5, sample_size = .8, feature_prop = .6) {
  #Step1: Generate
  feature_compiler <- c()
  
  for (i in 1:times) {
    cur_subsample <- sample(nrow(Xtr), sample_size*nrow(Xtr), replace = TRUE) #Bootstrap samples
    rand_feature <- sample(ncol(Xtr), feature_prop*ncol(Xtr), replace = FALSE) #random features
    x_sub <- x[cur_subsample,rand_feature]
    y_sub <- y[cur_subsample]
    cur_cv <- cv.glmnet(x_sub, y_sub,family = "binomial", standardize = FALSE,type.measure = 'auc', alpha = 1, nfolds = 3)
    cur_fit <- glmnet(x_sub, y_sub, family = "binomial", standardize = FALSE,
                      lambda = cur_cv$lambda.min, alpha = 1)
    features.sparse <- as.matrix(coef(cur_fit))
    features <- features.sparse[(which(features.sparse != 0)),]
    feature_compiler <- c(feature_compiler, features)
    
  }
  feature_split <- split(feature_compiler,names(feature_compiler))
  feature_imp <- abs(unlist(lapply(feature_split, sum))/times) ; feature_imp <- feature_imp[-1]
   ##Step 2: Select
  feature_compiler2 <- c()
  
  for (i in 1:times) {
    cur_subsample <- sample(nrow(Xtr), sample_size*nrow(Xtr), replace = TRUE) #Bootstrap samples
    rand_feature <- sample(attributes(feature_imp)$names, feature_prop*length(feature_imp), prob = feature_imp, replace = FALSE) #random features
    x_sub <- x[cur_subsample,rand_feature]
    y_sub <- y[cur_subsample]
    cur_cv <- cv.glmnet(x_sub, y_sub,family = "binomial", standardize = FALSE,type.measure = 'auc', alpha = 1, nfolds = 3)
    cur_fit <- glmnet(x_sub, y_sub, family = "binomial", standardize = FALSE,
                      lambda = cur_cv$lambda.min, alpha = 1)
    features.sparse <- as.matrix(coef(cur_fit))
    features <- features.sparse[(which(features.sparse != 0)),]
    feature_compiler2 <- c(feature_compiler2, features)
    
  }
  feature_split2 <- split(feature_compiler2,names(feature_compiler2))
  feature_imp2 <- unlist(lapply(feature_split2, sum))/times; feature_imp2 <- feature_imp2[-1]
  return(feature_imp2)
}
features <- randomLasso(times = 500, sample_size = .8)
# write.csv(sort(abs(features), decreasing = T), 'random_lasso_feat.csv')
str(features)

1/137

table(features)
sort(abs(features))
length(features)
features_stable <- sort(abs(features[(which(abs(features)>(30/137)))]),decreasing = T)

tmp <- sample(attributes(features)$names, 5)
Xtr[,tmp]

sample(c(1,2,3),1, prob = c(1,50,50))
#Bootstrapping
# Bootstrap 95% CI for regression coefficients
library(boot)
#non-scaled train data
names(features_stable)
names(dep_coef)
#random lasso features
Xtr <- data_mri_scaled_tr[,names(features_stable)] ; Xv <- data_mri_v[,names(features_stable)]
#lasso features
Xtr <- data_mri_scaled_tr[,names(dep_coef)] ; Xv <- data_mri_v[,names(dep_coef)]
Xtr <- data_mri_scaled_tr[,names(anx_coef)] ; Xv <- data_mri_v[,names(anx_coef)]



data_ols_tr <- cbind(dep = Ytr, Xtr)
lm(dep~., data = as.data.frame(data_ols_tr))
# function to obtain regression weights
bs <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula, data=d)
  return(coef(fit))
}
# bootstrapping with 1000 replications
results <- boot(data=as.data.frame(data_ols_tr), statistic=bs,
                R=1000, formula=dep~.)

results_anx <- boot(data=as.data.frame(data_ols_tr), statistic=bs,
                R=1000, formula=dep~.)
str(Xtr)
plot(results, index=1) # intercept
plot(results, index=3) # wt
plot(results, index=14) # disp

colnames(Xtr)
#Normal Lasso Bootstraps
boot.ci(results, type="bca", index=1) # intercept
boot.ci(results, type="bca", index=2) # "lh_parsopercularis_thickness"
boot.ci(results, type="bca", index=3) # 
boot.ci(results, type="bca", index=4) # lh_posteriorcingulate_thickness"
boot.ci(results, type="bca", index=5) # 
boot.ci(results, type="bca", index=6) # "Right.Hippocampus.x"
boot.ci(results, type="bca", index=7) # "rh_bankssts_area"
boot.ci(results, type="bca", index=8) # "Left.Cerebellum.White.Matter.y"
boot.ci(results, type="bca", index=9) # 
boot.ci(results, type="bca", index=10) # 
# boot.ci(results, type="bca", index=11) # lh_rostralanteriorcingulate_thickness
# boot.ci(results, type="bca", index=12) # lh_rostralanteriorcingulate_thickness
# boot.ci(results, type="bca", index=13) # lh_rostralanteriorcingulate_thickness
# boot.ci(results, type="bca", index=14) # lh_rostralanteriorcingulate_thickness
# boot.ci(results, type="bca", index=22) # lh_rostralanteriorcingulate_thickness

#Normal Anx Lasso
boot.ci(results_anx, type="bca", index=1) # intercept
boot.ci(results_anx, type="bca", index=2) # X4th.Ventricle.y
boot.ci(results_anx, type="bca", index=3) # 
boot.ci(results_anx, type="bca", index=4) # 
boot.ci(results_anx, type="bca", index=5) # rh_isthmuscingulate_area
boot.ci(results_anx, type="bca", index=6) # 
boot.ci(results_anx, type="bca", index=7) # "lh_transversetemporal_thickness"
boot.ci(results_anx, type="bca", index=8) # "Right.Cerebellum.Cortex.x" 
boot.ci(results_anx, type="bca", index=9) # "Left.Cerebellum.White.Matter.x" 
boot.ci(results_anx, type="bca", index=10) # 
boot.ci(results_anx, type="bca", index=11) # 

colnames(Xtr)


##
dep_lam <- cv.glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,type.measure = 'auc', alpha = 1)
depfit <- glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,
                 lambda = dep_lam$lambda.min, alpha = 1)

coef(depfit)

dep_coef <- as.matrix(depfit$beta)
dep_coef <- dep_coef[dep_coef!=0,]

sum(predict(depfit, newx = Xv, type = 'class') == Yv)/length(Yv)
dep_pred = predict(depfit, newx = Xv)
sum(Yv)/length(Yv)

result <- data.frame(MIND_ID = valid_ID, predict = predict(depfit, newx = Xv, type = 'class'), actual = Yv)
table(result[,-1])

#Test Case
Yts <- subset(disease_dsm_clean, (MIND_ID %in% test_ID))[,2]

sum(predict(depfit, newx = Xts, type = 'class') == Yts)/length(Yts)

coef(anxfit)

anx_coef <- as.matrix(anxfit$beta)
anx_coef <- anx_coef[anx_coef!=0,]

sum(predict(depfit, newx = Xts, type = 'class') == Yts)/length(Yts)
anx_pred = predict(anxfit, newx = Xv)
table(data.frame(predict = predict(depfit, newx = Xts, type = 'class'), actual = Yts))


#SVM
library(e1071)
for (i in 1:10){
  obj <- tune.svm(as.matrix(Xtr),as.factor(Ytr),kernel="polynomial",type="C-classification",
                  cost = 10^((-3):3),scale=FALSE)
  print(obj$best.parameters)
}
dep.model.svm<-svm(as.matrix(Xtr),as.factor(Ytr),cost=10,kernel="polynomial",
                   type ="C-classification",scale=FALSE)
pred<-predict(dep.model.svm,as.matrix(Xv))
sum(pred==Yv)/length(Yv)
table(pred,Yv)

#Test case
pred<-predict(dep.model.svm,as.matrix(Xts))
sum(pred==Yts)/length(Yts)
table(pred,Yts)
# #####Anxiety Binary Models
Yv <- subset(disease_dsm_clean, (MIND_ID %in% valid_ID))[,2]; Ytr <- subset(disease_dsm_clean, (MIND_ID %in% train_ID))[,2] #Anxiety indicator based on DSM

anx_lam <- cv.glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,type.measure = 'auc', alpha = 1)
anxfit <- glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,
                 lambda = anx_lam$lambda.min, alpha = 1)

coef(anxfit)
anx_coef <- as.matrix(anxfit$beta)
anx_coef[anx_coef!=0,]
sum(predict(anxfit, newx = Xv, type = 'class') == Yv)/length(Yv)
anx_pred = predict(anxfit, newx = Xv)

table(data.frame(predict = predict(anxfit, newx = Xv, type = 'class'), actual = Yv))

##TEST SET
Yts <- subset(disease_dsm_clean, (MIND_ID %in% test_ID))[,2]

sum(predict(anxfit, newx = Xts, type = 'class') == Yts)/length(Yts)

coef(anxfit)
anx_coef <- as.matrix(anxfit$beta)
anx_coef[anx_coef!=0,]
sum(predict(anxfit, newx = Xts, type = 'class') == Yts)/length(Yts)
anx_pred = predict(anxfit, newx = Xv)

table(data.frame(predict = predict(anxfit, newx = Xts, type = 'class'), actual = Yts))
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




####SVM
