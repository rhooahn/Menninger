#Load cleaned data
disease_dsm_clean = read.csv('master_list_clean.csv')
data_mri_clean = read.csv('mri_clean.csv')

HC_ID <- subset(data_mri_clean, (MIND_ID %in% 1:125))$MIND_ID #114 healthy controls post outlier removal
dep_ID <- subset(disease_dsm_clean, (depression == 1) & (anxiety == 0))$MIND_ID
anx_ID <- subset(disease_dsm_clean, (anxiety == 1) & (depression == 0))$MIND_ID

HC_ID_v <- sample(HC_ID, floor(1/6*115))
dep_ID_v <- sample(dep_ID, floor(1/6*length(dep_ID)))
anx_ID_v <- sample(anx_ID, floor(1/6*length(anx_ID)))

HC_ID_ts <- sample(setdiff(HC_ID, HC_ID_v), floor(1/6*115))
dep_ID_ts <- sample(setdiff(dep_ID,dep_ID_v), floor(1/6*length(dep_ID)))
anx_ID_ts <- sample(setdiff(anx_ID,anx_ID_v), floor(1/6*length(anx_ID)))

HC_ID_tr <- setdiff(HC_ID, c(HC_ID_v,HC_ID_ts))
dep_ID_tr <- setdiff(dep_ID, c(dep_ID_v,dep_ID_ts))
anx_ID_tr <- setdiff(anx_ID, c(anx_ID_v,anx_ID_ts))

valid_ID <- c(HC_ID_v, dep_ID_v)
train_ID <- c(HC_ID_tr, dep_ID_tr)

#initialize libraries
library(glmnet)
#intialize input data
data_mri_tr <- subset(data_mri_clean, (MIND_ID %in% train_ID))
data_mri_v <- subset(data_mri_clean, (MIND_ID %in% valid_ID))

data_mri_scaled_tr <- scale(data_mri_tr[,-1], center = T, scale = T)
data_mri_scaled_v <- scale(data_mri_v[,-1], center = T, scale = T)

Yv <- subset(disease_dsm_clean, (MIND_ID %in% valid_ID))[,3]; Ytr <- subset(disease_dsm_clean, (MIND_ID %in% train_ID))[,3] #Depression indicator based on DSM
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
dep_lam <- cv.glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,type.measure = 'auc', alpha = 1)
depfit <- glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,
                 lambda = dep_lam$lambda.min, alpha = 1)

coef(depfit)
dep_coef <- as.matrix(depfit$beta)
dep_coef[dep_coef!=0,]

sum(predict(depfit, newx = Xv, type = 'class') == Yv)/length(Yv)
dep_pred = predict(depfit, newx = Xv, type = 'class')
sum(Yv)/length(Yv)

table(data.frame(predict = predict(depfit, newx = Xv, type = 'class'), actual = Yv))

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


#df100
print(glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE, alpha = 1)) #Check lambda where df = 100
df100_features <- as.matrix(glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,
                                   alpha = 1, lambda = .001157)$beta)
write.csv(df100_features, "features_df100.csv")

# 
# write.csv(data.frame(predict = predict(fullfit, newx = Xv_full, type = 'class'), actual = Yv), "depress_full251016.csv")
# g_full_depression <- glmnet(Xtr_full, as.factor(Ytr),family = "binomial", standardize = FALSE, alpha = 1)


#key_predictor_plot(g_full_depression, alpha_min = .1) + ggtitle("Elastic net (alpha .5) regularization path - Binary Depression Model")

# #####Anxiety Binary Models
# Yv <- subset(disease_dsm_clean, (MIND_ID %in% valid_ID))[,2]; Ytr <- subset(disease_dsm_clean, !(MIND_ID %in% valid_ID))[,2] #Anxiety indicator based on DSM
# 
# anx_lam <- cv.glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,type.measure = 'auc', alpha = 1)
# anxfit <- glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,
#                  lambda = anx_lam$lambda.min, alpha = 1)
# 
# coef(anxfit)
# anx_coef <- as.matrix(anxfit$beta)
# anx_coef[anx_coef!=0,]
# sum(predict(anxfit, newx = Xv, type = 'class') == Yv)/length(Yv)
# anx_pred = predict(anxfit, newx = Xv)
# sum(Yv)/length(Yv)
# 
# table(data.frame(predict = predict(anxfit, newx = Xv, type = 'class'), actual = Yv))
# #df100
# print(glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE, alpha = 1)) #Check lambda where df = 100
# df100_features <- as.matrix(glmnet(Xtr, as.factor(Ytr),family = "binomial", standardize = FALSE,
#                                    alpha = 1, lambda = 0.009828)$beta)
# write.csv(df100_features, "features_anxiety_df100.csv")
# 
# #ROCR metrics
# library(ROCR)
# #Anxiety
# pred <- prediction(anx_pred, Yv)
# perf <- performance(pred, measure = "tpr", x.measure = "fpr")
# plot(perf, col=rainbow(10), main ="Anxiety ROC")
# abline(a=0, b= 1)
# 
# acc.perf = performance(pred, measure = 'acc')
# plot(acc.perf, main = "Anxiety Accuracy")
# #Depression
# pred <- prediction(dep_pred, Yv)
# perf <- performance(pred, measure = "tpr", x.measure = "fpr")
# plot(perf, col=rainbow(10), main = "Depression ROC")
# abline(a=0, b= 1)
# 
# acc.perf = performance(pred, measure = 'acc')
# plot(acc.perf, main = "Depression Accuracy")
# #Plotter function
# key_predictor_plot <- function(x, ..., num_predictors=10, alpha_min=0.4){
#   stopifnot(require("ggplot2"))
#   stopifnot(require("reshape2"))
#   stopifnot(require("gtools"))
#   
#   if(inherits(x, "glmnet")){
#     beta <- x$beta # Direct access of beta rather than coef() avoids intercept
#   } else if(inherits(x, "ncvreg")){
#     ## ncvreg includes an intercept term (unlike glmnet) so we omit it
#     beta <- if (length(x$penalty.factor) == NROW(x$beta)){
#       coef(x)
#     } else {
#       coef(x)[-1, , drop = FALSE] ## (intercept)
#     }
#   }
#   lambda <- x$lambda
#   
#   active = glmnet::nonzeroCoef(beta)
#   nactive = length(active)
#   beta <- as.matrix(beta[active, , drop=FALSE])
#   index = log(lambda)
#   xlab=expression(log(lambda))
#   
#   ## Identify the first num_predictors to enter the model
#   ## (or fewer, in the case when we add two variables simultaneously)
#   key_predictors <- names(which(beta[,max(which(colSums(beta != 0) <= num_predictors))] != 0))
#   if((length(key_predictors) > num_predictors) | (length(key_predictors) < 2)){
#     ## If we have lots of predictors (e.g., ridge)
#     ## just pull out the predictors with the largest |\hat{\beta}_j|
#     ## and carry on
#     key_predictors <- rownames(beta)[rank(abs(beta[,ncol(beta)]))][1:num_predictors]
#   }
#   
#   
#   ## Now we rearrange our data in a 'long' form preferred by ggplot2
#   mb <- melt(beta)
#   colnames(mb) <- c("Variable", "Step", "Beta")
#   mb$Variable <- as.character(mb$Variable)
#   mb <- cbind(mb, key_predictor=mb$Variable %in% key_predictors)
#   mb <- cbind(mb, LogLambda=rep(index, each=nrow(beta)))
#   mb <- cbind(mb, plot_color=ifelse(mb$key_predictor, mb$Variable, "Other"))
#   
#   ## And add some extra columns to use as plot parameters.
#   g <- ggplot(mb, aes(x=LogLambda, y=Beta, group=Variable))
#   if(!all(mb$key_predictor)){
#     g <- g + geom_line(aes(col=plot_color, alpha=key_predictor))
#     g <- g + scale_alpha_discrete(aes(alpha=key_predictor), guide="none", range=c(alpha_min, 1))
#   } else {
#     g <- g + geom_line(aes(col=plot_color))
#   }
#   
#   g <- g + xlab(expression(log(lambda))) + ylab(expression(hat(beta[j])))
#   g <- g + scale_color_discrete(guide=guide_legend(title="Predictor"), 
#                                 breaks=c(mixedsort(unique(mb$Variable)), "Other"))
#   
#   ## Conventionally, these are plotted with decreasing LogLambda
#   g <- g + scale_x_reverse()
#   
#   g 
# }
# 
# 


####SVM
