#Utilizing free surfer pull script, we now initialize baseline + neuroimaging models for Anxiety and Depression
data_dir = '/Users/rhuang/Documents/Classes/Fall_16/Stat_Research/R/data' #freesurfer_out.csv
script_dir = '/Users/rhuang/Documents/Classes/Fall_16/Stat_Research/R' #directory of this Rscript

##### 
setwd(data_dir)
#reading in data
data <- read.csv("freesurfer_out.csv", header = T)[,-1] #free-surfer
diag = load(file = 'menninger_clincial.rda')#Clinical data should be in data directory
clinical = get(diag); rm(dataset) 

#Formatting
data$MIND_ID <- as.numeric(substring(data$MIND_ID, 6, 9)) #format mind_ID to just the number so we can match with diagnosis data

#Begin experiment / Naming conventions fixed from here on out
# total_ID <- intersect(clinical$MIND_ID, data$MIND_ID)
# test_ID <- sample(total_ID, floor(1/6 * length(total_ID)))
# write.csv(test_ID, 'test_ID.csv')
test_ID <- read.csv('test_ID.csv')[,2] #load test ids saved randomly splitting
data_mri = subset(data, MIND_ID %in% clinical$MIND_ID & !(MIND_ID %in% test_ID)) #subsets to the patients that have both clinical and neuroimaging results
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
  
depression_indicator <- apply(disease_raw[,depression_dsm_names],1, sum, na.rm = TRUE)
depression_indicator[depression_indicator>0] <- 1

anxiety_indicator <- apply(disease_raw[,anxiety_dsm_names],1, sum, na.rm = TRUE)
anxiety_indicator[anxiety_indicator>0] <- 1

sum(anxiety_indicator) #190 of 369 with anxiety
sum(depression_indicator)  #208 of 369 with anxiety
sum((anxiety_indicator == 1)  & (depression_indicator == 1)) #132 of 369 have anxiety and depression

disease_dsm <- data.frame(anxiety = anxiety_indicator, depression = depression_indicator) #dsm9 grouped-disease matrix; depression and anxiety only for now

#Demographic data pulling
demo_features <- c("GENDER", "AGE","PI3_1", "PIMARRIED_1", "PI4_1")
demo_df <- data_clinical[, which(colnames(data_clinical) %in% demo_features)]
sapply(3:5, function(x) demo_df[,x] <<- as.factor(demo_df[,x]))


#Data pre-processing for model input
na_mri <- which(colSums(data_mri) == 0) ; data_mri <- data_mri[,-na_mri]#remove NA from MRI data
na_demo <- which(rowSums(is.na(demo_df))>0); data_mri <- data_mri[-na_demo,]; demo_df <- demo_df[-na_demo,] #remove NA rows from demo

demo_df$AGE <- scale(demo_df$AGE, center = T, scale = T)

#####
#initialize libraries
library(glmnet)
#Depression models
Xtr_demo <- model.matrix(~., demo_df)[,-1]
Xtr_full <- cbind(Xtr_demo, scale(data_mri[,-1], center = T, scale = T))
Ytr <- disease_dsm[-na_demo,2]
##### Data exploration/summary
summary(Xtr_full) #only finding is that there are some extreme outliers
data_mri$WM.hypointensities
str(Xtr_full)
#####
nrow(Xtr_full)
base_lam <-cv.glmnet(Xtr_demo, as.factor(Ytr), family = "binomial", #only penalize the continuous variables
          penalty.factor = c(0,1,rep(0,16)), standardize = FALSE, alpha = .5)$lambda.min

basefit <- glmnet(Xtr_demo, Ytr, family = "binomial", #only penalize the continuous variables
       penalty.factor = c(0,1,rep(0,16)), standardize = FALSE, lambda = base_lam, alpha = 1)

cv_full <- cv.glmnet(Xtr_full, as.factor(Ytr),
                     family = "binomial",standardize = FALSE, alpha = .5, 
                     type.measure = 'class',penalty.factor = c(0,1,rep(0,16), rep(1,692)))

plot(cv_full$lambda[1:10], cv_full$cvm[1:10])
cv_full$lambda.min
fullfit <- glmnet(Xtr_full, as.factor(Ytr),family = "binomial", standardize = FALSE,
                  lambda = cv_full$lambda.min,penalty.factor = c(0,1,rep(0,16), rep(1,692)), alpha = .5)
coef(fullfit)

g_full_depression <- glmnet(Xtr_full, as.factor(Ytr),family = "binomial", standardize = FALSE, alpha = .5)


key_predictor_plot(g_full_depression, alpha_min = .1)+ggtitle("Elastic net (alpha .5) regularization path - Binary Depression Model")







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