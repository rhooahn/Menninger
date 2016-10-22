require(randomForest)
set.seed(415)


initialize_data <- function(){
#Utilizing free surfer pull script, we now initialize baseline + neuroimaging models for Anxiety and Depression
data_dir = '/Users/kirylnovikau/Documents/stat444/jkstat450'
script_dir = '/Users/kirylnovikau/Documents/stat444/jkstat450'

setwd(data_dir)
#reading in data
data <- read.csv("freesurfer_out.csv", header = T)[,-1] #free-surfer
diag = load(file = 'menninger_clincial.rda')#diagnosis
clinical = get(diag); rm(dataset) 

#Formatting
data$MIND_ID <- as.numeric(substring(data$MIND_ID, 6, 9)) #format mind_ID to just the number so we can match with diagnosis data

#Begin experiment / Naming conventions fixed from here on out
# total_ID <- intersect(clinical$MIND_ID, data$MIND_ID)
# test_ID <- sample(total_ID, floor(1/6 * length(total_ID)))
# write.csv(test_ID, 'test_ID.csv')
test_ID <<- read.csv('test_ID.csv')[,2] #load test ids saved randomly splitting
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
        index <- c(unique(c(index1 + 2, index2 + 2)),1251,1252)
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

disease_dsm <<- data.frame(anxiety = anxiety_indicator, depression = depression_indicator) #dsm9 grouped-disease matrix; depression and anxiety only for now

#Demographic data pulling
demo_features <- c("GENDER", "AGE","PI3_1", "PIMARRIED_1", "PI4_1")
demo_df <- data_clinical[, which(colnames(data_clinical) %in% demo_features)]
sapply(3:5, function(x) demo_df[,x] <<- as.factor(demo_df[,x]))


#Data pre-processing for model input
na_mri <- which(colSums(data_mri) == 0) ; data_mri <- data_mri[,-na_mri]#remove NA from MRI data
na_demo <- which(rowSums(is.na(demo_df))>0); data_mri <- data_mri[-na_demo,]; demo_df <- demo_df[-na_demo,] #remove NA rows from demo

Xtr_demo <- model.matrix(~., demo_df)[,-1]
Xtr_full <<- cbind(Xtr_demo, data_mri[,-1])
Ytr <<- disease_raw[-na_demo,]
Ytr_bin <<- disease_dsm[-na_demo,]
}

initialize_data()

(rf_dep <- randomForest(x=Xtr_full, y=as.factor(Ytr_bin[,2]), importance = T, sampsize = 200, ntree=2000))
# Just predicting 1 for everything?
varImpPlot(rf_dep)

(rf_anx <- randomForest(x=Xtr_full, y=as.factor(Ytr_bin[,1]), importance = T, sampsize = 200, ntree=2000))
# Just predicting 1 for everything?
varImpPlot(rf_dep)

(rf_dep_reg <- randomForest(x=Xtr_full, y=Ytr[,45], importance = T, sampsize = 100, ntree=1000))
varImpPlot(rf_dep_reg)
summary(rf_dep_reg)

(rf_anx_reg <- randomForest(x=Xtr_full, y=Ytr[,46], importance = T, sampsize = 100, ntree=1000))
varImpPlot(rf_anx_reg)
summary(rf_anx_reg)
# Both of the regression models have a -% Var explained 
# so we would be better off just predicting the mean

ctrl = tune.control(sampling = c("fix"))
svm_dep <- best.svm(y = as.factor(Ytr_bin[,2]), x = Xtr_full, type = "C-classification", tunecontrol = ctrl)
summary(svm_dep)
# 303 Support vectors? Is this even predicting anything?

svm_anx <- best.svm(y = as.factor(Ytr_bin[,1]), x = Xtr_full, type = "C-classification", tunecontrol = ctrl)
summary(svm_anx)
# 304 Support vectors?

