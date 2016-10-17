#Utilizing free surfer pull script, we now initialize baseline + neuroimaging models for Anxiety and Depression
data_dir = '/Users/rhuang/Documents/Classes/Fall_16/Stat_Research/R/data'
script_dir = '/Users/rhuang/Documents/Classes/Fall_16/Stat_Research/R'

setwd(data_dir)
#reading in data
data <- read.csv("freesurfer_out.csv", header = T)[,-1] #free-surder
diag = load(file = 'menninger_clincial.rda')#diagnosis
clinical = get(diag); rm(dataset) 

#Formatting
data$MIND_ID <- as.numeric(substring(data$MIND_ID, 6, 9)) #format mind_ID to just the number so we can match with diagnosis data

#Begin experiment / Naming conventions fixed from here on out
test_ID <- read.csv('test_ID.csv')[,1] #load test ids saved randomly splitting
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

sum(anxiety_indicator) #224 of 369 with anxiety
sum(depression_indicator)  #250 of 369 with anxiety
sum((anxiety_indicator == 1)  & (depression_indicator == 1)) #156 of 369 have anxiety and depression

disease_dsm <- data.frame(anxiety = anxiety_indicator, depression = depression_indicator) #dsm9 grouped-disease matrix; depression and anxiety only for now

#Model creation


Xtr <- scale(data_mri[,-1], center = T, scale = T)
