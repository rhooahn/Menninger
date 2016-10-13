library(freesurfer)

files <- dir('/backup/crate/menninger/FREE_SURF/', pattern = 'Mind')

####################################################################################################
#### do this for all combinations of ("lh", "rh") and ("area", "volume", "thickness", "thicknessstd", 
#### "meancurv", "gauscurv", "foldind", "curvind")
#FOR CORTICAL
####################################################################################################
combos <- list(hemi = c("lh", "rh"), measure = c("area", "volume", "thickness", "thicknessstd", 
                                       "meancurv", "gauscurv", "foldind", "curvind"))
grid <- expand.grid(combos)
grid  <- data.frame(lapply(grid, as.character), stringsAsFactors=FALSE) #converts cols to char format

sample_dir <- "/home/rbh1/samples"
output_dir <- "/home/rbh1/output"

setwd(sample_dir)

sapply(1:nrow(grid), function(x) {
  sample_cur <- grid[x,] #index 1 is hemisphere, index 2 is measure
  file_name <- paste(sample_cur[[1]],"_",sample_cur[[2]],"_aparc", sep = '')
  
  aparcstats2table(subjects = files, outfile = file_name, sep = c("comma"), 
                   hemi = sample_cur[[1]], measure = sample_cur[[2]],  parc = "aparc",  subj_dir = '/backup/crate/menninger/FREE_SURF/',  
                   skip = TRUE)
})


####################################################################################################
#### do this for all combinations of measure = c("volume", "mean", "std")
#FOR SUBCORTICAL
####################################################################################################

grid2 <- expand.grid(c("volume", "mean", "std"))
grid2  <- data.frame(lapply(grid2, as.character), stringsAsFactors=FALSE) #converts cols to char format

sapply(1:nrow(grid2), function(x) {
  sample_cur <- grid2[x,] #index 1 is hemisphere, index 2 is measure
  file_name <- paste(sample_cur[[1]],"_subc", sep = '')
  
  asegstats2table(subjects = files, outfile = file_name,,
                  measure = sample_cur[[1]], sep = c("comma"), 
                  skip = TRUE, subj_dir = '/backup/crate/menninger/FREE_SURF/', 
                  verbose = TRUE)
})

#Have all individual samples, now merge them to one df

sample_files <- dir(sample_dir)
merged_df <- read.csv(sample_files[1]); colnames(merged_df)[1] <- "MIND_ID"
i <- sapply(merged_df, is.factor)
merged_df[i] <- lapply(merged_df[i], as.character)

sapply(2:length(sample_files), function(i) {
  new_df <-  read.csv(sample_files[i])
  colnames(new_df)[1] <- "MIND_ID"
  
  j <- sapply(new_df, is.factor)
  new_df[j] <- lapply(new_df[j], as.character)
merged_df <<- merge(merged_df, new_df, by= "MIND_ID") 
  })

setwd(output_dir)
write.csv(merged_df, 'freesurfer_out.csv')
