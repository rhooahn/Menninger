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
  file_name <- files[x]
  sample_cur <- grid[x,] #index 1 is hemisphere, index 2 is measure
  
  aparcstats2table(subjects = files, outfile = file_name, sep = c("comma"), 
                   hemi = sample_cur[1], measure = sample_cur[2],  parc = "aparc",  subj_dir = '/backup/crate/menninger/FREE_SURF/',  
                   skip = TRUE)
})


####################################################################################################
#### do this for all combinations of measure = c("volume", "mean", "std")
#FOR SUBCORTICAL
####################################################################################################

asegstats2table(subjects = files, outfile = 'youroutfile.txt',,
      measure = c("std"), sep = c("comma"), 
      skip = TRUE, subj_dir = '/backup/crate/menninger/FREE_SURF/', 
      verbose = TRUE)




