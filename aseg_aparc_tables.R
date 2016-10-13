library(freesurfer)

files <- dir('/backup/crate/menninger/FREE_SURF/', pattern = 'Mind')

####################################################################################################
#### do this for all combinations of ("lh", "rh") and ("area", "volume", "thickness", "thicknessstd", 
#### "meancurv", "gauscurv", "foldind", "curvind")
#FOR CORTICAL
####################################################################################################

aparcstats2table(subjects = files, outfile = 'youroutfile.txt', sep = c("comma"), 
hemi = "lh", measure = "area",  parc = "aparc",  subj_dir = '/backup/crate/menninger/FREE_SURF/',  
skip = TRUE) #skips over non-finished output


####################################################################################################
#### do this for all combinations of measure = c("volume", "mean", "std")
#FOR SUBCORTICAL
####################################################################################################

asegstats2table(subjects = files, outfile = 'youroutfile.txt',,
      measure = c("std"), sep = c("comma"), 
      skip = TRUE, subj_dir = '/backup/crate/menninger/FREE_SURF/', 
      verbose = TRUE)




