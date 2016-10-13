#Module Review
setwd("~/Documents/Classes/Fall_16/Stat_Research/Neurohacking_data/BRAINIX/DICOM/FLAIR")
require("oro.dicom")

slice = readDICOM("IM-0001-0011.dcm")
str(slice)
str(slice$img[1])

#DICOM images 
d = dim(slice$img[[1]])
#x,y locates the pixel in Euclidean space, z is a color mapping of each point.
image(1:d[1], 1:d[2], t(slice$img[[1]]), col = gray(0:64/64))

#histogram
hist(slice$img[[1]], breaks = 50, prob = T)
hdr = slice$hdr[[1]]
hdr[hdr$name == "PixelSpacing",]

#Full image
setwd("~/Documents/Classes/Fall_16/Stat_Research/Neurohacking_data/BRAINIX/DICOM")


full_image = readDICOM("T1/")

dim(full_image$img[[11]])

#Nifti format
nif_T1 = dicom2nifti(full_image)
d = dim(nif_T1)
d
#11th slice
image(1:d[1], 1:d[2], nif_T1[,,11], col = gray(0:64/64))

###NIFTI
setwd("~/Documents/Classes/Fall_16/Stat_Research/Neurohacking_data/BRAINIX/NIfTI")
fname = "Output_3D_File"
print({nii_T1 = nif_T1})

#imaging
image(1:d[1], 1:d[2], nif_T1[,,11], col = gray(0:64/64))

image(nif_T1, z=22, plot.type = "single")
length(nif_T1)
512*512*22


str(nif_T1[,,11])
str(nif_T1)
#Module Review
