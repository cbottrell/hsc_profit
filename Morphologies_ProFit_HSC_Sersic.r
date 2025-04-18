library(ProFuse) #do all the things
library(ProPane) #stack things
library(ProFound) #detect things
library(ProFit) #model things
library(Rfits) #work with FITS
library(foreach) #parallelise things

img_name = Sys.getenv('IMG_NAME')
universe = Sys.getenv('UNIVERSE')
simulation = Sys.getenv('SIMULATION')
snapnum = Sys.getenv('SNAPNUM')
subfindid = Sys.getenv('SUBFINDID')
camera = Sys.getenv('CAMERA')
band = Sys.getenv('BAND')
out_name = Sys.getenv('OUT_NAME')

hdulist = Rfits_read(img_name, header=TRUE, pointer=FALSE)
header = hdulist[[band]]$keyvalues
                      
# Stack bands for ProFound
stack = propaneStackFlatInVar(
    image_list = list(
        hdulist[['SUBARU_HSC.G']]$imDat, 
        hdulist[['SUBARU_HSC.R']]$imDat, 
        hdulist[['SUBARU_HSC.I']]$imDat, 
        hdulist[['SUBARU_HSC.Z']]$imDat, 
        hdulist[['SUBARU_HSC.Y']]$imDat
    ),
    sky_list = c(0,0,0,0,0),
    skyRMS_list = list(
        sqrt(hdulist[['SUBARU_HSC.G VARIANCE']]$imDat), 
        sqrt(hdulist[['SUBARU_HSC.R VARIANCE']]$imDat), 
        sqrt(hdulist[['SUBARU_HSC.I VARIANCE']]$imDat), 
        sqrt(hdulist[['SUBARU_HSC.Z VARIANCE']]$imDat), 
        sqrt(hdulist[['SUBARU_HSC.Y VARIANCE']]$imDat)
    ),
    magzero_in = c(22.5,22.5,22.5,22.5,22.5),
    magzero_out = 22.5)
                      
# ProFound segmantation and measurements using stacked image
regions_stack = profoundProFound(
    image=stack$image, sky=0, skyRMS=stack$skyRMS, skycut=1.5,
    redosky=FALSE, box=c(64,64), boundstats=TRUE, nearstats=TRUE,
    groupstats=TRUE, tolerance = 15, reltol = -10, watershed='ProFound',
    pixcut=10, sigma=1, SBdilate=2, verbose=FALSE, keepsegims=TRUE,
    dotot=TRUE, docol=TRUE
)

# Compute static g-Y colours using stacked segentation
regions_blue = profoundProFound(
    image = hdulist[['SUBARU_HSC.R']]$imDat,
    segim = regions_stack$segim, static_photom = TRUE,
)
regions_red = profoundProFound(
    image = hdulist[['SUBARU_HSC.Z']]$imDat,
    segim = regions_stack$segim, static_photom = TRUE,
)

# Add g-Y colour as new column in segstats for AutoMerge
segstats_tmp = regions_stack$segstats
segstats_tmp$col <- regions_blue$segstats$mag - regions_red$segstats$mag

groups = profoundAutoMerge(
    segim = regions_stack$segim_orig, segstats_tmp, 
    spur_lim = 2.0e-3, col_lim = c(0,0.8), Ncut = 1
)
segim_merged = profoundSegimKeep(
    segim=regions_stack$segim, segID_merge=groups$segID
)

regions_stack = profoundProFound(    
    image = stack$image,
    segim = segim_merged, static_photom = TRUE,
)

#regions_stack$segim = segim_merged

#all_fits = foreach(band = bands)%dopar%{ #loop over bands
    
result = profuseDoFit(
    image = hdulist[[band]]$imDat, 
    sigma = sqrt(hdulist[[paste(band,'VARIANCE')]]$imDat),
    segim = regions_stack$segim,
    boundstats=TRUE, rotstats=TRUE, 
    Ncomp = 1, nser_upper=8.,
    psf = hdulist[[paste(band,'PSF')]]$imDat,
    magzero = 22.5,
    tightcrop = FALSE, #so positions are relative to initial image
    seed = 666 #the seed of the beast
) 

model_im = profitMakeModel(
    result$finalmodel$modellist, 
    dim=result$Data$imagedim, 
    psf=result$Data$psf, 
    magzero = 22.5
)$z
                         
model_res = result$Data$image - model_im

result = list(
    best = result$parm,
    error = result$error,
    monitor = result$LD_last$Monitor,
    model_im = model_im,
    model_res = model_res,
    posterior = as.data.frame(result$LD_last$Posterior1),
    profound = result$profound$segstats
)

all_fits = list(result)
names(all_fits) = c(list(band))

# Add the common segim, segID and the ProFound segstats run on the stacked image to the front of the list
# Note we still include the per band ProFound (since this will be different), not sure what you want

all_fits = c(list(segim = regions_stack$segim),
             list(segID = regions_stack$segim[dim(regions_stack$segim)[1]/2, dim(regions_stack$segim)[2]/2]),
             list(profound_stack = regions_stack$segstats),
             all_fits
)

#Write out the FITS

Rfits_write(all_fits, filename = out_name, flatten = TRUE)



