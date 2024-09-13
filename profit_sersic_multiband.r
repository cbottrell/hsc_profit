library(ProFuse) #do all the things
library(ProPane) #stack things
library(ProFound) #detect things
library(ProFit) #model things
library(Rfits) #work with FITS
library(foreach) #parallelise things
library(doParallel) #parallelise things

registerDoParallel(cores=1)

img_name = Sys.getenv('IMG_NAME')
universe = Sys.getenv('UNIVERSE')
simulation = Sys.getenv('SIMULATION')
snapnum = Sys.getenv('SNAPNUM')
subfindid = Sys.getenv('SUBFINDID')
camera = Sys.getenv('CAMERA')
out_name = Sys.getenv('OUT_NAME')

hdulist = Rfits_read(img_name, header=TRUE, pointer=FALSE)
header = hdulist[['SUBARU_HSC.I']]$keyvalues

bands = list(
    'SUBARU_HSC.G',
    'SUBARU_HSC.R',
    'SUBARU_HSC.I',
    'SUBARU_HSC.Z',
    'SUBARU_HSC.Y'
)

# image_list = to_list(
#     for(band in bands) hdulist[[band]]$imDat)
# skyRMS_list = to_list(
#     for(band in bands) sqrt(hdulist[[paste(band,'VARIANCE')]],$imDat)
                      
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
                      
# ProFound regions
regions = profoundProFound(
    image = stack$image, sky = 0, skyRMS = stack$skyRMS,
    redosky = FALSE, box=62, boundstats=TRUE, nearstats=TRUE
)

all_fits = foreach(band = bands)%dopar%{ #loop over bands
    
    result = profuseDoFit(
        image = hdulist[[band]]$imDat, 
        sigma = sqrt(hdulist[[paste(band,'VARIANCE')]]$imDat),
        segim = regions$segim,
        Ncomp = 1, #just fit a free Sersic, can set to 1.5 and 2 for other version of this script, see ?profuseFound2Fit
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
                         
    return(result)
}

names(all_fits) = c(bands)

# Add the common segim, segID and the ProFound segstats run on the stacked image to the front of the list
# Note we still include the per band ProFound (since this will be different), not sure what you want

all_fits = c(list(segim = regions$segim),
             list(segID = regions$segim[dim(regions$segim)[1]/2, dim(regions$segim)[2]/2]),
             list(profound_stack = regions$segstats),
             all_fits
)

#Write out the FITS

Rfits_write(all_fits, filename = out_name, flatten = TRUE)



