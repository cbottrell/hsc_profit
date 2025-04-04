import os,sys
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from astropy.io import fits

def _hsccolour(rgb,ref_mag=26):
    rgb *= 10**(0.4*(22.5-ref_mag))
    u_min = -0.05
    u_max = 2. / 3.
    u_a = np.exp(10.)
    for i, x in enumerate(rgb):
        if i==2: x*=2 # make blue more blue
        x = np.arcsinh(u_a*x) / np.arcsinh(u_a)
        x = (x - u_min) / (u_max - u_min)
        rgb[i] = x
    return rgb

def colour_image(
    image_path, axes, 
    universe, simulation, snapnum, subfindid, camera, 
    bands='gri', level = 'realistic', # or 'idealized'
    zoom_factor = 1.,
):

    # file path (for idealized / realistic)
    if level=='idealized':
        filename = f'{image_path}/shalo_{snapnum:03}-{subfindid}_{camera}_photo.fits'
    else:
        filename = f'{image_path}/shalo_{snapnum:03}-{subfindid}_{camera}_HSC_GRIZY.fits'
        
    if not os.access(filename,0):
        sys.exit(f'File {filename} does not exist.')
        
    # filter names in FITS extension format    
    filters = np.array([f'Subaru_HSC.{band}' for band in bands.upper()])[::-1]

    # rgb cube and reference magnitudes (adjust as you like)
    rgb = []
    
    with fits.open(filename,mode='readonly') as hdul:
        header = hdul[0].header
        npix = header['NAXIS1']
        redshift = header['REDSHIFT']
        kpc_per_pixel = header['CDELT1']
        fov_kpc = header['FOVSIZE']
        for filt in filters:
            if level=='idealized':
                rgb+=[10**(0.4*(22.5-hdul[filt].data))]
                ref_mag=31
            else:
                rgb+=[hdul[filt].data]
                ref_mag=26
            
    # set buffer for zoom-in
    buf = int(npix*(1-1/zoom_factor)/2)
    rgb = np.array(rgb)
    rgb = _hsccolour(rgb,ref_mag).transpose(1,2,0)[
        buf:npix-buf,buf:npix-buf
    ]
    
    axes.imshow(rgb,origin='lower',interpolation='None',aspect='auto')
    return axes