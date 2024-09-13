from astropy.io import fits
import numpy as np
import pymysql
import morph
import os, sys, time
api_key = os.environ['TNG_API_KEY']
scratch_path = os.environ['MYSCRATCH']

universe = 'IllustrisTNG'
simulation = 'TNG50-1'
snapnum = 72
subfindid = 360895

virgotng_path = f'{scratch_path}/Simulations/virgotng'
sim_path = f'{virgotng_path}/data/{universe}/{simulation}'
img_path = f'{sim_path}/postprocessing/skirt_images_hsc/realistic/{snapnum:03}'
if not os.access(img_path,0): 
    os.system(f'mkdir -p {img_path}')

for camera in ['v0','v1','v2','v3']:

    out_name = f'{img_path}/profuse_{snapnum:03}-{subfindid}_{camera}_grizy.fits'
    if os.access(out_name,0):
        continue
    
    img_name = f'{img_path}/shalo_{snapnum:03}-{subfindid}_{camera}_HSC_GRIZY.fits'
    if not os.access(img_name,0):
        base_url = f'http://www.tng-project.org/api'
        img_url = f'{base_url}/{simulation}/snapshots/{snapnum}/subhalos/{subfindid}/skirt/skirt_images_hsc_realistic_{camera}.fits'
        sys_cmd = f'wget -nc -nv --content-disposition --header="API-Key:{api_key}" "{img_url}" -P {img_path}'
        os.system(sys_cmd)
    
    # Pass target image to R
    os.environ['IMG_NAME'] = img_name
    os.environ['UNIVERSE'] = universe
    os.environ['SIMULATION'] = simulation
    os.environ['SNAPNUM'] = str(snapnum)
    os.environ['SUBFINDID'] = str(subfindid)
    os.environ['CAMERA'] = camera
    os.environ['OUT_NAME'] = out_name
    
    os.system('Rscript profit_sersic_multiband.r')