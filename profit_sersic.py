from astropy.io import fits
from astropy.cosmology import Planck15 as cosmo
import numpy as np
import pymysql
import morph
import os, sys, time, copy
api_key = os.environ['TNG_API_KEY']
scratch_path = os.environ['MYSCRATCH']

universe = 'IllustrisTNG'
simulation = 'TNG50-1'

hostname = 'connor-db@icrar.org'
database = simulation.replace('-','_')
table = 'MorphologiesHSC_ProFit_Sersic'
db_commit=False

snapnum = 72 # 72 # 
subfindid = 360895 # 391507 # 391507 # 360895 
band = 'i'
camera = 'v0'

virgotng_path = f'{scratch_path}/Simulations/virgotng'
sim_path = f'{virgotng_path}/data/{universe}/{simulation}'
img_path = f'{sim_path}/postprocessing/skirt_images_hsc/realistic/{snapnum:03}'
if not os.access(img_path,0): 
    os.system(f'mkdir -p {img_path}')

# check if database entry exists
# update band list to fit/submit

rec = {
    'dbID': f'{snapnum}_{subfindid}',
    'SnapNum': snapnum,
    'SubfindID': subfindid,
    'Camera': camera,
    'Band': band,
    'ProcFlag': -1,
}

dbcmd = ' '.join([
    f'INSERT INTO {table} SET',
    f'dbID="{rec["dbID"]}",',
    f'SnapNum={rec["SnapNum"]},',
    f'SubfindID={rec["SubfindID"]},',
    f'Camera={rec["Camera"]},',
    f'Band={rec["Band"]},',
    f'ProcFlag={rec["ProcFlag"]}'
])

if db_commit:
    try:
        submit(dbcmd,database)
    except:
        # give warning and remove band from list
        print(f'DB {table} already contains {snapnum}-{subfindid}-{camera} in {band}-band. Skipping...\n')

            
img_name = f'{img_path}/shalo_{snapnum:03}-{subfindid}_{camera}_HSC_GRIZY.fits'
out_name = f'{img_path}/profuse_{snapnum:03}-{subfindid}_{camera}_{band}.fits'

if not os.access(out_name,0):

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
    os.environ['BAND'] = f'SUBARU_HSC.{band.capitalize()}'
    os.environ['OUT_NAME'] = out_name
    
    os.system('Rscript profit_sersic.r')