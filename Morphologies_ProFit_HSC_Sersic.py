from astropy.io import fits
from astropy.cosmology import Planck15 as cosmo
import numpy as np
import mymysql
import morph
import os, sys, time, copy, itertools

# global paths/credentials for all jobs
api_key = os.environ['TNG_API_KEY']
scratch_path = os.environ['MYSCRATCH']
software_path = os.environ['MYSOFTWARE']
cnf_path = '~/.mysql/ningaloo.cnf'

def ProFit_HSC_Sersic(
    universe, simulation, snapnum, subfindid,
    camera, band, img_path, 
    db_commit, database, table,
):

    # Primary key for database table
    databaseid = f'{snapnum:03}.{subfindid}.{camera}.{band}'

    # Initialization command for record on database
    dbcmd = ' '.join([
        f'INSERT INTO {table} SET',
        f'DatabaseID="{databaseid}",',
        f'SnapNum={snapnum},',
        f'SubfindID={subfindid},',
        f'Camera="{camera}",',
        f'Band="{band}",',
        f'ProcessFlag=-1'
    ])

    # Check if database entry exists, add new entry if not
    if db_commit:
        try:
            mymysql.submit(dbcmd,database=database,cnf_path=cnf_path)
        except:
            print(f'Table {table} already contains {databaseid}. Skipping...\n')
            return
                
    img_name = f'{img_path}/shalo_{snapnum:03}-{subfindid}_{camera}_HSC_GRIZY.fits'
    out_name = f'{img_path}/profuse_{snapnum:03}-{subfindid}_{camera}_{band}.fits'

    record = {'ProcessFlag': 0,}

    # Download image if not already available
    if not os.access(img_name,0):
        base_url = f'http://www.tng-project.org/api'
        img_url = f'{base_url}/{simulation}/snapshots/{snapnum}/subhalos/{subfindid}/skirt/skirt_images_hsc_realistic_{camera}.fits'
        sys_cmd = f'wget -nc -nv --content-disposition --header="API-Key:{api_key}" "{img_url}" -P {img_path}'
        os.system(sys_cmd)

    # Check that image is complete after download
    try: 
        hdul_in = fits.open(img_name,mode='readonly')
        header = hdul_in[f'SUBARU_HSC.{band}'].header  
        img = hdul_in[f'SUBARU_HSC.{band}'].data
        var = hdul_in[f'SUBARU_HSC.{band} VARIANCE'].data
        psf = hdul_in[f'SUBARU_HSC.{band} PSF'].data
    except:
        # Flag image download as unsuccessful
        record['ProcessFlag']=2**3
        if db_commit:
            dbcmd = mysql_table_update_cmd(record, table, databaseid)
            mymysql.submit(dbcmd, database=database, cnf_path=cnf_path)
        return  

    record['Redshift'] = float(header["redshift"])
    record['ApparentMagnitude'] = float(header["apmag"])
    record['RightAscension'] = float(header["ra"])
    record['Declination'] = float(header["dec"])

    # Run Pro-tools on image in R
    if not os.access(out_name,0):
        os.environ['IMG_NAME'] = img_name
        os.environ['UNIVERSE'] = universe
        os.environ['SIMULATION'] = simulation
        os.environ['SNAPNUM'] = str(snapnum)
        os.environ['SUBFINDID'] = str(subfindid)
        os.environ['CAMERA'] = camera
        os.environ['BAND'] = f'SUBARU_HSC.{band.capitalize()}'
        os.environ['OUT_NAME'] = out_name
        os.system('Rscript Morphologies_ProFit_HSC_Sersic.r')

    try: 
        hdul_out = fits.open(out_name, mode='readonly')
    except:
        # Flag if parametric fitting was unsuccessful
        record['ProcessFlag']+=2**0
        if db_commit: 
            dbcmd = mysql_table_update_cmd(record, table, databaseid)
            mymysql.submit(dbcmd, database=database, cnf_path=cnf_path)
        return
    
    # Long form column names (e.g. sersic.nser), convert to short
    columns = hdul_out[f'SUBARU_HSC.{band}.posterior'].columns.names
    cols = [col.replace('sersic.','') for col in columns]
    
    # maximum likelihood params
    r = hdul_out[f'SUBARU_HSC.{band}.best'].data
    record['Sersic_xcen'] = r[0]
    record['Sersic_ycen'] = r[1]
    record['Sersic_mag'] = r[2]
    record['Sersic_re'] = 10**r[3]
    record['Sersic_nser'] = 10**r[4]
    record['Sersic_ang'] = r[5]
    record['Sersic_axrat'] = 10**r[6]
    
    # posterior chain params
    r = hdul_out[f'SUBARU_HSC.{band}.posterior'].data
    
    for i,col in enumerate(cols):
        # these columns are in log units
        if col in ['re','nser','axrat']:
            chain = 10**r[f'sersic.{col}']
        else:
            chain = r[f'sersic.{col}']
        record[f'Sersic_{col}_med'] = np.median(chain)
        record[f'Sersic_{col}_std'] = np.std(chain,ddof=1)
        record[f'Sersic_{col}_p84'] = np.percentile(chain,84)-record[f'Sersic_{col}_med']
        record[f'Sersic_{col}_m16'] = np.percentile(chain,16)-record[f'Sersic_{col}_med']

    segim = hdul_out['segim'].data
    segid = hdul_out['segid'].data[0]
    model = hdul_out[f'SUBARU_HSC.{band}.model_im'].data
    residual = img-model
    
    # compute chi2nu for sersic fit (model dof = 7)
    model_ndof = 7
    mask = (segim==segid) | (segim==0)
    chi2nu = np.sum(residual[mask]**2/var[mask])/(np.sum(mask)-model_ndof)
    record['ReducedChiSquared'] = chi2nu
    
    # profound data (segmap, target id, and measurements)
    segim_primary = segim==segid # target object
    segim_secondary = (segim!=0) & (segim!=segid) # other non-sky objects
    cols = hdul_out[f'SUBARU_HSC.{band}.profound'].columns.names
    r = hdul_out[f'SUBARU_HSC.{band}.profound'].data
    rownum = np.argwhere(r['segID']==segid)[0][0]
    r = r[rownum]
    for col in cols:
        if col not in ['RAcen','Deccen','RAmax','Decmax']:
            record[f'ProFound_{col}'] = r[col]

    moments = morph.image_moments(img, segim, segid)
    pixscale = header['FOVARC']/header['NAXIS1']
    redshift = header['REDSHIFT']

    try:
        # Non-parametric morphologies using profound segmap
        npmorph = morph.nonparametric(
            data = img, variance = var, segmentation = segim, model = model, 
            residual = residual, primary_mask = segim_primary, secondary_mask = segim_secondary,
            moments = moments, psf = psf, pixel_scale = pixscale, redshift = redshift, cosmology = cosmo
        )
    except:
        # Flag that non-parametric calculations were unsuccessful
        record['ProcessFlag']+=2**1
        if db_commit:
            dbcmd = mysql_table_update_cmd(record, table, databaseid)
            mymysql.submit(dbcmd,database=database,cnf_path=cnf_path)
        return
    
    # elliptical and circular sizes
    record['RPetro_Elliptical'] = npmorph.rp_ellipse
    record['RPetro_Circular'] = npmorph.rp_circle
    record['RMax_Elliptical'] = npmorph.rmax_ellipse
    record['RMax_Circular'] = npmorph.rmax_circle
    record['R80_Elliptical'] = npmorph.r80_ellipse
    record['R80_Circular'] = npmorph.r80_circle
    record['R50_Elliptical'] = npmorph.r50_ellipse
    record['R50_Circular'] = npmorph.r50_circle
    record['R20_Elliptical'] = npmorph.r20_ellipse
    record['R20_Circular'] = npmorph.r20_circle
    
    # quantities derived from the image / segim only
    record['Asymmetry_xcen'] = npmorph.asymmetry_center[0]
    record['Asymmetry_ycen'] = npmorph.asymmetry_center[1]
    record['Asymmetry'] = npmorph.asymmetry_cas
    record['AsymmetryBkgDens'] = npmorph.bkg_asymmetry_density
    record['RMSAsymmetrySquared'] = npmorph.asymmetry_rms2
    record['OuterAsymmetry'] = npmorph.asymmetry_outer
    record['ShapeAsymmetry'] = npmorph.asymmetry_shape
    record['AsymmetryNoAperture_xcen'] = npmorph.asymmetry_center_nap[0]
    record['AsymmetryNoAperture_ycen'] = npmorph.asymmetry_center_nap[1]
    record['AsymmetryNoAperture'] = npmorph.asymmetry_cas_nap
    record['Concentration_Elliptical'] = npmorph.concentration_ellipse
    record['Concentration_Circular'] = npmorph.concentration_circle
    record['Smoothness'] = npmorph.smoothness
    record['SmoothnessBkgDens'] = npmorph.bkg_smoothness_density
    record['Gini'] = npmorph.gini
    record['M20'] = npmorph.m20
    record['SB1kpc'] = npmorph.sb_1kpc
    
    # model-dependent quantities computed from residuals
    record['ResidualAsymmetry_xcen'] = npmorph.asymmetry_residual_center[0]
    record['ResidualAsymmetry_ycen'] = npmorph.asymmetry_residual_center[1]
    record['ResidualAsymmetry'] = npmorph.asymmetry_residual
    record['ResidualAsymmetryNoAperture_xcen'] = npmorph.asymmetry_residual_center_nap[0]
    record['ResidualAsymmetryNoAperture_ycen'] = npmorph.asymmetry_residual_center_nap[1]
    record['ResidualAsymmetryNoAperture'] = npmorph.asymmetry_residual_nap

    if db_commit:
        try:
            dbcmd = mysql_table_update_cmd(record, table, databaseid)
            mymysql.submit(dbcmd,database=database,cnf_path=cnf_path)
        except:
            print(f'Could not add full record for {databaseid} to {table}. Skipping...\n')

    return
        

def mysql_table_update_cmd(
    record, table, databaseid,
):
    db_update = [
        f'UPDATE {table} SET'
    ]
    
    for key in record.keys():
        if type(record[key]) is str:
            db_update += [f'{key}="{record[key]}"']
        else:
            db_update += [f'{key}={record[key]}']
    
        if key != list(record.keys())[-1]:
            db_update += [',']
            
    db_update += [
        f'WHERE DatabaseID="{databaseid}"'
    ]
    return(' '.join(db_update))


def main_array():

    # Simulation identifiers
    universe = os.environ['UNIVERSE'] #'IllustrisTNG'
    simulation = os.environ['SIMULATION'] #'TNG50-1'

    # Mysql table identifiers
    database = os.environ['DATABASE'] # 'IllustrisTNG50_1'
    table = os.environ['TABLE'] # 'Morphologies_ProFit_HSC_Sersic']
    # Commit results to database?
    db_commit=True

    # Read in job array size, id from environment
    job_array_id = int(os.environ['JOB_ARRAY_INDEX'])
    job_array_size = int(os.environ['JOB_ARRAY_SIZE'])

    # Set limits of snapnum, choose subhalos
    snapmin,snapmax = 72,91
    snapnums = np.arange(snapmin,snapmax+1)
    logmstar_min = 9.
    logmstar_max = 99.

    # Camera angles and bands (independent tasks)
    cameras = ['v0','v1','v2','v3']
    bands = ['g','r','i','z','y']

    # Loop over snapshots, setting tasks for each job
    for snapnum in snapnums:

        # Set up paths to image directories
        virgotng_path = f'{scratch_path}/Simulations/virgotng'
        sim_path = f'{virgotng_path}/data/{universe}/{simulation}'
        img_path = f'{sim_path}/postprocessing/skirt_images_hsc/realistic/{snapnum:03}'
        if not os.access(img_path,0): 
            os.system(f'mkdir -p {img_path}')
        
        dbcmd = [
            f'SELECT SubfindID from Subhalos',
            f'WHERE SnapNum={snapnum}',
            f'AND SubhaloMassType_stars>={logmstar_min}',
            f'AND SubhaloMassType_stars<{logmstar_max}',
        ]
        subfindids = mymysql.query(
            ' '.join(dbcmd), database=database, cnf_path=cnf_path
        ).flatten()

        # All processing tasks for this snapnum
        iterables = np.array(list(itertools.product(subfindids,cameras,bands)))
        # Processing tasks for this job
        iterables = iterables[job_array_id::job_array_size]

        # Loop over tasks assigned to this job
        for iterable in iterables:
            subfindid = int(iterable[0])
            camera = iterable[1]
            band = iterable[2]

            ProFit_HSC_Sersic(
                universe=universe, simulation=simulation, snapnum=snapnum, subfindid=subfindid,
                camera=camera, band=band, img_path=img_path, db_commit=db_commit, 
                database=database, table=table,
            )


def main_mpi():

    from mpi4py import MPI

    # Read in job array size, id from environment
    # used to partition data for mpi processes
    array_size = int(os.environ['JOB_ARRAY_SIZE'])
    array_rank = int(os.environ['JOB_ARRAY_INDEX'])

    # For each array job, get mpi size and rank
    mpi_size = MPI.COMM_WORLD.Get_size()
    mpi_rank = MPI.COMM_WORLD.Get_rank()
    mpi_name = MPI.Get_processor_name()

    # Stagger jobs to avoid overloading SQL server
    time.sleep(array_rank*mpi_rank)

    # Simulation identifiers
    universe = os.environ['UNIVERSE'] 
    simulation = os.environ['SIMULATION'] 

    # Mysql table identifiers
    database = os.environ['DATABASE'] 
    table = os.environ['TABLE'] 
    # Commit results to database?
    db_commit=True

    # Set limits of snapnum, choose subhalos
    snapmin = 72
    snapmax = 91
    logmstar_min = 9.
    logmstar_max = 99.

    # Camera angles and bands (independent tasks)
    cameras = ['v0','v1','v2','v3']
    bands = ['g','r','i','z','y']
    
    dbcmd = ' '.join([
        f'SELECT SnapNum,SubfindID from Subhalos',
        f'WHERE SnapNum>={snapmin} and SnapNum<={snapmax}',
        f'AND SubhaloMassType_stars>={logmstar_min}',
        f'AND SubhaloMassType_stars<{logmstar_max}',
        f'ORDER BY snapnum,subfindid'
    ])
    df = mymysql.query_df(dbcmd, database=database, cnf_path=cnf_path)
    indices = np.array(df.index)
    
    # Create set of iterables for this job in array
    iters = np.array(list(itertools.product(indices,cameras,bands)))
    array_iters = iters[array_rank::array_size]

    # Get iterables for this specific mpi process 
    mpi_iters = array_iters[mpi_rank::mpi_size]
        
    # Loop over tasks assigned to this mpi process
    for iterable in mpi_iters:
        
        index = int(iterable[0])
        row = df.loc[index]
        snapnum = row['SnapNum']
        subfindid = row['SubfindID']
        camera = iterable[1]
        band = iterable[2]
    
        # Set up paths to image directories
        virgotng_path = f'{scratch_path}/Simulations/virgotng'
        sim_path = f'{virgotng_path}/data/{universe}/{simulation}'
        img_path = f'{sim_path}/postprocessing/skirt_images_hsc/realistic/{snapnum:03}'
        
        if not os.access(img_path,0): 
            os.system(f'mkdir -p {img_path}')
            
        ProFit_HSC_Sersic(
            universe=universe, simulation=simulation, snapnum=snapnum, subfindid=subfindid,
            camera=camera, band=band, img_path=img_path, db_commit=db_commit, 
            database=database, table=table,
        )

    
def main_single():

    # Simulation identifiers
    universe = 'IllustrisTNG'
    simulation = 'TNG50-1'

    # Mysql table identifiers
    database = 'IllustrisTNG50_1'
    table = 'Morphologies_ProFit_HSC_Sersic'
    # Commit results to database?
    db_commit=False

    snapnum=73
    subfindid=479757
    camera='v2'
    band='g'

    # Set up paths to image directories
    virgotng_path = f'{scratch_path}/Simulations/virgotng'
    sim_path = f'{virgotng_path}/data/{universe}/{simulation}'
    img_path = f'{sim_path}/postprocessing/skirt_images_hsc/realistic/{snapnum:03}'
    if not os.access(img_path,0): 
        os.system(f'mkdir -p {img_path}')

    ProFit_HSC_Sersic(
        universe=universe, simulation=simulation, snapnum=snapnum, subfindid=subfindid,
        camera=camera, band=band, img_path=img_path, db_commit=db_commit, 
        database=database, table=table,
    )
    

if __name__=='__main__':

    main_mpi()
