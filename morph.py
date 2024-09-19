import warnings
import numpy as np
import pandas as pd
import os,sys
import petrofit as pf
from petrofit.photometry import radial_photometry
import photutils
import skimage
import scipy.ndimage as ndi
import scipy.optimize as opt
from astropy.utils.exceptions import AstropyUserWarning
from astropy.cosmology import Planck15

def get_skirt_image(universe, simulation, snapnum, subfindid, camera, img_path, api_key):
    '''
    Get specified SKIRT image from virgotng.
    '''
    if not os.access(img_path,0):
        os.system(f'mkdir -p {img_path}')
        
    img_name = f'{img_path}/shalo_{snapnum:03}-{subfindid}_{camera}_HSC_GRIZY.fits'

    if not os.access(img_name,0):
    
        base_url = f'http://www.tng-project.org/api'
        img_url = f'{base_url}/{simulation}/snapshots/{snapnum}/subhalos/{subfindid}/skirt/skirt_images_hsc_realistic_{camera}.fits'
        sys_cmd = f'wget -nc -nv --content-disposition --header="API-Key:{api_key}" "{img_url}" -P {img_path}'
        os.system(sys_cmd)

    return img_name

def gaussian_kernel(sigma=3.,dims=(7,7)):
    '''
    Generate gaussian kernel for sep smoothing.
    '''
    # create kernel grid
    hw_x = int(dims[1]/2)
    hw_y = int(dims[0]/2)
    y, x = np.mgrid[-hw_y:hw_y+1, 
                    -hw_x:hw_x+1]
    # generate analytic Gaussian PSF
    kernel = np.exp(-(x**2 + y**2)/(2.0*sigma**2))
    # normalize
    kernel /= np.sum(kernel)
    return kernel

def image_moments(
    image,segmap,pflag,plot=False,
    pixelscale=0.168,vmin=18,vmax=30,
):
    '''
    Compute galaxy moments from galaxy pixels in [image]. 
    The target galaxy's flag in the segmentation image [segmap] should be [pflag]. 
    All other source flags are irrelevant. 
    The sky pixel flag should be 0. 
    The moments can be computed for all flagged sources in an image iteratively for each [pflag]. 
    Returns dictionary of image moments.
    
    References: 
    https://github.com/asgr/ProFound/blob/master/R/profoundSegim.R
    
    '''
    mask = (segmap==pflag)
    rows, cols = np.nonzero(mask)
    fluxes = image[rows,cols]
    # ignore negative fluxes
    fluxes[fluxes<0] = 0.
    
    m00 = np.sum(fluxes)
    x = cols + 0.5
    y = rows + 0.5
    xbar = np.sum(x * fluxes) / m00
    ybar = np.sum(y * fluxes) / m00
    x = x-xbar
    y = y-ybar
    sxx = np.sum(x**2 * fluxes)
    syy = np.sum(y**2 * fluxes)
    sxy = np.sum(x*y * fluxes)
    xvar = sxx / m00
    yvar = syy / m00
    covxy = sxy / m00

    u1 = -xvar-yvar
    u2 = xvar*yvar - covxy**2
    eigval_p = (-u1 + np.sqrt(u1**2-4*u2))/2
    eigval_m = (-u1 - np.sqrt(u1**2-4*u2))/2
    # add variance of uniform RV in quadrature
    semimaj = np.sqrt( eigval_p + 1./12 )
    semimin = np.sqrt( eigval_m + 1./12 )
    axrat = semimin/semimaj
    eigvec = (xvar-eigval_p)/covxy
    ang = 90-np.arctan(eigvec)*180/np.pi
    ang += (ang<0)*180
    
    if plot:
        import matplotlib.pyplot as plt
        _fig,_ax = plt.subplots(figsize=(5,5))
        from matplotlib.patches import Ellipse
        _ax.imshow(
            22.5-2.5*np.log10(image/pixelscale**2), origin='lower',
            interpolation='None',vmin=vmin,vmax=vmax, cmap='bone'
        )
        _ax.scatter(xbar,ybar,s=10)
        e = Ellipse(xy=[xbar,ybar],width=semimin,height=semimaj, 
                    edgecolor='red',facecolor='None',alpha=1., 
                    angle=ang,transform=_ax.transData)
        _ax.add_artist(e)
    return {
        'flux':m00,'xcen':xbar,'ycen':ybar, 'xvar':xvar,
        'yvar':yvar,'covxy':covxy,'semimaj':semimaj,
        'semimin':semimin,'axrat':axrat,'ang':ang,
    }

def make_segmap(data, variance, coords, detect_threshold=0.5, deblend_nthresh=32,
                filter_kernel=gaussian_kernel(5,(11,11)), deblend_cont=0.001,
                clean=True, clean_param=1.0, minarea=5):
    '''
    Make segmentation image, starting hot and getting colder.
    '''
    import sep
    sep.set_extract_pixstack(10000000)
    sep_flag = 0
    
    while True:
        # Make object list and segmentation image.
        obj,segmap = sep.extract(
            data, #.byteswap().newbyteorder(), 
            thresh=detect_threshold, err=np.sqrt(variance), mask=None, 
            minarea=minarea, filter_kernel=filter_kernel, filter_type='conv', 
            deblend_nthresh=deblend_nthresh, deblend_cont=deblend_cont, 
            clean=True, clean_param=1.0, segmentation_map=True)
        pflag = segmap[coords]
        
        # Set flag if detection threshold exceeds limit.
        if detect_threshold>3: sep_flag=1
    
        # Augment threshold if galaxy segmap is too big (touches image edge).
        edges = np.array([segmap[:,0],segmap[0,:],segmap[:,-1],segmap[-1,:]])
        if (pflag in edges.flatten()) or (np.sum(segmap==pflag)>0.8*data.size):
            detect_threshold+=0.5
            continue
        else:
            moments = image_moments(data,segmap,pflag,plot=False)
            break     
    return segmap,moments,pflag,sep_flag

def plot_segim_contours(ax,segim):
    '''Use opencv and segmap to make contour segim like ProFound and add them to axis.'''
    import cv2 as cv
    # Get the unique object labels from the image (excluding sky)
    labels = np.unique(segim)[1:]
    # Iterate over each object label
    for label in labels:
        # Create a mask for the current object
        mask = segim == label
        # Find the contours of the object
        contours, hierarchy = cv.findContours(mask.astype(np.uint8), cv.RETR_EXTERNAL, cv.CHAIN_APPROX_SIMPLE)
        # Draw the contours on the image
        for contour in contours:
            contour = contour.reshape(-1,2)
            contour = np.vstack([contour, contour[0]])
            ax.plot(contour[:,0], contour[:,1], linewidth=1)

def shift_img(img, shift_pix, order=1):
    shift_pix = shift_pix[::-1] 
    shifted_digit_image=shift(img, shift_pix, order = order)
    return shifted_digit_image
    
def rotate_image(img, rotate_pix, order =1):
    shift_pix = [-rotate_pix[0]*2, -rotate_pix[1]*2]
    shift_ = shift_img(img, shift_pix, order=order)
    rotate = np.flip(shift_)
    return rotate

def petrosian_radius(
    data,variance, segmentation, moments, primary_mask, secondary_mask,
    elliptical = True, eta = 0.2, max_apertures = 100):
    '''
    Compute Petrosian radius using PetroFit. If `elliptical`, use the axis ratio and position angle from image moments to define elliptical apertures. `eta` sets the Petrosian ratio.
    '''
    r_max = np.sqrt(data.shape[0]*data.shape[1])/np.sqrt(2)
    n_samples = min(int(2*r_max),max_apertures)
    r_list = pf.make_radius_list(max_pix=r_max, n=n_samples,log=True)

    masked_data = np.ones_like(data)*data
    masked_data[secondary_mask] = np.nan
    masked_stderr = np.ones_like(variance)*np.sqrt(variance)
    masked_stderr[secondary_mask] = np.nan
    skySigma = np.nanstd(data[segmentation==0])

    xc = moments['xcen']
    yc = moments['ycen']
    if not elliptical: 
        phi, q = 0, 1
    else:
        phi = moments['ang']
        q = moments['axrat']
    
    flux_arr, area_arr, error_arr = pf.photometry.radial_photometry(
        masked_data, (xc, yc), r_list, error=masked_stderr, 
        mask=~secondary_mask, elong=1./q, theta=phi*np.pi/180+np.pi/2,
        plot=False, vmin=-3*skySigma, vmax=3*skySigma, method='exact')

    petro = pf.Petrosian(
        r_list, area_arr, flux_arr, flux_err=error_arr)

    return petro.r_petrosian,r_list,flux_arr


def total_flux_fraction(radius,
    data, variance, segmentation, moments, primary_mask, secondary_mask,
    total_fraction, total_flux, elliptical = True):
    '''
    Helper function for calculating radius containing desired fraction
    of total flux (radius_total_flux_fraction). Returns difference between
    desired total flux fraction and total flux fraction in radius.
    '''
    assert (radius >= 0) & (total_fraction >= 0) & (total_fraction <= 1) & (total_flux > 0)
    if radius == 0:
        current_fraction = 0.0
    else:
        if elliptical:
            ap = photutils.aperture.EllipticalAperture(
                (moments['xcen'],moments['ycen']),a=radius, b=radius*moments['axrat'],
                theta=moments['ang']*np.pi/180+np.pi/2)
        else:
            ap = photutils.aperture.CircularAperture(
                (moments['xcen'],moments['ycen']),radius)

        # Force flux sum to be positive:
        ap_flux = np.abs(ap.do_photometry(data * primary_mask, method='exact')[0][0])
        current_fraction = ap_flux / total_flux

    # return value to be optimized
    return current_fraction - total_fraction

def radius_total_flux_fraction(    
    data, variance, segmentation, moments, primary_mask, secondary_mask,
    total_radius, total_fraction, elliptical = True):
    """
    Return the radius (in pixels) of an ellipse or circle that contains a 
    specified fraction of the total flux contained in the aperture defined
    by r_total. If elliptical, r_total is the semi-major axis. 
    """
    if elliptical:
        ap_total = photutils.aperture.EllipticalAperture(
            (moments['xcen'],moments['ycen']),a=total_radius, b=total_radius*moments['axrat'],
            theta=moments['ang']*np.pi/180+np.pi/2)
    else:
        ap_total = photutils.aperture.CircularAperture(
            (moments['xcen'],moments['ycen']),total_radius)

    total_flux = ap_total.do_photometry(data * primary_mask, method='exact')[0][0]

    # Find appropriate range for root finder
    npoints = 100
    r_grid = np.linspace(0.0, total_radius, num=npoints)
    i = 0  # initial value
    while True:
        assert i < npoints, 'Root not found within range.'
        r = r_grid[i]
        curval = total_flux_fraction(r,
            data, variance, segmentation, moments, primary_mask, secondary_mask,
            total_fraction, total_flux, elliptical)
        if curval <= 0:
            r_min = r
        elif curval > 0:
            r_max = r
            break
        i += 1
        
    r = opt.brentq(total_flux_fraction, r_min, r_max,
                   args=(data, variance, segmentation, moments, primary_mask, 
                         secondary_mask, total_fraction, total_flux, elliptical), 
                   xtol=1e-6)
    return r



# def sersic_ftot(sersic_Ie,sersic_Re,sersic_n,sersic_q):
#     '''
#     Compute total Sérsic model intensity given model params. Integrated to infinity. See Galfit 2010 paper: Peng+ 2010AJ....139.2097P (Equation 4).
#     '''
#     from scipy.special import gamma
#     # bn approx. (Graham & Driver 2005)
#     if sersic_n <= 0.36:
#         # MacArthur, Courteau, & Holtzman (2003)
#         bn = 0.01945 - 0.8902*sersic_n + 10.95*sersic_n**2 - 19.67*sersic_n**3 + 13.43*sersic_n**4
#     else:
#         # Ciotti & Bertin (1999)
#         bn = 2*sersic_n - 1./3 + 4/405/sersic_n + 46/25515/sersic_n**2 + 131/1148175/sersic_n**3 - 2194697/30690717750/sersic_n**4
#     ftot = sersic_Ie*2*np.pi*sersic_Re**2*sersic_n*np.exp(bn)/bn**(2*sersic_n)
#     ftot *= gamma(2*sersic_n)*sersic_q
#     return ftot

# def df_nullint(df):
#     '''Convert int to nullible int types.'''
#     dtypes = dict(df.dtypes)
#     for key in dtypes.keys():
#         if 'int' in str(dtypes[key]): 
#             dtypes[key] = pd.Int64Dtype()
#     return df.astype(dtypes) 

# def query_df(dbcmd,database,
#              host='nantai.ipmu.jp',user='bottrell',
#              cnf_path='/home/connor.bottrell/.mysql/nantai.cnf'):
#     '''Query database and return results as dataframe.'''
#     import pymysql
#     db = pymysql.connect(host=host,
#                          user=user,
#                          database=database,
#                          read_default_file=cnf_path)
#     df = pd.read_sql(dbcmd, con=db)
#     db.close()
#     df = df_nullint(df)
#     return df

# def query(dbcmd,database,
#           host='nantai.ipmu.jp',user='bottrell',
#           cnf_path='/home/connor.bottrell/.mysql/nantai.cnf'):
#     '''Query database and return results as numpy array.'''
#     import pymysql
#     import numpy as np
#     db = pymysql.connect(host=host,
#                          user=user,
#                          database=database,
#                          read_default_file=cnf_path)
#     c = db.cursor()
#     c.execute(dbcmd)
#     data = np.asarray(c.fetchall())
#     c.close()
#     db.close()
#     return data

# def submit(dbcmd,database,
#            host='nantai.ipmu.jp',user='bottrell',
#            cnf_path='/home/connor.bottrell/.mysql/nantai.cnf'):
#     '''Submit command to database.'''
#     import pymysql
#     import numpy as np
#     db = pymysql.connect(host=host,
#                          user=user,
#                          database=database,
#                          read_default_file=cnf_path,
#                          autocommit=True
#                         )
#     c = db.cursor()
#     c.execute(dbcmd)
#     c.close()
#     db.close()
#     return 

# def cbplot(nrows=1,ncols=1,sharex=False,sharey=False,
#            panelsize_x=5,panelsize_y=5):
#     '''My plotting style.'''
#     import matplotlib.pyplot as plt
#     fig,axarr = plt.subplots(nrows,ncols,figsize=(ncols*panelsize_x, nrows*panelsize_y), sharey=sharey,sharex=sharex,facecolor='white')
#     for ax in np.asarray([axarr]).flat:
#         ax.minorticks_on()
#         ax.tick_params(axis='both',which='major',direction='in',
#                        length=10,width=1,labelsize=18,right=1,top=1)
#         ax.tick_params(axis='both',which='minor',direction='in',
#                        length=5,width=0.5,right=1,top=1)
#         ax.tick_params(axis='both',which='major',pad=10)
#         for axis in ['top','bottom','left','right']:
#             ax.spines[axis].set_linewidth(1)
#     return fig,axarr

class nonparametric:
    '''
    Compute non-parametric morphologies
    '''
    def __init__(self, 
        data, variance, segmentation, model, residual, primary_mask, secondary_mask, moments, psf,
        pixel_scale=1, zeropoint = 22.5, petro_eta = 0.2, petro_extent_cas = 1.5, petro_extent_flux = 2.0,
        petro_fraction_cas = 0.25, petro_fraction_gini = 0.2, annulus_width = 1.0, 
        redshift = -99., cosmology = Planck15
    ):
        self.data = data
        self.variance = variance
        self.segmentation = segmentation
        self.model = model
        self.residual = residual
        self.primary_mask = primary_mask
        self.secondary_mask = secondary_mask
        self.moments = moments
        self.psf = psf
        self.pixel_scale = pixel_scale
        self.zeropoint = zeropoint
        self.petro_eta = petro_eta
        self.petro_extent_cas = petro_extent_cas
        self.petro_extent_flux = petro_extent_flux
        self.petro_fraction_cas = petro_fraction_cas
        self.petro_fraction_gini = petro_fraction_gini
        self.annulus_width = 1.0
        self.redshift = redshift
        self.cosmology = cosmology # Astropy cosmology object
        
        self.flag = 0
        
        self.moments_center = (self.moments['xcen'], self.moments['ycen'])
        
        # Some statistics for source and background
        
        self.sky_mean = np.nanmean(self.data[self.segmentation==0])
        self.sky_median = np.nanmedian(self.data[self.segmentation==0])
        self.sky_sigma = np.nanstd(self.data[self.segmentation==0])
        self.sky_npix = int(np.nansum(self.segmentation==0))
        self.target_npix = int(np.nansum(self.primary_mask))
        
        # Petrosian radii (elliptical / circular)
        
        self.rp_ellipse, r_list, flux_arr = self._petrosian_radius(
            elliptical=True)

        self.rp_circle, r_list, flux_arr = self._petrosian_radius(
            elliptical=False)
        
        # Half light radii (non-parametric)
        
        self.r50_ellipse = self._radius_total_flux_fraction( 
            self.petro_extent_cas*self.rp_circle,
            0.5, elliptical=True)
        
        self.r50_circle = self._radius_total_flux_fraction( 
            self.petro_extent_cas*self.rp_circle,
            0.5, elliptical=False)
        
        # Background asymmetry density using all sky pixels
        
        self.bkg_asymmetry_density = self._bkg_asymmetry_density(
            self.moments_center)
        
        # CAS Asymmetry (inside 1.5 x Petrosian radius)

        self.asymmetry_center = self._asymmetry_center(
            self.data, kind='cas')
        
        self.asymmetry_cas = self._asymmetry_function(
            self.asymmetry_center, self.data, kind='cas')
        
        # Outer Asymmetry
        
        self.rmax_circle = self._rmax_circ(
            self.asymmetry_center)
        
        self.rmax_ellipse = self._rmax_ellipse(
            self.asymmetry_center)
        
        self.asymmetry_outer = self._asymmetry_function(
            self.asymmetry_center, self.data, kind='outer')
        
        # Shape Asymmetry (use primary mask as data)
        
        self.asymmetry_shape = self._asymmetry_function(
            self.asymmetry_center, self.primary_mask, kind='shape')
        
        # CAS Asymmetry+ (all primary source pixels, no aperture)
        
        self.asymmetry_center_nap = self._asymmetry_center_no_aperture(
            self.data, kind='cas')
        
        self.asymmetry_cas_nap = self._asymmetry_function_no_aperture(
            self.asymmetry_center_nap, self.data, kind='cas')
        
        # Residual Asymmetry
        
        self.asymmetry_residual_center = self._asymmetry_center(
            self.data-self.model, kind='residual')
        
        self.asymmetry_residual = self._asymmetry_function(
            self.asymmetry_residual_center, self.data-self.model, kind='residual')
        
        # Residual Asymmetry+ (all primary source pixels, no aperture)
        
        self.asymmetry_residual_center_nap = self._asymmetry_center_no_aperture(
            self.data-self.model, kind='residual')
        
        self.asymmetry_residual_nap = self._asymmetry_function_no_aperture(
            self.asymmetry_residual_center_nap, self.data-self.model, kind='residual')
        
        # Square!!! of RMS Asymmetry (A_RMS^2)
        
        self.asymmetry_rms2 = self._asymmetry_function(
            self.asymmetry_center, self.data, kind='rms')
        
        # CAS Concentration
        
        self.r20_ellipse = self._radius_total_flux_fraction( 
            self.petro_extent_cas*self.rp_circle,
            0.20, elliptical=True)
        
        self.r20_circle = self._radius_total_flux_fraction( 
            self.petro_extent_cas*self.rp_circle,
            0.20, elliptical=False)
        
        self.r80_ellipse = self._radius_total_flux_fraction( 
            self.petro_extent_cas*self.rp_circle,
            0.80, elliptical=True)
        
        self.r80_circle = self._radius_total_flux_fraction( 
            self.petro_extent_cas*self.rp_circle,
            0.80, elliptical=False)
        
        self.concentration_ellipse = 5 * np.log10(self.r80_ellipse / self.r20_ellipse)
        
        self.concentration_circle = 5 * np.log10(self.r80_circle / self.r20_circle)
        
        # CAS Smoothness/Clumpiness
        
        self.bkg_smoothness_density = self._bkg_smoothness_density()
        
        self.smoothness = self._smoothness()
        
        # Gini Coefficient
        
        self.gini = self._gini()
        
        # M20
        
        self.m20 = self._m20()
        
        # Average surface brightness in 1 kpc
        
        self.sb_1kpc = self._sb_1kpc()

           
    def _asymmetry_function(
        self, center, data, kind = 'cas'):
        """
        SHAMELESSLY copied from statmorph.
        Helper function to determine the asymmetry and center of asymmetry.
        The idea is to minimize the output of this function.

        Parameters
        ----------
        center : tuple or array-like
            The (x,y) position of the center.
        data : array-like
            The 2D image.
        kind : {'cas', 'rms', 'outer', 'shape', 'residual'}
            Whether to calculate the traditional CAS asymmetry (default),
            RMS asymmetry, outer asymmetry, shape asymmetry.

        Returns
        -------
        asym : The asymmetry statistic for the given center.

        """
        
        data = np.float64(data)
        ny, nx = data.shape
        xc, yc = center

        if xc < 0 or xc >= nx or yc < 0 or yc >= ny:
            warnings.warn('[asym_center] Minimizer tried to exit bounds.',
                          AstropyUserWarning)
            asymmetry_flag = 1
            # Return large value to keep minimizer within range:
            return 100.0

        # Rotate around given center
        data_180 = skimage.transform.rotate(data, 180.0, center=center)

        # Apply symmetric mask
        mask = self.secondary_mask.copy() 
        mask_180 = skimage.transform.rotate(mask, 180.0, center=center)
        mask_180 = mask_180 >= 0.5  # convert back to bool
        mask_symmetric = mask | mask_180
        data = np.where(~mask_symmetric, data, 0.0)
        data_180 = np.where(~mask_symmetric, data_180, 0.0)

        # Create aperture for the chosen kind of asymmetry
        if kind == 'cas' or kind == 'rms' or kind == 'residual':
            r = self.petro_extent_cas * self.rp_circle
            ap = photutils.aperture.CircularAperture(center, r)      
        elif kind == 'outer':
            a_in = self.r50_ellipse
            a_out = self.rmax_ellipse
            b_out = a_out * self.moments['axrat']
            theta = self.moments['ang']*np.pi/180+np.pi/2
            assert (a_in > 0) & (a_out > 0)
            ap = photutils.aperture.EllipticalAnnulus(center, a_in, a_out, b_out, theta=theta)
        elif kind == 'shape':
            if np.isnan(self.rmax_circle) or (self.rmax_circle <= 0):
                warnings.warn('[shape_asym] Invalid rmax_circle value.',
                              AstropyUserWarning)
                return -99.0  # invalid
            ap = photutils.aperture.CircularAperture(center, self.rmax_circle)
        else:
            raise NotImplementedError('Asymmetry kind not understood:', kind)

        # Aperture area (in pixels)
        ap_area = ap.do_photometry(~mask_symmetric, method='exact')[0][0]
        
        if kind == 'residual':
            ap_abs_sum = ap.do_photometry(
                np.abs(self.model * ~mask_symmetric), method='exact')[0][0]
            ap_abs_diff = ap.do_photometry(np.abs(data_180-data), method='exact')[0][0]
        else:
            ap_abs_sum = ap.do_photometry(np.abs(data), method='exact')[0][0]
            ap_abs_diff = ap.do_photometry(np.abs(data_180-data), method='exact')[0][0]
    
        if ap_abs_sum == 0.0:
            warnings.warn('[asymmetry_function] Zero flux sum.',
                          AstropyUserWarning)
            # Return large value to get minimizer out of masked region:
            return 100.0

        if kind == 'shape':
            # The shape asymmetry of the background is zero
            asym = ap_abs_diff / ap_abs_sum
        elif kind == 'rms':
            # Apply eq. 27 from Sazonova et al. (2024)
            ap_sqr_sum = ap.do_photometry(data**2, method='exact')[0][0]
            ap_sqr_diff = ap.do_photometry((data_180-data)**2, method='exact')[0][0]

            asym = (ap_sqr_diff - 2*ap_area*self.sky_sigma**2) / (
                ap_sqr_sum - ap_area*self.sky_sigma**2)
        else:
            asym = (ap_abs_diff - ap_area*self.bkg_asymmetry_density) / ap_abs_sum

        return asym
    
    def _asymmetry_center(self, data, kind='cas'):
        """
        Find the position of the central pixel (relative to the
        "postage stamp" cutout) that minimizes the (CAS) asymmetry.
        """
        center = np.array([self.moments['xcen'], self.moments['ycen']])  # initial guess
        
        center_asym = opt.fmin(self._asymmetry_function, center,
                               args=(data, kind),
                               xtol=1e-6, disp=False)

        # If the asymmetry center drifted too far away from moment center
        # use moments center as asymmetry center
        dist = self.petro_extent_cas * self.rp_circle
        if np.linalg.norm(center_asym - center) > dist:
            center_asym = center
        return center_asym
    
    def _asymmetry_function_no_aperture(
        self, center, data, kind = 'cas'):
        """
        Same as previous asymmetry function but using all primary 
        source pixels and no aperture photometry. Only CAS and 
        residual asymmetries are supported.
        """
        
        data = np.float64(data)
        ny, nx = data.shape
        xc, yc = center

        if xc < 0 or xc >= nx or yc < 0 or yc >= ny:
            warnings.warn('[asym_center] Minimizer tried to exit bounds.',
                          AstropyUserWarning)
            # Return large value to keep minimizer within range:
            return 100.0

        # Rotate around given center
        data_180 = skimage.transform.rotate(data, 180.0, center=center)

        # Make secondary mask symmetric
        secondary = self.secondary_mask.copy()
        secondary_180 = skimage.transform.rotate(secondary, 180.0, center=center)
        secondary_180 = secondary_180 >= 0.5  # convert back to bool
        secondary_symmetric = secondary | secondary_180
        
        # Make primary mask symmetric
        primary = self.primary_mask.copy()
        primary_180 = skimage.transform.rotate(primary, 180.0, center=center)
        primary_180 = primary_180 >= 0.5  # convert back to bool
        primary_symmetric = primary | primary_180
        
        # Use pixels which are part of the primary and not secondary
        mask_symmetric = ~(primary_symmetric * ~secondary_symmetric)
        
        data = np.where(~mask_symmetric, data, 0.0)
        data_180 = np.where(~mask_symmetric, data_180, 0.0)
        
        # Primary target footprint area (in pixels)
        ap_area = np.sum(~mask_symmetric)

        if kind == 'residual':
            ap_abs_sum = np.sum(np.abs(self.model * ~mask_symmetric))
            ap_abs_diff = np.sum(np.abs(data_180-data))
        else:
            ap_abs_sum = np.sum(np.abs(data))
            ap_abs_diff = np.sum(np.abs(data_180-data))
    
        if ap_abs_sum == 0.0:
            warnings.warn('[asymmetry_function] Zero flux sum.',
                          AstropyUserWarning)
            # Return large value to get minimizer out of masked region:
            return 100.0

        asym = (ap_abs_diff - ap_area*self.bkg_asymmetry_density) / ap_abs_sum
                
        return asym
    
    def _asymmetry_center_no_aperture(self, data, kind='cas'):
        """
        Asymmetry center computed using all primary source pixels
        without aperture photometry. 
        """
        center = np.array([self.moments['xcen'], self.moments['ycen']])  # initial guess
        
        center_asym = opt.fmin(self._asymmetry_function_no_aperture, center,
                               args=(data, kind),
                               xtol=1e-6, disp=False)

        # If the asymmetry center drifted too far away from moment center
        # use moments center as asymmetry center
        dist = self.petro_extent_cas * self.rp_circle
        if np.linalg.norm(center_asym - center) > dist:
            center_asym = center
        return center_asym
    
    def _bkg_asymmetry_density(self, center):
        '''
        Estimate asymmetry background density from pixels
        that are not flagged as belonging to any source.
        '''
        data = np.float64(self.data)
        
        # Rotate around given center
        data_180 = skimage.transform.rotate(data, 180.0, center=center)
    
        # mask all pixels with source flags
        mask = self.primary_mask | self.secondary_mask
        mask_180 = skimage.transform.rotate(mask, 180.0, center=center)
        mask_180 = mask_180 >= 0.5  # convert back to bool
        mask_symmetric = mask | mask_180
        
        data = np.where(~mask_symmetric, self.data, 0.0)
        data_180 = np.where(~mask_symmetric, data_180, 0.0)
        
        # use full background to background estimate density
        bkg_asy_diff = np.sum(np.abs(data_180-data))
        bkg_asy_npix = np.sum(~mask_symmetric)
        
        return bkg_asy_diff/bkg_asy_npix
    
    def _smoothness(self):
        """
        Calculate smoothness (a.k.a. clumpiness) as defined in eq. (11)
        from Lotz et al. (2004). Note that the original definition by
        Conselice (2003) includes an additional factor of 10.
        """
        image = self.data * ~self.secondary_mask

        # Exclude central region during smoothness calculation:
        r_in = self.petro_fraction_cas * self.rp_circle
        r_out = self.petro_extent_cas * self.rp_circle
        ap = photutils.aperture.CircularAnnulus(self.asymmetry_center, r_in, r_out)

        boxcar_size = int(self.petro_fraction_cas * self.rp_circle)
        image_smooth = ndi.uniform_filter(image, size=boxcar_size)

        image_diff = image - image_smooth
        image_diff[image_diff < 0] = 0.0  # set negative pixels to zero

        ap_flux = ap.do_photometry(image, method='exact')[0][0]
        ap_diff = ap.do_photometry(image_diff, method='exact')[0][0]
        ap_area = ap.do_photometry(~self.secondary_mask, method='exact')[0][0]

        if ap_flux <= 0:
            warnings.warn('[smoothness] Nonpositive total flux.',
                          AstropyUserWarning)
            return -99.0  # invalid

        S = (ap_diff - ap.area*self.bkg_smoothness_density) / ap_flux

        if not np.isfinite(S):
            warnings.warn('Invalid smoothness.', AstropyUserWarning)
            return -99.0  # invalid

        return S
    
    def _bkg_smoothness_density(self):
        """
        Smoothness of the background. 
        Note the peculiar normalization (for reference only).
        """
        bkg = self.data.copy() 
        source_mask = self.primary_mask | self.secondary_mask
        bkg[source_mask] = 0.
        
        if bkg[~source_mask].size == 0:
            return 0.

        # If the smoothing "boxcar" is larger than the skybox itself,
        # this just sets all values equal to the mean:
        boxcar_size = int(self.petro_fraction_cas * self.rp_circle)
        bkg_smooth = ndi.uniform_filter(bkg, size=boxcar_size)

        bkg_diff = bkg - bkg_smooth
        bkg_diff[bkg_diff < 0] = 0. # set negative pixels to zero

        return np.sum(bkg_diff[~source_mask]) / float(bkg[~source_mask].size)
    
    def _rmax_circ(self, center):
        """
        Return distance from centre to the most distant pixel 
        on the primary target segment.
        """
        ny, nx = self.data.shape
        xc, yc = center

        # Distances from all pixels to the center
        ypos, xpos = np.mgrid[0:ny, 0:nx]
        distances = np.sqrt((ypos-yc)**2 + (xpos-xc)**2)

        # Only consider pixels within the segmap.
        return np.max(distances[self.primary_mask])
    
    def _rmax_ellipse(self, center):
        """
        Return the semimajor axis of the minimal ellipse (with fixed
        center, elongation and orientation) that contains the entire
        primary target segment.
        """
        ny, nx = self.data.shape
        xc, yc = center

        # Distances from all pixels to the center
        ypos, xpos = np.mgrid[0:ny, 0:nx]

        theta = self.moments['ang']*np.pi/180+np.pi/2
        y, x = np.mgrid[0:ny, 0:nx]

        xprime = (x-xc)*np.cos(theta) + (y-yc)*np.sin(theta)
        yprime = -(x-xc)*np.sin(theta) + (y-yc)*np.cos(theta)
        r_ellipse = np.sqrt(xprime**2 + (yprime/self.moments['axrat'])**2)

        # Only consider pixels within the segmap.
        return np.max(r_ellipse[self.primary_mask])
    
    def _petrosian_radius(self, elliptical = True, max_apertures = 100):
        '''
        Compute Petrosian radius using PetroFit. If `elliptical`, use the axis ratio and position angle from image moments to define elliptical apertures. `eta` sets the Petrosian ratio.
        '''
        r_max = np.sqrt(self.data.size)/np.sqrt(2)
        n_samples = min(int(2*r_max),max_apertures)
        r_list = pf.make_radius_list(max_pix=r_max, n=n_samples,log=True)

        masked_data = np.ones_like(self.data)*self.data
        masked_data[self.secondary_mask] = np.nan
        masked_stderr = np.ones_like(self.variance)*np.sqrt(self.variance)
        masked_stderr[self.secondary_mask] = np.nan
        skySigma = np.nanstd(self.data[self.segmentation==0])

        xc = self.moments['xcen']
        yc = self.moments['ycen']
        if not elliptical: 
            phi, q = 0, 1
        else:
            phi = self.moments['ang']
            q = self.moments['axrat']

        flux_arr, area_arr, error_arr = pf.photometry.radial_photometry(
            masked_data, (xc, yc), r_list, error=masked_stderr, 
            mask=~self.secondary_mask, elong=1./q, theta=phi*np.pi/180+np.pi/2,
            plot=False, vmin=-3*skySigma, vmax=3*skySigma, method='exact')

        petro = pf.Petrosian(
            r_list, area_arr, flux_arr, flux_err=error_arr)

        if np.isnan(petro.r_petrosian):
            return -99.,r_list,flux_arr
        else:
            return petro.r_petrosian,r_list,flux_arr

    def _total_flux_fraction(
        self, radius, total_fraction, total_flux, elliptical=True):
        '''
        Helper function for calculating radius containing desired fraction
        of total flux (radius_total_flux_fraction). Returns difference between
        desired total flux fraction and total flux fraction in radius.
        '''
        assert (radius >= 0) & (total_fraction >= 0) & (total_fraction <= 1) & (total_flux > 0)
        if radius == 0:
            current_fraction = 0.0
        else:
            if elliptical:
                ap = photutils.aperture.EllipticalAperture(
                    (self.moments['xcen'],self.moments['ycen']),a=radius, b=radius*self.moments['axrat'],
                    theta=self.moments['ang']*np.pi/180+np.pi/2)
            else:
                ap = photutils.aperture.CircularAperture(
                    (self.moments['xcen'],self.moments['ycen']),radius)

            # Force flux sum to be positive:
            ap_flux = np.abs(ap.do_photometry(self.data * self.primary_mask, method='exact')[0][0])
            current_fraction = ap_flux / total_flux

        # return value to be optimized
        return current_fraction - total_fraction

    def _radius_total_flux_fraction(
        self, total_radius, total_fraction, elliptical=True):
        """
        Return the radius (in pixels) of an ellipse or circle that contains a 
        specified fraction of the total flux contained in the aperture defined
        by r_total. If elliptical, r_total is the semi-major axis. 
        """
        if elliptical:
            ap_total = photutils.aperture.EllipticalAperture(
                (self.moments['xcen'],self.moments['ycen']),
                a=total_radius, b=total_radius*self.moments['axrat'],
                theta=self.moments['ang']*np.pi/180+np.pi/2)
        else:
            ap_total = photutils.aperture.CircularAperture(
                (self.moments['xcen'],self.moments['ycen']),total_radius)

        total_flux = ap_total.do_photometry(self.data * self.primary_mask, method='exact')[0][0]

        # Find appropriate range for root finder
        npoints = 100
        r_grid = np.linspace(0.0, total_radius, num=npoints)
        i = 0  # initial value
        while True:
            assert i < npoints, 'Root not found within range.'
            r = r_grid[i]
            curval = self._total_flux_fraction(
                r,total_fraction, total_flux, elliptical)
            if curval <= 0:
                r_min = r
            elif curval > 0:
                r_max = r
                break
            i += 1

        r = opt.brentq(self._total_flux_fraction, r_min, r_max,
                       args=(total_fraction, total_flux, elliptical), 
                       xtol=1e-6)
        return r
    
    def _gini(self):
        """
        Calculate the Gini coefficient as described in Lotz et al. (2004).
        Use user-provided segmentation map. 
        """
        image = self.data.flatten()
        segmap = self.primary_mask.flatten()

        sorted_pixelvals = np.sort(np.abs(image[segmap]))
        n = len(sorted_pixelvals)
        if n <= 1 or np.sum(sorted_pixelvals) == 0:
            warnings.warn('[gini] Not enough data for Gini calculation.',
                          AstropyUserWarning)
            return -99.0  # invalid

        indices = np.arange(1, n+1)  # start at i=1
        gini = (np.sum((2*indices-n-1) * sorted_pixelvals) /
                (float(n-1) * np.sum(sorted_pixelvals)))
        
        return gini
    
    def _m20(self):
        """
        Calculate the M_20 coefficient as described in Lotz et al. (2004).
        Use user-provided segmentation map.
        """

        # Use the same region as in the Gini calculation
        image = np.where(self.primary_mask, self.data, 0.0)
        image = np.float64(image)  # skimage wants double
        
        M00 = self.moments['flux']
        xc, yc = self.moments['xcen'], self.moments['ycen']
        Mxx = self.moments['xvar']
        Myy = self.moments['yvar']
        second_moment_tot = Mxx + Myy

        # Calculate threshold pixel value
        sorted_pixelvals = np.sort(image.flatten())
        flux_fraction = np.cumsum(sorted_pixelvals) / np.sum(sorted_pixelvals)
        sorted_pixelvals_20 = sorted_pixelvals[flux_fraction >= 0.8]
        if len(sorted_pixelvals_20) == 0:
            # This can happen when there are very few pixels.
            warnings.warn('[m20] Not enough data for M20 calculation.',
                          AstropyUserWarning)
            return -99.0  # invalid
        threshold = sorted_pixelvals_20[0]

        # Calculate second moment of the brightest pixels
        image_20 = np.where(image >= threshold, image, 0.0)
        # Calculate 2nd order moments of 20% brightest pixels using 1st order origin
        Mc_20 = skimage.measure.moments_central(image_20, center=(yc, xc), order=2)
        second_moment_20 = Mc_20[0, 2] + Mc_20[2, 0]

        if (second_moment_20 <= 0) | (second_moment_tot <= 0):
            warnings.warn('[m20] Negative second moment(s).',
                          AstropyUserWarning)
            m20 = -99.0  # invalid
        else:
            m20 = np.log10(second_moment_20 / second_moment_tot)

        return m20
    
    def _sb_1kpc(self):
        '''
        Average surface brightness (mag/arcsec2) within 1 kpc
        of moment centroid using a circular aperture.
        '''
        
        # proper units scaling (kpc/arcsecond)
        arcsec_per_kpc = self.cosmology.arcsec_per_kpc_proper(
            self.redshift).value
        
        # radius of circular aperture (pixels)
        radius = 0.5 * arcsec_per_kpc / self.pixel_scale 
        
        ap = photutils.aperture.CircularAperture(
            (self.moments['xcen'],self.moments['ycen']),radius)
        
        ap_sum = ap.do_photometry(self.data * self.primary_mask, method='exact')[0][0]
        ap_area = ap.do_photometry(self.primary_mask, method='exact')[0][0]
        
        SB_1kpc = self.zeropoint - 2.5 * np.log10( ap_sum / ap_area / self.pixel_scale**2 ) 
        
        return SB_1kpc
        

# npmorph = nonparametric(
#     data = target_data, variance = target_variance, 
#     segmentation = segmap, model = model_data, 
#     residual = residual_data, primary_mask = pri_mask,
#     secondary_mask = sec_mask, moments = moments, psf = target_psf,
#     pixel_scale = pixel_scale, redshift = redshift, cosmology = Planck15
# )
