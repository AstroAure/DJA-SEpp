from glob import glob
import numpy as np
from sourcextractor.config import *

import os


##### THESE NEXT LINES DO SOME MORE CONFIGURATION #######
fit_case = 'sersic' # this defines which model we want to fit
set_engine('levmar')
set_max_iterations(250)
use_iterative_fitting(True)
set_meta_iterations(3)
set_deblend_factor(0.95) # can be 0.95
set_meta_iteration_stop(0.01)  ## increase -> faster!

###################################################################################
################################ LOAD IMAGES ######################################

field = 'GDS'
list_of_IMG_names = glob(f"/home/aurelien/DAWN/DJA_SE++/image/{field}/cutout/*clear*sci*.fits")
list_of_WHT_names = glob(f"/home/aurelien/DAWN/DJA_SE++/image/{field}/cutout/*clear*wht*.fits")
list_of_PSF_names = glob(f"/home/aurelien/DAWN/DJA_SE++/psfex/{field}/*star_psf.psf")

# mag_zeropoint = {'F090W': 28.02,
#                  'F115W': 28.02,
#                  'F150W': 28.02,
#                  'F182M': 28.02,
#                  'F200W': 28.02,
#                  'F210M': 28.02,
#                  'F277W': 26.49,
#                  'F335M': 26.49,
#                  'F356W': 26.49,
#                  'F410M': 26.49,
#                  'F444W': 26.49}

mag_zeropoint = {'F090W': 28.90,
                 'F115W': 28.90,
                 'F150W': 28.90,
                 'F182M': 28.90,
                 'F200W': 28.90,
                 'F210M': 28.90,
                 'F277W': 28.90,
                 'F335M': 28.90,
                 'F356W': 28.90,
                 'F410M': 28.90,
                 'F444W': 28.90}

list_of_IMG_names = list( map( lambda x: x, list_of_IMG_names) )
list_of_WHT_names = list( map( lambda x: x, list_of_WHT_names) )
list_of_PSF_names = list( map( lambda x: x, list_of_PSF_names) )

# print to verify how they are loaded
imgroup = load_fits_images(
            images      = list_of_IMG_names,
            psfs        = list_of_PSF_names,
            weights     = list_of_WHT_names,
            weight_type = 'rms')


imgroup.split(ByKeyword('FILTER'))
mesgroup = MeasurementGroup(imgroup)
###################################################################################
###################################################################################



###################################################################################
########################### RUN APERTURE PHOTOMETRY ###############################
all_apertures = []
pix_diameter = 25 # pixels

# loop over every band in the measurement image group and measure aperture photometry in each
for band,img in mesgroup:
    all_apertures.extend(add_aperture_photometry(img, pix_diameter) )
add_output_column('APER', all_apertures)
###################################################################################
###################################################################################





###################################################################################
################################ DEFINE MODELS ####################################
if fit_case == "sersic" :
    # first get the source positions from the detection image
    x, y = get_pos_parameters()
    
    ## define Sersic model parameters: radius, angle, axis ratio, sersic index:
    
    rad = FreeParameter(lambda o: o.radius, Range(lambda v, o: (.01 * v, 50 * v), RangeType.EXPONENTIAL))
    angle = FreeParameter( lambda o: o.angle )
    
    ratio = FreeParameter( 0.5, Range( (0.01, 1.0), RangeType.LINEAR) )
    X_rd = DependentParameter( lambda q: np.log( (q-0.01)/(1.01-q) ), ratio )
    add_prior( X_rd, 0.03, 1.0 ) 
    
    sersic = FreeParameter(2.0, Range((0.3, 16.6), RangeType.LINEAR))
    X_sersic = DependentParameter(lambda nt: 1.4**(nt - 7.), sersic)
    add_prior(X_sersic, 0.0, 1.0)

    # we can also get the WCS coordinates, and the radius, angle and ratio in angular units
    # see https://sourcextractorplusplus.readthedocs.io/en/latest/config_api/model_fitting.html?highlight=get_world_parameters#config.model_fitting.get_world_parameters
    ra, dec, wc_rad, wc_angle, wc_ratio = get_world_parameters(x, y, rad, angle, ratio)
    
    # next, write these outputs 
    add_output_column('RA_MODEL', ra)
    add_output_column('DEC_MODEL', dec)
    add_output_column('RADIUS', wc_rad)
    add_output_column('AXRATIO', wc_ratio)
    add_output_column('ANGLE', wc_angle)
    add_output_column('SERSIC', sersic)
    add_output_column('X_SERSIC', X_sersic)
    
    ## The parameters above are fitted on all the images in the measurement group (all bands) in a S/N weighted average way (loosely speaking)
    ## But the flux is fitted on each band separately (but we can actually 
    
    i = 0
    flux = {}
    mag = {}
    # loop over each band in the measurement group
    for band,group in mesgroup: 
        print(band)
        
        flux[i] = get_flux_parameter()
        mag[i] = DependentParameter(lambda f, zp=mag_zeropoint[band]: -2.5 * np.log10(f) + zp, flux[i] )
        
        # add all compontents to the final model that is going to be fitted
        add_model(group, SersicModel(x, y, 
            flux[i], rad, ratio, angle, sersic))
        
        # save output columns for the parameters fitted in each band, in this case the flux and magnitude
        add_output_column(f'MAG_MODEL_{band}', mag[i])
        add_output_column(f'FLUX_MODEL_{band}', flux[i])
        i+=1


# we can print some info about the models        
print_model_fitting_info(mesgroup, show_params=True, prefix='')
