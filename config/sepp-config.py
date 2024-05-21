import sys
from glob import glob
import numpy as np
from sourcextractor.config import *

#### READ ARGUMENTS PASSED THROUGH --python-arg #########
def str2dict(str_dict):
    keys = []
    for e in str_dict.split(":")[:-1]:
        key = e.split(",")[-1]
        key = key.replace("{", "")
        key = key.replace(" ","")
        keys.append(key)
    values = []
    for e in str_dict.split(":")[1:]:
        val = []
        el = e.split(",")
        el = el if el[-1][-1]=="}" else el[:-1]
        for l in el:
            l = l.replace(" ", "")
            l = l.replace("[", "")
            l = l.replace("]", "")
            l = l.replace('}', "")
            val.append(l)
        if len(val)==1: val = val[0]
        values.append(val)
    dic = dict(zip(keys, values))
    return dic
args = str2dict(sys.argv[1])

##### THESE NEXT LINES DO SOME MORE CONFIGURATION #######
fit_case = args['fit_case']

set_engine('levmar')
set_max_iterations(500)
use_iterative_fitting(True)
set_meta_iterations(3)
set_deblend_factor(1.0) # can be 0.95
set_meta_iteration_stop(0.02)  ## increase -> faster!

###################################################################################
################################ LOAD IMAGES ######################################

list_of_IMG_names = args['list_of_IMG_names']
list_of_WHT_names = args['list_of_WHT_names']
list_of_PSF_names = args['list_of_PSF_names']

mag_zeropoint = {'F090W': 28.90,
                 'F115W': 28.90,
                 'F150W': 28.90,
                 'F182M': 28.90,
                 'F200W': 28.90,
                 'F210M': 28.90,
                 'F250M': 28.90,
                 'F277W': 28.90,
                 'F300M': 28.90,
                 'F335M': 28.90,
                 'F356W': 28.90,
                 'F410M': 28.90,
                 'F430M': 28.90,
                 'F444W': 28.90,
                 'F460M': 28.90,
                 'F480M': 28.90,}

filters_waveband = {'F090W': {'pivot': 0.901, 'band': 0.194},
                    'F115W': {'pivot': 1.154, 'band': 0.225},
                    'F150W': {'pivot': 1.501, 'band': 0.318},
                    'F182M': {'pivot': 1.845, 'band': 0.238},
                    'F200W': {'pivot': 1.990, 'band': 0.461},
                    'F210M': {'pivot': 2.093, 'band': 0.205},
                    'F250M': {'pivot': 2.503, 'band': 0.181},
                    'F277W': {'pivot': 2.786, 'band': 0.672},
                    'F300M': {'pivot': 2.996, 'band': 0.318},
                    'F335M': {'pivot': 3.365, 'band': 0.347},
                    'F356W': {'pivot': 3.563, 'band': 0.787},
                    'F410M': {'pivot': 4.092, 'band': 0.436},
                    'F430M': {'pivot': 4.280, 'band': 0.228},
                    'F444W': {'pivot': 4.421, 'band': 1.024},
                    'F460M': {'pivot': 4.624, 'band': 0.228},
                    'F480M': {'pivot': 4.834, 'band': 0.303},}

list_of_IMG_names = list( map( lambda x: x, list_of_IMG_names) )
list_of_WHT_names = list( map( lambda x: x, list_of_WHT_names) )
list_of_PSF_names = list( map( lambda x: x, list_of_PSF_names) )

# print to verify how they are loaded
imgroup = load_fits_images(
            images      = list_of_IMG_names,
            psfs        = list_of_PSF_names,
            weights     = list_of_WHT_names,
            weight_type = 'weight',
            weight_threshold=1.e-6)


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

### Sersic profile (detection mode)
if fit_case == "sersic_rg4" : 
    
    x,y = get_pos_parameters()
    
    rad = FreeParameter(lambda o: o.radius, Range(lambda v, o: (.0001, 1.5*v), RangeType.EXPONENTIAL))
    
    lrd=DependentParameter( lambda re: 1.015**(re - 10), rad )
    add_prior( lrd, 0.027/0.03,  0.5) 
    
    sersic = FreeParameter( 2.0, Range((0.3, 8.4), RangeType.LINEAR))
    X_sersic = DependentParameter( lambda n: np.log( (n-0.25)/(10-n) ), sersic )
    add_prior( X_sersic, -2.5, 1.5 )

    e1 = FreeParameter( 0.0, Range((-0.9999, 0.9999), RangeType.LINEAR))
    e2 = FreeParameter( 0.0, Range((-0.9999, 0.9999), RangeType.LINEAR))
    emod = DependentParameter( lambda x,y: np.sqrt( x*x + y*y ), e1, e2 )
    angle = DependentParameter( lambda e1,e2 : 0.5*np.arctan2( e1, e2 ), e1, e2 )
    ratio = DependentParameter( lambda e : np.abs(1-e)/(1+e), emod )
    add_prior( e1, 0.0, 0.25 )
    add_prior( e2, 0.0, 0.25 )

    ra, dec, wc_rad, wc_angle, wc_ratio = get_world_parameters(x, y, rad, angle, ratio)
    
    add_output_column('X_MODEL', x)
    add_output_column('Y_MODEL', y)
    add_output_column('RA_MODEL', ra)
    add_output_column('DEC_MODEL', dec)
    add_output_column('RADIUS', wc_rad)
    add_output_column('AXRATIO', wc_ratio)
    add_output_column('ANGLE', wc_angle)
    add_output_column('E1',e1)
    add_output_column('E2',e2)
    add_output_column('SERSIC', sersic)
    add_output_column('X_SERSIC', X_sersic)
    
    add_output_column('DET-IMG_RADIUS', rad)
    add_output_column('DET-IMG_AXRATIO', ratio)
    add_output_column('DET-IMG_ANGLE', angle)

    i = 0
    flux = {}
    mag = {}
    dx = {}
    xx={}
    dy = {}
    yy={}

    for band,group in mesgroup: 
        print(band)
        
        flux[i] = get_flux_parameter()
        mag[i] = DependentParameter(lambda f, zp=mag_zeropoint[band]: -2.5 * np.log10(f) + zp, flux[i] )
        
        # the centroid is fixed for all bands
        add_model(group, 
                  SersicModel(x, y, flux[i], rad, ratio, angle, sersic ) )

        add_output_column(f'MAG_MODEL_{band}', mag[i])
        add_output_column(f'FLUX_MODEL_{band}', flux[i])
        i+=1

### Bulge + Disc (detection mode)
if fit_case == "B+D":
    det_pix_scale = 0.04 # Pixel scale of detection image
    x,y = get_pos_parameters()

    r_b = FreeParameter(lambda o: o.radius, Range(lambda v,o: (0.001*v, 1.1*v), RangeType.EXPONENTIAL))
    r_d = FreeParameter(lambda o: o.radius, Range(lambda v,o: (0.01*v, 3*v), RangeType.EXPONENTIAL))
#     r_b = FreeParameter(lambda o: o.radius*0.5, Range(lambda v,o: (0.01*v, 2.1*v), RangeType.EXPONENTIAL))
#     r_d = FreeParameter(lambda o: o.radius*2.0, Range(lambda v,o: (0.01*v, 1.1*v), RangeType.EXPONENTIAL))
    
    lrd = DependentParameter( lambda y : np.log10(y), r_d )
    add_prior( lrd, 0.33+np.log10(0.1/det_pix_scale), 0.25 ) ## log10(rd) in pixels 
    
    rel_size = DependentParameter( lambda x,y : np.log10(y)-(1.14*np.log10(x)-1.2), r_b, r_d )
    add_prior( rel_size, 0.0, 0.4 )
    
    angle = FreeParameter( lambda o: o.angle )
    
    ratio_d = FreeParameter( 0.5, Range( (0.1, 1.0), RangeType.LINEAR) )
    X_rd = DependentParameter( lambda q: np.log( (q-0.01)/(1.01-q) ), ratio_d )
    add_prior( X_rd, 0.03, 1.0 )
    
    ratio_b = FreeParameter( 0.6, Range( (0.3, 1.0), RangeType.LINEAR) )
    X_rb = DependentParameter( lambda q: np.log( (q-0.01)/(1.01-q) ), ratio_b )
    add_prior( X_rb, 0.50, 1.1 )
    
    ra, dec, wc_rad_b, wc_angle, wc_ratio_b = get_world_parameters(x, y, r_b, angle, ratio_b)
    ra, dec, wc_rad_d, wc_angle, wc_ratio_d = get_world_parameters(x, y, r_d, angle, ratio_d)
    
    
    add_output_column('X_MODEL', x)
    add_output_column('Y_MODEL', y)
    add_output_column('RA_MODEL', ra)
    add_output_column('DEC_MODEL', dec)
    add_output_column('DISK_RADIUS_pix', r_d)
    add_output_column('BULGE_RADIUS_pix', r_b)
    add_output_column('DISK_RADIUS_deg', wc_rad_d)
    add_output_column('BULGE_RADIUS_deg', wc_rad_b)
    add_output_column('ANGLE_DetImg', angle)
    add_output_column('ANGLE', wc_angle)
    add_output_column('DISK_AXRATIO', wc_ratio_d)
    add_output_column('BULGE_AXRATIO', wc_ratio_b)
    add_output_column('DISK_AXRATIO_DetImg', ratio_d)
    add_output_column('BULGE_AXRATIO_DetImg', ratio_b)


        
    def func_X_BT_med(lambda_um):
        a = 0.59009933
        b = 0.35564443
        c = -1.99186112
        return a * np.log(b * lambda_um) + c

    def func_X_BT_wid(lambda_um):
        a = 0.85095641
        b = 0.68764073
        c = 2.87822699
        return a * np.log(b * lambda_um) + c
    
    i = 0
    flux1 = {}
    flux2 = {}
    mag1 = {}
    mag2 = {}
    flux = {}
    bt = {}
    X_bt= {}
    mag = {}

    for band,group in mesgroup:
        
        flux[i] = get_flux_parameter()
        mag[i] = DependentParameter(lambda f, zp=mag_zeropoint[band]: -2.5 * np.log10(f) + zp, flux[i] )
        
        bt[i] = FreeParameter(0.1, Range((0.00005,1.0), RangeType.LINEAR))
        
        flux1[i] = DependentParameter(lambda f, r: f*r, flux[i], bt[i] )
        flux2[i] = DependentParameter(lambda f, r: f*(1.0-r), flux[i], bt[i] )
        
        mag1[i] = DependentParameter(lambda f, zp=mag_zeropoint[band]: -2.5 * np.log10(f) + zp, flux1[i] )
        mag2[i] = DependentParameter(lambda f, zp=mag_zeropoint[band]: -2.5 * np.log10(f) + zp, flux2[i] )
        
        X_bt[i] = DependentParameter(lambda r: np.log( (r+0.01)/(1.01-r) ), bt[i] )
        add_prior( X_bt[i], func_X_BT_med(filters_waveband[band]['pivot']), func_X_BT_wid(filters_waveband[band]['pivot']) )
        
        add_model( group, ExponentialModel( x, y, flux2[i], r_d, ratio_d, angle) )
        add_model( group, DeVaucouleursModel( x, y, flux1[i], r_b, ratio_b, angle) )
        
        
        add_output_column(f'FLUX_MODEL_{band}',  flux[i])
        add_output_column(f'MAG_MODEL_{band}',  mag[i])
        
        add_output_column(f'FLUX_MODEL_BULGE_{band}',  flux1[i])
        add_output_column(f'MAG_MODEL_BULGE_{band}',  mag1[i])
        
        add_output_column(f'FLUX_MODEL_DISK_{band}',  flux2[i])
        add_output_column(f'MAG_MODEL_DISK_{band}',  mag2[i])
        
        add_output_column(f'B/T_{band}', bt[i])
        add_output_column(f'X_B/T_{band}', X_bt[i])
        i += 1

### Sersic profile (association mode)
if fit_case == "sersic_full_assoc" :
    # 1  : x           -> o.centroid_x
    # 2  : y           -> o.centroid_y
    # 3  : group_id    -> o.assoc_value_2
    # 4  : flux_mean   -> o.assoc_value_3
    # 5  : mag_mean    -> o.assoc_value_4
    # 6  : a_image     -> o.assoc_value_5
    # 7  : b_image     -> o.assoc_value_6
    # 8  : theta_image -> o.assoc_value_7
    # 9  : ax_ratio    -> o.assoc_value_8
    # 10 : bigsize     -> o.assoc_value_9

    coord_param_range = Range(lambda v, o: (v - 1*o.assoc_value_5, v + 1*o.assoc_value_5), RangeType.LINEAR)
    x = FreeParameter(lambda o: o.centroid_x, coord_param_range) 
    y = FreeParameter(lambda o: o.centroid_y, coord_param_range) 
    
    # rad = FreeParameter(lambda o: 1.3*o.assoc_value_5, Range(lambda v, o: (v*0.01, 5*v), RangeType.EXPONENTIAL))
    rad = FreeParameter(lambda o: 1.3*o.assoc_value_5, Range(lambda v, o: (v*0.01, 100*v), RangeType.EXPONENTIAL))
    
    lrd=DependentParameter( lambda re: 1.015**(re - 10), rad )
    add_prior( lrd, 0.027/0.03,  0.5) 

    if True:
        sersic = FreeParameter( 2.0, Range((0.3, 8.4), RangeType.LINEAR))
        X_sersic = DependentParameter( lambda n: np.log( (n-0.25)/(10-n) ), sersic )
        add_prior( X_sersic, -2.5, 1.5 )
    if False:
        sersic = FreeParameter( 2.0, Range((0.3, 5.0), RangeType.LINEAR))
        X_sersic = DependentParameter( lambda n: np.log( (n-0.25)/(6-n) ), sersic )
        add_prior( X_sersic, -2.5, 1.5 )

    e1 = FreeParameter( 0.0, Range((-0.9999, 0.9999), RangeType.LINEAR))
    e2 = FreeParameter( 0.0, Range((-0.9999, 0.9999), RangeType.LINEAR))
    emod = DependentParameter( lambda x,y: np.sqrt( x*x + y*y ), e1, e2 )
    angle = DependentParameter( lambda e1,e2 : 0.5*np.arctan2( e1, e2 ), e1, e2 )
    ratio = DependentParameter( lambda e : np.abs(1-e)/(1+e), emod )
    add_prior( e1, 0.0, 0.25 )
    add_prior( e2, 0.0, 0.25 )

    ra, dec, wc_rad, wc_angle, wc_ratio = get_world_parameters(x, y, rad, angle, ratio)
    
    add_output_column('X_MODEL', x)
    add_output_column('Y_MODEL', y)
    add_output_column('RA_MODEL', ra)
    add_output_column('DEC_MODEL', dec)
    add_output_column('RADIUS', wc_rad)
    add_output_column('AXRATIO', wc_ratio)
    add_output_column('ANGLE', wc_angle)
    add_output_column('E1',e1)
    add_output_column('E2',e2)
    add_output_column('SERSIC', sersic)
    add_output_column('X_SERSIC', X_sersic)
    
    add_output_column('DET-IMG_RADIUS', rad)
    add_output_column('DET-IMG_AXRATIO', ratio)
    add_output_column('DET-IMG_ANGLE', angle)

    i = 0
    flux = {}
    mag = {}
    dx = {}
    xx={}
    dy = {}
    yy={}

    for band,group in mesgroup: 
        print(band)
        
        flux[i] = FreeParameter(lambda o: o.assoc_value_3)
        # flux[i] = FreeParameter(lambda o, zp=mag_zeropoint[band]: 10**(0.4*(zp - o.assoc_value_4)))
        mag[i] = DependentParameter(lambda f, zp=mag_zeropoint[band]: -2.5 * np.log10(f) + zp, flux[i] )
        
        # the centroid is fixed for all bands
        add_model(group, 
                  SersicModel(x, y, flux[i], rad, ratio, angle, sersic ) )

        add_output_column(f'MAG_MODEL_{band}', mag[i])
        add_output_column(f'FLUX_MODEL_{band}', flux[i])
        i+=1

# we can print some info about the models        
print_model_fitting_info(mesgroup, show_params=True, prefix='')