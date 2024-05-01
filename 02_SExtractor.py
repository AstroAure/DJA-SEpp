import os
import glob
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
import astromatic_wrapper as aw
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.colors import SymLogNorm
from matplotlib.patches import Rectangle
import random as rd
from sklearn.cluster import DBSCAN
from sklearn.linear_model import RANSACRegressor

def run_sextractor(detect_img      : str, weight_img   : str, 
                   output_cat      : str, 
                   config_file     : str, params_file  : str,
                   dir_chckimg     : str,
                   detect_thresh   : float = 12.0,
                   analysis_thresh : float = 2.5,
                   verbose         : bool = False) -> None :
    os.makedirs("/".join(output_cat.split("/")[:-1]), exist_ok=True)
    os.makedirs(dir_chckimg, exist_ok=True)
    study_name = '.'.join(detect_img.split('/')[-1].split('.')[:-1])
    seg_img = f'{dir_chckimg}/{study_name}_seg.fits'
    verbose_type = 'NORMAL' if verbose else 'QUIET'
    os.environ['LD_LIBRARY_PATH'] = '/usr/lib' # Link 'libcfitsio.so.10'
    os.system(f'sex {detect_img} -c {config_file} \
            -CATALOG_NAME {output_cat} \
            -WEIGHT_IMAGE {weight_img} -WEIGHT_TYPE MAP_WEIGHT \
            -PARAMETERS_NAME {params_file} \
            -CHECKIMAGE_NAME {seg_img} \
            -DETECT_THRESH {detect_thresh} \
            -ANALYSIS_THRESH {analysis_thresh} \
            -VERBOSE_TYPE {verbose_type}')
    return output_cat
    
def show_vignets(data, grid_n):
    n = len(data)
    a, b = grid_n
    selected = rd.sample(range(n), a*b)
    fig, ax = plt.subplots(a,b,figsize=(4*b,4*a))
    for i in range(a):
        for j in range(b):
            ax[i,j].imshow(data['VIGNET'][selected[i*b+j]], 
                           origin='lower', cmap=colormaps['bone_r'], 
                           norm=SymLogNorm(linthresh=1e-2, linscale=0.1, vmin=-0.1, vmax=0.5))
            ax[i,j].set_title(selected[i*b+j])
    plt.show()

def plot_MuvMAG(data, star_selections=[], mag_bounds=(17,30), mu_bounds=(13,25), ax_custom=None):
    # star_selection = {<name> : {'label':<label>, 'color':<color>, 'flag':<array_of_indices>}}
    # MU_MAX = max pixel value in mag/arcsec² (MAG_ZEROPOINT - 2.5*log10(FLUX_MAX*PIXELSCALE²))
    # MAG_AUTO = total flux in mag (MAG_ZEROPOINT - 2.5*log10(FLUX_AUTO))
    ax = plt.subplots(figsize=(8,6))[1] if ax_custom is None else ax_custom
    ax.hexbin(data['MAG_AUTO'], data['MU_MAX'], bins='log', mincnt=1, extent=(mag_bounds[0],mag_bounds[1],mu_bounds[0],mu_bounds[1]), cmap='inferno', lw=0.01)
    # ax.plot(data['MAG_AUTO'], data['MU_MAX'], marker='o', ms=0.5, ls="", alpha=0.5, c='k', label=f'All [{len(data)}]', rasterized=True)
    for name in star_selections:
        selection = star_selections[name]
        ax.plot(data[selection['flag']]['MAG_AUTO'], data[selection['flag']]['MU_MAX'], 
                   marker='o', ms=1.0, ls="", c=selection['color'], label=f"{selection['label']} [{len(selection['flag'])}]",
                   rasterized=True)
    ax.set_xlabel('MAG_AUTO')
    ax.set_ylabel('MU_MAX')
    ax.set_xlim(mag_bounds[0],mag_bounds[1])
    ax.set_ylim(mu_bounds[0],mu_bounds[1])
    ax.legend(loc='best')
    if ax_custom is None: plt.show()

def hist_CLASS_STAR(data, star_selections=[], hist_bound=500, ax_custom=None):
    # star_selection = {<name> : {'label':<label>, 'color':<color>, 'flag':<array_of_indices>}}
    ax = plt.subplots(figsize=(8,6))[1] if ax_custom is None else ax_custom
    ax.hist(data['CLASS_STAR'], bins=100, color='k', label=f'All [{len(data)}]')
    for name in star_selections:
        selection = star_selections[name]
        ax.hist(data[selection['flag']]['CLASS_STAR'], bins=100, range=(0,1),
                color=selection['color'], label=f"{selection['label']} [{len(selection['flag'])}]")
    ax.set_ylabel('CLASS_STAR')
    ax.set_ylim(0,hist_bound)
    ax.legend(loc='best')
    if ax_custom is None: plt.show()

def plot_SNR_radius(data, star_selections=[], rad_bound=20, ax_custom=None):
    ax = plt.subplots(figsize=(8,6))[1] if ax_custom is None else ax_custom
    ax.hexbin(data['FLUX_RADIUS'], data['SNR_WIN'], bins='log', mincnt=1, yscale='log', extent=(0,rad_bound,0,6), cmap='inferno', lw=0.01)
    # ax.plot(data['FLUX_RADIUS'], data['SNR_WIN'], marker='o', ms=0.5, ls="", alpha=0.5, c='k', label=f'All [{len(data)}]', rasterized=True)
    for name in star_selections:
        selection = star_selections[name]
        ax.plot(data[selection['flag']]['FLUX_RADIUS'], data[selection['flag']]['SNR_WIN'],
                   marker='o', ms=1.0, ls="", c=selection['color'], label=f"{selection['label']} [{len(selection['flag'])}]",
                   rasterized=True)
    ax.set_xlabel('FLUX_RADIUS')
    ax.set_ylabel('SNR_WIN')
    ax.set_xlim(0, rad_bound)
    ax.set_yscale('log')
    ax.legend()
    if ax_custom is None: plt.show()

def plot_MuvMAG_manual(data, mag_bounds=(19,28), 
                       ax_custom=None, plot_mag_bounds=(18,30), plot_y_bounds=(-8,1)):
    """ 
    Plots MU_MAX-MAG_AUTO v. MAG_AUTO with
    guides to manually select the star line.
    """
    ax = plt.subplots(figsize=(8,6))[1] if ax_custom is None else ax_custom
    ax.hexbin(data['MAG_AUTO'], data['MU_MAX'], bins='log', mincnt=1, extent=(plot_mag_bounds[0],plot_mag_bounds[1],plot_y_bounds[0],plot_y_bounds[1]), cmap='inferno', lw=0.01)
    # ax.plot(data['MAG_AUTO'], data['MU_MAX']-data['MAG_AUTO'], marker='o', ls="", ms=0.5, c='k', rasterized=True)
    ax.axvline(mag_bounds[0], color='r', linewidth=0.5)
    ax.axvline(mag_bounds[1], color='r', linewidth=0.5)
    ax.set_xlabel('MAG_AUTO')
    ax.set_ylabel('MU_MAX - MAG_AUTO')
    ax.set_xlim(plot_mag_bounds[0], plot_mag_bounds[1])
    ax.set_ylim(plot_y_bounds[0], plot_y_bounds[1])
    ax.grid(axis='y', color='r', linewidth=0.5, linestyle=(0, (5, 5)), alpha=1)
    ax.set_axisbelow(True)
    if ax_custom is None: plt.show()

def find_star_line(data, eps_DBSCAN=0.1, y_max=-4, mag_fit=25, verbose=False, 
                   plot=False, plot_mag_bounds=(18,30), plot_y_bounds=(-8,1),
                   save=False, save_name=None):
    """
    Finds the star line on the MU_MAX v. MAG_AUTO plot.
    1. Remove the main galaxy cluster using DBSCAN (eps_DBSCAN parameter).
    2. RANSAC regression to find the star line on the remaining points,
        and with MU_MAX-MAG_AUTO < y_max
    3. Calculates the star line position at MAG_AUTO = mag_fit
    """
    # Remove galaxy cluster
    values = np.vstack([data['MAG_AUTO'], data['MU_MAX']-data['MAG_AUTO']])
    clustering = DBSCAN(eps=eps_DBSCAN).fit(values.T)
    lab, count = np.unique(clustering.labels_, return_counts=True)
    if verbose : print(f'DBSCAN clustering : {lab, count}')
    # Find star line
    cloud_idx = np.argmax(count[lab!=-1])
    A = data[clustering.labels_ != cloud_idx]
    B = A[A['MU_MAX']-A['MAG_AUTO'] < y_max]
    if plot | save:
        fig, ax = plt.subplots(figsize=(8,6))
        ax.plot(A['MAG_AUTO'], A['MU_MAX']-A['MAG_AUTO'], ms=0.5, ls="", marker='o', c='k', rasterized=True)
        ax.axhline(y_max, color='r', linewidth=0.5, linestyle='--')
        ax.set_xlim(plot_mag_bounds[0], plot_mag_bounds[1])
        ax.set_ylim(plot_y_bounds[0], plot_y_bounds[1])
        ax.set_xlabel('MAG_AUTO')
        ax.set_ylabel('MU_MAX - MAG_AUTO')
        ax.set_title("Star line detection after galaxy cluster removal")
        if save : fig.savefig(save_name, bbox_inches='tight', dpi=100)
        if plot : plt.show()
    ransac = RANSACRegressor()
    ransac.fit(B['MAG_AUTO'].reshape(-1,1), (B['MU_MAX']-B['MAG_AUTO']).reshape(-1,1))
    if (slope:=abs(ransac.estimator_.coef_[0][0])) > 0.1: 
        print(f"\33[33m WARNING \33[0m : Star line may be wrong (Slope = {slope:.3f})")
    c = ransac.predict([[mag_fit]])[0,0]
    if verbose: print(f"{'RANSAC Slope':<13}: {ransac.estimator_.coef_[0][0]:.3f}")
    if verbose: print(f"{'Star line':<13}: {c:.3f}")
    return c

def MUvMAG_star_selection(data, star_line, y_offsets=(-0.5,0.3), mag_bounds=(19,28), snr_min=1e2, 
                          plot=False, plot_mag_bounds=(18,30), plot_y_bounds=(-8,1),
                          save=False, save_name=None):
    """
    Selects the stars according to the MU_MAX v. MAG_AUTO box.
    star_line : Star line value of MU_MAX - MAG_AUTO
    y_offsets : Negative and positive offsets from star_line for the selection box
    mag_bounds : MAG_AUTO bounds for the selection box
    """
    if plot | save:
        fig, ax = plt.subplots(figsize=(8,6))
        ax.plot(data['MAG_AUTO'], data['MU_MAX']-data['MAG_AUTO'], ms=0.5, ls="", marker='o', c='k', rasterized=True)
        ax.add_patch(Rectangle((19, star_line + y_offsets[0]), 
                                mag_bounds[1]-mag_bounds[0], 
                                y_offsets[1]-y_offsets[0], 
                                fill=False, edgecolor='r', lw=1))
        ax.set_xlabel('MAG_AUTO')
        ax.set_ylabel('MU_MAX - MAG_AUTO')
        ax.set_xlim(plot_mag_bounds[0], plot_mag_bounds[1])
        ax.set_ylim(plot_y_bounds[0], plot_y_bounds[1])
        ax.set_title("Star detection with MUvMAG")
        if save : fig.savefig(save_name, bbox_inches='tight', dpi=100)
        if plot : plt.show()
    star_MUvMAG = np.where((data['MU_MAX'] < data['MAG_AUTO'] + star_line + y_offsets[1]) & \
                           (data['MU_MAX'] > data['MAG_AUTO'] + star_line + y_offsets[0]) & \
                           (data['MAG_AUTO'] < mag_bounds[1]) & \
                           (data['MAG_AUTO'] > mag_bounds[0]) & \
                           (data['SNR_WIN'] > snr_min))[0]
    return star_MUvMAG

def save_catalog(hdul, selection, outcat_name):
    hdul[2].data = hdul[2].data[selection]
    hdul[2].update_header()
    hdul.writeto(outcat_name, overwrite=True)

def extract_stars(detect_img        : str, 
                  weight_img        : str, 
                  output_cat        : str,
                  output_cat_star   : str, 
                  config_file       : str, 
                  params_file       : str,
                  dir_chckimg       : str,
                  detect_thresh     : float = 8.0,
                  analysis_thresh   : float = 2.5,
                  eps_DBSCAN        : float = 0.1, 
                  y_max             : float = -4.5, 
                  mag_fit           : float = 25,
                  y_offsets         : tuple = (-0.5,0.3),
                  mag_bounds        : tuple = (19,28),
                  snr_min           : float = 1e2,
                  plot_mag_bounds   : tuple = (18,30),
                  plot_mu_bounds    : tuple = (13,25), 
                  plot_y_bounds     : tuple = (-8,1),
                  plot_rad_bounds   : tuple = 20,
                  save_chckimg      : bool = True,
                  plot              : bool = False,
                  clean             : bool = True,
                  run_sex           : bool = True,
                  verbose           : bool = False) -> None :
    study_name = '.'.join(detect_img.split('/')[-1].split('.')[:-1])
    if run_sex: run_sextractor(detect_img, weight_img, output_cat, config_file, params_file, dir_chckimg, detect_thresh, analysis_thresh, verbose)
    hdul = fits.open(output_cat)
    data = hdul[2].data
    # data = aw.utils.ldac.get_table_from_ldac(output_cat, frame=1)
    star_line = find_star_line(data, eps_DBSCAN, y_max, mag_fit, verbose, plot, plot_mag_bounds, plot_y_bounds, save_chckimg, f"{dir_chckimg}/starDetect_{study_name}.png")
    star_MUvMAG = MUvMAG_star_selection(data, star_line, y_offsets, mag_bounds, snr_min, plot, plot_mag_bounds, plot_y_bounds, save_chckimg, f"{dir_chckimg}/starLine_{study_name}.png")
    if plot | save_chckimg:
        star_selections = {'MUvMAG' : {'label': 'Stars (MU v. MAG)', 'color': 'r', 'flag': star_MUvMAG}}
        fig, ax = plt.subplots(1,2,figsize=(12,6))
        plot_MuvMAG(data, star_selections, plot_mag_bounds, plot_mu_bounds, ax_custom=ax[0])
        plot_SNR_radius(data, star_selections, plot_rad_bounds, ax_custom=ax[1])
        if save_chckimg : fig.savefig(f"{dir_chckimg}/star_{study_name}.png", bbox_inches='tight', dpi=100)
        if plot : plt.show()
    if clean: os.remove(output_cat)
    if clean: os.remove(f"{dir_chckimg}/{study_name}_seg.fits")
    save_catalog(hdul, star_MUvMAG, output_cat_star)
    hdul.close()
    # aw.utils.ldac.save_table_as_ldac(Table(data[star_MUvMAG]), output_cat_star, overwrite=True)

def extract_stars_catalog(detect_img        : str, 
                          weight_img        : str,
                          output_cat        : str, 
                          output_cat_star   : str,
                          input_cat_star    : str, 
                          config_file       : str, 
                          params_file       : str,
                          dir_chckimg       : str,
                          detect_thresh     : float = 2.5,
                          analysis_thresh   : float = 2.5,
                          max_sep           : u.Quantity = 1.0*u.arcsec,
                          plot_mag_bounds   : tuple = (18,30),
                          plot_mu_bounds    : tuple = (13,25),
                          plot_rad_bounds   : tuple = 20, 
                          save_chckimg      : bool = True,
                          plot              : bool = False,
                          clean             : bool = True,
                          run_sex           : bool = True,
                          verbose           : bool = False) -> None:
    study_name = '.'.join(detect_img.split('/')[-1].split('.')[:-1])
    if run_sex: run_sextractor(detect_img, weight_img, output_cat, config_file, params_file, dir_chckimg, detect_thresh, analysis_thresh, verbose)
    if verbose: print ("Opening star catalogs")
    hdul = fits.open(output_cat, memmap=True, mode='denywrite') 
    data = hdul[2].data
    with fits.open(input_cat_star, memmap=True, mode='denywrite') as hdul_star: 
        data_star = hdul_star[2].data
    # data = aw.utils.ldac.get_table_from_ldac(output_cat, frame=1)
    # data_star = aw.utils.ldac.get_table_from_ldac(input_cat_star, frame=1)
    coord = SkyCoord(data['ALPHA_J2000'][:], data['DELTA_J2000'][:], unit='deg', frame='icrs')
    coord_star = SkyCoord(data_star['ALPHA_J2000'], data_star['DELTA_J2000'], unit='deg', frame='icrs')
    if verbose: print("Matching stars")
    idx, d2d, _ = coord_star.match_to_catalog_sky(coord)
    match = idx[d2d<max_sep]
    if plot | save_chckimg:
        print("Plotting")
        star_selections = {'MUvMAG' : {'label': 'Stars (MU v. MAG)', 'color': 'r', 'flag': match}}
        fig, ax = plt.subplots(1,2,figsize=(12,6))
        plot_MuvMAG(data, star_selections, plot_mag_bounds, plot_mu_bounds, ax_custom=ax[0])
        plot_SNR_radius(data, star_selections, plot_rad_bounds, ax_custom=ax[1])
        if save_chckimg : fig.savefig(f"{dir_chckimg}/star_{study_name}.png", bbox_inches='tight', dpi=100)
        if plot : plt.show()
    if clean: os.remove(output_cat)
    save_catalog(hdul, match, output_cat_star)
    hdul.close()
    # hdul_star.close()
    # aw.utils.ldac.save_table_as_ldac(Table(data[match]), output_cat_star, overwrite=True)


def main():
    # extract_stars(detect_img      = "/home/aurelien/DAWN/DJA_SE++/image/GDS/gds-grizli-v7.0-f444w-clear_drc_sci.fits", \
    #               weight_img      = "/home/aurelien/DAWN/DJA_SE++/image/GDS/gds-grizli-v7.0-f444w-clear_drc_wht.fits", \
    #               output_cat      = "/home/aurelien/DAWN/DJA_SE++/JADES-catalog/gdn-grizli-v7.0-f444w-clear_drc_cat.fits", \
    #               output_cat_star = "/home/aurelien/DAWN/DJA_SE++/JADES-catalog/gdn-grizli-v7.0-f444w-clear_drc_cat_star.fits", \
    #               config_file     = "/home/aurelien/DAWN/DJA_SE++/config/PSFEx-Cat-JWST.sex", \
    #               params_file     = "/home/aurelien/DAWN/DJA_SE++/config/PSFEx-Cat-JWST-SW.param", \
    #               dir_chckimg     = "/home/aurelien/DAWN/DJA_SE++/JADES-checkimages", \
    #               y_max = -4.5, \
    #               save_chckimg = True, plot = False, verbose = True)

    field = 'GDS'
    filter = 'f200w'
    extract_stars_catalog(detect_img      = glob.glob(f"/FlashStorage/image/{field}/*{filter}*sci.fits")[0], \
                                  weight_img      = glob.glob(f"/FlashStorage/image/{field}/*{filter}*wht.fits")[0], \
                                  output_cat      = f"/FlashStorage/catalog/{field}/{field}-{filter}-clear_drc_cat.fits", \
                                  output_cat_star = f"/FlashStorage/catalog/{field}/{field}-{filter}-clear_drc_cat_star.fits", \
                                  input_cat_star  = f"/FlashStorage/catalog/{field}/{field}_drc_cat_star.fits", \
                                  config_file     = "/home/ec2-user/DAWN/DJA-SEpp/config/PSFEx-Cat-JWST.sex", \
                                  params_file     = "/home/ec2-user/DAWN/DJA-SEpp/config/PSFEx-Cat-JWST-SW.param", \
                                  dir_chckimg     = f"/FlashStorage/checkimages/{'/'.join(field.split('/')[:-1])}", \
                                  detect_thresh = 2.5, \
                                  save_chckimg = False, plot = False, clean = False, verbose = True)
    
if __name__=="__main__":
    main()