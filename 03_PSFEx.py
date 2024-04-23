import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from astropy.io import fits

def run_psfex(cat_name, output_cat, dir_chckimg, dir_psfex, config_file, fwhm_range='1.0,6.0', min_snr=200, save_psf_png=False, verbose=False):
    os.makedirs(dir_chckimg, exist_ok=True)
    os.makedirs(dir_psfex, exist_ok=True)
    study_name = ".".join(cat_name.split("/")[-1].split(".")[:-1])
    output_chk = f"{dir_chckimg}/{study_name}"
    verbose_type = 'NORMAL' if verbose else 'QUIET'
    os.system(f'psfex {cat_name} -c {config_file} \
                -CHECKIMAGE_NAME {dir_chckimg}/chi2.fits,{dir_chckimg}/samp.fits,{dir_chckimg}/res.fits,{dir_chckimg}/snap.fits \
                -CHECKPLOT_NAME {output_chk}_selfwhm,{output_chk}_chi2,{output_chk}_counts \
                -XML_NAME {output_chk}_xml.xml \
                -OUTCAT_NAME {output_cat}_psf_cat.fits \
                -PSF_DIR {dir_psfex} \
                -SAMPLE_FWHMRANGE "{fwhm_range}" \
                -SAMPLE_MINSN {min_snr} \
                -VERBOSE_TYPE {verbose_type}')
    filename = f"{dir_psfex}/{study_name}_psf.psf"
    if save_psf_png:
        image = fits.open(filename)
        psf = image[1].data[0][0][0]
        image.close()
        fig, ax = plt.subplots(figsize=(6,6))
        ax.set_axis_off()
        ax.imshow(psf, origin='lower', cmap='bone_r', norm=SymLogNorm(linthresh=1e-5, linscale=0.8))
        fig.savefig(f"{output_chk}_psf.png", bbox_inches='tight', pad_inches=0, dpi=60)
        if verbose: plt.show()
    return filename
    
def compare_star(cat_name, cat_name_star, 
                 output_cat, output_cat_star, 
                 dir_chckimg, dir_psfex, config_file, 
                 fwhm_range='1.0,6.0', fwhm_range_star='1.0,10.0',
                 min_snr=200, min_snr_star=10,
                 verbose=False):
    run_psfex(cat_name, output_cat, dir_chckimg, dir_psfex, config_file, fwhm_range=fwhm_range, min_snr=min_snr, verbose=verbose)
    run_psfex(cat_name_star, output_cat_star, dir_chckimg, dir_psfex, config_file, fwhm_range=fwhm_range_star, min_snr=min_snr_star, verbose=verbose)
    study_name = ".".join(cat_name.split("/")[-1].split(".")[:-1])
    study_name_star = ".".join(cat_name_star.split("/")[-1].split(".")[:-1])
    image = fits.open(f"{dir_psfex}/{study_name}_psf.psf")
    image_star = fits.open(f"{dir_psfex}/{study_name_star}_psf.psf")
    psf = np.array(image[1].data[0][0][0])
    psf_star = np.array(image_star[1].data[0][0][0])
    fig, ax = plt.subplots(1,2,figsize=(12,6))
    ax[0].imshow(psf, origin='lower', cmap='bone_r', norm=SymLogNorm(linthresh=1e-5, linscale=0.8))
    ax[1].imshow(psf_star, origin='lower', cmap='bone_r', norm=SymLogNorm(linthresh=1e-5, linscale=0.8))
    ax[0].set_title(f"PSFEx auto-selection \n SAMPLE_FWHMRANGE = ({fwhm_range})")
    ax[1].set_title("MU_MAX v. MAG_AUTO selection")
    plt.show()
    return fig


def main():
    config_file = "/home/aurelien/DAWN/DJA_SE++/config/psfex_default.conf"
    cat_name = "/home/aurelien/DAWN/DJA_SE++/JADES-catalog/gds-grizli-v7.0-f444w-clear_drc_cat.fits"
    cat_name_star = "/home/aurelien/DAWN/DJA_SE++/JADES-catalog/gds-grizli-v7.0-f444w-clear_drc_star_cat.fits"
    output_cat = "/home/aurelien/DAWN/DJA_SE++/JADES-catalog/gds-grizli-v7.0-f444w-clear_drc_psf_cat.fits"
    output_cat_star = "/home/aurelien/DAWN/DJA_SE++/JADES-catalog/gds-grizli-v7.0-f444w-clear_drc_star_psf_cat.fits"
    dir_chckimg = "/home/aurelien/DAWN/DJA_SE++/JADES-checkimages"
    dir_psfex = "/home/aurelien/DAWN/DJA_SE++/JADES-psfex"
    compare_star(cat_name, cat_name_star, output_cat, output_cat_star, dir_chckimg, dir_psfex, config_file)

if __name__=='__main__':
    main()