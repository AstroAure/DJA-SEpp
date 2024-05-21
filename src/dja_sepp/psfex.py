import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from astropy.io import fits

# To work, this file requires you to have PSFEx installed on your machine.

def run_psfex(cat_name, 
              output_cat, 
              dir_chckimg, dir_psfex, 
              config_file, 
              fwhm_range='1.0,6.0', min_snr=200, 
              save_psf_png=False, verbose=False):
    """
    Run PSFEx on the FITS_LDAC catalog given, with vignets to create a PSF.

    cat_name : filename of the catalog
    output_cat : filename of the output catalog
    dir_chckimg : path to folder to save checkimages (will be created if non-existing)
    dir_psfex : path to folder to save the PSF (will be created if non-existing)
    config_file : filename of the PSFEx configuration file
    fwhm_range : range of FWHM (in pixels) used for PSFEx autoselect. 
                 Can be large if using a point-source only catalog
    min_snr : minimum SNR for PSFEx autoselect
    save_psf_png : save a png image of the PSF as checkimage
    verbose : verbose parameter
    """
    os.makedirs(dir_chckimg, exist_ok=True)
    os.makedirs(dir_psfex, exist_ok=True)
    study_name = ".".join(cat_name.split("/")[-1].split(".")[:-1])
    output_chk = f"{dir_chckimg}/{study_name}"
    verbose_type = 'NORMAL' if verbose else 'QUIET'
    subprocess.run(f'psfex {cat_name} -c {config_file} \
                   -CHECKIMAGE_NAME {dir_chckimg}/chi2.fits,{dir_chckimg}/samp.fits,{dir_chckimg}/res.fits,{dir_chckimg}/snap.fits \
                   -CHECKPLOT_NAME {output_chk}_selfwhm,{output_chk}_chi2,{output_chk}_counts \
                   -XML_NAME {output_chk}_xml.xml \
                   -OUTCAT_NAME {output_cat}_psf_cat.fits \
                   -PSF_DIR {dir_psfex} \
                   -SAMPLE_FWHMRANGE "{fwhm_range}" \
                   -SAMPLE_MINSN {min_snr} \
                   -VERBOSE_TYPE {verbose_type}',
                   shell=True)
    filename = f"{dir_psfex}/{study_name}_psf.psf"
    if save_psf_png:
        image = fits.open(filename)
        psf = image[1].data[0][0][0]
        image.close()
        fig, ax = plt.subplots(figsize=(6,6))
        ax.set_axis_off()
        ax.imshow(psf, origin='lower', cmap='bone_r', norm=SymLogNorm(linthresh=1e-5, linscale=0.8))
        fig.savefig(f"{filename.replace('.psf','.png')}", bbox_inches='tight', pad_inches=0, dpi=60)
        if verbose: plt.show()
    return filename
    
def compare_star(cat_name, cat_name_star, 
                 output_cat, output_cat_star, 
                 dir_chckimg, dir_psfex, config_file, 
                 fwhm_range='1.0,6.0', fwhm_range_star='1.0,10.0',
                 min_snr=200, min_snr_star=10,
                 verbose=False):
    """
    Function to compare the PSF between PSFEx standard auto-selection of point-like sources
    and a given point-like sources catalog.

    cat_name : filename of the raw SExtractor catalog
    cat_name_star : filename of the point-like sources catalog
    output_cat : filename of the output catalog
    output_cat_star : filename of the output point-like sources catalog
    dir_chckimg : path to folder to save checkimages (will be created if non-existing)
    dir_psfex : path to folder to save the PSF (will be created if non-existing)
    config_file : filename of the PSFEx configuration file
    fwhm_range : range of FWHM (in pixels) used for PSFEx autoselect
    fwhm_range_star : range of FWHM (in pixels) used for PSFEx autoselect (still used, with a large range)
    min_snr : minimum SNR for PSFEx autoselect
    min_snr_star : minimum SNR for PSFEx autoselect (still used, with a small value)
    verbose : verbose parameter
    """
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
    pass

if __name__=="__main__":
    main()