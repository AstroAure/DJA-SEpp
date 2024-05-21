import sys
import glob
import re
import dja_sepp

# python3 psf.py ceers-full-grizli-v7.2 /FlashStorage/DJA-SEpp /FlashStorage/DJA-SEpp/config aurelien-sepp

def main():
    field = sys.argv[1]
    home = sys.argv[2] if len(sys.argv)>2 else "/tmp"
    config_folder = sys.argv[3] if len(sys.argv)>3 else f"{home}/config"
    bucket = sys.argv[4] if len(sys.argv)>4 else 'aurelien-sepp'

    # F200W Point-like sources detection
    filter = 'f200w'
    dja_sepp.sextractor.extract_stars(detect_img      = glob.glob(f"{home}/fields/{field}/image/*{filter}*sci*.fits")[0], \
                                      weight_img      = glob.glob(f"{home}/fields/{field}/image/*{filter}*wht*.fits")[0], \
                                      output_cat      = f"{home}/fields/{field}/catalog/{field}_drc_cat.fits", \
                                      output_cat_star = f"{home}/fields/{field}/catalog/{field}_drc_cat_star.fits", \
                                      config_folder   = config_folder, \
                                      dir_chckimg     = f"{home}/fields/{field}/catalog/checkimages", \
                                      detect_thresh = 8.0, \
                                      y_max = -5.0, \
                                      save_chckimg = False, plot = False, clean = True, verbose = True)
    # Save star catalog to S3
    dja_sepp.s3.save_s3(f"{home}/fields/{field}/catalog/{field}_drc_cat_star.fits", bucket, f"{field}/catalog")
    # Filter list
    filter_list = [re.search('(f\d+\w+)', filename).group(1) for filename in glob.glob(f"{home}/fields/{field}/image/*clear*sci*")]
    filter_list.sort()
    # Run SExtractor, cross-match and PSFEx for each band
    for filter in filter_list:
        print(filter.upper())
        print("Running SExtractor")
        dja_sepp.sextractor.extract_stars_catalog(detect_img      = glob.glob(f"{home}/fields/{field}/image/*{filter}*sci.fits")[0], \
                                                  weight_img      = glob.glob(f"{home}/fields/{field}/image/*{filter}*wht.fits")[0], \
                                                  output_cat      = f"{home}/fields/{field}/catalog/{field}-{filter}-clear_drc_cat.fits", \
                                                  output_cat_star = f"{home}/fields/{field}/catalog/{field}-{filter}-clear_drc_cat_star.fits", \
                                                  input_cat_star  = f"{home}/fields/{field}/catalog/{field}_drc_cat_star.fits", \
                                                  config_folder   = config_folder, \
                                                  dir_chckimg     = f"{home}/fields/{field}/catalog/checkimages", \
                                                  detect_thresh = 2.5, \
                                                  save_chckimg = False, plot = False, clean = True, run_sex = True, verbose = True)
        print("Running PSFEX")
        psf = dja_sepp.psfex.run_psfex(cat_name    = f"{home}/fields/{field}/catalog/{field}-{filter}-clear_drc_cat_star.fits",
                                       output_cat  = f"{home}/fields/{field}/catalog/{field}-{filter}-clear_drc_cat_star_psf.fits",
                                       dir_chckimg = f"{home}/fields/{field}/catalog/checkimages",
                                       dir_psfex   = f"{home}/fields/{field}/psfex",
                                       config_file = f"{config_folder}/psfex_default.conf",
                                       fwhm_range  = '1.0, 10.0',
                                       save_psf_png = True,
                                       verbose = True)
        print("Save PSF to S3")
        dja_sepp.s3.save_s3(psf, 'aurelien-sepp', f"{field}/psfex")
        dja_sepp.s3.save_s3(psf.replace('.psf','.png'), 'aurelien-sepp', f"{field}/psfex")
        

if __name__=='__main__':
    main()