import os
import subprocess
import glob
import re
import sys
import boto3
from astropy.io import fits
import numpy as np
from . import utils

def find_images(generic_img : str,
                generic_wht : str,
                generic_psf : str,
                filter_list = None,
                verbose : bool = False):
    """Find the images, weight maps and PSFs according to generic names
    
    generic_img : generic filename for the images (can use *,? wildcards understandable by glob.glob)
    generic_wht : generic filename for the weight maps (can use *,? wildcards understandable by glob.glob)
    generic_psf : generic filename for the PSFs (can use *,? wildcards understandable by glob.glob)
    filter_list (optional) : list of filters to find in the file names.
        If set to None, will automatically be created by looking for f\d+\w+ Regex in images names
    verbose (default, False)
        
    Returns : list of images names, list of weight maps names, list of PSFs names"""
    images = glob.glob(generic_img)
    weights = glob.glob(generic_wht)
    psfs = glob.glob(generic_psf)
    if filter_list is None:
        filter_list = list(set([re.search('(f\d+\w+)', img).group(1) for img in images]))
        filter_list.sort()
    images = [[file for file in images if filter in file][0] for filter in filter_list]
    for i, img in enumerate(images):
        with fits.open(img, memmap=True) as hdul:
            if np.max(hdul[0].data)==0.0:
                blank_img = images.pop(i)
                blank_filter = filter_list.pop(i)
                print(f"Filter {blank_filter.upper()} is blank !! ({blank_img})")
    weights = [[file for file in weights if filter in file][0] for filter in filter_list]
    psfs = [[file for file in psfs if filter in file][0] for filter in filter_list]
    return images, weights, psfs

def run_sepp(detect_img : str,
             generic_img : str,
             generic_wht : str,
             generic_psf : str,
             output_catname : str,
             checkimg_path : str,
             python_config : str,
             filt : str,
             filter_list = None,
             fit_case : str = 'sersic_rg4',
             properties : str = 'PixelCentroid,WorldCentroid,SourceIDs,GroupInfo,GroupStamp,SourceFlags,NDetectedPixels,NCorePixel,AperturePhotometry,AutoPhotometry,FluxRadius,SNRRatio,ShapeParameters,FlexibleModelFitting',
             thread_count : int = 128,
             tile_memory_limit : int = 16384,
             tile_size : int = 1000,
             verbose : bool = True):
    """
    Run SourceXtractor++ on a given set of images
    
    detect_img : filename of the detection image to use
    generic_img : generic filename for the images (can use *,? wildcards understandable by glob.glob)
    generic_wht : generic filename for the weight maps (can use *,? wildcards understandable by glob.glob)
    generic_psf : generic filename for the PSFs (can use *,? wildcards understandable by glob.glob)
    output_catname : filename of the output catalog
    checkimg_path : path to the folder for the model and residual images
    python_config : filename of the Python configuration file
    filt : convolution filter filename to use for the detection
    filter_list (optional) : list of filters to find in the file names.
        If set to None, will automatically be created by looking for f\d+\w+ Regex in images names
    fit_case (default, sersic_rg4) : the model for fitting in SE++
    properties (default, complete) : SE++ properties to compute and export
    thread_count (default, 128) : number of threads to use for SE++
    tile_memory_limit (default, 16384) : maximum RAM allocated to SE++ (in MiB)
    tile_size (default, 1000) : size in px of the tiles created by SE++
    verbose (default, True)
    """
    images, weights, psfs = find_images(generic_img, generic_wht, generic_psf, filter_list, verbose)
    args = {'fit_case' : fit_case,
            'list_of_IMG_names' : images,
            'list_of_WHT_names' : weights,
            'list_of_PSF_names' : psfs
            }
    os.makedirs(checkimg_path, exist_ok=True)
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['MKL_DYNAMIC'] = 'FALSE'
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['OMP_DYNAMIC'] = 'FALSE'
    subprocess.run(f"sourcextractor++ \
                   --python-config-file {python_config} \
                   --python-arg '{args}' \
                   --output-catalog-filename {output_catname} \
                   --output-properties {properties} \
                   --output-flush-size 10 \
                   --detection-image {detect_img} \
                   --weight-type WEIGHT \
                   --weight-absolute 0 \
                   --segmentation-filter {filt} \
                   --background-cell-size 128 \
                   --segmentation-algorithm LUTZ \
                   --smoothing-box-size 5 \
                   --detection-threshold 0.80 \
                   --detection-minimum-area 7 \
                   --partition-corethreshold yes \
                   --core-threshold-value 1.5 \
                   --core-minimum-area 9 \
                   --grouping-algorithm split \
                   --partition-multithreshold yes \
                   --partition-minimum-area 18 \
                   --grouping-hard-limit 30 \
                   --partition-minimum-contrast 0.0001 \
                   --partition-threshold-count 42 \
                   --use-cleaning yes \
                   --cleaning-minimum-area 15 \
                   --model-fitting-iterations 550 \
                   --sampling-scale-factor 1 \
                   --check-image-model-fitting {checkimg_path}/model.fits \
                   --check-image-residual {checkimg_path}/resid.fits \
                   --log-file {output_catname.replace('.fits', '.log')} \
                   --log-level DEBUG \
                   --thread-count {thread_count} \
                   --tile-memory-limit {tile_memory_limit} \
                   --tile-size {tile_size} \
                   ",
                   shell=True)

def main_tile(tile, img_dir, sepp_dir, psf_dir, config_dir):
    # Python configuration file for SE++
    python_config=f'{config_dir}/sepp-config.py'
    # Detection image
    detect_img=glob.glob(f'{img_dir}/*ir*sci*-{tile}.fits')[0]
    # Filter to use on the detection image
    filt=f'{config_dir}/gauss_1.5_3x3.conv'
    # Output
    name = ".".join(detect_img.split("/")[-1].split(".")[:-1])
    output_catname=f'{sepp_dir}/{name}_sepp_cat.fits'
    checkimg_path=f'{sepp_dir}/checkimages'
    os.makedirs(checkimg_path, exist_ok=True)
    # Run SE++
    run_sepp(detect_img=detect_img,
             generic_img=f"{img_dir}/*clear*sci*-{tile}.fits",
             generic_wht=f"{img_dir}/*clear*wht*-{tile}.fits",
             generic_psf=f"{psf_dir}/*star_psf.psf",
             output_catname=output_catname,
             checkimg_path=checkimg_path,
             python_config=python_config,
             filt=filt,
             filter_list = None,
             fit_case='sersic_rg4',
             properties='PixelCentroid,WorldCentroid,SourceIDs,GroupInfo,GroupStamp,SourceFlags,NDetectedPixels,NCorePixel,AperturePhotometry,AutoPhotometry,FluxRadius,SNRRatio,ShapeParameters,FlexibleModelFitting',
             thread_count=128,
             tile_memory_limit=16384,
             tile_size=2000,
             verbose=True)
    # Save to S3
    utils.save_s3(output_catname, 'aurelien_sepp', 'sepp/GDS/tiles')

def main_run(img_dir, sepp_dir, psf_dir, config_dir):
    # Python configuration file for SE++
    python_config=f'{config_dir}/sepp-config.py'
    # Detection image
    detect_img=glob.glob(f'{img_dir}/*ir*sci*.fits')[0]
    # Filter to use on the detection image
    filt=f'{config_dir}/gauss_1.5_3x3.conv'
    # Output
    name = ".".join(detect_img.split("/")[-1].split(".")[:-1])
    output_catname=f'{sepp_dir}/{name}_sepp_cat.fits'
    checkimg_path=f'{sepp_dir}/checkimages'
    os.makedirs(checkimg_path, exist_ok=True)
    # Run SE++
    run_sepp(detect_img=detect_img,
             generic_img=f"{img_dir}/*clear*sci*.fits",
             generic_wht=f"{img_dir}/*clear*wht*.fits",
             generic_psf=f"{psf_dir}/*star_psf.psf",
             output_catname=output_catname,
             checkimg_path=checkimg_path,
             python_config=python_config,
             filt=filt,
             filter_list = None,
             fit_case='sersic_rg4',
             properties='PixelCentroid,WorldCentroid,SourceIDs,GroupInfo,GroupStamp,SourceFlags,NDetectedPixels,NCorePixel,AperturePhotometry,AutoPhotometry,FluxRadius,SNRRatio,ShapeParameters,FlexibleModelFitting',
             thread_count=128,
             tile_memory_limit=16384,
             tile_size=2000,
             verbose=True)

def main():
    tile = sys.argv[1]
    img_dir = sys.argv[2] if len(sys.argv)>2 else "/home/ec2-user/DAWN/DJA-SEpp/image/GDS/cutout/tiles"
    sepp_dir = sys.argv[3] if len(sys.argv)>3 else "/home/ec2-user/DAWN/DJA-SEpp/sepp/GDS/tiles"
    psf_dir = sys.argv[4] if len(sys.argv)>4 else "/home/ec2-user/DAWN/DJA-SEpp/psfex/GDS"
    config_dir = sys.argv[5] if len(sys.argv)>5 else "/home/ec2-user/DAWN/DJA-SEpp/config"
    main_tile(tile, img_dir, sepp_dir, psf_dir, config_dir)

if __name__=='__main__':
    main()