import os
import glob
import re
import sys

def run_sepp(python_config : str,
             args : dict,
             output_catname : str,
             properties : str,
             detect_img : str,
             filt : str,
             checkimg_path : str,
             thread_count : int,
             tile_memory_limit : int,
             tile_size : int):
    """
    Run SourceXtractor++ on a given set of images
    
    python_config : filename of the Python configuration file
    args : dictionary of arguments to give to the configuration file.
        It must contain at least :
            - fit_case : the model for fitting in SE++
            - list_of_IMG_names : list of images to use
            - list_of_WHT_names : list of weight images to use
            - list_of_PSF_names : list of PSF files to use
    output_catname : filename of the output catalog
    properties : SE++ properties to compute and export
    detect_img : filename of the detection image to use
    filt : convolution filter to use for the detection
    checkimg_path : path to the folder for the model and residual images
    thread_count : number of threads to use for SE++
    tile_memory_limit : maximum RAM allocated to SE++ (in MiB)
    tile_size : size in px of the tiles created by SE++
    """
    os.makedirs(checkimg_path, exist_ok=True)
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['MKL_DYNAMIC'] = 'FALSE'
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['OMP_DYNAMIC'] = 'FALSE'
    os.system(f"sourcextractor++ \
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
            --grouping-hard-limit 0 \
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
            ")

def main_tile(tile):
    # Python configuration file for SE++
    python_config='/home/ec2-user/DAWN/DJA-SEpp/config/sepp-config.py'
    # Field to run SE++ on (for file organization)
    field = 'GDS'
    # Detection image
    detect_img=glob.glob(f'/FlashStorage/image/{field}/tiles/*ir*sci*-{tile}.fits')[0]
    # Filter to use on the detection image
    filt='/home/ec2-user/DAWN/DJA-SEpp/config/gauss_1.5_3x3.conv'
    # Band filters to use for SE++ (here, automatically detected from the available images)
    filter_list = list(set([re.search('(f\d+\w+)', filename).group(1) for filename in glob.glob(f"/FlashStorage/image/{field}/tiles/*clear*sci*-{tile}.fits")]))
    filter_list.sort()
    print(f"Filters : {filter_list}")
    # Arguments to give to SE++ :
    args = {'fit_case' : 'sersic_rg4',
            'list_of_IMG_names' : [glob.glob(f"/FlashStorage/image/{field}/tiles/*{filter}*sci*-{tile}.fits")[0] for filter in filter_list],
            'list_of_WHT_names' : [glob.glob(f"/FlashStorage/image/{field}/tiles/*{filter}*wht*-{tile}.fits")[0] for filter in filter_list],
            'list_of_PSF_names' : [glob.glob(f"/home/ec2-user/DAWN/DJA-SEpp/psfex/{field}/*{filter}*star_psf.psf")[0] for filter in filter_list],
            }
    print(f"Images  : {args['list_of_IMG_names']}")
    print(f"Weights : {args['list_of_WHT_names']}")
    print(f"PSFs    : {args['list_of_PSF_names']}")
    # Output
    output_dir = "/FlashStorage/sepp/tiles"
    properties='PixelCentroid,WorldCentroid,SourceIDs,GroupInfo,GroupStamp,SourceFlags,NDetectedPixels,NCorePixel,AperturePhotometry,AutoPhotometry,FluxRadius,SNRRatio,ShapeParameters,FlexibleModelFitting'
    name = ".".join(detect_img.split("/")[-1].split(".")[:-1])
    output_catname=f'{output_dir}/{field}/{name}_sepp_cat.fits'
    checkimg_path=f'{output_dir}/{field}/checkimages'
    os.makedirs(checkimg_path, exist_ok=True)
    # Run SE++
    run_sepp(python_config=python_config, 
             args=args, 
             output_catname=output_catname,
             properties=properties,
             detect_img=detect_img,
             filt=filt,
             checkimg_path=checkimg_path,
             thread_count=128,
             tile_memory_limit=16384,
             tile_size=2000)

def main_cutout():
    ### Python configuration file for SE++
    python_config='/home/ec2-user/DAWN/DJA-SEpp/config/sepp-config.py'
    ### Field to run SE++ on (for file organization)
    field = 'GDS'
    ### Detection image
    detect_img=glob.glob(f'/FlashStorage/image/{field}/cutout/*ir*sci_cutout.fits')[0]
    ### Filter to use on the detection image
    filt='/home/ec2-user/DAWN/DJA-SEpp/config/gauss_1.5_3x3.conv'
    ### Band filters to use for SE++ (here, automatically detected from the available images)
    filter_list = list(set([re.search('(f\d+\w+)', filename).group(1) for filename in glob.glob(f"/FlashStorage/image/{field}/cutout/*clear*sci*.fits")]))
    filter_list.sort()
    print(f"Filters : {filter_list}")
    ### Arguments to give to SE++ :
    ###       - fit_case : the model for fitting in SE++ (keep 'sersic_rg4' for detection mode)
    ###       - list_of_IMG_names : list of images to use (here, automatically detected using 'filter_list')
    ###       - list_of_WHT_names : list of weight images to use (here, automatically detected using 'filter_list')
    ###       - list_of_PSF_names : list of PSF files to use (here, automatically detected using 'filter_list')
    args = {'fit_case' : 'sersic_rg4',
            'list_of_IMG_names' : [glob.glob(f"/FlashStorage/image/{field}/cutout/*{filter}*sci_cutout.fits")[0] for filter in filter_list],
            'list_of_WHT_names' : [glob.glob(f"/FlashStorage/image/{field}/cutout/*{filter}*wht_cutout.fits")[0] for filter in filter_list],
            'list_of_PSF_names' : [glob.glob(f"/home/ec2-user/DAWN/DJA-SEpp/psfex/{field}/*{filter}*star_psf.psf")[0] for filter in filter_list],
            }
    print(f"Images  : {args['list_of_IMG_names']}")
    print(f"Weights : {args['list_of_WHT_names']}")
    print(f"PSFs    : {args['list_of_PSF_names']}")
    ### Output
    output_dir = "/FlashStorage/sepp"
    properties='PixelCentroid,WorldCentroid,SourceIDs,GroupInfo,GroupStamp,SourceFlags,NDetectedPixels,NCorePixel,AperturePhotometry,AutoPhotometry,FluxRadius,SNRRatio,ShapeParameters,FlexibleModelFitting'
    name = ".".join(detect_img.split("/")[-1].split(".")[:-1])
    output_catname=f'{output_dir}/{field}/{name}_sepp_cat.fits'
    checkimg_path=f'{output_dir}/{field}/checkimages'
    ### Run SE++
    run_sepp(python_config=python_config, 
             args=args, 
             output_catname=output_catname,
             properties=properties,
             detect_img=detect_img,
             filt=filt,
             checkimg_path=checkimg_path,
             thread_count=128,
             tile_memory_limit=16384,
             tile_size=2000)

def main():
    main_tile(sys.argv[1])

if __name__=='__main__':
    main()