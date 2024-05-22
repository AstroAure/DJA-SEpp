import sys
import os
import glob
import dja_sepp

# python3 sepp.py ceers-full-grizli-v7.2 /FlashStorage/DJA-SEpp cutout /FlashStorage/DJA-SEpp/config aurelien-sepp B+D 32

def main():
    field = sys.argv[1]
    home = sys.argv[2] if len(sys.argv)>2 else "/tmp"
    img_path = sys.argv[3] if len(sys.argv)>3 else ""
    config_folder = sys.argv[4] if len(sys.argv)>4 else f"{home}/config"
    bucket = sys.argv[5] if len(sys.argv)>5 else 'aurelien-sepp'
    fit_case = sys.argv[6] if len(sys.argv)>6 else 'B+D' #'sersic_rg4', 'B+D'
    thread_count = sys.argv[7] if len(sys.argv)>7 else 32
    tile = f"tile-{sys.argv[8]}" if len(sys.argv)>8 else ""

    img_dir  = f"{home}/fields/{field}/image/{img_path}"
    sepp_dir = f"{home}/fields/{field}/sepp"
    psf_dir  = f"{home}/fields/{field}/psfex"

    # Python configuration file for SE++
    python_config = f'{config_folder}/sepp-config.py'
    # Detection image
    detect_img = glob.glob(f'{img_dir}/*ir*sci*{tile}*.fits')[0]
    # Filter to use on the detection image
    filt = f'{config_folder}/gauss_1.5_3x3.conv'
    # Output
    name = ".".join(detect_img.split("/")[-1].split(".")[:-1])
    output_catname = f'{sepp_dir}/{name}_sepp_cat.fits'
    checkimg_path = f'{sepp_dir}/checkimages'
    os.makedirs(checkimg_path, exist_ok=True)

    # Run SE++
    dja_sepp.sepp.run_sepp(detect_img=detect_img,
                           generic_img=f"{img_dir}/*clear*sci*{tile}*.fits",
                           generic_wht=f"{img_dir}/*clear*wht*{tile}*.fits",
                           generic_psf=f"{psf_dir}/*star_psf.psf",
                           output_catname=output_catname,
                           checkimg_path=checkimg_path,
                           python_config=python_config,
                           filt=filt,
                           filter_list=None,
                           fit_case=fit_case,
                           properties='PixelCentroid,WorldCentroid,SourceIDs,GroupInfo,GroupStamp,SourceFlags,NDetectedPixels,NCorePixel,AperturePhotometry,AutoPhotometry,FluxRadius,SNRRatio,ShapeParameters,FlexibleModelFitting',
                           thread_count=thread_count,
                           tile_memory_limit=16384,
                           tile_size=2000,
                           verbose=True)
    # Save to S3
    dja_sepp.s3.save_s3(output_catname, bucket, f"{field}/sepp")

if __name__=='__main__':
    main()