import sys
import dja_sepp

# python3 decompress.py ceers-full-grizli-v7.2 /FlashStorage/DJA-SEpp 0 aurelien-sepp

def main():
    field = sys.argv[1]
    home = sys.argv[2] if len(sys.argv)>2 else "/tmp"
    keep = sys.argv[3] if len(sys.argv)>3 else False
    bucket = sys.argv[4] if len(sys.argv)>4 else 'aurelien-sepp'
    
    files = dja_sepp.s3.find_files(bucket='grizli-v2', 
                                   path='JwstMosaics/v7', 
                                   regex=f"{field}.+((f\d+(w|m)-.*clear_drc)|ir).+(sci|wht).+")
    for i, file_name in enumerate(files):
        print(f"{i+1}/{len(files)}")
        dja_sepp.s3.decompress_save_to_S3(file_name,
                                        in_bucket='grizli-v2', in_path='JwstMosaics/v7',
                                        out_bucket=bucket, out_path=f"{field}/image",
                                        temp_folder=f"{home}/fields/{field}/image",
                                        deleting_file=not keep, verbose=True)

if __name__=='__main__':
    main()