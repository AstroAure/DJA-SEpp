import sys
import os
import boto3
import dja_sepp

# python3 download_tile.py ceers-full-grizli-v7.2 0 /FlashStorage/DJA-SEpp aurelien-sepp

def main():
    field = sys.argv[1]
    tile = sys.argv[2]
    home = sys.argv[3] if len(sys.argv)>3 else "/tmp"
    bucket = sys.argv[4] if len(sys.argv)>4 else 'aurelien-sepp'
    bucket_folder = sys.argv[5] if len(sys.argv)>5 else 'tiles'

    files = dja_sepp.s3.find_files(bucket=bucket, 
                                   path=f'{field}/image/{bucket_folder}', 
                                   regex=f"[^/]+((f\d+(w|m)-.*clear_drc)|ir).+(sci|wht).+(tile-{tile})\.fits")
    folder = f"{home}/fields/{field}/image/tiles"
    os.makedirs(folder, exist_ok=True)
    s3 = boto3.client('s3')
    for file in files:
        print(file)
        s3.download_file(bucket, f"{field}/image/{bucket_folder}/{file}", f"{folder}/{file}")

if __name__=='__main__':
    main()