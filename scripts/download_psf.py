import sys
import os
import boto3
import dja_sepp

# python3 download_psf.py ceers-full-grizli-v7.2 /FlashStorage/DJA-SEpp aurelien-sepp

def main():
    field = sys.argv[1]
    home = sys.argv[2] if len(sys.argv)>2 else "/tmp"
    bucket = sys.argv[3] if len(sys.argv)>3 else 'aurelien-sepp'

    files = dja_sepp.s3.find_files(bucket=bucket, 
                                   path=f'{field}/psfex', 
                                   regex=".+(f\d+(w|m)).+psf\.psf")
    folder = f"{home}/fields/{field}/psfex"
    os.makedirs(folder, exist_ok=True)
    s3 = boto3.client('s3')
    for file in files:
        print(file)
        s3.download_file(bucket, f"{field}/psfex/{file}", f"{folder}/{file}")

if __name__=='__main__':
    main()