import re
import os
import gzip
import shutil
import boto3

def find_files(bucket:str, path:str, regex:str) -> list[str]:
    """
    Find all files in a bucket matching a given regular expression.

    bucket : Bucket to search files in
    path : Path to search files in
    regex : Regular expression to match

    Returns : List of file names (without path)
    """
    files = []
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket)
    for obj in bucket.objects.filter(Prefix=f"{path}/"):
        if re.fullmatch(regex, obj.key[len(path)+1:]):
            files.append(obj.key[len(path)+1:])
    return files

def decompress_save(file_name:str, 
                    in_bucket:str, in_path:str,
                    out_path:str,
                    verbose=False):
    """
    Downloads a compressed .gz file from an S3 bucket, 
    and saves a decompressed version to your computer.

    file_name : Name of the compressed file
    in_bucket : Bucket where the compressed file is
    in_path : Path in the in_bucket to the compressed file
    out_path : Path in your computer to save the decompressed file
    """
    s3 = boto3.client('s3')
    os.makedirs(out_path, exist_ok=True)
    # Download compressed images from S3
    if verbose: print(f"Downloading : https://s3.amazonaws.com/{in_bucket}/{in_path}/{file_name}")
    s3.download_file(in_bucket, f"{in_path}/{file_name}", f"{out_path}/{file_name}")
    # Decompress images
    if verbose: print(f"Decompressing : {out_path}/{file_name}")
    with gzip.open(f"{out_path}/{file_name}", 'rb') as image_gz:
        with open(f"{out_path}/{file_name[:-3]}", 'wb') as image_fits:
            shutil.copyfileobj(image_gz, image_fits)
    # Delete compressed images
    if verbose: print(f"Removing : {out_path}/{file_name}")
    os.remove(f"{out_path}/{file_name}")

def decompress_save_to_S3(file_name:str, 
                          in_bucket:str, in_path:str,
                          out_bucket:str, out_path:str,
                          temp_folder:str,
                          deleting_file=True,
                          verbose=False):
    """
    Downloads a compressed .gz file from an S3 bucket, 
    and saves a decompressed version in another S3 bucket.

    file_name : Name of the compressed file
    in_bucket : Bucket where the compressed file is
    in_path : Path in the in_bucket to the compressed file
    out_bucket : Bucket to save the decompressed file to
    out_path : Path in the out_bucket to save the decompressed file
    temp_folder : Temporary folder to manipulate files
    """
    s3 = boto3.client('s3')
    # Download and decompress images from S3
    decompress_save(file_name=file_name, 
                    in_bucket=in_bucket, in_path=in_path,
                    out_path=temp_folder, verbose=verbose)
    #Save decompressed images to S3
    if verbose: print(f"Uploading : https://s3.amazonaws.com/{out_bucket}/{out_path}/{file_name[:-3]}")
    s3.upload_file(f"{temp_folder}/{file_name[:-3]}", out_bucket, f"{out_path}/{file_name[:-3]}")
    # Delete decompressed images
    if deleting_file:
        if verbose: print(f"Removing : {temp_folder}/{file_name[:-3]}")
        os.remove(f"{temp_folder}/{file_name[:-3]}")


def save_s3(file, bucket, path):
    """
    Save a file to S3 bucket (with the same name)
    
    file : filename of the file to save
    bucket : S3 bucket name
    path : path in the S3 bucket to save the file"""
    s3 = boto3.client('s3')
    name = file.split("/")[-1]
    s3.upload_file(file, bucket, f"{path}/{name}")