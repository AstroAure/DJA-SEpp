#!/usr/bin/env bash

FIELD='ceers-full-grizli-v7.2'
BASE='/home/aurelien/DAWN/DJA-SEpp'
BUCKET='aurelien-sepp'

# aws s3 sync s3://aurelien-sepp/$FIELD/psfex $BASE/fields/$FIELD/psfex --exclude="*" --include="*.psf" --dryrun
python3 download_psf.py $FIELD $BASE $BUCKET