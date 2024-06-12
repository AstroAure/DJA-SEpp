#!/usr/bin/env bash

FIELD=$1
BUCKET=$2
BASE=$3
DECOMPRESS=${4:-false}
TILE=$5
PSF=${6:-true}
BUCKT_FOLDER=${7:-'tiles'}

./setup.sh $BASE

# Decompress or download full images
if $DECOMPRESS
then
    python3 decompress.py $FIELD $BASE 1 $BUCKET
else
    # Download full or tile
    if [ -z "$TILE" ]
    then
        python3 download_full.py $FIELD $BASE $BUCKET
    else
        python3 download_tile.py $FIELD $TILE $BASE $BUCKET $BUCKT_FOLDER
    fi
fi

# Calculate PSF or download PSF
if $PSF
then
    python3 psf.py $FIELD $BASE $BASE/config $BUCKET
else
    python3 download_psf.py $FIELD $BASE $BUCKET
fi
