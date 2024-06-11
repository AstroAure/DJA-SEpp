#!/usr/bin/env bash

FIELD=$1
BUCKET=$2
HOME=$3
DECOMPRESS=${4:-false}
TILE=$5
PSF=${6:-true}

./setup.sh $HOME

# Decompress or download full images
if $DECOMPRESS
then
    python3 decompress.py $FIELD $HOME 1 $BUCKET
else
    # Download full or tile
    if [ -z "$TILE" ]
    then
        python3 download_full.py $FIELD $HOME $BUCKET
    else
        python3 download_tile.py $FIELD $TILE $HOME $BUCKET
    fi
fi

# Calculate PSF or download PSF
if $PSF
then
    python3 psf.py $FIELD $HOME $HOME/config $BUCKET
else
    python3 download_psf.py $FIELD $HOME $BUCKET
fi
