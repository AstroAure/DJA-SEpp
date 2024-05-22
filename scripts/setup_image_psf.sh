#!/bin/bash

FIELD=$1
BUCKET=$2
DECOMPRESS=${3:-false}
PSF=${4:-true}

./setup.sh

# Decompress or download full images
if $DECOMPRESS
then
    python3 decompress.py $FIELD /FlashStorage/DJA-SEpp 1 $BUCKET
else
    python3 download_full.py $FIELD /FlashStorage/DJA-SEpp $BUCKET
fi

# Calculate PSF or download PSF
if $PSF
then
    python3 psf.py $FIELD /FlashStorage/DJA-SEpp /FlashStorage/DJA-SEpp/config $BUCKET
else
    python3 download_psf.py $FIELD /FlashStorage/DJA-SEpp $BUCKET
fi
