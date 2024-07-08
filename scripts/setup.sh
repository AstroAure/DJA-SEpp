#!/bin/bash

BASE=$1

if [[ "$BASE" == *"FlashStorage"* ]]
then
    mount-storage /FlashStorage
fi
python3 -m pip install --upgrade dja_sepp
mkdir $BASE
cd $BASE
wget https://github.com/AstroAure/DJA-SEpp/archive/main.zip
unzip main.zip 'DJA-SEpp-main/config/*'
rm main.zip
mv DJA-SEpp-main/config/ config
rm -rf DJA-SEpp-main/