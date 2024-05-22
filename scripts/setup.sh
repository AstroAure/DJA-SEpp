#!/bin/bash

mount-storage /FlashStorage/
python3 -m pip install --upgrade dja_sepp
mkdir /FlashStorage/DJA-SEpp
cd /FlashStorage/DJA-SEpp
wget https://github.com/AstroAure/DJA-SEpp/archive/main.zip
unzip main.zip 'DJA-SEpp-main/config/*'
rm main.zip
mv DJA-SEpp-main/config/ config
rm -rf DJA-SEpp-main/