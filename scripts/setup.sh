#!/bin/bash

HOME=$1

if [[ "$HOME" == *"FlashStorage"* ]]
then
    mount-storage /FlashStorage/
fi
python3 -m pip install --upgrade dja_sepp
mkdir $HOME
cd $HOME
wget https://github.com/AstroAure/DJA-SEpp/archive/main.zip
unzip main.zip 'DJA-SEpp-main/config/*'
rm main.zip
mv DJA-SEpp-main/config/ config
rm -rf DJA-SEpp-main/