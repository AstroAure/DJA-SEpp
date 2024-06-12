#!/usr/bin/env bash
# WIP

FIELD=$1
BUCKET=$2
TILE=$3
FIT_CASE=${4:-'sersic_rg4'}
THREAD_COUNT=${5:-32}
BUCKET_FOLDER=${6:-'tiles'}
TEMPLATE=${7:-'lt-0f2b50c4559cd895e'}

# SCRIPT="#!/bin/bash
# cd /home/ec2-user
# echo 'Field : $FIELD' > /home/ec2-user/LOG.log
# echo 'Tile  : $TILE' >> /home/ec2-user/LOG.log
# echo 'Fit   : $FIT_CASE' >> /home/ec2-user/LOG.log
# echo '' >> /home/ec2-user/LOG.log
# echo 'Cloning git' >> /home/ec2-user/LOG.log
# git clone --branch package https://github.com/AstroAure/DJA-SEpp.git
# screen -S SEpp
# echo 'Running SE++' >> /home/ec2-user/LOG.log
# cd /home/ec2-user/DJA-SEpp/scripts
# chmod +x setup.sh
# chmod +x setup_image_psf.sh
# chmod +x run_tile.sh
# ./run_tile.sh $FIELD $BUCKET /home/ec2-user/RUN $TILE $FIT_CASE $THREAD_COUNT false false $BUCKET_FOLDER 1>> /home/ec2-user/LOG.log 2>> /home/ec2-user/LOG.log"
#shutdown now -h"

SCRIPT="#!/bin/bash
sudo -u ec2-user -i
cd /home/ec2-user
echo 'Field : $FIELD' > /home/ec2-user/LOG.log
echo 'Tile  : $TILE' >> /home/ec2-user/LOG.log
echo 'Fit   : $FIT_CASE' >> /home/ec2-user/LOG.log
echo '' >> /home/ec2-user/LOG.log
echo 'Cloning git' >> /home/ec2-user/LOG.log
git clone --branch package https://github.com/AstroAure/DJA-SEpp.git
echo 'Setting up' >> /home/ec2-user/LOG.log
cd DJA-SEpp/scripts
chmod +x setup.sh
./setup.sh /home/ec2-user/RUN
echo 'Downloading PSF files' >> /home/ec2-user/LOG.log
python3 -u /home/ec2-user/DJA-SEpp/scripts/download_psf.py $FIELD /home/ec2-user/RUN $BUCKET 1>> /home/ec2-user/LOG.log 2>> /home/ec2-user/LOG.log
echo 'Downloading tile images' >> /home/ec2-user/LOG.log
python3 -u /home/ec2-user/DJA-SEpp/scripts/download_tile.py $FIELD $TILE /home/ec2-user/RUN $BUCKET big-tiles 1>> /home/ec2-user/LOG.log 2>> /home/ec2-user/LOG.log"

echo "$SCRIPT"

# Launch instance and run code
aws ec2 run-instances  --launch-template LaunchTemplateId=$TEMPLATE \
                       --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=SEpp-auto-tile-$TILE}]" \
                       --user-data "$SCRIPT"
