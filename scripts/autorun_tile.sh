#!/usr/bin/env bash
# WIP

FIELD=$1
BUCKET=$2
TILE=$3
FIT_CASE=${4:-'sersic_rg4'}
THREAD_COUNT=${5:-32}
TEMPLATE=${6:-'lt-0f2b50c4559cd895e'}

# Create User Data script
# chmod +x run_tile.sh
# ./run_tile.sh $FIELD $BUCKET ~/RUN $TILE $FIT_CASE $THREAD_COUNT
# SCRIPT="#!/bin/bash
# cd /home/ec2-user
# git clone https://github.com/AstroAure/DJA-SEpp.git
# cd DJA-SEpp/scripts
# chmod +x setup.sh
# ./setup.sh /home/ec2-user/RUN
# python3 download_psf.py $FIELD /home/ec2-user/RUN $BUCKET
# python3 download_tile.py $FIELD $TILE /home/ec2-user/RUN $BUCKET
# python3 sepp.py $FIELD /home/ec2-user/RUN tiles /home/ec2-user/RUN/config $BUCKET $FIT_CASE $THREAD_COUNT $TILE"
#shutdown now -h"

SCRIPT="#!/bin/bash
cd /home/ec2-user
echo 'Field : $FIELD' > /home/ec2-user/LOG.log
echo 'Tile  : $TILE' >> /home/ec2-user/LOG.log
echo 'Fit   : $FIT_CASE' >> /home/ec2-user/LOG.log
echo '' >> /home/ec2-user/LOG.log
echo 'Cloning git' >> /home/ec2-user/LOG.log
git clone https://github.com/AstroAure/DJA-SEpp.git
echo 'Setting up' >> /home/ec2-user/LOG.log
cd DJA-SEpp/scripts
chmod +x setup.sh
./setup.sh /home/ec2-user/RUN
echo 'Downloading PSF files' >> /home/ec2-user/LOG.log
python3 -u /home/ec2-user/DJA-SEpp/scripts/download_psf.py $FIELD /home/ec2-user/RUN $BUCKET 1>> /home/ec2-user/LOG.log 2>> /home/ec2-user/LOG.log
echo 'Downloading tile images' >> /home/ec2-user/LOG.log
python3 -u /home/ec2-user/DJA-SEpp/scripts/download_tile.py $FIELD $TILE /home/ec2-user/RUN $BUCKET 1>> /home/ec2-user/LOG.log 2>> /home/ec2-user/LOG.log"

echo "$SCRIPT"
# eval "$SCRIPT"

# Launch instance and run code
aws ec2 run-instances  --launch-template LaunchTemplateId=$TEMPLATE \
                       --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=SEpp-auto-tile-$TILE}]" \
                       --user-data "$SCRIPT"
