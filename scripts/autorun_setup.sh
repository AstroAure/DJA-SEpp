#!/bin/bash

FIELD=$1
BUCKET=$2
DECOMPRESS=${3:-true}
PSF=${4:-true}
TILE_SIZE=${5:-5}
TILE_OVERLAP=${6:-0.5}
INSTANCE_TYPE=${7:-'c6id.4xlarge'}
TEMPLATE=${8:-'lt-0f2b50c4559cd895e'}

EC2_PATH='/home/ec2-user/miniconda3/envs/dawn-310/bin:/home/ec2-user/miniconda3/condabin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/home/ec2-user:/home/ec2-user/.local/bin:/home/ec2-user/bin'

SCREEN_SCRIPT="sudo -u ec2-user env "PATH=$EC2_PATH" ./setup_image_psf.sh $FIELD $BUCKET /FlashStorage $DECOMPRESS -1 $PSF
sudo -u ec2-user env "PATH=$EC2_PATH" python3 tile.py $FIELD $TILE_SIZE $TILE_OVERLAP /FlashStorage true none $BUCKET"

SCRIPT="#!/bin/bash
cd /home/ec2-user
sudo -u ec2-user git clone --branch package https://github.com/AstroAure/DJA-SEpp.git
cd /home/ec2-user/DJA-SEpp/scripts
sudo -u ec2-user screen -S SEpp -dm bash -c '$SCREEN_SCRIPT; sudo shutdown now -h'
"

echo "$SCRIPT"

# Launch instance and run code
aws ec2 run-instances  --launch-template LaunchTemplateId=$TEMPLATE \
                       --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=SEpp-setup}]" \
                       --cpu-options "CoreCount=8" \
                       --user-data "$SCRIPT" \
                       --instance-type $INSTANCE_TYPE
