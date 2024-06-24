#!/bin/bash
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
# sudo -u ec2-user -i
cd /home/ec2-user
echo 'Field : $FIELD' > /home/ec2-user/LOG.log
echo 'Tile  : $TILE' >> /home/ec2-user/LOG.log
echo 'Fit   : $FIT_CASE' >> /home/ec2-user/LOG.log
echo '' >> /home/ec2-user/LOG.log
echo 'Cloning git' >> /home/ec2-user/LOG.log
sudo -u ec2-user git clone --branch package https://github.com/AstroAure/DJA-SEpp.git
echo 'Running SE++' >> /home/ec2-user/LOG.log
cd /home/ec2-user/DJA-SEpp/scripts
# sudo -u ec2-user /home/ec2-user/miniconda3/envs/dawn-310/bin/python3 -m pip install dja-sepp 1>> /home/ec2-user/LOG2.log 2>> /home/ec2-user/LOG2.log
# sudo -u ec2-user /home/ec2-user/miniconda3/envs/dawn-310/bin/python3 download_psf.py $FIELD /home/ec2-user/RUN $BUCKET 1>> /home/ec2-user/LOG2.log 2>> /home/ec2-user/LOG2.log
sudo -u ec2-user env "PATH=/home/ec2-user/miniconda3/envs/dawn-310/bin" ./setup_image_psf.sh $FIELD $BUCKET /home/ec2-user/RUN false $TILE false $BUCKET_FOLDER  1>> /home/ec2-user/LOG2.log 2>> /home/ec2-user/LOG2.log
# sudo -u ec2-user env "PATH=/home/ec2-user/miniconda3/envs/dawn-310/bin" ./run_tile.sh $FIELD $BUCKET /home/ec2-user/RUN $TILE $FIT_CASE $THREAD_COUNT false false $BUCKET_FOLDER
# sudo -u ec2-user screen -S SEpp -dm bash -c 'pip install dja-sepp; bash' "

echo "$SCRIPT"

# Launch instance and run code
aws ec2 run-instances  --launch-template LaunchTemplateId=$TEMPLATE \
                       --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=SEpp-auto-tile-$TILE}]" \
                       --user-data "$SCRIPT"
