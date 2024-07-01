#!/bin/bash

FIELD=$1
BUCKET=$2
TILE=$3
FIT_CASE=${4:-'sersic_rg4'}
THREAD_COUNT=${5:-32}
BUCKET_FOLDER=${6:-'tiles'}
INSTANCE_TYPE=${7:-'c6a.4xlarge'}
TEMPLATE=${8:-'lt-0f2b50c4559cd895e'}

EC2_PATH='/home/ec2-user/miniconda3/envs/dawn-310/bin:/home/ec2-user/miniconda3/condabin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/home/ec2-user:/home/ec2-user/.local/bin:/home/ec2-user/bin'

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
sudo -u ec2-user screen -S SEpp -dm bash -c 'sudo -u ec2-user env "PATH=$EC2_PATH" ./run_tile.sh $FIELD $BUCKET /home/ec2-user/RUN $TILE $FIT_CASE $THREAD_COUNT false false $BUCKET_FOLDER; sudo shutdown now -h'
"

echo "$SCRIPT"

# Launch instance and run code
aws ec2 run-instances  --launch-template LaunchTemplateId=$TEMPLATE \
                       --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=SEpp-auto-tile-$TILE}]" \
                       --user-data "$SCRIPT" \
                       --instance-type $INSTANCE_TYPE
