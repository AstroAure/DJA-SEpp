#!/bin/bash
# WIP

FIELD=$1
BUCKET=$2
TILE=$3
FIT_CASE=${4:-'sersic_rg4'}
THREAD_COUNT=${5:-32}
TEMPLATE=${6:-'lt-0f2b50c4559cd895e'}

# Create User Data script
# > script_$TILE.sh
# git clone https://github.com/AstroAure/DJA-SEpp.git
# cd DJA_SEpp/scripts
# # chmod +x run_tile.sh
# # ./run_tile.sh $FIELD $BUCKET ~/RUN $TILE $FIT_CASE $THREAD_COUNT
# ./setup.sh ~/RUN
# python3 download_psf.py $FIELD ~/RUN $BUCKET
# python3 download_tile.py $FIELD $TILE ~/RUN $BUCKET
# python3 sepp.py $FIELD ~/RUN tiles ~/RUN/config $BUCKET $FIT_CASE $THREAD_COUNT $TILE
# shutdown -h

# Launch instance and run code
aws ec2 run-instances  --launch-template LaunchTemplateId=$TEMPLATE \
                       --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=SEpp-auto-tile}]' \
                       --user-data file://script_$TILE.sh
