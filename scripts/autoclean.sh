#!/bin/bash

FIELD=$1
BUCKET=${2:-'aurelien-sepp'}

aws s3 rm s3://$BUCKET/$FIELD/image/ --recursive --exclude="*" --include="*.fits"
aws s3 rm s3://$BUCKET/$FIELD/sepp/ --recursive --exclude="*tile-full*"