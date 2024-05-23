#!/bin/bash

TILE=$1

if [ -z "$TILE" ]
then
    echo 'False'
else
    echo $TILE
fi