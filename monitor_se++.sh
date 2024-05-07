#!/bin/bash

declare var=$(ps -C sourcextractor++ -o cmd,%cpu,%mem,vsz,rss,time)
echo "$var"
echo "$var" >> /home/ec2-user/DAWN/DJA-SEpp/sepp/GDS/benchmark/benchmark-$1.txt
sleep 10s
while [ ! -z "${var}" ]
do
    var=$(ps -C sourcextractor++ -o cmd,%cpu,%mem,vsz,rss,time --no-headers)
    echo "$var"
    echo "$var" >> /home/ec2-user/DAWN/DJA-SEpp/sepp/GDS/benchmark/benchmark-$1.txt
    sleep 10s
done
echo "DONE"