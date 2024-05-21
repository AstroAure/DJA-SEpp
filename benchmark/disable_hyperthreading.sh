#!/bin/bash
val=${1:-0}
for cpunum in $(cat /sys/devices/system/cpu/cpu*/topology/thread_siblings_list | cut -s -d, -f2- | tr ',' '\n' | sort -un)
do
	echo $val > /sys/devices/system/cpu/cpu$cpunum/online
done