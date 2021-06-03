#!/bin/bash

sfincsScan
for i in */; do cp norm.namelist species $i; done
scancel -u sbul

echo "Press any key to continue"
while [ true ] ; do
read -t 3 -n 1
if [ $? = 0 ] ; then
break ;
else
echo "waiting for the keypress"
fi
done

python update_subdir.py
for i in */; do cd $i; sfincsScan; cd ..; done
