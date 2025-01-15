#!/bin/bash

x=1
make cnc
while [ $x -le 100 ]
do
  echo ${x}
  ./cn_c $1 $2 
  x=$(( $x + 1 ))
  sleep 1
done