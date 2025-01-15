#!/bin/bash

x=1
make cnc
while [ $x -le 100 ]
do
  echo ${x}
  CILK_NWORKERS=$1 ./cn_c $2 $3 
  x=$(( $x + 1 ))
  sleep 1
done