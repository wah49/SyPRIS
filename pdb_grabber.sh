#!/bin/bash

for thing in `cat $1`
    do
       echo "$thing"".pdb.gz"
       wget http://www.rcsb.org/pdb/files/$thing.pdb.gz
       sleep 0.25
       gunzip $thing.pdb.gz
   done
