#!/bin/bash
for res in 10 40
do
    for i in {1..22}
    do
        # nohup python3 convertHDF5toSparseM.py --library FBD --chr $i --res $res &
        nohup python3 convertHDF5toSparseM.py --library FBP --chr $i --res $res &  
    done  
done

