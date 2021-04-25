#!/bin/bash
resolution=100
mkdir -p ../data/Hi-C/fetal_brain_Won_2016/raw_hdf5/${resolution}k
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77565/suppl/GSE77565_FBD_IC-heatmap-res-${resolution}k.hdf5.gz -P ../data/Hi-C/fetal_brain_Won_2016/raw_hdf5/${resolution}k
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77565/suppl/GSE77565_FBP_IC-heatmap-res-${resolution}k.hdf5.gz -P ../data/Hi-C/fetal_brain_Won_2016/raw_hdf5/${resolution}k

mv ../data/Hi-C/fetal_brain_Won_2016/raw_hdf5/${resolution}k/GSE77565_FBD_IC-heatmap-res-${resolution}k.hdf5.gz ../data/Hi-C/fetal_brain_Won_2016/raw_hdf5/${resolution}k/GSE77565_FBD_IC-heatmap-chr-${resolution}k.hdf5.gz
mv ../data/Hi-C/fetal_brain_Won_2016/raw_hdf5/${resolution}k/GSE77565_FBP_IC-heatmap-res-${resolution}k.hdf5.gz ../data/Hi-C/fetal_brain_Won_2016/raw_hdf5/${resolution}k/GSE77565_FBP_IC-heatmap-chr-${resolution}k.hdf5.gz

for resolution in 10 40
do

    mkdir -p ../data/Hi-C/fetal_brain_Won_2016/${resolution}k
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77565/suppl/GSE77565_FBD_IC-heatmap-chr-${resolution}k.hdf5.gz -P ../data/Hi-C/fetal_brain_Won_2016/raw_hdf5/${resolution}k
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77565/suppl/GSE77565_FBP_IC-heatmap-chr-${resolution}k.hdf5.gz -P ../data/Hi-C/fetal_brain_Won_2016/raw_hdf5/${resolution}k

done

# Decompressing
for resolution in 10 40 100
do
    echo $resolution
    gunzip ../data/Hi-C/fetal_brain_Won_2016/raw_hdf5/${resolution}k/GSE77565_FBD_IC-heatmap-chr-${resolution}k.hdf5.gz
    gunzip ../data/Hi-C/fetal_brain_Won_2016/raw_hdf5/${resolution}k/GSE77565_FBP_IC-heatmap-chr-${resolution}k.hdf5.gz

done
