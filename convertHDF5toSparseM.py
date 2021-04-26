import numpy as np
import scipy as sp
from scipy import sparse
import h5py
import sys
import csv
import itertools
import argparse
import os
parser = argparse.ArgumentParser()
parser.add_argument('--library', type=str, help='The brain library FBD/FBP')
parser.add_argument('--chr', default=0, type = str, help='Chromosome')
parser.add_argument('--res', default=0, type = str, help='resolution')


args = parser.parse_args()
library = args.library
chr = args.chr
res = args.res

inFile="../data/Hi-C/fetal_brain_Won_2016/raw_hdf5/" + res + "k/GSE77565_" + library + "_IC-heatmap-chr-" + res + "k.hdf5"
## inFile="/p/keles/yezheng/volumeA/Brain/GSE77565_FBP_IC-heatmap-chr-10k.hdf5"
inhdf=h5py.File(inFile, 'r')

## get the structure of contact matrix hdf5
list(inhdf.keys()) # look how many datasets are there
for i in inhdf.keys():
   print(i)

## key for chromosome 6 interactions are '(6, 6)'
## contact matrix for chromosome 6 is
inhdf['(' + chr + ', ' + chr + ')']
#for i,j in itertools.product(range(0,len(tmp)), range(0,len(tmp))):
#    print(tmp[i, j] == tmp[j,i]) #the matrix is a symmatric matrix, just as we thought

## to write out the matrix directly
folder = "../data/Hi-C/fetal_brain_Won_2016/bed/" + res + "k/" + library
if not os.path.exists(folder):
    os.mkdir(folder)
np.savetxt("../data/Hi-C/fetal_brain_Won_2016/bed/" + res + "k/" + library + "/GSE77565_"+library+"_chr" + chr + "." + res + "k.matrix", inhdf['(' + chr + ', ' + chr + ')'])

#tmp2 = open("../../data/Hi-C/fetal_brain_Won_2016/GSE77565_FBD_chr6.sparseMatrix2", 'w')
## to write the matrix out as sparse matrix in bed file format
with open("../data/Hi-C/fetal_brain_Won_2016/bed/"+ res + "k/" +library+"/GSE77565_"+library+"_chr"+chr+".bed", 'w') as f:
    writer = csv.writer(f, delimiter = '\t')
    resolution = int(res)*1000
    for (n, m), val in np.ndenumerate(np.matrix(inhdf['(' + chr + ', ' + chr + ')'])):
        if val > 0 and m>=n:
            writer.writerow(["chr" + chr, int(n * resolution + (resolution)/2),
                            "chr" + chr, int(m * resolution + (resolution)/2), val])



