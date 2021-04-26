#!/bin/bash
#for i in {1..22}
#do
#    nohup python3 generate_fragments.py --library FBD --chr $i &
#    nohup python3 generate_fragments.py --library FBP --chr $i &
#done
mkdir -p ../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fragment/10000
mkdir -p ../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fragment/40000

fithicPath="/z/Comp/lu_group/Software/fithic/fithic"

# prjPath="/z/Comp/lu_group/Members/khuang82/Hi-C_GWAS/data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fragment"
chromSizeFile="./chr_size/hg19.chrom.sizes.txt"
# perform 40k only
for resolution in 10000 40000
do
	for i in {1..22}
	do
		awk -v CHR=chr$i '$1 == CHR {print $0}' $chromSizeFile > $chromSizeFile.chr$i ## extract chromosome 6
		python3 $fithicPath/utils/createFitHiCFragments-fixedsize.py --chrLens "${chromSizeFile}.chr$i" --resolution "$resolution" --outFile "../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fragment/${resolution}/chr$i.fragments.FBD.gz"
		python3 $fithicPath/utils/createFitHiCFragments-fixedsize.py --chrLens "${chromSizeFile}.chr$i" --resolution "$resolution" --outFile "../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fragment/${resolution}/chr$i.fragments.FBP.gz"
	done
done

