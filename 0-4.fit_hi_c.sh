#!/bin/bash
fithicPath="/z/Comp/lu_group/Software/fithic/fithic"
# Only processed 40k
res=40000
resolution=40k
i=$1

for name in FBD FBP
do
    prjPath="../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fit_hi-c/${name}"
    mkdir -p $prjPath
    intx="../data/Hi-C/fetal_brain_Won_2016/bed/$resolution/$name/GSE77565_${name}_chr${i}.bed"
    frag="../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fragment/${res}/chr${i}.fragments.${name}.gz"
    fileName="${name}_chr${i}"
    gzip -c ${intx} > ${intx}.gz
    python3 $fithicPath/fithic.py -f ${frag} -i ${intx}.gz -o ${prjPath} -b 200 -L 20000 -U 10000000 -r ${res} -p 2 -l ${fileName}_FitHiC
    # python3 $fithicPath/fithic.py -f ../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fragment/${res}/${name}/chr${i}.fragments.${name}.gz -i ${intx}.gz -o ${prjPath} -b 200 -L 20000 -U 10000000 -r ${res} -p 2 -l ${fileName}_FitHiC
    # gunzip ../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/${name}/${name}_chr${i}_FitHiC.spline_pass2.res${res}.significances.txt.gz
done



