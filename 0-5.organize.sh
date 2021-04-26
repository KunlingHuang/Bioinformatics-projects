#!/bin/bash
for name in FBD FBP
do
    for i in {1..22}
    do
        # 10k
        res=10000
        cp /z/Comp/lu_group/Members/khuang82/Hi-C_GWAS/data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fit_hi-c/${name}/${name}_chr${i}_FitHiC.spline_pass2.res${res}.significances.txt ../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fit_hi-c_unzipped/${name}/

        # 40k
        res=40000
        # gunzip ../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fit_hi-c/${name}/${name}_chr${i}_FitHiC.spline_pass2.res${res}.significances.txt.gz
        cp ../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fit_hi-c/${name}/${name}_chr${i}_FitHiC.spline_pass2.res${res}.significances.txt ../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fit_hi-c_unzipped/${name}/
    done
done

# move the files
for name in FBD FBP
do
    for i in {1..22}
    do
        mv ../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fit_hi-c/${name}/${name}_chr${i}_FitHiC.spline_pass2.res${res}.significances.txt ../data/Hi-C/fetal_brain_Won_2016/bed_fit_hi-c/fit_hi-c_unzipped/FBD/
    done
done
