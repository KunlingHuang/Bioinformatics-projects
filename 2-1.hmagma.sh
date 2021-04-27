#!bin/bash
while getopts "c:g:" option
do
    case "${option}"
    in
	c) chr=${OPTARG};;
    g) gene=${OPTARG};;
    esac
done

for region in FBD FBP
do
    annot=..//data/H-MAGMA/input/annot_kl/${region}/chr${chr}/${gene}.prenatal.${region}.res40000.annot.txt
    ref=/z/Comp/lu_group/Members/khuang82/Hi-C_GWAS/data/H-MAGMA/input/1KG
    gwas=/z/Comp/lu_group/Members/khuang82/Hi-C_GWAS/results/GWAS_sumstat_H-MAGMA_compatible
    
    if [[ -f "$annot" ]]
    then
        /z/Comp/lu_group/Software/MAGMA/magma \
        --bfile $ref/g1000_eur \
        --pval $gwas/Grove_ASD_2019.1KG_rsid_update.txt use=SNP,P ncol=N \
        --gene-annot $annot \
        --out ../results/hmagma/${region}/chr${chr}/kl.${gene}.res40000.txt

    fi
done
