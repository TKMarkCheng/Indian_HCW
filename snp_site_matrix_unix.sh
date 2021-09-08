# target: receive multiple alignment.fasta, output = matrix of snp between each individual sequence
#!/bin/bash
# activate conda env
# source home/martha/miniconda3/etc/profile.d/conda.sh conda activate gupta_lab

###pipeline variable
indir=.
outdir=./output
refseq=./reference/MN908947.3.fasta

### end variables
mkdir ${outdir}

### aligned* searches for files starting with aligned---
for sample in *.fasta; do
    echo "perform function on $sample"

    ###snp-dists for future validation and cross reference
    snp-dists -c ${indir}/${sample} > ${outdir}/${sample%%.*}_snp_dist.csv

    ### snp-sites to generate the first .vcf
    snp-sites -v ${indir}/${sample} > ${outdir}/${sample%%.*}.vcf
    ### adjust the .vcf file 1) change all CHROM from 1 to MN908947.3, 2) split multiallelic variants
    echo "1 MN908947.3" >> ./chr_name_conv.txt

    bcftools annotate --rename-chrs chr_name_conv.txt ${outdir}/${sample%%.*}.vcf >> ${outdir}/${sample%%.*}_rename.vcf 

    python Split.py ${outdir}/${sample%%.*}_rename.vcf ${outdir}/${sample%%.*}_rename_split.vcf

    ### snpEff to convert to .ann.vcf
    java -Xmx8g -jar ~/snpEff/snpEff.jar -v MN908947.3 ${outdir}/${sample%%.*}_rename_split.vcf > ${outdir}/${sample%%.*}_rename_split.ann.vcf

    ### python script: python matrix.py [input.vcf] [output_mutation_matrix.xlsx]
    python matrix.py ${outdir}/${sample%%.*}_rename_split.ann.vcf ${outdir}/${sample%%.*}_matrix.xlsx

    echo "done ${sample}"
done