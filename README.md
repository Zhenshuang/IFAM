# IFAM
**I**ntegrating **F**unctional **A**nnotation information by genomic BLUP model with **M**ultiple random effects

## About
IFAM adopts a multiple random effects model to  integrate massive types of genomic functional annotation information such as enhancer, promoter, and transcription factor binding site et al. to improve the accuracy of genomic prediction for complex traits. In IFAM, functional annotation markers with similar contributions to phenotype variance are automatically merged to construct random effects, the association between markers and trait is utilized as another random effect. <br>
Details of IFAM could be found in our [IFAM manuscript](https:****).

## Tutorial for IFAM
In this softare, we use `**.R` scripts to make the usage of IFAM. 

### Input files and formats
* Genotype file: IFAM only accept the genotype in PLINK binary format, e.g. demo.fam, demo.bim and demo.bed, please see more details about these files at PLINK user manual. Users can convert any other format of genotype (e.g. VCF, HapMap, PED/MAP) to binary format by PLINK2.

* Phenotype file (e.g. `demo_data/pheno.txt`): this file includes the phenotypic records, the environmental covariates, fixed and random effects. The first column must be the individual id, the second column is phenotypic records (optional), header should be included in the file.

* Annotation files: a list of the annotation files, the file name must be "annotation_names.txt" eg. "enhancer.txt". Each annotation file has only one column representing the position of annotationed SNPs. If the annotation file is "*.bed" format, users can filtered all annotationed SNPs using "**.R" script.

* Annotation files: a list of the annotation files, the file name must be "annotation_names.txt" eg. "enhancer.txt". Each annotation file has only one column representing the position of annotationed SNPs, note that the genome version of the annotation file and genome file should be the same. If the annotation file is "*.bed" format, users can filtered all annotationed SNPs using bedtools, PLINK, vcftools, and other softwares.


For other parapeters, please use the commond "--help"

### Tutorial for Summary Analysis for Annotations (SAA)
```bash
# Set parameters
SAA=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/SAA.R
map=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/geno/pig.bim
map_format=bim  # or map3
anno=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/OCR.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/NFR.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/footprint.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/enhancer.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/any_motif.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/active_promoter.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/active_promoter_narrowPeak.txt
outPath=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/test/
output_prefix=test1

# Run SAA.R
Rscript ${SAA} --map ${map} --map_format ${map_format} --anno ${anno} --outPath ${outPath}\
        --output_prefix ${output_prefix}
````

### Tutorial for IFAM 
Please install [HIBLUP](https://www.hiblup.com/tutorials#running-hiblup) software in advance
```bash
# Set parameters
IFAM=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/IFAM.R
bfile=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/geno/pig
pheno=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/phe/BF.cv.phe.txt
anno=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/OCR.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/NFR.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/footprint.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/enhancer.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/any_motif.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/active_promoter.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/active_promoter_narrowPeak.txt
anno_spec=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/signals.txt
anno_GRM=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/OCR.GA,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/NFR.GA,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/footprint.GA,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/enhancer.GA,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/any_motif.GA,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/active_promoter.GA,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/active_promoter_narrowPeak.GA,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/signals.GA
weight=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/BF_cv_1_weight.txt
#VCfile=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/BF.vars
outPath=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/test/
output_prefix=test1
pheno_pos=3
randomMax=5
thread=2
VCmethod=AI
tmp_files=FALSE

# Run IFAM
Rscript ${IFAM} --bfile ${bfile} --pheno ${pheno} --anno ${anno} --anno_spec ${anno_spec}\
        --anno_GRM ${anno_GRM} --pheno_pos ${pheno_pos} --weight ${weight} --VCfile ${VCfile} --randomMax ${randomMax}\
        --VCmethod ${VCmethod} --thread ${thread} --tmp_files ${tmp_files} --outPath ${outPath}\
        --output_prefix ${output_prefix}
````

### Tutorial for Evaluating the Annotations (EA)
Please install [HIBLUP](https://www.hiblup.com/tutorials#running-hiblup) and [PLINK] softwares in advance
```bash
# Set parameters
EA=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/EA.R
bfile=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/geno/pig
pheno=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/phe/BF.cv.phe.txt
anno=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/OCR.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/NFR.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/footprint.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/enhancer.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/any_motif.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/active_promoter.txt,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/annotations/active_promoter_narrowPeak.txt
anno_GRM=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/OCR.GA,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/NFR.GA,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/footprint.GA,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/enhancer.GA,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/any_motif.GA,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/active_promoter.GA,/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/GRMs/active_promoter_narrowPeak.GA
Pruning=FALSE
indep_pairwise=1000,100,0.2
plink=plink
outPath=/share/home/yzbsl_tangzs/IFAM/IFAM_demodata/test/
output_prefix=test1
pheno_pos=3
thread=2
VCmethod=AI
tmp_files=FALSE

# Run EA
Rscript ${EA} --bfile ${bfile} --pheno ${pheno} --anno ${anno} --anno_GRM ${anno_GRM}\
        --Pruning ${Pruning} --indep_pairwise ${indep_pairwise} --plink ${plink} --pheno_pos ${pheno_pos} \
        --VCmethod ${VCmethod} --thread ${thread} --tmp_files ${tmp_files} --outPath ${outPath}\
        --output_prefix ${output_prefix}
````
 
# Citation
For IFAM:
...........   <br>
For HIBLUP software:
Lilin Yin, Haohao Zhang, Zhenshuang Tang, Dong Yin, Yuhua Fu, Xiaohui Yuan, Xinyun Li, Xiaolei Liu, Shuhong Zhao, HIBLUP: an integration of statistical models on the BLUP framework for efficient genetic evaluation using big genomic data, Nucleic Acids Research, Volume 51, Issue 8, 8 May 2023, Pages 3501–3512, https://doi.org/10.1093/nar/gkad074.

