# Havemeyer Demo
This repository contains a demonstration of how a variant processing pipeline can streamline the steps necessary to take raw FASTQ data to annotated VCFs. This tutorial is made up of three primary steps:
  1. FASTQ -> gVCF
  2. gVCF -> VCF
  3. VCF -> annotated VCF

## Set Up
*Conda environment*
First thing we need as a conda environment that holds a few necessary programs that will run portions of the pipeline. 
1st you'll have to download mamba if you don't already use conda/mamba

download Mambaforge installer (assuming Unix-like platform)
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
```
create a Snakemake environment which includes all Snakemake dependencies in
addition to miscellaneous modules
For older versions of mamba/conda
```
mamba create \
    -c conda-forge -c bioconda \
    -n snakemake \
    python=3.10 snakemake=7.19 pyaml xlsxwriter snpeff=5.2 bcftools=1.17
```
Creating environments can take considerable time depending on the speed of internet connection and HPC speed.

*Downloading Data*
[hold] fastq file + gvcfs of other horses + minivcf
*Downloading container*
The container contains many of the necessary reference files and programs that allow us to standardize genetic data processing. Due to the size of the files contained in the container this should take around 5 minutes. 
```
wget https://s3.msi.umn.edu/wags/horse/horse_wags.sif
```
Downloading pipeline*

```
git clone https://github.com/jonahcullen/wags.git
```
 
## Step 1: FASTQ to gVCF
activate snakemake?
python code

confirm sample CSV is there and formatted correctly
FQ_DIR=/scratch.global/marlo072/DemoTest2/AH4
PROC_DIR=/scratch.global/marlo072/DemoTest2/FASTQtogVCF
singularity exec --bind $PWD /scratch.global/marlo072/DemoTest2/horse_wags.sif \
    cp /home/refgen/horse/goldenPath/goldenPath.Ec_build-3.0_wMSY.dict $PWD
snakemake -s one_wag.smk \
    --use-singularity \
    --singularity-args "-B $PWD,$REF_DIR,$POP_VCF,$FQ_DIR,$PROC_DIR" \
    -c 8
    --configfile goldenPath_config.yaml \
    --keep-going

## Step 2: gVCFs to VCF
find -name '*goldenPath.g.vcf.gz' -type f -printf '%f\t%p\n' |sed 's/\.goldenPath.g.vcf.gz//' >> joint.list
#with change to config_prep.py script
python wags/wags/config_joint.py -r goldenPath -o . -b Outputs -g ../HavemeyerDemoMine/datafiles/Equus_caballus.EquCab3.0.104_sorted.gtf.gz --sif horse_wags.sif

python wags/wags/prep_joint.py -g joint.list -s snakemake -i horse_wags.sif -p msismall,msilarge -e marlo072@umn.edu -a durwa004 -o gVCFtoVCF -c goldenPath
##with my change to the config file that added the slurm profile
## Step 3: Annotating VCF
For the purposes of this tutorial we have created a small VCF to increase processing time. This VCF can be found in the datafiles directory. Currently this step uses VEP, SnpEff, and ANNOVAR to annotate every variant present in a VCF. The pipeline then creates a final report that compares all three annotators to each other to calculate agreement of impact predictions and variant classifications. In the future this step of the pipeline will create reports summarizing variants found in genes of interest, number of high, moderate or low impact variants and based on community outputs this step will use an agreed upon combination of tools for increased confidence in predictions. 

To activate the annotation pipeline on the sample VCF please run the following code

```
snakemake -s BandComp.smk \
  --configfile config.yaml \
  --use-singularity \
  --rerun-incomplete \
  -c 6 \
  --profile profile.go_compare
```
This should take approximately 10 minutes. 
Final outputs include a VCF annotated by all three annotators, summary files created by SnpEff and VEP for each chromosome, final comparison results of agreement between all three programs
 
