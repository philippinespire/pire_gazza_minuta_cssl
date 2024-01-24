# Gmi Data Processing Log

Log to track progress through capture bioinformatics pipeline for the Albatross and Contemporary *Gazza minuta* samples from Hamilo Cove & Basud River. First sequencing run.

---

## 1.  1st fastp

Raw data in `/home/e1garcia/shotgun_PIRE/Gmi/raw_fq_capture` (check Gazza minuta channel on Slack).  The root outdir for all analyses will be  `/home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta`. Both on Wahab/Turing (ODU HPCs).

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta

#runFASTP_1.sbatch <indir> <outdir>
sbatch runFASTP_1.sbatch /home/e1garcia/shotgun_PIRE/Gmi/raw_fq_capture fq_fp1
```

[Report](https://htmlpreview.github.io/?https://raw.githubusercontent.com/philippinespire/pire_cssl_data_processing/main/gazza_minuta/fq_fp1/1st_fastp_report.html?token=GHSAT0AAAAAABQHSGSTOR2F6QLQTSOR6EBQYTP7Z7A).

Potential issues:  
* % duplication - low for most Albatross, higher for contemporary
  * Alb: 30-40% (some in the 40s), Contemp: 60-70% (some ~30%)
* GC content - good
  * Alb: 50%, Contemp: 50%
* passing filter - most reads passed filters for both Albatross & contemporary
  * Alb: >90%, Contemp: ~95%
* % adapter - high for Albatross (but not as high as other species), low for contemporary
  * Alb: 55-70%, Contemp: 10%
* number of reads - good for Albatross, okay for contemporary but some libraries seemed to have failed completely
  * Alb: generally much higher # (>40 mil) w/ some very high (~500 mil), Contemp: ~10-20 mil w/some VERY low (only a couple thousand)

---

## 2. Clumpify

```
#on Turing
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta

enable_lmod
#runCLUMPIFY_r1r2.sbatch <indir> <outdir> <tempdir>
sbatch runCLUMPIFY_r1r2.sbatch fq_fp1 fq_fp1_clmp /scratch-lustre/r3clark
```

Checked that all files ran with `checkCLUMPIFY.R`. All ran (no RAM issues).

---

## 3. Run fastp2

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta

#runFASTP_2.sbatch <indir> <outdir>
sbatch runFASTP_2.sbatch fq_fp1_clmp fq_fp1_clmp_fp2
```

[Report](https://htmlpreview.github.io/?https://raw.githubusercontent.com/philippinespire/pire_cssl_data_processing/main/gazza_minuta/fq_fp1_clmp_fp2/2nd_fastp_report.html?token=GHSAT0AAAAAABQHSGSSVPH5E54MK42MSIKMYTP72SQ).

Potential issues:  
* % duplication - Albatross good, contemporary okay
  * Alb: ~5-10%, Contemp: ~30%
* GC content - good
*  Alb: 45%, Contemp: 50%
* passing filter - good
  * Alb: ~99%, Contemp: ~98%
* % adapter - good
  * Alb: <2.5%, Contemp: <1%
* number of reads - went down but Albatross still good, contemporary some are very low
  * Alb: ~5-25 mil with a few >100 mil, Contemp: about half 5-15 mil & about half >1 mil (and likely in the thousands)

---

## 4. Run fastq_screen

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta

#runFQSCRN_6.bash <indir> <outdir> <number of nodes to run simultaneously>
bash ../scripts/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 20
```

[Report](https://htmlpreview.github.io/?https://raw.githubusercontent.com/philippinespire/pire_cssl_data_processing/main/gazza_minuta/fq_fp1_clmp_fp2_fqscrn/fqscrn_report.html?token=GHSAT0AAAAAABQHSGSS6DXKLMN24OV6ZHZQYTP727Q).

Potential issues: 
* job 13 failed 
  * [out file](https://github.com/philippinespire/pire_cssl_data_processing/blob/main/gazza_minuta/logs/slurm-fqscrn.405009.13.out)
  * "No reads in GmA01013..., skipping"
  * Checked file and there are definitely reads here
* Also, looks like there is more contamination (20-40%) than with other species, esp. in Albatross.

Fixed errors: re-ran file that returned the "No reads in" error.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta

#runFQSCRN_6.bash <indir> <outdir> <number of nodes to run simultaneously> <fq file pattern to process>
bash ../scripts/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 1 GmA01013*r1.fq.gz
bash ../scripts/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 1 GmA01013*r2.fq.gz
```

Cleaned-up logs.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta

mkdir logs
mv *out logs
```

---

## 5. Re-pair fastq_screen paired end files

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta

#runREPAIR.sbatch <indir> <outdir> <threads>
sbatch runREPAIR.sbatch fq_fp1_clmp_fp2_fqscrn fq_fp1_clmp_fp2_fqscrn_repaired 40
```

[Report](https://htmlpreview.github.io/?https://raw.githubusercontent.com/philippinespire/pire_cssl_data_processing/main/gazza_minuta/fq_fp1_clmp_fp2_fqscrn_repaired/multiqc_postrepair.html?token=GHSAT0AAAAAABQHSGSS2UDCAZCIIUMJ7VRWYTP73VQ).

---

## 6. Rename files for dDocent HPC and put into mapping dir

Used decode file from Sharon Magnuson & Chris Bird.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta

cp /home/e1garcia/shotgun_PIRE/Gmi/shotgun_PIRE/Gmi/raw_fq_capture/Gmi_CaptureLibraries_SequenceNameDecode.tsv .
#removed _Ex1 extensions (kept in "original.tsv" copy)

#mkNewFileNames.bash <decode file name> <fqdir>
bash ../scripts/mkNewFileNames.bash Gmi_Capture_Libraries_SequenceNameDecode.tsv fq_fp1_clmp_fp2_fqscrn_repaired > decode_newnames.txt

#make old file names
cd fq_fp1_clmp_fp2_fqscrn_repaired
ls *.fq.gz > ../decode_oldnames.txt

#triple check that the old and new files are aligned
cd ..
module load parallel
bash
parallel --no-notice --link -kj6 "echo {1}, {2}" :::: decode_oldnames.txt decode_newnames.txt > decode_translation.csv
less -S decode_translation.csv
#names do NOT match --> some contemporary libraries that were in original decode file no longer present. Adjusted decode_newnames.txt file names so that everything matched.

#rename files and move to mapping dir
mkdir mkBAM
cd fq_fp1_clmp_fp2_fqscrn_repaired
parallel --no-notice --link -kj6 "mv {1} ../mkBAM/{2}" :::: ../decode_oldnames.txt ../decode_newnames.txt

#confirm success
ls ../mkBAM
```

---

## 7. Set up mapping dir and get reference genome

Pulled latest changes from dDocentHPC repo & copied `config.5.cssl` over.

```
#if you haven't already, you first need to clone the dDocentHPC.git repo
#cd /home/r3clark/PIRE/pire_cssl_data_processing/scripts
#git clone https://github.com/cbirdlab/dDocentHPC.git

cd /home/r3clark/PIRE/pire_cssl_data_processing/scripts/dDocentHPC
git pull

cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta
cp ../scripts/dDocentHPC/configs/config.5.cssl mkBAM
```

Because there is no whole genome reference for *G. minuta*, I am using the full "rwaw" reference fasta from the RAD data used to make probes.

```
#the destination reference fasta should be named as follows: reference.<assembly type>.<unique assembly info>.fasta
#<assembly type> is `ssl` for denovo assembled shotgun library or `rad` for denovo assembled rad library
#this naming is a little messy, but it makes the ref 100% tracable back to the source
#it is critical not to use `_` in name of reference for compatibility with ddocent and freebayes

mv /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkBAM/PIRE_GazzaMinuta.C.10.10.probes4development.fasta reference.rad.RAW-10-10.fasta
```

Updated the config file with the ref genome info.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkBAM
nano config.5.cssl
```

Inserted `<assembly type>` into the `Cutoff1` variable and `<unique assembly info>` into the `Cutoff2` variable.

```
----------mkREF: Settings for de novo assembly of the reference genome--------------------------------------------
PE              Type of reads for assembly (PE, SE, OL, RPE)        PE=ddRAD & ezRAD pairedend, non-overlapping reads; SE=singleend reads; OL=ddRAD & ezRAD overlapping reads, miseq; RPE=oregonRAD, restriction site + random shear
rad             Cutoff1 (integer)                                   Use unique reads that have at least this much coverage for making the reference     genome
RAW-10-10       Cutoff2 (integer)                                   Use unique reads that occur in at least this many individuals for making the reference genome
0.05            rainbow merge -r <percentile> (decimal 0-1)         Percentile-based minimum number of seqs to assemble in a precluster
0.95            rainbow merge -R <percentile> (decimal 0-1)         Percentile-based maximum number of seqs to assemble in a precluster
------------------------------------------------------------------------------------------------------------------

----------mkBAM: Settings for mapping the reads to the reference genome-------------------------------------------
Make sure the cutoffs above match the reference*fasta!
1		bwa mem -A Mapping_Match_Value (integer) 			bwa mem default is 1
4		bwa mem -B Mapping_MisMatch_Value (integer) 			bwa mem default is 4
6		bwa mem -O Mapping_GapOpen_Penalty (integer) 			bwa mem default is 6
30		bwa mem -T Mapping_Minimum_Alignment_Score (integer) 		bwa mem default is 30. Remove reads that have an alignment score less than this. don't go lower than 1 or else the resulting file will be huge. NOTE! in fltrBAM settings (below) there is an alignment score filter that uses a threshold relative to read length.  This -T setting here affects which reads the relative alignment score threshold will be applied to.
5		bwa mem -L Mapping_Clipping_Penalty (integer,integer) 		bwa mem default is 5
------------------------------------------------------------------------------------------------------------------
```

---

## 8. Map reads to reference - Filter Maps - Genotype Maps

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkBAM

#this script has to be run from dir with fq.gz files to be mapped and the ref genome
#this script is preconfigured to run mapping, filtering of the maps, and genotyping in 1 shot
sbatch ../dDocentHPC.sbatch config.5.cssl
```

**NOTE:** Because this was the first sequencing run (and we did not know we would do a 2nd at the time) we "finished" this pipeline all the way through making the "all sites" VCF. However, the filtered `.bam` files from this step were later merged with those from the second sequencing run and used for the official full dataset. See the Gmi README for more information.

---

## 9. Filter VCF Files

Pulled latest changes from fltrVCF and rad_haplotyper repos

```
#if you haven't already, you first need to clone the fltrVCF.git repo & the rad_haplotyper.git repo
#cd /home/r3clark/PIRE/pire_cssl_data_processing/scripts
#git clone https://github.com/cbirdlab/fltrVCF.git
#git clone https://github.com/cbirdlab/rad_haplotyper.git

cd /home/r3clark/PIRE/pire_cssl_data_processing/scripts/fltrVCF
git pull

cd /home/r3clark/PIRE/pire_cssl_data_processing/scripts/rad_haplotyper
git pull

cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta
mkdir filterVCF

cp ../scripts/fltrVCF/config_files/config.fltr.ind.cssl filterVCF
```

Because many *Gazza minuta* contemporary individuals had low numbers of reads, removed individuals with <100K sequences from VCF prior to running `fltrVCF.sbatch`. This removes individuals with the most amount of missing data beforehand, allowing us to retain more sites (losing fewer SNPs due to too much missing data, etc.).

```
#created list of individuals to remove (35 total, all contemporary)

cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/filterVCF

module load vcftools

vcftools --vcf ../mkBAM/TotalRawSNPs.rad.RAW-10.10.vcf --remove indvfewsequences.txt --recode --recode-INFO-all --out TotalRawSNPs.rad.RAW-10-10.noindvless100Kseq
```

Ran `fltrVCF.sbatch`.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/filterVCF

#before running, make sure the config file is updated with file paths and file extensions based on your species
#config file should ONLY run up to the second 07 filter (remove filters 18 & 17 from list of filters to run)
sbatch ../fltrVCF.sbatch config.fltr.ind.cssl

#troubleshooting will be necessary
```

---

## 10. Check for cryptic species

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta
mkdir pop_structure
cd pop_structure

#copy final VCF file made from fltrVCF step to `pop_structure` directory
cp ../filterVCF/*Fltr07.18.vcf .
```

Ran PCA w/PLINK. Instructions for installing PLINK with conda are [here](https://github.com/philippinespire/pire_cssl_data_processing/blob/main/scripts/popgen_analyses/README.md).

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/pop_structure

module load anaconda
conda activate popgen

plink --vcf Gmi.A.rad.RAW-10-10.Fltr07.18.vcf --allow-extra-chr --pca --out PIRE.Gmi.Ham.preHWE
conda deactivate
```

Made input files for ADMIXTURE with PLINK.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/pop_structure

module load anaconda
conda activate popgen

plink --vcf Gmi.A.rad.RAW-10-10.Fltr07.18.vcf --allow-extra-chr --make-bed --out PIRE.Gmi.Ham.preHWE

awk '{$1=0;print $0}' PIRE.Gmi.Ham.preHWE.bim > PIRE.Gmi.Ham.preHWE.bim.tmp
mv PIRE.Gmi.Ham.preHWE.bim.tmp PIRE.Gmi.Ham.preHWE.bim
conda deactivate
```

Ran ADMIXTURE (K = 1-5). Instructions for installing ADMIXTURE with conda are [here](https://github.com/philippinespire/pire_cssl_data_processing/blob/main/scripts/popgen_analyses/README.md).

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/pop_structure

module load anaconda
conda activate popgen

admixture PIRE.Gmi.Ham.preHWE.bed 1 --cv > PIRE.Gmi.Ham.preHWE.log1.out #run from 1-5
conda deactivate
```

Copied `*.eigenval`, `*.eigenvec` & `*.Q` files to local computer. Ran `pire_cssl_data_processing/scripts/popgen_analyses/pop_structure.R` on local computer to visualize PCA & ADMIXTURE results (figures in `/home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/pop_structure`).

---

## 11. Filter VCF file for HWE

**NOTE:** If PCA & ADMIXTURE results don't show cryptic structure, skip to running `fltrVCF.sbatch`.

PCA & ADMIXTURE showed cryptic structure. ABas & CBas all assigned to one deme ("A"). ~50% of AHam & CBat assigned to same deme ("A") as ABas & CBas and ~50% assigned to separate deme ("B"). Species IDs unknown at this point.

Adjusted popmap file to reflect new structure.

```
cd /home/PIRE/pire_cssl_data_processing/gazza_minuta/filterVCF
cp ../mkBAM/popmap.rad.RAW-10-10 ./popmap.rad.RAW-10-10.HWEsplit

#added -A or -B to end of pop assignment (second column) to assign individual to either group A or group B.
```

Ran `fltrVCF.sbatch`.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/filterVCF
cp config.fltr.ind.cssl ./config.fltr.ind.cssl.HWE

#before running, make sure the config file is updated with file paths and file extensions based on your species
#popmap path should point to popmap file (*.HWEsplit) just made (if cryptic structure detected)
#vcf path should point to vcf made at end of previous filtering run (the file PCA & ADMIXTURE was run with)
#config file should ONLY run filters 18 & 17 (in that order)
sbatch ../fltrVCF.sbatch config.fltr.ind.cssl.HWE

#troubleshooting will be necessary
```

---

## 12. Make VCF with monomorphic loci

Moved the files needed for genotyping from `mkBAM` to `mkVCF`

```
#run in scratch if need more space
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta
mkdir mKVCF_monomorphic
mv mkBAM/*bam* mkVCF_monomorphic
cp mkBAM/*fasta mkVCF_monomorphic
cp mkBAM/config.5.cssl mkVCF_monomorphic
```

Changed the config file so that the last setting (monomorphic) is set to yes and renamed it with the suffix `.monomorphic`

```
yes      freebayes    --report-monomorphic (no|yes)      Report even loci which appear to be monomorphic, and report allconsidered alleles,
```

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkVCF_monomorphic
mv config.5.cssl config.5.cssl.monomrphic
```

Genotyped

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkVCF_monomorphic
sbatch ../dDocentHPC_mkVCF.sbatch config.5.cssl.monomorphic
```

---

## 13. Filter VCF with monomorphic loci

Will filter for monomorphic & polymorphic loci separately, then merge the VCFs together for one "all sites" VCF. Again, probably best to do this in scratch.

Like with the original filtering step, because many *Gazza minuta* contemporary individuals had low numbers of reads, removed individuals with <100K sequences from VCF prior to running `fltrVCF.sbatch`.

```
#used same list as before (35 total, all contemporary)

cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkVCF_monomorphic

module load vcftools

vcftools --vcf ../mkBAM/TotalRawSNPs.rad.RAW-10.10.vcf --remove ../filterVCF/indvfewsequences.txt --recode --recode-INFO-all --out TotalRawSNPs.rad.RAW-10-10.noindvless100Kseq

mv TotalRawSNPs.rad.RAW-10.10.noindvless100Kseq.recode.vcf TotalRawSNPs.rad.RAW-10.10.noindvless100Kseq.vcf
```

Set-up filtering for monomorphic sites only.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkVCF_monomorphic
cp ../scripts/config.fltr.ind.cssl.mono .
```

Ran `fltrVCF.sbatch` for monomorphic sites.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkVCF_monomorphic

#before running, make sure the config file is updated with file paths and file extensions based on your species
#VCF file should be the TotalRawSNPs.noindvless100Kseq file made after the "make monomorphic VCF" step
#settings for filters 04, 14, 05, 16, 13 & 17 should match the settings used when filtering the original VCF file (step 9)
sbatch ../fltrVCF.sbatch config.fltr.ind.cssl.mono

#troubleshooting will  be necessary
```

Set-up filtering for polymorphic sites only.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkVCF_monomorphic
mkdir polymorphic_filter
cp ../../scripts/config.fltr.ind.cssl.poly .
```

Ran `fltrVCF.sbatch` for polymorphic sites.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkVCF_monomorphic/polymorphic_filter

#before running, make sure the config file is updated with file paths and file extensions based on your species
#VCF file should be the TotalRawSNPs.noindvless100Kseq file made after the "make monomorphic VCF" step
#popmap file should be file used in step 11, that accounts for any cryptic structure (*HWEsplit extension)
#settings should match the settings used when filtering the original VCF file (step 9)
sbatch ../../fltrVCF.sbatch config.fltr.ind.cssl.poly

#troubleshooting will be necessary
```

---

## 14. Merge monomorphic & polymorphic VCF files

Check monomorphic & polymorphic VCF file sto make sure that filtering removed the same individuals. If not, remove necessary individuals from files.

* mono.VCF: filtering removed ABas_001, ABas_004, ABas_008, ABas_028, ABas_030, ABas_033, ABas_034, ABas_035, ABas_045, ABas_047, ABas_049, ABas_055, ABas_058, ABas_059, ABas_060, AHam_001, AHam_015, AHam_074, CBas_013, CBas_015, CBas_017, CBas_019, CBas_034, CBas_038, CBas_041, CBas_068, CBas_069, CBas_082, CBas_092, CBas_095, CBat_002, CBat_013, CBat_036, CBat_047, CBat_058, CBat_062, CBat_070, CBat_071, CBat_093, CBat_094, CBat_095, ABas_017, AHam_002, CBas_014, CBas_049. Need to remove CBat_061 as well to match polymorphic VCF.
* poly.VCF: filtering removed ABas_001, ABas_004, ABas_008, ABas_028, ABas_030, ABas_033, ABas_034, ABas_035, ABas_045, ABas_047, ABas_049, ABas_055, ABas_058, ABas_059, ABas_060, AHam_001, AHam_015, AHam_074, CBas_013, CBas_015, CBas_017, CBas_019, CBas_034, CBas_038, CBas_041, CBas_068, CBas_069, CBas_082, CBas_092, CBas_095, CBat_002, CBat_013, CBat_036, CBat_047, CBat_058, CBat_062, CBat_070, CBat_071, CBat_093, CBat_094, CBat_095, ABas_017, AHam_002, CBas_014, CBas_049, CBat_061

Created `indv_missing.txt` in `mkVCF_monomorphic` directory. This is a list of all the individuals removed from either  file (total of 46 for Gmi). Used this list to make sure number of individuals matched in both filtered VCFs.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkVCF_monomorphic
mv polymorphic_filter/gmi.poly.rad.RAW-10.10.Fltr17.20.recode.vcf .

module load vcftools

vcftools --vcf gmi.mono.rad.RAW-10.10.Fltr17.11.recode.vcf --remove indv_missing.txt --recode --recode-INFO-all --out gmi.mono.rad.RAW-10.10.Fltr17.11.recode.nomissing
mv gmi.mono.rad.RAW-10.10.Fltr17.11.recode.nomissing.recode.vcf gmi.mono.rad.RAW-10.10.Fltr17.11.recode.nomissing.vcf

vcftools --vcf gmi.poly.rad.RAW-10.10.Fltr17.20.recode.vcf --remove indv_missing.txt --recode --recode-INFO-all --out gmi.poly.rad.RAW-10.10.Fltr17.20.recode.nomissing
mv gmi.poly.rad.RAW-10.10.Fltr17.20.recode.nomissing.recode.vcf gmi.poly.rad.RAW-10.10.Fltr17.20.recode.nomissing.vcf
```

Sorted each VCF file.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkVCF_monomorphic

module load vcftools

#sort monomorphic
vcf-sort gmi.mono.rad.RAW-10.10.Fltr17.11.recode.nomissing.vcf > gmi.mono.rad.RAW-10.10.Fltr17.11.recode.nomissing.sorted.vcf

#sort polymorphic
vcf-sort gmi.poly.rad.RAW-10.10.Fltr17.20.recode.nomissing.vcf > gmi.poly.rad.RAW-10.10.Fltr17.20.recode.nomissing.sorted.vcf
```

Zipped each VCF file.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkVCF_monomorphic

module load samtools/1.9

#zip monomorphic
bgzip -c gmi.mono.rad.RAW-10.10.Fltr17.11.recode.nomissing.sorted.vcf > gmi.mono.rad.RAW-10.10.Fltr17.11.recode.nomissing.sorted.vcf.gz

#zip polymorphic
bgzip -c  gmi.poly.rad.RAW-10.10.Fltr17.20.recode.nomissing.sorted.vcf >  gmi.poly.rad.RAW-10.10.Fltr17.20.recode.nomissing.sorted.vcf.gz
```

Indexed each VCF file.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkVCF_monomorphic

module load samtools/1.9

#index monomorphic
tabix  gmi.mono.rad.RAW-10.10.Fltr17.11.recode.nomissing.sorted.vcf.gz

#index polymorphic
tabix  gmi.poly.rad.RAW-10.10.Fltr17.20.recode.nomissing.sorted.vcf.gz
```

Merged files.

```
cd /home/r3clark/PIRE/pire_cssl_data_processing/gazza_minuta/mkVCF_monomorphic

module load container_env bcftools
module load samtools/1.9

crun bcftools concat --allow-overlaps  gmi.mono.rad.RAW-10.10.Fltr17.11.recode.nomissing.sorted.vcf.gz  gmi.poly.rad.RAW-10.10.Fltr17.20.recode.nomissing.sorted.vcf.gz -O z -o gmi.all.recode.nomissing.sorted.vcf.gz

tabix gmi.all.recode.nomissing.sorted.vcf.gz #index all sites VCF for downstream analyses
