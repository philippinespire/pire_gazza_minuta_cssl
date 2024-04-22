# Gmi Data Processing Log

Log to track progress through capture bioinformatics pipeline for the Albatross and Contemporary *Gazza minuta* samples from Hamilo Cove & Basud River. Second sequencing run (contains individuals sequenced in the first run as well as some new individuals).

---

## 0. Rename files for dDocent HPC

Raw data in `/home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run/raw_fq_capture` (check `Gazza-minuta` channel on Slack). Starting analyses in `/home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run/raw_fq_capture`.

Used decode file from Sharon Magnuson & Chris Bird.

```bash
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run/raw_fq_capture

salloc
bash

#check got back sequencing data for all individuals in decode file
ls | wc -l #494 files (2 additional files for README & decode.tsv = XX/2 = XX individuals (R&F)
wc -l Gmi_CaptureLibraries2_SequenceNameDecode.tsv #245 lines (1 additional line for header = XX individuals), checks out

#run renameFQGZ.bash first to make sure new names make sense
bash /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/siganus_spinus/raw_fq_capture/renameFQGZc.bash Gmi_CaptureLibraries2_SequenceNameDecode.tsv

#run renameFQGZ.bash again to actually rename files
#need to say "yes" 2X
bash /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/siganus_spinus/raw_fq_capture/renameFQGZc.bash Gmi_CaptureLibraries2_SequenceNameDecode.tsv rename

```

---

## 1.  Check data quality with fastqc

Ran [`Multi_FASTQC.sh`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/Multi_FASTQC.sh).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run/raw_fq_capture

#Multi_FastQC.sh "<indir>" "<file_extension>"
sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/Multi_FASTQC.sh "/home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run/raw_fq_capture" "fq.gz"
```

Potential issues:  
  * % duplication - 
    * Alb: 71.13%, Contemp: 59.79%
  * GC content - 
    * Alb: 43%, Contemp: 47%
  * number of reads - 
    * Alb: ~15 mil, Contemp: ~10-12 mil

---

## 2. 1st fastp

Ran [`runFASTP_1st_trim.sbatch`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/runFASTP_1st_trim.sbatch).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run

#runFASTP_1st_trim.sbatch <INDIR/full path to files> <OUTDIR/full path to desired outdir>
sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runFASTP_1st_trim.sbatch /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run/raw_fq_capture/ fq_fp1
```

Potential issues:  
  * % duplication - 
    * Alb: 37.40%, Contemp: 56.61%
  * GC content -
    * Alb: 41.84%, Contemp: 47.24%
  * passing filter - 
    * Alb: 97.94%, Contemp: 97.41%
  * % adapter - 
    * Alb: 67.56%, Contemp: 25.02%
  * number of reads - 
    * Alb: ~25 mil, Contemp: ~20-25 mil

---

## 3. Clumpify

Ran [`runCLUMPIFY_r1r2_array.bash`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/runCLUMPIFY_r1r2_array.bash).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run

#runCLUMPIFY_r1r2_array.bash <indir;fast1 files > <outdir> <tempdir> <max # of nodes to use at once>
bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runCLUMPIFY_r1r2_array.bash fq_fp1 fq_fp1_clmp /scratch/mmalabag 10

```

Ran [`checkClumpify_EG.R`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/checkClumpify_EG.R) to see if any failed.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run

salloc
module load container_env mapdamage2

#had to install tidyverse package first
crun R
install.packages("tidyverse") #said yes when prompted, when finished, exited & didn't save env

crun R < /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/checkClumpify_EG.R --no-save
#all files ran successfully
```

All files ran successfully

---

## 4. 2nd fastp

Ran [`runFASTP_2_cssl.sbatch`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/runFASTP_2_cssl.sbatch).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run

#runFASTP_2_cssl.sbatch <INDIR/full path to clumpified files> <OUTDIR/full path to desired outdir>
sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runFASTP_2_cssl.sbatch fq_fp1_clmp fq_fp1_clmp_fp2
```

Potential issues:  
  * % duplication - 
    * Alb: 10.98%, Contemp: 21.71%
  * GC content - 
    *  Alb: 42.00%, Contemp: 47.09%
  * passing filter - 
    * Alb: 98.79%, Contemp: 98.64%
  * % adapter - 
    * Alb: 1.18%, Contemp: 0.50%
  * number of reads - 
    * Alb: ~20 mil, Contemp: ~15-19 mil

---

## 5. Run fastq_screen

Ran [`runFQSCRN_6.bash`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/runFQSCRN_6.bash).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run/

#runFQSCRN_6.bash <indir> <outdir> <number of nodes running simultaneously>
bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_procesing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 20
```

Checked that all files were successfully completed.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run

#checked that all 5 output files from fastqc screen were created for each file (should be XX for each = XX R1 & XX R2)
ls fq_fp1_clmp_fp2_fqscrn/*tagged.fastq.gz | wc -l #488
ls fq_fp1_clmp_fp2_fqscrn/*tagged_filter.fastq.gz | wc -l #488 
ls fq_fp1_clmp_fp2_fqscrn/*screen.txt | wc -l #490
ls fq_fp1_clmp_fp2_fqscrn/*screen.png | wc -l #488
ls fq_fp1_clmp_fp2_fqscrn/*screen.html | wc -l #488

#checked all out files for any errors
grep 'error' slurm-fqscrn.*out #nothing
grep 'No reads in' slurm-fqscrn.*out #nothing
```

Everything looks good, no errors/missing files.

Ran [`runMultiQC.sbatch`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/runMULTIQC.sbatch) separately.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run

#runMULTIQC.sbatch <indir> <report name>
#do not use trailing / in paths
sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runMULTIQC.sbatch fq_fp1_clmp_fp2_fqscrn fastqc_screen_report
```

Potential issues:

  * one hit, one genome, no ID - 
    * Alb: 92%, Contemp: 90.5%
  * no one hit, one genome to any potential contaminators (bacteria, virus, human, etc) - 
    * Alb: 5%, Contemp: 2%
    
---

## 6. Re-pair fastq_screen paired end files

Ran [`runREPAIR.sbatch`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/runREPAIR.sbatch).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run

#runREPAIR.sbatch <indir> <outdir> <threads>
sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_procesing/runREPAIR.sbatch fq_fp1_clmp_fp2_fqscrn fq_fp1_clmp_fp2_fqscrn_repaired 40
```

Once finished, ran [`Multi_FASTQC.sh`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/Multi_FASTQC.sh) to assess quality.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run

#Multi_FastQC.sh "<indir>" "<file_extension>"
sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/Multi_FASTQC.sh "/home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run/fq_fp1_clmp_fp2_fqscrn_repaired" "fq.gz"
```

Potential issues:  
  * % duplication - 
    * Alb: 39.73%, Contemp: 19.58%
  * GC content - 
    * Alb: 40%, Contemp: 47%
  * number of reads - 
    * Alb: ~9 mil, Contemp: ~7-8.5 mil

---

## 7. Calculate the percent of reads lost in each step

Executed [`read_calculator_cssl.sh`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/read_calculator_cssl.sh).

```
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run

#read_calculator_cssl.sh "<Path to species home dir>" "<Path to dir with species raw files>"
sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/read_calculator.sh "/home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run" "raw_fq_capture"  
```

Reads lost:

  * fastp1 dropped 2.30% of the reads
  * 51.49% of reads were duplicates and were dropped by Clumpify
  * fastp2 dropped 1.27% of the reads after deduplication

Reads remaining:

  * Total reads remaining: 42.52%

---

## 8. Set up mapping dir and get reference genome

Make mapping directory and move `*fq.gz` files over.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run
mkdir mkBAM

mv fq_fp1_clmp_fp2_fqscrn_repaired/*fq.gz mkBAM
```

Pulled latest changes from dDocentHPC repo & copied `config.5.cssl` over.

```sh
#if you have cloned, just pull the latest changes
cd /pire_cssl_data_procesing/scripts/dDocentHPC
git pull

cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run/mkBAM

cp ../../scripts/dDocentHPC/configs/config.5.cssl .
```

Since this is the 2nd sequencing run and GMI does not have an assembled genome, the "raw" reference fasta that was used during the 1st sequencing run was copied and used here. 

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run/mkBAM

cp ../1st_sequencing_run/reference.rad.RAW-10-10.fasta

#the destination reference fasta should be named as follows: reference.<assembly type>.<unique assembly info>.fasta
#<assembly type> is `ssl` for denovo assembled shotgun library or `rad` for denovo assembled rad library
#this naming is a little messy, but it makes the ref 100% tracable back to the source
#it is critical not to use `_` in name of reference for compatibility with ddocent and freebayes
```

Updated the config file with the ref genome info. The nano from the 1st sequencing run was copied over from 1st sequencing run as well. 

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run/mkBAM

nano config.5.cssl
```

Inserted `<assembly type>` into the `Cutoff1` variable and `<unique assembly info>` into the `Cutoff2` variable.

```
----------mkREF: Settings for de novo assembly of the reference genome--------------------------------------------
PE              Type of reads for assembly (PE, SE, OL, RPE)          PE=ddRAD & ezRAD pairedend, non-overlapping reads; SE=singleend reads; OL=ddRAD & ezRAD overlapping reads, miseq; RPE=oregonRAD, restriction site + random shear
rad             Cutoff1 (integer)                                     Use unique reads that have at least this much coverage for making the reference genome
RAW-10-10       Cutoff2 (integer)                                     Use unique reads that occur in at least this many individuals for making the reference genome
0.05            rainbow merge -r <percentile> (decimal 0-1)           Percentile-based minimum number of seqs to assemble in a precluster
0.95            rainbow merge -R <percentile> (decimal 0-1)           Percentile-based maximum number of seqs to assemble in a precluster
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

## 9. Map reads to reference

Since mkVCF was already run during the 1st sequencing run, dDocentHPC.sbatch script was duplicated and had the mkVCF commented out in the script under the new name dDocentHPC_nomkVCF.sbatch.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/gazza_minuta/2nd_sequencing_run/mkBAM

#this script has to be run from dir with fq.gz files to be mapped and the ref genome
#this script is preconfigured to run mapping, filtering of the maps, and genotyping in 1 shot
sbatch ../../dDocentHPC_nomkVCF.sbatch config.5.cssl
```

Stopped here as the 1st and 2nd sequencing run was combined and ran through filtering and pop structure steps together (see main Gmi README for more information).
