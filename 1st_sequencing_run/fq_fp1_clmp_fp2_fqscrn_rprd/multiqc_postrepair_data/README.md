# correcting naming format for metadata summary

the names of these files are both non-standard format and dont match those from other directories in this project.  Consequently, I did the following to change th$

the root cause of this has been solved, and this is an anomly, so changing the metadata file read into R is the simplest solution

```bash
# commands run from this dir

# backup metadata file
cp multiqc_fastqc.txt multiqc_fastqc_backup.txt

# then I made seqid_sample.txt from the raw fqgz dir
cat <(echo -e 'seq_id\tSample') \
<(paste <(ls ../../raw_fq_capture/*1.fq.gz | sed -e 's/^.*\///' -e 's/_.*$//') <(ls ../../raw_fq_capture/*1.fq.gz | sed -e 's/^.*\///' -e 's/\(_L[1-9]\)_.*$/\1/')) > seqid_sample.txt

# then I replaced the sample in multiqc_fastqc_backup.txt with that in seqid_sample.txt
# there is a bit of code here to address the situation where there are no reads which did happen here.  May need to apply elsewhere, the sed null line.
paste <(sed 's/^Gmi\-...._...\-//' multiqc_fastqc_backup.txt | sed 's/\-.*$//' | sed 's/Sample.*$/seq_id/' ) \
<(cat multiqc_fastqc_backup.txt | cut -f2- | sed -e 's/ /@/g') | \
awk 'NR==FNR{a[$1]=$1OFS$2;next}{$1=a[$1];print}' OFS='\t' seqid_sample.txt - | \
grep -Pv "^\t" | \
sed 's/\tnull\t\(............\)\t\(...\)\t\(...\)\t/\tnull\t\1\t\2\t\3\tNA\t/' |
sed 's/@/ /g' | \
cut -f2- > multiqc_fastqc.txt
