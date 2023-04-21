# Gene Expression Cookbook
Helpful snippets on gene expression analysis


### Samtools manual  
If you have a ```bam``` file, it has to be sorted by coordinates and then indexed to enable most of the samtools functionality.  
```bash
samtools sort -@ 20 in.bam -o sorted.bam
#-@  number of processors
#-o  output file

samtools index sorted.bam
# creates sorted.bai file in the same folder
```
Bam file can contain alignments to the genome or to the transcriptome. In the former case it would have each sequencing read with a corresponding chromosome coordinate. In the latter case it would have transcrptomic coordinates, which is useful if you need to extract the coverage depth over a transcript or a genomic region.
```bash
# example of extracting coverage depth for a particular transcript from STAR output (in -quantMode TranscriptomeSAM)
samtools depth -a [in.bam|in.sam] -r 'rna-NM_007463.4' -o out.file
#-a  Output all positions (including those with zero depth)
#-aa Output absolutely all positions, including unused reference sequences.
#-o  Write output to FILE
#-r CHR:FROM-TO  Only report depth in a specified region. If aligned in transcriptome coordinates, the CHR is the name of the transcript.

# example of extracting reads specific to a single transcript. 
samtools view [in.bam|in.sam] transcriptName >output.sam
# requires position sorted bam or sam imput file
# transcript name should be in the format TRANSCRIPT:FROM-TO or just TRANSCRIPT
```

Create separate folders for each pair of sequencing files and transfer them to these folders.  
```bash
# cd to the directory containing *fastq.gz files (R1 and R2 pair)
for file in *_1.fastq.gz; do
  base="${file%_1.fastq.gz}"
  r1="$file"
  r2="$base""_2.fastq.gz"
  mkdir "$base"
  mv "$r1" "$base"
  mv "$r2" "$base"
done
```
Illumina adapter trimming in batch  
```bash
cutadapt -j 20 -m 75 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -o trimmed_1.fq.gz -p trimmed_2.fq.gz read.1.fq.gz read.2.fq.gz
# -j      - number of threads
# -m      - discard read pair if any of the mates if shorter than 50 nucleotides after adapter trimming

#automate with bash for loop if needed. In the folder containing sample subfolders (R1 and R1 files in each subfolder) run this:
for dir in */; do
file1=$(find "$dir" -name "*_1.fq.gz");
file2=$(find "$dir" -name "*_2.fq.gz");
echo "$(cutadapt -j 20 -m 75 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -o "$dir""trimmed_1.fq.gz" -p "$dir""trimmed_2.fq.gz" "$file1" "$file2")";
done
```
Count reads per gene in sam files
```bash
grep -ve '^@' aligned.sam | cut -f 3 | sort -bgr | uniq -c
#### grep parameters ###
# e   - exclude lines with the pattern (in this case exclude lines that start with @ symbol, which marks headers in sam)
# v   - take regular expressions as pattern
#### cut parameters ###
# f   - select a column/field
#### sort parameters, easy to remember as 'bigger' ###
# g   - general numeric sort
# b   - ignore leading blanks
# r   - sort from smaller to larger
### uniq parameters ###
# c   - count the identical lines (until a different line is met). Returns the line and the count.

```

Remove non-templated nucleotide addition from the 5' end. 
```
# a line from bowtie-1 output (native format)
# notice 0:G>T field which indicates a mismatch in the first position of the read. It is very common in ribosome profiling library preparations due to non-templated nucleotide addition during the reverse transcription phase.
VH01141:12:AACG2FYM5:1:1101:55042:1189_TAGGGTTA 1:N:0:GGATGT	+	NR_003287.4_28S	3838	TTCCCTACCTACTATCCAGCGAAACCACAGCC	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;	0	0:G>T
# this is how it should look like after removing the mismatch
VH01141:12:AACG2FYM5:1:1101:55042:1189_TAGGGTTA 1:N:0:GGATGT	+	NR_003287.4_28S	3839	TCCCTACCTACTATCCAGCGAAACCACAGCC	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;	0
```

```bash
awk -F"\t" '{if ($8 ~ "0:[A-Z]>[A-Z]$") {print $1"\t"$2"\t"$3"\t"$4+1"\t"substr($5,2)"\t"substr($6,2)"\t"$7} else {print $0} }' aligned_reads.bwt >aligned_reads_corrected.bwt


```
