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
#-r CHR:FROM-TO  Only report depth in specified region.
```
