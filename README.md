
# Description

Tool for integrating cutadapt into the gatk workflow instead of the default MarkIlluminaAdapters. 

Create an unmapped queryname sorted bam using fastqToSam from the picard toolkit or similar tools. Then trim the SamToFastq ouput using cutadapt (emitting all reads in the same order! not omitting 0 length reads ) to output the lossy trim. Then use this tool to annotate the XT tags with the location of possible adapter. Then finally use FastqToSam on the annotated unaligned bam to mark the bases as low quality(basequality=2).

# Installation

A basic `setup.py` will do

```
python setup.py
```

# Example

# run FastqToSam to creat a queryname sorted unaligned.bam 
# first run SamToFastq and cutadapt on the resulting fastqs
```sh
cat in_unaligned_queryname_sorted.bam| marktrimming.py --input - --output - --fastq cutadapt_trimmed_1.fastq.gz --fastq cutadapt_trimmed_2.fastq.gz --ubam-out | 