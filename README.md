
# Description

Tool for integrating cutadapt into the gatk workflow instead of the default MarkIlluminaAdapters. 

This aims to integrate cutadapts flexable but rigid trimming for the lossless bam workflow of the gatk. By doing so it allows more workflow fexability for example trimming non-illumina prepkit adapters (read "Twist UMI Adapter Systems"). 

# Installation

A basic `setup.py` will do

```
python setup.py
```

# Example

Run picard FastqToSam (or fgBio fastqToSam + sortSam -so "queryname") to create a queryname_sorted_unaligned.bam from the initial raw fastqs.  

Then dump those reads with picard SamToFastq to a file and use the outputs as input for cutadapt or any other trimming program retaining 0 length reads and file ordering. 

Finally run this tool to integrate the cutadapt results into the bam fixing any desyncing issues with the umis or other added barcodes to the sam tags. As shown below for pipeline integration.

```sh
cat queryname_sorted_unaligned.bam | marktrimming.py --input - --output - --fastq cutadapt_trimmed_1.fastq.gz --fastq cutadapt_trimmed_2.fastq.gz --ubam-out | samtools view -b > queryname_sorted_unaligned_adaptermarked.bam
```

This can be used with the picard toold SamToFastq marking adapter trimmed reads as low quality for bwa input. Finally producing a an alignment file with picards MergeBamAlignment.

# debugging

```bash
(ml marktrimming && python3 setup.py install --prefix  tests/install/  && export PYTHONPATH='./tests/install/lib/:'"$PYTHONPATH" && python ./build/lib/marktrimming/marktrimming.py --help)
```