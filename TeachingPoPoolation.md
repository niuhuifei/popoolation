

# Introduction #

This is a Quick Guide intented to demonstrate the functionality of PoPoolation while requiring as little 'user-time' as possible!

We provide a very small prefiltered sample data set. If time is not an issue and you prefer a more thorough introduction using unfiltered data we recommend to use the extensive Walkthrough instead: http://code.google.com/p/popoolation/wiki/PoPOOLationWalkthrough

However the sample data set is the area flanking the Cyp6g1 gene in Drosophila melanogaster. Cyp6g1 has recently sweeped in D.mel as it confers resistance to DDT. This is a very complex genomic area where several transposons have been inserted (Accord, P-elements etc) and several copy number variations have been described.

# Data #

Download the small sample data set: http://popoolation.googlecode.com/files/teaching-data.zip

Unzip the folder

# Software #

Install the following tools:
  * IGV: http://www.broadinstitute.org/igv/
  * bwa: http://bio-bwa.sourceforge.net/bwa.shtml
  * samtools: http://samtools.sourceforge.net/
  * Perl: http://www.perl.org/
  * R: http://www.r-project.org/
  * PoPoolation: http://code.google.com/p/popoolation/


# Quick Tour #

## Lazy Man's Quick Guide ##

copy the following shell script into the unzipped folder 'data' containing the data of the quick guide:

http://popoolation.googlecode.com/files/popool-teaching.sh

enter the folder in the command line and type:

```
sh popool-teaching.sh <your-path-to-popoolation>
```

eg:
```
sh popool-teaching.sh ~/dev/popoolation
```

**NOTE**: do not provide the character "/" in the end of the path!

  * Check out the resulting files (like the sorted-genewise.pi or the cyp6g1.pdf)
  * Continue with loading the data into the IGV

## Normal Quick Guide ##

enter the command line and change directory to the folder 'data' mentioned above.

### Trim reads ###

As a first step trim the reads by quality as shown here

```
perl <local-popoolation-installation>/basic-pipeline/trim-fastq.pl --input1 read_1.fastq --input2 read_2.fastq --output trim --quality-threshold 20 --min-length 50
```

### Prepare the reference genome ###
The reference genome needs to be prepared:

```
mkdir wg
mv dmel-2R-chromosome-r5.22.fasta wg
awk '{print $1}' wg/dmel-2R-chromosome-r5.22.fasta > wg/dmel-2R-short.fa
bwa index wg/dmel-2R-short.fa
```

### Map the trimmed reads to the reference genome ###
```
bwa aln wg/dmel-2R-short.fa trim_1 > trim_1.sai
bwa aln wg/dmel-2R-short.fa trim_2 > trim_2.sai
bwa sampe wg/dmel-2R-short.fa trim_1.sai trim_2.sai trim_1 trim_2 > maped.sam
```

**NOTE**: As this guide is intented to be as fast as possible this are definitely not the optimal parameters for mapping of the reads; Instead the following parameters would be prefererable : bwa aln -l 100 -o 2 -d 12 -e 12 -n 0.01 wg/dmel-2R-short.fa trim\_1 > trim\_1.sai

### Create a pileup file ###

Extract reads with a mapping quality of at least 20 (unambiguously mapped reads) and create a sorted bam file.

```
samtools view -q 20 -bS maped.sam| samtools sort - maped.sort
```

Create a pileup file

```
samtools pileup maped.sort.bam > cyp6g1.pileup
```

### Calculate Tajima's Pi using a sliding window approach ###
```
perl <local-popoolation-installation>/Variance-sliding.pl --measure pi --input cyp6g1.pileup --min-count 2 --min-qual 20 --min-coverage 4 --max-coverage 70 --pool-size 500 --window-size 1000 --step-size 1000 --output cyp6g1.varslid.pi --region 2R:7800000-8300000
```

### Create a small overview ###
This step creates a pdf showing a quick overview of pi for chromosome 2R. This step is optional.

```
perl <local-popoolation-installation>/Visualise-output.pl --input cyp6g1.varslid.pi --output cyp6g1.pdf --chromosomes "2R" --ylab pi
```

### Prepare the output for visualization using the IGV ###

Convert the output of Variance-sliding.pl into a wiggle file:

```
perl <local-popoolation-installation>/VarSliding2Wiggle.pl --input cyp6g1.varslid.pi --output cyp6g1.pi.wig --trackname "nat-pop-pi"
```

Index the bam file:

```
samtools index maped.sort.bam
```


### Calculate a genewise pi ###

```
perl <local-popoolation-installation>/Variance-at-position.pl --measure pi --pileup cyp6g1.pileup --gtf cyp6g1.gtf --output genewise.pi --pool-size 500 --min-count 2 --min-coverage 4 --max-coverage 70 --min-qual 20
```

Sort the output showing the genes with the smallest Tajima's Pi on top.

```
sort -k 4,4n genewise.pi> sorted-genewise.pi
```

view the output
```
less sorted-genewise.pi
```

# Display the results in the IGV #

  * Open the IGV

## Import a genome ##

  * Press Import Genome...
  * Enter Name: Dmel-2R
  * Enter sequence file: dmel-2R-chromosome-[r5](https://code.google.com/p/popoolation/source/detail?r=5).22.fasta
  * Press Save

## Load the data ##

  * Press File, Load from File..
  * select the file: maped.sort.bam
  * Similarly load the file: cyp6g1.gtf
  * Finally load: cyp6g1.pi.wig

## View the results ##

  * Zoom in on the whole region: 2R:7800000-8300000
  * Zoom in on Cyp6g1: CG8453-RA

# Compare with our results #
Download the following archive and compare your results with these screenshots: http://popoolation.googlecode.com/files/cyp6g1-screen.zip

# Concluding remarks #

The Cyp6g1 gene is sweeped in D.melanogaster as it confers resistance to DDT. Tajima's Pi should be very low around Cyp6g1. However there is a small peak of Tajima's Pi directly in the Cyp6g1 gene.
The Cyp6g1 gene has also been duplicated in natural populations of D.melanogaster, thus the peak of Tajiama's Pi is not caused by real SNPs in the population rather copy number variation and copy variable sites (CVS) are responsible for this 'artificial' peak. It may thus be necessary to adjust the maximum coverage to a more conservative value of for example '50'.


**In summary this example shows that both the minimum and the maximum coverage needs to be choosen carefully**

Enjoy using PoPoolation,

yours Robert