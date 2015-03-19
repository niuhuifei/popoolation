



# Data set #

| file | read number|read count| read length|
|:-----|:-----------|:---------|:-----------|
| s\_7\_1\_sequence\_Jul2009.fastq|1 |17957720|75|
| s\_7\_2\_sequence\_Jul2009.fastq|2 |17957720|75|

**NOTE**: Data available at short read archiv http://trace.ddbj.nig.ac.jp/DRASearch/submission?acc=SRA023610

**NOTE**: A very small sample data set can alternatively be found in the Quick Guide: http://code.google.com/p/popoolation/wiki/TeachingPoPoolation

# Requirements #
see [Manual](Manual.md)


# Walkthrough #

## Download PoPoolation ##

PoPoolation may be obtained directly from the subversion repository. Just go to the directory where you want to install PoPoolation in the command line and enter the command:
```
svn checkout http://popoolation.googlecode.com/svn/trunk/ popoolation
```

To update your copy of PoPoolation with the latest improvements enter your PoPoolation directory and enter the command:
```
svn update
```

Alternatively PoPoolation may be downloaded from the project main page: http://code.google.com/p/popoolation/

**However we recommend to use subversion as bugfixes will be immediately available in the repository**

## Trimming of the reads ##
```
perl trim-fastq.pl --input1 s_7_1_sequence_Jul2009.fastq --input2 s_7_2_sequence_Jul2009.fastq --quality-threshold 20 --min-length 50 --output dmel_trimed
```
### Resulting trimming statistics ###

```
FINISHED: end statistics
Read-pairs processed: 17957720
Read-pairs trimmed in pairs: 13495995
Read-pairs trimmed as singles: 3802238


FIRST READ STATISTICS
First reads passing: 15656946
5p poly-N sequences trimmed: 16377
3p poly-N sequences trimmed: 61990
Reads discarded during 'remaining N filtering': 0
Reads discarded during length filtering: 2300774
Count sequences trimed during quality filtering: 8551460

Read length distribution first read
length	count
50	110528
51	112155
52	115546
53	119517
54	120622
55	123213
56	127909
57	132907
58	139580
59	143970
60	151528
61	166721
62	178000
63	185352
64	192283
65	202482
66	221528
67	247456
68	270715
69	312185
70	364457
71	461386
72	830395
73	1220619
74	9405892


SECOND READ STATISTICS
Second reads passing: 15137282
5p poly-N sequences trimmed: 99060
3p poly-N sequences trimmed: 35878
Reads discarded during 'remaining N filtering': 0
Reads discarded during length filtering: 2820438
Count sequences trimed during quality filtering: 9051269

Read length distribution second read
length	count
50	103225
51	102568
52	108397
53	117810
54	130272
55	146265
56	153847
57	166468
58	182299
59	191046
60	199036
61	192029
62	223088
63	224462
64	190667
65	183919
66	190513
67	202838
68	230031
69	269701
70	360601
71	460748
72	773657
73	1181002
74	8852793

```


## Mapping of reads (BWA) ##

### Prepare reference sequence ###
First obtain a reference genome of D. melanogaster from http://flybase.org/

We used: dmel-all-chromosome-[r5](https://code.google.com/p/popoolation/source/detail?r=5).22.fasta.gz

Remove everything after the first whitespace from the reference genome. This is a precautioniary measure as some mappers and software downstream of mapping have difficulties with fasta-ids containing whitespace.
```
awk '{print $1}' dmel-all-chromosome-r5.22.fasta > dmel-short-header.fa
```

Index the reference sequence
```
bwa index dmel-short-header.fa
```
### Mapping using bwa ###

Assuming the reference sequence is in the folder 'wg' the following command will map the reads which have been trimmed in a pair to the reference:

```
bwa aln -t 3 -o 2 -d 12 -e 12 -l 100 -n 0.01 wg/dmel-short-header.fa dmel_trimed_1 > dmel_trimed_1.sai
bwa aln -t 3 -o 2 -d 12 -e 12 -l 100 -n 0.01 wg/dmel-short-header.fa dmel_trimed_2 > dmel_trimed_2.sai
```

Converting mapping results to a sam file
```
bwa sampe wg/dmel-short-header.fa dmel_trimed_1.sai dmel_trimed_2.sai dmel_trimed_1 dmel_trimed_2 > dmel.sam
```



## Filter reads by mapping quality and convert to a pileup file ##

### Convert sam-file into a sorted bam-file ###

Filter by a mapping quality of 20 and convert the sam file into a sorted bam file.
Filtering by a mapping qualiy of 20 removes the ambiguously mapped reads.
```
samtools view -q 20 -b -S dmel.sam|samtools sort - dmel.sort
```
### Crosscheck ###
```
samtools flagstat dmel.sort.bam
```

The result should be this
```
22467169 in total
0 QC failure
0 duplicates
22467156 mapped (100.00%)
22467169 paired in sequencing
11234128 read1
11233041 read2
22324803 properly paired (99.37%)
22439809 with itself and mate mapped
27347 singletons (0.12%)
73897 with mate mapped to a different chr
73897 with mate mapped to a different chr (mapQ>=5)
```

### Convert the sorted bam file into a pileup file ###
Convert the bam file into a pileup file:
```
samtools pileup dmel.sort.bam > dmel.pileup
```

## Run PoPoolation ##

### Run Variance-sliding.pl ###

Run the script Variance-sliding.pl with the dme.pileup file requesting Tajima's Pi. We use a window size and a step size of 10000.

```
perl Variance-sliding.pl --input dmel.pileup --output dmel.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 500
```

We assume a population size of 500 individuals. Furthermore we require a minimum allele count of 2, a minimum base quality of 20, a minimum coverage of 4 and a maximum coverage of 400.

### Visualise output of Variance-sliding.pl ###

#### Create an overview ####

Create the overview:
```
perl Visualise-output.pl --input dmel.pi --output dmel.pi.pdf --ylab pi --chromosomes "X 2L 2R 3L 3R 4"
```

and check out the result:

http://popoolation.googlecode.com/files/dmel.pi.pdf

#### Use IGV ####
The IGV may be downloaded from: http://www.broadinstitute.org/igv/

First convert the Drosophila melanogaster pi-file into a wiggle file
```
perl ~/dev/PopGenTools/VarSliding2Wiggle.pl --input dmel.pi --trackname "dmel Pi" --output dmel.pi.wig
```

than index the bam file:
```
samtools index dmel.sort.bam
```

Than:
  1. open the IGV.
  1. switch to the genome D. melanogaster reference 5.22
  1. load the sorted bam file:  dmel.sort.bam
  1. load the wiggle file: dmel.pi.wig

and check out the results

http://popoolation.googlecode.com/files/igv_ex1.pdf

http://popoolation.googlecode.com/files/igv_ex2.pdf




### Run Variance-at-position ###

#### Obtain and prepare a annotation ####
Go to the FlyBase homepage (http://flybase.org/) and get the annotation: dmel-all-[r5](https://code.google.com/p/popoolation/source/detail?r=5).22.gff.gz

Unzip the file.

Filter for exons and convert it into a gtf file:
```
cat dmel-all-r5.22.gff| awk '$2=="FlyBase" && $3=="exon"'| perl -pe 's/ID=([^:;]+)([^;]+)?;.*/gene_id "$1"; transcript_id "$1:1";/'> exons.gtf
```

**Note:** the script Variance-at-position.pl will group all exons with the same gene\_id; the transcript\_id is not considered

**Note:** UCSC (the Tables section) already provides gtf formated annotation -> if you use those and the proper reference sequence no reformating would be necessary!

Run Variance-at-position.pl
```
perl Variance-at-position.pl --pool-size 500 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 4 --pileup dmel.pileup --gtf exons.gtf --output dmel.genes.pi --measure pi
```
We set the pool size to 500, request pi as Population Genetic estimator, use the previously generated pileup file and the reformated exon annotation.


Output will be something like the following:
```
FBgn0042083	7	0.961	0.002524220
FBgn0027066	9	0.970	0.001810036
FBgn0033100	3	0.989	0.001972413
FBgn0033101	2	0.664	0.001424038
CG9438	8	0.900	0.003091130
FBgn0085421	14	0.936	0.003365280
```
  * column 1: the gene id (or id of the region)
  * column 2: the number of SNPs found in the gene
  * column 3: covered fraction; how much of the gene is sufficiently covered by the pileup file; values 0-1; 0..not a single base of the gene has the required coverage, 1.. all bases of the gene have the required coverage
  * column 4: the measure; here Tajima's pi

**NOTE:** This script is loading the whole annotation into the memory and internally converts it into a per-base hash. For very large gtf files this can be quite memory demanding. In case you do not have enough memory you may split your gtf file and analyse them separately. As long as you are not splitting genes, this will not affect the  outcome.
For example the following command will retrieve only entries from chromosome 2R;
```
cat exons.gtf| awk '$1=="2R"' > 2R-exons.gtf
```