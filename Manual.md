

# Requirements #

  * Perl 5.8 or higher
  * R 2.7 or higher
  * a short read aligner (e.g.: bwa)
  * samtools
  * Mauve (http://asap.ahabs.wisc.edu/mauve/)

# With which data can I use PoPoolation #

The software is designed for next generation sequencing data of pooled genomic DNA. Moreover at least a partial assembly of the used species should be available.

We do not recommend using the software with cDNA or unigenes as these sequences do not contain introns which will negatively affect the mapability of genomic reads. Furthermore we do not recommend to use the software on pooled RNA samples as RNA is subject to many modifications which will distort the population genetic measures, like differential gene expression, RNA editing, alternative splicing, allele specific gene expression.

# Download PoPoolation #

PoPoolation may be obtained directly from the subversion repository. Just go to the directory where you want to install PoPoolation in the command line and enter the command:
```
svn checkout http://popoolation.googlecode.com/svn/trunk/ popoolation
```

To update your copy of PoPoolation with the latest improvements enter your PoPoolation directory and enter the command:
```
svn update
```

Alternatively PoPoolation may be downloaded from the project main page: http://code.google.com/p/popoolation/

However we recommend to use subversion as bugfixes will be immediately available in the repository.



# PoPoolation scripts #

## General note ##

Help for all PoPoolation scripts can be requested using the '--help' option.
For example
```
perl trim-fastq.pl --help
```

The main scripts of PoPoolation contain Unit tests which can be run using the '--test' option
For example:
```
perl trim-fastq.pl --test
```

## Variance-sliding.pl ##

Calculates Tajima's Pi or Wattersons Theta or Tajima's D along chromosomes using a sliding window approach.
The variable sites (SNPs) used for calculating the requested Population Genetic measure may be output to a separate file. The user can choose a step size, a window size, the minimum base quality, the minium allele count, a minimum coverage and a maximum coverage.

## Variance-at-position.pl ##

Calculates Tajima's Pi or Wattersons Theta or Tajima's D for target regions. The target regions have to be provided as gtf file. Genes frequently consist of several exons, the script groups all regions which have the same "gene\_id". The script is also able to handle overlapping regions. If the overlapping regions belong to the same gene (e.g.: sometimes an exon is reported for every transcripts belonging to the same gene), the overlapping region is only considered once. If the overlapping regions belong to different genes they are treated ndependently.
If requested the SNPs used for calculating the Population Genetic measure may be output to distinct file.
Again the user may choose a minimum base quality, a minimum allele count, a minimum coverage and a maximum coverage.

## Visualise-output.pl ##

Converts the output of Variance-sliding.pl into a pdf. Thus a overview of the general pattern of polymorphism for the provided chromosomes (contigs) may be quickly created.

## VarSliding2Wiggle.pl ##

Converts the output of Variance-sliding.pl into a wiggle file (http://genome.ucsc.edu/goldenPath/help/wiggle.html). Wiggle files are for example accepted by the Integrative Genomics Viewer (http://www.broadinstitute.org/igv/) or by the UCSC Genome Browser (http://genome.ucsc.edu/).

## VarSliding2Flybase.pl ##

Converts the output of Variance-sliding.pl into a format which is compatible with the FlyBase Genome browser.

## mauve-parser.pl ##

Parses the output of the Mauve Genome Alignment Software (http://asap.ahabs.wisc.edu/mauve/).
To determine the ancestral allele it is necessary to align three species where one acts as an outgroup, this may be accomplised with Mauve. This script converts a Mauve multiple genome alignment into a single tab delimited table where one genome (the first) acts as the reference.

## calculate-divergence.pl ##

Calculates the divergence of the species along chromosomes using a sliding window approach. This scripts requires the parsed Mauve output (or any other Multiple alignment software as long as the above mentioned tab delimited table is produced; For more information use: perl calculate-divergence.pl --help).

## basic-pipeline/trim-fastq.pl ##

This script trims reads in the fastq file format and outputs again a fastq file. Reads may either be trimmed as singel end reads or paired end reads. In case of paired end reads three output files are produced one for the first read, one for the second read and one for all reads which lost their mate because of not fulfilling the requirements (i.e.: they have become single end reads). First a trimming of the character 'N' is performed and secondly a quality trimming using  a modified Mott algorithm as implemented by Phred and described here (http://www.phrap.org/phredphrap/phred.html) is performed. Reads which do not have the required minimum length are discarded. If only one read of a pair has the required minimum length the paired end will be broken and the read will be reported as single end (SE).

## basic-pipeline/mask-sam-indelregions.pl ##

Indels frequently cause misalignments and thus may lead to false positive SNPs. This script masks the region surrounding an indel in the sam file. The sam file may be converted into a pileup file (samtools) and used with PoPoolation

## basic-pipeline/identify-genomic-indel-regions.pl ##

Identifies indel regions from a pileup file and outputs the coordinates into a gtf file. This gtf file may for example be used to filter regions surrounding the indel from a pileup file (see basic-pipeline/filter-pileup-by-gtf.pl).


## basic-pipeline/filter-pileup-by-gtf.pl ##
Allows to filter a pileup file by the regions specified in a gtf file. The specified regions may either be kept or discarded. For example the output of RepeatMasker may be used directly with this script to remove repetitive regions. Or as another example exons/CDS sequences may be removed before calculating Tajima's D. Or another examples whole regions surrounding indels may be removed (see basic-pipeline/identify-genomic-indel-regions.pl). The pileup file may subsequently be use with PoPoolation

## basic-pipeline/convert-fastq.pl ##
Convert between the different quality encodings of a fastq file. For example convert a illumina encoded fastq file to a sanger encoded fastq file


## basic-pipeline/subsample-pileup.pl ##

May be used to create a uniform coverage in a pileup file. The bases of a pileup file may be randomly sampled to reduce the coverage to the given threshold. Two methods may be used for subsampling: either random drawing without replacement or a simple down scaling.

# Calculating Tajima's D #

Popoolation offers two options for calculating Tajima's D. For both options we recommend to subsample the coverage to a uniform coverage.

  * Classic Tajima's D; minimum count needs to be set to 1 (--dissable-corrections)
  * Corrected Tajima's D; minimum count needs to be set to 2; The work of Achatz (2008; http://www.ncbi.nlm.nih.gov/pubmed/18562660) is the basis of our correction method; Requirement: Poolsize > 3 x minimum coverage

## Detailed description of Tajima's D ##


As a basis for the equation for calculating Tajima's D we are using the results from Achaz (2008; http://www.ncbi.nlm.nih.gov/pubmed/18562660) for individual sequencing when singletons are omitted (detailed equations please see http://popoolation.googlecode.com/files/correction_equations.pdf). The difficulty arises when calculating the variance of Tajima's D.
The following simplifying assumptions are made in our variance approximation:
  * The locally fluctuating coverage is replaced by the minimum coverage. This makes the variance estimator larger, and therefore leads to conservative estimates of Tajima’s D.
  * The random number of different individuals sequenced (nbase) under a given coverage is replaced by its expected value. This assumption should not affect the results much, if the pool size is large compared to the coverage, sequencing the same individual more than once is uncommon under these circumstances. Furthermore the number of different individuals sequenced will have a low variance.  As we are working with the minimum coverage, nbase will be biased downwards—tending to give a conservative estimate oft he variance.
  * Different subsets from the pool are sequenced at different positions. Their coalescent histories will be correlated but not identical. As the classical formula for Tajima’s D are for a single sample sharing a common coalescent history, there is more independence in the data than assumed with the classical formula—which should again make the variance approximation more conservative.
Summing up, the approximate variance in our equations provides a conservative approximation, and the values for Tajima’s D will tend to be smaller than those that would be expected for an experiment based on individual sequencing of a single sample.

For a correlation between the classical Tajima's D and the corrected Tajima's D using the equation described above please see (x: classical, y: corrected). Correlation has been done with real data from _Drosophila_ with a coverage of 12, a window size of 500 and a minimum count of 1) ![http://popoolation.googlecode.com/files/correlation_classic_correctedTajimasD.pnghttp://popoolation.googlecode.com/files/correlation_classic_correctedTajimasD.png](http://popoolation.googlecode.com/files/correlation_classic_correctedTajimasD.pnghttp://popoolation.googlecode.com/files/correlation_classic_correctedTajimasD.png)

# After PoPoolation #

What to do after you successfully used PoPoolation, having identified several hundreds (or thousands) of SNPs with unusual low (or high) Tajima's D or Watterson's Theta. You may use Gowinda to test whether your SNPs show an enrichment for any GO category (or any gene set): http://code.google.com/p/gowinda/