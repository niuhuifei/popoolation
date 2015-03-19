# Requirements #

see [Manual](Manual.md)

# Data #

Download the full genomes of _D. melanogaster_, _D. simulans_ and _D. yakuba_ from FlyBase (http://flybase.org/).
used files:
  * dmel-all-chromosome-[r5](https://code.google.com/p/popoolation/source/detail?r=5).29.fasta.gz
  * dsim-all-chromosome-[r1](https://code.google.com/p/popoolation/source/detail?r=1).3.fasta.gz
  * dyak-all-chromosome-[r1](https://code.google.com/p/popoolation/source/detail?r=1).3.fasta.gz

Unzip the files

# Walkthrough #

## Prepare reference genomes ##

Many downstream software tools have problems with headers containing any whitespace. As a cautionary measure it is recommended to keep only the not-whitespace characters of the fasta headers.

```
awk '{print $1}' dmel-all-chromosome-r5.29.fasta > dmel-short.fasta
awk '{print $1}' dsim-all-chromosome-r1.3.fasta > dsim-short.fasta
awk '{print $1}' dyak-all-chromosome-r1.3.fasta > dyak-short.fasta
```

## Create multiple genome alignments using mauve ##

Start Mauve by double clicking.

In case Mauve complains about insufficient memory you may assign more memory by using (progressiveMauve has to be added to the path):

```
java -Xmx2000m -jar Mauve.jar
```

  * Click File and Align with progressiveMauve
  * Click Add sequence and select dmel-short.fasta
  * Click Add sequence and select dsim-short.fasta
  * Click Add sequence and select dyak-short.fasta
**Note:** the order of the sequences is important. The first should be the reference sequence, i.e the genome for which the divergence will be calculated and the last should be the outgroup (D. yakuba)
  * Click output file and enter: dmel-dsim-dyak
  * Press Align

## Parse the Mauve output ##
The output of Mauve is in the multi fasta alignment format.
To calculate the divergence of D.melanogaster using a sliding window approach it is necessary to reformat the output.
```
perl mauve-parser.pl --ref-input dmel-short.fasta --input mel-sim-yak --output multi-align-dmel-anchored --crosscheck
```

This will create an output like this:
```
2L      10588   A       A       A
2L      10589   T       -       -
2L      10590   T       -       T
2L      10591   T       T       T
2L      10592   T       T       T
```

  * column 1: the chromosome id of _D. melanogaster_
  * column 2: the position in the respective chromosome
  * column 3: the reference character of _D. melanogaster_
  * column 4: the reference character of _D. simulans_
  * column 5: the reference character of _D. yakuba_

## Calculate the divergence using a sliding window approach ##

```
perl calculate-dxy.pl --input multi-align-dmel-anchored  --window-size 50000 --step-size 50000 --output divergence
```

## Reformat for the script Visualize-output.pl ##

```
awk 'BEGIN{OFS="\t"}{print $1,$2,"-","-",$3}' divergence > divergence.reformated
```

## Create a graphical overview ##
```
perl ~/dev/PopGenTools/Visualise-output.pl --input divergence.reformated --output dmel.div5k.pdf --ylab dxy --chromosomes "X 2L 2R 3L 3R 4"
```

The result should be like this:
http://popoolation.googlecode.com/files/dmel.div5k.pdf