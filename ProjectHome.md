


# Introduction #

Currently next generation sequencing (NGS) technologies are mainly used to sequence individuals. However, the high coverage required and the resulting costs may be
prohibitive for population scale studies. Sequencing pools of individuals instead may often be more cost effective and more accurate than sequencing individuals. PoPoolation is a pipeline for analysing pooled next generation sequencing data. PoPoolation builds upon open source tools (bwa, samtools) and uses standard file formats (gtf, sam, pileup) to ensure a wide compatibility. Currently PoPoolation allows to calculate Tajima’s Pi, Watterson’s Theta and Tajima’s D for reference sequences using a sliding window approach. Alternatively these population genetic estimators may be calculated for a set of genes (provided as gtf). One of the main challenges in population genomics is to identify regions of intererest on a genome wide scale. We believe that PoPoolation will greatly aid this task by allowing a fast and user friendly analysis of NGS data from DNA pools.

**Note a major novel version has been released: ReleaseNotes**

# Usage #

  * Quick walkthrough: TeachingPoPoolation
  * Detailed walkthrough: [PoPOOLationWalkthrough](PoPOOLationWalkthrough.md)
  * Manual: [Manual](Manual.md)
  * Slides from a course http://drrobertkofler.wikispaces.com/PoPoolationGenomics

# How to cite PoPoolation #

Please cite the following two paper

  * Kofler R, et al. (2011) PoPoolation: A Toolbox for Population Genetic Analysis of Next Generation Sequencing Data from Pooled Individuals. PLoS ONE http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0015925

You may also be interested in our Pool-seq review (Nature Reviews Genetics) where we provide some recommendations for the analysis of Pool-seq data:

  * http://www.nature.com/nrg/journal/vaop/ncurrent/abs/nrg3803.html



# Partner projects #

**Gowinda**: unbiased analysis of gene set enrichement (e.g: Gene Ontology) for Genome Wide Association Studies. Gowinda may thus be used for biological interpretation of the results of PoPoolation and PoPoolation2:  http://code.google.com/p/gowinda/

**PoPoolation2**: Allows analyzing the population frequencies of SNPs from two or more populations. It may be used to identify differentiation between populations or to analyze data from genome wide association studies.
http://code.google.com/p/popoolation2/

**PoPoolation TE**: A quick and simple pipeline for the analysis of transposable element insertion frequencies in populations from pooled next generation sequencing data. PoPoolation TE identifies TE insertions that are present in the reference genome as well as novel TE insertions and estimates their population frequencies. This also allows for an comparision of TE insertion frequencies between different populations http://code.google.com/p/popoolationte.

**PoPoolation DB**: A user-friendly web-based database for the retrieval of natural variation in _Drosophila melanogaster_
http://www.popoolation.at/pgt/

# Authors #

  * Robert Kofler http://drrobertkofler.wikispaces.com/
  * Christian Schlötterer http://i122server.vu-wien.ac.at/pop/lab_members/christian_schloetterer.html


