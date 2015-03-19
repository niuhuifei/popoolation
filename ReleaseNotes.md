# 1.2.2 #

Calculation of Tajima's D has been updated (equation: http://popoolation.googlecode.com/files/correction_equations.pdf)
  * Novel requirement for corrected Tajima's D: Poolsize > 3 x minimum coverage
  * Novel requirement for corrected Tajima's D: minimum count must be 2!
  * Novel recommendation for corrected Tajima's D. We recommend subsampling of the pileup to a uniform coverage when calculating Tajima's D
  * nbase=expected\_value(minimum coverage) => conservative (for details see below)

**Details about calculating Tajima's D**

As a basis for the equation for calculating Tajima's D we are using the results from Achaz (2008; http://www.ncbi.nlm.nih.gov/pubmed/18562660) for individual sequencing when singletons are omitted. The difficulty arises when calculating the variance of Tajima's D.
The following simplifying assumptions are made in our variance approximation:
  * The locally fluctuating coverage is replaced by the minimum coverage. This makes the variance estimator larger, and therefore leads to conservative estimates of Tajima’s D.
  * The random number of different individuals sequenced (nbase) under a given coverage is replaced by its expected value. This assumption should not affect the results much, if the pool size is large compared to the coverage, sequencing the same individual more than once is uncommon under these circumstances. Furthermore the number of different individuals sequenced will have a low variance.  As we are working with the minimum coverage, nbase will be biased downwards—tending to give a conservative estimate oft he variance.
  * Different subsets from the pool are sequenced at different positions. Their coalescent histories will be correlated but not identical. As the classical formula for Tajima’s D are for a single sample sharing a common coalescent history, there is more independence in the data than assumed with the classical formula—which should again make the variance approximation more conservative.
Summing up, the approximate variance in our equations provides a conservative approximation, and the values for Tajima’s D will tend to be smaller than those that would be expected for an experiment based on individual sequencing of a single sample.

For a correlation between the classical Tajima's D and the corrected Tajima's D using the equation described above please see (x: classical, y: corrected). Correlation has been done with real data from _Drosophila_ with a coverage of 12, a window size of 500 and a minimum count of 1) ![http://popoolation.googlecode.com/files/correlation_classic_correctedTajimasD.png](http://popoolation.googlecode.com/files/correlation_classic_correctedTajimasD.png)


# 1.2.1 #
This is a major novel release of PoPoolation which has several new features and a major bug fix. Most notably there was a problem with calculating Tajima's D.

  * Fixed a bug in calculating Tajima's D. For this novel equations have been developed to account for the idiosyncracies of pooling (Pool size, minor allele count). This novel method is aiming to closely reproduce the classical Tajima's D values. The novel equations can be found here: http://popoolation.googlecode.com/files/correction_equations.pdf
  * PoPoolation now allows to calculate the classical Tajima's D, Tajima's Pi, and Watterson's Theta (use `--dissable-corrections`). Since no corrections can be applied to the classical measures a minimum allele count of 1 has to be used. This feature allows to explore the effect of the correction factors, by comparing corrected and uncorrected measures.
  * a new script has been added that allows to randomly sub-sample bases of a pileup file to achieve a uniform coverage throughout the whole genome. This may be useful to account for coverage effects. This is highly recommended when for example calculating the classical Tajima's D
  * As a major novel feature PoPoolation now allows to calculate the synonymous and non-synonymous Pi, Theta, Tajima's D, both corrected and uncorrected. The synonymous and non-synonymous measures may be calculated using a sliding window approach or for every feature (gene, exon) separately.