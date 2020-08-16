
## Characterizing the between-host transmission of SARS-CoV-2 minor-variants in a linked cluster of infections

-----

### Project Overview

The influence of genetic diversity accumulated within an infected
individual on the global evolutionary dynamics of SARS-CoV-2 is unknown.
Characterizing this dynamic requires an estimate of the amount of
genetic variation transferred between individuals during an infection.
Here, we leverage a unique dataset to measure this bottleneck.
Deep-sequencing of RNA collected from a linked cluster of infections on
an isolated fishing vessel enables the observation of the transmission
of minor variants between individuals. We propose using this data to
provide the highest resolution estimate of the inter-host transmission
bottleneck for SARS-CoV-2 to date.

### Methods

#### Cohort

There was an outbreak on a fishing boat in Seattle. Individuals
(120/122) were screened for SARS-CoV-2 RNA by RT-qPCR and seroconversion
by the Abbot IgG assay before departure and after the return of the
vessel. Although all were negative after the initial testing, an
outbreak on the boat occurred. A total of 104/122 crew members had viral
RNA or seroconverted upon return to shore (84% Attack Rate).

#### Sequencing Data

From the 104 reportedly infected crew members, there are 39 metagenomic
RNA-sequencing runs available. Libraries were sequenced on a 1x75 bp
Illumina NextSeq run. A median of 509,551 sequencing reads were obtained
for each sample. The data is publicly available at NCBI BioProject
[`PRJNA610428`](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=2&WebEnv=MCID_5f344f743a76fab9ab22e42b&o=acc_s%3Aa).

#### Variant Calling

To measure intra-host variation in each sample, the data was processed
through a variant calling pipeline. The raw sequencing reads were
processed using `Fastp` to trim adaptor sequences, low-quality bases,
filter low-quality reads, and trim poly-A tails. Processed reads were
subsequently aligned to a reference (NC\_045512) using either `BWA` or
`STAR`, and variants were called using either `Varscan2` or `Lofreq`.

### Analysis

#### Included Samples

The paper has a supplementary table (Supplementary Table 2.) that
includes the strain names for each individual and the GISAID ID for the
consensus, but not the SRA accession. I used this information to search
the Bioproject
([`PRJNA610428`](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=2&WebEnv=MCID_5f344f743a76fab9ab22e42b&o=acc_s%3Aa))
for the relevant samples. This Bioproject contains all of the samples
sequenced by the UW Virology Department (~1000 Runs). Of the 39 samples
referenced in the paper, there are 28 that had corresponding Runs
deposited in this project and 11 that did not. As a quasi-control, I
randomly sampled an equal number of sequencing runs from the rest of the
UW bioproject that were done on the same sequencer, were also
single-ended, and were within one standard deviation of the mean library
size for the transmission cluster Runs.

#### Quality Contol

Data for the quality control portion of the analysis was collated by the
tool `multiqc`. It consists of statistics produced by `samtools`,
`STAR`, `fastp`, and `fastQC.` I’ve also included the coverage over each
sample in this part of the analyis.

**The percentage of reads that mapped to
SARS-CoV-2.**

<img src="Bottleneck-Analysis_files/figure-gfm/Percentage Mapped-1.png" style="display: block; margin: auto;" />

There is a tremendous amount of variation in the proportion of reads
that correctly map. It’s fairly well distributed between the control and
cluster samples, but on average the cluster samples have a lower mapping
rate. Interestingly, there is a bigger effect of aligner choice (`STAR`
v. `BWA`) on mapping rate in the cluster than the control. I wonder if
this speaks to some underlying issues with sample quality that isn’t
well summarized by mapping rate?

**The mean coverage over the SARS-CoV-2 reference genome for each
sample.**

<img src="Bottleneck-Analysis_files/figure-gfm/Average Coverage-1.png" style="display: block; margin: auto;" />

Coverage is more consistent than the mapping rate between different
aligners. The same relationship holds that the controls have higher
coverage than the cluster samples on average.

**The ratio of transitions to transversions in each
sample.**

<img src="Bottleneck-Analysis_files/figure-gfm/Ts/Tv Ratio-1.png" style="display: block; margin: auto;" />

In theory, transitions are much more common than transversions
biologically, however Illumina sequencing errors tend to bais towards
transversions. Above, I plotted the ratio transitions (Ts) to
transversions (Tv) for each sample and colored them by the experiment
(cluster v. control). The cluster samples seem to have a higher ratio of
transitions. I’m not sure what to make of this yet. This is just for BWA
and Lofreq. From here on, I’ll just use that aligner-caller combo.

**The concordance between the pipeline for annotated
SNPs.**

<img src="Bottleneck-Analysis_files/figure-gfm/Pipeline Comparison, -1.png" style="display: block; margin: auto;" /><img src="Bottleneck-Analysis_files/figure-gfm/Pipeline Comparison, -2.png" style="display: block; margin: auto;" />

The difference in concordance between all variants (first plot) and just
variants that have an allele frequency above 2% (second plot) is huge.
It’s worth taking this into consideration in down stream analyses and
perhaps filtering out alleles with frequencies less than 2%.

#### Mutational Spectra

**The ratio of Synonymous/Non-Synonymous/Missense mutations - all
variants and only minor
variants.**

<img src="Bottleneck-Analysis_files/figure-gfm/Effect Distribution-1.png" style="display: block; margin: auto;" />

The ratio of mutational effects is very similar between the four
combinations and consistent with other studies that I’ve seen.

**The distribution of nucleotide changes - all variants and only minor
variants.**

<img src="Bottleneck-Analysis_files/figure-gfm/Mutational Spectra-1.png" style="display: block; margin: auto;" /><img src="Bottleneck-Analysis_files/figure-gfm/Mutational Spectra-2.png" style="display: block; margin: auto;" />

Similar to the Ts/Tv plot above, I wanted to see the mutational spectra
for minor and major variants. Overall, there seems to be a slight bias
towards C\>T transtions, which is also consistent with other studies,
but it’s not as apparent in only the minor variants.

**The distribution of nucleotide changes across genes, normalized to
size - all variants and only minor
variants.**

<img src="Bottleneck-Analysis_files/figure-gfm/Gene Distribution-1.png" style="display: block; margin: auto;" />

#### Shared Variation

**Dirstribution of minor
variants.**

<img src="Bottleneck-Analysis_files/figure-gfm/Minor SNP Distribution-1.png" style="display: block; margin: auto;" />

**Upset plot of shared Minor
SNPs.**

<img src="Bottleneck-Analysis_files/figure-gfm/Minor SNP Uspet-1.png" style="display: block; margin: auto;" />

\*\* Correlation between frequency of overlap and allele frequency in
minor SNPs.
\*\*

<img src="Bottleneck-Analysis_files/figure-gfm/Correlation between AF and Overlap-1.png" style="display: block; margin: auto;" />

    ##      [,1] [,2]
    ## [1,]   71  217
    ## [2,]  127  593

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  
    ## p-value = 0.01387
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  1.080774 2.147583
    ## sample estimates:
    ## odds ratio 
    ##   1.527071

![](Bottleneck-Analysis_files/figure-gfm/Pairwise%20Comparison-1.png)<!-- -->

    ## `summarise()` ungrouping output (override with `.groups` argument)

    ## # A tibble: 2 x 3
    ##   Experiment  mean median
    ##   <chr>      <dbl>  <dbl>
    ## 1 Cluster    11.7      10
    ## 2 Control     9.02      7
