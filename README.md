snpEffect
=========

Determine the effect of a SNP, given GFF annotation.  An `R` function which uses [BioConductor](http://www.bioconductor.org).

The main functions are `snpEffect()` and `indelEffect()`.  These take a `data.frame` of SNPs/indels, a collection of GFF annotation as formatted by BioConductor's `import.gff()`, and a genome sequence as formatted by BioConductor's `read.DNAStringSet()`.  These latter two functions can be initiated with names of standard GFF and Fasta files.  The SNPs on the other hand are in a non-standard format, but it shouldn't be too difficult to get them into the format required.  If I find the time I would love to accept a BioConductor object filled from a VCF file for the SNPs.

I originally wrote this to save me a bunch of time because I couldn't get [SnpEff](http://snpeff.sourceforge.net) to work for my set of SNPs.  As a result, it is strongly skewed toward the original formulation of my problem and still has a number of blind spots.  For example, it doesn't always do a good job of dealing with the '-' strand in coding regions.

### Input and output formats

SNPs are described in a simple `data.frame` format, and SNP effects are returned as an annotated version of the input `data.frame`.

### Example

````R
source("snpEffect.R")
snp = read.delim("snps.txt")  # dataframe format to be described
genome = read.DNAStringSet("genome.fa")
gff = import.gff3("annotation.gff")
snp.effects = snpEffect(snp, gff, genome)
````
