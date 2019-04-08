# Athena
#### A toolkit for exploring structural variation mutation rates and dosage sensitivity

Copyright (c) 2019, [Ryan L. Collins](mailto:rlcollins@g.harvard.edu) and the Talkowski Laboratory.  
Distributed under terms of the [MIT License](/LICENSE) (see `LICENSE`).  

---  

## Table of Contents

#### Getting started
  * [Run Athena with Docker](https://github.com/talkowski-lab/athena#run-from-docker)
  * [Manual Installation](https://github.com/talkowski-lab/athena#manual-installation)
  * [Invoking Athena](https://github.com/talkowski-lab/athena#invoking-athena)

#### The Athena workflow
_Mutation rate modeling_ 
  1. [Filter SVs for mutation rate training](https://github.com/talkowski-lab/athena#step-1)  
  2. [Calculate SV size & spacing distributions](https://github.com/talkowski-lab/athena#step-2) (Optional)   
  3. [Create training bins](https://github.com/talkowski-lab/athena#step-3)   
  4. [Annotate training bins](https://github.com/talkowski-lab/athena#step-4)  
  5. [Collapse redundant annotations](https://github.com/talkowski-lab/athena#step-5)  

_Dosage sensitivity modeling_  

#### Other
  * [A note on design](https://github.com/talkowski-lab/athena#a-note-on-design)
  * [About the name](https://github.com/talkowski-lab/athena#about-the-name)


---  

## Run from Docker

The recommended way to run Athena is from its dedicated Docker container. This will handle all dependencies and installation for you, and ensure you are running the latest version.

```
$ docker pull talkowski/athena
$ docker run --rm -it talkowski/athena
```

## Manual installation

If you would prefer to install Athena on your own system, you can do so with `pip`.

```
$ git clone https://github.com/talkowski-lab/athena.git
$ cd athena
$ pip install -e .
```

## Invoking Athena

Athena is called from the command line:
```
$ athena --help
Usage: athena [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  annotate    Annotate bins
  count-sv    Intersect SV and bins
  decomp      Decompose bin annotations
  make-bins   Create bins
  query       Mutation rate lookup
  vcf-filter  Filter an input VCF
  vcf-stats   Get SV size & spacing
```

Athena has numerous subcommands. Specify `--help` with any subcommand to see a list of options available. For example:
```
$ athena vcf-filter --help
Usage: athena vcf-filter [OPTIONS] VCF OUT

  Filter an input VCF

Options:
  --include-chroms TEXT  Chromosomes to include (comma-separated) [default:
                         include all chromosomes]
  --exclude-chroms TEXT  Chromosomes to exclude (comma-separated) [default:
                         exclude no chromosomes]
  --svtypes TEXT         SV classes to include (comma-separated) [default: all
                         SVs]
  -x, --blacklist TEXT   BED file of regions to exclude, based on SV overlap
                         [default: None]
  --minAF FLOAT          Minimum allowed allele frequency [default: 0]
  --maxAF FLOAT          Maximum allowed allele frequency [default: 1.0]
  --minAC INTEGER        Minimum allowed allele count [default: 0]
  --maxAC INTEGER        Maximum allowed allele count [default: None]
  --vcf-filters TEXT     VCF FILTER statuses to be included [default: PASS]
  --minQUAL INTEGER      Minimum allowed QUAL score [default: 0]
  --maxQUAL INTEGER      Maximum allowed QUAL score [default: None]
  --pHWE FLOAT           Minimum Hardy-Weinberg equilibrium p-value to include
                         [default: Do not filter on HWE]
  --keep-infos TEXT      INFO fields to retain in output VCF (comma-
                         separated). Will always retain the following: END,
                         CHR2, SVTYPE, SVLEN. Specifying "ALL" will keep all
                         INFO. Note: this does _not_ filter records based on
                         INFO. [default: Keep no other INFO]
  -z, --bgzip            Compress output with bgzip
  --help                 Show this message and exit.
```

--- 

# The Athena workflow

### Step 1
**Filter SVs for mutation rate training**  
First, SVs must be filtered to an appropriate subset to train a mutation rate model.  

This is accomplished with `athena vcf-filter`, which offers many options to filter an input VCF. These options are enumerated in the `--help` text.  

Conceptually, this filtering step aims to accomplish a few points:
1. Restrict to rare variants, as these will have more recently arisen in the population  
2. Restrict to a subset of near-neutral variants, as these will limit the degree to which selection confounds mutation rate estimates  
3. Restrict to high-quality variants to reduce the influence of purely technical factors in mutation rate estimates  

For instance, to generate a training set of deletions from the [gnomAD-SV v2 dataset](https://gnomad.broadinstitute.org/downloads) (`gnomad_v2_sv.sites.vcf.gz`) for a deletion mutation rate model, you could run the following:  
```
$ athena vcf-filter -z \
  	--exclude-chroms X,Y \
  	--svtypes DEL \
  	--blacklist data/athena.SV_selection_blacklist.v1.GRCh37.bed.gz \
  	--maxAF 0.01 \
  	--minAC 1 \
  	--minQUAL 100 \
  	--pHWE 0.01 \
  	gnomad_v2_sv.sites.vcf.gz \
  	athena_training_deletions.vcf.gz
```
The above command will return a set of 129,033 rare, high-quality autosomal deletions from gnomAD-SV for training a deletion mutation rate model downstream.  

See `data/` for a description of the file provided to `--blacklist`.  

--- 

### Step 2 (Optional)
**Calculate SV size & spacing distributions**
Selecting an appropriate resolution for the mutation rate model is important, as it determines the accuracy of your results and the file sizes of the model.  

To aid in this process, you can compute SV size & spacing distributions with `athena vcf-stats`.  

For example, you can invoke `athena vcf-stats` on the training set [generated above](https://github.com/talkowski-lab/athena#step-1) as follows:
```
$ athena vcf-stats athena_training_deletions.vcf.gz
```

This will print two distributions directly to the console:
```
SV spacing quantiles:
---------------------
25.0%: 0 bp
50.0%: 3,138 bp
75.0%: 10,951 bp
90.0%: 22,031 bp
95.0%: 30,908 bp
99.0%: 52,568 bp
99.9%: 89,697 bp

SV size quantiles:
------------------
25.0%: 139 bp
50.0%: 642 bp
75.0%: 2,843 bp
90.0%: 7,303 bp
95.0%: 12,000 bp
99.0%: 41,716 bp
99.9%: 150,958 bp
```

From these data, we can see that about 99% of all deletions in the training set are smaller than 42kb, which means that we probably do not need to extend the training of our 2D model beyond 40kb to capture almost all of the informative 2D signal.  

We can also observe that roughly half of all deletions are no farther than a few kilobases away from the next-nearest deletion, which means that bin sizes of several kilobases should result in reasonably dense training data.  

Based on this logic along with some approximate rounding, we can decide on the two key parameters used in the remaining steps of the mutation rate model:
```
1D bin size: 4kb
2D maximum bin-pair distance: 40kb
```

--- 

### Step 3
**Create training bins**  
The next step is to segment the genome into sequential, uniform bins.  

These bins will only be used for training the mutation rate model. Once trained, you can extend the mutation rate model beyond this subset of training bins.  

Following [the example from step 2](https://github.com/talkowski-lab/athena#step-2) (above), we can generate 4kb bins for all autosomes using `athena make-bins`.  

We can also introduce any number of `--blacklist` BED files, which will exclude bins overlapping the blacklist(s) from being included in model training. We can also apply a buffer around the elements in the blacklists. For this example, we will exclude all bins within ±4kb of the intervals used for SV blacklisting ([see step 1](https://github.com/talkowski-lab/athena#step-1), above) as well as all bins overlapping any unalignable (i.e., N-masked) regions of the reference genome.  

We can create these bins as follows:
```
$ athena make-bins -z \
  	-x data/athena.SV_selection_blacklist.v1.GRCh37.bed.gz \
  	-x GRCh37.Nmask.bed.gz \
  	--buffer 4000 \
  	--exclude-chroms X,Y,M \
  	GRCh37.genome \
  	4000 \
  	athena_training_bins.bed.gz
```

Where:
  * `GRCh37.Nmask.bed.gz` is a BED file containing all N-masked regions of the reference genome. This file is available from numerous sources, including [UCSC](https://genome.ucsc.edu).  
  * `GRCh37.genome` is a two-column, tab-delimited file [per BEDTools spec](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html).  

--- 

### Step 4
**Annotate training bins**  
Once the genome has been segmented into sequential, uniform bins (see the [example in step 3](https://github.com/talkowski-lab/athena#step-3)), we next must annotate these bins with any features to be considered in mutation rate modeling.  

Athena has a flexible interface to apply multiple annotations directly from various sources, invoked as `athena annotate-bins`.  

Currently, supported data sources are as follows:  

| Track format | Source | Athena flag | Available actions |  
| :--- | :--- | :--- | :--- |   
| Any file compatible with BEDTools [`coverage`](https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html) or [`intersect`](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) (e.g., BED, VCF, GFF, etc.) | Local | `--track` | `--count`, `--count-unique`, `--coverage` |  
| [BigWig](https://genome.ucsc.edu/goldenPath/help/bigWig.html) | Local *or remote<sup>1</sup>* | `--track` | `--map-mean`, `--map-sum`, `--map-min`, `--map-max` |  
| [UCSC Genome Browser Tables](https://genome.ucsc.edu/cgi-bin/hgTables)<sup>2</sup> | Hosted by UCSC | `--ucsc-track` | `--count`, `--count-unique`, `--coverage` |  
| [UCSC-Hosted BigWig Tracks](https://genome.ucsc.edu/cgi-bin/hgTables)<sup>2</sup> | Hosted by UCSC | `--ucsc-track` | `--map-mean`, `--map-sum`, `--map-min`, `--map-max` |  
| [FASTA](https://zhanglab.ccmb.med.umich.edu/FASTA/) | Local | `--fasta` | _None required<sup>3</sup>_ |  


#### Notes:
 1. If specifying a remote BigWig file, pass the full URL of the remote-hosted file to `--track`.  
 2. If any UCSC tracks are requested, the UCSC reference build must also be specified with `--ucsc-ref`. Athena will automatically handle necessary conversions between GRC & UCSC contig nomenclature (e.g., `chr1` vs `1`).  
 3. Specifying `--fasta` adds a `pct_gc` annotation, which is the GC content of each bin.  

The precise annotations added at this stage are up to the user.  

For example, we could annotate the bins from [step 3 (above)](https://github.com/talkowski-lab/athena#step-3) with the following six tracks:
 1. Counts per bin vs. a custom local annotation file (`my_local_annotation.bed`)
 2. Average ovary expression level per bin from ENCODE (accession [ENCFF250KXQ](https://www.encodeproject.org/experiments/ENCFF250KXQ/))
 3. Maximum ovary chromatin accessibility score per bin from ENCODE (accession [ENCFF416KSV](https://www.encodeproject.org/experiments/ENCFF416KSV/))  
 4. Coverage per bin by segmental duplications from UCSC (table `genomicSuperDups`)  
 5. Mean mapability of 100mers from UCSC (table `wgEncodeCrgMapabilityAlign24mer`, which links to a remote bigWig file)  
 6. GC content from the GRCh37 reference (`GRCh37.fa`).  

We would add these six annotations to the bins from step 3 as follows:
```
$ athena annotate-bins -z \
    -t my_local_annotation.bed -a count -n my_annotation \
    -t https://www.encodeproject.org/files/ENCFF250KXQ/@@download/ENCFF250KXQ.bigWig -a map-mean -n ovary_sense_strand_expression \
    -t https://www.encodeproject.org/files/ENCFF416KSV/@@download/ENCFF416KSV.bigWig -a map-max -n ovary_peak_dnase_hypersensitivity \
    -u genomicSuperDups -a coverage -n segdup_coverage \
    -u wgEncodeCrgMapabilityAlign100mer -a map-mean -n mapability_100mers \
    --ucsc-ref hg19 \
    --fasta GRCh37.fa \
    -z \
    athena_training_bins.bed.gz \
    athena_training_bins.annotated.bed.gz
```

Don't worry about correlations between various annotation tracks -- this correlation structure will be handled in [step 5, below](https://github.com/talkowski-lab/athena#step-5).  

--- 

### Step 5  
**Collapse redundant annotations**  
Many genomic anntations are redudant (e.g., microsatellites and low sequence uniqueness, or GC content and protein-coding exons). 

To control for the underlying correlation structure of the annotations included in [step 4 (above)](https://github.com/talkowski-lab/athena#step-4), the next step in the Athena workflow is to perform Eigenvector decomposition of the annotation matrix.  

This can be accomplished with `athena eigen-bins`. For example, we could decompose the annotations from [the example in step 4]((https://github.com/talkowski-lab/athena#step-4) into the top three Eigenfeatures as follows:  
```
athena eigen-bins \
  --eigenfeatures 3 \
  --stats eigenfeature_stats.txt \
  -z \
  athena_training_bins.annotated.bed.gz \
  athena_training_bins.annotated.decomped.bed.gz
```  

If you include the `--stats` option, Athena will also generate a log file with the summary of all Eigenfeatures, including total variance explained by each, and the Spearman correlation of raw annotation tracks significantly correlated with each Eigenfeature. This log file can be useful for determining _what_ each Eigenfeature represents (roughly speaking).  

---  

### A note on design
This package was designed with [the gnomAD-SV callset](https://gnomad.broadinstitute.org/downloads) in mind. To that end, it assumes input data follows gnomAD-SV formatting standards. This may cause issues for alternative styles of SV representation, or different metadata labels. If using non-gnomAD data with Athena, please compare your VCF formatting standards, and the `INFO` field in particular.  

You can read more about the gnomAD-SV dataset [in the corresponding preprint](https://broad.io/gnomad_sv).

### About the name
This package is named after [Athena](https://en.wikipedia.org/wiki/Athena), the Greek goddess of wisdom, strategy, tactics, and mathematics. She was selected as the namesake for this package given that it relies on understanding the features that influence structural variation mutation rates (_wisdom_), incorporating those features into a statistical model (_mathematics_), and using these models to infer which components of the genome are vulnerable to changes in copy number (a kind of genomic _tactics_/_strategy_).