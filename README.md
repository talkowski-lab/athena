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
  1. [Filter SVs for mutation rate training](https://github.com/talkowski-lab/athena#step-1-filter-svs-for-mutation-rate-training)  
  2. [Calculate SV size & spacing distributions](https://github.com/talkowski-lab/athena#step-2-optional-calculate-sv-size--spacing-distributions) (_optional_)   
  3. [Create 1D training bins](https://github.com/talkowski-lab/athena#step-3-create-1d-bins-genome-wide)   
  4. [Annotate 1D training bins](https://github.com/talkowski-lab/athena#step-4-annotate-1d-training-bins)  
  5. [Collapse redundant 1D annotations](https://github.com/talkowski-lab/athena#step-5-collapse-correlated-1d-annotations)  

_Dosage sensitivity modeling_  

#### Other
  * [A note on design](https://github.com/talkowski-lab/athena#a-note-on-design)
  * [About the name](https://github.com/talkowski-lab/athena#about-the-name)


---  

## Run from Docker

The recommended way to run Athena is from its dedicated Docker container. This will handle all dependencies and installation for you, and ensure you are running the latest version.

```
$ docker pull rlcollins/athena
$ docker run --rm -it rlcollins/athena
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
  annotate-bins  Annotate bins
  count-sv       Intersect SV and bins
  eigen-bins     Eigendecomposition of annotations
  feature-hists  Plot bin annotation distributions
  make-bins      Create sequential bins
  query          Mutation rate lookup
  vcf-filter     Filter an input VCF
  vcf-stats      Get SV size & spacing
```

Athena has numerous subcommands. Specify `--help` with any subcommand to see a list of options available.  

--- 

# The Athena workflow

### Step 1: Filter SVs for mutation rate training  

First, SVs must be filtered to an appropriate subset to train a mutation rate model.  

This is accomplished with `athena vcf-filter`, which offers many options to filter an input VCF. These options are enumerated in the `--help` text.  

Conceptually, this filtering step aims to accomplish a few points:
1. Restrict to rare variants, as these will have more recently arisen in the population  
2. Restrict to a subset of near-neutral variants, as these will limit the degree to which selection confounds mutation rate estimates  
3. Restrict to high-quality variants to reduce the influence of purely technical factors in mutation rate estimates  

For instance, to generate a training set of deletions from the [gnomAD-SV v2.1 dataset](https://gnomad.broadinstitute.org/downloads) (`gnomad_v2.1_sv.sites.vcf.gz`) for a deletion mutation rate model, you could run the following:  

```
# First, download the gnomAD-SV v2.1 VCF
$ wget https://storage.googleapis.com/gnomad-public/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz
$ wget https://storage.googleapis.com/gnomad-public/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz.tbi

# Second, filter the gnomad-SV VCF with Athena
$ athena vcf-filter -z \
    --exclude-chroms X,Y,M \
    --svtypes DEL \
    --blacklist data/athena.SV_selection_blacklist.v1.GRCh37.bed.gz \
    --maxAF 0.001 \
    --minAC 1 \
    --minAN 20402 \
    --minQUAL 100 \
    --pHWE 0.01 \
    gnomad_v2.1_sv.sites.vcf.gz \
    athena_training_deletions.vcf.gz
```

The above command will return a set of 96,464 rare, high-quality autosomal deletions from gnomAD-SV for training a deletion mutation rate model downstream.  

See `data/` for a description of the file provided to `--blacklist`.  

--- 

### Step 2 (Optional): Calculate SV size & spacing distributions  

Selecting an appropriate resolution for the mutation rate model is important, as it determines the accuracy of your results and the file sizes of the model.  

To aid in this process, you can compute SV size & spacing distributions with `athena vcf-stats`.  

For example, you can invoke `athena vcf-stats` on the training set [generated above](https://github.com/talkowski-lab/athena#step-1-filter-svs-for-mutation-rate-training) as follows:
```
$ athena vcf-stats athena_training_deletions.vcf.gz
```

This will print two distributions directly to the console:
```
SV spacing quantiles:
---------------------
25.0%: 143 bp
50.0%: 5,556 bp
75.0%: 15,561 bp
90.0%: 29,410 bp
95.0%: 40,475 bp
99.0%: 67,920 bp
99.9%: 114,766 bp

SV size quantiles:
------------------
25.0%: 117 bp
50.0%: 631 bp
75.0%: 3,132 bp
90.0%: 8,060 bp
95.0%: 13,127 bp
99.0%: 46,221 bp
99.9%: 156,913 bp
```

From the first distribution, `SV spacing`, we can observe that over half of all deletions are no farther than 6 kilobases away from the next-nearest deletion, which means that bin sizes of several kilobases should result in reasonably dense training data.  

And based on the second distribution, `SV size`, we can also see that 99% of all deletions in the training set are smaller than 47kb, which means that we probably do not need to extend the training of our 2D model beyond \~50kb to capture almost all of the informative 2D signal.  

Likewise, we can see that 99.9% of deletions are smaller than 156kb, so we can safely restrict computing 2D mutation rates to a maximum distance of 160kb while only failing to model <0.1% of all deletions.

Based on this logic along with some approximate rounding, we can decide on the three key parameters used in the remaining steps of the mutation rate model:
```
1D bin size: 5kb
2D maximum bin-pair distance (training): 50kb
2D maximum bin-pair distance (full model): 160kb
```

--- 

### Step 3: Create 1D bins genome-wide  

The next step is to segment the genome into sequential, uniform bins.  

Using the bin size determined in [the example from step 2](https://github.com/talkowski-lab/athena#step-2-optional-calculate-sv-size--spacing-distributions) (above), we can generate sequential 5kb bins for all autosomes using `athena make-bins`.  

By default, sequential bins will be created from telomere to telomere for every chromosome. However, we can introduce any number of `-x`/`--blacklist-all` BED files to exclude bins overlapping the blacklist(s).  In this example, we will exclude bins overlapping unalignable (_i.e._, N-masked) regions.  

We can simultaneously define a smaller subset of bins for training the mutation rate model. Once trained, you can extend the mutation rate model beyond these training bins. 

To generate a subset of training bins, we can pass any number of `--training-blacklist` BED files, which will apply further blacklisting to the bins generated above. In this example, we will exclude training bins overlapping the intervals used for SV blacklisting ([see step 1](https://github.com/talkowski-lab/athena#step-1-filter-svs-for-mutation-rate-training), above).  

Finally, we have the option to apply a buffer around the elements in the blacklists. For this example, we will conservatively exclude all bins within ±5kb of any interval from the blacklists.  

We can create these bins as follows:  
```
$ athena make-bins -z \
    --blacklist-all data/GRCh37.Nmask.bed.gz \
    --blacklist-training data/athena.SV_selection_blacklist.v1.GRCh37.bed.gz \
    --buffer 5000 \
    --exclude-chroms X,Y,M \
    --training-bins-out GRCh37.autosomes.5kb_bins.training.bed.gz \
    data/GRCh37.genome \
    5000 \
    GRCh37.autosomes.5kb_bins.all.bed.gz
```

Where:
  * `GRCh37.Nmask.bed.gz` is a BED file containing all N-masked regions of the reference genome. This file is available from numerous sources, including [UCSC](https://genome.ucsc.edu). A copy for GRCh37 is [packaged with Athena](https://github.com/talkowski-lab/athena/tree/master/data)  
  * `GRCh37.genome` is a two-column, tab-delimited file [per the BEDTools specifications](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html). A copy for GRCh37 is [packaged with Athena](https://github.com/talkowski-lab/athena/tree/master/data).    

--- 

### Step 4: Annotate 1D training bins  

Once the genome has been segmented into sequential, uniform bins (see the [example in step 3](https://github.com/talkowski-lab/athena#step-3-create-1d-bins-genome-wide)), we next must annotate these bins with any features to be considered in mutation rate modeling.  

Athena has a flexible interface to apply multiple annotations directly from various sources, invoked as `athena annotate-bins`.  

Supported data sources are listed in below:  

| Track format | Source | Athena flag | Available actions |  
| :--- | :--- | :--- | :--- |   
| Any file compatible with BEDTools [`coverage`](https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html) or [`intersect`](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) (_e.g._, BED, VCF, GFF, etc.) | Local | `--track` | `--count`, `--count-unique`, `--coverage`, `--any-overlap` |  
| [BigWig](https://genome.ucsc.edu/goldenPath/help/bigWig.html) | Local *or remote<sup>1</sup>* | `--track` | `--map-mean`, `--map-sum`, `--map-min`, `--map-max` |  
| [UCSC Genome Browser Tables](https://genome.ucsc.edu/cgi-bin/hgTables)<sup>2</sup> | Hosted by UCSC | `--ucsc-track` | `--count`, `--count-unique`, `--coverage`, `--any-overlap` |  
| [UCSC-Hosted BigWig Tracks](https://genome.ucsc.edu/cgi-bin/hgTables)<sup>2</sup> | Hosted by UCSC | `--ucsc-track` | `--map-mean`, `--map-sum`, `--map-min`, `--map-max` |  
| [FASTA](https://zhanglab.ccmb.med.umich.edu/FASTA/) | Local | `--fasta` | _None required<sup>3</sup>_ |  

#### Notes:
 1. If specifying a remote BigWig file, pass the full URL of the remote-hosted file to `--track`.  
 2. If any UCSC tracks are requested, the UCSC reference build must also be specified with `--ucsc-ref`. Athena will automatically handle necessary conversions between GRC & UCSC contig nomenclature (_e.g._, `chr1` vs `1`).  
 3. Specifying `--fasta` adds a `pct_gc` annotation, which is the GC content of each bin. Further passing a trinucleotide SNV mutation rate lookup table as `--snv-mutrate` will add a `snv_mu` annotation, with is the cumulative SNV mutation rate for the DNA sequence of each bin. An example trinucleotide SNV mutation rate lookup table can be found in `data/snv_mutation_rates.Samocha_2014.tsv.gz`.  

The exact annotations added at this stage are for the user to decide.  

For example, we could annotate the bins from [step 3 (above)](https://github.com/talkowski-lab/athena#step-3-create-1d-bins-genome-wide) with the following six tracks:
 1. Counts per bin vs. a custom local annotation file (`my_local_annotation.bed`)
 2. Average ovary expression level per bin from ENCODE (accession [ENCFF250KXQ](https://www.encodeproject.org/experiments/ENCFF250KXQ/))
 3. Maximum ovary chromatin accessibility score per bin from ENCODE (accession [ENCFF416KSV](https://www.encodeproject.org/experiments/ENCFF416KSV/))  
 4. Coverage per bin by segmental duplications from UCSC (table `genomicSuperDups`)  
 5. Mean mapability of 100mers from UCSC (table `wgEncodeCrgMapabilityAlign24mer`, which links to a remote bigWig file)  
 6. GC content from the GRCh37 reference ([`Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz`](https://grch37.ensembl.org/Homo_sapiens/Info/Index)).  

We would add these six annotations to the bins from step 3 as follows:
```
# First, download and unzip the reference assembly fasta file 
# (note: BEDTools currently does not support compressed reference fasta files)
$ wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
$ gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
$ samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa

# Second, run athena annotate-bins
$ athena annotate-bins -z \
    -t my_local_annotation.bed -a count -n my_annotation \
    -t https://www.encodeproject.org/files/ENCFF250KXQ/@@download/ENCFF250KXQ.bigWig -a map-mean -n ovary_sense_strand_expression \
    -t https://www.encodeproject.org/files/ENCFF416KSV/@@download/ENCFF416KSV.bigWig -a map-max -n ovary_peak_dnase_hypersensitivity \
    -u genomicSuperDups -a coverage -n segdup_coverage \
    -u wgEncodeCrgMapabilityAlign100mer -a map-mean -n mapability_100mers \
    --ucsc-ref hg19 \
    --fasta Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    GRCh37.autosomes.5kb_bins.all.bed.gz \
    GRCh37.autosomes.5kb_bins.all.annotated.bed.gz
```

Note that we need to annotate _all_ bins (rather than just the training bins) in order to expand the trained mutation rate model genome-wide.  

Don't worry about correlations between various annotation tracks, either. This correlation structure will be handled in [step 5, below](https://github.com/talkowski-lab/athena#step-5-collapse-correlated-1d-annotations).  

#### Customizing UCSC Genome Browser track queries:  
Athena supports a limited range of conditional filtering and column manipulation for UCSC tracks.  

By default, Athena will extract the minimal BED coordinate information for each track (`chrom`, `start`, `end`).  

Additional columns or column filters can be specified by adding `:` after the track name (`-u`/`--ucsc-track`).  

The syntax for these options follows a set of conventions:
  * Track name must be specified only once, and directly after `-u`/`--ucsc-track`  
  * *All* desired columns must be listed as a column-delimited list after the track name, separated by `:`  
  * Column filtering must be specified directly after the desired column using one of a limited set of operators. Valid operators include `=`, `!=`, `<`, `>`, `<=`, and `>=`.  

If fewer than three columns are specified, `chrom`, `chromStart`, and `chromEnd` will always be appended.  

If no filtering is specified for a given column, no filters will be applied.  

The exact syntax will vary based on columns available for each track, which can be checked using the [UCSC Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables).  

**Warning #1**  
Be careful to wrap modified UCSC track names in quotes, because some operators (especially `>`) can be interpreted as shell flow control characters.

**Example Use Case #1**  
We could restrict the segmental duplications annotation from the example above to only those with high sequence identity as follows:
```
$ athena annotate-bins -z \
    ... 
    -u "genomicSuperDups:fracMatch>0.99" -a coverage -n segdup_coverage_high_identity \
    ...
```

**Example Use Case #2**  
We could calculate the max sequence identity among segmental duplications overlapping each bin as follows:
```
$ athena annotate-bins -z \
    ... 
    -u "genomicSuperDups:fracMatch" -a map-max -n max_segdup_identity \
    ...
```

**Warning #2 & Example Use Cases #3-4**  
It is recommended to check the column names for each UCSC track before querying, as some UCSC tracks have non-standard coordinate field designations.  

For instance, the RepeatMasker track uses `genoName`/`genoStart`/`genoEnd` instead of `chrom`/`chromStart`/`chromEnd`.  

Thus, we could count all RepeatMasker elements per bin as follows:
```
$ athena annotate-bins -z \
    ... 
    -u "rmsk:genoName,genoStart,genoEnd" -a count -n repeatmasker_count \
    ...
```

And we could further restrict this count to only _SINE_ elements as follows:
```
$ athena annotate-bins -z \
    ... 
    -u "rmsk:genoName,genoStart,genoEnd,repClass=SINE" -a count -n repeatmasker_sine_count \
    ...
```

--- 

### Step 5: Collapse correlated 1D annotations  

Many genomic anntations are correlated (_e.g._, microsatellites and low sequence uniqueness, or GC content and protein-coding exons). 

To control for the underlying correlation structure of the annotations included in [step 4 (above)](https://github.com/talkowski-lab/athena#step-4-annotate-1d-training-bins), the next step in the Athena workflow is to perform Eigendecomposition of the binwise annotations.  

This can be accomplished with `athena eigen-bins`.  

For example, we could decompose the annotations from [the example in step 4](https://github.com/talkowski-lab/athena#step-4-annotate-1d-training-bins) into the top three Eigenfeatures as follows:  
```
$ athena eigen-bins \
    --eigenfeatures 3 \
    --stats eigenfeature_stats.txt \
    -z \
    GRCh37.autosomes.5kb_bins.all.annotated.bed.gz \
    GRCh37.autosomes.5kb_bins.all.annotated.decomped.bed.gz
```  

If you include the `--stats` option, Athena will also generate a log file with the summary of all Eigenfeatures, including total variance explained by each, and the Spearman correlation of raw annotation tracks significantly correlated with each Eigenfeature.  

This log file can be useful for determining _what_ each Eigenfeature represents (roughly speaking).  

#### Transforming annotations prior to Eigendecomposition:  

An important guideline is that `athena eigen-bins` assumes all annotations approximately follow a Gaussian (_i.e._, normal) distribution.  

Therefore, it is recommended to perform a cursory exploration of each annotation prior to annotation Eigendecomposition. Athena provides a simple utility for exploring each annotation via `athena feature-hists`.

For example, we could visualize the distributions of the six annotations applied in [the example in step 4](https://github.com/talkowski-lab/athena#step-4-annotate-1d-training-bins) as follows:  
```
$ athena feature-hists \
    GRCh37.autosomes.5kb_bins.all.annotated.bed.gz \
    ./athena_feature_hist
```  
This will generate a small `.png` file with a histogram of each feature for all bins.  

Upon manually reviewing each histogram, we can observe a few non-normal annotations:  
 * `mapability_100mers` is left-skewed; and
 * `ovary_peak_dnase_hypersensitivity` and `ovary_sense_strand_expression` are right-skewed.

Our Eigendecomposition from the example above will more accurately capture the variance among bins if we can transform these annotations to appear approximately Gaussian prior to Eigendecomposition.  

Fortunately, Athena can perform a few simple transformations of specified annotation(s), as summarized below:  

| Transformation | Athena flag | Details |  
| :--- | :--- | :--- |  
| log<sub>10</sub> | `--log-transform` | Transforms values as `log10(x + E)`, where `E = max(x) / 1000` |  
| Square root | `--sqrt-transform` | Transforms values as square roots |  
| Exponential | `--exp-transform` | Transforms values as `e ^ x` |  
| Square | `--square-transform` | Transforms values as `x ^ 2` |  
| Box-Cox | `--boxcox-transform` | Box-Cox power transformation of `x + E`, where `E = max(x) / 1000` |  

Note that both `athena eigen-bins` and `athena feature-hists` accept these transformation options, allowing you to re-plot features after transformation to check the distributions of your transformed data as it will be used during Eigendecomposition.  

#### Eigendecomposition with transformed annotations  

Putting everything together, the final `athena eigen-bins` call with necessary annotation transformations would be as follows:  
```
$ athena eigen-bins \
    --eigenfeatures 3 \
    --log-transform ovary_peak_dnase_hypersensitivity \
    --log-transform ovary_sense_strand_expression \
    --boxcox-transform mapability_100mers \
    --stats eigenfeature_stats.txt \
    -z \
    GRCh37.autosomes.5kb_bins.all.annotated.bed.gz \
    GRCh37.autosomes.5kb_bins.all.annotated.decomped.bed.gz
```

---  

### Step 6: Count SVs per training bin  

After annotating all bins with Eigenfeatures as described in [step 5 (above)](https://github.com/talkowski-lab/athena#step-5-collapse-correlated-1d-annotations), we next need to intersect the set of SVs and bins to be used for training the mutation rate model.  

This can be accomplished with `athena count-sv`; for example:  
```
$ athena count-sv -z \
    --sv-format vcf \
    --comparison interval \
    GRCh37.autosomes.5kb_bins.training.bed.gz \
    athena_training_deletions.vcf.gz \
    GRCh37.autosomes.5kb_bins.training.wSV.bed.gz
```

---  

### Step 7: Train 1D mutation rate model  

Training the 1D mutation rate model reqiuires the following inputs:  
 1. A BED file of training bins (see [step 3](https://github.com/talkowski-lab/athena#step-3-create-1d-bins-genome-wide)) with counts of SVs overlapping each bin (see [step 6](https://github.com/talkowski-lab/athena#step-6-count-svs-per-training-bin)); and 
 2. A BED of bins that:
    * includes all training bins; and
    * has been annotated with the desired covariates to include in the mutation rate model (see [step 4](https://github.com/talkowski-lab/athena#step-4-annotate-1d-training-bins)); and
    * for which the correlation structure between annotations has been collapsed (see [step 5](https://github.com/talkowski-lab/athena#step-5-collapse-correlated-1d-annotations)).  

Given these two inputs, you can train a 1D mutation rate model with `athena train-mu`.  

Following the example dataset from above, we could create a mutation rate model as follows:  
```
$ athena train-mu \
    GRCh37.autosomes.5kb_bins.training.wSV.bed.gz \
    GRCh37.autosomes.5kb_bins.all.annotated.decomped.bed.gz
```

**Please note that this functionality is still in early development**

---  

### Step 8: Create 2D bin-pairs genome-wide  

After training the 1D mutation rate model, the next component of the Athena workflow is to build the 2D mutation rate model.  

The first step of the 2D modeling process is to create bin-pairs genome-wide.  

Bin pairing can be accomplished with `athena pair-bins`, which follows many similar parameters as `athena make-bins` (as described in [step 3, above](https://github.com/talkowski-lab/athena#step-3-create-1d-bins-genome-wide)).  

However, there are two key differences in `athena pair-bins` as compared to `athena make-bins`:
  1. Users must specify a maximum bin-to-bin distance to be considered during pairing by using the `--max-dist-all` argument; and
  2. Blacklists are applied to the outer span of the distance between the two bins. For more information, see the [`bedtools pairtobed` documentation](https://daler.github.io/pybedtools/autodocs/pybedtools.bedtool.BedTool.pair_to_bed.html).  

Finally, a smaller distance for pairs to be used in model training can be specified using `--max-dist-training`.  

For an example, using the genome-wide bin BED file output from [step 3](https://github.com/talkowski-lab/athena#step-3-create-1d-bins-genome-wide)), and based on the parameters estimated in [step 2](https://github.com/talkowski-lab/athena#step-2-optional-calculate-sv-size--spacing-distributions) of 160kb maximum distance between any two bins and 50kb maximum distance between any two bins used for training, we could 

### A note on design  

This package was designed with canonical CNVs from [the gnomAD-SV callset](https://gnomad.broadinstitute.org/downloads) in mind.  

To that end, it assumes input data follows gnomAD-SV formatting standards. This may cause issues for alternative styles of SV representation, for SV types other than canonical CNVs, or different metadata labels.  

If using non-gnomAD data with Athena, please compare your VCF formatting standards, and the `INFO` field in particular.  

You can read more about the gnomAD-SV dataset [in the corresponding preprint](https://broad.io/gnomad_sv).  

### A note on parallellization  

Some steps (especially [bin annotation](https://github.com/talkowski-lab/athena#step-4-annotate-1d-training-bins)) can be computationally demanding. In practice, it is recommended to parallelize by chromosome where possible.  

### About the name  

This package is named after [Athena](https://en.wikipedia.org/wiki/Athena), the Greek goddess of wisdom, strategy, tactics, and mathematics. She was selected as the namesake for this package given that it relies on understanding the features that influence structural variation mutation rates (_wisdom_), incorporating those features into a statistical model (_mathematics_), and using these models to infer which components of the genome are vulnerable to changes in copy number (a kind of genomic _tactics_/_strategy_).