# Athena
#### A toolkit for exploring structural variation mutation rates and dosage sensitivity

Copyright (c) 2019, [Ryan L. Collins](mailto:rlcollins@g.harvard.edu) and the Talkowski Laboratory.  
Distributed under terms of the [MIT License](/LICENSE) (see `LICENSE`).  

---  

## Table of Contents

#### Getting started
  * [Run Athena with Docker](https://github.com/talkowski-lab/athena#run-from-docker)
  * [Manual Installation](https://github.com/talkowski-lab/athena#manual-installation)

#### The Athena workflow
_Mutation rate modeling_ 
  1. [Filter SVs for mutation rate training](https://github.com/talkowski-lab/athena#step-1)  

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

For instance, from the [gnomAD-SV v2 dataset](https://gnomad.broadinstitute.org/downloads) (`gnomad_v2_sv.sites.vcf.gz`), you could run the following:  
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

We can also observe that roughly half of all deletions are no farther than a few kilobases away from the next-nearest deletion, which means that bin sizes of several kilobases should result in adequately dense data.  

Based on this logic along with some approximate rounding, we can decide on the two key parameters used in the remaining steps of the mutation rate model:
```
1D bin size: 4kb
2D maximum bin-pair distance: 40kb
```

--- 

### A note on design
This package was designed with [the gnomAD-SV callset](https://gnomad.broadinstitute.org/downloads) in mind. To that end, it assumes input data follows gnomAD-SV formatting standards. This may cause issues for alternative styles of SV representation, or different metadata labels. If using non-gnomAD data with Athena, please compare your VCF formatting standards, and the `INFO` field in particular.  

You can read more about the gnomAD-SV dataset [in the corresponding preprint](https://broad.io/gnomad_sv).

### About the name
This package is named after [Athena](https://en.wikipedia.org/wiki/Athena), the Greek goddess of wisdom, strategy, tactics, and mathematics. She was selected as the namesake for this package given that it relies on understanding the features that influence structural variation mutation rates (_wisdom_), incorporating those features into a statistical model (_mathematics_), and using these models to infer which components of the genome are vulnerable to changes in copy number (a kind of genomic _tactics_/_strategy_).