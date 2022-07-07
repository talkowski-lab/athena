# Athena
#### A toolkit for exploring structural variation mutation rates and dosage sensitivity

Copyright (c) 2019-Present, [Ryan L. Collins](mailto:rlcollins@g.harvard.edu) and the Talkowski Laboratory.  
Distributed under terms of the [MIT License](/LICENSE) (see `LICENSE`).  

#### Note: the functionality of this repo is incomplete and is under active development.  

---  

## Table of Contents

#### Getting started
  * [Run Athena with Docker](https://github.com/talkowski-lab/athena#run-from-docker)
  * [Manual Installation](https://github.com/talkowski-lab/athena#manual-installation)
  * [Invoking Athena](https://github.com/talkowski-lab/athena#invoking-athena)

#### Other
  * [A note on design](https://github.com/talkowski-lab/athena#a-note-on-design)
  * [About the name](https://github.com/talkowski-lab/athena#about-the-name)


---  

## Run from Docker

The recommended way to run Athena is from its dedicated Docker container hosted on Google Container Registry. This will handle all dependencies and installation for you, and ensure you are running the latest version.

```
$ docker pull us.gcr.io/broad-dsmap/athena
$ docker run --rm -it us.gcr.io/broad-dsmap/athena
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
  annotate-bins          Annotate bins
  annotate-pairs         Annotate pairs
  breakpoint-confidence  Annotate breakpoint uncertainty
  count-sv               Intersect SV and 1D bins or 2D bin-pairs
  eigen-bins             Eigendecomposition of annotations
  feature-hists          Plot bin annotation distributions
  feature-stats          Compute feature distributions
  make-bins              Create sequential bins
  mu-predict             Predict mutation rates with a trained model
  mu-query               Query a mutation rate matrix
  mu-train               Train mutation rate model
  pair-bins              Create pairs of bins
  slice-remote           Localize slices of remote genomic data
  transform              Transform one or more annotations
  vcf-filter             Filter an input VCF
  vcf-stats              Get SV size & spacing
```

Athena has numerous subcommands. Specify `--help` with any subcommand to see a list of options available.  

### A note on design  

This package was designed with canonical CNVs from [the gnomAD-SV callset](https://gnomad.broadinstitute.org/downloads) in mind.  

To that end, it assumes input data follows gnomAD-SV formatting standards. This may cause issues for alternative styles of SV representation, for SV types other than canonical CNVs, or different metadata labels.  

If using non-gnomAD data with Athena, please compare your VCF formatting standards, and the `INFO` field in particular.  

You can read more about the gnomAD-SV dataset [in the corresponding preprint](https://broad.io/gnomad_sv).  

### About the name  

This package is named after [Athena](https://en.wikipedia.org/wiki/Athena), the Greek goddess of wisdom, strategy, tactics, and mathematics. She was selected as the namesake for this package given that it relies on understanding the features that influence structural variation mutation rates (_wisdom_), incorporating those features into a statistical model (_mathematics_), and using these models to infer which components of the genome are vulnerable to changes in copy number (a kind of genomic _tactics_/_strategy_).
