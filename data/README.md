# Athena
#### A toolkit for exploring structural variation mutation rates and dosage sensitivity

Copyright (c) 2019, [Ryan L. Collins](mailto:rlcollins@g.harvard.edu) and the Talkowski Laboratory.  
Distributed under terms of the [MIT License](/LICENSE) (see `LICENSE`).  

---  

## Resources

The Athena distribution comes packaged with supplementary resources, described below.  

### Athena mutation rate modeling resources  

**SV selection blacklist:** `athena.SV_selection_blacklist.v1.GRCh37.bed.gz`  
This is a BED file of regions under mutational constraint. This file can be used as a blacklist in `athena filter-vcf` using the `-x`/`--blacklist` option to remove SVs likely under strong selection. This is an important step during mutation rate modeling; without excluding SVs under strong selection, the mutation rate model will be partially confounded by strong selection.    

This file was compiled by taking the union of two datasets: (1) all exons from the [Gencode v19](https://www.gencodegenes.org/human/release_19.html) canonical transcripts of protein-coding genes in the top half of all genes (n=9,857 genes) based on constraint against loss-of-function point mutations [(i.e., LOEUF)](https://gnomad.broadinstitute.org/downloads) from the [gnomAD v2.1 SNV and indel analyses](https://broad.io/gnomad_lof), and (2) a list of all mutationally constrained DNA elements [per GERP++](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001025) (GERP RS-score ≥ 3).  

The nonredundant union of these two sources comprises 253,628 distinct intervals over a total of 44.5Mb, or roughly \~1.5% of the alignable primary GRCh37 reference assembly.  

### GRCh37 resources  

**N-mask:** `GRCh37.Nmask.bed.gz`  
This is a BED file of the N-masked regions of the GRCh37 reference assembly, as provided by the [UCSC Genome Browser](http://genome.ucsc.edu).  

**Genome file:** `GRCh37.genome`  
This is a [BEDTools-style](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html) genome file, which lists all contigs in the GRCh37 primary assembly and their lengths.  

### Other resources  

***SNV mutation rates per trinucleotide context*** `snv_mutation_rates.Samocha_2014.tsv.gz`  
This is a tsv enumerating the relative mutation rates of all single-nucleotide changes given all possible trinucleotide contexts, as presented in [Samocha _et al._, _Nat. Genet._ (2014)](https://www.ncbi.nlm.nih.gov/pubmed/25086666).  
