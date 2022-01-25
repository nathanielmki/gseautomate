# gseautomate


## Introduction

gseautomate acts as a functional wrapper around the GSEApy python package, extending it's functionality. Utilizing Pandas, users are able to supply any gene matrix file, assuming that the gene column is in ENSEMBL format (required for preprocessing).

This utility will iterate through each column of your supplied file together with the gene ID, and submits it to GSEApy. It generates a simple named output (col_1, etc) underneath the library directory from GSEA that it was ran against.

There are a few bugs with user-provided library and organism IDs, for the moment they need to be changed manually, and the organism-specific library sourced from ENRICHR.

## Requirements

To properly use this software you'll need to acquire the required library files from the Enrichr [website](https://maayanlab.cloud/Enrichr/).

From the bottom of the page, select the organism specific database to work with (FlyEnrichr, etc), follow the `Libraries` link, and download the desired library.

Then move the library file into the `~/.cache/gseapy` directory. If it does not exist, please create it.

## Usage

An example execution:

```python3 ~/Documents/Github/mdibl/gseautomate/run_gseautomate.py -i jcoffman.meta.wt-Q2_full.DESeq2_out.tsv -cl 11```

* Use the `-i` argument for the supplied file, and `-cl` for the number of columns to iterate through (n+1 more than you need, as index column counts).

* The only `required` input is the actual file itself. No `-cl` argument will cause gseautomate to iterate through all columns of your file.
