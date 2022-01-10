# gseautomate


## Introduction

gseautomate acts as a functional wrapper around the GSEApy python package, extending it's functionality. Utilizing Pandas, users are able to supply any gene matrix file, assuming that the gene column is in ENSEMBL format (required for preprocessing).

This utility will iterate through each column of your supplied file together with the gene ID, and submits it to GSEApy. It generates a simple named output (col_1, etc) underneath the library directory from GSEA that it was ran against.

There are a few bugs with user-provided library and organism IDs, for the moment they need to be changed manually, and the organism-specific library sourced from ENRICHR.

## Usage

An example execution:

```python3 ~/Documents/Github/mdibl/gseautomate/run_gseautomate.py -i jcoffman.meta.wt-Q2_full.DESeq2_out.tsv -cl 11```

* Use the -i argument for the supplied file, and -cl for the number of columns to iterate through (n+1 more than you need, as index column counts).

* The only `required` input is the actual file itself. No -cl argument will cause gseautomate to iterate through all columns of your file.
