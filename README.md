# gseautomate


## Introduction

gseautomate is a wrapper around the GSEApy python package, adding to it's functionality. Utilizing Pandas, users are able to supply any gene matrix file, assuming that the gene column is in ENSEMBL format (required for preprocessing).

This utility will iterate through each column of your supplied file together with the gene ID, and submits it to GSEApy. It generates a simple named output (col_1, etc) underneath the library directory from GSEA that it was ran against.

There are a few bugs with user-provided library and organism IDs, for the moment they need to be changed manually, and the organism-specific library sourced from ENRICHR.

## Requirements

To install the required python libraries, create a virtual python environment using `venv`, activate it, and execute `pip install -r requirements.txt`. 

Next you'll need to acquire the required library files from the Enrichr [website](https://maayanlab.cloud/Enrichr/).

From the bottom of the page, select the organism specific database to work with (FlyEnrichr, etc), follow the `Libraries` link, and download the desired library.

#### Note: The below requirement is under active development to be deprecated once the path/output dir can be given as an argument

Then move the library file into the `~/.cache/gseapy` directory. If it does not exist, please create it.

Lastly, you need to update the path located in the preRank function to point to your specific enrichr library, and the name of the output dir that you want created:

```
gmt_dict = gp.parser.gsea_gmt_parser(
             '/home/nathanielmki/.cache/gseapy/enrichr.KEGG_2019.gmt', organism='Fish')
```
```
 pre_res = gp.prerank(rnk=rnk, gene_sets=gmt_dict, processes=4,
                             permutation_num=100, outdir='./KEGG_2019/'+dirname, format='png', seed=6)
```
## Usage

An example execution:

```python3 ~/Documents/Github/mdibl/gseautomate/run_gseautomate.py -i jcoffman.meta.wt-Q2_full.DESeq2_out.tsv -cl 11```

* Use the `-i` argument for the supplied file, and `-cl` for the number of columns to iterate through (n+1 more than you need, as index column counts).

* The only `required` input is the actual file itself. No `-cl` argument will cause gseautomate to iterate through all columns of your file.
