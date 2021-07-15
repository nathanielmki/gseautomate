import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import argparse
import json
import sys
import os
import csv
import re
import string
from gseapy.parser import Biomart

def processID(infile, gene_symbol, delim='\t'):
    # Load in gene list from data frame
    df = pd.read_csv(infile, sep='\t')
    # Reduce gene ID to gene symbol, removing ENSEMBL flag (precedes ID by a _ )
    df['gene_symbol'] = df['NAME'].str.split('_').str[1]
    # Shift last column (gene_symbol) to be first
    cols = list(df.columns)
    cols = [cols[-1]] + cols[:-1]
    df = df[cols]
    # Drop redundant second column containing ENSEMBL IDs and gene symbols
    df_final = df.drop(columns=['NAME'])
    #print(df_final)
    df_final.to_csv(gene_symbol, sep='\t', encoding='utf-8')
    
    return gene_symbol

# def biomartConversion(gene_symbol, outfile, delim='\t'):
#     df = pd.read_csv(gene_symbol, sep='\t')
#     bm = Biomart()
#     # View validated marts
#     marts = bm.get_marts()
#     # View validated dataset
#     datasets = bm.get_datasets(mart='ENSEMBL_MART_ENSEMBL')
#     # View validated attributes
#     attrs = bm.get_attributes(dataset='drerio_gene_ensembl')
#     # View validated filters
#     filters = bm.get_filters(dataset='drerio_gene_ensembl')
#     # Pull out column 1 as list to feed to Biomart
#     col_one_list = df['NAME'].tolist()
#     #print(col_one_list)
#     # Query results
#     queries = col_one_list
#     results = bm.query(dataset='drerio_gene_ensembl', attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'go_id'],
#                         filters={'ensemble_gene_id': queries})
#     #results.to_csv(outfile, sep='\t', encoding='utf-8')
#     print(results)

def inputVerification(parsed_args):
    # Check if input file exists, throw warning if it does not
    infile = parsed_args.infile
    if not os.path.isfile(parsed_args.infile):
        print("The input file %s does not exist, quitting." %
              parsed_args.infile)
        sys.exit()

    # Auto-generate an output file name based upon checked input
    if parsed_args.outfile == '':
        outfile = infile.split('.')[-2] + '.biomart.txt'
    else:
        outfile = parsed_args.outfile

    return infile, outfile


def main():
    # Definitions for input and output files
    parser = argparse.ArgumentParser(
        description='Load gene list for processing.')

    parser.add_argument("-i", "--infile", required=True, dest="infile",
                        help="Input gene list, acquired from DESeq2, PCA, etc")

    parser.add_argument("-o", "--outfile", dest="outfile",
                        help="Output filename, should end in '.biomart.txt'")

    infile, outfile = inputVerification(parser.parse_args())

    gene_symbol = processID(gene_symbol)

    processID(infile, gene_symbol)

    biomartConversion(gene_symbol, outfile)


if __name__ == '__main__':
    main()