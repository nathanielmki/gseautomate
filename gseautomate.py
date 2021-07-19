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


def processID(infile, outfile, delim='\t'):

    # Load in gene list from data frame
    df = pd.read_csv(infile, sep='\t')

    # Set Index Name
    df.rename(columns={"": "NAME"}, inplace=True)

    # Reduce gene ID to gene symbol, removing ENSEMBL flag (preceded ID by a _ )
    df['gene_symbol'] = df['NAME'].str.split('_').str[1]

    # Shift last column (gene_symbol) to be first
    cols = list(df.columns)
    cols = [cols[-1]] + cols[:-1]
    df = df[cols]

    # Drop redundant second column containing ENSEMBL IDs and gene symbols
    df_toPreRank = df.drop(columns=['NAME'])
    return df_toPreRank


def preRank(df_toPreRank):

    # Dump gene column as Dataframe
    df_gene_symbol = df_toPreRank['gene_symbol']

    # Iterate from 1 to n, dump column i
    # Take dump of gene column and join with column i as dataframe
    # Stack as larger dataframe with n objects
    # Each object is of gene_symbol and PCi
    # [[gene_column, PC1], [gene_column, PC2], [i, j]]
    # nested-nested list of arrays
    # Return entire [[gene_column, PC1], [gene_column, PC2], [i, j]]

    # TODO: Update for loop to work with user input for PC limit

    for col in df_toPreRank.columns[1:11]:
        df_pc = df_toPreRank[col]
        frames = pd.concat([df_gene_symbol, df_pc], axis=1, join='inner')
        #frames.rename(columns=frames.iloc[0]).drop(frames.index[0])
        #print(reduced)

    # To rank, iterate over entire object
        for frame in frames:
            # TODO: Strip header in preparation for Prerank submission
            rnk = frames
            print(rnk)
            # TODO: write to new folder for each PC (output is being overwritten)
            pre_res = gp.prerank(rnk=rnk, gene_sets='KEGG_2016', processes=4,
                             permutation_num=100, outdir='test/prerank_report_kegg', format='png', seed=6)
            print(pre_res)
        return(pre_res)

    # df_pc = df_toPreRank.iloc[:,[0,1]]
    # col_num = len(df_toPreRank.columns[1:])

    # # split dataframes in a list
    # split_df = []

    # for i in range(0, col_num, 1):
    #     split_df.append(df_toPreRank.iloc[:, [i,i+1]])
    #     print(split_df)

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
    # TODO: Implement this feature
    parser.add_argument("-pc", "--pclimit",
                        help="Number of PCs to work with", type=int)

    infile, outfile = inputVerification(parser.parse_args())

    df_toPreRank = processID(infile, outfile)

    preRank(df_toPreRank)


if __name__ == '__main__':
    main()
