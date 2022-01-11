from gseapy.parser import Biomart
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import argparse
import sys
import os
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Processes input file (in this case, PCA loading scores)
# Removes ENSEMBL Identifier and reduces to gene ID


def processID(infile, delim='\t'):

    # Load in gene list from data frame
    df = pd.read_csv(infile, sep='\t')

    # Give column 0 an intermediary ID
    df.rename(columns={list(df)[0]: "NAME"}, inplace=True)

    # Reduce gene ID to gene symbol, removing ENSEMBL flag (preceded ID by a _ )
    df['gene_symbol'] = df['NAME'].str.split('_').str[1]

    # Shift last column (gene_symbol) to be first
    cols = list(df.columns)
    cols = [cols[-1]] + cols[:-1]
    df = df[cols]

    # Drop redundant second column containing ENSEMBL IDs combined with gene symbols
    df_toPreRank = df.drop(columns=['NAME'])
    return df_toPreRank

# Takes processed and converted dataframe as input,
# running it through the Prerank function provided by gseapy


def preRank(df_toPreRank, col_limit, library, organism):

    # Dump gene column as Dataframe
    df_gene_symbol = df_toPreRank['gene_symbol']

    # Iterate from 1 to n, dump column i
    # Take dump of gene column and join with column i as dataframe
    # Stack as larger dataframe with n objects
    # Each object is of gene_symbol and PCi
    # [[gene_column, PC1], [gene_column, PC2], [i, j]]
    # nested-nested list of arrays
    # Return entire [[gene_column, PC1], [gene_column, PC2], [i, j]]

    frames = []

    column_headers = list(df_toPreRank.columns.values)
    print(column_headers)

    for col in df_toPreRank.columns[1:col_limit]:
        df_col = df_toPreRank[col]
        frames.append(pd.concat([df_gene_symbol, df_col], axis=1, join='inner'))

    # To rank, iterate over entire object
    i = 1
    for frame in frames:
        
        rnk = frame

        # Remove any column in df with na value
        rnk = frame.dropna()

        # Convert all gene_symbol to uppercase
        rnk = frame['gene_symbol'].str.upper()

        # Remove duplicate IDs, keep highest value
        rnk = frame.groupby('gene_symbol', as_index=False).max()

        print(rnk)
        dirname = column_headers[i]
        i += 1

        # Add support for defining organism dataset to be used
        gmt_dict = gp.parser.gsea_gmt_parser(
             '/Users/nathanielmaki/.cache/gseapy/enrichr.GO_Biological_Process_2018.gmt', organism='Fish')
        
        #gmt_dict = gp.parser.gsea_gmt_parser(library, organism)

        pre_res = gp.prerank(rnk=rnk, gene_sets=gmt_dict, processes=4,
                             permutation_num=100, outdir='./GO_Biological_Process_2018/'+dirname, format='png', seed=6)
        print(pre_res)
    return pre_res

# Verifies that the initial input file exists, if not throw warning


def inputVerification(parsed_args):
    infile = parsed_args.infile
    if not os.path.isfile(parsed_args.infile):
        print("The input file %s does not exist, quitting." %
              parsed_args.infile)
        sys.exit()
    return infile


def main():
    # Definitions for input and output files
    parser = argparse.ArgumentParser(
        description='Load gene list for processing.')

    parser.add_argument("-i", "--infile", required=True, dest="infile",
                        help="Input gene list, acquired from DESeq2, PCA, etc")

    #parser.add_argument("-o", "--outfile", dest="outfile",
                        #help="Output filename")

    parser.add_argument("-cl", "--col_limit", dest="col_limit",
                        help="Number of columns to parse through", type=int)

    parser.add_argument("-lib", "--library", dest="library",
                        help="Enrichr library to pull from")

    parser.add_argument("-org", "--organism", dest="organism",
                        help="Organism to pull enrichr library for")

    infile = inputVerification(parser.parse_args())

    df_toPreRank = processID(infile)

    args = parser.parse_args()

    preRank(df_toPreRank, args.col_limit, args.library, args.organism)


if __name__ == '__main__':
    main()
