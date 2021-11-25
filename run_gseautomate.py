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

# Takes processed dataframe as input, converts to alternative
# gene symbol for downstream GSEA run


def geneConversion(df_toPreRank, gdb, delim='\t'):

    # For gene_id in df_gene_symbol,
    # Check if gene_id exists in col 1 of gdb
    # if true, replace with value in col 3 of gdb
    # write to new df
    # if false, store in missing_df, print
    # return converted gene list for submission to Prerank

    # Load in dataframe from processID
    df_gene_symbol = df_toPreRank['gene_symbol']

    # Load in dataframe from gdb
    df_gdb = pd.read_csv(gdb, sep='\t')
    #print(df_gdb)

    #for gene_symbol in df_gene_symbol:


# Takes processed and converted dataframe as input,
# running it through the Prerank function provided by gseapy
def preRank(df_toPreRank, pc_limit):

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
    #pc_limit = 11
    frames = []
    for col in df_toPreRank.columns[1:pc_limit]:
        df_pc = df_toPreRank[col]
        frames.append(pd.concat([df_gene_symbol, df_pc], axis=1, join='inner'))

    # To rank, iterate over entire object
    i = 1
    for frame in frames:
        rnk = frame
        # Convert all gene_symbol to uppercase
        rnk = frame['gene_symbol'].str.upper()
        # Remove any column in df with na value
        rnk = frame.dropna()
        # Remove duplicate IDs, keep highest value
        rnk = frame.groupby('gene_symbol', as_index=False).max()
        print(rnk)
        dirname = 'PC_%d' % (i,)
        i += 1

        # Set enrichr database, can choose from Human, Mouse, Yeast, Fly, Fish, Worm)
        gss = gp.get_library_name(organism='Fish')

        # Add support for defining organism dataset to be used
        gmt_dict = gp.parser.gsea_gmt_parser(
            '/Users/nathanielmaki/.cache/gseapy/enrichr.GO_Biological_Process_2018.gmt', organism='Fish')

        # TODO: Update to larger default image number (20 to 50)
        pre_res = gp.prerank(rnk=rnk, gene_sets=gmt_dict, processes=4,
                             permutation_num=100, outdir='./GO_Biological_Process_2018/'+dirname, format='png', seed=6)
        print(pre_res)
    return pre_res

# Verifies that the initial input file exists, if not throw warning
def inputVerification(parsed_args):
    # Check if input file exists, throw warning if it does not
    infile = parsed_args.infile
    if not os.path.isfile(parsed_args.infile):
        print("The input file %s does not exist, quitting." %
              parsed_args.infile)
        sys.exit()

    # # Check if gene db file exists, throw warning if it does not
    # gdb = parsed_args.gdb
    # if not os.path.isfile(parsed_args.gdb):
    #     print("The gene database file %s is missing, quitting." %
    #           parsed_args.gdb)
    #     sys.exit()

    # Auto-generate an output file name based upon checked input
    if parsed_args.outfile == '':
        outfile = infile.split('.')[-2] + '.biomart.txt'
    else:
        outfile = parsed_args.outfile

    # return infile, gdb, outfile
    return infile, outfile


def main():
    # Definitions for input and output files
    parser = argparse.ArgumentParser(
        description='Load gene list for processing.')

    parser.add_argument("-i", "--infile", required=True, dest="infile",
                        help="Input gene list, acquired from DESeq2, PCA, etc")

    parser.add_argument("-o", "--outfile", dest="outfile",
                        help="Output filename, should end in '.biomart.txt'")

    parser.add_argument("-pc", "--pc_limit", dest="pc_limit",
                        help="Number of PCs to work with", type=int)

    # TODO: Implement this feature
    parser.add_argument("-gdb", "--gene_database", dest="gdb",
                        help="Database to convert genes against")

    #infile, gdb, outfile = inputVerification(parser.parse_args())

    infile, outfile = inputVerification(parser.parse_args())

    df_toPreRank = processID(infile, outfile)

    #geneConversion(df_toPreRank, gdb)

    args = parser.parse_args()

    preRank(df_toPreRank, args.pc_limit)


if __name__ == '__main__':
    main()
