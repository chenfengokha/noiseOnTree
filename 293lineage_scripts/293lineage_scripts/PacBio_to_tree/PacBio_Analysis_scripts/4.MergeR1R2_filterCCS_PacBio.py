import os
import sys
import pandas as pd
import argparse
## import self define modules
sys.path.append("/mnt/data/home/lzz/project/Reconstruct_lineage/scripts/NewReTree")
from modules.SequenceParser import formFASTA


def getCCSreadDF(CCS_result_path):
    CCS_fa = open(CCS_result_path, 'r').read()
    readName = []
    readSeq = []
    for record in formFASTA(CCS_fa):
        name, seq = record
        readName.append(name)
        readSeq.append(seq)
    # build dataframe
    ccs_df = pd.DataFrame({"name":readName, "sequence":readSeq})
    ccs_df["name"] = ">" + ccs_df["name"]
    return ccs_df

def Parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--R1', type=str, nargs='?', required=True, help="Input R1 CCS filter FASTA file")
    parser.add_argument('-2', '--R2', type=str, nargs='?', required=True, help="Input R2 CCS filter FASTA file")
    parser.add_argument('-o', '--union', type=str, nargs='?', help="Output Union CCS filter FASTA file")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    options = Parsers()
    # input file path
    R1_ccs_path = options.R1
    R2_ccs_path = options.R2
    # output path
    union_ccs_path = options.union
    # run
    R1_ccs_df = getCCSreadDF(R1_ccs_path)
    R2_ccs_df = getCCSreadDF(R2_ccs_path)
    # merge two df
    R1_R2_ccsUnion_df = pd.merge(R1_ccs_df, R2_ccs_df, on=["name", "sequence"], how="outer")
    # write out results
    R1_R2_ccsUnion_df.to_csv(union_ccs_path, sep='\n', index=0, header=0)