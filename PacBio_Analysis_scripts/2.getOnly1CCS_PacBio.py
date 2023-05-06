import sys
from collections import defaultdict
import os
import multiprocessing
import argparse

## import self define modules
sys.path.append("/mnt/data/home/lzz/project/Reconstruct_lineage/scripts/NewReTree")
from modules.SequenceParser import formFASTQ


def AvgQual(qual):
    return sum([ord(i) for i in qual]) / len(qual)

def getCCSresultAvgQual(ccs_result_fq):
    qual_list = []
    with open(ccs_result_fq, 'r') as infq:
        for record in formFASTQ(infq):
            _, _, _, qual = record
            avg_qual = AvgQual(qual)
            qual_list.append(avg_qual)
    return(qual_list)

# def getCutoffAvgQual(qual_list, cutoff=0.95):
#     cutoff_index = round(len(qual_list) * cutoff)
#     qual_list.sort(reverse=True)
#     return qual_list[cutoff_index]

# def filterCCSresult(ccs_result_fq, output_path, cutoff=0.95):
#     qual_list = getCCSresultAvgQual(ccs_result_fq)
#     cutoff_qual = getCutoffAvgQual(qual_list)
#     with open(ccs_result_fq, 'r') as inf, open(output_path, 'w') as outf:
#         for record in formFASTQ(inf):
#             header, seq, ind, qual = record
#             avg_qual = AvgQual(qual)
#             if avg_qual >= cutoff_qual:
#                 outf.write("{}\n{}\n{}\n{}\n".format(*record))

def getFQdict(ccs_fq_path):
    zw_fq_dict = defaultdict(list)
    with open(ccs_fq_path, 'r') as infq:
        for record in formFASTQ(infq):
            header,_,_,_ = record
            zwID = header.split("/")[1]
            zw_fq_dict[zwID].append(record)
    return zw_fq_dict

def getHighestQualCCS(fq_record_list):
    max_qual = 0
    max_qual_ccs = ""
    for each in fq_record_list:
        _, _, _, qual = each
        avg_qual = AvgQual(qual)
        if avg_qual > max_qual:
            max_qual = avg_qual
            max_qual_ccs = each
    return max_qual_ccs

def getOnlyCCS(out_path, zw_fq_dict):
    outf = open(out_path, 'w')
    for k, v in zw_fq_dict.items():
        if len(v) == 1:
            out_record = v[0]
            outf.write("{}\n{}\n{}\n{}\n".format(*out_record))
        else:
            out_record = getHighestQualCCS(v)
            outf.write("{}\n{}\n{}\n{}\n".format(*out_record))
    outf.close()

def Parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=str, nargs='?', required=True, help="Input RAW CCS FASTQ file")
    parser.add_argument('-o', '--outfile', type=str, nargs='?', help="Output only one CCS reads for each zwID")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    options = Parsers()
    # input path
    raw_CCS_fq_path = options.infile
    # output path
    out_CCS_fq_path = options.outfile
    ## run
    zw_FQ_dict = getFQdict(raw_CCS_fq_path)
    getOnlyCCS(out_CCS_fq_path, zw_FQ_dict)
