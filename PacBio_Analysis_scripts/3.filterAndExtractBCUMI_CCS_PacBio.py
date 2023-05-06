import sys
import os
import regex
from collections import defaultdict
import operator
import multiprocessing
import argparse
## import self define modules
sys.path.append("/mnt/data/home/lzz/project/Reconstruct_lineage/scripts/NewReTree")
from modules.SequenceParser import convert2string, fqTOfa, formFASTA, seqco, seqre, formFASTQ

# some function

def extractBC_UMI(read, P2, P4, mismatch=2):
    ## extract BC UMI barcode pattern
    BC_UMI_fwd = "(?e)(?P<p2>%s){e<=%d}(?P<raw_seq>.*)(?P<UMI>.{10})(?P<BC>.{16})(?P<p4>%s){e<=%d}" % (P2, mismatch, P4, mismatch)
    BC_UMI_rev = "(?e)(?P<p4>%s){e<=%d}(?P<BC>.{16})(?P<UMI>.{10})(?P<raw_seq>.*)(?P<p2>%s){e<=%d}" % (seqre(seqco(P4)), mismatch, seqre(seqco(P2)), mismatch)
    # trim P2 and P4, extract BC and UMI
    fwd_match = regex.search(BC_UMI_fwd, read)
    rev_match = regex.search(BC_UMI_rev, read)
    if fwd_match:
        BC = seqre(seqco(fwd_match.groupdict()["BC"]))
        UMI = seqre(seqco(fwd_match.groupdict()["UMI"]))
        raw_seq = fwd_match.groupdict()["raw_seq"]
        return BC, UMI, raw_seq
    elif rev_match:
        BC = rev_match.groupdict()["BC"]
        UMI = rev_match.groupdict()["UMI"]
        raw_seq = seqre(seqco(rev_match.groupdict()["raw_seq"]))
        return BC, UMI, raw_seq
    else:
        return None

def trimFPRP(read, FP, RP, mismatch=2):
    ## trim primers pattern
    FP_pattern = "(?e)(?P<FP>%s){e<=%d}" % (FP, mismatch)
    RP_pattern = "(?e)(?P<RP>%s){e<=%d}" % (RP, mismatch)
    ## search 
    FP_match = regex.search(FP_pattern, read)
    RP_match = regex.search(RP_pattern, read)
    if FP_match and RP_match:
        return "both", read[FP_match.end(): RP_match.start()]
    elif FP_match:
        return "fwd", read[FP_match.end():]
    elif RP_match:
        return "rev", read[:RP_match.start()]
    else:
        return "noMatch", read

def diffbase(seq1, seq2):
    return sum(map(operator.ne, seq1, seq2))
    
def trimPolyA(seq, windowSize=6, mismatch=1):
    PolyA = "A" * windowSize
    seq = seq.upper()
    cutPos = -1
    minDiff = -1
    for i in range(len(seq)-windowSize+1):
        slideW = seq[i:i+windowSize]
        if slideW.startswith("AAA"):
            diffNum = diffbase(slideW, PolyA)
            if diffNum <= mismatch:
                if minDiff == -1 and cutPos == -1:
                    minDiff = diffNum
                    cutPos = i
                else:
                    if diffNum <= minDiff:
                        minDiff = diffNum
                        cutPos = cutPos
                        break
    if cutPos != -1:
        return seq[:cutPos]
    else:
        return seq     

    
def filterCCS(header, read, P2, P4, FP, RP, polyAwin=6, pAmismatch=1, mismatch=2):
    ## get zw ID
    zwID = header.split('/')[1]
    ## trim primers pattern
    FP_pattern = "(?e)(?P<FP>%s){e<=%d}" % (FP, mismatch)
    RP_pattern = "(?e)(?P<RP>%s){e<=%d}" % (RP, mismatch)
    # 1. trim P2 and P4, extract BC and UMI
    trimP2P4_result = extractBC_UMI(read, P2, P4, mismatch)
    if trimP2P4_result:
        BC, UMI, raw_seq = trimP2P4_result
        noPolyA_seq = trimPolyA(raw_seq, windowSize=polyAwin, mismatch=pAmismatch)
        rType, final_seq = trimFPRP(noPolyA_seq, FP, RP, mismatch)
        #rType, final_seq = trimFPRP(raw_seq, FP, RP, mismatch)
        if rType == "both":
            name = "zwID={} BC={} UMI={} P2P4 FPRP".format(zwID, BC, UMI)
            return "pass", name, final_seq, zwID, BC, UMI
        elif rType == "fwd":
            name = "zwID={} BC={} UMI={} P2P4 FP".format(zwID, BC, UMI)
            return "pass", name, final_seq, zwID, BC, UMI
        elif rType == "rev":
            name = "zwID={} BC={} UMI={} P2P4 RP".format(zwID, BC, UMI)
            return "pass", name, final_seq, zwID, BC, UMI
        elif rType == "noMatch":
            name = "zwID={} BC={} UMI={} P2P4 noMatch".format(zwID, BC, UMI)
            return "pass", name, final_seq, zwID, BC, UMI
    else:
        return "nope", header.strip(), read.strip(), zwID, "noBC", "noUMI"

def mycallback(results):
    filre, header, seq, zwid, BC, UMI = results
    if filre == "pass" and seq != "":
        outFA.write(">{}\n{}\n".format(header, seq))
        outFA.flush()
    elif filre == "nope":
        outNoPass.write(">{}\n{}\n".format(header, seq))
        outNoPass.flush()
    outSta.write("{}\t{}\t{}\t{}\n".format(zwid, BC, UMI, filre))
    outSta.flush()


def Parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=str, nargs='?', required=True, help="Input the FASTQ/FASTA format file")
    parser.add_argument('-o', '--outfile', type=str, nargs='?', required=True, help="Output the FASTA format file")
    parser.add_argument('-s', '--statistics', type=str, nargs='?', required=True, help="Output the statistics file about the filter")
    parser.add_argument('--nopass', type=str, nargs='?', required=True, help="Output the reads can not pass filter, FASTA format file")
    parser.add_argument('--format', type=str, nargs='?', required=True, default= 'FASTQ', choices=['FASTQ', 'FASTA'], help="Format of input file, FASTQ or FASTA")
    parser.add_argument('--FP', type=str, nargs='?', required=True, help="Forward primer of PCR amplification")
    parser.add_argument('--RP', type=str, nargs='?', required=True, help="Reverse primer of PCR amplification")
    parser.add_argument('--inFP', type=str, nargs='?', required=True, help="Flanking forward primer of target sites")
    parser.add_argument('--inRP', type=str, nargs='?', required=True, help="Flanking reverse primer of target sites")
    parser.add_argument('--pAw', type=int, nargs='?', required=True, default=6, help="PolyA window size, use to trim polyA")
    parser.add_argument('--pMs', type=int, nargs='?', required=True, default=1, help="PolyA mismatch permit")
    parser.add_argument('--mismatch', type=int, nargs='?', required=True, default=2, help="Primers mismatch permit")
    parser.add_argument('--numcpu', type=int, nargs='?', required=True, help="Number of cpu needed")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    options = Parsers()
    ## input file path
    inputFile = open(options.infile, 'r')
    ## output file path
    outFA = open(options.outfile, 'w')
    outNoPass = open(options.nopass, 'w')
    outSta = open(options.statistics, 'w')
    # write statistics header
    outSta.write("{}\t{}\t{}\t{}\n".format("zwID", "BC", "UMI", "pass_filter"))
    ## some input parameters
    P2 = options.FP
    P4 = options.RP
    inFP = options.inFP
    inRP = options.inRP
    polyAwin = options.pAw
    pAmismatch = options.pMs
    mismatch = options.mismatch
    ## get file handle
    if options.format == "FASTQ":
        file_handle = fqTOfa(inputFile)
    elif options.format == "FASTA":
        file_handle = formFASTA(inputFile.read())
    ## run filter 
    pfq = multiprocessing.Pool(options.numcpu)
    for record in file_handle:
        header, raw_seq = record
        prf = pfq.apply_async(filterCCS, args=(header, raw_seq, P2, P4, inFP, inRP, polyAwin, pAmismatch, mismatch, ), callback=mycallback)
    print("waite")
    pfq.close()
    pfq.join()
    print("done")
    ## close and save files
    inputFile.close()
    outFA.close()
    outNoPass.close()
    outSta.close()






