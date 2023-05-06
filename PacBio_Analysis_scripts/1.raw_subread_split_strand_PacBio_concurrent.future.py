#import multiprocessing
from concurrent.futures import ProcessPoolExecutor as Pool 
import os
import sys
import regex
import argparse


def runFilter(raw_line, mismatch=2):
    # define match pattern
    polyA_pattern = "(?e)(AAA)(?P<pA>AAAAAAAAAA){e<=%d}" % mismatch
    polyT_pattern = "(?e)(TTT)(?P<pA>TTTTTTTTTT){e<=%d}" % mismatch
    # fuzz match
    pA_match = regex.search(polyA_pattern, raw_line)
    pT_match = regex.search(polyT_pattern, raw_line)
    if pA_match and pT_match:
        return "other", raw_line
    elif pA_match:
        return "R1", raw_line
    elif pT_match:
        return "R2", raw_line
    else:
        return "other", raw_line
    
def mycallback(future):
    #print(future)
    r_type, line = future.result()
    if r_type == "R1":
        R1_sam.write(line)
        R1_sam.flush()
    elif r_type == "R2":
        R2_sam.write(line)
        R2_sam.flush()
    elif r_type == "other":
        other_sam.write(line)
        other_sam.flush()


def Parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--raw', type=str, nargs='?', required=True, help="Input raw PacBio subreads SAM file")
    parser.add_argument('-1', '--R1', type=str, nargs='?', required=True, help="Output R1 PacBio subreads SAM file")
    parser.add_argument('-2', '--R2', type=str, nargs='?', required=True, help="Output R2 PacBio subreads SAM file")
    parser.add_argument('--other', type=str, nargs='?', required=True, help="Output other PacBio subreads SAM file")
    parser.add_argument('--numcpu', type=int, nargs='?', required=True, help="Number of cpu needed")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    options = Parsers()
    # input file path
    raw_sam_file = options.raw
    # output file path
    R1_sam_file = options.R1
    R2_sam_file = options.R2
    other_sam_file = options.other
    ##
    ## input raw file 
    raw_sam = open(raw_sam_file, 'r')
    ## output results
    R1_sam = open(R1_sam_file, 'w')
    R2_sam = open(R2_sam_file, 'w')
    other_sam = open(other_sam_file, 'w')
    ## get header
    ## NOTE: the header is only subreads form
    header = ""
    i = 1
    while i <= 5:
        header += next(raw_sam)
        i += 1
    ## write header to output file
    R1_sam.write(header)
    R1_sam.flush()
    R2_sam.write(header)
    R2_sam.flush()
    other_sam.write(header)
    other_sam.flush()
    ## start run
    filterPool = pool = Pool(max_workers=options.numcpu)
    for rawline in raw_sam:
        #prf = filterPool.apply_async(runFilter, args=(rawline,), callback=mycallback)
        future = filterPool.submit(runFilter, rawline)
        future.add_done_callback(mycallback)
    # print('waite')
    # filterPool.close()
    # filterPool.join()
    # print('done')
    ## close files
    raw_sam.close()
    R1_sam.close()
    R2_sam.close()
    other_sam.close()
