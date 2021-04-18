import os
import sys
import argparse
import tempfile
import shutil
import multiprocessing
import time
## import self define modules
sys.path.append("/mnt/data/home/lzz/project/Reconstruct_lineage/scripts/NewReTree")
from modules.FIFO import FIFO
from modules.SequenceParser import getReference


def Parsers():
    parser = argparse.ArgumentParser(description="Align reads to reference use MUSCLE.")
    parser.add_argument('-r', '--reference', type=str, nargs='?', required=True, help="References file")
    parser.add_argument('-s', '--seq', type=str, nargs='?',required=True, help="Sequence file")
    parser.add_argument('-o', '--out', type=str, nargs='?', required=True, help="Out put FASTA format alignment file")
    parser.add_argument('-c', '--cpu', type=int,  nargs='?',required=True, help="The Number of cpu needed")
    args = parser.parse_args()
    return args

def getReads(seq):
    while True:
        name = next(seq)
        record = next(seq)
        seq_record = name + record
        yield seq_record
        if not name:
            break

def getpath(filename):
    path = os.path.dirname(filename)
    return path

### run align with muscle
def runMuscle(ref, seq, outdir):
    alignResult = ""
    with FIFO(dir=outdir, name="tempin.fa") as temp_in, FIFO(dir=outdir, name="tempout.fa") as temp_out:
        with open(temp_in.filename, 'w') as tempin:
            tempin.write(ref + "\n" + seq + "\n")
        os.system("muscle -in %s -out %s"%(temp_in.filename,temp_out.filename))
        with open(temp_out.filename, 'r') as tempout:
            for each in tempout:
                alignResult += each        
    return alignResult

def mycallback(results):
    alignf.write(results)
    alignf.flush()

if __name__ == "__main__":
    start = time.time()
    options = Parsers()
    seq = open(options.seq,'r')
    ref = open(options.reference,'r')
    outdir = getpath(options.out)
    reference = getReference(ref)
    alignf = open(options.out, 'w')
    pfq = multiprocessing.Pool(options.cpu)
    results = []
    for each in getReads(seq):
        prf = pfq.apply_async(runMuscle, args=(reference, each, outdir,), callback=mycallback)
        #results.append(prf)
    print('waite')
    pfq.close()
    pfq.join()
    print('done')
    #for i in results:
        #alignf.write(i.get())
    alignf.close()
    end = time.time()
    print(end - start)



