import os
import sys
import tempfile
import shutil
from collections import OrderedDict, defaultdict
import subprocess
import argparse
import multiprocessing
sys.path.append("/mnt/data/home/lzz/project/Reconstruct_lineage/scripts/NewReTree")
from modules.SequenceParser import convert2string, fqTOfa, formFASTA, seqco, seqre
from modules.FIFO import FIFO

def group_reads(fa_file):
    BC_UMI_seq = defaultdict(lambda: defaultdict(list))
    for record in formFASTA(fa_file):
        name, seq = record
        BC = name.split("BC=")[1].split(" ")[0]
        UMI = name.split("UMI=")[1].split(" ")[0]
        BC_UMI_seq[BC][UMI].append(seq)
    return BC_UMI_seq

def call_consensus(alignment, minproportion = 0.5):
    pos_base = OrderedDict()
    record_size = 0
    for faObj in formFASTA(alignment):
        name, seq = faObj
        reads_count = 1
        record_size += reads_count
        index = 1
        for base in seq:
            index += 1
            if index not in pos_base.keys():
                pos_base[index] = defaultdict(int)
            if base != "-":
                pos_base[index][base] += reads_count
    
    con_seq = ""
    for pos, count in pos_base.items():
        total = float(sum(count.values()))
        if total / record_size >= minproportion:
            maxnum = 0
            maxbase = ""
            for base, num in count.items():
                if num > maxnum:
                    maxnum = num
                    maxbase = base
            con_seq += maxbase
    return con_seq, record_size
            


def runMuscleMultipalign(barcode, umi, seqObj_list, temp_dir):
    con_seq = ""
    read_num = 0
    ind = 1
    with FIFO(dir=temp_dir, name="tempin.fa") as temp_in, FIFO(dir=temp_dir, name="tempout.fa") as temp_out:
        with open(temp_in.filename, 'w') as tempin:
            for each in seqObj_list:
                tempin.write(">{}\n{}\n".format(umi + str(ind), each))
                ind += 1
        CMD = "muscle -in %s -out %s" % (temp_in.filename, temp_out.filename)
        p = subprocess.Popen(CMD,stdout=subprocess.PIPE,shell=True)
        _, error = p.communicate()
        with open(temp_out.filename, 'r') as tempout:
            con_seq, read_num = call_consensus(tempout.read())
    return barcode, umi, con_seq, read_num


def mycallback(results):
    barcode, umi, con_sequence, reads_num = results
    new_header = ">con_seq BC={} UMI={} zwCCS.num={} length={}".format( barcode, umi, reads_num, len(con_sequence))
    new_record = "{}\n{}\n".format(new_header, con_sequence)
    out_fq.write(new_record)
    out_fq.flush()

def Parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=str, nargs='?', required=True, help="Input the FASTA format file")
    #parser.add_argument('--readtype', type=str, nargs='?', required=True, default= 'R1', choices=['R1', 'R2'], help="Direction of FASTQ file, choices R1 or R2")
    parser.add_argument('--numcpu', type=int, nargs='?', required=True, help="Number of cpu needed")
    parser.add_argument('-o', '--outfile', type=str, nargs='?', help="The out file of results")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    options = Parsers()

    in_fq = open(options.infile, 'r').read()
    out_fq = open(options.outfile, 'w')
    outdir = os.path.dirname(options.outfile)
    #direction = options.readtype

    bc_umi_reads = group_reads(in_fq)

    pfq = multiprocessing.Pool(options.numcpu)
    for bc, umi_dict in bc_umi_reads.items():
        for k, v in umi_dict.items():
            prf = pfq.apply_async(runMuscleMultipalign, args=(bc, k, v, outdir,), callback=mycallback)
    
    print('waite')
    pfq.close()
    pfq.join()
    print('done')
    out_fq.close()

    