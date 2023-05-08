import sys
import os
import optparse
import multiprocessing
import subprocess
from collections import OrderedDict, defaultdict
import pandas as pd
## import self define modules
# sys.path.append("/mnt/data/home/lzz/project/Reconstruct_lineage/scripts/NewReTree")
# from modules.FIFO import FIFO

def seq_dist(name1, name2, seq1, seq2, transM):
    seq_len = max(len(seq1), len(seq2))
    all_dist = 0
    for i,j in zip(seq1, seq2):
        if i == j:
            sub_dist = 0
        elif transM.loc[i,j] == "":
            sub_dist = 190 - int(transM.loc[j,i])
        else:
            sub_dist = 190 - int(transM.loc[i,j])
        all_dist += sub_dist
    return name1, name2, all_dist / float(seq_len)

def get_seqs(file_path):
    protein_seq = OrderedDict()
    with open(file_path, 'r') as inf:
        next(inf)
        for line in inf:
            spl = line.strip().split(" ")
            name = spl[0]
            pro_seq = spl[-1]
            protein_seq[name] = pro_seq
    return protein_seq


###
# ## run neighbor to get the treefile
# ###
def exc_Neighbor(dir, dist_file):
    if not os.path.exists(dir):
        os.makedirs(dir)
    os.chdir(dir)
    args=["L","Y"]
    params = '\n'.join(args)
    with FIFO(dir=dir, name="neighbor_args.txt") as Argfile:
        with open(Argfile.filename, 'w') as Args:
            all_cmd = "{}\n{}\n".format(dist_file,params)
            Args.write(all_cmd)
        neighborCMD = "neighbor < " + Argfile.filename
        p = subprocess.Popen(neighborCMD,stdout=subprocess.PIPE,shell=True)
        _, error = p.communicate()
        return error



## parameters pass from command line
def Parsers():
    usage = "Usage: python %prog -p [protein] -m [trans matrix] -o [output plot positions]"
    parser = optparse.OptionParser(usage=usage)
    group = optparse.OptionGroup(parser, "Construct tree from IQtree input pseudo protein sequences", "This args will define the input output")
    group.add_option("-p", "--protein", action = "store", dest = "protein_file", type = str, help = "Protein sequence file")
    group.add_option("-m", "--matrix", action = "store", dest = "transMatrix", type = str, help = "Transform matrix file")
    group.add_option("-o", "--outfile", action = "store", dest = "outfile", type = str, help = "Output path of output files")
    group.add_option("-d", "--dir", action = "store", dest = "out_dir", type = str, help = "Output path of output directory")
    group.add_option("-n", "--cpu", action = "store", dest = "ncpu", type = int, help = "Number of CPU to use")
    parser.add_option_group(group)
    options, argvs = parser.parse_args()
    return options, argvs

def main():
    options, _ = Parsers()
    ### input files
    pro_seqs = get_seqs(options.protein_file)
    all_seq_names = list(pro_seqs.keys())
    transMX = pd.read_csv(options.transMatrix,sep='\t', header=0, index_col=0, keep_default_na=False)
    ##
    processPool = multiprocessing.Pool(options.ncpu)
    results = []
    for n1 in all_seq_names:
        for n2 in all_seq_names:
            dist_re = processPool.apply_async(seq_dist, args = (n1, n2, pro_seqs[n1], pro_seqs[n2], transMX, ))
            results.append(dist_re)
    print('waite')
    processPool.close()
    processPool.join()
    print('done')
    # extract results
    all_seq_dist = defaultdict(dict)
    for sub_re in results:
        sub_dist = sub_re.get()
        seq1_name, seq2_name, compDist = sub_dist
        all_seq_dist[seq1_name][seq2_name] = compDist
    # output results
    with open(options.outfile, 'w') as outf:
        outf.write("\t" + str(len(all_seq_names)) + "\n")
        for each_name in all_seq_names:
            outf.write("{}\t{}\n".format(each_name, "\t".join([str(all_seq_dist[each_name][i]) for i in all_seq_names])))
    ## exec neighbor
    exc_Neighbor(options.dir, options.outfile)
    
if __name__ == "__main__":
    main()
