import os
import subprocess
import argparse
import sys
import time
from ete3 import Tree

def Parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--prefix", type=str, help="Give the output prefix", required=True)
    parser.add_argument("-d","--directory", type=str, help="Set up an directory to run mix", required=True)
    parser.add_argument("-i","--input", type=str, help="The sequence of iqtree input", required=True)
    args = parser.parse_args()
    return args

def exc_IQtree(dir, sequence, outpre):
    if not os.path.exists(dir):
        os.makedirs(dir)
    os.chdir(dir)
    iqtrCMD = "iqtree -s {} -m LG+FO -pre {} -nt AUTO".format(sequence, outpre)
    p = subprocess.Popen(iqtrCMD,stdout=subprocess.PIPE,shell=True)
    _, error = p.communicate()
    return error

def changeTree(inTree, outTree):
    rawTree = Tree(inTree)
    for node in rawTree.iter_descendants():
        node.dist = 1
    rawTree.write(outfile = outTree)

if __name__ == "__main__":
    options = Parsers()
    start = time.time()
    ## run iqtree
    iqtree_error = exc_IQtree(options.directory, options.input, options.prefix)
    assert iqtree_error == None, "Runing iqtree Error occured on {}: {}".format(options.prefix, iqtree_error)
    ## change distance to 1
    rawTree = os.path.join(options.directory, options.prefix + ".treefile")
    newTree = os.path.join(options.directory, options.prefix + ".nwk")
    changeTree(rawTree, newTree)
    end = time.time()
    print("{}\t{}\n".format(options.prefix,end - start))
    
