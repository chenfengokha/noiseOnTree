import os
import sys
import operator
import argparse
from collections import defaultdict
import networkx as nx
from itertools import combinations
from copy import deepcopy
sys.path.append("/mnt/data/home/lzz/project/Reconstruct_lineage/scripts/NewReTree/modules")
from SequenceParser import formFASTA

def filterBCs(filter_bc_path):
    bcs = []
    with open(filter_bc_path, 'r') as inf:
        for line in inf:
            bcs.append(line.strip())
    return bcs
# get raw bc and umi info
def getRawBCumi(raw_fa_path):
    rawFA = open(raw_fa_path, 'r').read()
    raw_BC_umi_dict = defaultdict(list)
    for fa_read in formFASTA(rawFA):
        rawname, rawseq = fa_read
        rawBC = rawname.split("BC=")[1].split(" ")[0]
        rawUMI = rawname.split("UMI=")[1].split(" ")[0]
        raw_BC_umi_dict[rawBC].append(rawUMI)
    return raw_BC_umi_dict

# adjust BC , get the rawBC to adjustBC map
def hDistance(seq1, seq2):
    return sum(map(operator.ne, seq1, seq2))

def getFilterBCcluster(rawBClist, filterList):
    BCmap = {}
    for BC in rawBClist:
        for filterBC in filterList:
            if hDistance(BC, filterBC) <= 1:
                BCmap[BC] = filterBC
                #print(BC)
                break
    return BCmap

def changeBCinBCumi(rawBCumi_dict, BCmap):
    adjustBC_rawUMI = defaultdict(list)
    for rawBC, adjustBC in BCmap.items():
        rawUMI_list = rawBCumi_dict[rawBC]
        adjustBC_rawUMI[adjustBC].extend(rawUMI_list)
    return adjustBC_rawUMI

def countList(umiList):
    unique = set(umiList)
    count_dict = {}
    for i in unique:
        count_dict[i] = umiList.count(i)
    return count_dict, list(unique)

def adjustUMI(umiList, dist = 1):
    umi_count, umi_set_list = countList(umiList)
    # 1. compare all umis find singleton and creat Graph to get center
    umiG = nx.Graph()
    singleton_umi = deepcopy(umi_set_list)
    umi_map = {}
    for u1, u2 in combinations(umi_set_list, 2):
        if hDistance(u1, u2) <= 1:
            umiG.add_edges_from([(u1,u2)])
            if u1 in singleton_umi:
                singleton_umi.remove(u1)
            if u2 in singleton_umi:
                singleton_umi.remove(u2)
    if singleton_umi:
        for eachSgl in singleton_umi:
            umi_map[eachSgl] = eachSgl
    if list(nx.connected_components(umiG)):
        for eachG in nx.connected_component_subgraphs(umiG):
            centerUmilist = nx.center(eachG)
            if len(centerUmilist) == 1:
                for rawUmi in list(nx.connected_components(eachG))[0]:
                    umi_map[rawUmi] = centerUmilist[0]
            else:
                max_umiCount = 0
                centerUmi = ""
                for cenUmi in centerUmilist:
                    if umi_count[cenUmi] > max_umiCount:
                        max_umiCount = umi_count[cenUmi]
                        centerUmi = cenUmi
                for rawUmi in list(nx.connected_components(eachG))[0]:
                    umi_map[rawUmi] = centerUmi
    return umi_map

def getAdjustBC_UMImap_dict(adjustBCrawUMI):
    adjustBCandUMI = {}
    for adBC, rawUMI_list in adjustBCrawUMI.items():
        umiMap_dict = adjustUMI(rawUMI_list)
        adjustBCandUMI[adBC] = umiMap_dict
    return adjustBCandUMI

def changeBCandUMIinRawFA(raw_fa_path, adjust_fa_path, adBCumiMap, BCmap):
    rawBCfa = open(raw_fa_path, 'r').read()
    adBCumifa = open(adjust_fa_path, 'w')
    for record in formFASTA(rawBCfa):
        name, seq = record
        rawBC = name.split("BC=")[1].split(" ")[0]
        if rawBC in BCmap:
            adBC = BCmap[rawBC]
            rawUMI = name.split("UMI=")[1].split(" ")[0]
            adUMI = adBCumiMap[adBC][rawUMI]
            # construct new name
            prefix = name.split(rawBC)[0]
            suffix = name.split(rawUMI)[1]
            new_name = prefix + adBC + " UMI=" + adUMI + suffix
            adBCumifa.write(">{}\n{}\n".format(new_name, seq))
    adBCumifa.close()



def Parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--rawFA', type=str, nargs='?', required=True, help="Input raw fasta file")
    parser.add_argument('-o', '--outFA', type=str, nargs='?', required=True, help="Output adjust BC and UMI file")
    #parser.add_argument('-f', '--filterBC', type=str, nargs='?', required=True, help="filter BCs from scRNA-seq")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    options = Parsers()
    # 1. get filter BC list
    #filter_BC_list = filterBCs(options.filterBC)
    # 2. get rawBC umi dict
    raw_BCumi = getRawBCumi(options.rawFA)
    raw_BC_list = list(raw_BCumi.keys())
    # 3. adjust BC
    #BCmap = getFilterBCcluster(raw_BC_list, filter_BC_list)
    BCmap = adjustUMI(raw_BC_list, dist=1)
    # 4. change rawBC in BCumi dict
    adjustBC_rawUMI_dict = changeBCinBCumi(raw_BCumi, BCmap)
    # 5. adjust rawUMI in each BC
    adjustBC_adUMImap = getAdjustBC_UMImap_dict(adjustBC_rawUMI_dict)
    # 6. change raw fasta file
    changeBCandUMIinRawFA(options.rawFA, options.outFA, adjustBC_adUMImap, BCmap)
