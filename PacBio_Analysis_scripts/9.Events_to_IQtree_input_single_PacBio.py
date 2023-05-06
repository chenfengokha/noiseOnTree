import sys
import os
import math
import optparse
import pickle
from collections import OrderedDict, defaultdict
import pandas as pd
import subprocess
import time
from ete3 import Tree


def identNONE(edit_list, start):
    if edit_list:
        _, las_end, _, _ = edit_list[-1]
        none_event = (las_end, start, 0, "NONE")
    else:
        none_event = ("1", start, 0, "NONE")
    return none_event

def identEdit(editString):
    start = int(editString.split('+')[1])
    editLen = int(editString.split('+')[0][0:-1])
    end = editLen + int(start)
    if editString.split('+')[0][-1] == "I":
        color = "INSERTION"
        end = start + 1
    elif editString.split('+')[0][-1] == "D":
        color = "DELETION"
    else:
        color = "NONE"
    return start, end, editLen, color

def editAnotation(editline, lenght=447):
    spl = editline.strip().split('\t')
    edit_events = []
    only_edit = sorted(set(spl[0:]), key=spl[0:].index)
    if only_edit == ["NONE"]:
        return [(1, lenght, 0, "NONE")]
    else: 
        if "NONE" in only_edit:
            only_edit.remove("NONE")
        checkList = []
        #print(only_edit)
        for each in only_edit:
            if "&" not in each and "S" not in each and each not in checkList:
                    start, end, editLen, color = identEdit(each)
                    checkList.append(each)
                    edition = (start, end, editLen, color)
                    none_event = identNONE(edit_events, start)
                    edit_events.extend([none_event, edition])
                    if only_edit.index(each) == only_edit.index(only_edit[-1]) and end < lenght:
                        edit_events.append((end, lenght, 0, "NONE"))
            else:
                event_end = 0
                for i in each.split("&"):
                    if "S" not in i and i not in checkList:
                        checkList.append(i)
                        start, end, editLen, color = identEdit(i)
                        event_end = int(end)
                        edition = (start, end, editLen, color)
                        none_event = identNONE(edit_events, start)
                        edit_events.extend([none_event,edition])
                if only_edit.index(each) == only_edit.index(only_edit[-1]) and event_end < lenght:
                    edit_events.append((event_end, lenght, 0, "NONE"))
        return edit_events

def countEvents(event_file):
    next(event_file)
    event_count = OrderedDict()
    for line in event_file:
        spl = line.strip().split('\t')
        #cellNum = int(spl[0].split("_")[1])
        for each in set(spl[1:]):
            if each != "NONE":
                if each in event_count.keys():
                    event_count[each] += 1
                else:
                    event_count[each] = 1
    return event_count

def getEventSet(eventLine):
    eventSet = set()
    for each in eventLine.strip().split('\t')[1:]:
        if "&" in each:
            for every in each.split("&"):
                eventSet.add(every)
        else:
            eventSet.add(each)
    return eventSet

def countEventsFromLine(lineList):
    event_count = OrderedDict()
    for line in lineList:
        eventSet = getEventSet(line)
        for each in eventSet:
            if each != "NONE":
                if each in event_count.keys():
                    event_count[each] += 1
                else:
                    event_count[each] = 1
    return event_count


def calculateWeights(event_count, matrix):
    max_value = event_count[max(event_count,key=event_count.get)]
    # if "NONE" in event_count.keys():
    #     del event_count["NONE"]
    transM = pd.read_csv(matrix, sep='\t', header=0, index_col=0, keep_default_na=False)
    transDict = {}
    for i in transM.index:
        for j in transM.index:
            if transM.loc[i,j] == "" or transM.loc[i,j] == 0 or transM.loc[i,j] == "0":
                continue
            else:
                transDict[int(transM.loc[i,j])] = i+j
    indexRange = list(range(190,0, -1))
    weights = OrderedDict()
    for k, v in event_count.items():
        weight_index = int(round((math.log(v) / math.log(max_value))*(len(transDict) - 1)))
        weight = transDict[indexRange[weight_index]]
        weights[k] = weight
    return weights

def getFilterBClist(filterBC_file):
    if filterBC_file == "All":
        return None
    else:
        BClist = [i.strip() for i in open(filterBC_file, 'r')]
        return BClist

def transTOprotein(eventLine_list, weights, fakePro_file, position_file, alleleInfo_file, edit_fre_file, ref_len=448):
    eventTOunique = defaultdict(list)
    # define WT is or no
    haveWT = False
    for line in eventLine_list:
        spl = line.strip().split('\t')
        BC = spl[0]
        edits = '\t'.join(spl[1:])
        eventTOunique[edits].append(BC)
    # write headers
    fakePro_file.write('%d\t%d\n' % (len(eventTOunique.keys()), len(weights.keys())))
    position_file.write('\t'.join(["y","start", "end", "editLen","event"]) + '\n')
    alleleInfo_file.write("{}\t{}\t{}\t{}\n".format("NodeName", "BC", "cellNum", "\t".join(["target" + str(i) for i in range(1,14)])))
    edit_fre_file.write("{}\t{}\t{}\t{}\t{}\n".format("Position", "Insertion", "InserPercent", "Deletion", "DeletPercent"))
    ## start to get fake protein
    wildtype = ""
    for key, value in weights.items():
        if key == "NONE":
            wildtype += value[0]
        else:
            wildtype += value[1]
    ## build a dict to count the edit frequency in every base
    pos = list(range(1, ref_len))
    Edit_fre = {}
    for p in pos:
        Edit_fre[p] = [0,0]
    ## count cell numbers
    total = 0
    ## get the fake protein sequences
    ind = 1
    for k,v in eventTOunique.items():
        total += len(v)
        fake_pro = ""
        nodeName = "WT"
        for each in weights.keys():
            if each in k:
                fake_pro += weights[each][0]
            else:
                fake_pro += weights[each][1]
        if fake_pro == wildtype:
            haveWT = True
            nodeName = nodeName
        else:
            nodeName = "N" + str(ind) + "_" + str(len(v))
            ind += 1
        ## write out the morphological file to construct tree in iqtree
        pddName = nodeName + " "*(10 - len(nodeName))
        fakePro_file.write(pddName + fake_pro + '\n')
        ## write out the frequency info of each node(alleles)
        for eachBC in v:
            alleleInfo_file.write("{}\t{}\t{}\t{}\n".format(pddName.strip(), eachBC, len(v), k))
        ## write out the plot position of each node
        edit_pos = editAnotation(k)
        for plot in edit_pos:
            editAno = plot
            position_file.write('{}\t{}\t{}\t{}\t{}\n'.format(nodeName, *editAno))
        ## fill in the Edit_fre dict
        for eachEdit in getEventSet(k):
            if "D" in eachEdit:
                start = int(eachEdit.split("+")[1])
                Dlen = int(eachEdit.split("D")[0])
                for ep in range(start, start+Dlen):
                    Edit_fre[ep][1] += len(v)
            if "I" in eachEdit:
                start = int(eachEdit.split("+")[1])
                Edit_fre[start][0] += len(v)
    ## write out frequency
    for positi, freq in Edit_fre.items():
        edit_fre_file.write("{}\t{}\t{}\t{}\t{}\n".format(positi, freq[0], freq[0]/total, freq[1], freq[1]/total))
    ## close and save result
    fakePro_file.close()
    position_file.close()
    alleleInfo_file.close()
    edit_fre_file.close()
    ## retrun haveWT
    return haveWT

###
## run iqtree to get the treefile
###
def exc_IQtree(dir, sequence, outpre, haveWT=False):
    if not os.path.exists(dir):
        os.makedirs(dir)
    os.chdir(dir)
    if haveWT:
        iqtrCMD = "iqtree -s {} -m LG+FO -pre {} -o WT -nt AUTO -czb".format(sequence, outpre)
    else:
        iqtrCMD = "iqtree -s {} -m LG+FO -pre {} -nt AUTO -czb".format(sequence, outpre)
    p = subprocess.Popen(iqtrCMD,stdout=subprocess.PIPE,shell=True)
    _, error = p.communicate()
    return error

def changeTree(inTree, outTree):
    rawTree = Tree(inTree)
    for node in rawTree.iter_descendants():
        node.dist = 1
    rawTree.write(outfile = outTree)

## parameters pass from command line
def Parsers():
    usage = "Usage: python %prog -e [edit events] -o [output morphological] -s [output stastics] -p [output plot positions]"
    parser = optparse.OptionParser(usage=usage)
    group = optparse.OptionGroup(parser, "Transpot edit events to morphological format, and calculate alleles frequency, output the plot position of each event", "This args will define the input output")
    group.add_option("-e", "--events", action = "store", dest = "event_file", type = str, help = "Edit events file")
    group.add_option("-m", "--matrix", action = "store", dest = "transMatrix", type = str, help = "Transform matrix file")
    group.add_option("-d", "--outdir", action = "store", dest = "outdir", type = str, help = "Output directory of output files")
    group.add_option("-n", "--name", action= "store", dest="name", type = str, help = "Output prefix of output files")
    group.add_option("-l", '--len', action="store", dest="refLen", type=int, default=448, help="Reference length, use to plot")
    group.add_option("-f", '--filter', action="store", dest="filterBC", type=str, default="All", help="The BC which wans to keep in the tree, pass a file contain each BC each line, if is All, keep all BCs")
    parser.add_option_group(group)
    options, argvs = parser.parse_args()
    return options, argvs


def main():
    options, _ = Parsers()
    ### input files
    event_file = open(options.event_file, 'r')
    if options.filterBC == "All":
        BClist = False
    else:
        BClist = getFilterBClist(options.filterBC)
    #print(BClist)
    ## get eventLineList
    next(event_file)
    eventLineList = []
    if BClist:
        for line in event_file:
            BC = line.split('\t')[0]
            #print(BC)
            if BC in BClist:
                eventLineList.append(line)
    else:
        for line in event_file:
            eventLineList.append(line)
    #print(eventLineList)
    #eventcountDic = open(options.eventcount, 'rb')
    ### output files
    protein = open(os.path.join(options.outdir, options.name + ".IQtreeInput.txt"), 'w')
    position_file = open(os.path.join(options.outdir, options.name + ".Plots.txt"), 'w')
    alleleInfos = open(os.path.join(options.outdir, options.name + ".AllelesInfo.txt"), 'w')
    editfre_file = open(os.path.join(options.outdir, options.name + ".EditFre.txt"), 'w')
    event_count = countEventsFromLine(eventLineList)
    #print(event_count)
    weights = calculateWeights(event_count, options.transMatrix)
    haveWT = transTOprotein(eventLineList, weights, protein, position_file, alleleInfos, editfre_file, options.refLen)
    ## run iqtree
    iqtree_error = exc_IQtree(options.outdir, os.path.join(options.outdir, options.name + ".IQtreeInput.txt"), options.name, haveWT)
    assert iqtree_error == None, "Runing iqtree Error occured on {}: {}".format(options.name, iqtree_error)
    ## change distance to 1
    rawTree = os.path.join(options.outdir, options.name + ".treefile")
    newTree = os.path.join(options.outdir, options.name + ".nwk")
    changeTree(rawTree, newTree)
    # ## rm the no need files
    # for f in [".bionj", ".ckp.gz", ".iqtree", ".log", ".mldist", ".parstree"]:
    #     file = options.name + f
    #     if os.path.exists(file):
    #         os.remove(os.path.join(options.outdir, file))

if __name__ == "__main__":
    main()

