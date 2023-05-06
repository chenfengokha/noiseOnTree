import sys
import os
import math
import optparse
from collections import OrderedDict, defaultdict
import subprocess
import time
from ete3 import Tree
## import self define modules
sys.path.append("/mnt/data/home/lzz/project/Reconstruct_lineage/scripts/NewReTree")
from modules.FIFO import FIFO

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
                for i in each.split("&"):
                    if "S" not in i and i not in checkList:
                        checkList.append(i)
                        start, end, editLen, color = identEdit(i)
                        edition = (start, end, editLen, color)
                        none_event = identNONE(edit_events, start)
                        edit_events.extend([none_event,edition])
                if only_edit.index(each) == only_edit.index(only_edit[-1]) and list(edit_events[-1])[1] < lenght:
                    edit_events.append((list(edit_events[-1])[1], lenght, 0, "NONE"))
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

def calculateWeights(event_count, weights_file):
    max_value = event_count[max(event_count,key=event_count.get)]
    if "NONE" in event_count.keys():
        del event_count["NONE"]
    characterArray = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    weights = []
    for v in event_count.values():
        weight = int(round((math.log(v) / math.log(max_value))*(len(characterArray) - 1)))
        weights.append(characterArray[weight])
    weights_file.write(''.join(weights))
    weights_file.close()

def getFilterBClist(filterBC_file):
    if filterBC_file == "All":
        return None
    else:
        BClist = [i.strip() for i in open(filterBC_file, 'r')]
        return BClist

def transTObinary(eventLine_list, event_count, binary_events, position_file, alleleInfo_file, edit_fre_file, ref_len=448):
    if "NONE" in event_count.keys():
        del event_count["NONE"]
    eventTOunique = defaultdict(list)
    # define WT is or no
    haveWT = False
    for line in eventLine_list:
        spl = line.strip().split('\t')
        BC = spl[0]
        edits = '\t'.join(spl[1:])
        eventTOunique[edits].append(BC)
    # write headers
    binary_events.write('%d\t%d\n' % (len(eventTOunique.keys()), len(event_count.keys())))
    position_file.write('\t'.join(["y","start", "end", "editLen","event"]) + '\n')
    alleleInfo_file.write("{}\t{}\t{}\t{}\n".format("NodeName", "BC", "cellNum", "\t".join(["target" + str(i) for i in range(1,14)])))
    edit_fre_file.write("{}\t{}\t{}\t{}\t{}\n".format("Position", "Insertion", "InserPercent", "Deletion", "DeletPercent"))
    ## start to get binary seq
    wildtype = "0" * len(event_count.keys())
    ## build a dict to count the edit frequency in every base
    pos = list(range(1, ref_len))
    Edit_fre = {}
    for p in pos:
        Edit_fre[p] = [0,0]
    ## count cell numbers
    total = 0
    ## get the binary seq
    ind = 1
    for k,v in eventTOunique.items():
        total += len(v)
        binary = ""
        nodeName = "WT"
        for each in event_count.keys():
            if each in k:
                binary += "1"
            else:
                binary += "0"
        if binary == wildtype:
            haveWT = True
            nodeName = nodeName
        else:
            nodeName = "N" + str(ind) + "_" + str(len(v))
            ind += 1
        ## write out the morphological file to construct tree in MIX
        pddName = nodeName + " "*(10 - len(nodeName))
        binary_events.write(pddName + binary + '\n')
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
    binary_events.close()
    position_file.close()
    alleleInfo_file.close()
    edit_fre_file.close()
    ## retrun haveWT
    return haveWT, wildtype

###
## run MIX to get the treefile
###
def exc_MIX(dir, binary, weight, haveWT, wildtype):
    if not os.path.exists(dir):
        os.makedirs(dir)
    os.chdir(dir)
    if not haveWT:
        args=["P","W","4","5","Y"]
        params = '\n'.join(args)
        with FIFO(dir=dir, name="mix_args.txt") as Argfile:
            with open(Argfile.filename, 'w') as Args:
                all_cmd = "{}\n{}\n{}\n".format(binary, params, weight)
                Args.write(all_cmd)
            mixCMD = "mix < " + Argfile.filename
            p = subprocess.Popen(mixCMD,stdout=subprocess.PIPE,shell=True)
            _, error = p.communicate()
            return(error)
    else:
        args=["P","A","W","4","5","Y"]
        with open(os.path.join(dir, "ancestor.txt"), 'w') as ancestor:
            ancestor.write(wildtype + "\n")
        params = '\n'.join(args)
        with FIFO(dir=dir, name="mix_args.txt") as Argfile:
            with open(Argfile.filename, 'w') as Args:
                all_cmd = "{}\n{}\n{}\n{}\n".format(binary, params, weight, os.path.join(dir, "ancestor.txt"))
                Args.write(all_cmd)
            mixCMD = "mix < " + Argfile.filename
            p = subprocess.Popen(mixCMD,stdout=subprocess.PIPE,shell=True)
            _, error = p.communicate()
            return error


def rename_file(dir, name, program):
    os.chdir(dir)
    os.rename("outtree", name + program + "Tree")
    os.rename("outfile", name + program + "Outfile")
    return os.path.join(dir, name + program + "Tree")
    

def run_consense(dir, mix_outtree):
    if not os.path.exists(dir):
        os.makedirs(dir)
    os.chdir(dir)
    with FIFO(dir=dir, name="consense_args.txt") as conArg:
        with open(conArg.filename, 'w') as conArgs:
            conArgs.write('{}\n{}\n'.format(mix_outtree, "Y"))
        conCMD = "consense < " + conArg.filename
        q = subprocess.Popen(conCMD, stdout=subprocess.PIPE, shell=True)
        _, con_error = q.communicate()
        return con_error

## parameters pass from command line
def Parsers():
    usage = "Usage: python %prog -e [edit events] -o [output morphological] -s [output stastics] -p [output plot positions]"
    parser = optparse.OptionParser(usage=usage)
    group = optparse.OptionGroup(parser, "Transpot edit events to morphological format, and calculate alleles frequency, output the plot position of each event", "This args will define the input output")
    group.add_option("-e", "--events", action = "store", dest = "event_file", type = str, help = "Edit events file")
    #group.add_option("-m", "--matrix", action = "store", dest = "transMatrix", type = str, help = "Transform matrix file")
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
    binary_file = open(os.path.join(options.outdir, options.name + ".MixInput.txt"), 'w')
    mix_weight_file = open(os.path.join(options.outdir, options.name + ".MixWeight.txt"), 'w')
    position_file = open(os.path.join(options.outdir, options.name + ".MixPlots.txt"), 'w')
    alleleInfos = open(os.path.join(options.outdir, options.name + ".MixAllelesInfo.txt"), 'w')
    editfre_file = open(os.path.join(options.outdir, options.name + ".MixEditFre.txt"), 'w')
    event_count = countEventsFromLine(eventLineList)
    #print(event_count)
    calculateWeights(event_count, mix_weight_file)
    haveWT, wildtype = transTObinary(eventLineList, event_count, binary_file, position_file, alleleInfos, editfre_file, options.refLen)
    ## run mix
    mix_error = exc_MIX(options.outdir, os.path.join(options.outdir, options.name + ".MixInput.txt"), os.path.join(options.outdir, options.name + ".MixWeight.txt"), haveWT, wildtype)
    assert mix_error == None, "Runing mix. Error occured on {}: {}".format(options.name, mix_error)
    mix_tree = rename_file(options.outdir, options.name, "MIX")
    con_error = run_consense(options.outdir, mix_tree)
    assert con_error == None, "Runing consense. Error occured on {}: {}".format(options.name, con_error)
    con_tree = rename_file(options.outdir, options.name, "consense")

if __name__ == "__main__":
    main()