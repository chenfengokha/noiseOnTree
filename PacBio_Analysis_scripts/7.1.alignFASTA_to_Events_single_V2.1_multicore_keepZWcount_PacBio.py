import collections
import sys
import os
import multiprocessing
import argparse
from collections import defaultdict
import pickle
## import self define modules
sys.path.append("/mnt/data/home/lzz/project/Reconstruct_lineage/scripts/NewReTree")
from modules.SequenceParser import formFASTA
from modules.Alignment import Alignment
### build a class to hanndle alignment object
### the composition of the Alignment class is every Position of reference, the base in the same
### position of reference and reads.

### filter the Aligment List
def filterEnds(eventlist, minMatchBase = 8, matchPercent = 0.6):
    startIndex = -1
    endIndex = -1
    for Findex in range(len(eventlist)):
        Fevent = eventlist[Findex]
        if startIndex < 0 and \
            Fevent.indicator == "M" and \
            Fevent.matchPct() > matchPercent and \
            (len(Fevent.readbase) - Fevent.refbase.count("N")) >= minMatchBase:
            startIndex = Findex
    for Rindex in range((len(eventlist) - 1), -1, -1):
        Revent = eventlist[Rindex]
        if endIndex < 0 and \
            Revent.indicator == "M" and \
            Revent.matchPct() > matchPercent and \
            (len(Revent.readbase) - Revent.refbase.count("N")) >= minMatchBase:
            endIndex = Rindex
    return eventlist[startIndex: endIndex + 1]

### drop out the reads which have great deletions


### Find the Scar("S") events 
def findScar(refseq, readseq, minScarSize = 3, maxScarSize = 10, minScarProp = 0.75):
    if len(refseq) != len(readseq):
        raise Exception("Error: The lengths of reference and read is different, read: {}, reference:{}".format(len(readseq),len(refseq)))
    scarList = []
    #windowScale = (minScarSize, maxScarSize)
    i = 0   
    while i <= len(refseq) - minScarSize:
        maxPosition = -1
        maxMismatches = 0
        maxMismatchProp = 0.0
        windowSize = 1
        if i <= len(refseq) - maxScarSize:
            windowScale = list(range(minScarSize, maxScarSize + 1))
        else: 
            windowScale = list(range(minScarSize, len(refseq)-i+1))
        for j in windowScale:
            miscount = 0
            for b1, b2 in zip(refseq[i: i+j].upper(), readseq[i: i+j].upper()):
                if b1 != b2 and b1 != "-" and b2 != "-":
                    miscount += 1
            misProp = miscount / j
            if miscount >= maxMismatches and misProp >= minScarProp:
                maxPosition = i
                maxMismatchProp = misProp
                windowSize = j
                maxMismatches = miscount
        if maxMismatches >= minScarSize:
            scarList.append(Alignment(maxPosition, refseq[maxPosition: maxPosition + windowSize], readseq[maxPosition: maxPosition + windowSize], "S"))
        i += windowSize
    if scarList:
        return scarList
    else:
        return None


### Analysis the editevents. "D" represent deletion; "I" represent insertion; "M" represent match or mismatch
### Iterate the aligment, for each base, can be a Alignment Class, return the list of edit events
### refseq = reference read.  readseq = the read record
def analyzeEdit(readname,refseq, readseq, miniMatchEnd = 8):
    refpos = 1
    editEvents = []
    for line in zip(refseq,readseq):
        refbase, readbase = line
        if (refbase, readbase) == (refbase, "-"):         
            if editEvents:
                newEvent = editEvents.pop().combine(Alignment(refpos, refbase, readbase, "D"))
                editEvents.extend(newEvent)
            else:
                editEvents.append(Alignment(refpos, refbase, readbase, "D"))
            refpos += 1
        elif (refbase, readbase) == ("-", readbase):
            if editEvents:
                newEvent = editEvents.pop().combine(Alignment(refpos, refbase, readbase, "I"))
                editEvents.extend(newEvent)
            else:
                editEvents.append(Alignment(refpos, refbase, readbase, "I"))
        else:
            if editEvents:
                newEvent = editEvents.pop().combine(Alignment(refpos, refbase, readbase, "M"))
                editEvents.extend(newEvent)
            else:
                editEvents.append(Alignment(refpos, refbase, readbase, "M"))            
            refpos += 1
    #editEvents = filterEnds(editEvents,miniMatchEnd)
    scarList = findScar(refseq, readseq)
    if scarList:
        editEvents.extend(scarList)
    editEvents.append(readname)
    return editEvents

# def findalledit(align_fa, miniMend = 8):
#     editDict = collections.OrderedDict()
#     all_seq = formFASTA(align_fa)
#     while True:
#         try:
#             _, refseq = next(all_seq)
#             readname, readseq = next(all_seq)
#             editDict[readname] = analyzeEdit(refseq, readseq, miniMend)
#         except StopIteration:
#             break
#     return editDict

def overlap(start1, end1, start2, end2):
    return end1 >= start2 and end2 >= start1

def callToEvents(start, end, editlist):
    targetEvents = ""
    candidates = []
    for each in editlist:
        if each.indicator == "I" and overlap(start, end, each.refpos, each.refpos):
            candidates.append(each)
        elif each.indicator == "D" and overlap(start, end, each.refpos, each.refpos + len(each.refbase) -1):
            candidates.append(each)
        elif each.indicator == "S" and overlap(start, end, each.refpos, each.refpos + len(each.readbase) -1):
            candidates.append(each)
    if len(candidates) == 0:
        targetEvents = "NONE"
    elif len(candidates) == 1:
        targetEvents = candidates[-1].toEditString()
    else:
        targetEvents = "&".join([i.toEditString() for i in candidates])
    return targetEvents

def outputEdit(events, targetPosition):
    eventList = []
    for tu in targetPosition:
        start, end = tu
        eventList.append(callToEvents(start, end, events))
    return '\t'.join(eventList)

### filter out have number of threshold scars
# def filterScar(events, num):
#     Scount = 0
#     for s in events:
#         if s.indicator == "S":
#             Scount += 1
#     return Scount < num

# def filterDeletion(events, threshold, mode):
#     Dlen = []
#     eCount = 0
#     for i in events:
#         if i.indicator != "M":
#             eCount += 1
#             if i.indicator == "D":
#                 Dlen.append(len(i.readbase))
#     if Dlen:
#         maxD = max(Dlen)
#         if mode == 'single':
#             if eCount == 1:
#                 return maxD < threshold
#             else: return True
#         elif mode == 'multi':
#             return maxD < threshold
#     else: return True

# def filterEdit(editDict, maxscar, maxdele, delemode):
#     new_editDict = {}
#     for k,v in editDict.items():
#         if maxscar > 0 and maxdele >0:
#             if filterScar(v, maxscar) and filterDeletion(v, maxdele, delemode):
#                 new_editDict[k] = v
#         elif maxscar > 0 and maxdele == 0:
#             if filterScar(v, maxscar):
#                 new_editDict[k] = v
#         elif maxscar == 0 and maxdele > 0:
#             if filterDeletion(v, maxdele, delemode):
#                 new_editDict[k] = v
#         else:
#             new_editDict[k] = v
#     return new_editDict

def Parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--align', type=str, nargs='?', required=True, help="Input the FASTA format align file")
    parser.add_argument('-t', '--targetpos', type=str, nargs='?', required=True, help="Positions of each target in reference, including start end")
    parser.add_argument('-n', '--name', type=str,  nargs='?',required=True, help="The prefix of the output file")
    parser.add_argument('--dir', type=str,  nargs='?',required=True, help="Output diretory ")
    #parser.add_argument('--len', type=int,  nargs='?',required=True, default=448, help="Reference length, use to plot")
    parser.add_argument('--cpu', type=int,  nargs='?',default=1, help="Number of cpu wants to use")
    args = parser.parse_args()
    return args


def main():
    options = Parsers()
    
    alignFa = open(options.align, 'r').read()
    targetPos = open(options.targetpos, 'r')
    ### output file
    output = open(os.path.join(options.dir, options.name + ".Events.txt"), 'w')
    # editfre = open(os.path.join(options.dir, options.name + ".EditFrequency.txt"), 'w')
    # eventcount = open(os.path.join(options.dir, options.name + ".EventCounts.txt"), 'w')
    # eventcountDic = open(os.path.join(options.dir, options.name + ".EventCounts.pkl"), "wb")

    next(targetPos)
    tarPos = []
    header = []
    ### parser target positions 
    for eachline in targetPos:
        sple = eachline.strip().split('\t')
        header.append(sple[0])
        tarPos.append((int(sple[1]),int(sple[2])))
    ### build a dict to count the edit frequency in every base
    # pos = list(range(1, options.len))
    # Edit_fre = {}
    # for p in pos:
    #     Edit_fre[p] = [0,0]
    ### write the output header
    #editfre.write("{}\t{}\t{}\t{}\t{}\n".format("Position", "Insertion", "InserPercent", "Deletion", "DeletPercent"))
    output.write("{}\t{}\t{}\t{}\n".format("BC", "UMI", "ZWcount", '\t'.join(header)))
    #editfre.flush()
    output.flush()
    ### get the output data
    processPool = multiprocessing.Pool(options.cpu)
    results = []
    all_seq = formFASTA(alignFa)
    while True:
        try:
            _, refseq = next(all_seq)
            readname, readseq = next(all_seq)
            editresult = processPool.apply_async(analyzeEdit, args = (readname,refseq, readseq,))
            results.append(editresult)
        except StopIteration:
            break
    print('waite')
    processPool.close()
    processPool.join()
    print('done')
    editEventDict = collections.OrderedDict()
    for each_result in results:
        rawEditlist = each_result.get()
        editEventDict[rawEditlist[-1]] = rawEditlist[:-1]
    #EventCounts = defaultdict(int)
    #total = 0
    for key, value in editEventDict.items():
        #frequency = int(key.strip().split('_')[-1])
        BC = key.split("BC=")[1].split(" ")[0]
        UMI = key.split("UMI=")[1].split(" ")[0]
        zwCount = key.split("zwCCS.num=")[1].split(" ")[0]
        stringline = outputEdit(value, tarPos)
        ### write out edit events
        output.write("{}\t{}\t{}\t{}\n".format(BC, UMI, zwCount, stringline))
        #total += 1
        # for align in value:
        #     EventCounts[align.toEditString()] += 1
        #     if align.indicator == "D":
        #         for i in range(align.refpos, align.refpos + len(align.readbase)):
        #             Edit_fre[i][1] += 1
        #     elif align.indicator == "I":
        #         Edit_fre[align.refpos][0] += 1
    # ### write out the edit frequency
    # for k, v in Edit_fre.items():
    #     editfre.write("{}\t{}\t{}\t{}\t{}\n".format(k, v[0], v[0]/total, v[1], v[1]/total))
    ### write out the Edit events counts
    # eventcount.write("{}\t{}\n".format("Editevent", "counts"))
    # clean_EventCounts = {}
    # for eve, num in EventCounts.items():
    #     if "M" not in eve:
    #         clean_EventCounts[eve] = num
    #         eventcount.write("{}\t{}\n".format(eve, num))
    # pickle.dump(clean_EventCounts,eventcountDic)


    #editfre.close()
    output.close()
    #eventcount.close()
    #eventcountDic.close()
if __name__ == "__main__":
    main()


    


        

