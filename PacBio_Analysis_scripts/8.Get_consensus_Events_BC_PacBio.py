import random
import os
import sys
import argparse
from collections import defaultdict

def getBCumiDict(rawEventFile):
    BCumi_dict = defaultdict(lambda: defaultdict(int))
    next(rawEventFile)
    for line in rawEventFile:
        spl = line.strip().split("\t")
        BC = spl[0]
        UMI = spl[1]
        events = "\t".join(spl[2:])
        BCumi_dict[BC][events] += 1
    return BCumi_dict

def getMaxVkey(countDict):
    res = []
    for k, v in countDict.items():
        if v == max(countDict.values()):
            res.append(k)
    return res

def getConEvents(EventsCountDict):
    ## 1. choice the max count event
    maxCountEvents = getMaxVkey(EventsCountDict)
    # when only one maxCounEvents, that is the consensus event
    if len(maxCountEvents) == 1:
        return maxCountEvents[0]
    # else random choice
    else:
        # random choice an event no S
        noSevent = [i for i in maxCountEvents if "S" not in i]
        if noSevent:
            conEvent = random.choice(noSevent)
            return conEvent
        # if all have S, random choice one
        else:
            conEvent = random.choice(maxCountEvents)
            return conEvent

def Parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=str, nargs='?', required=True, help="Input the RAW Event file")
    parser.add_argument('-o', '--outfile', type=str, nargs='?', help="Output BC consensus Events file")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    options = Parsers()
    ## input file
    inRawEvent = open(options.infile, 'r')
    ## output file
    outBCconEvent = open(options.outfile, 'w')
    outBCconEvent.write("{}\t{}\n".format("BC", "\t".join(["target" + str(i) for i in range(1, 14)])))
    ## run get consensus events
    BC_umiCountEvents = getBCumiDict(inRawEvent)
    for BC, eventCountDict in BC_umiCountEvents.items():
        conEvent = getConEvents(eventCountDict)
        outBCconEvent.write("{}\t{}\n".format(BC, conEvent))
    ## close and save results
    inRawEvent.close()
    outBCconEvent.close()








