import os
import sys
import argparse
from collections import defaultdict, Counter
from copy import deepcopy

def getFilterBClist(filterBC_file):
    if filterBC_file == "All":
        return None
    else:
        BClist = [i.strip() for i in open(filterBC_file, 'r')]
        return BClist


def getBCumiDict(rawEventFile, filterBC):
    BCumi_dict = defaultdict(lambda: defaultdict(list))
    next(rawEventFile)
    for line in rawEventFile:
        spl = line.strip().split("\t")
        BC = spl[0]
        UMI = spl[1]
        zwCount = spl[2]
        events = "\t".join(spl[3:])
        if filterBC != "no":
            if BC in filterBC:
                BCumi_dict[BC][UMI] = [events, zwCount]
        else:
            BCumi_dict[BC][UMI] = [events, zwCount]
    return BCumi_dict

def getMaxVkey(countDict):
    res = []
    for k, v in countDict.items():
        if v == max(countDict.values()):
            res.append(k)
    return res

### score edit 
def scoreUMIandEdit(umi_edit_dict):
    ## 1. count event in every target
    target_edit_count = defaultdict(int)
    for umi, infoList in umi_edit_dict.items():
        eString = infoList[0]
        for event in set(list(eString.strip().split('\t'))):
            if "&" in event:
                for each in event.split("&"):
                    target_edit_count[each] += 1
            else:
                target_edit_count[event] += 1
    ## 2. score the umi events base on the count
    umi_score_dict = {}
    edit_score_dict = {}
    for u, il in umi_edit_dict.items():
        score = 0
        es = il[0]
        for e in set(list(es.strip().split('\t'))):
            if e == "NONE":
                score += 0
            elif "&" in e:
                for ev in e.split("&"):
                    score += target_edit_count[ev]
            else:
                score += target_edit_count[e]
        umi_score_dict[u] = score
        edit_score_dict[es] = score
    return umi_score_dict, edit_score_dict

def getEventCount(umiInfo_dict):
    eventCount_ = defaultdict(int)
    for umi, infl in umiInfo_dict.items():
        es = infl[0]
        eventCount_[es] += 1
    return eventCount_

def checkInfo(conEvent, checkEvent):
    isCon = False
    discard = False
    if checkEvent != conEvent:
        for ec_, e in zip(conEvent.strip().split('\t'), checkEvent.strip().split('\t')):
            if e != "NONE" and e not in ec_:
                discard = True
    else:
        isCon = True
    if isCon == True:
        return "con"
    else:
        if discard == True:
            return "dis"
        else:
            return "anc"

def getRawConEventCount_comSCOREandUMI(BCumiInfo_dict):
    bc_conE_count = defaultdict(int)
    for B_, u_info in BCumiInfo_dict.items():
        e_l = [u_info[i][0] for i in u_info.keys()]
        # 1. score edit
        _, e_s = scoreUMIandEdit(u_info)
        # 2. count edit
        e_c = getEventCount(u_info)
        # 3. change edit base on score
        for e1_ in set(e_l):
            for e_ in e_l:
                if checkInfo(e1_, e_) == "anc":
                    e_c[e1_] += 1
        # 4. get max count edit
        mce_l = getMaxVkey(e_c)
        # 5. if more than 1 max count edit select the high score
        maxS_e = {}
        for mce_ in mce_l:
            maxS_e[mce_] = e_s[mce_]
        maxSC_l = getMaxVkey(maxS_e)
        for each_e_ in maxSC_l:
            bc_conE_count[each_e_] += 1
    return bc_conE_count

def getBCconEvent_comSCOREandUMI(BCumiInfo_dict):
    BC_raw_con_count_comSCOREandUMI = getRawConEventCount_comSCOREandUMI(BCumiInfo_dict)
    BC_conEvent = {}
    moreInfo_dict = defaultdict(lambda: defaultdict(dict))
    BC_info = defaultdict(dict)
    for BC, umiInfo_dict in BCumiInfo_dict.items():
        # get some raw info
        BC_info[BC]["numOfUMI"] = len(umiInfo_dict.keys())
        event_list = [umiInfo_dict[j][0] for j in umiInfo_dict.keys()]
        # get score info
        umi_score, edit_score = scoreUMIandEdit(umiInfo_dict)
        # count edit
        raw_edit_count = getEventCount(umiInfo_dict)
        com_edit_count = deepcopy(raw_edit_count)
        raw_maxUMI_edit_list = getMaxVkey(raw_edit_count)
        BC_info[BC]["UMImax"] = raw_edit_count[raw_maxUMI_edit_list[0]]
        BC_info[BC]["numOfUMImax"] = len(raw_maxUMI_edit_list)
        # chage the count of edit
        max_score_edit = getMaxVkey(edit_score)
        BC_info[BC]["numOfSCOREmax"] = len(max_score_edit)
        ###
        have_inher = False
        for check_edit in set(event_list):
            for each_edit in event_list:
                if checkInfo(check_edit, each_edit) == "anc":
                    have_inher = True
                    com_edit_count[check_edit] += 1
        # get inher info 
        if have_inher:
            BC_info[BC]["haveInher"] = "Y"
        else:
            BC_info[BC]["haveInher"] = "N"
        # get max count edit
        maxCount_com_edit = getMaxVkey(com_edit_count)
        BC_info[BC]["numOfcomSCOREandUMImax"] = len(maxCount_com_edit)
        BC_info[BC]["comSCOREandUMImax"] = com_edit_count[maxCount_com_edit[0]]
        ## comprehensive thinking to get the con event
        # 1. if more than 1 max count edit select the high score
        comMax_SCORE = {}
        for raw_max_com in maxCount_com_edit:
            comMax_SCORE[raw_max_com] = edit_score[raw_max_com]
        comMax_SCOREmax_edit = getMaxVkey(comMax_SCORE)
        # 2. if more than 1 max count and max score select the most frequency in the tree
        conEvent = ""
        comCount = 0
        for maxCom_edit in comMax_SCOREmax_edit:
            if BC_raw_con_count_comSCOREandUMI[maxCom_edit] > comCount:
                comCount = BC_raw_con_count_comSCOREandUMI[maxCom_edit]
                conEvent = maxCom_edit
        BC_conEvent[BC] = conEvent
        #####
        # get moreinfos
        for umi, info in umiInfo_dict.items():
            moreInfo_dict[BC][umi]["events"] = info[0]
            moreInfo_dict[BC][umi]["zwCount"] = info[1]
            moreInfo_dict[BC][umi]["rawUMIcount"] = raw_edit_count[info[0]]
            moreInfo_dict[BC][umi]["comUMIcount"] = com_edit_count[info[0]]
            moreInfo_dict[BC][umi]["SCORE"] = umi_score[umi]
            moreInfo_dict[BC][umi]["checkInfo"] = checkInfo(conEvent, info[0])
    return BC_conEvent, moreInfo_dict, BC_info

def Parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=str, nargs='?', required=True, help="Input the RAW Event file")
    parser.add_argument('-f', '--filter', type=str, nargs='?', required=True, help="The BC which wans to keep in the tree, pass a file contain each BC each line,\nif no filter file, use no")
    parser.add_argument('--outDir', type=str, nargs='?', help="Output results PATH")
    parser.add_argument('-n', '--name', type=str, nargs='?', help="Output files prefix")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    options = Parsers()
    ## input file
    inRawEvent = open(options.infile, 'r')
    ## output files
    outComSCOREandUMI = open(os.path.join(options.outDir, options.name + "_comSCOREandUMIconEvents.txt"), 'w')
    outMoreInfos = open(os.path.join(options.outDir, options.name + "_moreInfos.txt"), 'w')
    outStaInfos = open(os.path.join(options.outDir, options.name + "_BCstaInfos.txt"), 'w')
    ## run
    if options.filter == "no":
        BClist = options.filter
    else:
        BClist = getFilterBClist(options.filter)
    BC_umi_info_dict = getBCumiDict(inRawEvent, BClist)
    conEvent_dict, moreInfo_dict, Sta_dict=getBCconEvent_comSCOREandUMI(BC_umi_info_dict)
    ## write out results
    #1. write con event results
    # write header
    con_header = "{}\t{}\n".format("BC", "\t".join(["target" + str(i) for i in range(1, 14)]))
    outComSCOREandUMI.write(con_header)
    # write events
    for BC, conEvent in conEvent_dict.items():
        outComSCOREandUMI.write("{}\t{}\n".format(BC, conEvent))
    outComSCOREandUMI.close()
    #2. write out more infos
    info_list = ["BC", "UMI", "zwCount", "rawUMIcount", "comUMIcount", "SCORE", "checkInfo"]
    key_list = ["zwCount", "rawUMIcount", "comUMIcount", "SCORE", "checkInfo", "events"]
    moreInfo_header = "{}\t{}\n".format("\t".join(info_list), "\t".join(["target" + str(i) for i in range(1, 14)]))
    outMoreInfos.write(moreInfo_header)
    for BC_, umi_info_dict in moreInfo_dict.items():
        for umi_, info_dict in umi_info_dict.items():
            outMoreInfos.write("{}\t{}\t{}\n".format(BC_, umi_, "\t".join([str(info_dict[i]) for i in key_list])))
    outMoreInfos.close()
    #3. write out Sta infos
    outStaInfos.write("{}\n".format("\t".join(["BC", "numOfUMI", "UMIcountMax", "numOfUMImax", "numOfSCOREmax", "numOfcomSCOREandUMImax", "comSCOREandUMImax", "haveInher"])))
    key_list2 = ["numOfUMI", "UMImax", "numOfUMImax", "numOfSCOREmax", "numOfcomSCOREandUMImax", "comSCOREandUMImax", "haveInher"]
    for bc, sta_info_dict in Sta_dict.items():
        outStaInfos.write("{}\t{}\n".format(bc, "\t".join([str(sta_info_dict[l]) for l in key_list2])))
    outStaInfos.close()




