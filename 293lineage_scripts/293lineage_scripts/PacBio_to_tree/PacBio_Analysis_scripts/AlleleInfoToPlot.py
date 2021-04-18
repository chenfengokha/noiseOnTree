import sys
import os
import argparse

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

def readAlleleInfoToEditDict(alleleInfo_path):
    EditDict = {}
    with open(alleleInfo_path, 'r') as alleleInfo:
        next(alleleInfo)
        for line in alleleInfo:
            spl = line.strip().split("\t")
            name = spl[0].strip()
            if name not in EditDict.keys():
                EditLine = "\t".join(spl[3:])
                EditDict[name] = EditLine
    return EditDict

def EditDictToPlotFile(EditDict, outPlot_path):
    with open(outPlot_path, 'w') as outPlot:
        outPlot.write('\t'.join(["y","start", "end", "editLen","event"]) + '\n')
        for name, editline in EditDict.items():
            edit_pos = editAnotation(editline)
            for plot in edit_pos:
                outPlot.write("{}\t{}\t{}\t{}\t{}\n".format(name, *plot))

def Parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inf', type=str, nargs='?', required=True, help="Input AlleleInfo file")
    parser.add_argument('-o', '--out', type=str, nargs='?', required=True, help="Output plots file")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    options = Parsers()
    # run
    nodeEditDict = readAlleleInfoToEditDict(options.inf)
    EditDictToPlotFile(nodeEditDict, options.out)