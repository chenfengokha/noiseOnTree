### defined a class pase alignment file

class Alignment:
    def __init__(self, refpos, refbase, readbase, indicator):
        self.refpos, self.refbase, self.readbase, self.indicator = (refpos, refbase, readbase, indicator)
    def combine(self,nextAlign):
        if self.indicator == nextAlign.indicator:
            return [Alignment(self.refpos, self.refbase + nextAlign.refbase, self.readbase + nextAlign.readbase, self.indicator)]
        else:
            return [self, nextAlign]   
    def toEditString(self):
        editEvent = "{}{}+{}".format(len(self.readbase), self.indicator, self.refpos)
        if self.indicator == "I" or self.indicator == "S":
            editEvent = editEvent + "+" + self.readbase
        return editEvent
    def matchPct(self):
        matchCount = 0
        for t in zip(self.refbase, self.readbase):
            if t[0] == t[1]:
                matchCount += 1
        percent = matchCount / len(self.readbase)
        return percent

class Alignment2:
    def __init__(self, refpos, refbase, readbase, indicator):
        self.start, self.refbase, self.readbase, self.indicator = (refpos, refbase, readbase, indicator)
        if self.indicator == "I":
            self.end = self.start
        else:
            self.end = self.start + len(readbase) - 1
        self.editLen = len(self.readbase)
        if self.indicator == "I" or self.indicator == "S":
            self.editseq = self.readbase
        else:
            self.editseq = ""
    def combine(self,nextAlign):
        if self.indicator == nextAlign.indicator:
            return [Alignment2(self.start, self.refbase + nextAlign.refbase, self.readbase + nextAlign.readbase, self.indicator)]
        else:
            return [self, nextAlign]   
    def toEditString(self):
        editEvent = "({},{},{},{},{})".format(self.start, self.end, self.indicator, self.editLen, self.editseq)
        return editEvent
    def matchPct(self):
        matchCount = 0
        for t in zip(self.refbase, self.readbase):
            if t[0] == t[1]:
                matchCount += 1
        percent = matchCount / len(self.readbase)
        return percent
    @classmethod
    def buildFromString(cls, editString):
        spl = list(editString.split(","))
        start = int(spl[0][1:])
        end = int(spl[1])
        indicator = spl[2]
        editLen = int(spl[3][:-1])
        if indicator == "D":
            editSeq = "N" * (end - start + 1)
        else:
            editSeq = spl[4]
        return Alignment2(start, "", editSeq, indicator)