### some functions to parse sequence


### give a reference file handle, return the reference id and sequence
def getReference(ref):
    ref_id = next(ref)
    ref_seq = next(ref)
    reference = ref_id + ref_seq
    return reference

### give a FASTA file handle, return an iterator of the FASTA sequence object (id + sequence)
def formFASTA(fa):
    for each in fa.split('>')[1:]:
        faobj = (each.split('\n')[0] , "".join(each.split('\n')[1:]))
        yield faobj

def convert2string(line):
    if type(line) == str:
        return line
    else:
        return line.decode("utf-8", "ignore")

### give a FASTA file handle, return an iterator of the FASTQ sequence object (id + sequence + ind + quality)
def formFASTQ(fq_file_handle):
    while True:
        try:
            header = convert2string(next(fq_file_handle)).strip()
            seq = convert2string(next(fq_file_handle)).strip()
            ind = next(fq_file_handle).strip()
            qul = next(fq_file_handle).strip()
            yield header, seq, ind, qul
        except StopIteration:
            break

### give a FASTQ file handle, return an iterator just the (id + sequence)
def fqTOfa(fq_file_handle):
    while True:
        try:
            header = convert2string(next(fq_file_handle))
            seq = convert2string(next(fq_file_handle))
            next(fq_file_handle)
            next(fq_file_handle)
            yield header, seq
        except StopIteration:
            break

### get the comlement sequence
def seqco(sequence):
    sequence = sequence.upper()
    comlement_Dict = {"A":"T", "T":"A", "C":"G", "G":"C"}
    com_seq = ""
    for i in sequence:
        com_seq += comlement_Dict[i]
    return com_seq.upper()

### get the reversed sequence
def seqre(sequence):
    return sequence.upper()[::-1]
