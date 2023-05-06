import os
import sys
import glob
import time
## load config file
configfile: "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/scripts/PacBio_analysis_pipline/PacBio_Analysis_config_CF2_mt10_adjust_with10x_rawBC.yaml"

# input raw BAM file
raw_file_path = config["raw_subreads_path"]

# filter BC base on 10x cellrangers
filterBC_dir = config["filterBC_dir"]

# reference 
bowtie_index = config["bowtie_index"]
ref_seq = config["ref_seq"]
tarPos = config["tarPos"]
iqtree_transM = config["iqtree_transMatrix"]

# some hand write scripts
split_strand = config["split_strand"]
HigQulCCS = config["HigQulCCS"]
FilterCCS = config["FilterCCS"]
MergeR1R2CCS = config["MergeR1R2CCS"]
BCumi_con = config["BCumi_con"]
Align_toRef = config["Align_toRef"]
Call_Events = config["Call_Events"]
BC_conEvents = config["BC_conEvents"]
Events_toIqtreeInput = config["Events_toIqtreeInput"]
Run_Iqtree = config["Run_Iqtree"]

AdjustBC = config["AdjustBC"]

## add stastics
stainfo = config["stainfo"]
plotSta = config["plotSta"]
plotTOP100Edits = config["plotTOP100Edits"]
plotEditLine = config["plotEditLine"]
plotAllEdits = config["plotAllEdits"]
plotTreeANDnodeCellNum = config["plotTreeANDnodeCellNum"]


## some parameters
FP = config["P2"]
RP = config["P4"]
inFP = config["inFP"]
inRP = config["inRP"]
primerMis = config["primerMismatch"]
pAwinS = config["polyAwindowSize"]
pAmis = config["polyAmismatch"]
ref_len = config["ref_len"]

# input raw data path list
def get_sample_name(sample_path, suffix=".subreads.bam"):
    return [file.split(".subreads.bam")[0] for file in os.listdir(sample_path) if file.endswith(suffix)]

# sample_name use to define outputs 
#sample_name = get_sample_name(raw_file_path)
sample_name = ["CF2"]
# output directory
## addition out
outDir_raw = config["outDir"]
sample_out_paths = [os.path.join(outDir_raw, j) for j in sample_name]
for each in sample_out_paths:
    if not os.path.exists(each):
        os.makedirs(each)

outDir = config["outDir"]

rule all:
    input: 
        expand(outDir + "/{sample}/mt10/{sample}.adjust_staPlots.pdf", sample = sample_name),
        # expand(outDir + "/{sample}/mt10/{sample}.union.IQtreeInput.txt", sample = sample_name),
        # expand(outDir + "/{sample}/mt10/{sample}.union.Plots.txt", sample = sample_name),
        # expand(outDir + "/{sample}/mt10/{sample}.union.AllelesInfo.txt", sample = sample_name),
        # expand(outDir + "/{sample}/mt10/{sample}.union.EditFre.txt", sample = sample_name)

# rule RAW_BAMtoSAM:
#     input:
#         raw_file_path + "/{sample}.subreads.bam"
#     output:
#         outDir + "/{sample}.subreads.sam"
#     threads: 40
#     shell:
#         "samtools view -h -S -O sam -@ {threads} {input} -o {output}"

# rule split_strand:
#     input:
#         rules.RAW_BAMtoSAM.output
#     output:
#         outDir + "/{sample}/mt10/{sample}.subreads.R1.sam",
#         outDir + "/{sample}/mt10/{sample}.subreads.R2.sam",
#         outDir + "/{sample}/mt10/{sample}.subreads.Other.sam"
#     threads: 20
#     message: "Split the raw SAM files to R1, R2 and Other, base on the read contain polyA(R1), polyT(R2), both and neight is Other"
#     shell:
#         "python {split_strand} -r {input} -1 {output[0]} -2 {output[1]} --other {output[2]} --numcpu {threads}"

# rule SAMtoBAM_beforCCS:
#     input:
#         raw_file_path + "/{sample}.subreads.R1.sam",
#         raw_file_path + "/{sample}.subreads.R2.sam",
#         raw_file_path + "/{sample}.subreads.Other.sam"
#     output:
#         temp(outDir_raw + "/{sample}/mt10/{sample}.subreads.R1.bam"),
#         temp(outDir_raw + "/{sample}/mt10/{sample}.subreads.R2.bam"),
#         #temp(outDir_raw + "/{sample}/mt10/{sample}.subreads.other.bam")
#     threads: 20
#     shell:
#         """
#         samtools view -@ {threads} -b -o {output[0]} {input[0]}
#         samtools view -@ {threads} -b -o {output[1]} {input[1]}
#         """
#         #samtools view -@ {threads} -b -o {output[2]} {input[2]}

# ================================================================================================================

# rule runCCS:
#     input:
#         outDir_raw + "/{sample}/mt10/{sample}.subreads.R1.bam",
#         outDir_raw + "/{sample}/mt10/{sample}.subreads.R2.bam"
#     output:
#         outDir + "/{sample}/mt10/{sample}.R1.ccs.bam",
#         outDir + "/{sample}/mt10/{sample}.R1.ccs.report.txt",
#         outDir + "/{sample}/mt10/{sample}.R2.ccs.bam",
#         outDir + "/{sample}/mt10/{sample}.R2.ccs.report.txt",
#         # outDir + "/{sample}/mt10/{sample}.other.ccs.bam",
#         # outDir + "/{sample}/mt10/{sample}.other.ccs.report.txt",
#     threads: 20
#     shell:
#         """
#         ccs --reportFile={output[1]} --minPasses=10 -j 20 {input[0]} {output[0]}
#         ccs --reportFile={output[3]} --minPasses=10 -j 20 {input[1]} {output[2]}
#         """
#         #         ccs --reportFile={output[5]} --minPasses=3 -j 15 {input[2]} {output[4]}
# rule ccsBAMtoFASTQ:
#     input:
#         rules.runCCS.output[0],
#         rules.runCCS.output[2]
#         #rules.runCCS.output[4]
#     output:
#         outDir + "/{sample}/mt10/{sample}.R1.ccs.raw.fastq",
#         outDir + "/{sample}/mt10/{sample}.R2.ccs.raw.fastq",
#         #outDir + "/{sample}/mt10/{sample}.other.ccs.raw.fastq"
#     shell:
#         """
#         bedtools bamtofastq -i {input[0]} -fq {output[0]}
#         bedtools bamtofastq -i {input[1]} -fq {output[1]}
#         """
#         # bedtools bamtofastq -i {input[2]} -fq {output[2]}

# rule getNoGenomeCCS_FASTQ:
#     input:
#         rules.ccsBAMtoFASTQ.output
#     output:
#         temp(outDir + "/{sample}/mt10/{sample}.R1.ccs.noG.fastq"),
#         temp(outDir + "/{sample}/mt10/{sample}.R1.ccs.unmap.sam"),
#         temp(outDir + "/{sample}/mt10/{sample}.R2.ccs.noG.fastq"),
#         temp(outDir + "/{sample}/mt10/{sample}.R2.ccs.unmap.sam")
#         # temp(outDir + "/{sample}/mt10/{sample}.other.ccs.noG.fastq"),
#         # temp(outDir + "/{sample}/mt10/{sample}.other.ccs.unmap.sam")
#     threads: 20
#     message: "Align to the human genome, exclude the nonspecific amplification"
#     shell:
#         """
#         bowtie2 --threads {threads} -x {bowtie_index} --un {output[0]} -q {input[0]} -S {output[1]}
#         bowtie2 --threads {threads} -x {bowtie_index} --un {output[2]} -q {input[1]} -S {output[3]}
#         """
# #        bowtie2 --threads {threads} -x {bowtie_index} --un {output[4]} -q {input[2]} -S {output[5]}

# rule getHightestQualCCSforEachZwID:
#     input:
#         rules.getNoGenomeCCS_FASTQ.output[0],
#         rules.getNoGenomeCCS_FASTQ.output[2]
#         #rules.getNoGenomeCCS_FASTQ.output[4]
#     output:
#         outDir + "/{sample}/mt10/{sample}.R1.ccs.clean.fastq",
#         outDir + "/{sample}/mt10/{sample}.R2.ccs.clean.fastq"
#         # outDir + "/{sample}/mt10/{sample}.other.ccs.clean.fastq"
#     shell:
#         """
#         python {HigQulCCS} -i {input[0]} -o {output[0]}
#         python {HigQulCCS} -i {input[1]} -o {output[1]}
#         """

# rule filterCCSandExtractBCUMI:
#     input:
#         rules.getHightestQualCCSforEachZwID.output
#     output:
#         outDir + "/{sample}/mt10/{sample}.R1.filter.fa",
#         outDir + "/{sample}/mt10/{sample}.R1.staDF.txt",
#         outDir + "/{sample}/mt10/{sample}.R1.noPass.fa",
#         outDir + "/{sample}/mt10/{sample}.R2.filter.fa",
#         outDir + "/{sample}/mt10/{sample}.R2.staDF.txt",
#         outDir + "/{sample}/mt10/{sample}.R2.noPass.fa"
#     threads: 20
#     shell:
#         """
#         python {FilterCCS} -i {input[0]} -o {output[0]} -s {output[1]} --nopass {output[2]} \
#             --format FASTQ --FP {FP} --RP {RP} --inFP {inFP} --inRP {inRP} --pAw {pAwinS} --pMs {pAmis} --mismatch {primerMis} --numcpu {threads}
#         python {FilterCCS} -i {input[1]} -o {output[3]} -s {output[4]} --nopass {output[5]} \
#             --format FASTQ --FP {FP} --RP {RP} --inFP {inFP} --inRP {inRP} --pAw {pAwinS} --pMs {pAmis} --mismatch {primerMis} --numcpu {threads}
#         """

# rule MergeR1R2:
#     input:
#         outDir + "/{sample}/mt10/{sample}.R1.filter.fa",
#         outDir + "/{sample}/mt10/{sample}.R2.filter.fa"
#         #rules.filterCCSandExtractBCUMI.output[6]
#     output:
#         outDir + "/{sample}/mt10/{sample}.union.filter.fa"
#     shell:
#         "python {MergeR1R2CCS} -1 {input[0]} -2 {input[1]} -o {output}"

rule AdjustBCandUMI:
    input:
        "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results" + "/{sample}/{sample}.union.filter.fa"
    output:
        outDir + "/{sample}/{sample}.union.adjust.fa"
    params:
        filterBC_file = filterBC_dir + "/{sample}_mt10_BC.txt"
    shell:
        "python {AdjustBC} -i {input[0]} -o {output[0]} -f {params.filterBC_file}"

rule getUMIconSeqEachBC:
    input:
        rules.AdjustBCandUMI.output
    output:
        outDir + "/{sample}/mt10/{sample}.union.adjust.UMIcon.fa"
    threads: 20
    shell:
        "python {BCumi_con} -i {input} -o {output} --numcpu {threads}"
    
rule AlignToRef:
    input:
        rules.getUMIconSeqEachBC.output
    output:
        outDir + "/{sample}/mt10/{sample}.union.adjust.UMIcon.align"
    threads: 20
    shell:
        "python {Align_toRef} -r {ref_seq} -s {input} -o {output} -c {threads}"

rule CallEditEventsFromAlignment:
    input:
        rules.AlignToRef.output
    output:
        outDir + "/{sample}/mt10/{sample}.union.adjust.Events.txt"
    threads: 20
    params: 
        output_dir=outDir + "/{sample}/mt10",
        name="{sample}.union.adjust"
    shell:
        "python {Call_Events} -a {input} -t {tarPos} -n {params.name} --dir {params.output_dir} --cpu {threads}"


rule getBCconsensusEditEvents:
    input:
        rules.CallEditEventsFromAlignment.output
    output:
        outDir + "/{sample}/mt10/{sample}.union.adjust_comSCOREandUMIconEvents.txt",
        outDir + "/{sample}/mt10/{sample}.union.adjust_moreInfos.txt",
        outDir + "/{sample}/mt10/{sample}.union.adjust_BCstaInfos.txt"
    params:
        output_dir = outDir + "/{sample}/mt10",
        outprefix = "{sample}.union.adjust",
        filterBC_file = filterBC_dir + "/{sample}_mt10_BC.txt"
    shell:
        "python {BC_conEvents} -i {input} -n {params.outprefix} --outDir {params.output_dir} -f {params.filterBC_file}"



rule eventsToIqtreeInput:
    input:
        rules.getBCconsensusEditEvents.output
    output:
        outDir + "/{sample}/mt10/{sample}.union.adjust.IQtreeInput.txt",
        outDir + "/{sample}/mt10/{sample}.union.adjust.Plots.txt",
        outDir + "/{sample}/mt10/{sample}.union.adjust.AllelesInfo.txt",
        outDir + "/{sample}/mt10/{sample}.union.adjust.EditFre.txt",
        outDir + "/{sample}/mt10/{sample}.union.adjust.nwk"
    params:
        output_dir=outDir + "/{sample}/mt10",
        outprefix = "{sample}.union.adjust",
        filterBC_file = filterBC_dir + "/{sample}_mt10_BC.txt"
    shell:
        "python {Events_toIqtreeInput} -e {input[0]} -m {iqtree_transM} -d {params.output_dir} -n {params.outprefix} -l {ref_len} -f {params.filterBC_file}"

rule allele_staInfo:
    input:
        outDir + "/{sample}/mt10/{sample}.union.adjust.AllelesInfo.txt"
    output:
        outDir + "/{sample}/mt10/{sample}.adjust.Statistics.txt",
        outDir + "/{sample}/mt10/{sample}.adjust.EventInfo.txt"
    params: 
        out_name = "{sample}.adjust",
        output_dir = outDir + "/{sample}/mt10"
    shell:
        "python {stainfo} -i {input[0]} -n {params.out_name} -o {params.output_dir}"

rule plot_StaInfo:
    input:
        outDir + "/{sample}/mt10/{sample}.adjust.Statistics.txt",
        outDir + "/{sample}/mt10/{sample}.union.adjust.AllelesInfo.txt",
        outDir + "/{sample}/mt10/{sample}.union.adjust.Plots.txt",
        outDir + "/{sample}/mt10/{sample}.union.adjust.EditFre.txt",
        outDir + "/{sample}/mt10/{sample}.union.adjust.nwk"

    output:
        outDir + "/{sample}/mt10/{sample}.adjust_staPlots.pdf",
        outDir + "/{sample}/mt10/{sample}.adjust_TOP100_edits.pdf",
        outDir + "/{sample}/mt10/{sample}.adjust_edit_freq_line.pdf",
        outDir + "/{sample}/mt10/{sample}.adjust_all_edits.pdf",
        outDir + "/{sample}/mt10/{sample}.adjust_tree_nodeCellNum.pdf",

    params:
        out_name = "{sample}.adjust",
        output_dir = outDir + "/{sample}/mt10"
    shell:
        """
        Rscript {plotSta} -i {input[0]} -o {output[0]}
        Rscript {plotTOP100Edits} --ai {input[1]} --pl {input[2]} -n {params.out_name} --pdf {output[1]}
        Rscript {plotEditLine} --EditFreq {input[3]} --pdf {output[2]}
        Rscript {plotAllEdits} --pl {input[2]} --pdf {output[3]}
        Rscript {plotTreeANDnodeCellNum} --ai {input[1]} --tree {input[4]} --pdf {output[4]}
        """

## print some messages when success
onsuccess:
    print("Workflow finished, no error.")

