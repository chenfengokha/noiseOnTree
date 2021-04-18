import os
import sys
import glob
import time
## load config file
configfile: "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/scripts/PacBio_analysis_pipline/PacBio_Analysis_sample_config.yaml"

# input raw BAM file
raw_file_path = config["raw_subreads_path"]

# filter BC base on 10x cellrangers
filterBC = config["filterBC"]
exp_mx_dir = config["exp_mx_dir"]

# reference 
bowtie_index = config["bowtie_index"]
ref_seq = config["ref_seq"]
tarPos = config["tarPos"]
iqtree_transM = config["iqtree_transMatrix"]

# some scripts
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
# R scripts
getPair_dist = config["getPair_dist"]
getExpMatrix = config["getExpMatrix"]
getTreeVSexpDist = config["getTreeVSexpDist"]
analysis_share = config["analysis_share"]



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
sample_name = get_sample_name(raw_file_path)

# output directory
outDir = config["outDir"]
sample_out_paths = [os.path.join(outDir, j) for j in sample_name]
for each in sample_out_paths:
    if not os.path.exists(each):
        os.makedirs(each)

rule all:
    input: 
        expand(outDir + "/{sample}/{sample}.union.nwk", sample = sample_name),
        # expand(outDir + "/{sample}/{sample}.union.IQtreeInput.txt", sample = sample_name),
        # expand(outDir + "/{sample}/{sample}.union.Plots.txt", sample = sample_name),
        # expand(outDir + "/{sample}/{sample}.union.AllelesInfo.txt", sample = sample_name),
        # expand(outDir + "/{sample}/{sample}.union.EditFre.txt", sample = sample_name)

rule RAW_BAMtoSAM:
    input:
        raw_file_path + "/{sample}.subreads.bam"
    output:
        outDir + "/{sample}.subreads.sam"
    threads: 20
    shell:
        "samtools view -h -S -O sam -@ {threads} {input} -o {output}"

rule split_strand:
    input:
        rules.RAW_BAMtoSAM.output
    output:
        outDir + "/{sample}/{sample}.subreads.R1.sam",
        outDir + "/{sample}/{sample}.subreads.R2.sam",
        outDir + "/{sample}/{sample}.subreads.Other.sam"
    threads: 20
    message: "Split the raw SAM files to R1, R2 and Other, base on the read contain polyA(R1), polyT(R2), both and neight is Other"
    shell:
        "python {split_strand} -r {input} -1 {output[0]} -2 {output[1]} --other {output[2]} --numcpu {threads}"

rule SAMtoBAM_beforCCS:
    input:
        rules.split_strand.output[0],
        rules.split_strand.output[1]
    output:
        temp(outDir + "/{sample}/{sample}.subreads.R1.bam"),
        temp(outDir + "/{sample}/{sample}.subreads.R2.bam")
    threads: 20
    shell:
        """
        samtools view -@ {threads} -b -o {output[0]} {input[0]}
        samtools view -@ {threads} -b -o {output[1]} {input[1]}
        """

rule runCCS:
    input:
        rules.SAMtoBAM_beforCCS.output
    output:
        outDir + "/{sample}/{sample}.R1.ccs.bam",
        outDir + "/{sample}/{sample}.R1.ccs.report.txt",
        outDir + "/{sample}/{sample}.R2.ccs.bam",
        outDir + "/{sample}/{sample}.R2.ccs.report.txt"
    #threads: 20
    shell:
        """
        ccs --reportFile={output[1]} --minPasses=3 -j 15 {input[0]} {output[0]}
        ccs --reportFile={output[3]} --minPasses=3 -j 15 {input[1]} {output[2]}
        """
rule ccsBAMtoFASTQ:
    input:
        rules.runCCS.output[0],
        rules.runCCS.output[2]
    output:
        outDir + "/{sample}/{sample}.R1.ccs.raw.fastq",
        outDir + "/{sample}/{sample}.R2.ccs.raw.fastq"
    shell:
        """
        bedtools bamtofastq -i {input[0]} -fq {output[0]}
        bedtools bamtofastq -i {input[1]} -fq {output[1]}
        """

rule getNoGenomeCCS_FASTQ:
    input:
        rules.ccsBAMtoFASTQ.output
    output:
        temp(outDir + "/{sample}/{sample}.R1.ccs.noG.fastq"),
        temp(outDir + "/{sample}/{sample}.R1.ccs.unmap.sam"),
        temp(outDir + "/{sample}/{sample}.R2.ccs.noG.fastq"),
        temp(outDir + "/{sample}/{sample}.R2.ccs.unmap.sam")
    threads: 20
    message: "Align to the human genome, exclude the nonspecific amplification"
    shell:
        """
        bowtie2 --threads {threads} -x {bowtie_index} --un {output[0]} -q {input[0]} -S {output[1]}
        bowtie2 --threads {threads} -x {bowtie_index} --un {output[2]} -q {input[1]} -S {output[3]}
        """

rule getHightestQualCCSforEachZwID:
    input:
        rules.getNoGenomeCCS_FASTQ.output[0],
        rules.getNoGenomeCCS_FASTQ.output[2]
    output:
        outDir + "/{sample}/{sample}.R1.ccs.clean.fastq",
        outDir + "/{sample}/{sample}.R2.ccs.clean.fastq"
    shell:
        """
        python {HigQulCCS} -i {input[0]} -o {output[0]}
        python {HigQulCCS} -i {input[1]} -o {output[1]}
        """

rule filterCCSandExtractBCUMI:
    input:
        rules.getHightestQualCCSforEachZwID.output
    output:
        outDir + "/{sample}/{sample}.R1.filter.fa",
        outDir + "/{sample}/{sample}.R1.staDF.txt",
        outDir + "/{sample}/{sample}.R1.noPass.fa",
        outDir + "/{sample}/{sample}.R2.filter.fa",
        outDir + "/{sample}/{sample}.R2.staDF.txt",
        outDir + "/{sample}/{sample}.R2.noPass.fa"
    threads: 10
    shell:
        """
        python {FilterCCS} -i {input[0]} -o {output[0]} -s {output[1]} --nopass {output[2]} \
            --format FASTQ --FP {FP} --RP {RP} --inFP {inFP} --inRP {inRP} --pAw {pAwinS} --pMs {pAmis} --mismatch {primerMis} --numcpu {threads}
        python {FilterCCS} -i {input[1]} -o {output[3]} -s {output[4]} --nopass {output[5]} \
            --format FASTQ --FP {FP} --RP {RP} --inFP {inFP} --inRP {inRP} --pAw {pAwinS} --pMs {pAmis} --mismatch {primerMis} --numcpu {threads}
        """

rule MergeR1R2:
    input:
        rules.filterCCSandExtractBCUMI.output[0],
        rules.filterCCSandExtractBCUMI.output[3]
    output:
        outDir + "/{sample}/{sample}.union.filter.fa"
    shell:
        "python {MergeR1R2CCS} -1 {input[0]} -2 {input[1]} -o {output}"

rule AdjustBCandUMI:
    input:
        rules.MergeR1R2.output
    output:
        outDir + "/{sample}/{sample}.adjust.fa"
    params:
        filterBC_file=filterBC_dir + "/{sample}.filterBCs.txt"
    shell:
        "python {AdjustBC} -i {input[0]} -o {output[0]} -f {params.filterBC_file}"


rule getUMIconSeqEachBC:
    input:
        rules.AdjustBCandUMI.output
    output:
        outDir + "/{sample}/{sample}.union.UMIcon.fa"
    threads: 10
    shell:
        "python {BCumi_con} -i {input} -o {output} --numcpu {threads}"
    
rule AlignToRef:
    input:
        rules.getUMIconSeqEachBC.output
    output:
        outDir + "/{sample}/{sample}.union.UMIcon.align"
    threads: 10
    shell:
        "python {Align_toRef} -r {ref_seq} -s {input} -o {output} -c {threads}"

rule CallEditEventsFromAlignment:
    input:
        rules.AlignToRef.output
    output:
        outDir + "/{sample}/{sample}.union.Events.txt"
    threads: 10
    params: 
        output_dir=outDir + "/{sample}",
        name="{sample}.union"
    shell:
        "python {Call_Events} -a {input} -t {tarPos} -n {params.name} --dir {params.output_dir} --cpu {threads}"

rule getBCconsensusEditEvents:
    input:
        rules.CallEditEventsFromAlignment.output
    output:
        outDir + "/{sample}/{sample}_comSCOREandUMIconEvents.txt",
        outDir + "/{sample}/{sample}_moreInfos.txt",
        outDir + "/{sample}/{sample}_BCstaInfos.txt"
    params:
        output_dir = outDir + "/{sample}",
        outprefix = "{sample}",
        filterBC_file=filterBC_dir + "/{sample}.filterBCs.txt"
    shell:
        "python {BC_conEvents} -i {input} -n {params.outprefix} --outDir {params.output_dir} -f {params.filterBC_file}"

rule eventsToIqtreeInput:
    input:
        rules.getBCconsensusEditEvents.output
    output:
        outDir + "/{sample}/{sample}.nwk",
        outDir + "/{sample}/{sample}.AllelesInfo.txt"
    params:
        output_dir = outDir + "/{sample}",
        outprefix = "{sample}",
        filterBC_file = filterBC_dir + "/{sample}.filterBCs.txt"
    threads:
        40
    shell:
        """
        python {Events_toIqtreeInput} -e {input[0]} -m {iqtree_transM} -d {params.output_dir} -n {params.outprefix} -l {ref_len} -f {params.filterBC_file}
        """

rule getTreeNodePairDist:
    input:
        rules.eventsToIqtreeInput.output
    output:
        outDir + "/{sample}/{sample}.oneCellNode.pair_Dist.txt"
    threads:
        40
    shell:
        """
        Rscript {getPair_dist} -i {input[0]} -o {output[0]}
        """

rule getExpMx:
    input:
        outDir + "/{sample}/{sample}.AllelesInfo.txt",
    output:
        outDir + "/{sample}/{sample}_oneCellexp.Rds",
        outDir + "/{sample}/{sample}_alltreeExp.Rds"
    params:
        filterBC_file = filterBC_dir + "/{sample}.filterBCs.txt",
        exp_mx = exp_mx_dir + "/{sample}_allCell_exp.Rds"
    shell:
        """
        Rscript {getExpMatrix} -i {params.exp_mx} -f {params.filterBC_file} -a {input[0]} --outOneCell {output[0]} --outAllTree {output[1]}
        """

rule getTreeVSexpDist:
    input:
        outDir + "/{sample}/{sample}_oneCellexp.Rds",
        outDir + "/{sample}/{sample}.oneCellNode.pair_Dist.txt"
    output:
        outDir + "/{sample}/{sample}_treeVSexpDist_result.Rds"
    shell:
        """
        Rscript {getTreeVSexpDist} --inTreeDist {input[1]} --inExpMx {input[0]} -o {output[0]}
        """

rule analysis_ANC_share:
    input:
        outDir + "/{sample}/{sample}.AllelesInfo.txt",
        outDir + "/{sample}/{sample}_moreInfos.txt",
        outDir + "/{sample}/{sample}.oneCellNode.pair_Dist.txt",
        outDir + "/{sample}/{sample}_treeVSexpDist_result.Rds"
    output:
        outDir + "/{sample}/{sample}_ANC_share_analysis_results.Rds"
    threads:
        60
    shell:
        "Rscript {analysis_share} --allele {input[0]} --moreInfo {input[1]} --pairDist {input[2]} -o {output}"


## print some messages when success
onsuccess:
    print("Workflow finished, no error.")

