#######################################
## STARna: Snakefile for RNA mapping ##
#######################################

import os
#conda_prefix = str(os.environ["CONDA_PREFIX"])

import sys
#sys.path.insert(1, conda_prefix+"/lib/python"+str(sys.version_info[0])+"."+str(sys.version_info[1])+"/site-packages")

from typing import List
import pathlib
import re
import numpy
import pandas as pd
import math
from itertools import combinations


# Define general variables
genome_fasta = str(config["genome_fasta"])


### working directory
home_dir = os.path.join(config["output_directory"],"")
shell('mkdir -p {home_dir}')
workdir: home_dir



### Other variables
if (eval(str(config["paired_end"])) == True):
    paired = True
else:
    paired = False



# get the unique samples names and other variables
if not (os.path.exists(config["fastq_directory"])):
    os.system("printf '\033[1;31m\\n!!! *fastq_directory* does not exist !!!\\n\\n\033[0m'")

FILENAMES = next(os.walk(config["fastq_directory"]))[2]
RUNNAMES = [re.sub(rf"{config['fastq_suffix']}$", "", i) for i in FILENAMES]
SAMPLENAMES = numpy.sort(numpy.unique([re.sub(rf"{config['read_suffix'][0]}|{config['read_suffix'][1]}.*$", "", i) for i in RUNNAMES]))


### Generation of global wildcard_constraints
# Function to handle the values for the wilcards
def constraint_to(values: List[str]) -> str:
    """
    From a list, return a regular expression allowing each
    value and not other.
    ex: ["a", "b", "v"] -> (a|b|v)
    """
    if isinstance(values, str):
            raise ValueError("constraint_to(): Expected a list, got str instead")
    return "({})".format("|".join(values))

wildcard_constraints:
    SAMPLE = constraint_to(SAMPLENAMES),
    RUNS = constraint_to(RUNNAMES)


# ruleorder: fastQC_filtered_BAM > normalized_bigWig > raw_bigWig

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================
# Function to run all funtions
rule AAA_initialization:
    input:
        multiqc_bam_report = os.path.join("03_quality_controls/multiQC_bam.html"),
        merged_counts_unstranded = os.path.join("04_gene_counts/unstranded_merged_gene_counts.txt"),
        merged_counts_plus = os.path.join("04_gene_counts/strand.forward_merged_gene_counts.txt"),
        merged_counts_minus = os.path.join("04_gene_counts/strand.reverse_merged_gene_counts.txt")
    shell:
        """
        printf '\033[1;36mPipeline ended!\\n\033[0m'
        """

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================

# Generate STAR index if required [OPTIONAL]  ---------------------------------------------
if not os.path.exists(os.path.join(config["genome_directory"], "transcriptInfo.tab")):
    rule generate_genome_index:
        input:
            genome = ancient(genome_fasta),
            gtf = ancient(config["gtf"])
        output:
            transcriptInfo = os.path.join(config["genome_directory"], "transcriptInfo.tab")
        params:
            genome_dir = config["genome_directory"],
            overhang = config["sjdbOverhang"],
            extra_param = config["extra_param_index"]
        threads:
            workflow.cores
        shell:
            """
            printf '\033[1;36mGenerating the genome index...\\n\033[0m'
            ${{CONDA_PREFIX}}/bin/STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {params.genome_dir} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf} \
-           --sjdbOverhang {params.extra_param}

            printf '\033[1;36mGenome index done.\\n\033[0m'
            """



# cutdapat -------------------------------------------------------------------------------
if paired:
    rule cutadapt_PE:
        input:
            R1 = os.path.join(config["fastq_directory"], "".join(["{SAMPLE}", config['read_suffix'][0], config['fastq_suffix']])),
            R2 = os.path.join(config["fastq_directory"], "".join(["{SAMPLE}", config['read_suffix'][1], config['fastq_suffix']]))
        output:
            R1_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][0], "_trimmed.fastq.gz"])),
            R2_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][1], "_trimmed.fastq.gz"]))
        params:
            sample = "{SAMPLE}",
            opts = str(config["cutadapt_trimm_options"]),
            fw_adapter_sequence = str(config["fw_adapter_sequence"]),
            rv_adapter_sequence = str(config["fw_adapter_sequence"])
        log:
            out = "01_trimmed_fastq/logs/cutadapt.{SAMPLE}.out",
            err = "01_trimmed_fastq/logs/cutadapt.{SAMPLE}.err"
        threads:
            workflow.cores
        priority: 50
        shell:
            """
            printf '\033[1;36m{params.sample}: reads trimming...\\n\033[0m'
            mkdir -p 01_trimmed_fastq/logs/

            ${{CONDA_PREFIX}}/bin/cutadapt \
            -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 25 \
            -a {params.fw_adapter_sequence} -A {params.rv_adapter_sequence} {params.opts} \
            -o {output.R1_trimm} -p {output.R2_trimm} {input.R1} {input.R2} > {log.out} 2> {log.err}
            """
else:
    rule cutadapt_SE:
        input:
            R1 = os.path.join(config["fastq_directory"], "".join(["{SAMPLE}", config['read_suffix'][0], config['fastq_suffix']]))
        output:
            R1_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][0], "_trimmed.fastq.gz"]))
        params:
            sample = "{SAMPLE}",
            opts = str(config["cutadapt_trimm_options"]),
            fw_adapter_sequence = str(config["fw_adapter_sequence"]),
            rv_adapter_sequence = str(config["fw_adapter_sequence"])
        log:
            out = "01_trimmed_fastq/logs/cutadapt.{SAMPLE}.out",
            err = "01_trimmed_fastq/logs/cutadapt.{SAMPLE}.err"
        priority: 50
        threads:
            workflow.cores
        shell:
            """
            printf '\033[1;36m{params.sample}: reads trimming...\\n\033[0m'
            mkdir -p 01_trimmed_fastq/logs

            $CONDA_PREFIX/bin/cutadapt \
            -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 25 \
            -a {params.fw_adapter_sequence} {params.opts} \
            -o {output.R1_trimm} {input.R1} > {log.out} 2> {log.err}
            """



# STAR mapping -----------------------------------------------------------------------------
if (paired == True):
    rule STAR_PE:
        input:
            R1_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][0], "_trimmed.fastq.gz"])),
            R2_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][1], "_trimmed.fastq.gz"])),
            transcriptInfo = ancient(os.path.join(config["genome_directory"], "transcriptInfo.tab"))
        output:
            mapping_results = expand(os.path.join("02_mapping/{SAMPLE}_{ext}"), ext = ["ReadsPerGene.out.tab", "SJ.out.tab"], allow_missing=True),
            bam = os.path.join("02_mapping/{SAMPLE}_Aligned.sortedByCoord.out.bam"),
            bai = os.path.join("02_mapping/{SAMPLE}_Aligned.sortedByCoord.out.bai"),
            temp_folder = temp(directory("02_mapping/{SAMPLE}__STARtmp")),
            flagstat = os.path.join("02_mapping/flagstat/{SAMPLE}_flagstat.txt")
        params:
            genome_dir = config["genome_directory"],
            sample = "{SAMPLE}",
            output_prefix = "02_mapping/{SAMPLE}_",
            readFilesCommand = config["readFilesCommand"],
            star_extra_options = config["star_extra_options"]
        threads:
            max(math.floor(workflow.cores/int(config["n_parallel_samples"])), 1)
        log:
            log_final = os.path.join("02_mapping/{SAMPLE}_Log.final.out"),
            log_out = os.path.join("02_mapping/{SAMPLE}_Log.out"),
            log_progress = os.path.join("02_mapping/{SAMPLE}_Log.progress.out")
        shell:
            """
            printf '\033[1;36m{params.sample}: STAR mapping and gene counts...\\n\033[0m'

            $CONDA_PREFIX/bin/STAR \
            --runThreadN {threads} \
            --genomeDir {params.genome_dir} \
            --readFilesIn {input.R1_trimm} {input.R2_trimm} \
            --readFilesCommand {params.readFilesCommand} \
            --outFileNamePrefix {params.output_prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts {params.star_extra_options}

            printf '\033[1;36m{params.sample}: Indexing bam...\\n\033[0m'
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam} {output.bai}

            printf '\033[1;36m{params.sample}: Computing flagstat...\\n\033[0m'
            $CONDA_PREFIX/bin/samtools flagstat {output.bam} > {output.flagstat}
            """

else:
    rule STAR_SE:
        input:
            R1_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][0], "_trimmed.fastq.gz"])),
            transcriptInfo = ancient(os.path.join(config["genome_directory"], "transcriptInfo.tab"))
        output:
            mapping_results = expand(os.path.join("02_mapping/{SAMPLE}_{ext}"), ext = ["ReadsPerGene.out.tab", "SJ.out.tab"], allow_missing=True),
            bam = os.path.join("02_mapping/{SAMPLE}_Aligned.sortedByCoord.out.bam"),
            bai = os.path.join("02_mapping/{SAMPLE}_Aligned.sortedByCoord.out.bai"),
            temp_folder = temp(directory("02_mapping/{SAMPLE}__STARtmp")),
            flagstat = os.path.join("02_mapping/flagstat/{SAMPLE}_flagstat.txt")
        params:
            genome_dir = config["genome_directory"],
            sample = "{SAMPLE}",
            output_prefix = "02_mapping/{SAMPLE}_",
            readFilesCommand = config["readFilesCommand"],
            star_extra_options = config["star_extra_options"]
        threads:
            max(math.floor(workflow.cores/int(config["n_parallel_samples"])), 1)
        log:
            log_final = os.path.join("02_mapping/{SAMPLE}_Log.final.out"),
            log_out = os.path.join("02_mapping/{SAMPLE}_Log.out"),
            log_progress = os.path.join("02_mapping/{SAMPLE}_Log.progress.out")
        shell:
            """
            printf '\033[1;36m{params.sample}: STAR mapping and gene counts...\\n\033[0m'

            $CONDA_PREFIX/bin/STAR \
            --runThreadN {threads} \
            --genomeDir {params.genome_dir} \
            --readFilesIn {input.R1_trimm} \
            --readFilesCommand {params.readFilesCommand} \
            --outFileNamePrefix {params.output_prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts {params.star_extra_options}

            printf '\033[1;36m{params.sample}: Indexing bam...\\n\033[0m'
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam} {output.bai}

            printf '\033[1;36m{params.sample}: Computing flagstat...\\n\033[0m'
            $CONDA_PREFIX/bin/samtools flagstat {output.bam} > {output.flagstat}
            """




# fastqc on bams -----------------------------------------------------------------------------
rule fastQC_bams:
    input:
        bam = os.path.join("02_mapping/{SAMPLE}_Aligned.sortedByCoord.out.bam")
    output:
        fastqc_html = os.path.join("03_quality_controls/fastqc/{SAMPLE}_Aligned.sortedByCoord.out_fastqc.html"),
        fastqc_zip = os.path.join("03_quality_controls/fastqc/{SAMPLE}_Aligned.sortedByCoord.out_fastqc.zip")
    params:
        sample = "{SAMPLE}"
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    shell:
        """
        printf '\033[1;36m{params.sample}: fastQC on bam...\\n\033[0m'
        $CONDA_PREFIX/bin/fastqc -t {threads} --outdir 03_quality_controls/fastqc {input.bam}
        """


rule multiQC_bam:
    input:
        flagstat = expand(os.path.join("02_mapping/flagstat/{sample}_flagstat.txt"), sample = SAMPLENAMES),
        fastqc = expand(os.path.join("03_quality_controls/fastqc/{sample}_Aligned.sortedByCoord.out_fastqc.zip"), sample = SAMPLENAMES)
    output:
        multiqc_bam_report = os.path.join("03_quality_controls/multiQC_bam.html")
    params:
        out_directory = "03_quality_controls",
        multiqc_bam_report_name = "multiQC_bam.html"
    log:
        out = "03_quality_controls/log/multiQC_bam_filtered.out",
        err = "03_quality_controls/log/multiQC_bam_filtered.err"
    threads: 1
    shell:
        """
        printf '\033[1;36mPerforming multiQC on BAMs/flagstat/STAR metrics...\\n\033[0m'

        mkdir -p {params.out_directory}

        $CONDA_PREFIX/bin/multiqc -f \
        -o {params.out_directory} \
        -n {params.multiqc_bam_report_name} \
        --dirs 02_mapping 02_mapping/flagstat 03_quality_controls/fastqc > {log.err} 2> {log.out}
        """



# merging gene counts -----------------------------------------------------------------------------
rule merge_counts:
    input:
        count_table = expand(os.path.join("02_mapping/{sample}_ReadsPerGene.out.tab"), sample = SAMPLENAMES),
    output:
        merged_counts_unstranded = os.path.join("04_gene_counts/unstranded_merged_gene_counts.txt"),
        merged_counts_plus = os.path.join("04_gene_counts/strand.forward_merged_gene_counts.txt"),
        merged_counts_minus = os.path.join("04_gene_counts/strand.reverse_merged_gene_counts.txt")
    params:
        colnames = ['geneID'] + list(SAMPLENAMES)
    threads:
        workflow.cores
    priority: 100
    run:
        shell("printf '\033[1;36mMerging gene counts...\\n\033[0m'")

        import pandas as pd

        # Dictionary to hold dataframes
        dfs = {}

        # Read all files
        for i, f in enumerate(list(input.count_table), start=1):
            df = pd.read_csv(f, sep="\t", skiprows=4, header=None)  # skip header rows, no headers
            df.columns = ["geneID", "unstranded", "forward", "reverse"]           # assign headers
            dfs[f"df{i}"] = df

        # Merge function
        def merge_on_col(dfs, colname):
            merged = None
            for name, df in dfs.items():
                sub = df[["geneID", colname]].copy()
                sub = sub.rename(columns={colname: f"{colname}_{name}"})  # avoid duplicate names
                if merged is None:
                    merged = sub
                else:
                    merged = pd.merge(merged, sub, on="geneID", how="outer")
            return merged

        # Create the 3 merged tables
        merged_1_2 = merge_on_col(dfs, "unstranded")
        merged_1_2.columns = params.colnames

        merged_1_3 = merge_on_col(dfs, "forward")
        merged_1_3.columns = params.colnames

        merged_1_4 = merge_on_col(dfs, "reverse")
        merged_1_4.columns = params.colnames

        # Save results if needed
        merged_1_2.to_csv(str(output.merged_counts_unstranded), sep="\t", index=False)
        merged_1_3.to_csv(str(output.merged_counts_plus), sep="\t", index=False)
        merged_1_4.to_csv(str(output.merged_counts_minus), sep="\t", index=False)



# ------------------------------------------------------------------------------
#                                 END pipeline
# ------------------------------------------------------------------------------
