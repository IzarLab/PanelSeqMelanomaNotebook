# Snakefile — WES/panel pre-processing → deduplicated BAM (streaming, tidy outputs)

configfile: "config.yaml"

SAMPLES     = config["samples"]
IN_DIR      = config["input_dir"].rstrip("/")
OUT_DIR     = config["output_dir"].rstrip("/")
REF_FASTA   = config["reference"]["fasta"]
REF_FAI     = config["reference"]["fai"]
REF_DICT    = config["reference"]["dict"]
GATK        = config["tools"]["gatk"]
BWA         = config["tools"]["bwa"]
SAMTOOLS    = config["tools"]["samtools"]

# Final targets
rule all:
    input:
        expand(f"{OUT_DIR}/{{sample}}_dedup.bam", sample=SAMPLES),
        expand(f"{OUT_DIR}/{{sample}}_dedup.bai",  sample=SAMPLES)

###############################################################################
# 1) FastqToSam → uBAM
###############################################################################
rule fastq_to_sam:
    input:
        read1 = lambda wc: f"{IN_DIR}/{wc.sample}/{wc.sample}_R1.fastq.gz",
        read2 = lambda wc: f"{IN_DIR}/{wc.sample}/{wc.sample}_R2.fastq.gz"
    output:
        bam = temp(f"{OUT_DIR}/{{sample}}_unaligned_reads.bam")
    params:
        gatk = GATK,
        rg   = "rg0013"
    threads: 2
    log:
        f"{OUT_DIR}/{{sample}}/logs/{{sample}}.fastq_to_sam.log"
    shell:
        r"""
        mkdir -p {OUT_DIR}/{wildcards.sample}/logs {OUT_DIR}/{wildcards.sample}/metrics
        {params.gatk} FastqToSam \
          F1={input.read1} \
          F2={input.read2} \
          O={output.bam} \
          SAMPLE_NAME='{wildcards.sample}' \
          RG='{params.rg}' \
          2> {log}
        """

###############################################################################
# 2) RevertSam
###############################################################################
rule revert_sam:
    input:
        bam = f"{OUT_DIR}/{{sample}}_unaligned_reads.bam"
    output:
        bam = temp(f"{OUT_DIR}/{{sample}}_unaligned_reads_revert.bam")
    params:
        gatk = GATK
    threads: 2
    log:
        f"{OUT_DIR}/{{sample}}/logs/{{sample}}.revert_sam.log"
    shell:
        r"""
        mkdir -p {OUT_DIR}/{wildcards.sample}/logs {OUT_DIR}/{wildcards.sample}/metrics
        {params.gatk} RevertSam \
          -I {input.bam} \
          -O {output.bam} \
          2> {log}
        """

###############################################################################
# 3) MarkIlluminaAdapters
###############################################################################
rule mark_illumina_adapters:
    input:
        bam = f"{OUT_DIR}/{{sample}}_unaligned_reads_revert.bam"
    output:
        bam     = temp(f"{OUT_DIR}/{{sample}}_markilluminaadapters.bam"),
        metrics = f"{OUT_DIR}/{{sample}}/metrics/{{sample}}_markilluminaadapters_metrics.txt"
    params:
        gatk = GATK
    threads: 2
    log:
        f"{OUT_DIR}/{{sample}}/logs/{{sample}}.mark_illumina_adapters.log"
    shell:
        r"""
        mkdir -p {OUT_DIR}/{wildcards.sample}/logs {OUT_DIR}/{wildcards.sample}/metrics
        {params.gatk} MarkIlluminaAdapters \
          -I {input.bam} \
          -O {output.bam} \
          -M {output.metrics} \
          2> {log}
        """

###############################################################################
# 4) Stream: SamToFastq → bwa mem → MergeBamAlignment
###############################################################################
rule align_and_merge:
    input:
        unmapped_bam = f"{OUT_DIR}/{{sample}}_unaligned_reads_revert.bam",
        marked_bam   = f"{OUT_DIR}/{{sample}}_markilluminaadapters.bam",
        ref          = REF_FASTA,
        fai          = REF_FAI,
        dict         = REF_DICT
    output:
        bam = temp(f"{OUT_DIR}/{{sample}}_mergebamalignment.bam"),
        bai = temp(f"{OUT_DIR}/{{sample}}_mergebamalignment.bai")
    params:
        gatk = GATK,
        bwa  = BWA,
        ref  = REF_FASTA
    threads: 7
    log:
        f"{OUT_DIR}/{{sample}}/logs/{{sample}}.align_and_merge.log"
    shell:
        r"""
        set -o pipefail
        mkdir -p {OUT_DIR}/{wildcards.sample}/logs {OUT_DIR}/{wildcards.sample}/metrics
        {params.gatk} SamToFastq \
          -I {input.marked_bam} \
          -FASTQ /dev/stdout \
          -CLIPPING_ATTRIBUTE XT \
          -CLIPPING_ACTION 2 \
          -INTERLEAVE true \
          -NON_PF true \
        | {params.bwa} mem -M -t {threads} -p {params.ref} - \
        | {params.gatk} MergeBamAlignment \
            -R {params.ref} \
            -UNMAPPED_BAM {input.unmapped_bam} \
            -ALIGNED_BAM /dev/stdin \
            -O {output.bam} \
            -CREATE_INDEX true \
            -ADD_MATE_CIGAR true \
            -CLIP_ADAPTERS false \
            -CLIP_OVERLAPPING_READS true \
            -INCLUDE_SECONDARY_ALIGNMENTS true \
            -MAX_INSERTIONS_OR_DELETIONS -1 \
            -PRIMARY_ALIGNMENT_STRATEGY MostDistant \
            -ATTRIBUTES_TO_RETAIN XS \
          2> {log}
        """

###############################################################################
# 5) MarkDuplicates (final deliverables)
###############################################################################
rule mark_duplicates:
    input:
        bam = f"{OUT_DIR}/{{sample}}_mergebamalignment.bam",
        bai = f"{OUT_DIR}/{{sample}}_mergebamalignment.bai"
    output:
        bam     = f"{OUT_DIR}/{{sample}}_dedup.bam",
        bai     = f"{OUT_DIR}/{{sample}}_dedup.bai",
        metrics = f"{OUT_DIR}/{{sample}}/metrics/{{sample}}_dedup.metrics.txt"
    params:
        gatk = GATK,
        java_mem = "16g",
        tmp_java = f"{OUT_DIR}/tmp_java_{{sample}}"
    threads: 4
    log:
        f"{OUT_DIR}/{{sample}}/logs/{{sample}}.markduplicates.log"
    shell:
        r"""
        set -euo pipefail

        # freshen outputs & dirs
        rm -rf "{output.bam}" "{output.bai}"
        mkdir -p "{OUT_DIR}/{wildcards.sample}/logs" "{OUT_DIR}/{wildcards.sample}/metrics" "{params.tmp_java}"

        # sanity check input BAM
        samtools quickcheck -v "{input.bam}" || exit 1

        # run GATK MarkDuplicates (non-Spark), create index, keep seq dups
        "{params.gatk}" \
          --java-options "-Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp_java}" \
          MarkDuplicates \
            -I "{input.bam}" \
            -O "{output.bam}" \
            -M "{output.metrics}" \
            --CREATE_INDEX true \
            --REMOVE_SEQUENCING_DUPLICATES false \
          2> "{log}"
        """
