# Snakefile

import pandas as pd

# 설정 파일 로드
configfile: "/home/rkawk/Snakemake/config/config.yaml"

# 샘플 정보 로드
try:
    samples_df = pd.read_csv(config["paths"]["samples_tsv"], sep="\t").set_index("sample_id", drop=False)
    SAMPLES = samples_df["sample_id"].tolist()
except FileNotFoundError:
    print("WARNING: samples.tsv not found or path in config.yaml is incorrect. Assuming SAMPLES = ['control_1']")
    SAMPLES = ["control_1"]


# rule all: 최종적으로 생성하고자 하는 모든 파일을 명시
rule all:
    input:
        # 통합된 카운트 매트릭스
        "results/featurecounts/all_samples.counts.tsv",
        expand("results/star/{sample}/Log.final.out", sample=SAMPLES)

rule star_genome_index:
    input:
        fasta = config["paths"]["genome_fasta"],
        gtf = config["paths"]["genome_gtf"]
    output:
        # STAR 인덱스 생성 완료를 나타내는 임의의 파일
        touch("results/star_index/SA")
    params:
        index_dir = "results/star_index/",
        sjdbOverhang = config["params_star"]["sjdbOverhang"]
    log:
        "logs/star_index/genome_generate.log"
    threads: config["threads"].get("star_index", 8)
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeDir {params.index_dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.sjdbOverhang} > {log} 2>&1"

rule star_align:
    input:
        r1 = lambda wildcards: pd.read_csv(config["paths"]["samples_tsv"], sep="\t").set_index("sample_id").loc[wildcards.sample, "fastq1"],
        r2 = lambda wildcards: pd.read_csv(config["paths"]["samples_tsv"], sep="\t").set_index("sample_id").loc[wildcards.sample, "fastq2"],
        index_dir_sentinel = rules.star_genome_index.output,
        gtf = config["paths"]["genome_gtf"]
    output:
        bam = "results/star/{sample}/Aligned.sortedByCoord.out.bam",
        log_final = "results/star/{sample}/Log.final.out",
        log_progress = "results/star/{sample}/Log.progress.out",
        log_stdout = "results/star/{sample}/Log.out",
        sj_out = "results/star/{sample}/SJ.out.tab"
    params:
        index_dir = "results/star_index/",
        prefix = lambda wildcards: f"results/star/{wildcards.sample}/",
        outSAMtype = "BAM SortedByCoordinate",
        outFilterMultimapNmax = config["params_star"]["outFilterMultimapNmax"],
        twopassMode = config["params_star"]["twopassMode"],
        extra_star_params = config["params_star"].get("extra_star_params", "")
    log:
        "logs/star/{sample}.log"
    threads: config["threads"].get("star_align", 10)
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runThreadN {threads} "
        "--genomeDir {params.index_dir} "
        "--readFilesIn {input.r1} {input.r2} "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMtype {params.outSAMtype} "
        "--outFilterMultimapNmax {params.outFilterMultimapNmax} "
        "--twopassMode {params.twopassMode} "
        "{params.extra_star_params} "
        "> {output.log_stdout} 2> {log}"

rule feature_counts:
    input:
        bam = rules.star_align.output.bam,
        gtf = config["paths"]["genome_gtf"]
    output:
        counts = "results/featurecounts/{sample}.counts.txt",
        summary = "results/featurecounts/{sample}.counts.txt.summary"
    params:
        strand = config["params_featurecounts"]["strand_specificity"],
        feature_type = config["params_featurecounts"]["feature_type"],
        attribute_type = config["params_featurecounts"]["attribute_type"],
        paired_end_opt = "-p" if config["params_featurecounts"]["is_paired_end"] else "",
        extra_fc_params = config["params_featurecounts"].get("extra_fc_params", "-O --minOverlap 1 --primary")
    log:
        "logs/featurecounts/{sample}.log"
    threads: config["threads"].get("featurecounts", 8)
    conda:
        "../envs/subread.yaml"
    shell:
        "featureCounts -T {threads} "
        "-a {input.gtf} "
        "-o {output.counts} "
        "-s {params.strand} "
        "-t {params.feature_type} "
        "-g {params.attribute_type} "
        "{params.paired_end_opt} "
        "{params.extra_fc_params} "
        "{input.bam} > {log} 2>&1"

rule merge_counts:
    input:
        expand("results/featurecounts/{sample}.counts.txt", sample=SAMPLES)
    output:
        merged_counts = "results/featurecounts/all_samples.counts.tsv"
    params:
        gene_id_col_name = config["params_featurecounts"]["gene_id_col_name"],
        count_col_index = int(config["params_featurecounts"]["count_col_index"])
    log:
        "logs/merge_counts/merge.log"
    threads: 1
    conda:
        "../envs/python_pandas.yaml"
    script:
        "../workflow/scripts/merge_featurecounts.py"
