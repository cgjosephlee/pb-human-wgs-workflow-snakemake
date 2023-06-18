# localrules: bgzip_vcf, tabix_vcf


rule bgzip_vcf:
    input: f"samples/{sample}/{{prefix}}.vcf"
    output: f"samples/{sample}/{{prefix}}.vcf.gz"
    log: f"samples/{sample}/logs/bgzip/{{prefix}}.log"
    benchmark: f"samples/{sample}/benchmarks/bgzip/{{prefix}}.tsv"
    threads: 2
    conda: "envs/htslib.yaml"
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


rule tabix_vcf:
    input: f"samples/{sample}/{{prefix}}.vcf.gz"
    output: f"samples/{sample}/{{prefix}}.vcf.gz.tbi"
    log: f"samples/{sample}/logs/tabix/index/{{prefix}}.log"
    benchmark: f"samples/{sample}/benchmarks/tabix/index/{{prefix}}.tsv"
    params: "-p vcf"
    conda: "envs/htslib.yaml"
    shell: "tabix {params} {input} > {log} 2>&1"


rule samtools_index_bam:
    input: f"samples/{sample}/{{prefix}}.bam"
    output: f"samples/{sample}/{{prefix}}.bam.bai"
    log: f"samples/{sample}/logs/samtools/index/{{prefix}}.log"
    benchmark: f"samples/{sample}/logs/samtools/index/{{prefix}}.tsv"
    threads: 4
    conda: "envs/samtools.yaml"
    shell: "(samtools index -@ 3 {input}) > {log} 2>&1"
