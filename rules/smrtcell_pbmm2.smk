ruleorder: pbmm2_align_ubam > pbmm2_align_fastq


rule pbmm2_align_ubam:
    input:
        reference = config['ref']['fasta'],
        ref_index = config['ref']['index'],
        query = lambda wildcards: ubam_dict[wildcards.sample][wildcards.movie]
    output:
        bam = f"samples/{{sample}}/aligned/{{movie}}.{ref}.bam",
        bai = f"samples/{{sample}}/aligned/{{movie}}.{ref}.bam.bai"
    log: f"samples/{{sample}}/logs/pbmm2/align/{{movie}}.{ref}.log"
    benchmark: f"samples/{{sample}}/benchmarks/pbmm2/align/{{movie}}.{ref}.tsv"
    params:
        sample = lambda wildcards: wildcards.sample,
        preset = "CCS",
        extra = "--sort --unmapped",
        loglevel = "INFO"
    threads: 28
    conda: "envs/pbmm2.yaml"
    shell:
        """
        (pbmm2 align --num-threads {threads} \
            --preset {params.preset} \
            --sample {params.sample} \
            --log-level {params.loglevel} \
            {params.extra} \
            {input.reference} \
            {input.query} \
            {output.bam}) > {log} 2>&1
        """


rule pbmm2_align_fastq:
    input:
        reference = config['ref']['fasta'],
        ref_index = config['ref']['index'],
        query = lambda wildcards: fastq_dict[wildcards.sample][wildcards.movie]
    output:
        bam = f"samples/{{sample}}/aligned/{{movie}}.{ref}.bam",
        bai = f"samples/{{sample}}/aligned/{{movie}}.{ref}.bam.bai"
    log: f"samples/{{sample}}/logs/pbmm2/align/{{movie}}.{ref}.log"
    benchmark: f"samples/{{sample}}/benchmarks/pbmm2/align/{{movie}}.{ref}.tsv"
    params:
        sample = lambda wildcards: wildcards.sample,
        preset = "CCS",
        extra = "--sort --unmapped",
        loglevel = "INFO"
    threads: 28
    conda: "envs/pbmm2.yaml"
    shell:
        """
        (pbmm2 align --num-threads {threads} \
            --preset {params.preset} \
            --sample {params.sample} \
            --log-level {params.loglevel} \
            {params.extra} \
            {input.reference} \
            {input.query} \
            {output.bam}) > {log} 2>&1
        """
