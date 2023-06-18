ruleorder: samtools_fasta > seqtk_fastq_to_fasta


rule samtools_fasta:
    input: lambda wildcards: ubam_dict[wildcards.sample][wildcards.movie]
    output: temp("samples/{sample}/jellyfish/{movie}.fasta")
    log: "samples/{sample}/logs/samtools/fasta/{movie}.log"
    benchmark: "samples/{sample}/benchmarks/samtools/fasta/{movie}.tsv"
    threads: 4
    conda: "envs/samtools.yaml"
    shell: "(samtools fasta -@ 3 {input} > {output}) > {log} 2>&1"


rule seqtk_fastq_to_fasta:
    input: lambda wildcards: fastq_dict[wildcards.sample][wildcards.movie]
    output: temp("samples/{sample}/jellyfish/{movie}.fasta")
    log: "samples/{sample}/logs/seqtk/seq/{movie}.log"
    benchmark: "samples/{sample}/benchmarks/seqtk/seq/{movie}.tsv"
    conda: "envs/seqtk.yaml"
    shell: "(seqtk seq -A {input} > {output}) > {log} 2>&1"


rule jellyfish_count:
    input: "samples/{sample}/jellyfish/{movie}.fasta"
    output: "samples/{sample}/jellyfish/{movie}.jf"
    log: "samples/{sample}/logs/jellyfish/count/{movie}.log"
    benchmark: "samples/{sample}/benchmarks/jellyfish/count/{movie}.tsv"
    params:
        kmer_length = config['kmer_length'],
        size = 1000000000,
        extra = "--canonical --disk"
    threads: 28
    conda: "envs/jellyfish.yaml"
    shell:
        """
        (jellyfish count {params.extra} \
            --mer-len={params.kmer_length} \
            --size={params.size} \
            --threads={threads} \
            --output={output} \
            {input}) > {log} 2>&1
        """


rule dump_modimers:
    input: "samples/{sample}/jellyfish/{movie}.jf"
    output: "samples/{sample}/jellyfish/{movie}.modimers.tsv.gz"
    log: "samples/{sample}/logs/dump_modimers/{movie}.log"
    benchmark: "samples/{sample}/benchmarks/dump_modimers/{movie}.tsv"
    conda: "envs/jellyfish.yaml"
    threads: 2
    shell:
        """
        (jellyfish dump -c -t {input} \
            | PYTHONHASHSEED=0 python workflow/scripts/modimer.py -N 5003 /dev/stdin \
            | sort | gzip > {output}) > {log} 2>&1
        """
