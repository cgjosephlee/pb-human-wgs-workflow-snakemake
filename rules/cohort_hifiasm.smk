ruleorder: samtools_fasta > seqtk_fastq_to_fasta
# localrules: asm_stats, gfa2fa


rule samtools_fasta:
    input: lambda wildcards: ubam_dict[wildcards.sample][wildcards.movie]
    output: temp(f"cohorts/{cohort}/fasta/{{sample}}/{{movie}}.fasta")
    log: f"cohorts/{cohort}/logs/samtools/fasta/{{sample}}/{{movie}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/samtools/fasta/{{sample}}/{{movie}}.tsv"
    threads: 4
    conda: "envs/samtools.yaml"
    shell: "(samtools fasta -@ 3 {input} > {output}) > {log} 2>&1"


rule seqtk_fastq_to_fasta:
    input: lambda wildcards: fastq_dict[wildcards.sample][wildcards.movie]
    output: temp(f"cohorts/{cohort}/fasta/{{sample}}/{{movie}}.fasta")
    log: f"cohorts/{cohort}/logs/seqtk/seq/{{sample}}/{{movie}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/seqtk/seq/{{sample}}/{{movie}}.tsv"
    conda: "envs/seqtk.yaml"
    shell: "(seqtk seq -A {input} > {output}) > {log} 2>&1"


rule yak_count:
    input: lambda wildcards: expand(f"cohorts/{cohort}/fasta/{wildcards.sample}/{{movie}}.fasta", movie=ubam_fastq_dict[wildcards.sample])
    output: temp(f"cohorts/{cohort}/yak/{{sample}}.yak")
    log: f"cohorts/{cohort}/logs/yak/{{sample}}.yak.count.log"
    benchmark: f"cohorts/{cohort}/benchmarks/yak/{{sample}}.yak.tsv"
    conda: "envs/yak.yaml"
    threads: 28
    shell: "(yak count -t {threads} -o {output} {input}) > {log} 2>&1"


rule hifiasm_assemble:
    input: 
        fasta = lambda wildcards: expand(f"cohorts/{cohort}/fasta/{wildcards.sample}/{{movie}}.fasta", movie=ubam_fastq_dict[wildcards.sample]),
        parent1_yak = lambda wildcards: f"cohorts/{cohort}/yak/{trio_dict[wildcards.sample]['parent1']}.yak",
        parent2_yak = lambda wildcards: f"cohorts/{cohort}/yak/{trio_dict[wildcards.sample]['parent2']}.yak"
    output:
        temp(f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.hap1.p_ctg.gfa"),
        f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.hap1.p_ctg.lowQ.bed",
        f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.hap1.p_ctg.noseq.gfa",
        temp(f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.hap2.p_ctg.gfa"),
        f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.hap2.p_ctg.lowQ.bed",
        f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.hap2.p_ctg.noseq.gfa",
        temp(f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.p_utg.gfa"),
        f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.p_utg.lowQ.bed",
        f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.p_utg.noseq.gfa",
        temp(f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.r_utg.gfa"),
        f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.r_utg.lowQ.bed",
        f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.r_utg.noseq.gfa",
        temp(f"cohorts/{cohort}/hifiasm/{{sample}}.asm.ec.bin"),
        temp(f"cohorts/{cohort}/hifiasm/{{sample}}.asm.ovlp.reverse.bin"),
        temp(f"cohorts/{cohort}/hifiasm/{{sample}}.asm.ovlp.source.bin")
    log: f"cohorts/{cohort}/logs/hifiasm/{{sample}}.hifiasm.log"
    benchmark: f"cohorts/{cohort}/benchmarks/hifiasm/{{sample}}.hifiasm.tsv"
    conda: "envs/hifiasm.yaml"
    params: 
        prefix = f"cohorts/{cohort}/hifiasm/{{sample}}.asm",
        parent1 = lambda wildcards: trio_dict[wildcards.sample]['parent1'],
        parent2 = lambda wildcards: trio_dict[wildcards.sample]['parent2'],
        extra = "-c1 -d1"
    threads: 56
    shell:
            """
            (
                hifiasm -o {params.prefix} -t {threads} {params.extra} \
                    -1 {input.parent1_yak} -2 {input.parent2_yak} {input.fasta} \
                && (echo -e "hap1\t{params.parent1}\nhap2\t{params.parent2}" > cohorts/{cohort}/hifiasm/{wildcards.sample}.asm.key.txt) \
            ) > {log} 2>&1
            """
    

rule gfa2fa:
    input: f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.{{infix}}.gfa"
    output: f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta"
    log: f"cohorts/{cohort}/logs/gfa2fa/{{sample}}.asm.dip.{{infix}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/gfa2fa/{{sample}}.asm.dip.{{infix}}.tsv"
    conda: "envs/gfatools.yaml"
    shell: "(gfatools gfa2fa {input} > {output}) 2> {log}"


rule bgzip_fasta:
    input: f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta"
    output: f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta.gz"
    log: f"cohorts/{cohort}/logs/bgzip/{{sample}}.asm.dip.{{infix}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/bgzip/{{sample}}.asm.dip.{{infix}}.tsv"
    threads: 4
    conda: "envs/htslib.yaml"
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


rule yak_trioeval:
    input: 
        parent1_yak = lambda wildcards: f"cohorts/{cohort}/yak/{trio_dict[wildcards.sample]['parent1']}.yak",
        parent2_yak = lambda wildcards: f"cohorts/{cohort}/yak/{trio_dict[wildcards.sample]['parent2']}.yak",
        fasta = f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta.gz"
    output: f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta.trioeval.txt"
    log: f"cohorts/{cohort}/logs/yak/{{sample}}.asm.dip.{{infix}}.fasta.trioeval.log"
    conda: "envs/yak.yaml"
    threads: 16
    shell: "(yak trioeval -t {threads} {input.parent1_yak} {input.parent2_yak} {input.fasta} > {output}) > {log} 2>&1"


rule asm_stats:
    input: f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta.gz"
    output: f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta.stats.txt"
    log: f"cohorts/{cohort}/logs/asm_stats/{{sample}}.asm.dip.{{infix}}.fasta.log"
    benchmark: f"cohorts/{cohort}/benchmarks/asm_stats/{{sample}}.asm.dip.{{infix}}.fasta.tsv"
    conda: "envs/k8.yaml"
    shell: f"(k8 workflow/scripts/calN50/calN50.js -f {config['ref']['index']} {{input}} > {{output}}) > {{log}} 2>&1"


rule align_hifiasm:
    input:
        target = config['ref']['fasta'],
        query = [f"cohorts/{cohort}/hifiasm/{{sample}}.asm.dip.{infix}.fasta.gz"
                 for infix in ["hap1.p_ctg", "hap2.p_ctg"]]
    output: f"cohorts/{cohort}/hifiasm/{{sample}}.asm.{ref}.bam"
    log: f"cohorts/{cohort}/logs/align_hifiasm/{{sample}}.asm.{ref}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/align_hifiasm/{{sample}}.asm.{ref}.tsv"
    params:
        minimap2_args = "-L --secondary=no --eqx -ax asm5",
        minimap2_threads = 12,
        readgroup = f"@RG\\tID:{{sample}}_hifiasm\\tSM:{{sample}}",
        samtools_threads = 3,
        samtools_mem = "8G"
    threads: 16  # minimap2 + samtools(+3)
    conda: "envs/align_hifiasm.yaml"
    shell:
        """
        (minimap2 -t {params.minimap2_threads} {params.minimap2_args} -R '{params.readgroup}' {input.target} {input.query} \
            | samtools sort -@ {params.samtools_threads} -T $TMPDIR -m {params.samtools_mem} > {output}) > {log} 2>&1
        """


rule samtools_index_bam:
    input: f"cohorts/{cohort}/hifiasm/{{prefix}}.bam"
    output: f"cohorts/{cohort}/hifiasm/{{prefix}}.bam.bai"
    log: f"cohorts/{cohort}/logs/samtools/index/{{prefix}}.log"
    benchmark: f"cohorts/{cohort}/logs/samtools/index/{{prefix}}.tsv"
    threads: 4
    conda: "envs/samtools.yaml"
    shell: "(samtools index -@ 3 {input}) > {log} 2>&1"
