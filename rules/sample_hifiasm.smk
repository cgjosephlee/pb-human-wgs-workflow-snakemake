ruleorder: samtools_fasta > seqtk_fastq_to_fasta
localrules: asm_stats, htsbox_bcftools_stats, gfa2fa


rule samtools_fasta:
    input: lambda wildcards: ubam_dict[wildcards.movie]
    output: temp(f"samples/{sample}/fasta/{{movie}}.fasta")
    log: f"samples/{sample}/logs/samtools/fasta/{{movie}}.log"
    benchmark: f"samples/{sample}/benchmarks/samtools/fasta/{{movie}}.tsv"
    threads: 4
    conda: "envs/samtools.yaml"
    shell: "(samtools fasta -@ 3 {input} > {output}) > {log} 2>&1"


rule seqtk_fastq_to_fasta:
    input: lambda wildcards: fastq_dict[wildcards.movie]
    output: temp(f"samples/{sample}/fasta/{{movie}}.fasta")
    log: f"samples/{sample}/logs/seqtk/seq/{{movie}}.log"
    benchmark: f"samples/{sample}/benchmarks/seqtk/seq/{{movie}}.tsv"
    conda: "envs/seqtk.yaml"
    shell: "(seqtk seq -A {input} > {output}) > {log} 2>&1"


rule hifiasm_assemble:
    input: expand(f"samples/{sample}/fasta/{{movie}}.fasta", movie=movies)
    output:
        temp(f"samples/{sample}/hifiasm/{sample}.asm.bp.hap1.p_ctg.gfa"),
        f"samples/{sample}/hifiasm/{sample}.asm.bp.hap1.p_ctg.lowQ.bed",
        f"samples/{sample}/hifiasm/{sample}.asm.bp.hap1.p_ctg.noseq.gfa",
        temp(f"samples/{sample}/hifiasm/{sample}.asm.bp.hap2.p_ctg.gfa"),
        f"samples/{sample}/hifiasm/{sample}.asm.bp.hap2.p_ctg.lowQ.bed",
        f"samples/{sample}/hifiasm/{sample}.asm.bp.hap2.p_ctg.noseq.gfa",
        temp(f"samples/{sample}/hifiasm/{sample}.asm.bp.p_ctg.gfa"),
        temp(f"samples/{sample}/hifiasm/{sample}.asm.bp.p_utg.gfa"),
        temp(f"samples/{sample}/hifiasm/{sample}.asm.bp.r_utg.gfa"),
        temp(f"samples/{sample}/hifiasm/{sample}.asm.ec.bin"),
        temp(f"samples/{sample}/hifiasm/{sample}.asm.ovlp.reverse.bin"),
        temp(f"samples/{sample}/hifiasm/{sample}.asm.ovlp.source.bin")
    log: f"samples/{sample}/logs/hifiasm.log"
    benchmark: f"samples/{sample}/benchmarks/hifiasm.tsv"
    conda: "envs/hifiasm.yaml"
    params: prefix = f"samples/{sample}/hifiasm/{sample}.asm"
    threads: 56
    shell: "(hifiasm -o {params.prefix} -t {threads} {input}) > {log} 2>&1"


rule gfa2fa:
    input: f"samples/{sample}/hifiasm/{sample}.asm.bp.{{infix}}.gfa"
    output: f"samples/{sample}/hifiasm/{sample}.asm.bp.{{infix}}.fasta"
    log: f"samples/{sample}/logs/gfa2fa/{sample}.asm.bp.{{infix}}.log"
    benchmark: f"samples/{sample}/benchmarks/gfa2fa/{sample}.asm.bp.{{infix}}.tsv"
    conda: "envs/gfatools.yaml"
    shell: "(gfatools gfa2fa {input} > {output}) 2> {log}"


rule bgzip_fasta:
    input: f"samples/{sample}/hifiasm/{sample}.asm.bp.{{infix}}.fasta"
    output: f"samples/{sample}/hifiasm/{sample}.asm.bp.{{infix}}.fasta.gz"
    log: f"samples/{sample}/logs/bgzip/{sample}.asm.bp.{{infix}}.log"
    benchmark: f"samples/{sample}/benchmarks/bgzip/{sample}.asm.bp.{{infix}}.tsv"
    threads: 4
    conda: "envs/htslib.yaml"
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


rule asm_stats:
    input: f"samples/{sample}/hifiasm/{sample}.asm.bp.{{infix}}.fasta.gz"
    output: f"samples/{sample}/hifiasm/{sample}.asm.bp.{{infix}}.fasta.stats.txt"
    log: f"samples/{sample}/logs/asm_stats/{sample}.asm.bp.{{infix}}.fasta.log"
    benchmark: f"samples/{sample}/benchmarks/asm_stats/{sample}.asm.bp.{{infix}}.fasta.tsv"
    conda: "envs/k8.yaml"
    shell: f"(k8 workflow/scripts/calN50/calN50.js -f {config['ref']['index']} {{input}} > {{output}}) > {{log}} 2>&1"


rule align_hifiasm:
    input:
        target = config['ref']['fasta'],
        query = [f"samples/{sample}/hifiasm/{sample}.asm.bp.{infix}.fasta.gz"
                 for infix in ["hap1.p_ctg", "hap2.p_ctg"]]
    output: f"samples/{sample}/hifiasm/{sample}.asm.{ref}.bam"
    log: f"samples/{sample}/logs/align_hifiasm/{sample}.asm.{ref}.log"
    benchmark: f"samples/{sample}/benchmarks/align_hifiasm/{sample}.asm.{ref}.tsv"
    params:
        minimap2_args = "-L --secondary=no --eqx -ax asm5",
        minimap2_threads = 12,
        readgroup = f"@RG\\tID:{sample}_hifiasm\\tSM:{sample}",
        samtools_threads = 3,
        samtools_mem = "8G"
    threads: 16  # minimap2 + samtools(+3)
    conda: "envs/align_hifiasm.yaml"
    shell:
        """
        (minimap2 -t {params.minimap2_threads} {params.minimap2_args} -R '{params.readgroup}' {input.target} {input.query} \
            | samtools sort -@ {params.samtools_threads} -T $TMPDIR -m {params.samtools_mem} > {output}) > {log} 2>&1
        """


rule htsbox:
    input:
        bam = f"samples/{sample}/hifiasm/{sample}.asm.{ref}.bam",
        bai = f"samples/{sample}/hifiasm/{sample}.asm.{ref}.bam.bai",
        reference = config['ref']['fasta']
    output: f"samples/{sample}/hifiasm/{sample}.asm.{ref}.htsbox.vcf"
    log: f"samples/{sample}/logs/htsbox/{sample}.asm.log"
    benchmark: f"samples/{sample}/benchmarks/htsbox/{sample}.asm.tsv"
    params: '-q20'
    conda: "envs/htsbox.yaml"
    shell: "(htsbox pileup {params} -c -f {input.reference} {input.bam} > {output})> {log} 2>&1"


rule htsbox_bcftools_stats:
    input: f"samples/{sample}/hifiasm/{sample}.asm.{ref}.htsbox.vcf.gz"
    output: f"samples/{sample}/hifiasm/{sample}.asm.{ref}.htsbox.vcf.stats.txt"
    log: f"samples/{sample}/logs/bcftools/stats/{sample}.asm.{ref}.htsbox.vcf.log"
    benchmark: f"samples/{sample}/benchmarks/bcftools/stats/{sample}.asm.{ref}.htsbox.vcf.tsv"
    params: f"--fasta-ref {config['ref']['fasta']} -s samples/{sample}/hifiasm/{sample}.asm.{ref}.bam"
    threads: 4
    conda: "envs/bcftools.yaml"
    shell: "(bcftools stats --threads 3 {params} {input} > {output}) > {log} 2>&1"
