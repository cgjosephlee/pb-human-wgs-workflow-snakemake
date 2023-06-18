ruleorder: trgt_genotype > bgzip_vcf


rule trgt_genotype:
    input:
        reference = config['ref']['fasta'],
        bam = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam",
        bai = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam.bai",
        bed = config['ref']['trgt_bed'],
    output:
        vcf = f"samples/{sample}/trgt/{sample}.{ref}.trgt.vcf.gz",
        bam = f"samples/{sample}/trgt/{sample}.{ref}.trgt.spanning.bam",
    log: f"samples/{sample}/logs/trgt/genotype.log"
    benchmark: f"samples/{sample}/benchmarks/trgt/genotype.tsv"
    conda: "envs/trgt.yaml"
    params:
        prefix = f"samples/{sample}/trgt/{sample}.{ref}.trgt"
    threads: 28
    message: "Executing {rule}: Genotyping tandem repeat regions from {input.bed} in {input.bam}."
    shell:
        """
        (trgt --threads {threads} \
            --genome {input.reference} \
            --repeats {input.bed} \
            --reads {input.bam} \
            --output-prefix {params.prefix}) > {log} 2>&1
        """


rule trgt_coverage_dropouts:
    input:
        bam = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam",
        bai = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam.bai",
        bed = config['ref']['trgt_bed']
    output: f"samples/{sample}/trgt/{sample}.{ref}.trgt.dropouts.txt"
    log: f"samples/{sample}/logs/trgt/{sample}.dropouts.log"
    benchmark: f"samples/{sample}/benchmarks/trgt/{sample}.dropouts.tsv"
    conda: "envs/tandem-genotypes.yaml"
    message: "Executing {rule}: Identify coverage dropouts in {input.bed} regions in {input.bam}."
    shell: "(python3 workflow/scripts/check_trgt_coverage.py {input.bed} {input.bam} > {output}) > {log} 2>&1"
