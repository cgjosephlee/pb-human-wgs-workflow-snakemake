def list_svsigs(region):
    """Given a region, return all cohort svsig files from this region."""
    return [f"samples/{sample}/pbsv/svsig/{movie}.{ref}.{region}.svsig.gz" for movie in movies]


# localrules: bcftools_concat_pbsv_vcf


rule pbsv_discover:
    input:
        bam = f"samples/{sample}/aligned/{{movie}}.{ref}.bam",
        bai = f"samples/{sample}/aligned/{{movie}}.{ref}.bam.bai",
        tr_bed = config['ref']['tr_bed']
    output: f"samples/{sample}/pbsv/svsig/{{movie}}.{ref}.{{region}}.svsig.gz"
    log: f"samples/{sample}/logs/pbsv/discover/{{movie}}.{ref}.{{region}}.log"
    benchmark: f"samples/{sample}/benchmarks/pbsv/discover/{{movie}}.{ref}.{{region}}.tsv"
    params:
        extra = "--hifi",
        region = lambda wildcards: wildcards.region,
        loglevel = "INFO"
    conda: "envs/pbsv.yaml"
    shell:
        """
        (pbsv discover {params.extra} \
            --log-level {params.loglevel} \
            --region {wildcards.region} \
            --tandem-repeats {input.tr_bed} \
            {input.bam} {output}) > {log} 2>&1
        """


rule pbsv_call:
    input:
        svsigs = lambda wildcards: list_svsigs(wildcards.region),
        reference = config['ref']['fasta']
    output: temp(f"samples/{sample}/pbsv/{sample}.{ref}.chrom_vcfs/{sample}.{ref}.{{region}}.pbsv.vcf")
    log: f"samples/{sample}/logs/pbsv/call/{sample}.{ref}.{{region}}.log"
    benchmark: f"samples/{sample}/benchmarks/pbsv/call/{sample}.{ref}.{{region}}.tsv"
    params:
        region = lambda wildcards: wildcards.region,
        extra = "--hifi -m 20 " + config['pbsv_call_extra'],
        loglevel = "INFO"
    threads: 8
    conda: "envs/pbsv.yaml"
    shell:
        """
        (pbsv call {params.extra} \
            --log-level {params.loglevel} \
            --num-threads {threads} \
            {input.reference} {input.svsigs} {output}) > {log} 2>&1
        """


rule bcftools_concat_pbsv_vcf:
    input:
        calls = expand(f"samples/{sample}/pbsv/{sample}.{ref}.chrom_vcfs/{sample}.{ref}.{{region}}.pbsv.vcf.gz", region=all_chroms),
        indices = expand(f"samples/{sample}/pbsv/{sample}.{ref}.chrom_vcfs/{sample}.{ref}.{{region}}.pbsv.vcf.gz.tbi", region=all_chroms)
    output: f"samples/{sample}/pbsv/{sample}.{ref}.pbsv.vcf"
    log: f"samples/{sample}/logs/bcftools/concat/{sample}.{ref}.pbsv.vcf.log"
    benchmark: f"samples/{sample}/benchmarks/bcftools/concat/{sample}.{ref}.pbsv.vcf.tsv"
    conda: "envs/bcftools.yaml"
    shell: "(bcftools concat -a -o {output} {input.calls}) > {log} 2>&1"


# TODO: cleanup pbsv intermediates
