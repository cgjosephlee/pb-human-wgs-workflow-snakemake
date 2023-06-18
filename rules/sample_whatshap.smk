# localrules: split_deepvariant_vcf, whatshap_bcftools_concat


rule split_deepvariant_vcf:
    input: f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz",
    output: temp(f"samples/{sample}/whatshap/{sample}.{ref}.regions/{sample}.{ref}.{{region}}.deepvariant.vcf")
    log: f"samples/{sample}/logs/tabix/query/{sample}.{ref}.{{region}}.deepvariant.vcf.log"
    benchmark: f"samples/{sample}/benchmarks/tabix/query/{sample}.{ref}.{{region}}.deepvariant.vcf.tsv"
    params: region = lambda wildcards: wildcards.region, extra = '-h'
    conda: "envs/htslib.yaml"
    shell: "tabix {params.extra} {input} {params.region} > {output} 2> {log}"


rule whatshap_phase:
    input:
        reference = config['ref']['fasta'],
        vcf = f"samples/{sample}/whatshap/{sample}.{ref}.regions/{sample}.{ref}.{{chromosome}}.deepvariant.vcf.gz",
        tbi = f"samples/{sample}/whatshap/{sample}.{ref}.regions/{sample}.{ref}.{{chromosome}}.deepvariant.vcf.gz.tbi",
        phaseinput = abams,
        phaseinputindex = [f"{x}.bai" for x in abams]
    output: temp(f"samples/{sample}/whatshap/{sample}.{ref}.regions/{sample}.{ref}.{{chromosome}}.deepvariant.phased.vcf.gz")
    log: f"samples/{sample}/logs/whatshap/phase/{sample}.{ref}.{{chromosome}}.log"
    benchmark: f"samples/{sample}/benchmarks/whatshap/phase/{sample}.{ref}.{{chromosome}}.tsv"
    params: chromosome = lambda wildcards: wildcards.chromosome, extra = "--indels"
    conda: "envs/whatshap.yaml"
    shell:
        """
        (whatshap phase {params.extra} \
            --chromosome {wildcards.chromosome} \
            --output {output} \
            --reference {input.reference} \
            {input.vcf} \
            {input.phaseinput}) > {log} 2>&1
        """


rule whatshap_bcftools_concat:
    input:
        calls = expand(f"samples/{sample}/whatshap/{sample}.{ref}.regions/{sample}.{ref}.{{region}}.deepvariant.phased.vcf.gz", region=all_chroms),
        indices = expand(f"samples/{sample}/whatshap/{sample}.{ref}.regions/{sample}.{ref}.{{region}}.deepvariant.phased.vcf.gz.tbi", region=all_chroms)
    output: f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz"
    log: f"samples/{sample}/logs/bcftools/concat/{sample}.{ref}.whatshap.log"
    benchmark: f"samples/{sample}/benchmarks/bcftools/concat/{sample}.{ref}.whatshap.tsv"
    params: "-a -Oz"
    conda: "envs/bcftools.yaml"
    shell: "bcftools concat {params} -o {output} {input.calls} > {log} 2>&1"


rule whatshap_stats:
    input:
        vcf = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz",
        tbi = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz.tbi",
        chr_lengths = config['ref']['chr_lengths']
    output:
        gtf = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.gtf",
        tsv = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.tsv",
        blocklist = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.blocklist"
    log: f"samples/{sample}/logs/whatshap/stats/{sample}.{ref}.log"
    benchmark: f"samples/{sample}/benchmarks/whatshap/stats/{sample}.{ref}.tsv"
    conda: "envs/whatshap.yaml"
    shell:
        """
        (whatshap stats \
            --gtf {output.gtf} \
            --tsv {output.tsv} \
            --block-list {output.blocklist} \
            --chr-lengths {input.chr_lengths} \
            {input.vcf}) > {log} 2>&1
        """


rule whatshap_haplotag:
    input:
        reference = config['ref']['fasta'],
        vcf = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz",
        tbi = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz.tbi",
        bam = lambda wildcards: abam_dict[wildcards.movie]
    output: temp(f"samples/{sample}/whatshap/{sample}.{ref}.{{movie}}.deepvariant.haplotagged.bam")
    log: f"samples/{sample}/logs/whatshap/haplotag/{sample}.{ref}.{{movie}}.log"
    benchmark: f"samples/{sample}/benchmarks/whatshap/haplotag/{sample}.{ref}.{{movie}}.tsv"
    params: "--tag-supplementary"
    threads: 4
    conda: "envs/whatshap.yaml"
    shell:
        """
        (whatshap haplotag {params} \
            --output-threads {threads} \
            --output {output} \
            --reference {input.reference} \
            {input.vcf} {input.bam}) > {log} 2>&1
        """


rule merge_haplotagged_bams:
    input: expand(f"samples/{sample}/whatshap/{sample}.{ref}.{{movie}}.deepvariant.haplotagged.bam", movie=movies)
    output: f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam"
    log: f"samples/{sample}/logs/samtools/merge/{sample}.{ref}.haplotag.log"
    benchmark: f"samples/{sample}/benchmarks/samtools/merge/{sample}.{ref}.haplotag.tsv"
    threads: 8
    conda: "envs/samtools.yaml"
    shell:
        """
        if [[ "{input}" =~ " " ]]
        then
            (samtools merge -@ 7 {output} {input}) > {log} 2>&1
        else
            mv {input} {output}
        fi
        """
