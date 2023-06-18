# localrules: deepvariant_bcftools_stats, deepvariant_bcftools_roh


shards = [f"{x:05}" for x in range(config['N_SHARDS'])]

# variables affected by config['cpu_only']
deepvariant_version = config['deepvariant_cpu_version'] if config['cpu_only'] else config['deepvariant_gpu_version']
call_variants_threads = 256 if config['cpu_only'] else 8
call_variants_extra = os.environ.get('DEEPVARIANT_CPU_EXTRA', '') if config['cpu_only'] else os.environ['DEEPVARIANT_GPU_EXTRA']
call_variants_partition =  os.environ['PARTITION'] if config['cpu_only'] else os.environ['DEEPVARIANT_GPU_PARTITION']


rule deepvariant_make_examples:
    input:
        bams = abams,
        bais = [f"{x}.bai" for x in abams],
        reference = config['ref']['fasta']
    output:
        tfrecord = temp(f"samples/{sample}/deepvariant/examples/examples.tfrecord-{{shard}}-of-{config['N_SHARDS']:05}.gz"),
        nonvariant_site_tfrecord = temp(f"samples/{sample}/deepvariant/examples/gvcf.tfrecord-{{shard}}-of-{config['N_SHARDS']:05}.gz")
    log: f"samples/{sample}/logs/deepvariant/make_examples/{sample}.{ref}.{{shard}}-of-{config['N_SHARDS']:05}.log"
    benchmark: f"samples/{sample}/benchmarks/deepvariant/make_examples/{sample}.{ref}.{{shard}}-of-{config['N_SHARDS']:05}.tsv"
    container: f"docker://google/deepvariant:{deepvariant_version}"
    params:
        vsc_min_fraction_indels = "0.12",
        shard = lambda wildcards: wildcards.shard,
        reads = ','.join(abams)
    resources: 
        extra = os.environ.get('DEEPVARIANT_AVX2_CONSTRAINT', '')
    shell:
        f"""
        (/opt/deepvariant/bin/make_examples \
            --norealign_reads \
            --vsc_min_fraction_indels {{params.vsc_min_fraction_indels}} \
            --pileup_image_width 199 \
            --track_ref_reads \
            --phase_reads \
            --partition_size=25000 \
            --max_reads_per_partition=600 \
            --alt_aligned_pileup=diff_channels \
            --add_hp_channel \
            --sort_by_haplotypes \
            --parse_sam_aux_fields \
            --min_mapping_quality=1 \
            --mode calling \
            --ref {{input.reference}} \
            --reads {{params.reads}} \
            --examples samples/{sample}/deepvariant/examples/examples.tfrecord@{config['N_SHARDS']}.gz \
            --gvcf samples/{sample}/deepvariant/examples/gvcf.tfrecord@{config['N_SHARDS']}.gz \
            --task {{wildcards.shard}}) > {{log}} 2>&1
        """


rule deepvariant_call_variants:
    input: expand(f"samples/{sample}/deepvariant/examples/examples.tfrecord-{{shard}}-of-{config['N_SHARDS']:05}.gz", shard=shards)
    output: temp(f"samples/{sample}/deepvariant/{sample}.{ref}.call_variants_output.tfrecord.gz")
    log: f"samples/{sample}/logs/deepvariant/call_variants/{sample}.{ref}.log"
    benchmark: f"samples/{sample}/benchmarks/deepvariant/call_variants/{sample}.{ref}.tsv"
    container: f"docker://google/deepvariant:{deepvariant_version}"
    params: model = "/opt/models/pacbio/model.ckpt"
    threads: call_variants_threads
    resources:
        partition = call_variants_partition,
        extra = call_variants_extra
    shell:
        f"""
        (/opt/deepvariant/bin/call_variants \
            --outfile {{output}} \
            --examples samples/{sample}/deepvariant/examples/examples.tfrecord@{config['N_SHARDS']}.gz \
            --checkpoint {{params.model}}) > {{log}} 2>&1
        """


rule deepvariant_postprocess_variants:
    input:
        tfrecord = f"samples/{sample}/deepvariant/{sample}.{ref}.call_variants_output.tfrecord.gz",
        nonvariant_site_tfrecord = expand(f"samples/{sample}/deepvariant/examples/gvcf.tfrecord-{{shard:05}}-of-{config['N_SHARDS']:05}.gz",
                                          shard=range(config['N_SHARDS'])),
        reference = config['ref']['fasta']
    output:
        vcf = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz",
        vcf_index = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz.tbi",
        gvcf = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.g.vcf.gz",
        gvcf_index = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.g.vcf.gz.tbi",
        report = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.visual_report.html"
    log: f"samples/{sample}/logs/deepvariant/postprocess_variants/{sample}.{ref}.log"
    benchmark: f"samples/{sample}/benchmarks/deepvariant/postprocess_variants/{sample}.{ref}.tsv"
    container: f"docker://google/deepvariant:{deepvariant_version}"
    threads: 6
    resources: 
        extra = os.environ.get('DEEPVARIANT_AVX2_CONSTRAINT','')
    shell:
        f"""
        (/opt/deepvariant/bin/postprocess_variants \
            --ref {{input.reference}} \
            --infile {{input.tfrecord}} \
            --outfile {{output.vcf}} \
            --nonvariant_site_tfrecord_path samples/{sample}/deepvariant/examples/gvcf.tfrecord@{config['N_SHARDS']}.gz \
            --gvcf_outfile {{output.gvcf}}) > {{log}} 2>&1
        """


rule deepvariant_bcftools_stats:
    input: f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz"
    output: f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.stats.txt"
    log: f"samples/{sample}/logs/bcftools/stats/{sample}.{ref}.deepvariant.vcf.log"
    benchmark: f"samples/{sample}/benchmarks/bcftools/stats/{sample}.{ref}.deepvariant.vcf.tsv"
    params: f"--fasta-ref {config['ref']['fasta']} --apply-filters PASS -s {sample}"
    threads: 4
    conda: "envs/bcftools.yaml"
    shell: "(bcftools stats --threads 3 {params} {input} > {output}) > {log} 2>&1"


rule deepvariant_bcftools_roh:
    input: f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz"
    output: f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.roh.bed"
    log: f"samples/{sample}/logs/bcftools/stats/{sample}.{ref}.deepvariant.vcf.log"
    benchmark: f"samples/{sample}/benchmarks/bcftools/stats/{sample}.{ref}.deepvariant.vcf.tsv"
    params: default_allele_frequency = 0.4
    conda: "envs/bcftools.yaml"
    shell:
        """
        (echo -e "#chr\tstart\tend\tqual" > {output}
        bcftools roh --AF-dflt {params.default_allele_frequency} {input} \
        | awk -v OFS='\t' '$1=="RG" {{ print $3, $4, $5, $8 }}' \
        >> {output}) > {log} 2>&1
        """
