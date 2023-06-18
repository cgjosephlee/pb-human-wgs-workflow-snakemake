# localrules: slivar_tsv


rule bcftools_norm:
    input:
        vcf = slivar_input,
        tbi = slivar_input + ".tbi"
    output: temp(f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.norm.bcf")
    log: f"cohorts/{cohort}/logs/bcftools/norm/{cohort}.{ref}.deepvariant.phased.vcf.log"
    benchmark: f"cohorts/{cohort}/benchmarks/bcftools/norm/{cohort}.{ref}.deepvariant.phased.vcf.tsv"
    params: f"--multiallelics - --output-type b --fasta-ref {config['ref']['fasta']}"
    conda: "envs/bcftools.yaml"
    shell: "(bcftools norm {params} {input.vcf} | bcftools sort --output-type b -o {output}) > {log} 2>&1"


slivar_filters = [
        f"""--info 'variant.FILTER==\"PASS\" \
                && INFO.gnomad_af <= {config['max_gnomad_af']} \
                && INFO.hprc_af <= {config['max_hprc_af']} \
                && INFO.gnomad_nhomalt <= {config['max_gnomad_nhomalt']} \
                && INFO.hprc_nhomalt <= {config['max_hprc_nhomalt']}'""",
        "--family-expr 'recessive:fam.every(segregating_recessive)'",
        "--family-expr 'x_recessive:(variant.CHROM == \"chrX\") && fam.every(segregating_recessive_x)'",
        f"""--family-expr 'dominant:fam.every(segregating_dominant) \
                       && INFO.gnomad_ac <= {config['max_gnomad_ac']} \
                       && INFO.hprc_ac <= {config['max_hprc_ac']}'""",
        f"""--family-expr 'x_dominant:(variant.CHROM == \"chrX\") \
                       && fam.every(segregating_dominant_x) \
                       && INFO.gnomad_ac <= {config['max_gnomad_ac']} \
                       && INFO.hprc_ac <= {config['max_hprc_ac']}'""",
        f"--sample-expr 'comphet_side:sample.het && sample.GQ > {config['min_gq']}'",
]


rule slivar_small_variant:
    input:
        bcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.norm.bcf",
        csi = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.norm.bcf.csi",
        ped = f"cohorts/{cohort}/{cohort}.ped",
        gnomad_af = {config['ref']['gnomad_gnotate']},
        hprc_af = {config['ref']['hprc_dv_gnotate']},
        js = config['slivar_js'],
        gff = config['ref']['ensembl_gff'],
        ref = config['ref']['fasta']
    output: f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.vcf"
    log: f"cohorts/{cohort}/logs/slivar/filter/{cohort}.{ref}.deepvariant.phased.slivar.vcf.log"
    benchmark: f"cohorts/{cohort}/benchmarks/slivar/filter/{cohort}.{ref}.deepvariant.phased.slivar.tsv"
    params: filters = slivar_filters
    threads: 12
    conda: "envs/slivar.yaml"
    shell:
        """
        (pslivar --processes {threads} \
            --fasta {input.ref}\
            --pass-only \
            --js {input.js} \
            {params.filters} \
            --gnotate {input.gnomad_af} \
            --gnotate {input.hprc_af} \
            --vcf {input.bcf} \
            --ped {input.ped} \
            | bcftools csq -l -s - --ncsq 40 \
                -g {input.gff} -f {input.ref} - -o {output}) > {log} 2>&1
        """


skip_list = [
    'non_coding_transcript',
    'intron',
    'non_coding',
    'upstream_gene',
    'downstream_gene',
    'non_coding_transcript_exon',
    'NMD_transcript',
    '5_prime_UTR',
    '3_prime_UTR'
    ]


rule slivar_compound_hets:
    input: 
        vcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.vcf.gz",
        tbi = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.vcf.gz.tbi",
        ped = f"cohorts/{cohort}/{cohort}.ped"
    output:
        vcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.compound-hets.vcf"
    log: f"cohorts/{cohort}/logs/slivar/compound-hets/{cohort}.{ref}.deepvariant.phased.slivar.compound-hets.vcf.log"
    benchmark: f"cohorts/{cohort}/benchmarks/slivar/compound-hets/{cohort}.{ref}.deepvariant.phased.slivar.compound-hets.vcf.tsv"
    params: skip = ",".join(skip_list)
    conda: "envs/slivar.yaml"
    shell:
        """
        (slivar compound-hets \
            --skip {params.skip} \
            --vcf {input.vcf} \
            --sample-field comphet_side \
            --ped {input.ped} \
            --allow-non-trios \
            | python3 workflow/scripts/add_comphet_phase.py \
            > {output.vcf}) > {log} 2>&1
        """


info_fields = [
    'gnomad_af',
    'hprc_af',
    'gnomad_nhomalt',
    'hprc_nhomalt',
    'gnomad_ac',
    'hprc_ac'
    ]


rule slivar_tsv:
    input:
        filt_vcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.vcf.gz",
        comphet_vcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.compound-hets.vcf.gz",
        ped = f"cohorts/{cohort}/{cohort}.ped",
        lof_lookup = config['lof_lookup'],
        clinvar_lookup = config['clinvar_lookup'],
        phrank_lookup = f"cohorts/{cohort}/{cohort}_phrank.tsv"
    output:
        filt_tsv = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.tsv",
        comphet_tsv = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.compound-hets.tsv"
    log: f"cohorts/{cohort}/logs/slivar/tsv/{cohort}.{ref}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/slivar/tsv/{cohort}.{ref}.tsv"
    params: info = "".join([f"--info-field {x} " for x in info_fields])
    conda: "envs/slivar.yaml"
    shell:
        """
        (slivar tsv \
            {params.info} \
            --sample-field dominant \
            --sample-field x_dominant \
            --sample-field recessive \
            --sample-field x_recessive \
            --csq-field BCSQ \
            --gene-description {input.lof_lookup} \
            --gene-description {input.clinvar_lookup} \
            --gene-description {input.phrank_lookup} \
            --ped {input.ped} \
            --out /dev/stdout \
            {input.filt_vcf} \
            | sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
            > {output.filt_tsv}
        slivar tsv \
            {params.info} \
            --sample-field slivar_comphet \
            --info-field slivar_comphet \
            --csq-field BCSQ \
            --gene-description {input.lof_lookup} \
            --gene-description {input.clinvar_lookup} \
            --gene-description {input.phrank_lookup} \
            --ped {input.ped} \
            --out /dev/stdout \
            {input.comphet_vcf} \
            | sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
            > {output.comphet_tsv}) > {log} 2>&1
        """
