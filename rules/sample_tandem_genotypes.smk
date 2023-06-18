# localrules: tandem_genotypes_absolute_count, tandem_genotypes_plot, tandem_repeat_coverage_dropouts


rule last_align:
    input:
        [f"{config['ref']['last_index']}.{suffix}"
         for suffix in ['bck', 'des', 'prj', 'sds', 'ssp', 'suf', 'tis']],
        bam = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam",
        bai = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam.bai",
        bed = config['ref']['tg_bed'],
        score_matrix = config['score_matrix']
    output: temp(f"samples/{sample}/tandem-genotypes/{sample}.maf.gz")
    log: f"samples/{sample}/logs/last/align/{sample}.log"
    benchmark: f"samples/{sample}/benchmarks/last/align/{sample}.tsv"
    conda: "envs/last.yaml"
    params: 
        last_index = config['ref']['last_index'],
        extra = "-C2"
    threads: 28
    shell:
        """
        (samtools view -@3 -bL {input.bed} {input.bam} | samtools fasta \
         | lastal -P20 -p {input.score_matrix} {params.extra} {params.last_index} - \
         | last-split | bgzip > {output}) > {log} 2>&1
        """


rule tandem_genotypes:
    input:
        maf = f"samples/{sample}/tandem-genotypes/{sample}.maf.gz",
        repeats = config['ref']['tg_list']
    output: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.txt"
    log: f"samples/{sample}/logs/tandem-genotypes/{sample}.log"
    benchmark: f"samples/{sample}/benchmarks/tandem-genotypes/{sample}.tsv"
    conda: "envs/tandem-genotypes.yaml"
    shell: "(tandem-genotypes {input.repeats} {input.maf} > {output}) > {log} 2>&1"


rule tandem_genotypes_absolute_count:
    input: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.txt"
    output: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.absolute.txt"
    log: f"samples/{sample}/logs/tandem-genotypes/{sample}.absolute.log"
    benchmark: f"samples/{sample}/benchmarks/tandem-genotypes/{sample}.absolute.tsv"
    shell:
        """
        (awk -v OFS='\t' \
            '$0 ~ /^#/ {{print $0 " modified by adding reference repeat count"}}
            $0 !~ /^#/ {{
                ref_count=int(($3-$2)/length($4));
                num_fwd=split($7, fwd, ",");
                num_rev=split($8, rev, ",");
                new_fwd=result=fwd[1] + ref_count;
                for (i=2; i<=num_fwd; i++)
                    new_fwd = new_fwd "," fwd[i] + ref_count;
                new_rev=rev[1] + ref_count;
                for (i=2; i<=num_rev; i++)
                    new_rev = new_rev "," rev[i] + ref_count;
                print $1, $2, $3, $4, $5, $6, new_fwd, new_rev;
            }}' {input} > {output} \
        ) > {log} 2>&1
        """


rule tandem_genotypes_plot:
    input: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.txt"
    output: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.pdf"
    log: f"samples/{sample}/logs/tandem-genotypes/plot/{sample}.log"
    benchmark: f"samples/{sample}/benchmarks/tandem-genotypes/plot/{sample}.tsv"
    conda: "envs/tandem-genotypes.yaml"
    params: top_N_plots = 100
    shell: "(tandem-genotypes-plot -n {params.top_N_plots} {input} {output}) > {log} 2>&1"


rule tandem_repeat_coverage_dropouts:
    input:
        bam = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam",
        bai = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam.bai",
        bed = config['ref']['tg_bed']
    output: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.dropouts.txt"
    log: f"samples/{sample}/logs/tandem-genotypes/{sample}.dropouts.log"
    benchmark: f"samples/{sample}/benchmarks/tandem-genotypes/{sample}.dropouts.tsv"
    conda: "envs/tandem-genotypes.yaml"
    shell: "(python3 workflow/scripts/check_tandem_repeat_coverage.py {input.bed} {input.bam} > {output}) > {log} 2>&1"
