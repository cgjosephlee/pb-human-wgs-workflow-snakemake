# localrules: infer_sex_from_coverage, calculate_m2_ratio, calculate_gc_coverage


rule infer_sex_from_coverage:
    input: f"samples/{{sample}}/mosdepth/{{movie}}.{ref}.mosdepth.summary.txt"
    output: f"samples/{{sample}}/mosdepth/{{movie}}.{ref}.mosdepth.inferred_sex.txt"
    log: f"samples/{{sample}}/logs/quality_control/infer_sex_from_coverage.{{movie}}.{ref}.log"
    benchmark: f"samples/{{sample}}/benchmarks/quality_control/infer_sex_from_coverage.{{movie}}.{ref}.tsv"
    conda: "envs/pandas.yaml"
    shell: "(python3 workflow/scripts/infer_sex_from_coverage.py {input} > {output}) > {log} 2>&1"


rule calculate_m2_ratio:
    input: f"samples/{{sample}}/mosdepth/{{movie}}.{ref}.mosdepth.summary.txt"
    output: f"samples/{{sample}}/mosdepth/{{movie}}.{ref}.mosdepth.M2_ratio.txt"
    log: f"samples/{{sample}}/logs/quality_control/calculate_M2_ratio.{{movie}}.{ref}.log"
    benchmark: f"samples/{{sample}}/benchmarks/quality_control/calculate_M2_ratio.{{movie}}.{ref}.tsv"
    conda: "envs/pandas.yaml"
    shell: "(python3 workflow/scripts/calculate_M2_ratio.py {input} > {output}) > {log} 2>&1"


rule calculate_gc_coverage:
    input:
        mosdepth_regions = f"samples/{{sample}}/mosdepth/{{movie}}.{ref}.regions.bed.gz",
        ref = config['ref']['fasta']
    output: f"samples/{{sample}}/mosdepth/{{movie}}.{ref}.gc_coverage.summary.txt"
    log: f"samples/{{sample}}/logs/quality_control/calculate_gc_coverage.{{movie}}.{ref}.log"
    benchmark: f"samples/{{sample}}/benchmarks/quality_control/calculate_gc_coverage.{{movie}}.{ref}.tsv"
    conda: "envs/gc_coverage.yaml"
    shell:
        """
        (bedtools nuc -fi {input.ref} -bed {input.mosdepth_regions} \
            | awk '($11==0) {{ print 0.05*int($6/0.05) "\t" $4; }}' \
            | sort -k1,1g \
            | datamash -g1 q1 2 median 2 q3 2 count 2 \
            | tr '\t' ',' \
            | awk 'BEGIN {{ print "#gc_perc,q1,median,q3,count"; }} {{ print $0; }}' > {output}) > {log} 2>&1
        """
