# localrules: calculate_sample_gc_coverage


rule calculate_sample_gc_coverage:
    input:
        mosdepth_regions = f"samples/{{sample}}/mosdepth/{{sample}}.{ref}.deepvariant.haplotagged.regions.bed.gz",
        ref = config['ref']['fasta']
    output: f"samples/{{sample}}/mosdepth/{{sample}}.{ref}.gc_coverage.summary.txt"
    log: f"samples/{{sample}}/logs/quality_control/calculate_gc_coverage.{{sample}}.{ref}.log"
    benchmark: f"samples/{{sample}}/benchmarks/quality_control/calculate_gc_coverage.{{sample}}.{ref}.tsv"
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
