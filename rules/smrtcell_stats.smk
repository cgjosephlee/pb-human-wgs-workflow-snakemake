ruleorder: smrtcell_stats_ubam > smrtcell_stats_fastq
# localrules: smrtcell_summary_stats


rule smrtcell_stats_ubam:
    input: lambda wildcards: ubam_dict[wildcards.sample][wildcards.movie]
    output: "samples/{sample}/smrtcell_stats/{movie}.read_length_and_quality.tsv"
    log: "samples/{sample}/logs/smrtcell_stats/{movie}.log"
    benchmark: "samples/{sample}/benchmarks/smrtcell_stats/{movie}.tsv"
    conda: "envs/smrtcell_stats.yaml"
    shell: "(python3 workflow/scripts/extract_read_length_and_qual.py {input} > {output}) > {log} 2>&1"


rule smrtcell_stats_fastq:
    input: lambda wildcards: fastq_dict[wildcards.sample][wildcards.movie]
    output: "samples/{sample}/smrtcell_stats/{movie}.read_length_and_quality.tsv"
    log: "samples/{sample}/logs/smrtcell_stats/{movie}.log"
    benchmark: "samples/{sample}/benchmarks/smrtcell_stats/{movie}.tsv"
    conda: "envs/smrtcell_stats.yaml"
    shell: "(python3 workflow/scripts/extract_read_length_and_qual.py {input} > {output}) > {log} 2>&1"


rule smrtcell_summary_stats:
    input: "samples/{sample}/smrtcell_stats/{movie}.read_length_and_quality.tsv"
    output:
        rlsummary = "samples/{sample}/smrtcell_stats/{movie}.read_length_summary.tsv",
        rqsummary = "samples/{sample}/smrtcell_stats/{movie}.read_quality_summary.tsv"
    log: "samples/{sample}/logs/smrtcell_stats/{movie}.summary.log"
    benchmark: "samples/{sample}/benchmarks/smrtcell_stats/{movie}.summary.tsv"
    conda: "envs/smrtcell_stats.yaml"
    shell:
        """
        (awk '{{ b=int($2/1000); b=(b>39?39:b); print 1000*b "\t" $2; }}' {input} |
            sort -k1,1g | datamash -g 1 count 1 sum 2 |
            awk 'BEGIN {{ for(i=0;i<=39;i++) {{ print 1000*i"\t0\t0"; }} }} {{ print; }}' |
            sort -k1,1g | datamash -g 1 sum 2 sum 3 > {output.rlsummary}) 2> {log}

        (awk '{{ print ($3>50?50:$3) "\t" $2; }}' {input} |
            sort -k1,1g | datamash -g 1 count 1 sum 2 |
            awk 'BEGIN {{ for(i=0;i<=60;i++) {{ print i"\t0\t0"; }} }} {{ print; }}' |
            sort -k1,1g | datamash -g 1 sum 2 sum 3 > {output.rqsummary}) 2>> {log}
        """
