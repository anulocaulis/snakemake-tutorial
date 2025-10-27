
# define functions. This is a key step to fix Lowry-nitrogen snakefile,
# because inputs for some rules can't be determined during initialization.
def get_bwa_map_inputs(wildcards):
    return config["samples"][wildcards.sample]

def get_read_group(wildcards):
    return f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}"

# map reads to reference genome with bwa mem
rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_inputs # function defined above to get sample-specific input correctly from wildcards
    threads: 8
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg=get_read_group # function defined above to get read group string
    log: 
        "logs/bwa_map/{sample}.log"
    benchmark:
        "benchmarks/{sample}_bwa_map.benchmark.txt"
    shell:
        "bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb - > {output} 2> {log}"

# sort BAM files with samtools sort
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample}"
        " -o {output} {input}"

# index sorted BAM files
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index {input}"
# after adding this rule, we can construct the dag with:
# snakemake --dag dot sorted_reads/A.bam | dot -Tpng > dag.png