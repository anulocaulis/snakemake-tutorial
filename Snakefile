# add config file
configfile: "config.yaml"

# include rules from other files
include: "mapping/rules_mapping.smk"

# Define samples
# moved to config.yaml

# rule all
rule all:
    input:
        "plots/quals.svg"


# call genomic variants with bcftools mpileup
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    params:
        mutation_rate=config["mutation_prior"]
    log:
        "logs/bcftools_call/all.log"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv -P {params.mutation_rate} - > {output} 2> {log}"

#### done for today at the end of step 5 from
#### https://hbctraining.github.io/Intro-to-NGS-hpc-O2/lessons/05b_variant_calling.html
#### https://snakemake.readthedocs.io/en/v8.1.1/tutorial/basics.html

# Using custom scripts
rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"

# I added a new dir called scripts and put the plot-quals.py file in there
# The content of plot-quals.py generates histograms of variant quality scores
# using matplotlib and pysam, and calling the input and output files via
# snakemake.input and snakemake.output

# add a target rule to the top of the Snakefile
# then run the entire workflow with:
# snakemake --use-conda --cores all plots/quals.svg

    
# done with basic tutorial. on to advanced tutorial next.
# I will update the above rules as needed based on the advanced tutorial.
# added config files that create params that can get called in rules, input functions
# to get sample-specific inputs, and read group functions to get sample-specific read groups
# also added logs to rules to capture stdout and stderr
# adding temporary files to reduce storage use now
# also added protected output to prevent accidental deletion of important files