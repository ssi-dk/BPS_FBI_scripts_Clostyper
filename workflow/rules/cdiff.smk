"""
----------------------------
Friedrich-Loeffler-Institut (https://www.fli.de/), IBIZ
date: April, 27, 2022
Author: Mostafa Abdel-Glil (mostafa.abdel-glil@fli.de)
-------------------------------
"""
import os
import tempfile

#configfile: "config.yaml"
#working_dir=config['working_dir']
#snakemake_folder=config['snakemake_folder']

#directories
raw_data_dir= working_dir + "1_input_samples/"
fasta_dir= working_dir + "1_input_samples/"
results_per_sample_dir= working_dir + "2_results_per_sample/"
results_all_samples_dir= working_dir + "3_results_all_samples/"
assembly_dir_all= results_all_samples_dir + "assembly_fasta/"
reports_dir= working_dir + "4_all_reports/"

#samples
SAMPLES, = glob_wildcards( raw_data_dir + "{sample}_1.fastq.gz")
GENOMES, = glob_wildcards( fasta_dir + "{genome}.fasta")
SAMPLES = set(SAMPLES)
GENOMES = set(GENOMES)
DATASET =  [SAMPLES, GENOMES]

#scripts_paths
bin_dir= snakemake_folder + "../scripts/",
envs_folder= snakemake_folder + "../envs/"

#variables
#if "c_difficile_db" in config:
#    AMR_db_abricate= config["AMR_db_abricate"]
#else:
#    AMR_db_abricate="ncbi"

if "c_difficile_db" in config:
    reference= config["reference"]
    transposon_db_dir= config["c_difficile_db"]["transposons"]["database_dir"]
    transposon_db_name= config["c_difficile_db"]["transposons"]["database_name"]
    transposon_coverage= config["c_difficile_db"]["transposons"]["coverage"]
    plasmids_db_dir= config["c_difficile_db"]["plasmids"]["database_dir"]
    plasmids_db_name= config["c_difficile_db"]["plasmids"]["database_name"]
    plasmids_coverage= config["c_difficile_db"]["plasmids"]["coverage"]
    phages_db_dir= config["c_difficile_db"]["phages"]["database_dir"]
    phages_db_name= config["c_difficile_db"]["phages"]["database_name"]
    phages_coverage= config["c_difficile_db"]["phages"]["coverage"]
    res_mutations_db_dir= config["c_difficile_db"]["resitance_mutations"]["database_dir"]
    res_mutations_db_name= config["c_difficile_db"]["resitance_mutations"]["database_name"]
    res_mutations_aminoacid= config["c_difficile_db"]["resitance_mutations"]["aminoacid_position"]
    res_mutations_nucleotide= config["c_difficile_db"]["resitance_mutations"]["nucleotide_position"]
    pubmlst_loci= config["c_difficile_db"]["pubmlst_loci"]
else:
    c_difficile_db=""

#Rules output
cdiff_trsanposons= expand( results_per_sample_dir + "{sample}/cdiff_trsanposons/transposons.tab", sample=SAMPLES)
cdiff_plasmids= expand( results_per_sample_dir + "{sample}/cdiff_plasmids/plasmids.tab", sample=SAMPLES)
cdiff_phages= expand( results_per_sample_dir + "{sample}/cdiff_phages/phages.tab", sample=SAMPLES)
cdiff_resistance_mutations = expand( results_per_sample_dir + "{sample}/cdiff_resistance_mutations/aminoacid_at_mutated_positions.tab", sample=SAMPLES)

#master rule
#rule all:
#    input: cdiff_trsanposons, cdiff_plasmids, cdiff_phages, cdiff_resistance_mutations


"""
trsanposons
"""
rule trsanposons_search_with_srst2:
    input:
        r1 = raw_data_dir + "{sample}_1.fastq.gz",
        r2 = raw_data_dir + "{sample}_2.fastq.gz",
    output:
        srst2_summary = results_per_sample_dir + "{sample}/cdiff_trsanposons/transposons.tab",
        srst2_fullsummary = results_per_sample_dir + "{sample}/cdiff_trsanposons/transposons_full.tab",
        srst2_raw_pileup = temp(results_per_sample_dir + "{sample}/cdiff_trsanposons/srst2__{sample}.transposons_clustered.pileup"),
        srst2_raw_bam = temp(results_per_sample_dir + "{sample}/cdiff_trsanposons/srst2__{sample}.transposons_clustered.sorted.bam"),
    threads: 16
    conda:
        envs_folder + "srst2.yaml"
    params:
        trsanposons_dir = directory (results_per_sample_dir + "{sample}/cdiff_trsanposons/"),
        srst2_raw_summary = results_per_sample_dir + "{sample}/cdiff_trsanposons/srst2__genes__transposons_clustered__results.txt",
        srst2_raw_summary_full = results_per_sample_dir + "{sample}/cdiff_trsanposons/srst2__fullgenes__transposons_clustered__results.txt",
    shell:
        "srst2 --input_pe {input.r1} {input.r2} --output {params.trsanposons_dir}srst2 --log --gene_db {transposon_db_dir}{transposon_db_name} --threads {threads} --min_coverage {transposon_coverage} --stop_after 1000000 --report_new_consensus --report_all_consensus "
        "&& if [[ -f {params.srst2_raw_summary} ]]; then  mv {params.srst2_raw_summary} {output.srst2_summary}; else touch {output.srst2_summary}; fi " #incase of negative results create empty file
        "&& if [[ -f {params.srst2_raw_summary_full} ]]; then mv {params.srst2_raw_summary_full} {output.srst2_fullsummary}; else touch {output.srst2_fullsummary}; fi"
"""
collect all transposons results
"""
rule compile_transposons:
    input:
        trsanposons_results_per_sample = expand( results_per_sample_dir + "{sample}/cdiff_trsanposons/srst2__genes__transposons_clustered__results.txt", sample=SAMPLES),
    output:
        trsanposons_results_all_sample = results_all_samples_dir + "transposons_allsamples.tab",
    threads: 16
    conda:
        envs_folder + "srst2.yaml"
    shell:
        "srst2 --prev_output {input.trsanposons_results_per_sample} --output {output.trsanposons_results_all_sample}"

"""
plasmids
"""
rule plasmids_search_with_srst2:
    input:
        r1 = raw_data_dir + "{sample}_1.fastq.gz",
        r2 = raw_data_dir + "{sample}_2.fastq.gz",
    output:
        srst2_summary = results_per_sample_dir + "{sample}/cdiff_plasmids/plasmids.tab",
        srst2_fullsummary = results_per_sample_dir + "{sample}/cdiff_plasmids/plasmids_full.tab",
        srst2_raw_pileup = temp(results_per_sample_dir + "{sample}/cdiff_plasmids/srst2__{sample}.plasmids_clustered.pileup"),
        srst2_raw_bam = temp(results_per_sample_dir + "{sample}/cdiff_plasmids/srst2__{sample}.plasmids_clustered.sorted.bam"),
    threads: 16
    conda:
        envs_folder + "srst2.yaml"
    params:
        plasmids_dir = directory (results_per_sample_dir + "{sample}/cdiff_plasmids/"),
        srst2_raw_summary = results_per_sample_dir + "{sample}/cdiff_plasmids/srst2__genes__plasmids_clustered__results.txt",
        srst2_raw_summary_full = results_per_sample_dir + "{sample}/cdiff_plasmids/srst2__fullgenes__plasmids_clustered__results.txt",
    shell:
        "srst2 --input_pe {input.r1} {input.r2} --output {params.plasmids_dir}srst2 --log --gene_db {plasmids_db_dir}{plasmids_db_name} --threads {threads} --min_coverage {plasmids_coverage} --stop_after 1000000 --report_new_consensus --report_all_consensus "
        "&& if [[ -f {params.srst2_raw_summary} ]]; then  mv {params.srst2_raw_summary} {output.srst2_summary}; else touch {output.srst2_summary}; fi " #incase of negative results create empty file
        "&& if [[ -f {params.srst2_raw_summary_full} ]]; then mv {params.srst2_raw_summary_full} {output.srst2_fullsummary}; else touch {output.srst2_fullsummary}; fi"

"""
phages
"""
rule phages_search_with_srst2:
    input:
        r1 = raw_data_dir + "{sample}_1.fastq.gz",
        r2 = raw_data_dir + "{sample}_2.fastq.gz",
    output:
        srst2_summary = results_per_sample_dir + "{sample}/cdiff_phages/phages.tab",
        srst2_fullsummary = results_per_sample_dir + "{sample}/cdiff_phages/phages_full.tab",
        srst2_raw_pileup = temp(results_per_sample_dir + "{sample}/cdiff_phages/srst2__{sample}.phages_clustered.pileup"),
        srst2_raw_bam = temp(results_per_sample_dir + "{sample}/cdiff_phages/srst2__{sample}.phages_clustered.sorted.bam"),
    threads: 16
    conda:
        envs_folder + "srst2.yaml"
    params:
        phages_dir = directory (results_per_sample_dir + "{sample}/cdiff_phages/"),
        srst2_raw_summary = results_per_sample_dir + "{sample}/cdiff_phages/srst2__genes__phages_clustered__results.txt",
        srst2_raw_summary_full = results_per_sample_dir + "{sample}/cdiff_phages/srst2__fullgenes__phages_clustered__results.txt",
    shell:
        "srst2 --input_pe {input.r1} {input.r2} --output {params.phages_dir}srst2 --log --gene_db {phages_db_dir}{phages_db_name} --threads {threads} --min_coverage {phages_coverage} --stop_after 1000000 --report_new_consensus --report_all_consensus "
        "&& if [[ -f {params.srst2_raw_summary} ]]; then  mv {params.srst2_raw_summary} {output.srst2_summary}; else touch {output.srst2_summary}; fi " #incase of negative results create empty file
        "&& if [[ -f {params.srst2_raw_summary_full} ]]; then mv {params.srst2_raw_summary_full} {output.srst2_fullsummary}; else touch {output.srst2_fullsummary}; fi"

"""
resistance mutations
"""
rule resistance_mutations_positions_calling:
    input:
        r1 = raw_data_dir + "{sample}_1.fastq.gz",
        r2 = raw_data_dir + "{sample}_2.fastq.gz",
    output:
        aminoacid_change_summary = results_per_sample_dir + "{sample}/cdiff_resistance_mutations/aminoacid_at_mutated_positions.tab",
        nucleotide_change_summary = results_per_sample_dir + "{sample}/cdiff_resistance_mutations/nucleotide_at_mutated_positions.tab",
        sam = temp(results_per_sample_dir + "{sample}/cdiff_resistance_mutations/{sample}.sam"),
        unsortedbam = temp(results_per_sample_dir + "{sample}/cdiff_resistance_mutations/{sample}.unsorted.bam"),
        sortedbam = temp(results_per_sample_dir + "{sample}/cdiff_resistance_mutations/{sample}.sorted.bam"),
        vcf = temp(results_per_sample_dir + "{sample}/cdiff_resistance_mutations/{sample}.vcf.gz"),
        consensus = temp(results_per_sample_dir + "{sample}/cdiff_resistance_mutations/{sample}.consensus.fa"),
        consensus_aa = temp(results_per_sample_dir + "{sample}/cdiff_resistance_mutations/{sample}.consensus.faa"),
    threads: 16
    conda:
        envs_folder + "res_mutations.yaml" #samtools 1.9, bcftools 1.9, seqkit v0.16.1, bowtie2 2.3.0, tabix 1.4-6-g10081c4
    params:
        resistance_dir = directory (results_per_sample_dir + "{sample}/cdiff_resistance_mutations/"),
        REF= res_mutations_db_dir + res_mutations_db_name,
        aminoacid_positions= res_mutations_db_dir + res_mutations_aminoacid,
        nucleotide_positions= res_mutations_db_dir + res_mutations_nucleotide
    shell:
        #make a Bowtie index file  && #index with samtools faidx
        "bowtie2-build --quiet {params.REF} {params.REF}"
        " && samtools faidx {params.REF} "
        #align the short reads with mode very-sensitive-local
        " && mkdir -p {params.resistance_dir}"
        " && bowtie2 --quiet -1 {input.r1} -2 {input.r2} -S {output.sam} -q --very-sensitive-local --no-unal -a -x {params.REF} --threads {threads} -u 1000000"
        #make the BAM file
        " && samtools view -@ {threads} -b -O BAM -o {output.unsortedbam} -q 1 -S {output.sam}"
        " && samtools sort -@ {threads} -O BAM {output.unsortedbam} -o {output.sortedbam}"
        #produce a vcf file from the BAM file & index it
        #https://www.biostars.org/p/367960/
        " && bcftools mpileup -Ou -f {params.REF} {output.sortedbam} | bcftools call -Ou --variants-only --multiallelic-caller | bcftools norm -f {params.REF} -Oz -o {output.vcf} "
        " && tabix  {output.vcf}"
        #create a consensus sequence (ref seq seeded with sample variants)
        " && bcftools consensus --fasta-ref {params.REF} --haplotype 1 {output.vcf} > {output.consensus}"
        #get the aminoacid at the specific positions.
        #seqkit grep  -n -p gyrA consensus.fa | seqkit subseq -r 127:129 | seqkit translate --seq-type dna --transl-table 11
        " && seqkit translate --seq-type dna --transl-table 11 {output.consensus}  --out-file {output.consensus_aa}"
        " && seqkit subseq --bed {params.aminoacid_positions} {output.consensus_aa} | seqkit fx2tab > {output.aminoacid_change_summary}"
        " && bcftools mpileup -T {params.nucleotide_positions} -Ou -f {params.REF}  {output.sortedbam} | bcftools call -Ou -vm | bcftools norm -f {params.REF} -Oz | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' --print-header -o {output.nucleotide_change_summary}"
#===================================================
