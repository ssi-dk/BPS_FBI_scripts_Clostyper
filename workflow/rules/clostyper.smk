"""----------------------------
Friedrich-Loeffler-Institut (https://www.fli.de/), IBIZ
date: April, 06, 2022
Author: Mostafa Abdel-Glil (mostafa.abdel-glil@fli.de)
-------------------------------


# TODO: URGENT
* change `fa` with `seqkit`
* the Rule fastp -> read the size of the reference genome to report depth correctly
* produce the interactive tree with shiptv
* Abricate summary figure in the Rmarkdown will fail if multiple blast hits exist for the same ref gene
#others
* isolates.txt, use python code to prepare it nicely
* make the log files more nice and clear
*add seqkit to the snippy conda #!/usr/bin/env python
* add gotree to the fasttree env


Assembly
* assembler_options =config["assembler_options"] ---snake variable and rule
* -opts {params.assembler} options does not wotŕk in the bash --- bash
 * shovill --outdir $OUTDIR --cpus $THREADS --R1 $READ1 --R2 $READ2 --force --ram 10 $assembler_options $assembler_options_additional

Prokka
* --gram pos can not be used, no `signalp` in conda
"""
import os
import tempfile
#configfile: "config.yaml"
#defineFolders
#working_dir=config['working_dir']
#snakemake_folder=config['snakemake_folder']
#raw_data_dir=config['raw_data_dir']
#fasta_dir=config['fasta_dir']
#results_dir=config['results_dir']
raw_data_dir= working_dir + "1_input_samples/"
fasta_dir= working_dir + "1_input_samples/"
#results_dir= working_dir + "2_pipeline_results/"
results_per_sample_dir= working_dir + "2_results_per_sample/"
results_all_samples_dir= working_dir + "3_results_all_samples/"
assembly_dir_all= results_all_samples_dir + "assembly_fasta/"
reports_dir= working_dir + "4_all_reports/"
temporary_todelete= results_all_samples_dir + "temporary_todelete/"

#samples
SAMPLES, = glob_wildcards( raw_data_dir + "{sample}_1.fastq.gz")
GENOMES, = glob_wildcards( fasta_dir + "{genome}.fasta")
SAMPLES = set(SAMPLES)
GENOMES = set(GENOMES)
DATASET =  [SAMPLES, GENOMES]
#scripts_paths
bin_dir= snakemake_folder + "../scripts/",
#lib_dir= snakemake_folder + "lib/"

#nullarbor_bin= snakemake_folder + "lib/nullarbor/bin/"
envs_folder= snakemake_folder + "../envs/"
#envs_folder="lib/envs/"
#variables
if "AMR_db_abricate" in config:
    AMR_db_abricate= config["AMR_db_abricate"]
else:
    AMR_db_abricate="ncbi"
if "VF_db_abricate" in config:
    VF_db_abricate= config["VF_db_abricate"]
else:
    VF_db_abricate="vfdb"
if "prokka_genustag" in config:
    prokka_genustag= config["prokka_genustag"]
else:
    prokka_genustag= ""
if "prokka_params" in config:
    prokka_params= config["prokka_params"]
else:
    prokka_params= ""
#mlst_parameters
if "MLST_options" in config:
    MLST_options= config["MLST_options"]
else:
    MLST_options= ""
#abricate_parameters
if "abricate_vir_options" in config:
    abricate_vir_options= config["abricate_vir_options"]
else:
    abricate_vir_options= ""
if "abricate_res_options" in config:
    abricate_res_options= config["abricate_res_options"]
else:
    abricate_res_options= ""
if "snippy_options" in config:
    snippy_options= config["snippy_options"]
else:
    snippy_options= ""
if "assembler_options" in config:
    assembler_options= config["assembler_options"]
else:
    assembler_options= ""
if "kraken2_opts" in config:
    kraken2_opts= config["kraken2_opts"]
else:
    kraken2_opts= ""
if "shovill_assembler" in config:
    shovill_assembler= config["shovill_assembler"]
else:
    shovill_assembler= "spades"

if "shovill_options" in config:
    shovill_options= config["shovill_options"]
else:
    shovill_options= ""

#assembler_options= ""
taxoner= config["taxoner"]
reference= config["reference"]

ASSEMBLY=[expand( assembly_dir_all + "{sample}.fasta", sample=SAMPLES), expand( assembly_dir_all + "{genome}.fasta", genome=GENOMES)],
KRAKEN= expand(results_per_sample_dir + "{sample}/kraken.reads.tab", sample=SAMPLES),
KRAKEN2ASSEMBLIES=[ expand(results_per_sample_dir + "{sample}/kraken.contigs.tab", sample=SAMPLES), expand(results_per_sample_dir + "{genome}/kraken.contigs.tab", genome=GENOMES) ],
FASTASTATS = results_all_samples_dir + "denovo.tab",
REFERENCE=results_all_samples_dir + "ref.fa"
SNPANALYSIS=[expand( results_per_sample_dir + "{sample}/snippy/snps.tab", sample=SAMPLES),
expand( results_per_sample_dir + "{genome}/snippy/snps.tab", genome=GENOMES), results_all_samples_dir  + "snippycore_results/core.aln", results_all_samples_dir + "snippycore_results/distances.tab"]
SAMPLENAMES= results_all_samples_dir + "isolates.txt"
FASTP= expand( results_per_sample_dir + "{sample}/fastp.html", sample=SAMPLES)
VIRULENCE= results_all_samples_dir + "virulence.vfdb.tab"
RESISTANCE =[results_all_samples_dir + "resistance.alldatabases.tab",
results_all_samples_dir + "resistance.ncbi.tab", results_all_samples_dir + "resistance.card.tab", results_all_samples_dir + "resistance.resfinder.tab",  results_all_samples_dir + "resistance.argannot.tab", results_all_samples_dir + "resistance.megares.tab"]
PROKKA= [expand(results_per_sample_dir + "{sample}/prokka/{sample}.gff", sample=SAMPLES),
expand( results_per_sample_dir + "{genome}/prokka/{genome}.gff", genome=GENOMES)]
PHYLOGENY= [results_all_samples_dir + "core.newick", results_all_samples_dir + "core-midpoint.newick"]
PANGENOME= [results_all_samples_dir + "roary_results/pan_genome_reference.fa",
results_all_samples_dir + "roary_results/pangenome_frequency.png", results_all_samples_dir + "roary_results/pan.svg"]
MLST = results_all_samples_dir + "mlst.tab"
REPORT = reports_dir + "clostyper_report.html"

#master rule
#rule all:
#    input: KRAKEN2ASSEMBLIES,# ASSEMBLY, FASTASTATS, REFERENCE, SNPANALYSIS, SAMPLENAMES, FASTP, VIRULENCE, RESISTANCE, PROKKA, PHYLOGENY, PANGENOME, MLST, REPORT


"""
collect assemblies
"""
rule collect_assemblies:
    input:
        assemblies = fasta_dir + "{genome}.fasta",
#        yield_na = lib_dir + "yield.tab"
    output:
        contig = assembly_dir_all + "{genome}.fasta",
#        contig_yield = results_dir + "{genome}/{genome}.yield.tab",
    conda:
        envs_folder + "bioawk.yaml" #spades, sickle
    shell:
        "bash {bin_dir}fa_rename.sh {input.assemblies} {output.contig}"
        #" && cp -v -f {input.yield_na} {output.contig_yield}" #ln -s -f ##* Prokka does not like contigs ID > 37, '--centre X --compliant' is not the way always
"""
assembly
"""
rule assembly:
    input:
        r1 = raw_data_dir + "{sample}_1.fastq.gz",
        r2 = raw_data_dir + "{sample}_2.fastq.gz",
    output:
#        contig = results_dir + "{sample}/contigs.fa",
        contigs_fasta = assembly_dir_all + "{sample}.fasta",
    threads: 16 #increasing threads produces errors with spades
    conda:
        envs_folder + "shovill.yaml" #spades, sickle
    params:
        #shovill_assembler = shovill_assembler,
        assembly_dir = directory (results_per_sample_dir + "{sample}")
    shell:
        "shovill --assembler {shovill_assembler} --R1 {input.r1} --R2 {input.r2} --outdir {params.assembly_dir}/shovill  --cpus {threads} --force {shovill_options}"
        "&& cp {params.assembly_dir}/shovill/contigs.fa {output.contigs_fasta} "
        #"&& ln -sf ./shovill/contigs.fa {output.contig}" #-t trimmomatic -f true
"""
TODO kraken
"""
rule kraken:
    input:
        r1 = raw_data_dir + "{sample}_1.fastq.gz",
        r2 = raw_data_dir + "{sample}_2.fastq.gz",
        db=config["kraken"],
        #taxoner= config["taxoner"]
    output:
        kraken_tab= results_per_sample_dir + "{sample}/kraken.reads.tab"
    threads: 32
    conda:
      envs_folder + "kraken2.yaml"
    shell:
        "kraken2 --db {input.db} --threads {threads} --gzip-compressed  {input.r1} {input.r2} --output -  --report {output.kraken_tab} {kraken2_opts} 2>&1 | sed 's/^/[kraken] /'"

rule kraken_assemblies:
    input:
        assemblies = assembly_dir_all + "{sample}.fasta",
        db=config["kraken"],
    output:
        kraken_tab= results_per_sample_dir + "{sample}/kraken.contigs.tab",
    threads: 32
    conda:
      envs_folder + "kraken2.yaml"
    shell:
        "kraken2 --db {input.db} --threads {threads}  {input.assemblies} --output -  --report {output.kraken_tab} {kraken2_opts} 2>&1 | sed 's/^/[kraken] /'"

"""
denovo fasta
"""
rule fasta_denovo:
    input:
        fasta = expand( assembly_dir_all + "{sample}.fasta", sample=SAMPLES),
        fasta_assembly = expand( assembly_dir_all + "{sample}.fasta", sample=GENOMES),
    output:
        denovo = results_all_samples_dir + "denovo.tab",
    conda:
        envs_folder + "perl.yaml"
    shell:
        """var={results_all_samples_dir} && {bin_dir}fa --minsize 500 -e -t {input.fasta} {input.fasta_assembly} | sed "s#$var##g" | tee -a {output.denovo}"""
"""
reference
"""
rule reference:
    output:
        ref_fasta= results_all_samples_dir + "ref.fa",
    conda:
        envs_folder + "reference.yaml"
    shell:
        "any2fasta -q -u {reference} > {output.ref_fasta} && \
        samtools faidx {output.ref_fasta}"
"""
isolates.txt #revision
"""
names =expand( "{sample}", sample=SAMPLES),
names_genomes =expand( "{sample}", sample=GENOMES),
isolates = [names, names_genomes]
rule isolates_list:
    output:
        isolates_list= results_all_samples_dir + "isolates.txt",
        tempfile=temp(results_all_samples_dir + "tempfile",)
    shell:
        '''
        echo -ne "{isolates}" > {output.tempfile} && bash {bin_dir}fixnames.sh {output.tempfile} > {output.isolates_list}
        '''
#import sys
#sys.stdout = open('results/isolates.txt', 'w')
#for w in isolates:
#print(w, sep="\n")
#print (*names, sep='\n')
#print (*names_genomes, sep='\n')

"""
fastp
"""
rule fastp:
    input:
        r1 = raw_data_dir + "{sample}_1.fastq.gz",
        r2 = raw_data_dir + "{sample}_2.fastq.gz",
        ref_fasta= results_all_samples_dir + "ref.fa",
    output:
        fastp_html = results_per_sample_dir + "{sample}/fastp.html",
        fastp_json = results_per_sample_dir + "{sample}/fastp.json",
        #pandoc_md = temp(results_dir + "{sample}/{sample}.pandoc.md"),
        #before_filtr_stats= temp(results_dir + "{sample}/{sample}.stats.txt"),
        #fastp_stats= results_dir + "{sample}/{sample}.fastp.tab",
    threads: 16
    conda:
      envs_folder + "fastp.yaml" #fastp, pandoc
    params:
      report_title="{sample}",
    shell:
        "fastp -i {input.r1}  -I {input.r2} --disable_quality_filtering --disable_adapter_trimming --disable_length_filtering --disable_trim_poly_g --json {output.fastp_json} --html {output.fastp_html} --thread {threads} --verbose --report_title {params.report_title} 2>&1 | sed 's/^/[fastp] /' "
        #" && pandoc -s -r html {output.fastp_html} -o {output.pandoc_md}"
        #" && grep -m 1 -A 10 'Before filtering' {output.pandoc_md} | tail -n5 > {output.before_filtr_stats}"
        #" && sh {bin_dir}/fastp_stats.sh {output.before_filtr_stats} > {output.fastp_stats}"

"""
virulence
"""
rule virulence:
    input:
        contig = assembly_dir_all + "{sample}.fasta",
    output:
        virulence = results_per_sample_dir + "{sample}/virulence.vfdb.tab",
    threads: 8
    conda:
        envs_folder + "abricate.yaml"
    params:
        VF_db_abricate = VF_db_abricate,
    shell:
        "abricate --nopath --db {params.VF_db_abricate} --threads {threads} {abricate_vir_options} {input.contig} > {output.virulence}"

"""
resistance
"""
rule resistance:
    input:
        contig = assembly_dir_all + "{sample}.fasta",
    output:
        resistance_all = results_per_sample_dir + "{sample}/resistance.alldatabases.tab",
        resistance_ncbi = results_per_sample_dir + "{sample}/resistance.ncbi.tab",
        resistance_card = results_per_sample_dir + "{sample}/resistance.card.tab",
        resistance_resfinder = results_per_sample_dir + "{sample}/resistance.resfinder.tab",
        resistance_argannot = results_per_sample_dir + "{sample}/resistance.argannot.tab",
        resistance_megares = results_per_sample_dir + "{sample}/resistance.megares.tab",
    threads: 8
    conda:
        envs_folder + "abricate.yaml"
    params:
        AMR_db_abricate = AMR_db_abricate,
    shell:
        "abricate --nopath --db {params.AMR_db_abricate} --threads {threads} {abricate_res_options} {input.contig} > {output.resistance_all}"
        "&& abricate --nopath --db card --threads {threads} {abricate_res_options} {input.contig} > {output.resistance_card}"
        "&& abricate --nopath --db ncbi --threads {threads} {abricate_res_options} {input.contig} > {output.resistance_ncbi}"
        "&& abricate --nopath --db resfinder --threads {threads} {abricate_res_options} {input.contig} > {output.resistance_resfinder}"
        "&& abricate --nopath --db argannot --threads {threads} {abricate_res_options} {input.contig} > {output.resistance_argannot}"
        "&& abricate --nopath --db megares --threads {threads} {abricate_res_options} {input.contig} > {output.resistance_megares}"

"""
resistance and virulence summary
"""

rule summary_res_vir: #run roary
    input:
        #bin_dir=config["bin_dir"],
        resistance_all= expand(results_per_sample_dir + "{sample}/resistance.alldatabases.tab",sample=SAMPLES),
        resistance_ncbi = expand(results_per_sample_dir + "{sample}/resistance.ncbi.tab",sample=SAMPLES),
        resistance_card = expand(results_per_sample_dir + "{sample}/resistance.card.tab",sample=SAMPLES),
        resistance_resfinder = expand(results_per_sample_dir + "{sample}/resistance.resfinder.tab",sample=SAMPLES),
        resistance_argannot = expand(results_per_sample_dir + "{sample}/resistance.argannot.tab",sample=SAMPLES),
        resistance_megares = expand(results_per_sample_dir + "{sample}/resistance.megares.tab", sample=SAMPLES),
        virulence = expand(results_per_sample_dir + "{sample}/virulence.vfdb.tab", sample=SAMPLES),
        resistance_all_genome= expand(results_per_sample_dir + "{genome}/resistance.alldatabases.tab",genome=GENOMES),
        resistance_ncbi_genome = expand(results_per_sample_dir + "{genome}/resistance.ncbi.tab",genome=GENOMES),
        resistance_card_genome = expand(results_per_sample_dir + "{genome}/resistance.card.tab",genome=GENOMES),
        resistance_resfinder_genome = expand(results_per_sample_dir + "{genome}/resistance.resfinder.tab",genome=GENOMES),
        resistance_argannot_genome = expand(results_per_sample_dir + "{genome}/resistance.argannot.tab",genome=GENOMES),
        resistance_megares_genome = expand(results_per_sample_dir + "{genome}/resistance.megares.tab", genome=GENOMES),
        virulence_genome = expand(results_per_sample_dir + "{genome}/virulence.vfdb.tab", genome=GENOMES),
    output:
        sum_resistance_all = results_all_samples_dir + "resistance.alldatabases.tab",
        sum_resistance_ncbi = results_all_samples_dir + "resistance.ncbi.tab",
        sum_resistance_card = results_all_samples_dir + "resistance.card.tab",
        sum_resistance_resfinder = results_all_samples_dir + "resistance.resfinder.tab",
        sum_resistance_argannot = results_all_samples_dir + "resistance.argannot.tab",
        sum_resistance_megares = results_all_samples_dir + "resistance.megares.tab",
        sum_virulence = results_all_samples_dir + "virulence.vfdb.tab",
    conda:
        envs_folder + "abricate.yaml"
    params:
        AMR_db_abricate = AMR_db_abricate,
    shell:
        "abricate --summary {input.resistance_all} {input.resistance_all_genome} > {output.sum_resistance_all}"
        "&& abricate --summary {input.resistance_ncbi} {input.resistance_ncbi_genome} > {output.sum_resistance_ncbi}"
        "&& abricate --summary {input.resistance_card} {input.resistance_card_genome} > {output.sum_resistance_card}"
        "&& abricate --summary {input.resistance_resfinder} {input.resistance_resfinder_genome} > {output.sum_resistance_resfinder}"
        "&& abricate --summary {input.resistance_argannot} {input.resistance_argannot_genome} > {output.sum_resistance_argannot}"
        "&& abricate --summary {input.resistance_megares} {input.resistance_megares_genome} > {output.sum_resistance_megares}"
        "&& abricate --summary {input.virulence} {input.virulence_genome} > {output.sum_virulence}"
"""
Prokka
"""
rule prokka:
    input:
        contig= assembly_dir_all + "{sample}.fasta",
    output:
        prokka_gff= results_per_sample_dir + "{sample}/prokka/{sample}.gff",
        prokka_gbk= results_per_sample_dir + "{sample}/prokka/{sample}.gbk",
    threads: 32
    conda: envs_folder + "prokka_1.14.5.yaml"
    params:
        prokka_outdir= results_per_sample_dir + "{sample}/prokka",
        prefix="{sample}",
        strain_name="{sample}",
        locustag= "{sample}",
        sqn= results_per_sample_dir + "{sample}/prokka/{sample}.sqn",
        fna= results_per_sample_dir + "{sample}/prokka/{sample}.fna",
        fsa= results_per_sample_dir + "{sample}/prokka/{sample}.fsa",
        tbl= results_per_sample_dir + "{sample}/prokka/{sample}.tbl",
    shell:
        """prokka --cpus {threads} --force --kingdom Bacteria --outdir {params.prokka_outdir}/ --prefix {params.prefix} --genus {prokka_genustag} --gcode 11 --locustag {params.locustag} --strain {params.strain_name} {prokka_params} {input.contig} 2>&1 | sed 's/^/[prokka] /' && rm -f {params.sqn} {params.fsa} {params.fna} {params.tbl}""" #--centre Institute --compliant to generate clean contig names
rule prokka_ln:
    input:
        prokka_gff= results_per_sample_dir + "{sample}/prokka/{sample}.gff",
        prokka_gbk= results_per_sample_dir + "{sample}/prokka/{sample}.gbk",
    output:
        prokka_gff= results_per_sample_dir + "{sample}/contigs.gff",
        prokka_gbk= results_per_sample_dir + "{sample}/contigs.gbk",
    params:
        prokka_outdir= results_per_sample_dir + "{sample}/prokka",
    shell:
        "ln -s -f $(realpath {input.prokka_gff}) {output.prokka_gff}; ln -s -f $(realpath {input.prokka_gbk}) {output.prokka_gbk}"
"""
snippy
"""
#export LD_LIBRARY_PATH=/usr/local/lib:/usr/lib:/usr/local/lib64:/usr/lib64
rule snippy:
    input:
        r1 = raw_data_dir + "{sample}_1.fastq.gz",
        r2 = raw_data_dir + "{sample}_2.fastq.gz",
    output:
        snippy_snps= results_per_sample_dir + "{sample}/snippy/snps.tab",
        #snippy_folder=  results_dir + "{sample}",
    threads: 32
    conda: envs_folder + "snippy.yaml"
    params:
        snippy_outdir= results_per_sample_dir + "{sample}/snippy",
        #snippy_outdir_cp= results_dir + "{sample}/",
    shell:
        "snippy --force  --cpus {threads} --ram {threads} --outdir {params.snippy_outdir} --ref {reference} --R1 {input.r1} --R2 {input.r2} {snippy_options} 2>&1 | sed 's/^/[snippy] /'"
rule snippy_ln:
    input:
        snippy_snps= results_per_sample_dir + "{sample}/snippy/snps.tab",
    output:
        snippy_snps_tab=  temp(results_per_sample_dir + "{sample}/snps.tab"),
        snippy_snps_aligned= temp(results_per_sample_dir + "{sample}/snps.aligned.fa"),
        #snippy_snps_rawvcf= results_dir + "{sample}/snps.raw.vcf",
        snippy_snps_vcf= temp(results_per_sample_dir + "{sample}/snps.vcf"),
        #snippy_snps_bam= results_dir + "{sample}/snps.bam",
        #snippy_snps_bai= results_dir + "{sample}/snps.bam.bai",
        #snippy_outdir_ref= temp(directory(results_dir + "{sample}/snippy/reference/")),
        snippy_snps_log= temp(results_per_sample_dir + "{sample}/snps.log"),
    conda: envs_folder + "snippy.yaml"
    params:
        snippy_outdir= results_per_sample_dir + "{sample}/snippy",
        #snippy_outdir_ref= temp(results_dir + "{sample}/snippy/reference"),
    shell:
        "ln -s -f $(realpath {params.snippy_outdir}/snps.tab) {output.snippy_snps_tab}; ln -s -f $(realpath {params.snippy_outdir}/snps.aligned.fa) {output.snippy_snps_aligned} ;  ln -s -f $(realpath {params.snippy_outdir}/snps.vcf) {output.snippy_snps_vcf};  ln -s -f $(realpath {params.snippy_outdir}/snps.log) {output.snippy_snps_log} && rm -f {params.snippy_outdir}/reference/{{ref.fa,ref.gff,ref.txt}} " #ln -s -f $(realpath {params.snippy_outdir}/snps.bam) {output.snippy_snps_bam}; ln -s -f $(realpath {params.snippy_outdir}/snps.bam.bai) {output.snippy_snps_bai}; ln -s -f $(realpath {params.snippy_outdir}/snps.raw.vcf) {output.snippy_snps_rawvcf};
rule snippy_assemblies:
    input:
        contig = fasta_dir + "{genome}.fasta"
    output:
        snippy_snps= results_per_sample_dir + "{genome}/snippy/snps.tab",
    threads: 32
    conda: envs_folder + "snippy.yaml"
    params:
        snippy_outdir= results_per_sample_dir + "{genome}/snippy",
    shell:
        "snippy --force  --cpus {threads} --ram {threads} --outdir {params.snippy_outdir} --ref {reference} --ctgs {input.contig} {snippy_options} 2>&1 | sed 's/^/[snippy] /'"

def snippy_folders(wildcards):
    return expand(results_per_sample_dir + "{sample}", sample=SAMPLES)
def snippy_folders_assemblies(wildcards):
    return expand(results_per_sample_dir + "{genome}", genome=GENOMES)

rule snippy_core:
    input:
        #snippy_folders=  expand( results_dir + "{sample}", sample=SAMPLES),
        #snippy_folders_assemblies=  expand( results_dir + "{genome}", genome=GENOMES),
        a=expand( results_per_sample_dir + "{sample}/snps.tab", sample=SAMPLES),
        b=expand( results_per_sample_dir + "{sample}/snps.aligned.fa", sample=SAMPLES),
        #c=expand( results_dir + "{sample}/snps.raw.vcf", sample=SAMPLES),
        d=expand( results_per_sample_dir + "{sample}/snps.vcf", sample=SAMPLES),
        #e=expand( results_dir + "{sample}/snps.bam", sample=SAMPLES),
        #f=expand( results_dir + "{sample}/snps.bam.bai", sample=SAMPLES),
        g=expand( results_per_sample_dir + "{sample}/snps.log", sample=SAMPLES),
        h=expand( results_per_sample_dir + "{genome}/snps.tab", genome=GENOMES),
        i=expand( results_per_sample_dir + "{genome}/snps.aligned.fa", genome=GENOMES),
        #j=expand( results_dir + "{genome}/snps.raw.vcf", genome=GENOMES),
        k=expand( results_per_sample_dir + "{genome}/snps.vcf", genome=GENOMES),
        #l=expand( results_dir + "{genome}/snps.bam", genome=GENOMES),
        #m=expand( results_dir + "{genome}/snps.bam.bai", genome=GENOMES),
        n=expand( results_per_sample_dir + "{genome}/snps.log", genome=GENOMES),
    output:
        snippycorenoRef= results_all_samples_dir  + "snippycore_results/core-noRef.aln",
        snippycorenoRefsnpsites= results_all_samples_dir  + "snippycore_results/core-noRef-snpsites.aln",
        snippycore= results_all_samples_dir  + "snippycore_results/core.aln",
        snippycore_text= results_all_samples_dir  + "snippycore_results/core.txt",
    conda: envs_folder + "snippy.yaml"
    params:
        snippycoreoutdir= results_all_samples_dir + "snippycore_results/core",
        snippy_folders = snippy_folders,
        snippy_folders_assemblies = snippy_folders_assemblies,
        #snippy_outdir_cp= results_dir + "{sample}/",
    shell:
        """
        var={{results_per_samples_dir}}Reference && snippy-core  --ref {reference} {params.snippy_folders} $(echo {params.snippy_folders_assemblies} | sed "s#$var##g") --prefix {params.snippycoreoutdir}  2>&1 | sed 's/^/[snippy-core] /' &&
        seqkit grep --quiet -v -p Reference {output.snippycore} -o {output.snippycorenoRef}  &&
        snp-sites -o {output.snippycorenoRefsnpsites} {output.snippycorenoRef}
        """
"""
Core genome phylogeny using FastTree
"""
rule FastTree: #run fasttree
    input:
        snippycore= results_all_samples_dir  + "snippycore_results/core-noRef.aln",
    output:
        phylogeny_tree= results_all_samples_dir + "core.newick",
        phylogeny_svg= results_all_samples_dir + "core.svg"
    #threads: 64
    #benchmark: benchmarks_folder + "FastTree.txt"
    conda: envs_folder + "FastTree.yaml" #needs revision
    shell:
        "FastTree -nt -gtr {input.snippycore} > {output.phylogeny_tree} "
        "&& nw_display -S -s -w 1024 -l 'font-size:12;font-family:sans-serif;' -i 'opacity:0' -b 'opacity:0' -v 16 {output.phylogeny_tree}  > {output.phylogeny_svg}"
#seqkit grep -v Reference core.aln > Core-noRef.aln
"""
reroot tree
"""
rule gotree: #run fasttree
    input:
        phylogeny_tree= results_all_samples_dir + "core.newick",
    output:
        phylogeny_tree_reroot= results_all_samples_dir + "core-midpoint.newick",
    #threads: 64
    #benchmark: benchmarks_folder + "FastTree.txt"
    conda: envs_folder + "gotree.yaml" #needs revision
    shell:
        "gotree reroot midpoint -i {input.phylogeny_tree} -o {output.phylogeny_tree_reroot} "
#seqkit grep -v Reference core.aln > Core-noRef.aln
"""
Roary
"""
#export PERL5LIB=$PERL5LIB:/home/software/Roary/lib/
rule Roary: #run roary
    input:
        #bin_dir=config["bin_dir"],
        prokka_gff= expand(results_per_sample_dir + "{sample}/prokka/{sample}.gff", sample=SAMPLES),
        prokka_assemblies_gff= expand(results_per_sample_dir + "{genome}/prokka/{genome}.gff", genome=GENOMES),
    output:
        tmp_dir=temp(directory(temporary_todelete)),
        Roary_pangenome_fa= results_all_samples_dir + "roary_results/pan_genome_reference.fa",
        Roary_pangenome_summary= results_all_samples_dir + "roary_results/summary_statistics.txt",
        #Roary_aln= results_dir + "roary/core_gene_alignment.aln",
        roary_presenceabsence= results_all_samples_dir + "roary_results/gene_presence_absence.csv",
        roary_acc= results_all_samples_dir + "roary_results/accessory_binary_genes.fa.newick",
    threads: 64
    log: results_all_samples_dir + "roary.log"
    conda: envs_folder + "roary.yaml"
    params:
        options=config["roary_params"],
        Roary_dir=results_all_samples_dir + "roary_results"
    shell:
        "bash {bin_dir}fixRoaryOutDirError.sh {params.Roary_dir} {output.tmp_dir}"
        " && roary -p {threads} -f {params.Roary_dir} {params.options} {input.prokka_gff} {input.prokka_assemblies_gff} 2>&1 | sed 's/^/[roary] /' | tee -a {log}" #-e, core genes alignment using PRANK, -r, Rplots ,
        # sometimes you may need to set -i (minimum percentage identity for blastp) and -s (dont split paralogs), according to the organism
        #" && python3 {bin_dir}roary_plots.py {params.Roary_dir}/accessory_binary_genes.fa.newick {params.Roary_dir}/gene_presence_absence.csv"
        #" && mv -t {params.Roary_dir} pangenome_frequency.png pangenome_matrix.png pangenome_pie.png"
rule Roary_plots:
    input:
        roary_presenceabsence= results_all_samples_dir + "roary_results/gene_presence_absence.csv",
        roary_acc= results_all_samples_dir + "roary_results/accessory_binary_genes.fa.newick",
    output:
        Roary_dir=results_all_samples_dir + "roary_results/pangenome_frequency.png",
        Roary_pie=results_all_samples_dir + "roary_results/pangenome_pie.png",
    conda: envs_folder + "roary_plots.yaml"
    params:
        #options=config["roary_params"],
        Roary_dir=results_all_samples_dir + "roary_results"
    shell:
        " python3.5 {bin_dir}roary_plots.py {params.Roary_dir}/accessory_binary_genes.fa.newick {params.Roary_dir}/gene_presence_absence.csv"
        " && mv -t {params.Roary_dir} pangenome_frequency.png pangenome_matrix.png pangenome_pie.png"
#source ~/.bash_profile SVG.pm
rule Roary_svg: #run roary
    input:
        roary_pangenome_fa= results_all_samples_dir + "roary_results/pan_genome_reference.fa",
        roary_presenceabsence= results_all_samples_dir + "roary_results/gene_presence_absence.csv",
        roary_acc= results_all_samples_dir + "roary_results/accessory_binary_genes.fa.newick",
    output:
        acc_svg= results_all_samples_dir + "roary_results/acc.svg",
        pan_svg= results_all_samples_dir + "roary_results/pan.svg",
    #threads: 64
    conda: envs_folder + "perl.yaml"
    params:
        options=config["roary_params"],
        Roary_dir=results_all_samples_dir + "roary_results"
    shell:
        " nw_display -S -s -w 1024 -l 'font-size:12;font-family:sans-serif;' -i 'opacity:0' -b 'opacity:0' -v 16 {input.roary_acc} > {output.acc_svg}"
        " && {bin_dir}roary2svg.pl {input.roary_presenceabsence} > {output.pan_svg}"

"""
mlst
"""
rule MLST:
    input:
        contigs= expand(assembly_dir_all + "{sample}.fasta", sample=SAMPLES),
        contigs_assemblies= expand(assembly_dir_all + "{genome}.fasta", genome=GENOMES),
        #ref_fasta= results_dir + "ref.fa",
    output:
        mlstresults= results_all_samples_dir + "mlst.tab",
    conda:
        envs_folder + "mlst.yaml"
    shell:
        """var={results_all_samples_dir} && mlst --quiet --nopath {MLST_options} {input.contigs} {input.contigs_assemblies} | sed "s#$var##g" | tee -a {output.mlstresults}"""   #{input.ref_fasta}
"""
snpdists
"""
rule snpdists:
    input:
        snippycore= results_all_samples_dir  + "snippycore_results/core-noRef.aln",
    output:
        snpdistsresults= results_all_samples_dir  + "snippycore_results/distances.tab",
    conda:
        envs_folder + "snpdists.yaml"
    shell:
        "snp-dists -b {input.snippycore} > {output.snpdistsresults}"

"""
create report
"""
rule make_report:
    input:
        snpdistsresults= results_all_samples_dir + "snippycore_results/distances.tab",
        isolates_list= results_all_samples_dir + "isolates.txt",
        fastp_json = expand(results_per_sample_dir + "{sample}/fastp.json", sample=SAMPLES),
        kraken_tab= expand(results_per_sample_dir + "{sample}/kraken.reads.tab", sample=SAMPLES),
        kraken_tab_contigs= expand(results_per_sample_dir + "{sample}/kraken.contigs.tab", sample=SAMPLES),
        kraken_tab_assemblies= expand(results_per_sample_dir + "{genome}/kraken.contigs.tab", genome=GENOMES),
        assembledata = results_all_samples_dir + "denovo.tab",
        sum_resistance_all = results_all_samples_dir + "resistance.alldatabases.tab",
        sum_resistance_ncbi = results_all_samples_dir + "resistance.ncbi.tab",
        sum_resistance_card = results_all_samples_dir + "resistance.card.tab",
        sum_resistance_resfinder = results_all_samples_dir + "resistance.resfinder.tab",
        sum_resistance_argannot = results_all_samples_dir + "resistance.argannot.tab",
        sum_resistance_megares = results_all_samples_dir + "resistance.megares.tab",
        sum_virulence = results_all_samples_dir + "virulence.vfdb.tab",
        mlstresults= results_all_samples_dir + "mlst.tab",
        snippycore_text= results_all_samples_dir  + "snippycore_results/core.txt",
        phylogeny_tree= results_all_samples_dir + "core-midpoint.newick",
        Roary_pangenome_summary= results_all_samples_dir + "roary_results/summary_statistics.txt",
        Roary_pie=results_all_samples_dir + "roary_results/pangenome_pie.png",
        pan_svg= results_all_samples_dir + "roary_results/pan.svg",
        acc_svg= results_all_samples_dir + "roary_results/acc.svg",
        Roary_dir=results_all_samples_dir + "roary_results/pangenome_frequency.png",

        #kraken_tab_assemblies= expand(results_dir + "{genome}/{genome}.kraken.tab", genome=GENOMES), #localres=expand(kraken_fastq_dir+"{sample}.report.txt", sample=SAMPLES),#expand(kraken_fastq_dir + "{sample}.report.txt",sample=SAMPLES),
    output:
        clostyper_report = reports_dir + "clostyper_report.html",#final result as input for MQC
#    conda:
#        envs_folder + "report.yaml"
    params:
        results_dir = results_all_samples_dir,
        reports_dir = reports_dir,
    script:
        "../scripts/clostyper_report.Rmd"
    #shell:
        #"""        Rscript --vanilla -e 'rmarkdown::render("/home/mostafa.abdel/aProjects/gitProjects/ClosTyper/bin/clostyper_report.Rmd", output_file="/home/mostafa.abdel/aProjects/gitProjects/ClosTyper/123/4_all_reports/clostyper_report.html", quiet=TRUE, knit_root_dir = "/home/mostafa.abdel/aProjects/gitProjects/ClosTyper", params = list(rmd="/home/mostafa.abdel/aProjects/gitProjects/ClosTyper/#.snakemake/scripts/tmp0p75jyvg.clostyper_report.Rmd"))'
        #"""
