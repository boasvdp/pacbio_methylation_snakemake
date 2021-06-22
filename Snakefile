configfile: 'config.yaml'

IDS = ['WT', 'LMA', 'LME', 'LMH', 'dHsdS']

rule all:
  input:
    expand("unicycler_out/{sample}", sample=IDS),    
    expand("pbmm2_out/{sample}_pbmm2.bam", sample=IDS),
    expand("ipdsummary_out/{sample}_modifications.gff", sample=IDS),
    expand("multimotifmaker_out/{sample}_motifs.csv", sample=IDS),

#rule flye:
#  input:
#    nanopore = "nanopore_fastq/{sample}.fastq.gz",
#    pacbio = "pacbio_fa/{sample}.fa"
#  output:
#    directory("flye_out/{sample}")
#  conda:
#    "envs/flye.yaml"
#  log:
#    "logs/flye/{sample}.log"
#  params:
#    genome_size = config['flye']['genome_size']
#  threads: 16
#  shell:
#    """
#    flye -o {output} --pacbio-raw {input.pacbio} --nano-raw {input.nanopore} --genome-size {params.genome_size} --threads {threads}
#    """

rule unicycler:
  input:
    nanopore = "nanopore_fastq/{sample}.fastq.gz",
    pacbio = "pacbio_fa/{sample}.fa",
    fw = "illumina_fastq/{sample}_R1.fastq.gz",
    rv = "illumina_fastq/{sample}_R2.fastq.gz",
    startgenes = "startgenes.fasta"
  output:
    directory("unicycler_out/{sample}")
  conda:
    "envs/unicycler.yaml"
  log:
    "logs/unicycler/{sample}.log"
  params:
    general = config['unicycler']['general']
  threads: 16
  shell:
    """
    unicycler --long {input.pacbio} --long {input.nanopore} -1 {input.fw} -2 {input.rv} -o {output} --threads {threads} --start_genes {input.startgenes} {params.general} 2>&1>{log}
    """

rule pbmm2:
  input:
    unicycler_dir = "unicycler_out/{sample}",
    pacbio = "pacbio_bam/{sample}.bam"
  output:
    "pbmm2_out/{sample}_pbmm2.bam"
  conda:
    "envs/pbmm2.yaml"
  log:
    "logs/pbmm2/{sample}.log"
  params:
    general = config['pbmm2']['general']
  threads: 16
  shell:
    """
    pbmm2 align {params.general} -j {threads} {input.unicycler_dir}/assembly.fasta {input.pacbio} {output} 2>&1>{log}
    """

rule ipdsummary:
  input:
    bam = "pbmm2_out/{sample}_pbmm2.bam",
    unicycler_dir = "unicycler_out/{sample}"
  output:
    gff = "ipdsummary_out/{sample}_modifications.gff",
    csv = "ipdsummary_out/{sample}_kinetics.csv"
  log:
    "logs/ipdsummary/{sample}.log"
  params:
    general = config['ipdsummary']['general'],
    mods = config['ipdsummary']['mods']
  threads: 16
  shell:
    """
    samtools faidx {input.unicycler_dir}/assembly.fasta
    smrtlink/install/smrtlink-release_10.1.0.119588/bundles/smrttools/install/smrttools-release_10.1.0.119588/private/pacbio/python3pkgs/kineticstools-py3/binwrap/ipdSummary \
 --reference {input.unicycler_dir}/assembly.fasta --identify {params.mods} --gff {output.gff} --csv {output.csv} -j {threads} {params.general} {input.bam} 2>&1>{log}
    """

rule multimotifmaker:
  input:
    unicycler_dir = "unicycler_out/{sample}",
    gff = "ipdsummary_out/{sample}_modifications.gff",
  output:
    "multimotifmaker_out/{sample}_motifs.csv"
  log:
    "logs/multimotifmaker/{sample}.log"
  threads: 16
  shell:
    """
    java -jar ./tools/MultiMotifMaker.jar find --fasta {input.unicycler_dir}/assembly.fasta --gff {input.gff} --output {output} --layer 1 --thread {threads} 2>&1>{log}
    """
