configfile: 'config.yaml'

IDS = ['WT', 'LMA', 'LME', 'LMH', 'dHsdS']

rule all:
  input:
    expand("unicycler_out/{sample}", sample=IDS),    
    expand("pbmm2_out/{sample}.bam", sample=IDS),
    expand("ipdsummary_out/{sample}_modifications.gff", sample=IDS),
    expand("multimotifmaker_out/{sample}_motifs.csv", sample=IDS),

rule unicycler:
  input:
    fw = "illumina_fastq/{sample}_R1.fastq.gz",
    rv = "illumina_fastq/{sample}_R2.fastq.gz",
    pacbio = "pacbio_fa/{sample}.fa"
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
    unicycler -1 {input.fw} -2 {input.rv} -l {input.pacbio} -o {output} --threads {threads} {params.general} 2>&1>{log}
    """

rule pbmm2:
  input:
    unicycler_dir = "unicycler_out/{sample}",
    pacbio = "pacbio_fa/{sample}.fa"
  output:
    "pbmm2_out/{sample}.bam"
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
    bam = "pbmm2_out/{sample}.bam",
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
