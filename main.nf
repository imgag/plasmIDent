#!/usr/bin/env nextflow

/* 
===============================================================================
                              P L A S M I D E N T  
===============================================================================
Nextflow pipeline for resistance plasmid identification and annotation using 
Nanopore reads and bacterial genome assemblies
-------------------------------------------------------------------------------
@ Author
Caspar Groß <post@caspar.bio>
-------------------------------------------------------------------------------
@ Documentation
https://github.com/caspargross/hybridassembly/README.md
------------------------------------------------------------------------------
*/

// Check special input parameters

process filter_reads {
// Subsample large datasets to given target coverage
// Filtering with focus on high qual to keep enough short reads for plasmids
    tag { id }

    input:
    tuple val(id), path(assembly), path(lr)

    output:
    tuple val(id), path(assembly), path('reads_filtered.fastq'), emit: filtered

    script:
    if (params.noSubsampling)
        """
        zcat -f ${lr} > reads_filtered.fastq
        """
    else
        """
        ${params.env}
        len=\$(grep -v '>' ${assembly} | wc -c)
        nbases=\$(expr \$len * ${params.mappingCov})
        filtlong -t \$nbases --length_weight 0 ${lr} > reads_filtered.fastq
        """
}

process save_plasmids {
// Save plasmids as a separate fasta file
    tag { id + ":" + contigName }
    publishDir({ "${params.outDir}/${id}/plasmids/" }, mode: 'copy')

    input:
    tuple val(id), path(lr), val(contigName), val(length), val(sequence)
   
    output:
    path("${contigName}.fasta")

    script:
    """
    echo ">${contigName} len=${length}" > ${contigName}.fasta
    echo ${sequence} >> ${contigName}.fasta

    """
}

process pad_plasmids {
// Add prefix and suffix with sequence from oppsig end to each plasmid
    tag { id + ":" + contigName }

    input: 
    tuple val(id), path(lr), val(contigName), val(length), val(sequence)

    output: 
    tuple val(id), path("${id}_${contigName}_padded.fasta"), path(lr), val(contigName), emit: padded
    
    shell:
    '''
    echo '>!{contigName}' >  !{id}_!{contigName}_padded.fasta

    echo !{sequence} | awk '{print \
        substr($1, length($1)-(!{params.seqPadding} - 1), length($1))\
        $1 \
        substr($1, 1, !{params.seqPadding})\
        }' >> !{id}_!{contigName}_padded.fasta

    '''
}

process combine_padded_contigs {
// Recombines padded contigs into a single fasta
    tag { id + ":" + contigName }
    publishDir({ "${params.outDir}/${id}/alignment/" }, mode: 'copy')

    input:
    tuple val(id), path(assembly), val(lr), val(contigName)

    output:
    tuple val(id), path("${id}_padded.fasta"), val(lr), val('padded'), emit: padded_assembly

    script:
    """
    cat \$(echo ${assembly} | tr -d '[],') > ${id}_padded.fasta 
    """
}

process map_longreads {
// Use minimap2 to align longreads to padded contigs
    publishDir({ "${params.outDir}/${id}/alignment/" }, mode: 'copy')
    tag { id }

    input:
    tuple val(id), path(assembly), path(lr), val(type)

    output:
    tuple val(id), path(assembly), val(type), path("${id}_${type}_lr.bam"), path("${id}_${type}_lr.bam.bai"), emit: bam

    script:
    """
    ${params.env}
    minimap2 -Y -P -ax map-ont -t ${task.cpus} ${assembly} ${lr} \
    | samtools sort | samtools view -b -F 4 -o  ${id}_${type}_lr.bam 
    samtools index ${id}_${type}_lr.bam ${id}_${type}_lr.bam.bai
    """
}

process find_ovlp_reads {
// Creates circos file from bam, uses R script to find overlapping reads
    tag { id + ":" + contig_name }

    input:
    tuple val(id), path(lr), val(contig_name), val(length), val(seq), path(assembly), val(type), path(bam), path(bai)

    output:
    tuple val(id), val(contig_name), val(length), path('reads.txt'), path('ovlp.txt'), path('cov_ovlp.txt'), optional: true, emit: reads

    script:
    """
    ${params.env}
    bedtools bamtobed -i ${bam} > reads.bed
    echo -e ${contig_name}'\\t'\$(expr ${params.seqPadding} - 10)'\\t'\$(expr ${params.seqPadding} + 10) > breaks.bed
    echo -e ${contig_name}'\\t'\$(expr ${length} + ${params.seqPadding} - 10 )'\\t'\$(expr ${length} + ${params.seqPadding} + 10) >> breaks.bed
    samtools view -L breaks.bed -b ${bam} > region.bam
    intersectBed -wa -a reads.bed -b breaks.bed > ovlp.bed
    awk '{print \$4}' ovlp.bed | sort | uniq -D | uniq > readID.txt
    samtools view -H region.bam > ovlp.sam 
    samtools view region.bam | grep -f readID.txt >> ovlp.sam || true
    samtools view -b ovlp.sam > ovlp.bam
    samtools index ovlp.bam
    
    bedtools bamtobed -i ovlp.bam > ovlp_extracted.bed
     
    mosdepth -t ${task.cpus} -F 4 -n -b ${params.covWindow} ${contig_name} ovlp.bam
    gunzip -c ${contig_name}.regions.bed.gz > cov_ovlp.bed
    
    03_prepare_bed.R ovlp_extracted.bed ${params.seqPadding} ovlp.txt FALSE ${contig_name} ${length}
    03_prepare_bed.R cov_ovlp.bed ${params.seqPadding} cov_ovlp.txt TRUE
    03_prepare_bed.R reads.bed ${params.seqPadding} reads.txt FALSE ${contig_name} ${length}
    """
}

process identify_resistance_genes {
// Find antibiotic resistance genes in the CARD database
    publishDir({ "${params.outDir}/${id}/resistances" }, mode: 'copy')
    tag { id }

    input:
    tuple val(id), path(assembly), path(lr)
    
    output:
    tuple val(id), path("${id}_rgi.txt"), emit: rgi

    script:
    """
    ${params.env}
    rgi main -i ${assembly} -n ${task.cpus} -o ${id}_rgi
    """
}

process format_data_rgi {
// Converts gff file to circos readable format    
    tag { id }

    input:
    tuple val(id), path(rgi)

    output:
    tuple val(id), path('rgi.txt'), path('rgi_span.txt'), emit: circos_rgi

    script:
    """
    ${params.env}
    02_create_rgi_circos.R ${rgi}
    """
}

process mos_depth {
// Calculate coverage depth
    publishDir({ "${params.outDir}/${id}/coverage" }, mode: 'copy')
    tag { id }

    input:
    tuple val(id), path(assembly), val(type), path(aln_lr), path(aln_lr_idx)

    output:
    tuple val(id), path("${id}_cov_${type}.bed.gz"), val(type), emit: cov_bed

    script:
    """
    ${params.env}
    mosdepth -t ${task.cpus} -F 4  -n -b ${params.covWindow} ${id} ${aln_lr} 
    mv ${id}.regions.bed.gz ${id}_cov_${type}.bed.gz
    """
}

process format_data_cov {
// Formats coverage data for use in circos
    tag { id }
    
    input:
    tuple val(id), path(bed), val(type)

    output:
    tuple val(id), path('cov.txt'), val(type), emit: cov_formatted

    script:
    if (type == 'padded')
        """
        ${params.env}
        gunzip -c ${bed} > cov.bed
        03_prepare_bed.R cov.bed ${params.seqPadding} cov.txt TRUE
        """
    else
        """
        ${params.env}
        gunzip -c ${bed} > cov.bed
        03_prepare_bed.R cov.bed 0 cov.txt FALSE 
        """
}

process calcGC {
// Calculate gc conten
    publishDir({ "${params.outDir}/${id}/gc" }, mode: 'copy')
    tag { id }

    input:
    tuple val(id), path(assembly), val(lr), val(type)
    
    output:
    tuple val(id), path('gc1000.txt'), path(assembly), emit: table_gc
    tuple val(id), path('gc50.txt'), path('gc1000.txt'), path('gcskew50.txt'), path('gcskew1000.txt'), path('gcskewsum50.txt'), path('gcskewsum1000.txt'), emit: circos_gc

    script:
    """
    ${params.env}
    01_calculate_GC.R ${assembly} ${params.seqPadding}
    """
}

process glimmer {
// Predict gene positions with glimmer3
    publishDir({ "${params.outDir}/${id}/genes" }, mode: 'copy')
    tag { id }

    input:
    tuple val(id), path(assembly), path(lr)

    output:
    tuple val(id), path("${id}.predict"), emit: genes_glimmer
    path("${id}.detail")

    script:
    """
    ${params.env}
    long-orfs -n -t 1.15 ${assembly} ${id}.longorfs
    extract -t ${assembly} ${id}.longorfs > ${id}.train
    build-icm -r ${id}.icm < ${id}.train
    glimmer3 -o50 -g110 -t30 ${assembly} ${id}.icm ${id}
    """
}

process format_glimmer {
// Format predicted glimmer genes for circos
    tag { id }

    input:
    tuple val(id), path(genes)

    output: 
    tuple val(id), path('genes.txt'), emit: circos_genes

    script:
    """
    ${params.env}
    05_convert_glimmer.R ${genes}
    """
}

process circos{
// Use the combined data to create circular plots
    publishDir({ "${params.outDir}/${id}/plots" }, mode: 'copy')
    tag { id + ":" + contigID }

    input:
    tuple val(id), val(contigID), val(length), path(reads), path(ovlp), path(cov_ovlp), path(gc50), path(gc1000), path(gcskew50), path(gcskew1000), path(gcskewsum50), path(gcskewsum1000), path(cov), val(type), path(rgi), path(rgi_span), path(genes)

    output:
    path("${id}_${contigID}_plasmid.*")

    script:
    """
    ${params.env}
    echo "chr	-	${contigID}	1	0	${length}	chr1	color=lblue" > contig.txt
    ln -s ${workflow.projectDir}/conf/circos//* .
    circos
    mv circos.png ${id}_${contigID}_plasmid.png
    mv circos.svg ${id}_${contigID}_plasmid.svg
    """
}

process table{
// Create table with contig informations
    publishDir({ "${params.outDir}/${id}/" }, mode: 'copy')
    tag { id }

    input:
    tuple val(id), path(gc), path(assembly), path(cov), val(type), path(rgi)

    output:
    path("${id}_summary.csv")

    script:
    """
    ${params.env}
    04_summary_table.R ${assembly} ${rgi} ${cov} ${gc} ${params.seqPadding}
    mv contig_summary.txt ${id}_summary.csv
    """
}

workflow {
    if (params.help) exit 0, helpMessage()
    if (params.version) exit 0, pipelineMessage()
    if (!params.input) exit 0, helpMessage()

    def inputFile = workflow.profile.contains('test') ? file("${workflow.projectDir}/" + params.input) : file(params.input)
    def samples = getFiles(params.input)

    startMessage(inputFile)
    runParamCheck()

    def samples_filtered = filter_reads(samples).filtered
    def sample_channels = samples_filtered.multiMap {
        samples_rgi: it
        samples_glimmer: it
        samples_map: it
        samples_split: it
    }

    def contig_candidates = sample_channels.samples_split
        .splitFasta(record: [id: true, seqString: true])
        .map {
            def id = it[0]
            def lr = it[2]
            def contigName = it[1]['id']
            def length = it[1]['seqString'].length()
            def sequence = it[1]['seqString']
            [id, lr, contigName, length, sequence]
        }
        .filter { it[3] < params.maxLength }
        .filter { it[3] > params.minLength }

    def contig_channels = contig_candidates.multiMap {
        contigs: it
        contigs_2: it
        contigs_3: it
    }

    save_plasmids(contig_channels.contigs_3)

    def contigs_padded = pad_plasmids(contig_channels.contigs_2).padded
    def assembly_padded = combine_padded_contigs(contigs_padded.groupTuple()).padded_assembly
    def padded_channels = assembly_padded.multiMap {
        map_padded: it
        gc_padded: it
    }

    def to_mapping = sample_channels.samples_map
        .map { [it[0], it[1], it[2], 'normal'] }
        .mix(padded_channels.map_padded.map { [it[0], it[1], it[2][0], it[3]] })

    def bam_lr = map_longreads(to_mapping).bam
    def bam_channels = bam_lr.multiMap {
        bam_cov: it
        bam_ovlp: it
    }

    def circos_reads = find_ovlp_reads(contig_channels.contigs.combine(bam_channels.bam_ovlp.filter { it[2] == 'padded' }, by: 0)).reads

    def from_rgi = identify_resistance_genes(sample_channels.samples_rgi).rgi
    def rgi_channels = from_rgi.multiMap {
        rgi_txt: it
        table_data_rgi: it
    }
    def circos_data_rgi = format_data_rgi(rgi_channels.rgi_txt).circos_rgi

    def cov_formated = format_data_cov(mos_depth(bam_channels.bam_cov).cov_bed).cov_formatted
    def circos_data_cov = cov_formated.filter { it[2] == 'padded' }
    def table_data_cov = cov_formated.filter { it[2] != 'padded' }

    def calc_gc = calcGC(padded_channels.gc_padded)
    def table_data_gc = calc_gc.table_gc
    def circos_data_gc = calc_gc.circos_gc

    def circos_data_genes = format_glimmer(glimmer(sample_channels.samples_glimmer).genes_glimmer).circos_genes

    def circos_data = circos_data_gc
        .join(circos_data_cov)
        .join(circos_data_rgi)
        .join(circos_data_genes)

    def combined_data = circos_reads.combine(circos_data, by: 0)

    def table_data = table_data_gc
        .join(table_data_cov)
        .join(rgi_channels.table_data_rgi)

    circos(combined_data)
    table(table_data)
}


/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/
def runParamCheck() {
  if (params.minLength < params.seqPadding) {
	log.info "Minimum contig length cannot be shorter then seqPadding. Please adjust parameters and restart"
  }
}

def getFiles(tsvFile) {
  // Extracts Read Files from TSV
    def inputFile
  if (workflow.profile.contains('test')) {
      inputFile = file("${workflow.projectDir}/" + tsvFile)
  } else {
      inputFile = file(tsvFile)
  }
  log.info "------------------------------"
  Channel.fromPath(inputFile)
      .ifEmpty {exit 1,   "Cannot find path file ${tsvFile}"}
      .splitCsv(sep:'\t')
      .map { row ->
            [row[0], returnFile(row[1]), returnFile(row[2])]
            }   
}

def returnFile(it) {
// Return file if it exists and is readable
  if (workflow.profile.contains('test')) {
    return(file("$workflow.projectDir/" + it))
  }
  if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
  if (!file(it).canRead()) exit 1, "Cannot read file in TSV file: ${it}"
  return file(it)
}


def helpMessage() {
  // Display help message
  // pipelineMessage()
  log.info "  Usage:"
  log.info "       nextflow run caspargross/plasmident --input <file.csv> [options] "
  log.info "    --input <file.tsv>"
  log.info "       TSV file containing paths to files (id | assembly | longread)"
  log.info "  Parameters: "
  log.info "    --outDir "
  log.info "    Output locattion (Default: current working directory"
  log.info "    --maxLength <bases> (Default: 500000)"
  log.info "    Contigs larger then maxLength will not be considered a putative plasmid"
  log.info "    --seqPadding <bases> (Default: 2000)"
  log.info "    Length of recycled sequences at contig edges for long read mapping."
  log.info "    --covWindow <bases> (Default: 50)"
  log.info "    Moving window size for coverage calculation"
  log.info "    --mappingCov <coverage> (Default: 50)"
  log.info "    Target coverage for long read sampling"
  log.info "    --noSubsampling"
  log.info "    Skips the read subsampling step. Use when read coverage is not uniform."
  log.info "    --version"
  log.info "      Displays pipeline version"
  log.info "    --help"
  log.info "      Displays this help"
  log.info "           "
  log.info "  Profiles:"
  log.info "    -profile local "
  log.info "    Pipeline runs with locally installed conda environments (found in env/ folder)"
  log.info "    -profile test "
  log.info "    Runs complete pipeline on small included test dataset"
  log.info "    -profile localtest "
  log.info "    Runs test profile with locally installed conda environments"


}

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def minimalInformationMessage(inputFile) {
  // Minimal information message
  log.info "Command Line  : " + workflow.commandLine
  log.info "Input file    : " + inputFile
  log.info "Profile       : " + workflow.profile
  log.info "Project Dir   : " + workflow.projectDir
  log.info "Launch Dir    : " + workflow.launchDir
  log.info "Work Dir      : " + workflow.workDir
  log.info "Cont Engine   : " + workflow.containerEngine
  log.info "Out Dir       : " + params.outDir
  log.info "Align. Overlp.: " + params.seqPadding
  log.info "Cov. window   : " + params.covWindow
  log.info "Max Plasm. Len: " + params.maxLength
  log.info "Min Plasm. Len: " + params.minLength
  log.info "Target cov.   : " + params.mappingCov
  log.info "read sampling : " + !params.noSubsampling
  log.info "Containers    : " + workflow.container 
}

def pipelineMessage() {
  // Display hybridAssembly info  message
  log.info "PlasmIdent Pipeline ~  version ${workflow.manifest.version} - revision " + grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def startMessage(inputFile) {
  // Display start message
  asciiArt()
  pipelineMessage()
        minimalInformationMessage(inputFile)
}

def asciiArt() {
    log.info "       _                    ___________           _   "
    log.info "      | |                  |_   _|  _  \\         | |  "
    log.info " _ __ | | __ _ ___ _ __ ___  | | | | | |___ _ __ | |_ "
    log.info "| '_ \\| |/ _` / __| '_ ` _ \\ | | | | | / _ \\ '_ \\| __|"
    log.info "| |_) | | (_| \\__ \\ | | | | || |_| |/ /  __/ | | | |_ "
    log.info "| .__/|_|\\__,_|___/_| |_| |_\\___/|___/ \\___|_| |_|\\__|"
    log.info "| |                                                   "
    log.info "|_|                                                   "
}
