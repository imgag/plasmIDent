#!/usr/bin/env/ nextflow

/* 
===============================================================================
                              P L A S M I D E N T  
===============================================================================
Nextflow pipeline for resistance plasmid identification and annotation using 
Nanopore reads and bacterial genome assemblies
-------------------------------------------------------------------------------
@ Author
Caspar Groß <mail@caspar.one>
-------------------------------------------------------------------------------
@ Documentation
https://github.com/caspargross/hybridassembly/README.md
------------------------------------------------------------------------------
*/

// Check special input parameters
if (params.help) exit 0, helpMessage()
if (params.version) exit 0, pipelineMessage()
if (!params.input) exit 0, helpMessage()

// Setup
samples = getFiles(params.input)
env = params.env
startMessage()
runParamCheck()

process filter_reads {
// Subsample large datasets to given target coverage
// Filtering with focus on high qual to keep enough short reads for plasmids
    tag{id}

    input:
    set id, assembly, lr from samples 

    output:
    set id,  assembly, file('reads_filtered.fastq') into samples_filtered

    script:
    if (params.noSubsampling)
        """
        zcat -f ${lr} > reads_filtered.fastq
        """
    else
        """
        ${env}
        len=\$(grep -v '>' ${assembly} | wc -c)
        nbases=\$(expr \$len * ${params.mappingCov})
        filtlong -t \$nbases --length_weight 0 ${lr} > reads_filtered.fastq
        """
}

// Duplicate channel
samples_filtered.into{samples_rgi; samples_glimmer; samples_map; samples_split}

// Split into contigs and filter for length channel
samples_split
    .map{[
        it[0],
        it[1],
        it[2]
        ]}
    .splitFasta(record: [id: true, seqString: true])
    .map{
        def id = it[0]
        def lr = it[2]
        //def contigName = it[1]['id'].replaceAll('_', '-')
        def contigName = it[1]['id']
        def length = it[1]['seqString'].length()
        def sequence = it[1]['seqString']
        [id, lr, contigName, length, sequence]
       }
    .filter{it[3] < params.maxLength}
    .filter{it[3] > params.minLength}
  //.view()
    .into{contigs; contigs_2; contigs_3}

process save_plasmids {
// Save plasmids as a separate fasta file
    tag{id + ":" + contigName} 
    publishDir "${params.outDir}/${id}/plasmids/", mode: 'copy'

    input:
    set id, lr, contigName, length, sequence from contigs_3
   
    output:
    file("${contigName}.fasta")

    script:
    """
    echo ">${contigName} len=${length}" > ${contigName}.fasta
    echo ${sequence} >> ${contigName}.fasta

    """
}

process pad_plasmids {
// Add prefix and suffix with sequence from oppsig end to each plasmid
    tag{id + ":" + contigName}

    input: 
    set id, lr, contigName, length, sequence from contigs_2

    output: 
    set id, file("${id}_${contigName}_padded.fasta"), lr, contigName into contigs_padded
    
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
    tag{id + ":" + contigName}
    publishDir "${params.outDir}/${id}/alignment/", mode: 'copy'

    input:
    set id, assembly, lr, contigName from contigs_padded.groupTuple()

    output:
    set id, file("${id}_padded.fasta"), lr, val("padded") into assembly_padded

    script:
    """
    cat \$(echo ${assembly} | tr -d '[],') > ${id}_padded.fasta 
    """
}

assembly_padded.into{map_padded; gc_padded}

// Mix channel with padded and normal contigs
samples_map
  //.view()
    .map{[it[0], 
        it[1], 
        it[2], 
        'normal']}
    .mix(map_padded
        .map{[it[0], 
            it[1], 
            it[2][0],
            it[3]]})
    .set{to_mapping}

process map_longreads {
// Use minimap2 to align longreads to padded contigs
    publishDir "${params.outDir}/${id}/alignment/", mode: 'copy'
    tag{id}

    input:
    set id, assembly, lr, type from to_mapping

    output:
    set id, assembly, type, file("${id}_${type}_lr.bam"), file("${id}_${type}_lr.bam.bai") into bam_lr

    script:
    """
    ${env}
    minimap2 -Y -P -ax map-ont -t ${task.cpus} ${assembly} ${lr} \
    | samtools sort | samtools view -b -F 4 -o  ${id}_${type}_lr.bam 
    samtools index ${id}_${type}_lr.bam ${id}_${type}_lr.bam.bai
    """
}

// Distribute bamfiles for coverage and read overlap identification
bam_cov = Channel.create()
bam_ovlp = Channel.create()
bam_lr.into{bam_cov; bam_ovlp}

process find_ovlp_reads {
// Creates circos file from bam, uses R script to find overlapping reads
    tag{id + ":" + contig_name}

    input:
    set id, lr, contig_name, length, seq, file(assembly), type, bam, bai from contigs.combine(bam_ovlp.filter{it[2] == 'padded'}, by : 0)

    output:
    set id, contig_name, length, file("reads.txt"), file("ovlp.txt"), file("cov_ovlp.txt") optional true into circos_reads 

    script:
    """
    ${env}
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
    publishDir "${params.outDir}/${id}/resistances", mode: 'copy'
    tag{id}

    input:
    set id, assembly, lr from samples_rgi
    
    output:
    set id, file("${id}_rgi.txt") into from_rgi

    script:
    """
    ${env}
    rgi main -i ${assembly} -n ${task.cpus} -o ${id}_rgi
    """
}

from_rgi.into{rgi_txt; table_data_rgi}

process format_data_rgi {
// Converts gff file to circos readable format    
    tag{id}

    input:
    set id, rgi from rgi_txt

    output:
    set id, file("rgi.txt"), file("rgi_span.txt") into circos_data_rgi

    script:
    """
    ${env}
    02_create_rgi_circos.R ${rgi}
    """
}

process mos_depth {
// Calculate coverage depth
    publishDir "${params.outDir}/${id}/coverage", mode: 'copy'
    tag{id}

    input:
    set id, assembly, type, file(aln_lr), file(aln_lr_idx) from bam_cov

    output:
    file("${id}_cov_${type}.bed.gz")
    set id, file("${id}_cov_${type}.bed.gz"), type into cov_bed

    script:
    """
    ${env}
    mosdepth -t ${task.cpus} -F 4  -n -b ${params.covWindow} ${id} ${aln_lr} 
    mv ${id}.regions.bed.gz ${id}_cov_${type}.bed.gz
    """
}

process format_data_cov {
// Formats coverage data for use in circos
    tag{id}
    
    input:
    set id, bed, type from cov_bed

    output:
    set id, file("cov.txt"), type into cov_formated

    script:
    if (type == "padded")
        """
        ${env}
        gunzip -c ${bed} > cov.bed
        03_prepare_bed.R cov.bed ${params.seqPadding} cov.txt TRUE
        """
    else
        """
        ${env}
        gunzip -c ${bed} > cov.bed
        03_prepare_bed.R cov.bed 0 cov.txt FALSE 
        """
}

// Distribute coverage file for circos (padded)  and summary table (normal)
circos_data_cov = Channel.create()
table_data_cov = Channel.create()
cov_formated.choice(circos_data_cov, table_data_cov) { it[2] == 'padded' ? 0 : 1 }

process calcGC {
// Calculate gc conten
    publishDir "${params.outDir}/${id}/gc", mode: 'copy'
    tag{id}

    input:
    set id, assembly, lr, type from gc_padded
    
    output:
    set id, file('gc1000.txt'), assembly into table_data_gc
    set id, file('gc50.txt'), file('gc1000.txt'), file('gcskew50.txt'), file('gcskew1000.txt'), file('gcskewsum50.txt'), file('gcskewsum1000.txt') into circos_data_gc

    script:
    """
    ${env}
    01_calculate_GC.R ${assembly} ${params.seqPadding}
    """
}

process glimmer {
// Predict gene positions with glimmer3
    publishDir "${params.outDir}/${id}/genes", mode: 'copy'
    tag{id}

    input:
    set id, assembly, lr from samples_glimmer

    output:
    set id, file("${id}.predict") into genes_glimmer
    file("${id}.detail")

    script:
    """
    ${env}
    long-orfs -n -t 1.15 ${assembly} ${id}.longorfs
    extract -t ${assembly} ${id}.longorfs > ${id}.train
    build-icm -r ${id}.icm < ${id}.train
    glimmer3 -o50 -g110 -t30 ${assembly} ${id}.icm ${id}
    """
}

process format_glimmer {
// Format predicted glimmer genes for circos
    tag{id}

    input:
    set id, genes from genes_glimmer

    output: 
    set id, file("genes.txt") into circos_data_genes

    script:
    """
    ${env}
    05_convert_glimmer.R ${genes}
    """
}

// Combine all finished circos data based on the id
circos_data_gc
   .join(circos_data_cov)
       .join(circos_data_rgi)
           .join(circos_data_genes)
           .set{circos_data}

// Combine contig data with sample wide circos data
combined_data = circos_reads.combine(circos_data, by: 0)

// Combine all table data based on id
table_data_gc
    .join(table_data_cov)
        .join(table_data_rgi)
        .set{table_data}

process circos{
// Use the combined data to create circular plots
    publishDir "${params.outDir}/${id}/plots", mode: 'copy'
    tag{id + ":" + contigID}

    input:
    set id, contigID, length, file(reads), file(ovlp), file(cov_ovlp), file(gc50), file(gc1000), file(gcskew50), file(gcskew1000), file(gcskewsum50), file(gcskewsum1000), file(cov), type, file(rgi), file(rgi_span), file(genes) from combined_data

    output:
    file("${id}_${contigID}_plasmid.*")

    script:
    """
    ${env}
    echo "chr	-	${contigID}	1	0	${length}	chr1	color=lblue" > contig.txt
    ln -s ${workflow.projectDir}/conf/circos//* .
    circos
    mv circos.png ${id}_${contigID}_plasmid.png
    mv circos.svg ${id}_${contigID}_plasmid.svg
    """
}

process table{
// Create table with contig informations
    publishDir "${params.outDir}/${id}/", mode: 'copy'
    tag{id}

    input:
    set id, gc, assembly, cov, type, rgi from table_data

    output:
    file("${id}_summary.csv")

    script:
    """
    ${env}
    04_summary_table.R ${assembly} ${rgi} ${cov} ${gc} ${params.seqPadding}
    mv contig_summary.txt ${id}_summary.csv
    """
}


/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/
def runParamCheck() {
  if (params.minLength > params.seqPadding) {
	log.info "Minimum contig length cannot be shorter then seqPadding. Please adjust parameters and restart"
  }
}

def getFiles(tsvFile) {
  // Extracts Read Files from TSV
  if (workflow.profile in ['test', 'localtest'] ) {
      inputFile = file("$workflow.projectDir/" + tsvFile)
  } else {
      inputFile = file(tsvFile)
  }
  log.info "------------------------------"
  Channel.fromPath(inputFile)
      .ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
      .splitCsv(sep:'\t')
      .map { row ->
            [id:row[0], assembly:returnFile(row[1]), lr:returnFile(row[2])]
            }   
}

def returnFile(it) {
// Return file if it exists and is readable
  if (workflow.profile in ['test', 'localtest'] ) {
      return(file("$workflow.projectDir/data/" + it))
  }
  if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
  if (!file(it).canRead()) exit 1, "Cannot read file in TSV file: ${it}"
  return file(it)
}


def helpMessage() {
  // Display help message
  // this.pipelineMessage()
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

def minimalInformationMessage() {
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
  log.info "PlasmIdent Pipeline ~  version ${workflow.manifest.version} - revision " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def startMessage() {
  // Display start message
  this.asciiArt()
  this.pipelineMessage()
  this.minimalInformationMessage()
}

workflow.onComplete {
  // Display complete message
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
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
