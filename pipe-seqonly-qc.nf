#!/usr/bin/env nextFlow

// set variables
seq = params.sequencing
runfolder = params.runfolder
basedir = params.basedir
metaID = params.metaid
OUTDIR = params.outdir
FQDIR = params.fqdir
QCDIR = params.qcdir
FQCDIR = params.fqcdir
CTGQC = params.ctgqc
demux = params.demux
b2farg= params.bcl2fastq_arg

// Read and process sample sheet
sheet = file(params.sheet)

// create new file for reading into channels that provide sample info!
newsheet = file("${basedir}/samplesheet.nf.seqonly.csv")

// Read and process sample sheet
allLines = sheet.readLines()
writeB = false // if next lines has sample info
newsheet.text=""     

for ( line in allLines ) {
    if ( writeB ) {
	newsheet.append(line + "\n")
    }
    if (line.contains("[Data]")) {
	writeB = true
    }
}

println "============================="
println ">>> sc-rna-10x pipeline for multiple projects / run >>>"
println ""
println "> INPUT: "
println "> sequnecing		: $seq "
println "> runfolder		: $runfolder "
println "> sample-sheet		: $sheet "
println "> run-meta-id		: $metaID "
println "> basedir		: $basedir "
println "> bcl2fastq args	: $b2farg "
println ""
println "> OUTPUT: "
println "> output-dir		: $OUTDIR "
println "> fastq-dir		: $FQDIR "
println "> qc-dir		: $QCDIR "
println "> fastqc-dir		: $FQCDIR "
println "> ctg-qc-dir		: $CTGQC "
println "============================="

// sample info
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Project ) }
    .unique()
    .tap{infoSamples}
    .set{ fastqc_ch;  }

println " > Samples to process: "
infoSamples.subscribe{ println "Sample: $it" }


// bcl2fastq
process bcl2fastq {

    tag "$metaID"

    input:
    val sheet 

    output:
    val "x" into fastqc_go
    val "x" into md5_ch

    when:
    demux == "y"
        
    """
    
    bcl2fastq -R $runfolder \\
              --sample-sheet $sheet \\
              --no-lane-splitting  \\
              -r 1 \\
              -p $task.cpus  \\
              -w 1  \\
              --output-dir $FQDIR \\
	      $b2farg
   
     """
}

// fastqc 
process fastqc {

	tag "${sid}-${projid}"

	input:
	val y from fastqc_go
	set sid, projid from fastqc_ch

        output:
        val sid into multiqc_fastqc

	when:
	demux == "y"
	
	"""

        mkdir -p ${QCDIR}
        mkdir -p ${FQCDIR}

	read1=\$(echo ${FQDIR}/$projid/$sid*R1*fastq.gz)
   	read2=\$(echo ${FQDIR}/$projid/$sid*R2*fastq.gz)

    	# Check if fastq is not found with pattern above (due to sample fastq are not put in sample id folder
    	if [[ \$read1 == *"*R1*"* ]]; then
       	   read1=\$(echo ${FQDIR}/${projid}/${sid}/${sid}*R1*fastq.gz)
       	   read2=\$(echo ${FQDIR}/${projid}/${sid}/${sid}*R2*fastq.gz)
    	fi

	fastqc -t $task.cpus --outdir $FQCDIR \$read1
	fastqc -t $task.cpus --outdir $FQCDIR \$read2
 	
	"""
    
}

process multiqc {

    tag "$metaID"

    input:
    val y from multiqc_fastqc.collect()

    """
    
    cd $OUTDIR
    multiqc -f ${OUTDIR} ${runfolder}/ctg-interop --outdir ${QCDIR}/multiqc/ -n multiqc_seqonly_${metaID}_report.html

    mkdir -p ${CTGQC}

    cp -r ${QCDIR} ${CTGQC}/

    """

}

process md5sum_fastq {

    tag "$metaID"

    input:
    val x from md5_ch

    """
    cd $FQDIR
    find -type f \\( -name "*fastq.gz" \\) -exec md5sum '{}' \\; > ctg-md5.fastq.txt
    
    """
}


