# ctg-seqonly

The ctg-seqonly pipeline performs demultiplexing and QC of NGS Illumina runs

Based on Nextflow + Singularity


## Input files


The following files must be in the runfolder to start pipeline successfully.

1. Samplesheet (CTG_SampleSheet.seqonly.csv)

### Samplesheet requirements:

!! NOTE: 
- `Sample_Name` does not have to be included in sheet - better to exclude. If it is included it MUST be identical to`Sample_ID`

- `ProjectID` has to be added (somewhere in the samplesheet above [Data]). The driver will use grep 'ProjectID' from the samplesheet, and take the value after the comma (`metaid=$(grep "ProjectID" $sheet | cut -f2 -d"," | tr -d '\n\r')`). So it has to be of the format "ProjectID,2021_0XX".

- `Adapter` under  `[Settings]` can be omitted.
- `[Reads]` Settings can also be omitted.

Illumina IEM based. The Bold/Italic field below must be correct! Other fields not in bold does not have to be changed for the pipeline to work, but can be changed if wanted.

[Header]         
***ProjectID,2021_0XX***     
IEMFileVersion,5  
Investigator Name,X  
Experiment Name,X  
Date,YYYY-MM-DD  
Workflow,GenerateFASTQ  
Application,NovaSeq FASTQ Only  
Instrument Type,NovaSeq  
Assay,Nextera XT  
Index Adapters,"Nextera XT v2 Index Kit A"  
Chemistry,Amplicon  
  
[Reads]  
***26***  
***78***  
  
[Settings]  
Adapter,***CTGTCTCTTATACACATCT***  
AdapterRead2,***CTGTCTCTTATACACATCT***  

[Data]  
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description  
***S1***,***S1***,,,***N702***,***CGTACTAG***,,,***2021_024***,  
***S2***,***S2***,,,***N706***,***TAGGCATG***,,,***2021_024***,  

 
 
### Samplesheet template 


```
[Header]
ProjectID,2021_0XX
IEMFileVersion,5  
Investigator Name,X  
Experiment Name,X  
Date,YYYY-MM-DD  
Workflow,GenerateFASTQ  
Application,NovaSeq FASTQ Only  
Instrument Type,NovaSeq  
Assay,Nextera XT  
Index Adapters,"Nextera XT v2 Index Kit A"  
Chemistry,Amplicon  
  
[Reads]  
26  
78  
  
[Settings]  
Adapter,CTGTCTCTTATACACATCT
AdapterRead2,CTGTCTCTTATACACATCT
[Data]  
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description  
S1,S1,,,N702,CGTACTAG,,,2021_024,  
S2,S2,,,N706,TAGGCATG,,,2021_024,  
```
## USAGE

1. Clone and build the Singularity container for this pipeline: https://github.com/perllb/ctg-seqonly/tree/main/container. Add the path to the .sif in the nextflow.config `container = ` parameter under process{}
2. Edit your samplesheet to match the example samplesheet. See section `SampleSheet` below
3. Edit the nextflow.config file to fit your project and system. 
4. Run pipeline (from runfolder)
```
nohup nextflow run pipe-seqonly-qc.nf > log.pipe-seqonly-qc.txt &
```

## Run with Driver
.. 1-3 from above  
4. Edit the driver folder declaration to fit your system  
5. Add cloned git repo folder to your system PATH  
6. Start driver from runfolder  

```
seqonly-driver 
```
- This will use default values (e.g. /path/to/runfolder/CTG_SampleSheet.csv)

## Usage Seqonly-Driver
```
Usage: seqonly-driver [ -i META_ID ] [ -s SAMPLESHEET ] [ -b BCL2FASTQ ARGUMENTS ] [ -r RESUME ] [ -d DEMUX-OFF ] [ -h HELP ] 


Optional arguments: 
META-ID       -i : Set 'meta-id' for runfolder (e.g. 210330-seqonly). Default: Takes date of runfolder (before first _) and adds '-seqonly' as suffix 
SAMPLESHEET   -s : Set samplesheet used for run (Default: CTG_SampleSheet.csv) 
BCL2FASTQ arg -b : String with bcl2fastq argument. e.g. '--minimum-trimmed-read-length 20 --mask-short-adapter-reads 20' 
RESUME        -r : Set if to resume nf-pipeline
DEMUX-OFF     -d : Set flag to skip bcl2fastq (then fastq must be in FQDIR) 
HELP          -h : print help message
```

### Pipeline steps:

* `Demultiplexing` (bcl2fastq): Converts raw basecalls to fastq, and demultiplex samples based on index (https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf).
* `FastQC`: FastQC calculates quality metrics on raw sequencing reads (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). MultiQC summarizes FastQC reports into one document (https://multiqc.info/).
* `multiQC`: Compile fastQC and cellranger count metrics in multiqc report (https://multiqc.info/)
* `md5sum`: md5sum of all fastq files


### Output:
* ctg-PROJ_ID-output
    * `qc`: Quality control output. 
        * fastqc output (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
        * multiqc output: Summarizing FastQC output and demultiplexing (https://multiqc.info/)
    * `fastq`: Contains raw fastq files from cellranger mkfastq.
    

  
### Container  
https://github.com/perllb/ctg-seqonly/tree/main/container  
