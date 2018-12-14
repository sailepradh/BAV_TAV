## Codes for the generation of HiC mapping and HiCap interaction files

#### *HiC Alignment*

This folder consists of scripts for generating alignments files (bam) of the promoter mediated HiC from three different biological samples.

* Few points to remember in these scripts
    * The first script merges the individual fastq.gz lane files to give the combined R1 and R2 fastq.gz files. After that the HiCUP truncator script is run.

    * The second script maps the above created truncated fastq.gz files to reference genome with Bowtie2. This is kind time consuming steps .

    * The third script processes the above mapped bam files .

More information about the HiCUP are provided on the following site [HiCUP documentation](https://www.bioinformatics.babraham.ac.uk/projects/hicup/read_the_docs/html/index.html)

#### *HiCAP ProximityDetector interaction calling*

The HiCaptools designs probes and gives the proximity statistics of the fragments.
ProbeDesigner and ProximityDetector then takes the bam files and gives the probes and rest of genomes.

* Note while running HiCaptools in uppmax *
    * Gitclone the [HiCaptools](https://github.com/sahlenlab/HiCapTools)
    * module load cmake/3.7.2  gcc/4.9.2
    * Compile HiCaptools using the shell script  *./buildHiCapTools.sh*
    *  *export LD_LIBRARY_PATH=/path/to/HiCapTools/bamtools/:$LD_LIBRARY_PATH*
    * The executable hicaptools are in bin folder
    * ./bin/HiCapTools -h

Something is not right in here. I guess the earlier HiCap was run with different probe than the new ones. Currently I am testing with the earlier probe set but must be checked meticulously.

More scripts and codes will be updated as we go forward with the project
