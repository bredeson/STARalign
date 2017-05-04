STARalign.sh v0.1.1 

Usage: STARalign.sh [--workdir STR] [--RNA-seq-dir STR] [--ref-genome STR]
       [--runtime HH:MM:SS] [--sj-support INT] [--help] [-h]

Pipeline to align RNA-seq reads to a genome using STAR with splice junctions.

Required arguments:

	--workdir STR       path for working directory
	--RNA-seq-dir STR   path for directory with RNA-seq files
	--ref-genome STR    path for genome fasta file

Optional arguments:

	--runtime HH:MM:SS  runtime for job submissions [12:0:0]
	--sj-support INT    minimum support for splice junctions [20]
	--debug INT         debug mode to run test data on local number of threads
	-h, --help          show this help message and exit

Notes:
   1. This script assumes a standard UNIX/Linux install. In addition, this script
      assumes access to an SGE job-scheduling system that can be accessed by 'qsub'.
   2. This script assumes that fastq files end in the format [1,2][fastq format] and 
      are the only fastq files in the RNA-seq-dir and are all the same fastq format.
      Accepted fastq formats are .fastq, .fq, .fastq.gz, and .fq.gz. These files
      should already be adapter trimmed.
   3. It is highly recommended that absolute paths be passed via the required
      flags, as those paths may be used in a submission to the cluster.

Test data:
   1. Human RNA-seq data is available from https://github.com/griffithlab/rnaseq_tutorial/wiki/RNAseq-Data 
   at http://genome.wustl.edu/pub/rnaseq/data/practical.tar.
   2. The human genome can be downloaded from the 1000 Genomes Project at ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/.
