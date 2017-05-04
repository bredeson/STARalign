#!/bin/bash

###############################################################################
# Setup
###############################################################################

# Set hard and optional variables
# ============================================================
RUNTIME='12:0:0';
SJSUPPORT='20';
VERSION='0.1.1'; # Updated 05/03/17 before Jessen
DATE=`date +%F`;
COMMAND="$0 $*";
THREADS=0; 
# Load modules
# ============================================================
## SUGGESTION: let the user choose the versions of star and samtools. This
##   allows greater flexibility with later versions/releases. From a developer
##   perspective, your code will be more portable to systems that do not have
##   modules if you collect tool paths from the $PATH variable (via `which`)
##   at the start of execution.

MODULE_STAR="star/2.5.2b";
MODULE_SAMTOOLS="samtools/0.1.19";
module rm star;
module rm samtools;
module load ${MODULE_STAR};
module load ${MODULE_SAMTOOLS};

# Error function
# ============================================================
function error () {
    printf "[%s] ERROR: " `basename $0` >&2;
    echo "$2." >&2;
    exit $1;
}

# Help usage function
# ============================================================
function usage () {
    printf "\n" >&2;
    printf "%s v%s \n" `basename $0` $VERSION >&2;
    printf "\n" >&2;
    printf "Usage: %s [--workdir STR] [--RNA-seq-dir STR] [--ref-genome STR]\n" `basename $0` >&2;
    printf "       [--runtime HH:MM:SS] [--sj-support INT] [--threads INT] [--help] [-h]\n" >&2;
    printf "\n" >&2;
    printf "This script runs a pipeline to align RNA-seq reads to a genome using STAR.\n" >&2;
    printf "\n" >&2;
    printf "Required arguments:\n" >&2;
    printf "       --workdir STR       path for working directory\n" >&2;
    printf "       --RNA-seq-dir STR   path for directory with RNA-seq files\n" >&2;
    printf "       --ref-genome STR    path for genome fasta file\n" >&2;
    printf "\n" >&2;
    printf "Optional arguments:\n" >&2;
    printf "       --runtime HH:MM:SS  runtime for job submissions [12:0:0]\n" >&2;
    printf "       --sj-support INT    minimum support for splice junctions [20]\n" >&2;
    printf "       --threads INT       run pipeline locally with INT number of threads [SGE]\n" >&2;
    printf "       -h, --help          show this help message and exit\n" >&2;
    printf "\n" >&2;
    printf "Notes:\n" >&2;
    printf "   1. This script assumes a standard UNIX/Linux install. In addition, this script\n" >&2;
    printf "      assumes access to an SGE job-scheduling system that can be accessed by \'qsub\'.\n" >&2;
    printf "   2. This script assumes that fastq files end in the format [1,2][fastq format] and \n" >&2;
    printf "      are the only fastq files in the RNA-seq-dir and are all the same fastq format.\n" >&2;
    printf "      Accepted fastq formats are .fastq, .fq, .fastq.gz, and .fq.gz. These files\n" >&2;
    printf "      should already be adapter trimmed.\n" >&2;
    printf "   3. It is highly recommended that absolute paths be passed via the required\n" >&2;
    printf "      flags, as those paths may be used in a submission to the cluster.\n" >&2;
    printf "\n" >&2;
}

# Check that external tools are accessible
# ============================================================
CAT=`which cat`;
ZCAT=`which zcat`;
PERL=`which perl`;
QBATCH=`which qbatch`;
STAR=`which STAR`;
SAMTOOLS=`which samtools`;

if [[ -z "${PERL}" || ! -x "${PERL}" ]]; then
    error 127 "perl not in PATH env variable or not executable";

elif [[ -z "${CAT}" || ! -x "${CAT}" ]]; then
    error 127 "cat not in PATH env variable or not executable";

elif [[ -z "${ZCAT}" || ! -x "${ZCAT}" ]]; then
    error 127 "zcat not in PATH env variable or not executable";

elif [[ -z "${QBATCH}" || ! -x "${QBATCH}" ]]; then
    error 127 "qbatch not in PATH env variable or not executable";

elif [[ -z "${STAR}" || ! -x "${STAR}" ]]; then
    error 127 "star not in PATH env variable or not executable";

elif [[ -z "${SAMTOOLS}" || ! -x "${SAMTOOLS}" ]]; then
    error 127 "samtools not in PATH env variable or not executable";

fi

## Retrieve options on the command line and check for errors
# ============================================================
HELP_MESSAGE=;
while [[ -n $@ ]]; do
    case "$1" in
        '--workdir') shift; WORKDIR=$1;;
        '--RNA-seq-dir') shift; RNADIR=$1;;
        '--ref-genome') shift; GENOME=$1;;
        '--runtime') shift; RUNTIME=$1;;
        '--sj-support') shift; SJSUPPORT=$1;;
        '--threads') THREADS=$1;;
        '--help') HELP_MESSAGE=1;;
        '-h') HELP_MESSAGE=1;;
        -*) usage; error 2 "Invalid option: ${1}";;
        *) break;;
    esac;
    shift;
done

if [[ -n "${HELP_MESSAGE}" ]]; then
    usage;
    exit 1;

elif [[ -z "${WORKDIR}" || ! -d "${WORKDIR}" ]]; then
    usage; error 1 "WORKDIR not defined or does not exist";

elif [[ -z "${RNADIR}" || ! -d "${RNADIR}" ]]; then
    usage; error 1 "RNADIR not defined or does not exist";

elif [[ -z "${GENOME}" || ! -f "${GENOME}" ]]; then
    usage; error 1 "GENOME not defined or does not exist";

else
    printf "[%s] Starting %s %s %s %s %s %s\n" `basename $0` `date`     >&2;
    printf "[%s] Command-line: $COMMAND\n" `basename $0` >&2;
    printf "[%s] Version: $VERSION\n" `basename $0` >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "WORKDIR"   $WORKDIR   >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "RNADIR"    $RNADIR    >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "GENOME"    $GENOME    >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "RUNTIME"   $RUNTIME   >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "SJSUPPORT" $SJSUPPORT >&2;
fi

# Determine fastq format
# ============================================================
# what if the user has files with different suffixes (some .fq and some .fastq files?)
# You should put in the Usage page that you are assuming only one suffix type
# and that the user should consistently name their files (and/or compress them).
# also check explicitly and throw error if not consistent.

GZ=;
UNGZ=;
CATPATH=;
if [[ "$(find ${RNADIR}/* -maxdepth 0 -name '*.fastq.gz' | wc -l)" -ne 0 ]]; then
    FQFORMAT=".fastq.gz";
    CATPATH=$ZCAT;
    GZ=1;

elif [[ "$(find ${RNADIR}/* -maxdepth 0 -name '*.fq.gz' | wc -l)" -ne 0 ]]; then
    FQFORMAT=".fq.gz";
    CATPATH=$ZCAT;
    GZ=1;
fi
if [[ "$(find ${RNADIR}/* -maxdepth 0 -name '*.fastq' | wc -l)" -ne 0 ]]; then
    FQFORMAT=".fastq";
    CATPATH=$CAT;
    UNGZ=1;

elif [[ "$(find ${RNADIR}/* -maxdepth 0 -name '*.fq' | wc -l)" -ne 0 ]]; then
    FQFORMAT=".fq";
    CATPATH=$CAT;
    UNGZ=1;
fi

if [[ -n $GZ && -n $UNGZ ]]; then
    error 1 "Mixed compressed and uncompressed fastq files detected. Please standardize your files to either all compressed or uncompressed."

# we do not want to check explicitly that $THREADS (formerly $DEBUG) is 
# defined, as it will throw an error whenever we want to use the cluster
# ($THREADS is unset in that case)

elif [[ -z $GZ && -z $UNGZ ]]; then
    usage; error 1 "Fastq files could not be found in RNADIR";
fi


# Determine if fastq files are single, pair, or both
# ============================================================
# Best to start with the most restrictive case first. If there are no 
# second-end reads, it cannot be "pair" or "both". But the presence
# of second-end reads does not immediately inform us whether we have
# "pair" or "both".

if [[ "$(ls ${RNADIR}/*0${FQFORMAT} | wc -l)" -ne 0 ]]; then
    READEND="single";
fi

if [[ "$(ls ${RNADIR}/*1${FQFORMAT} | wc -l)" -ne 0 || "$(ls ${RNADIR}/*2${FQFORMAT} | wc -l)" -ne 0 ]]; then
    if [[ $READEND == "single" ]]; then
	READEND="both";
    else
	READEND="paired";
    fi
fi

if [[ $READEND == "pair" || $READEND == "both" ]]; then
    if [[ "$(ls ${RNADIR}/*1${FQFORMAT} | wc -l)" != "$(ls ${RNADIR}/*2${FQFORMAT} | wc -l)" ]]; then
	error 1 "Unequal number of paired files detected";
    fi
fi


# Determine max read length - 1
# ============================================================
READLEN=`for f in $RNADIR/*1$FQFORMAT; do $CATPATH ${f} | head -n1000; done | awk '{if (NR % 4 == 2) { print length($1)}}' | sort -k1,1nr | head -1`;


printf "[%s] PARAM: %s = %s\n" `basename $0` "READLEN"   $READLEN   >&2;
printf "[%s] PARAM: %s = %s\n" `basename $0` "READEND"   $READEND   >&2;
printf "[%s] PARAM: %s = %s\n" `basename $0` "FQFORMAT"  $FQFORMAT  >&2;

# Make directories
# ============================================================
mkdir -p $WORKDIR/Index/1st_Index;
mkdir -p $WORKDIR/Index/2nd_Index;
mkdir -p $WORKDIR/Split_Fastq;
mkdir -p $WORKDIR/RNA_Mapping/1st_Mapping;
mkdir -p $WORKDIR/RNA_Mapping/2nd_Mapping;
mkdir -p $WORKDIR/Final_Bam;


###############################################################################
# First Round:
###############################################################################

# Split transcriptome files
# ============================================================
printf "[%s] Splitting transcriptome files \n" `basename $0` >&2;
if [[ -f "${WORKDIR}/Split_Files.submit" && -f "${WORKDIR}/Split_Files.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else
    for f in $RNADIR/*$FQFORMAT; do 
	sample=`basename $f $FQFORMAT`; 
	echo "$CATPATH $f | split -dl 4000000 -a 4 - $WORKDIR/Split_Fastq/${sample}_";
    done >$WORKDIR/Split_Files.submit;
    
    $QBATCH submit -T $THREADS -W -n Split_Files.submit $WORKDIR/Split_Files.submit $WORKDIR/BATCH_Split;

    if [[ "$(grep 'done' ${WORKDIR}/BATCH_Split/Split_Files.submit.e* | wc -l)" == "$(find ${WORKDIR}/BATCH_Split/* -maxdepth 0 -name 'Split_Files.submit.e*' | wc -l)" ]]; then
	$QBATCH touch $WORKDIR/Split_Files.done;
    else
	error 1 "Splitting transcriptome files failed; please identify error and restart";
    fi
fi

# Index genome
# ============================================================
printf "[%s] FIRST ROUND \n" `basename $0` >&2;
printf "[%s] Indexing genome \n" `basename $0` >&2;
if [[ -f "${WORKDIR}/1st_Index.submit" && -f "${WORKDIR}/1st_Index.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else
    echo "cd $WORKDIR && $STAR --runMode genomeGenerate --runThreadN 64 --genomeDir $WORKDIR/Index/1st_Index --genomeFastaFiles $GENOME --limitGenomeGenerateRAM=225000000000" >$WORKDIR/1st_Index.submit;
    $QBATCH submit -T $THREADS -W -R 7 -p 32 -t $RUNTIME -n Index_1st.submit $WORKDIR/1st_Index.submit $WORKDIR/BATCH_1st_Index;

    if [[ "$(grep 'done' ${WORKDIR}/BATCH_1st_Index/1st_Index.submit.e* | wc -l)" == "$(find ${WORKDIR}/BATCH_1st_Index/* -maxdepth 0 -name '1st_Index.submit.e*' | wc -l)" ]]; then
	$QBATCH touch $WORKDIR/1st_Index.done;
    else
	error 1 "Indexing genome failed; please identify error and restart";
    fi
fi

# Align split transcriptome files
# ============================================================
printf "[%s] Aligning split transcriptome files \n" `basename $0` >&2;
if [[ -f "$WORKDIR/1st_Align.submit" && -f "$WORKDIR/1st_Align.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else

    # For paired reads
    # ============================================================
    # The use of the >> (append) redirection requires us to first
    # make sure to clear any previous .batch files or the commands
    # will be duplicated and clobber eachother when run simultaneously
    # on the cluster.
    printf '' >$WORKDIR/1st_Align.submit;  # clear previous runs

    if [[ $READEND == "pair" || $READEND == "both" ]]; then
	# grep -v 1_ is not a very specific pattern. There may be unexpected 
	# side-effects when using this expression on something containing
	# that substring, eg. replicate1_muscle.2_0001, etc. use rev | cut -c7- | rev
	# to remove the _[0-9][0-9][0-9][0-9] at the end:
	ls $WORKDIR/Split_Fastq/*1_[0-9][0-9][0-9][0-9] | rev | cut -c7- | rev | sort | uniq -c | awk '{print $2"\t"$1-1}' | \
	while read sample N; do 
	    for n in `seq 0 $N`; do 
		n=`printf %04d $n`; 
		echo "$STAR --genomeDir $WORKDIR/Index/1st_Index --readFilesIn $WORKDIR/Split_Fastq/${sample}1_${n} $WORKDIR/Split_Fastq/${sample}2_${n} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $WORKDIR/RNA_Mapping/1st_Mapping/${sample}_${n} --runThreadN 32 --chimOutType SeparateSAMold --chimSegmentMin 20 --chimJunctionOverhangMin 20 --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No"; 
	    done; 
	done >$WORKDIR/1st_Align.submit;
    fi

    # For single reads
    # ============================================================
    if [[ $READEND == "single" || $READEND == "both" ]]; then
	# Are we expecting single files to have the same naming convetion as paired?
	ls $WORKDIR/Split_Fastq/*0_[0-9][0-9][0-9][0-9] | rev | cut -c7- | rev | sort | uniq -c | awk '{print $2"\t"$1-1}' | \
        while read sample N; do 
	    for n in `seq 0 $N`; 
	    do n=`printf %04d $n`;
		echo "$STAR --genomeDir $WORKDIR/Index/1st_Index --readFilesIn $WORKDIR/Split_Fastq/${sample}1_${n} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $WORKDIR/RNA_Mapping/1st_Mapping/${sample}_${n} --runThreadN 32 --chimOutType SeparateSAMold --chimSegmentMin 20 --chimJunctionOverhangMin 20 --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No"; 
	    done; 
	done >>$WORKDIR/1st_Align.submit;
    fi

    $QBATCH submit -T $THREADS -W -t $RUNTIME -p 16 -R 5 -S 4 -n Align_1.submit $WORKDIR/1st_Align.submit $WORKDIR/BATCH_1st_Align;

    if [[ "$(grep -o 'done' ${WORKDIR}/BATCH_1st_Align/1st_Align.submit.e* | wc -l)" == "$(find ${WORKDIR}/Split_Fastq/* -maxdepth 0 -name '*1_*' | wc -l)" ]]; then
	$QBATCH touch $WORKDIR/1st_Align.done;
    else
	error 1 "Aligning split transcriptome files failed; please identify error and restart";
    fi
fi

# Identify splice junctions
# ============================================================
printf "[%s] Identifying splice junctions \n" `basename $0` >&2;
if [[ -f "${WORKDIR}/sjdb.out" && -f "${WORKDIR}/Find_SJ.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else
    cat $WORKDIR/RNA_Mapping/1st_Mapping/*SJ.out.tab | awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if ($7>0) print $1,$2,$3,strChar[$4],$5,$7}' | $PERL -ane '$d{"$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]"} += $F[5]; END{map{print("$_\t$d{$_}\n")}keys(%d)}' | sort -k1,1 -k2,2n | awk -v sjs="$SJSUPPORT" '{if($6>sjs) print $0}' >$WORKDIR/sjdb.out;
    if [[ "$(find ${WORKDIR}/* -maxdepth 0 -type f -name 'sjdb.out' -size +1c | wc -l)" -ne 0 ]]; then
	$QBATCH touch $WORKDIR/Find_SJ.done;
    else
	error 1 "Identifying splice junctions failed; please identify error and restart";
    fi
fi


###############################################################################
# Second Round:
###############################################################################

# Index genome
# ============================================================
printf "[%s] SECOND ROUND \n" `basename $0` >&2;
printf "[%s] Indexing genome \n" `basename $0` >&2;
if [[ -f "${WORKDIR}/2nd_Index.submit" && -f "${WORKDIR}/2nd_Index.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else
    echo "cd $WORKDIR && $STAR --runMode genomeGenerate --runThreadN 64 --genomeDir $WORKDIR/Index/2nd_Index --genomeFastaFiles $GENOME --sjdbFileChrStartEnd $WORKDIR/sjdb.out --limitGenomeGenerateRAM=225000000000 --sjdbOverhang $READLEN" >$WORKDIR/2nd_Index.submit;

    $QBATCH submit -T $THREADS -W -R 7 -p 32 -t $RUNTIME -n Index_2nd.submit $WORKDIR/2nd_Index.submit $WORKDIR/BATCH_2nd_Index;

    if [[ "$(grep 'done' ${WORKDIR}/BATCH_2nd_Index/2nd_Index.submit.e* | wc -l)" != "$(find ${WORKDIR}/BATCH_2nd_Index/* -maxdepth 0 -name '2nd_Index.submit.e*' | wc -l)" ]]; then
	error 1 "Indexing genome failed; please identify error and restart";
    else
	$QBATCH touch $WORKDIR/2nd_Index.done;
    fi
fi

# Align split transcriptome files
# ============================================================
printf "[%s] Aligning split transcriptome files \n" `basename $0` >&2;
if [[ -f "$WORKDIR/2nd_Align.submit" && -f "$WORKDIR/2nd_Align.done" ]]; then
    printf "[%s] ... Already done. \n" `basename $0` >&2;
else

    # For paired reads
    # ============================================================
    printf '' >2nd_Align.submit;
    if [[ $READEND == "pair" || $READEND == "both" ]]; then
	# using grep -v 1_ is problematic...
	ls $WORKDIR/Split_Fastq/*1_[0-9][0-9][0-9][0-9] | rev | cut -c7- | rev | sort | uniq -c | awk '{print $2"\t"$1-1}' | \
        while read sample N; do 
	    for n in `seq 0 $N`; do
		n=`printf %04d $n`; 
		echo "$STAR --genomeDir $WORKDIR/Index/2nd_Index --readFilesIn $WORKDIR/Split_Fastq/${sample}1_${n} $WORKDIR/Split_Fastq/${sample}2_${n} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $WORKDIR/RNA_Mapping/2nd_Mapping/${sample}_${n} --runThreadN 32 --chimOutType SeparateSAMold --chimSegmentMin 20 --sjdbFileChrStartEnd $WORKDIR/sjdb.out --chimJunctionOverhangMin 20 --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No"; 
	    done; 
	done >$WORKDIR/2nd_Align.submit;
    fi

    # For single reads
    # ============================================================
    if [[ $READEND == "single" || $READEND == "both" ]]; then
	ls $WORKDIR/Split_Fastq/*0_[0-9][0-9][0-9][0-9] | rev | cut -c7- | rev | sort | uniq -c | awk '{print $2"\t"$1-1}' | \
	while read sample N; do 
	    for n in `seq 0 $N`; do 
		n=`printf %04d $n`; 
		echo "$STAR --genomeDir $WORKDIR/Index/2nd_Index --readFilesIn $WORKDIR/Split_Fastq/${sample}1_${n} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $WORKDIR/RNA_Mapping/2nd_Mapping/${sample}_${n} --runThreadN 32 --chimOutType SeparateSAMold --chimSegmentMin 20 --sjdbFileChrStartEnd $WORKDIR/sjdb.out --chimJunctionOverhangMin 20 --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No"; 
	    done; 
	done >>$WORKDIR/2nd_Align.submit;
    fi

    $QBATCH submit -T $THREADS -W -t $RUNTIME -p 16 -R 5 -S 4 -n Align_2.submit $WORKDIR/2nd_Align.submit $WORKDIR/BATCH_2nd_Align;

    if [[ "$(grep -o 'done' ${WORKDIR}/BATCH_2nd_Align/2nd_Align.submit.e* | wc -l)" == "$(find ${WORKDIR}/Split_Fastq/* -maxdepth 0 -name '*1_*' | wc -l)" ]]; then
	$QBATCH touch $WORKDIR/2nd_Align.done;
    else
	error 1 "Aligning split transcriptome files failed; please identify error and restart";
    fi
fi


###############################################################################
# Finish:
###############################################################################

# Merge second round bam files
# ============================================================  
printf "[%s] Merging final bam files \n" `basename $0` >&2;
if [[ "$(find ${WORKDIR}/Final_Bam/* -maxdepth 0 -name '*bam' -size +1c | wc -l)" != "$(find ${RNADIR}/* -maxdepth 0 -name '*1${FQFORMAT}' | wc -l)" ]]; then
    ls $WORKDIR/Split_Fastq/*[012]_[0-9][0-9][0-9][0-9] | rev | cut -c8- | rev | sort | uniq -c | awk '{print $2"\t"$1-1}' | \
    while read sample; do 
	$SAMTOOLS merge -f -@ 10 $WORKDIR/Final_Bam/${sample}.bam $WORKDIR/RNA_Mapping/2nd_Mapping/${sample}*Aligned.sortedByCoord.out.bam; 
    done;
else
    error 1 "Merging final bam files failed; please identify error and restart";
fi

if [[ "$(find ${WORKDIR}/Final_Bam/* -maxdepth 0 -name '*bam' -size +1c | wc -l)" -ne 0 ]]; then
    printf "[%s] Finished %s %s %s %s %s %s\n" `basename $0` `date` >&2;
else
    error 1 "Merging final bam files failed; please identify error and restart";

exit 0;
