#!/bin/bash


echo "WARNING: use https://github.com/CSB5/misc-scripts/blob/master/bwamem_pe.sh instead"  1>&2
exit 1

set -o pipefail

# Decontaminate FastQ (single-end, paired-end, optionally gzipped) by
# mapping against a contamination source

echoinfo() {
    echo "INFO: $@" 1>&2
}
echofatal() {
    echo "FATAL: $@" 1>&2
}

# defaults
# illumina quality encoding
illumina=0
# number of threads to use
threads=4
# keep temp directory
keep_temp=0
# default picard directory; override with env var PICARDDIR
PICARDDIR_DEFAULT="/mnt/software/stow/picard-1.74/bin/"

usage() {
    # keep in sync with arg parsing below
cat <<EOF
$(basename $0): decontaminate fastq

Performs a mapping of given SR/PE reads (gzip supported) with BWA
against given source of contamination and produces a BAM file with
contaminated reads and new fastq file(s) with clean reads (qualities
will be Sanger encoded).

Prerequisites: BWA, Picard and samtools. Point the PICARDDIR
environment variable to your picard installation (will use
$PICARDDIR_DEFAULT otherwise)

  Mandatory options:
    -f | --fastq1    : Input fastq[.gz] file
    -r | --ref       : Reference (contamination source) fasta file
    -o | --outprefix : Output prefix
  Optional:
    -g | --fastq2    : Fastq[.gz], second in pair (optional)
    -t | --threads   : Number of threads to use (default=$threads)
    -I | --illumina  : Phred qualities are ASCII64, ie. Illumina 1.3-1.7 instead of Sanger (check with FastQC)
    -T | --tmpdir    : Use this as temp directory instead of automatically determined one
    -K | --keep      : Keep temp directory
    -B | --reusebam  : Reuse already created BAM, which contains unmapped reads and reads mapped against the contaminant
                       Still needs original fastq files for auto setting of output filenames and determining SR or PE.
    -h | --help      : Display this help
EOF
}



# check for required programs
#
for bin in java samtools bwa; do
    if ! which $bin >/dev/null 2>&1; then
        echofatal "Couldn't find $bin. make sure it's in your path"
        exit 1
    fi
done
test -z "$PICARDDIR" && export PICARDDIR=$PICARDDIR_DEFAULT
# check for picard needed for samtofastq adding (needed for GATK)
picard_samtofastq_jar=${PICARDDIR}/SamToFastq.jar
if [ ! -s $picard_samtofastq_jar ]; then
    echofatal "Couldn't find Picard's $(basename picard_samtofastq_jar). Please set PICARDDIR to your Picard installation"
    exit 1
fi


# parse arguments
#
reusebam=""
reffa=""
fastq1=""
fastq2=""
while [ "$1" != "" ]; do
    case $1 in
        -f | --fastq1 )
            shift
            fastq1=$1
            test -e $fastq1 || exit 1
            ;;
        -g | --fastq2 )
            shift
            fastq2=$1
            test -e $fastq2 || exit 1
            ;;
        -h | --help )
            usage
            exit
            ;;
        -I | --illumina )
            illumina=1
            ;;
        -K | --keep )
            keep_temp=1
            ;;
        -o | --outprefix )           
            shift
            outprefix=$1
            ;;
        -r | --ref )           
            shift
            reffa=$1
            test -e $reffa || exit 1
            ;;
        -t | --threads )
            shift
            threads=$1
            ;;
        -T | --tmpdir )
            shift
            tmpdir=$1
            ;;
        -B | --reusebam )
            shift
            reusebam=$1
            test -e $reusebam || exit 1
            ;;
        * ) 
            echofatal "unknown argument \"$1\""
            usage
            exit 1
    esac
    shift
done

# check arguments
# 
if [ -z "$fastq1" ]; then
    echofatal "fastq file \"$fastq1\" missing"
    usage
    exit 1
fi
if [ -z "$reffa" ] && [ -z "$reusebam" ]; then
    echofatal "reference fasta file \"$reffa\" missing"
    usage
    exit 1
fi
if [ -z "$outprefix" ]; then
    echofatal "output prefix missing"
    usage
    exit 1
fi

if [ -z "$tmpdir" ]; then
    tmpdir=$(mktemp --tmpdir -d "$(basename $0).XXXXXX")
fi
echoinfo "Using temp dir $tmpdir"

# set first batch of output file names that we need early on
contbam=${outprefix}_cont.bam
fastq_clean_1=${outprefix}_1.fastq
if [ ! -z "$fastq2" ]; then
    fastq_clean_2=${outprefix}_2.fastq
fi

for f in $fastq_clean_1 ${fastq_clean_1}.gz $contbam; do
    if [ -s $f ]; then
        echofatal "refusing to overwrite already existing file $f"
        exit 1
    fi
done


#
# Now after all this paranoia checking and user friendliness, can we
# please get on with it?
#

# index reference if necessary
if [ -n "$reffa" ]; then
    test -s ${reffa}.pac || bwa index $reffa || exit 1
fi
allbam=$tmpdir/all.bam
if [ -n "$reusebam" ]; then
    allbam=$reusebam
fi
bwalog=$tmpdir/bwa.log


# STEP 1: bwa aln
# --------------------------------------------------------------------

# -q 3: remove q3 in accordance with illumina guidelines
bwa_aln_extra_args="-t $threads -q 3"
if [ $illumina -eq 1 ]; then
    # -I: if phred qualities are ascii64, ie. Illumina 1.3-1.7
    bwa_aln_extra_args="$bwa_aln_extra_args -I"
fi

# bwa aln for each fastq. skip if sai already exists
sais=""
for fastq in $fastq1 $fastq2; do
    sai=$tmpdir/$(basename $fastq .gz | sed -e 's,.fastq$,,' | sed -e 's,.txt$,,').sai
    sais="$sais $sai"
    if [ -s "$sai" ]; then
        echoinfo "Reusing already existing $sai"
        continue
    fi
    if [ -s "$allbam" ]; then
        echoinfo "Skipping alignment step (bwa aln) and reusing already existing $allbam"
        continue
    fi
#cat <<EOF
    if ! bwa aln $bwa_aln_extra_args -f $sai $reffa $fastq >> $bwalog 2>&1; then
        echofatal "bwa aln failed. See $bwalog"
        exit 1
    fi
#EOF
done
echoinfo "BWA aln done"



# STEP 2: bwa sam[sp]e
# --------------------------------------------------------------------

# SR/PE specific options to bwa
bwa_samse_extra_args=""
bwa_sampe_extra_args="-s"
# -s disable Smith-Waterman for the unmapped mate
if [ -z "$fastq2" ]; then
    args="samse $bwa_samse_extra_args $reffa $sais $fastq1"
    # remove unmapped reads from single end mapping
else
    args="sampe $bwa_sampe_extra_args $reffa $sais $fastq1 $fastq2"
    # keep only reads mapped in proper pair
fi

if [ -s "$allbam" ]; then
    echoinfo "Skipping sam conversion and reusing already existing $allbam"
else
#cat <<EOF
    if ! bwa $args 2>>$bwalog | samtools view -bS - > $allbam; then
        echofatal "bwa or samtools failed. see $bwalog"
        exit 1
    fi
#EOF
fi
echoinfo "BWA sam[sp]e done"



# STEP 3: split into contaminated reads (BAM) and clean FastQ
# --------------------------------------------------------------------

# In case of PE reads, we consider clean if neither read nor mate are
# mapped you can give -F/-f only once to samtools). Another advantage
# of using awk is that you can avoid a nasty BWA bug, which results in
# inconsistent mate-pair info. See http://www.biostars.org/p/60100/
# and http://www.biostars.org/p/60001/
#
# Use picard for bam2fastq because it's extra pedantic. Good
# alternative would be
# http://www.hudsonalpha.org/gsl/software/bam2fastq.php which has the
# (dis)advantage of simply skipping unpaired mate-pairs. With
# "LENIENT" Picard will ignore errors like "Error parsing text SAM
# file. MAPQ should be 0 for unmapped read." for cicular chromosomes.
#
# All communication via fifo to save space and only save BAM/gz.
# Disadvantage: can't catch errors directly, only via log.

# setup fifos
fastq_clean_1_fifo=$tmpdir/$(basename $fastq_clean_1 fastq)fifo
if [ ! -z "$fastq2" ]; then
    fastq_clean_2_fifo=$tmpdir/$(basename $fastq_clean_2 fastq)fifo
fi
contbam_fifo=$tmpdir/$(basename $contbam bam)fifo
fifos="$fastq_clean_1_fifo $fastq_clean_2_fifo $contbam_fifo"
for fifo in $fifos; do
	test -e $fifo && rm $fifo
	mkfifo $fifo
done

execlog=${outprefix}_decont-exec.log
samtools view -bS - < $contbam_fifo > $contbam 2>>$execlog &
gzip < $fastq_clean_1_fifo > ${fastq_clean_1}.gz 2>>$execlog &
if [ -z "$fastq2" ]; then
# && [ ! -s "$allbam" ]; then
    echoinfo "SE mode"
else
    echoinfo "PR mode"
    gzip < $fastq_clean_2_fifo > ${fastq_clean_2}.gz 2>>$execlog &
    samtofastq_pe_arg="SECOND_END_FASTQ=$fastq_clean_2_fifo"
fi


echoinfo "Splitting into contaminated BAM and clean fastq..."
samtools view -h $allbam 2>>$execlog | \
    awk -v fifo=$contbam_fifo 'BEGIN {FS="\t"; OFS="\t"}
{if (/^@/ && substr($2, 3, 1)==":") {print >> fifo; next}
if (( (!and($2, 0x1) && and($2, 0x4)) || (and($2, 0x1) && (and($2, 0x4) || and($2, 0x8))) ) && $3=="*") {print} else {print >> fifo}}' | \
    java -jar $picard_samtofastq_jar \
        VERBOSITY=ERROR QUIET=TRUE \
        INCLUDE_NON_PF_READS=true VALIDATION_STRINGENCY=LENIENT \
        INPUT=/dev/stdin FASTQ=$fastq_clean_1_fifo $samtofastq_pe_arg 2>>$execlog

rm $fifos
grep -v '[samopen] SAM header is present:' $execlog 1>&2

if [ $keep_temp -ne 1 ]; then
    test -d $tmpdir && rm -rf $tmpdir
else
    echoinfo "Keeping $tmpdir"
fi
echoinfo "Successful exit"
