#!/bin/bash

# Decontaminate FastQ
#
# Largely based on bwa_unique.sh

# defaults
illumina=0
threads=4
keep_temp=0
PICARDDIR_DEFAULT="/mnt/software/stow/picard-1.74/bin/"

usage() {
    # keep in sync with arg parsing below
cat <<EOF
$(basename $0): wrapper decontaminating a fastq file

Performs a mapping of given SR/PE reads (gzip supported) with BWA
against given reference / source of contamination and produces a BAM
file with contaminated reads and a new fastq file with clean reads
(qualities will be Sanger encoded).

Prerequisites: BWA, Picard and samtools. Point the PICARDDIR
environment variable to your picard installation (will use
$PICARDIR_DEFAULT otherwise)

  Mandatory options:
    -f | --fastq1    : Input fastq[.gz] file
    -r | --ref       : Reference fasta file
    -o | --outprefix : Output prefix
  Optional:
    -h | --help      : Display this help
    -g | --fastq2    : Fastq[.gz], second in pair (optional)
    -k | --keep      : Keep temp directory
    -t | --threads   : Number of threads to use (default=$threads)
         --illumina  : Phred qualities are ASCII64, ie. Illumina 1.3-1.7 instead of Sanger (check with FastQC)
         --tmpdir    : Use this as temp directory instead of automatically determined one
EOF
}



# check for required programs
#
for bin in java samtools; do
    if ! which $bin >/dev/null 2>&1; then
        echo "FATAL: couldn't find $bin. make sure it's in your path" 1>&2
        exit 1
    fi
done
test -z "$PICARDDIR" && export PICARDDIR=$PICARDDIR_DEFAULT
# check for picard needed for samtofastq adding (needed for GATK)
picard_samtofastq_jar=${PICARDDIR}/SamToFastq.jar
if [ ! -s $picard_samtofastq_jar ]; then
    echo "FATAL: couldn't find Picard's $(basename picard_samtofastq_jar). Please set PICARDDIR to your Picard installation" 1>&2
    exit 1
fi


# parse arguments
#
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
        --illumina )
            illumina=1
            ;;
        --keep )
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
        --tmpdir )
            shift
            tmpdir=$1
            ;;
        * ) 
            echo "FATAL: unknown argument \"$1\"" 1>&2
            usage
            exit 1
    esac
    shift
done

# check arguments
# 
if [ -z "$fastq1" ]; then
    echo "FATAL: fastq file \"$fastq1\" missing" 1>&2
    echo
    usage
    exit 1
fi
if [ -z "$reffa" ]; then
    echo "FATAL: reference fasta file \"$reffa\" missing" 1>&2
    usage
    exit 1
fi
if [ -z "$outprefix" ]; then
    echo "FATAL: output prefix missing" 1>&2
    usage
    exit 1
fi

if [ -z "$tmpdir" ]; then
    tmpdir=$(mktemp --tmpdir -d "$(basename $0).XXXXXX")
fi

allbam=$tmpdir/all.bam
cleanbam=$tmpdir/clean.bam
contbam=${outprefix}_cont.bam

fastq_clean_s1=${outprefix}_s1.fastq
if [ ! -z $fastq2 ]; then
    fastq_clean_s2=${outprefix}_s2.fastq
fi

if [ -s $fastq_clean_s1 ] || [ -s ${fastq_clean_s1}.gz ]; then
    echo "ERROR: refusing to overwrite already existing $fastq_clean_s1" 1>&2
    exit 1
fi


# index reference if necessary
test -s ${reffa}.pac || bwa index $reffa || exit 1


#
# Now after all this paranoia checking and user friendliness, can we
# please get on with it?
#

# remove q3 in accordance with illumina guidelines
bwa_aln_extra_args="-t $threads -q 3"
if [ $illumina -eq 1 ]; then
    bwa_aln_extra_args="$bwa_aln_extra_args -I"
    # -I: if phred qualities are ascii64, ie. Illumina 1.3-1.7
fi

# bwa aln for each fastq. skip if sai already exists
sais=""
for fastq in $fastq1 $fastq2; do
    sai=$tmpdir/$(basename $fastq .gz | sed -e 's,.fastq$,,' | sed -e 's,.txt$,,').sai
    sais="$sais $sai"
    if [ -s "$sai" ]; then
        echo "INFO: reusing already existing $sai" 1>&2
        continue
    fi
    if [ -s "$allbam" ]; then
        echo "INFO: skipping alignment step (bwa aln) and reusing already existing $allbam" 1>&2
        continue
    fi
#cat <<EOF
    bwa aln $bwa_aln_extra_args -f $sai $reffa $fastq || exit 1;
#EOF
done


# SR/PE specific options to bwa
bwa_samse_extra_args=""
bwa_sampe_extra_args="-s"
# -s disable Smith-Waterman for the unmapped mate
if [ -z $fastq2 ]; then
    args="samse $bwa_samse_extra_args $reffa $sais $fastq1"
    # remove unmapped reads from single end mapping
else
    args="sampe $bwa_sampe_extra_args $reffa $sais $fastq1 $fastq2"
    # keep only reads mapped in proper pair
fi


# -f INT: Only output alignments with all bits in INT present in the FLAG field
# -F INT: Skip alignments with bits present in INT [0]
# 0x4 segment unmapped
# 0x8 next segment in the template unmapped
#
# In case of PE reads, we consider clean if neither read nor mate are
# mapped read 0x4 and mate 0x8 unmapped. both required. therefore -f
# 12 Note: -f 12 is not the same as -f 4 and -f 8 (let alone 0x12)
# Numbers in clean and cont do not necessarily tally.
#
# Contaminated are then reads that map or where mate maps:
# -F 8 -F 4



# Note: numbers don't add up (not sure why)
#
# grep ^Total SRR341581_1_fastqc/fastqc_data.txt 
# Total Sequences16615077
# echo 16615077*2 | bc -l
# 33230154
# samtools view -c  -F 12  /tmp/decont.sh.4zE4mU/all.bam  
# 38
# samtools view -c  -F 8  /tmp/decont.sh.4zE4mU/all.bam  
# 1880
# samtools view -c  -F 4  /tmp/decont.sh.4zE4mU/all.bam  
# 1880
# samtools view -c  -F 8 -F 4  /tmp/decont.sh.4zE4mU/all.bam  
# 1880
# samtools view -c  -f 12  /tmp/decont.sh.4zE4mU/all.bam  
# 33226432
# samtools view -c  -f 4 -f 8  /tmp/decont.sh.4zE4mU/all.bam  
# 33228274
# echo 33228274+1880 | bc 
# 33230154


if [ -s "$allbam" ]; then
    echo "INFO: skipping sam conversion and reusing already existing $allbam" 1>&2
else
#cat <<EOF
    bwa $args | samtools view -bS - > $allbam || exit 1
#EOF
fi

#samtools index $allbam
if [ -z $fastq2 ]; then
    clean_flag='-f 4'
    cont_flag='-F 4'
else
    clean_flag='-f 12'
    cont_flag='-F 8 -F 4'
fi
#cat <<EOF
samtools view -b $clean_flag $allbam > $cleanbam
samtools view -b $cont_flag $allbam > $contbam
#EOF
echo "WARN: samtools filtering step untested on SE read mapping" 1>&2


if [ ! -z $fastq2 ]; then
    pe_arg="SECOND_END_FASTQ=$fastq_clean_s2"
fi

# samtofastq has no gzip support. could use /dev/stdout and | gzip for
#first fastq but wouldn't be able to catch errors then

#cat <<EOF
java -jar $picard_samtofastq_jar \
        INPUT=$cleanbam FASTQ=$fastq_clean_s1 $pe_arg INCLUDE_NON_PF_READS=true || exit 1
gzip ${outprefix}_s[12].fastq
#EOF

if [ $keep_temp -ne 1 ]; then
    test -d $tmpdir && rm -rf $tmpdir
else
    echo "Keeping $tmpdir"
fi
