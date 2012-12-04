bam="$1"
test -n "$bam" || exit 1
samtools view -H $bam | \
    awk '/^@SQ/ {for (i=2; i<=NF; i++) {if ($i ~ /^SN:/) {print substr($i, 4, length($i))}}}'
