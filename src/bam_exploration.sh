total=$(samtools view -c results/alignment/ERR458495Aligned.sortedByCoord.out.bam)
mapq0=$(samtools view -c results/alignment/ERR458495Aligned.sortedByCoord.out.bam)
mapq_ge1=$(samtools view -c -q 1 results/alignment/ERR458495Aligned.sortedByCoord.out.bam)

echo "Total: $total"
echo "MAPQ 0: $((total - mapq_ge1))"
