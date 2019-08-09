# Get read counts from fastq and bam files
# generated at each pipeline step
# to plot and see where reads are lost

# Get fastq file lines
for p in $(ls 1_fastq/); do 
	for f in $(ls 1_fastq/$p); do 
		sample=$(echo $f | cut -f1 -d ".")
		lines=$(( $(zcat 1_fastq/$p/$f | wc -l) / 4 ))
		echo -e $sample","$p","$lines","$(readlink -f 1_fastq/$p/$f)
	done
 done > fastq_reads.txt


# Get bam file lines
for p in $(ls 2_bam/); do
	for f in $(ls 2_bam/$p | grep -v bai); do
		sample=$(echo $f | cut -f1 -d "_")
		lines=$(samtools view -c 2_bam/$p/$f)
		echo -e $sample","$p","$lines","$(readlink -f 2_bam/$p/$f)
	done
done > bam_reads.txt


cat fastq_reads.txt bam_reads.txt > all_reads.txt

