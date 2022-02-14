### quality control with fastqc
for i in *.fastq.gz; do 
fastqc $i
done

###cut adapter sequences from reads with cutadapt
#single read data:
for i in *.fastq.gz; do 
trim_galore --nextera $i
done

#paired-end data:
for read1 in read1*.fastq.gz; do 
read2=$(echo $read1| sed 's/read1/read2/g')
trim_galore --paired --nextera $read1 $read2
done


###alignment to reference genome with hisat2
wget <urlgenomehisat2>
tar zxvf genome.tar.gz

#single read data:
for i in *.fq.gz; do
hisat2 -q -x genome_tran -U $i -S $i.sam
done

#paired-end data:
for read1 in read1*.fq.gz; do 
read2=$(echo $read1| sed 's/read1/read2/g' | sed 's/val_1/val_2/g')
hisat2 -q -x genome_tran -1 $read1 -2 $read2 -S $read1.sam
done


###convert sam to bam
for i in *.sam;do
samtools view -bS -o $i.bam $i
done

###sort bam files
for i in *.bam;do
samtools sort -o sorted.$i $i
done

###get read counts matrix with featurecounts
featureCounts -p -M -t exon -g gene_id -a Homo_sapiens.GRCh37.75.gtf -o featcounts5.txt *.bam
