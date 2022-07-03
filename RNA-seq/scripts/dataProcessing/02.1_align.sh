### map each sample to the mouse reference genome
## enable quantification
## index the produced BAM files

star=/home/ibarra01/STAR-2.6.0c/bin/Linux_x86_64_static/STAR

for i in $(seq 16433 16442)
do
  dir=/user01/group_folders/Personal/Ximena/somitogenesis2020/RNA-seq/data/
  sample=do$i
  mkdir $dir"STAR/CONTROL/"$sample

  echo "Processing sample "$sample

  $star --genomeDir /user01/group_folders/Personal/Ximena/REFERENCE/STAR_MOUSE \
  --readFilesIn $dir"FASTq/"$sample".fq.gz" --readFilesCommand zcat \
  --outFileNamePrefix $dir"STAR/CONTROL/"$sample"/" --runThreadN 14 --outFilterMismatchNmax 1 \
  --outFilterMatchNminOverLread 0.5 --outFilterScoreMinOverLread 0.5 \
  --sjdbGTFfile /user01/group_folders/Personal/Ximena/REFERENCE/Mus_musculus.GRCm38.96.gtf \
  --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
  --outSAMstrandField intronMotif --quantMode GeneCounts --outReadsUnmapped None
  
  samtools index $dir"STAR/CONTROL/"$sample"/Aligned.sortedByCoord.out.bam"