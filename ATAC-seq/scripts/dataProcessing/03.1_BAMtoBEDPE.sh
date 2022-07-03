## convert BAM files into BEDPE

dir=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA

bams=$(ls -d -1 $dir/* | grep bam | grep -v bai)

for b in $bams
do
  sample="${b##*/}"
  sample="${sample%%.*}"
  echo Processing sample $sample
  samtools sort -n --threads 20 $b | bedtools bamtobed -bedpe -mate1 -i stdin > $dir/$sample.noDUPs.GQ.bedpe
done


