### merge BAM files for each sample
## there are eleven files for each sample, coming from the 11 lanes used for sequencing.

# /ark02/repository/do26181-do26258

out=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA

for i in $(seq 26181 26258)
do
  dir=/ark02/repository/do$i
  
  files=$(ls -d -1 $dir/* | grep bam)
  samtools merge --threads 20 $out/do$i.bam $files
done