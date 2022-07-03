### count number of mapped reads in final (good-quality) BAM file

dir=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA

files=$(ls -d -1 $dir/* | grep flagstat)

for file in $files
do
  sample="${file##*/}"
  echo "${sample%%.*}"
  grep QC $file
done
