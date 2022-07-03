### count number of mapped reads from Picard's output

# /ark02/repository/do26181-do26258

dir=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics

files=$(ls -d -1 $dir/*)

for file in $files
do
  grep -A1 LIBRARY $file
done