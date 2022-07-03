### get library size of FASTQs from FASTQC files

# /ark02/repository/do26181-do26258

out=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/scripts/dataProcessing

for i in $(seq 26181 26258)
do
  dir=/ark02/repository/do$i
  perl $out/01.4_countTotalReads_fromFASTQC.pl -dir $dir
done



