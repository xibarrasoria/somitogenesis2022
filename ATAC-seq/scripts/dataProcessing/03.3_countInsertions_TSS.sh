## use bedtools to count the number of fragments overlapping 2kb fragments centered at the TSS of expressed genes

dir=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq

insertions=$(ls -d -1 $dir/data/BWA/* | grep insertionSites.bed)
intervals=$dir/data/TSScoords_expressedGenes_flanking1kb.bed 

for ins in $insertions
do
  sample="${ins##*/}"
  sample=${sample/.noDUPs.GQ.insertionSites.bed/}
  echo Processing sample $sample
  bedtools coverage -a $intervals -b $ins -d > $dir/TSS/$sample".TSScount.bed"
done
