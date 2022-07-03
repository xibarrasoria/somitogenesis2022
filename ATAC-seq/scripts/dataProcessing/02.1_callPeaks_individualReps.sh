## use MACS2 to call peaks in each sample
## this serves as a proxy for the signal-to-noise ratio

dir=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq

bams=$(ls -d -1 $dir/data/BWA/* | grep bam | grep -v bai)

for b in $bams
do
  sample="${b##*/}"
  echo Processing sample "${sample%%.*}"
  macs2 callpeak -f BAMPE -g mm --keep-dup all --broad -t $b --outdir $dir/peaks/individualReps/ -n $sample
done


