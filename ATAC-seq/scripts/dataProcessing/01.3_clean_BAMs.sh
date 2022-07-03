java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e1_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e1_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e1_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e1_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e2_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e2_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e2_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e2_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e3_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e3_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e3_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e3_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e5_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e5_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e5_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e5_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e6_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e6_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e6_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e9_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e9_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e9_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e9_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e13_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e13_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e13_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e13_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e14_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e14_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e14_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e15_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e15_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e15_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e16_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e16_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e16_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e16_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e17_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e17_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e17_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e17_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e19_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e19_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e19_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e19_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e20_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e20_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e20_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e20_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e21_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e21_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e21_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e21_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e22_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e22_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e22_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e22_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e23_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e23_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e23_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e23_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e24_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e24_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e24_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e24_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e25_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e25_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e25_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e25_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e25_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e25_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e25_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e25_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e25_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e25_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e26_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e26_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e26_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e26_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e27_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e27_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e27_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e27_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e28_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e28_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e28_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e28_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e29_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e29_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e29_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e29_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e30_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e30_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e30_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e31_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e31_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e31_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e31_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e32_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e32_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e32_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e32_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e33_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e33_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e33_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e33_SIII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SI-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SI-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e34_SI-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SI-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SI-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SI-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SI-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SI-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SI-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SI-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e34_SII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SII-2.noDUPs.bam

java -jar $picard MarkDuplicates I=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SIII-2.bam O=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SIII-2.noDUPs.bam M=/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/DUPSmetrics/e34_SIII-2.markDup_metrics.txt REMOVE_DUPLICATES=TRUE
samtools view -bh -q 30 -f 0x02 -F 0x800 -L /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/mm10_usedChr.bed /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SIII-2.noDUPs.bam > /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SIII-2.noDUPs.GQ.bam
samtools index /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SIII-2.noDUPs.GQ.bam
samtools flagstat /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SIII-2.noDUPs.GQ.bam >> /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SIII-2.flagstat
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SIII-2.bam
rm /user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2020/ATAC-seq/data/BWA/e34_SIII-2.noDUPs.bam

