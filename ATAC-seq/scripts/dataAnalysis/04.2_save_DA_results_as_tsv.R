## Generate text files with differential accessibility results
# from ATAC-seq/scripts/04_differentialAccessibilty.Rmd
# report subset of peaks most significant, to avoid very large files

library(csaw)

# dir <- "/Users/ibarra01/OneDrive - CRUK Cambridge Institute/github/somitogenesis2020/"
dir <- "/Users/xi629080/Documents/somitogenesis2022/"

## DA between somite trios
# per stage
somiteTrios.perStage <- readRDS(paste0(dir, "ATAC-seq/results/04_diffAccessibility_somiteTrios_perStage.Rds"))

for(contrast in names(somiteTrios.perStage)){
  tmp <- as.data.frame(somiteTrios.perStage[[contrast]])
  tmp <- tmp[tmp$PValue < 0.1,]
  write.table(tmp, paste0(dir, "ATAC-seq/results/04_diffAccessibility_somiteTrios_perStage.", contrast, ".tsv"),
              quote = FALSE, sep="\t", row.names = FALSE)
}

# average across stages
somiteTrios.average <- readRDS(paste0(dir, "ATAC-seq/results/04_diffAccessibility_somiteTrios_average.Rds"))

for(contrast in names(somiteTrios.average)){
  tmp <- as.data.frame(somiteTrios.average[[contrast]])
  tmp <- tmp[tmp$PValue < 0.1,]
  write.table(tmp, paste0(dir, "ATAC-seq/results/04_diffAccessibility_somiteTrios_average.", contrast, ".tsv"),
              quote = FALSE, sep="\t", row.names = FALSE)
}



## DA between stages
# per stage
stage.perSomite <- readRDS(paste0(dir, "ATAC-seq/results/04_diffAccessibility_stages_perSomite.Rds"))

for(contrast in names(stage.perSomite)){
  tmp <- as.data.frame(stage.perSomite[[contrast]])
  tmp <- tmp[tmp$PValue < 0.1,]
  write.table(tmp, paste0(dir, "ATAC-seq/results/04_diffAccessibility_stages_perSomite.", contrast, ".tsv"),
              quote = FALSE, sep="\t", row.names = FALSE)
}

# average across stages
stage.average <- readRDS(paste0(dir, "ATAC-seq/results/04_diffAccessibility_stages_average.Rds"))
tmp <- as.data.frame(stage.average)
tmp <- tmp[tmp$PValue < 0.1,]
write.table(tmp, paste0(dir, "ATAC-seq/results/04_diffAccessibility_stages_average.tsv"),
            quote = FALSE, sep="\t", row.names = FALSE)

