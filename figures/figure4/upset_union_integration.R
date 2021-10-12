library(ComplexHeatmap)
library(dplyr)
library(openxlsx)

# we need the whole set of DEGs (1046 vs 745 oder so)
ago21_specifics <- read.xlsx("../../data/TableS1_RNA-seq.xlsx",
                             sheet = "Ago2&1_KO specific DEGs",
                             startRow=3)

# we need to load all tabs in the excel sheet
analyses = c("H3K27me3 ChIP-seq", "ATAC-seq", "diffTF-derived", "KLF4 ChIP-seq", "CTCF ChIP-seq")
for (status in c("UP", "DOWN")) {
  ago21_filtered <- filter(ago21_specifics, Status == status)

  lt = sapply(analyses, function(analysis) {

    df <- read.xlsx("../../data/TableS5_Ago2and1_DEGs_Integration.xlsx",
                    sheet = analysis,
                    startRow=3)
    filtered <- filter(df, Status == status)
    return(filtered$Gene.name)
  })
  mat <- make_comb_mat(lt, mode = "union", universal_set = ago21_filtered$Gene.name)
  pdf(paste0("analyses_overlap_", status, "s.pdf"), width=8, height=3)
  UpSet(mat, comb_order = order(comb_size(mat), decreasing=TRUE))
  dev.off()
}

