library(dplyr)

# read .tsv files
# Vertebrates
ampCit <- read.csv("~/biohazard/data/arcCen_ampZal_ampCit/sbst3_arcCen_ampZal_ampCit_20241121_ampCit.tsv", sep="\t", header=TRUE)
ampZal <- read.csv("~/biohazard/data/arcCen_ampZal_ampCit/sbst3_arcCen_ampZal_ampCit_20241121_ampZal.tsv", sep="\t", header=TRUE)

# Animals
batBro <- read.csv("~/biohazard/data/gigPla_batSep_batBro/mut3_gigPla_batSep_batBro_20240610_batBro.tsv", sep="\t", header=TRUE)
batSep <- read.csv("~/biohazard/data/gigPla_batSep_batBro/mut3_gigPla_batSep_batBro_20240610_batSep.tsv", sep="\t", header=TRUE)
mytEdu <- read.csv("~/biohazard/data/mytTro_mytEdu_mytGal/mytTro_mytEdu_mytGal_20240513_mytEdu.tsv", sep="\t", header=TRUE)
mytGal <- read.csv("~/biohazard/data/mytTro_mytEdu_mytGal/mytTro_mytEdu_mytGal_20240513_mytGal.tsv", sep="\t", header=TRUE)
oikAlb <- read.csv("~/biohazard/data/oikDio_oikAlb_oikVan/oikDio_oikAlb_oikVan_20240513_oikAlb.tsv", sep="\t", header=TRUE)
oikVan <- read.csv("~/biohazard/data/oikDio_oikAlb_oikVan/oikDio_oikAlb_oikVan_20240513_oikVan.tsv", sep="\t", header=TRUE)

# Algae
ostLuc <- read.csv("~/biohazard/data/ostMed_ostTau_ostLuc/mut3_ostMed_ostTau_ostLuc_20240423_ostLuc.tsv", sep="\t", header=TRUE)
ostTau <- read.csv("~/biohazard/data/ostMed_ostTau_ostLuc/mut3_ostMed_ostTau_ostLuc_20240423_ostTau.tsv", sep="\t", header=TRUE)
proBov <- read.csv("~/biohazard/data/proCut_proCif_proBov/mut3_proCut_proCif_proBov_20240423_proBov.tsv", sep="\t", header=TRUE)
proCif <- read.csv("~/biohazard/data/proCut_proCif_proBov/mut3_proCut_proCif_proBov_20240423_proCif.tsv", sep="\t", header=TRUE)
proBovCh <- read.csv("~/biohazard/data/proCutCh_proCifCh_proBovCh/mut3_proCutCh_proCifCh_proBovCh_20240903_proBovCh.tsv", sep="\t", header=TRUE)
proCifCh <- read.csv("~/biohazard/data/proCutCh_proCifCh_proBovCh/mut3_proCutCh_proCifCh_proBovCh_20240903_proCifCh.tsv", sep="\t", header=TRUE)
ulvCom <- read.csv("~/biohazard/data/ulvPro_ulvMut_ulvCom/mut3_ulvPro_ulvMut_ulvCom_20240610_ulvCom.tsv", sep="\t", header=TRUE)
ulvMut <- read.csv("~/biohazard/data/ulvPro_ulvMut_ulvCom/mut3_ulvPro_ulvMut_ulvCom_20240610_ulvMut.tsv", sep="\t", header=TRUE)
ulvComCh <- read.csv("~/biohazard/data/ulvProCh_ulvMutCh_ulvComCh/mut3_ulvProCh_ulvMutCh_ulvComCh_20240826_ulvComCh.tsv", sep="\t", header=TRUE)
ulvMutCh <- read.csv("~/biohazard/data/ulvProCh_ulvMutCh_ulvComCh/mut3_ulvProCh_ulvMutCh_ulvComCh_20240826_ulvMutCh.tsv", sep="\t", header=TRUE)
walHed <- read.csv("~/biohazard/data/walMel_walIch_walHed/mut3_walMel_walIch_walHed_20240611_walHed.tsv", sep="\t", header=TRUE)
walIch <- read.csv("~/biohazard/data/walMel_walIch_walHed/mut3_walMel_walIch_walHed_20240611_walIch.tsv", sep="\t", header=TRUE)

# Fungi
aspChe <- read.csv("~/biohazard/data/aspCos_aspChe_aspCri/mut3_aspCos_aspChe_aspCri_20240617_aspChe.tsv", sep="\t", header=TRUE)
aspCri <- read.csv("~/biohazard/data/aspCos_aspChe_aspCri/mut3_aspCos_aspChe_aspCri_20240617_aspCri.tsv", sep="\t", header=TRUE)

# Bacteria
parMar <- read.csv("~/biohazard/data/cyaGra_parMar_proMar/mut3_cyaGra_parMar_proMar_20240911_parMar.tsv", sep="\t", header=TRUE)
proMar <- read.csv("~/biohazard/data/cyaGra_parMar_proMar/mut3_cyaGra_parMar_proMar_20240911_proMar.tsv", sep="\t", header=TRUE)
theSic <- read.csv("~/biohazard/data/pseAzo_theVes_theSic/mut3_pseAzo_theVes_theSic_20240910_theSic.tsv", sep="\t", header=TRUE)
theVes <- read.csv("~/biohazard/data/pseAzo_theVes_theSic/mut3_pseAzo_theVes_theSic_20240910_theVes.tsv", sep="\t", header=TRUE)

readFile <- function(filePath) {
  data <- read.csv(filePath, sep="\t", header=TRUE)
  df <- data.frame(data)
  df <- cbind(df, SbstPerOrig=as.numeric((df[,"mutNum"]/df[,"totalRootNum"])*100))
  df <- df[, c("mutType","SbstPerOrig")]
}

return df

data_list <- list(df1, df2, df3)
# combine multiple data frames by "mutType"
merged_df <- Reduce(function(x, y) merge(x, y, by = "mutType"), data_list)
# 列名を適切に変更（オプション）
colnames(merged_df) <- c("mutType", "Sample1", "Sample2", "Sample3")
pca_matrix <- as.matrix(merged_df[,-1])  # mutType列を除外
rownames(pca_matrix) <- merged_df$mutType  # mutTypeを行名に設定
pca_matrix <- t(pca_matrix)

pca_result <- prcomp(pca_matrix, scale. = TRUE)

summary(pca_result)

#　あとは可視化
