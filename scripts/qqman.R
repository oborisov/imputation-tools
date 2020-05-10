library(qqman)
library(data.table)
# phenotypes full names
sumstats=fread(snakemake@input[[1]])
colnames(sumstats)[1:3]=c("CHR", "BP", "SNP")
# adjusting CHR and p
sumstats[, CHR := as.numeric(CHR)]
sumstats=sumstats[!is.na(sumstats$P) & P > 0]
# producing Manhattan plot
jpeg(paste0(snakemake@input[[1]], "_manh.jpeg"), width = 8, height = 4, units = "in", res = 300)
print(manhattan(rbind(sumstats[P<5e-2], sumstats[P>5e-2][seq(1,nrow(sumstats[P>5e-2]),10)]), chr="CHR", bp="BP", p="P", snp="SNP", annotatePval = 1, annotateTop = T))
dev.off()
# producing Q-Q plot
