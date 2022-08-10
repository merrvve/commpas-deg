library(dplyr)
meta2$fage50.<- ifelse(meta2$age<50,"young","old")
meta2$fage50. <- as.factor(meta2$fage50.)
tb1<-table(meta2$fage50.,meta2$gender)
chisq.test(tb1)
tb1<-table(meta2$fage50.,meta2$race)
fisher.test(tb1)
tb1<-table(meta2$fage50.,meta2$iss)
fisher.test(tb1)
tb1<-table(meta2$fage50.,meta2$age)
t.test(tb1)
library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData=counts,colData=meta2,design =~age+fage50.)
rld<-vst(dds,blind=TRUE)
plotPCA(rld,intgroup="fage50.")
plotPCA(rld,intgroup="gender")
plotPCA(rld,intgroup="iss")
plotPCA(rld,intgroup="race")
rld_mat<-assay(rld)
rld_cor<-cor(rld_mat)
library(pheatmap)
pheatmap(rld_cor)
pheatmap(rld_cor)
dds<-DESeq(dds)
plotDispEsts(dds)
contrast<-c("fage50.","young","old")
res1<-results(dds,contrast=contrast,alpha=0.05)
summary(res1)
sig_res1<- data.frame(res1) %>% filter(padj<0.05&abs(log2FoldChange)>0.58)
res2<-results(dds,contrast=contrast,alpha=0.05,lfcThreshold = 0.58)
normalized_counts<-counts(dds,normalized=TRUE)
normalized_counts<-data.frame(normalized_counts)
write.csv(normalized_counts,"c:/users/merve/documents/commpass-deg2/results/norm_counts.csv")
write.csv(res2_sig,"c:/users/merve/documents/commpass-deg2/results/sig2.csv")
write.csv(sig_res1,"c:/users/merve/documents/commpass-deg2/results/sig1.csv")
top20<-head(arrange(sig_res1,padj),20)
top20norm<-filter(normalized_counts,rownames(normalized_counts) %in% rownames(top20))
library(tidyr)
top20norm$gene<-rownames(top20norm)
top20gather<-gather(top20norm,colnames(top20norm)[1:634],key="samplename",value="normcount")

library(ggplot2)
meta2$samplename<-rownames(meta2)
joined<-inner_join(meta2,top20gather,by="samplename")
View(joined)
library(ggplot2)

write.csv(joined,"c:/users/merve/documents/commpass-deg2/results/top20toplot.csv")
library(EnsDb.Hsapiens.v86)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= joined$gene, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
View(geneIDs)
colnames(geneIDs)[2]="gene"
joined<-inner_join(joined,geneIDs,by="gene")

ggplot(joined) +
geom_point(aes(x = SYMBOL, y = normcount, color = fage50.)) +
scale_y_log10() +
xlab("Genes") +
ylab("log10 Normalized Counts") +
ggtitle("Top 20 Significant DE Genes") +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
theme(plot.title = element_text(hjust = 0.5))
colnames(joined)[6]="AgeCat"
ggplot(joined) +
geom_point(aes(x = SYMBOL, y = normcount, color = fage50.)) +
scale_y_log10() +
xlab("Genes") +
ylab("log10 Normalized Counts") +
ggtitle("Top 20 Significant DE Genes") +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
theme(plot.title = element_text(hjust = 0.5))
colnames(joined)[6]="AgeCat40"
colnames(joined)[8]="AgeCat"
ggplot(joined) +
geom_point(aes(x = SYMBOL, y = normcount, color = AgeCat)) +
scale_y_log10() +
xlab("Genes") +
ylab("log10 Normalized Counts") +
ggtitle("Top 20 Significant DE Genes") +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
theme(plot.title = element_text(hjust = 0.5))
View(top20norm)
res1<-data.frame(res1)
res1 <- mutate(res1, thresh=padj<0.05 & abs(log2FoldChange)>0.58)
ggplot(res1) +
geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = thresh)) +
ggtitle("Young Patients DEG") +
xlab("log2 fold change") +
ylab("-log10 adjusted p-value") +
#scale_y_continuous(limits = c(0,50)) +
theme(legend.position = "none",
plot.title = element_text(size = rel(1.5), hjust = 0.5),
axis.title = element_text(size = rel(1.25)))
mutate(arrange(res1,padj),genelabel="")
res1<-mutate(arrange(res1,padj),genelabel="")
res1$gene<-rownames(res1)
res1<-left_join(res1,geneIDs,by="gene")
res1$genelabel[1:10]<-res1$SYMBOL[1:10]
library(ggrepel)
ggplot(res1, aes(x = log2FoldChange, y = -log10(padj))) +
geom_point(aes(colour = thresh)) +
geom_text_repel(aes(label = genelabel)) +
ggtitle("Young Patients DEG") +
xlab("log2 fold change") +
ylab("-log10 adjusted p-value") +
theme(legend.position = "none",
plot.title = element_text(size = rel(1.5), hjust = 0.5),
axis.title = element_text(size = rel(1.25)))


