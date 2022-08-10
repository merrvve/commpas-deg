library(EnsDb.Hsapiens.v86)
sig_deg <- read.csv("c:/users/merve/documents/commpass-deg2/results/sig1.csv", row.names=1)
library(dplyr)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(sig_deg), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
write.csv(geneIDs,"c:/users/merve/documents/commpass-deg2/results/sig-deg-genesymbolv86.csv")
sig_deg$geneid<-rownames(sig_deg)
sig_deg<- left_join(sig_deg,geneIDs,by=c("geneid"="GENEID"))
write.csv(sig_deg,"c:/users/merve/documents/commpass-deg2/results/sig-deg-withsymbol.csv")
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
countsp <- read.csv("c:/users/merve/documents/commpass-deg2/results/norm_counts.csv", row.names=1)
all_genes<-rownames(countsp)
ego <- enrichGO(gene = sig_deg$geneid,
                universe = all_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)
save(ego, file="c:/users/merve/documents/commpass-deg2/data/ego.rda")
dotplot(ego, showCategory=50, label_format=70, font.size=9)
foldchanges <- sig_deg$log2FoldChange
names(foldchanges)<-sig_deg$gene
library(ggnewscale)
cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange=foldchanges,
vertex.label.font=6)
library(org.Hs.eg.db)
library(AnnotationDbi)
sig_deg$entrezid<- mapIds(org.Hs.eg.db, keys = sig_deg$geneid,
column = "ENTREZID", keytype = "ENSEMBL")
res_entrez <- dplyr::filter(sig_deg, entrezid != "NA")
## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrezid) == F), ]
foldchanges <- res_entrez$log2FoldChange
names(foldchanges) <- res_entrez$entrezid
foldchanges <- sort(foldchanges, decreasing = TRUE)
head(foldchanges)
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
organism = "hsa", # supported organisms listed below
nPerm = 1000, # default number permutations
minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
pvalueCutoff = 0.05, # padj cutoff value
verbose = FALSE)
library(clusterProfiler)
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
organism = "hsa", # supported organisms listed below
nPerm = 1000, # default number permutations
minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
pvalueCutoff = 0.05, # padj cutoff value
verbose = FALSE)
gseaKEGG@result
View(gseaKEGG@result)
write.csv(gseaKEGG@result,file="c:/users/merve/documents/commpass-deg2/results/kegg.csv")
library(pathview)
pathview(gene.data = foldchanges,
pathway.id = "hsa04060",
species = "hsa",
limit = list(gene = 2, # value gives the max/min limit for foldchanges
cpd = 1))
gseaGO <- gseGO(geneList = foldchanges,
OrgDb = org.Hs.eg.db,
ont = 'BP',
nPerm = 1000,
minGSSize = 20,
pvalueCutoff = 0.05,
verbose = FALSE)
View(gseaGO@result)
write.csv(gseaGO@result,file="c:/users/merve/documents/commpass-deg2/results/gseago.csv")
