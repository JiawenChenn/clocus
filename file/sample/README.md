# Step-by-step guide

## Step 1: install ```clocus```
```
install_github("JiawenChenn/clocus/clocus")
library(clocus)
```
## Step 2: Download sample data of gene C1QTNF4.
Download all th files in ``clocus\file``.
1. gencode.v19.withdir.txt
2. C1QTNF4.ld
3. C1QTNF4.txt
4. C1QTNF4_snp.txt

## Step 3: Read all files into R

```
library(data.table)
all_gene<-fread("gencode.v19.withdir.txt",header=T)
```
target GENE is C1QTNF4
```
C1QTNF4<-fread("C1QTNF4.txt",header=T)
```
LD file
```
C1QTNF4_ld<-fread("C1QTNF4.ld",header=T)
C1QTNF4<-C1QTNF4[order(C1QTNF4$POS)]
C1QTNF4_ld<-C1QTNF4_ld[order(C1QTNF4_ld$BP_B)]
C1QTNF4<-C1QTNF4[C1QTNF4$SNP %in% C1QTNF4_ld$SNP_B]
```
Combine the LD dataset with GWAS dataset =>gene_GWAS
```
C1QTNF4_p<-data.frame(P=C1QTNF4$P,r2=C1QTNF4_ld$R2,snp=C1QTNF4$SNP,pos=C1QTNF4$POS)
top_p =9.702e-11
top_pos=47380340
top_snp="rs3740688"
all_gene_11<-all_gene[all_gene$chr==11]
```
## Step 4: Generate upper part.
```
g1<-clocus_upper(C1QTNF4_p,grid.show = element_blank(),top_p = top_p,top_pos = top_pos,top_snp = top_snp)
```
## Step 5: Generate lower part.
```
g2<-clocus_lower(C1QTNF4_p,all_gene_11,x_title = "chromosome 11 position (Mb)")
```
## Step 5: Combine two figures

```
grid.arrange(g1,g2,nrow = 2)
```
