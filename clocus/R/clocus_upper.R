#' A clocus_upper Funtion.
#'
#' This function allows you to plot th upper plot in the example.
#' @param gene_GWAS The Gene GWAS file.
#' @param grid.show Show grid or not, default is False.
#' @param top_pos Position of top SNP you want to highlight
#' @param top_p P value of top snp
#' @param top_snp Name of top snp
#' @param color Color used to plot the points. 10 levels.
#' @keywords upper
#' @export
#' @examples
#' all_gene<-fread("gencode.v19.withdir.txt",header=T)
#' C1QTNF4<-fread("C1QTNF4.txt",header=T)
#' #LD file
#' C1QTNF4_ld<-fread("C1QTNF4.ld",header=T)
#' C1QTNF4<-C1QTNF4[order(C1QTNF4$POS)]
#' C1QTNF4_ld<-C1QTNF4_ld[order(C1QTNF4_ld$BP_B)]
#' C1QTNF4<-C1QTNF4[C1QTNF4$SNP %in% C1QTNF4_ld$SNP_B]
#' #combine the LD dataset with GWAS dataset =>gene_GWAS
#' C1QTNF4_p<-data.frame(P=C1QTNF4$P,r2=C1QTNF4_ld$R2,snp=C1QTNF4$SNP,pos=C1QTNF4$POS)
#' top_p =9.702e-11
#' top_pos=47380340
#' top_snp="rs3740688"
#' all_gene_11<-all_gene[all_gene$chr==11]
#' g1<-clocus_upper(C1QTNF4_p,grid.show = element_blank(),top_p = top_p,top_pos = top_pos,top_snp = top_snp)
#' g1


clocus_upper<-function(gene_GWAS, grid.show= element_line(color = "grey94"),top_pos,top_p,top_snp,
                       color=c("#c1d8e6","#537ea3","#d0b8d9","#a283a3","#739e6f",
                               "#b3d49b","#ffc887","#e3987f","#ffb5b5","#bf4d4d")){
  ##order the gene using r2
  library(gggenes)
  library(ggplot2)
  library(ggfittext)
  library(gplots)
  library(data.table)
  library(ggrepel)
  library(latex2exp)
  library(dplyr)
  gene_GWAS$rank<-""
  for(i in 1:10){
    gene_GWAS$rank[gene_GWAS$r2<=(i/10) & gene_GWAS$r2>=((i-1)/10)] = paste("(0.",i-1,",0.",i,"]",sep="")
  }
  gene_GWAS$rank[gene_GWAS$r2>=0 & gene_GWAS$r2<=0.1] ="[0.0,0.1]"
  gene_GWAS$rank[gene_GWAS$r2>0.9 &gene_GWAS$r2<=1]   ="(0.9,1]"
  gene_GWAS$rank<-as.factor(gene_GWAS$rank)
  gene_GWAS$rank<-factor(gene_GWAS$rank,levels = c("[0.0,0.1]","(0.1,0.2]","(0.2,0.3]","(0.3,0.4]","(0.4,0.5]",
                                                   "(0.5,0.6]","(0.6,0.7]","(0.7,0.8]","(0.8,0.9]","(0.9,1]"))

  color_match<-data.frame(locuscolor = color,
                          level_col= c("[0.0,0.1]","(0.1,0.2]","(0.2,0.3]","(0.3,0.4]","(0.4,0.5]",
                                       "(0.5,0.6]","(0.6,0.7]","(0.7,0.8]","(0.8,0.9]","(0.9,1]"))
  rank0<-gene_GWAS$rank[!duplicated(gene_GWAS$rank)]
  locusPalette<-as.character(merge(as.data.frame(rank0),color_match,by.x = "rank0",by.y="level_col",sort=F)$locuscolor)

  g1= gene_GWAS %>%
    ggplot(aes(x=pos/1000000,y=-log10(P),color=rank))+geom_point()+ scale_color_manual("LD",values =locusPalette)+
    geom_point(aes(x=top_pos/1000000,y=-log10(top_p)),fill="red3",shape = 23,color = "black",size=2.7)+
    geom_hline(yintercept = -log10(5e-08),color="snow3",linetype="dashed")+
    annotate("text",x=top_pos/1000000,y=-log10(top_p)+0.7,label=top_snp,color="black",size=4)+
    xlim(min(gene_GWAS$pos/1000000),max(gene_GWAS$pos/1000000))+ylab(TeX("-log_{10}(P)"))+
    theme(legend.position = c(0.94,0.7),panel.border = element_blank(),panel.background = element_blank(),
          panel.grid.minor = grid.show,panel.grid.major = grid.show,
          axis.line = element_line(colour = "grey34"),axis.title.x = element_blank(),
          legend.background = element_blank(),legend.box.background =element_blank(),legend.key = element_blank(),plot.title = element_text(hjust = 0.5))
  return(g1)
}
