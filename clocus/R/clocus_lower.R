#' A clocus_lower Funtion.
#'
#' This function allows you to plot th lower plot in the example.
#' @param gene_GWAS The Gene GWAS file.
#' @param gene_chr Chromosome specific gene list file
#' @param gene_everyline The maximum number of genes in every line, the default number is 6.
#' @param grid.show Show grid or not, default is False.
#' @param top_snp Name of top snp
#' @param omit_gene The genes you do not want to show in the figure. Default is a empty list.
#' @keywords lower
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
#' g2<-clocus_lower(C1QTNF4_p,all_gene_11,x_title = "chromosome 11 position (Mb)")
#' g2

clocus_lower<-function(gene_GWAS,gene_chr,gene_everyline=6,x_title,grid.show=element_blank(),omit_gene=c()){
    library(dplyr)
    library(gggenes)
    library(ggplot2)
    library(ggfittext)
    library(gplots)
    library(data.table)
    library(ggrepel)
    library(latex2exp)
  gene_list<-gene_chr[gene_chr$start>=min(gene_GWAS$pos) & gene_chr$end <=max(gene_GWAS$pos),]
  if(length(omit_gene)>0){
    gene_list<-gene_list[-which(gene_list$gene_name%in%omit_gene),]
  }
  gene_line<-nrow(gene_list)%/%gene_everyline# every line includes 6 gene(defult)
  last_line<-nrow(gene_list)-gene_everyline*gene_line
  gene_list$line<-0
  gene_list$line[order(gene_list$start)<=(gene_line+1)]<-1:(gene_line+1)
  for(i in 1:gene_line){
    gene_next_index<-which(gene_list$line==i)
    end = gene_list$end[gene_list$line==i]
    start = gene_list$start[gene_list$line==i]
    for(j in 1:(gene_everyline-1)){
      gene_len<-nchar(as.character(gene_list$gene_name[gene_next_index]))
      gene_next_index =which(gene_list$start/1000000 > (start/1000000 +gene_len*0.02) & gene_list$line ==0)
      if(length(gene_next_index)>=1){
        gene_next_index<-which(gene_list$start==min(gene_list$start[gene_next_index]))
        gene_list$line[gene_next_index]=i
        start=gene_list$start[gene_next_index]
      }
    }
  }
  gene_list$line[gene_list$line==0]=gene_line+1
  gene_list$gene_name[gene_list$direction=="+"]<-paste(gene_list$gene_name[gene_list$direction=="+"],">",sep="")
  gene_list$gene_name[gene_list$direction=="-"]<-paste("<",gene_list$gene_name[gene_list$direction=="-"],sep="")

  #randomly order the gene
  # gene_list$line<-c(rep(1:gene_line,gene_everyline),rep(gene_line+1,(nrow(gene_list)-gene_line*gene_everyline)))
  # gene_list$line<-sample(gene_list$line,nrow(gene_list))
  #
  # #keep only first n line and your target gene(default is 4)
  # gene_list$line[gene_list$gene_name %in% target_gene]<-sample(1:line_show,length(target_gene))
  # gene_list<-gene_list[gene_list$line<=line_show]
  # #check whether your target gene exists
  gene_list$line<-format(as.character(gene_list$line),width=5)
  #####

  g2<-ggplot(gene_list, aes(xmin = start/1000000, xmax = end/1000000, y =line,label = gene_name)) +
    geom_gene_arrow(fill="#c1d8e6",color="#c1d8e6",arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"),arrow_body_height=unit(2,"mm")) +
    ylab("")+xlab(x_title)+xlim(min(gene_GWAS$pos/1000000),max(gene_GWAS$pos/1000000))+
    theme(axis.ticks =element_line(color = "white"),axis.text = element_text(color="white"),legend.position ="None" ,
          panel.background = element_blank(),panel.grid.major = grid.show,panel.grid.minor =grid.show)+
    geom_text_repel(aes(x=gene_list$start/1000000,y=gene_list$line,label = gene_list$gene_name),
                    show.legend=F,size = 3,nudge_y=0.4,direction = "x",arrow=NULL, segment.size = 0,hjust=0.5)

  return(g2)
}
