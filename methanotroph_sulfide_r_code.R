#RNA-seq analysis of sulfide toxicity in Methylococcus.capsulatus
#Auther: Sichong Pei
#packages used
library(limma)
library(edgeR)
library(rtracklayer)
library(Glimma)
library(magrittr)
library(gplots)
library(ggplot2)
library(ggfortify)
library(factoextra)
library(reshape2)
library(ggpubr)
library(WGCNA)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(readxl)
setwd("C:/Users/Sichong Pei/Box Sync/Methane project/rnaseq")

#Raw featurecount and headline file-------------------------------------------------------
targets <- readTargets("Targets_Final.txt") 
targets$Group <- factor(targets$Group,levels=c("WT","Shock","P01","P05","P075"))
targets$is_h2s<-factor(rep(c(0,1,1,1,1),each=3,levels=c(0,1)))
targets$col <- as.numeric(targets$Group)
#Reorder by number
targets<-targets[order(targets$Group),]

#Read fate process Feature count matrix import------------------------------------------------------------
ReadFateBWA <- dplyr::select(targets, Assigned:Unassigned_NoFeatures)
ReadFateBWA <- ReadFateBWA/targets$Total *100
ReadFateBWA$Sample <- targets$Sample
ReadFateBWA$Group<-targets$Group
print(paste0('Does total count sum to 100?: ',sum(rowSums(ReadFateBWA[,1:4]))/dim(ReadFateBWA)[1]==100))
summary(ReadFateBWA)
# Melt to long dataframe
df0=melt(ReadFateBWA,variable.name='Fate',value.name='pct')
df0$Fate <- factor(df0$Fate, levels = unique(df0$Fate), ordered = TRUE)
df0$Group <- factor(df0$Group, levels = unique(df0$Group), ordered = TRUE)
df0$Sample<- factor(df0$Sample, levels = paste0('SP',c(1:15)), ordered = TRUE)
df0$Triplicate<-factor(rep(rep(c(1,2,3),each=4),5),levels=c(1,2,3))

ggplot(data=df0, aes(x=Sample, y=pct, fill=Fate)) +
  geom_bar(stat="identity", position='stack') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Read Fates (Alignment: BWA-MEM)") + ylab("Read%")#+facet_grid(~Group)
#BAM stat
bam_stat=data.frame(readxl::read_excel('csv/bam_stat.xlsx',sheet = 1))
ggplot(data=bam_stat,aes(x=Sample,y=coverage))+
  geom_bar(stat='identity')

#Read featurecount file----------------------------------
file_list=paste0(targets$Basename,'_fc.txt')
y <- readDGE(file_list, path = "featureCounts/", columns = c(1,7), 
             labels = targets$Sample, comment.char = "#", header = TRUE)

#all.equal(targets$Assigned, y$samples$lib.size) #<=make sure size of reads match.
#median(rowSums(y$counts))  #median of gene counts throughout all samples


#Reorganize gff file----------------------------------------------
#Parse tag name
mca_gff=readGFF('GCF_000008325.1_ASM832v1_genomic.gff',v=3)
old_gff=mca_gff[mca_gff$type=='gene',]
oldtag_df=data.frame(geneID=old_gff$locus_tag,oldgeneID=old_gff$old_locus_tag)
rm(old_gff)

mca_gff=mca_gff[mca_gff$type=='CDS',]
notation_df=data.frame(ID=mca_gff$ID,geneID=mca_gff$locus_tag,gene_name=mca_gff$gene,product=mca_gff$product);rm(mca_gff)
for (i in 1:dim(notation_df)[1]){
  cds_id=notation_df$ID[i];gene_id=notation_df$geneID[i]
  #old_name=oldtag_df$oldgeneID[oldtag_df$geneID==gene_id]
  #if (length(old_name)<1){old_name=gene_id}
  #if (is.na(old_name)==TRUE){old_name=gene_id}
  rownames(y$counts)[rownames(y$counts) == cds_id] <- gene_id#old_name
}
notation_df$oldID=oldtag_df$oldgeneID[match(notation_df$geneID,oldtag_df$geneID)]

y$genes <- data.frame(GeneID = rownames(y$counts));
max_counted_gene <- apply(y$counts, 2, max)/y$samples$lib.size * 100    #Calculate maximum counted genes.
head(y$counts)

#Normalize by TMM##################################################
y<-calcNormFactors(y,method="TMM")
y$samples$group=targets$Group
y$samples$H2S_Conc=rep(c(0,0.75,0.1,0.5,0.75),each=3)
y$samples$is_H2S=rep(c(0,1,1,1,1),each=3)
y$samples$H2S_Time=rep(c(0,1,18,20,21),each=3)
write.csv(y$counts,file='csv/raw_count.csv')

#Normalized read distribution plot################################
pseudocount=cpm(y,log=TRUE,prior.count = 1)
pseudocount_long=melt(pseudocount,variable.name='Sample',value.name='TMM')
ggplot(data=pseudocount_long,aes(x=Sample,y=TMM))+
  geom_violin()


#Exact test for wt and shock 
et_wt_vs_shock_tag$gene_name=notation_df$gene_name[match(et_wt_vs_shock_tag$GeneID,notation_df$geneID)]
et_wt_vs_shock_tag$gene_name[is.na(et_wt_vs_shock_tag$gene_name)]=notation_df$oldID[match(et_wt_vs_shock_tag$GeneID[is.na(et_wt_vs_shock_tag$gene_name)],notation_df$geneID)]
et_wt_vs_shock_result=plot_volcano(et_wt_vs_shock_tag)
et_wt_vs_shock_result[[2]]
write.csv(et_wt_vs_shock_result[[1]],file='csv/wt_vs_shock_new.csv')

#Compare 0% and 0.1%
h2s_0v01=exactTest(y,pair=c('WT','P01'))
h2s_0v01_tag=data.frame(topTags(h2s_0v01,n=Inf))
h2s_0v01_tag$gene_name=notation_df$gene_name[match(h2s_0v01_tag$GeneID,notation_df$geneID)]
h2s_0v01_tag$gene_name[is.na(h2s_0v01_tag$gene_name)]=notation_df$oldID[match(h2s_0v01_tag$GeneID[is.na(h2s_0v01_tag$gene_name)],notation_df$geneID)]
h2s_0v01_result=plot_volcano(h2s_0v01_tag)
#h2s_0v01_result[[2]]
write.csv(h2s_0v01_result[[1]],file='csv/h2s_0v01_comparison.csv')

#Compare 0% and 0.5%
h2s_0v05=exactTest(y,pair=c('WT','P05'))
h2s_0v05_tag=data.frame(topTags(h2s_0v05,n=Inf))
h2s_0v05_tag$gene_name=notation_df$gene_name[match(h2s_0v05_tag$GeneID,notation_df$geneID)]
h2s_0v05_tag$gene_name[is.na(h2s_0v05_tag$gene_name)]=notation_df$oldID[match(h2s_0v05_tag$GeneID[is.na(h2s_0v05_tag$gene_name)],notation_df$geneID)]
h2s_0v05_result=plot_volcano(h2s_0v05_tag)
#h2s_0v05_result[[2]]
write.csv(h2s_0v05_result[[1]],file='csv/h2s_0v05_comparison.csv')

#Compare 0% and 0.75%
h2s_0v075=exactTest(y,pair=c('WT','P075'))
h2s_0v075_tag=data.frame(topTags(h2s_0v075,n=Inf))
h2s_0v075_tag$gene_name=notation_df$gene_name[match(h2s_0v075_tag$GeneID,notation_df$geneID)]
h2s_0v075_tag$gene_name[is.na(h2s_0v075_tag$gene_name)]=notation_df$oldID[match(h2s_0v075_tag$GeneID[is.na(h2s_0v075_tag$gene_name)],notation_df$geneID)]
h2s_0v075_result=plot_volcano(h2s_0v075_tag)
#h2s_0v075_result[[2]]
write.csv(h2s_0v075_result[[1]],file='csv/h2s_0v075_comparison.csv')


#Venn graph of 0.1% 0.5% and 0.75% DGE#####################################################
p01 <- readxl::read_excel('csv/h2s_conc_comparison.xlsx',sheet = 1)%>% filter(abs(logFC)>1) %>% pull(GeneID)#A
p05 <- readxl::read_excel('csv/h2s_conc_comparison.xlsx',sheet = 2)%>% filter(abs(logFC)>1) %>% pull(GeneID)#B
p07 <- readxl::read_excel('csv/h2s_conc_comparison.xlsx',sheet = 3)%>% filter(abs(logFC)>1) %>% pull(GeneID)#C
h2s_shock <- readxl::read_excel('csv/et_wt_vs_shock_new.xlsx',sheet=2)%>% filter(abs(logFC)>1) %>% pull(GeneID)#S

fc01 <- readxl::read_excel('csv/h2s_conc_comparison.xlsx',sheet = 1)%>% filter(abs(logFC)>1) %>% pull(logFC)#A
fc05 <- readxl::read_excel('csv/h2s_conc_comparison.xlsx',sheet = 2)%>% filter(abs(logFC)>1) %>% pull(logFC)#B
fc07 <- readxl::read_excel('csv/h2s_conc_comparison.xlsx',sheet = 3)%>% filter(abs(logFC)>1) %>% pull(logFC)#C
fc_shock <- readxl::read_excel('csv/et_wt_vs_shock_new.xlsx',sheet=2)%>% filter(abs(logFC)>1) %>% pull(logFC)#S

p01_up=p01[fc01>0];p01_down=p01[fc01<0]
p05_up=p05[fc05>0];p05_down=p05[fc05<0]
p07_up=p07[fc07>0];p07_down=p07[fc07<0]
h2s_shock_up=h2s_shock[fc_shock>0];h2s_shock_down=h2s_shock[fc_shock<0]

h2s_comparison=list('0.1%'=p01,
                    '0.5%'=p05,
                    '0.75%'=p07,
                    shock=h2s_shock)

h2s_up_comparison=list('0.1%'=p01_up,
                       '0.5%'=p05_up,
                       '0.75%'=p07_up,
                       shock=h2s_shock_up)

h2s_down_comparison=list('0.1%'=p01_down,
                         '0.5%'=p05_down,
                         '0.75%'=p07_down,
                         shock=h2s_shock_down)

library(ggvenn)
ggvenn(h2s_comparison)
ggvenn(h2s_up_comparison)
ggvenn(h2s_down_comparison)
#install.packages('venneuler')
#library(venneuler)
#install.packages('eulerr')
library(eulerr)
h2s_euler=c('0.1%'=length(p01),
            '0.5%'=length(p05),
            '0.75%'=length(p07),
            'Shock'=length(h2s_shock),
            '0.1%&0.5%'=length(intersect(p01,p05)),
            '0.1%&0.75%'=length(intersect(p01,p07)),
            '0.1%&Shock'=length(intersect(p01,h2s_shock)),
            '0.5%&0.75%'=length(intersect(p05,p07)),
            '0.5%&Shock'=length(intersect(p05,h2s_shock)),
            '0.75%&Shock'=length(intersect(p07,h2s_shock)),
            '0.1%&0.5%&0.75%'=length(intersect(intersect(p01,p05),p07)),
            '0.1%&0.5%&Shock'=length(intersect(intersect(p01,p05),h2s_shock)),
            '0.1%&0.75%&Shock'=length(intersect(intersect(p01,p07),h2s_shock)),
            '0.5%&0.75%&Shock'=length(intersect(intersect(p05,p07),h2s_shock)),
            '0.1%&0.5%&0.75%&Shock'=length(intersect(intersect(intersect(p01,p05),p07),h2s_shock)))

a=intersect(p01,p05)
b=intersect(p01,p07)
c=intersect(p01,h2s_shock)
d=intersect(p05,p07)
e=intersect(p05,h2s_shock)
f=intersect(p07,h2s_shock)
g=intersect(intersect(p01,p05),p07)
h=intersect(intersect(p01,p05),h2s_shock)
j=intersect(intersect(p01,p07),h2s_shock)
k=intersect(intersect(p05,p07),h2s_shock)
l=intersect(intersect(intersect(p01,p05),p07),h2s_shock)
write.csv(data.frame(l),'csv/l.csv')


h2s_euler_up=c(A=length(p01_up),
               B=length(p05_up),
               C=length(p07_up),
               S=length(h2s_shock_up),
               'A&B'=length(intersect(p01_up,p05_up)),
               'A&C'=length(intersect(p01_up,p07_up)),
               'A&S'=length(intersect(p01_up,h2s_shock_up)),
               'B&C'=length(intersect(p05_up,p07_up)),
               'B&S'=length(intersect(p05_up,h2s_shock_up)),
               'C&S'=length(intersect(p07_up,h2s_shock_up)),
               'A&B&C'=length(intersect(intersect(p01_up,p05_up),p07_up)),
               'A&B&S'=length(intersect(intersect(p01_up,p05_up),h2s_shock_up)),
               'A&C&S'=length(intersect(intersect(p01_up,p07_up),h2s_shock_up)),
               'B&C&S'=length(intersect(intersect(p05_up,p07_up),h2s_shock_up)),
               'A&B&C&S'=length(intersect(intersect(intersect(p01_up,p05_up),p07_up),h2s_shock_up)))

intersect(intersect(intersect(p01,p05),p07),h2s_shock)


euler_up_plt=euler(h2s_euler_up,shape='ellipse') 

euler_plt=euler(h2s_euler,
                input='union') 
plot(euler_plt,
     fills =list(fill=c(viridis::plasma(n = 4))),
     alpha=0.8,
     edges=list(lty = 1),
     quantities=TRUE)

plot(euler_up_plt,
     quantities=TRUE,
     fill='transparent')


#plot(venneuler(h2s_comparison))
#vd <- venneuler(c(A=10, B=6, "A&B"=6));plot(vd)
# install.packages('VennDiagram')
# library(VennDiagram) 
# venn.diagram(h2s_comparison,scale=TRUE,'venn_diagram.tiff')


#Volcano plot of top tags#######################################################
#Function
plot_volcano<-function(top,fc_limit=log2(2),fdr_limit=0.1){
  #Input y: top-table (matrix like)
  topdf=top;topdf$dfex='NO';topdf$if_name=NA
  if(!('adj.P.Val' %in% colnames(topdf))){topdf$adj.P.Val=topdf$FDR}
  topdf$dfex[topdf$logFC>fc_limit & topdf$adj.P.Val<fdr_limit]<-'UP'
  topdf$dfex[topdf$logFC<(-fc_limit) & topdf$adj.P.Val<fdr_limit]<-'DOWN'
  sig_id=(abs(topdf$logFC)>1.5 & topdf$adj.P.Val<0.05)
  topdf$if_name[sig_id]<-topdf$gene_name[sig_id]
  topdf$product=notation_df$product[match(topdf$Gene,notation_df$geneID)]
  
  mycolors <- c("blue", "red", "gray")
  names(mycolors) <- c("DOWN", "UP", "NO")
  p=ggplot(data=topdf,aes(x=logFC,y=-log10(adj.P.Val),col=dfex))+
    geom_point()+
    theme_minimal()+
    geom_label_repel(
      aes(label=if_name),
      size=4,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    ) +
    geom_vline(xintercept=c(-fc_limit, fc_limit), col="red")+
    geom_hline(yintercept=-log10(fdr_limit), col="red")+
    scale_colour_manual(values=mycolors)
  return(list(topdf,p))}


# Find most/least variable gene--------------------------------------------------
ps_iqr=apply(X = pseudocount,MARGIN = 1,IQR)
ps_iqr_sorted=sort(ps_iqr)
iqr_df=data.frame(iqr=ps_iqr_sorted)
iqr_df$name=notation_df$gene_name[match(rownames(iqr_df),notation_df$geneID)]
write.csv(iqr_df,file='csv/iqr_sorted.csv')


#select most variable 1000 genes
pseudocount_top1000=pseudocount[tail(names(ps_iqr_sorted),1000),]

pca_res=prcomp(t(pseudocount_top1000),scale=TRUE)
fviz_eig(pca_res,ncp=10)
fviz_pca_ind(pca_res,
             col.ind="cos2",
             repel=TRUE)

pseudocount_top1000_adapt=pseudocount[tail(names(ps_iqr_sorted),1000),c(1:3,7:15)]
pca_res_adapt=prcomp(t(pseudocount_top1000_adapt),scale=TRUE)
fviz_eig(pca_res_adapt)
fviz_pca_ind(pca_res_adapt,
             col.ind="cos2",
             repel=TRUE)

pscount_adapt=pseudocount[,c(1:3,7:15)]
pca_full_adapt=prcomp(t(pscount_adapt),scale=TRUE)
fviz_pca_ind(pca_full_adapt,
             col.ind="cos2",
             repel=TRUE)




#Read pathway from KEGG-----------------------------------------------------------
library(KEGGREST)
library(ggpmisc)
pathways.list <- keggList("pathway", "mca")
pathway.codes <- sub("path:", "", names(pathways.list))
#Get genes in each pathway  (old locus tag name)
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             pw2 = notation_df$geneID[match(pw2,notation_df$oldID)]
                             return(pw2[!is.na(pw2)])
                           }
)
head(genes.by.pathway)

#TCA cycle mca00020
tca_shock=pseudocount[genes.by.pathway$mca00020,1:6]
write.csv(tca_shock,'csv/tca_shock_count.csv')


#Enriched pathways in shock
shock_gene=readxl::read_excel('csv/et_wt_vs_shock_new.xlsx',sheet=2) %>% pull(GeneID)
pathway_code=read.csv('csv/pathway_df.csv',header = T)
shock_pathway_count=c()
for (i in 1:length(genes.by.pathway)){
  shock_pathway_count=c(shock_pathway_count,sum(shock_gene %in% genes.by.pathway[[i]]))
}
names(shock_pathway_count)=pathway_code$pathway.name
shock_pathway_count=shock_pathway_count[shock_pathway_count>5]
#shock_pathway_count=sort(shock_pathway_count,decreasing = TRUE)[1:10]
shock_path_df=data.frame(
  Pathway=names(shock_pathway_count),
  Count=shock_pathway_count
)
library(plotly)
pie=plot_ly(shock_path_df,
            labels=~Pathway,
            values= ~Count,
            type='pie',
            textposition='inside',
            textinfo='label+percent+value',
            insidetextfont=list(color = '#FFFFFF'),
            hoverinfo = 'text') %>%
  layout(title = 'Top DEGs in H2S shock',
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
pie


#Heatmap of H2S sensitive gene#################################
y_adapt=DGEList(counts=y$counts[,-(4:6)],lib.size=colSums(y$counts[,-(4:6)]),samples=y$samples[-(4:6),])
feature_heatmap=function(y_adapt,feature){
  #y_adapt:dgelist of adapting exp
  #feature: array of gene name to plot
  ps_count=cpm(y_adapt$counts,log=TRUE) %>% t() %>% scale() %>% t()
  ps_count=ps_count[feature,]
  row_annotation=data.frame(group=rep(c('0%','0.1%','0.5%','0.75%'),each=3))
  rownames(row_annotation)=colnames(ps_count)
  ann_colors = list(group = c('0%' = "#FFCC00", '0.1%' = "#FF9900",'0.5%'='#FF6600','0.75%'='#FF3300'))
  pheatmap(t(ps_count),
           border_color = "NA",
           scale='column',
           show_rownames = F,
           annotation_row = row_annotation,
           annotation_colors = ann_colors,
           cluster_cols = F,
           cluster_rows = F,
           clustering_method="complete",
           color = colorRampPalette(c("blue", "black", "red"))(50))
  
}
df1=read.csv('csv/df1_conc_sig_df.csv')
df2=read.csv('csv/df2_conc_sig_df.csv')
df3=read.csv('csv/df3_conc_sig_df.csv')
df4=read.csv('csv/df4_conc_sig_df.csv')
ddf=rbind(df1[,c(1,3,4,5,6)],df2[,c(1,4,5,6,7)],df3[,c(1,5,6,7,8)],df4[,c(1,3,4,5,6)])

trend_gene_index=unique(c(df1$X,df2$X,df3$X,df4$X))
write.csv(trend_gene_index,file='csv/h2s_sensitive_gene_index.csv')
#Heatmap of concentration-sensitive genes
h2s_gene_index=read.csv('csv/h2s_sensitive_gene_index.csv')[,2]
h2s_count=pseudocount[h2s_gene_index,-(4:6)] %>% t() %>% scale() %>% t()
pheatmap(h2s_count,
         main='H2S concentration sensitive genes',
         border_color = "black",
         show_rownames = TRUE,
         scale='row',
         #annotation_row = gene_path_id,
         #annotation_col = col_annotation,
         #annotation_colors = ann_colors,
         cluster_cols = F,
         cluster_rows = T,
         clustering_method="complete",
         color = colorRampPalette(c("blue", "black", "red"))(50))


#Metabolome heatmap visualization---------------------------------------------------------------
cyto=data.frame(read_excel('csv/04-05-2021_GC-MS-profiling.xlsx',sheet=1))
rownames(cyto)=cyto[,1]
cyto=cyto[,-1]
cyto[is.na(cyto)]=0
cyto=cyto[rowSums(cyto==0)<5,]
cyto_average=data.frame(Control=rowMeans(cyto[,1:3]),
                        Shock=rowMeans(cyto[,4:6]),
                        S01=rowMeans(cyto[,7:9]),
                        S05=rowMeans(cyto[,10:12]),
                        S07=rowMeans(cyto[,13:15]))
cyto_average=cyto_average[rowSums(cyto_average==0)==0,]
cyto_fc=data.frame(shock=cyto_average$Shock/cyto_average$Control,
                   S01=cyto_average$S01/cyto_average$Control,
                   S05=cyto_average$S05/cyto_average$Control,
                   S07=cyto_average$S07/cyto_average$Control);
rownames(cyto_fc)=rownames(cyto_average)
cyto_fc=log2(cyto_fc)

ext=data.frame(read_excel('csv/04-05-2021_GC-MS-profiling.xlsx',sheet=2))
rownames(ext)=ext[,1]
ext=ext[,-1]
ext[is.na(ext)]=0
ext=ext[rowSums(ext==0)<=3,]
ext_average=data.frame(Control=rowMeans(ext[,1:3]),
                       Shock=rowMeans(ext[,4:6]),
                       S01=rowMeans(ext[,7:9]),
                       S05=rowMeans(ext[,10:12]),
                       S07=rowMeans(ext[,13:15]))
ext_average=ext_average[rowSums(ext_average==0)==0,]
ext_fc=data.frame(shock=ext_average$Shock/ext_average$Control,
                  S01=ext_average$S01/ext_average$Control,
                  S05=ext_average$S05/ext_average$Control,
                  S07=ext_average$S07/ext_average$Control)
rownames(ext_fc)=rownames(ext_average)
ext_fc=log2(ext_fc)


paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(cyto_fc), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(cyto_fc)/paletteLength, max(cyto_fc), length.out=floor(paletteLength/2)))
hm_cyto=pheatmap(cyto_fc,
                 scale='none',
                 border_color = 'NA',
                 show_rownames = TRUE,
                 cluster_cols = FALSE,
                 cluster_rows=TRUE,
                 color = myColor,
                 breaks= myBreaks)

cyto_average_avg=data.frame(compound=rownames(cyto),conc=rowMeans(cyto_average))
cyto_average_df=cyto_average_avg[hm_cyto$tree_row$order,]
cyto_average_df$compound=factor(cyto_average_df$compound,levels=rev(cyto_average_df$compound))
ggplot(data=cyto_average_df,aes(x=factor(compound),y=conc))+
  geom_bar(stat='identity',aes(fill=conc))+
  scale_y_log10()+
  coord_flip()+
  theme_bw()


paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(ext_fc), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(ext_fc)/paletteLength, max(ext_fc), length.out=floor(paletteLength/2)))
hm_ext=pheatmap(ext_fc,
                scale='none',
                border_color = 'NA',
                show_rownames = TRUE,
                cluster_cols = FALSE,
                cluster_rows=TRUE,
                color = myColor,
                breaks=myBreaks)

ext_average_avg=data.frame(compound=rownames(ext),conc=rowMeans(ext_average))
ext_average_df=ext_average_avg[hm_ext$tree_row$order,]
ext_average_df$compound=factor(ext_average_df$compound,levels=rev(ext_average_df$compound))
ggplot(data=ext_average_df,aes(x=factor(compound),y=conc))+
  geom_bar(stat='identity',aes(fill=conc))+
  scale_y_log10()+
  coord_flip()+
  theme_bw()


setwd("C:/Users/Sichong Pei/Box Sync/Methane project/rnaseq")
go=data.frame(readxl::read_excel('csv/Enriched_control_vs_shock.xlsx',sheet = 1))
go_match=go[,4]
n=length(go_match)

n_genes=apply(as.matrix(go$Matches,ncol=1),
              1,
              function(x) length(unlist(strsplit(x,'//'))))
go$number=n_genes
go=go[,-4]

ps_control=apply(pseudocount[,1:3],1,mean)
ps_shock=apply(pseudocount[,4:6],1,mean)
ps_data=data.frame(control=ps_control,shock=ps_shock)

#transcript abundance genomic mapping
plot_range=1:1000
ggplot(data=ps_data[plot_range,],aes(x=plot_range))+
  geom_line(aes(y=control))+
  geom_line(aes(y=shock),color='Red')+
  labs(x='Loci number',y='Normalized count')
