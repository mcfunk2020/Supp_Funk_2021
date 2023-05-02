---
title: "20220412_CapC16_MouseStrains_postChicago"
author: "A. Feldmann"
date: "27 January 2020"
output: html_document
editor_options: 
  chunk_output_type: console
---

```{r setup, include=FALSE}

base='C:/Users/angelika/OneDrive/Documents/bioinfo/'
scorecut=5
readcut=1
mode='score'
workdir=paste0(base,"analysis/CaptureC/MyData/CapC20/postchicago_fullnewpipeline_cutoff",scorecut,"scores_",readcut,"reads_mode.",mode)
dir.create(workdir)

source(paste0(base,'src/Rscripts/functions.r'))
setwd(workdir)

library(devtools)
library(argparser)
library(Chicago)
library(Hmisc)
library(MASS)
library(Delaporte)
require(Delaporte)
library(GenomicRanges)
library(Hmisc)
library(pheatmap)
library(RColorBrewer)

##Make hind file:
# hind=read.delim('../designDir2/DpnII.mm10.rmap',stringsAsFactors=FALSE,header=FALSE)
# hind.=bed2gr(hind[,1:3])
# hind.$id=hind$V4
# hind=hind.
# rm(hind.)
# saveRDS(hind,'../designDir2/DpnII.mm10.rmap.rds')
hind=readRDS('../designDir2/DpnII.mm10.rmap.rds')

ints=read.delim('ints_all.txt')

# ##read in genes table and convert to GRanges object (just once):
# # genes=read.csv(paste0(base,'analysis/MouseGenomes/rnaseq/masked_SNP.SV/ResultTables/refSeq.geneBodyUnion.ignore.strand.TRUE.inter.feature.FALSE.Peaks.Castaneus_vs_Bl6_annotated_q0.005_fe1.5_expr_annotated.csv'))
# # gr=bed2gr(genes[,c('Chr','Start','End','Strand')])
# # values(gr)=genes[,7:ncol(genes)]
# # names(values(gr))[1]='refseq'
# # names(gr)=gr$refseq
# # genes=gr
# # saveRDS(genes, paste0(base,'analysis/MouseGenomes/rnaseq/masked_SNP.SV/ResultTables/refSeq.geneBodyUnion.ignore.strand.TRUE.inter.feature.FALSE.Peaks.Castaneus_vs_Bl6_annotated_q0.005_fe1.5_expr_annotated.rds'))
# genes=readRDS(paste0(base,'analysis/MouseGenomes/rnaseq/masked_SNP.SV/ResultTables/refSeq.geneBodyUnion.ignore.strand.TRUE.inter.feature.FALSE.Peaks.Castaneus_vs_Bl6_annotated_q0.005_fe1.5_expr_annotated.rds'))

########################################################################
##baited genes, add dpn_id from baitmap
########################################################################

# ##just once:

# baitmap=read.delim('../designDir2/My.CapC16.20.baitmap',header=FALSE,stringsAsFactors=FALSE)
# gr=bed2gr(baitmap[,1:3])
# gr$dpn_id=baitmap$V4
# gr$genenames=baitmap$V5
# saveRDS(gr,'../designDir2/baited_genes.rds')

baited_genes=readRDS('../designDir2/baited_genes.rds')
names(baited_genes)=baited_genes$genenames

##grl should contain atac data received from Maja Funk 9boutros' lab)
atac=read.csv(paste0(base,'analysis/Organoids_Boutros_Funk/atac/annotated_differential_peaks_edgeR_exvivo.csv'))
atacOld=atac[atac$log2FoldChange>= log2(1.3) & atac$padj<0.1,]
atacYoung=atac[atac$log2FoldChange<= -log2(1.3) & atac$padj<0.1,]

refseq=readRDS(paste0(base,'annotations/Mmusculus.UCSC.mm10.RefSeq.20150821.transcripts.rds'))
refseqtss=resize(getTSSgr(refseq),1001,fix='center')

grl=GRangesList(bed2gr(atac[,c('seqnames','start','end')]), bed2gr(atacOld[,c('seqnames','start','end')]), bed2gr(atacYoung[,c('seqnames','start','end')]), stripGR(refseqtss))
  

names(grl)=c('ATAC_all','ATAC_old','ATAC_young','RefseqTSS')


oldAtacPeak=GRanges('chr16',IRanges(10471865, 10472410))
ciitaTSS1=GRanges('chr16',IRanges(10480059,10480059))

```


##Create list with all interactions:

```{r}

filedir='../results'

types=list.files(pattern='cd',path=filedir)
types=strsplit_string(types,s1='cd_')
types

L=list()
for (type in types){
  load(paste('../results/cd',type,sep='_'))
  L=c(L,list(cd@x))
}
names(L)=types

ids=unique(L[[1]]$baitID)
length(ids)==length(baited_genes)

```

##Create a list with all replicate interactions:

```{r}

types=list.files(pattern='cd',path=filedir)
types
types=strsplit_string(types[grep('rep',types)],s1='cd_')
types

LL=list()
for (type in types){
  load(paste(filedir,'/cd_',type,sep=''))
  a=cd@x
  # a$intID=paste(a$baitID,a$otherEndID,sep=';')
  LL=c(LL,list(a))
}

names(LL)=types

```

##################################################
##################################################
##1. Some general stats
##################################################
##################################################

```{r}
setwd(workdir)

library('BSgenome')
library(BSgenome.Mmusculus.UCSC.mm10.masked)
library(GenomicRanges)

##create tables with interesting parameters and append to a list, then plot this list
L_interesting_params=list()
for(sample in names(L)){
  print(sample)
  L_interesting_params=c(L_interesting_params,list(getInterestingParams(sample=sample,hind=hind,baited_genes=baited_genes)))
}
names(L_interesting_params)=names(L)

plotInterestingParams(L_interesting_params)


```


######################################################
######################################################
##2. Initial visualization of the data               
######################################################
######################################################

############################################
##Plot line plots of random selected baits:
############################################

```{r}


pdf(paste('Lineplots_intervalsIncluded.pdf',sep=''), width=21, height=10,onefile = TRUE)

for (id in ids){
  
  par(mfrow=c(2,3))
  print(id)
  name=paste(baited_genes[baited_genes$dpn_id==id,]$genenames,baited_genes[baited_genes$dpn_id==id,]$inttype,baited_genes[baited_genes$dpn_id==id,]$genetype,sep=';')
  
  colintervals=c('black','purple','darkgreen','orange')
  col=rep(c('blue','lightblue','green','darkred','red','orange'))
    
 
  k=75
  zoom=50000
  ylim=c(0,400)
  plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name=name, col=col,hind=hind,intervals=grl,colintervals = colintervals)
  
  k=101
  zoom=100000
  ylim=c(0,200)
  plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name=name, col=col,hind=hind,intervals=grl,colintervals = colintervals)
  
  zoom=200000
  k=155
  ylim=c(0,100)
  plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name=name, col=col,hind=hind,intervals=grl,colintervals = colintervals)
   
  zoom=500000
  k=301
  ylim=c(0,50)
  plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name=name, col=col,hind=hind,intervals=grl,colintervals = colintervals)

  zoom=750000
  k=355
  ylim=c(0,50)
  plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name=name, col=col,hind=hind,intervals=grl,colintervals = colintervals)

  
  zoom=1000000
  k=401
  ylim=c(0,20)
  plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name=name, col=col,hind=hind,intervals=grl,colintervals = colintervals)
  
  # ylim=c(0,5)
  # zoom=10000000
  # k=350
  # plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name=name,intervals=grl[c(1:4,9:13)],colintervals = colintervals, col=col,hind=hind)
  #  
}
dev.off()


##Keep specific intervals for Ciita:

pdf('Lineplots_Ciita_fixedIntervals_noPeaks.pdf')

colintervals=c('black','purple','darkgreen','orange')
  col=rep(c('blue','lightblue','green','darkred','red','orange'))
   
par(mfrow=c(3,2))

for (id in baited_genes[grep('Ciita',baited_genes$genenames)]$dpn_id){
  
   name=paste(baited_genes[baited_genes$dpn_id==id,]$genenames,baited_genes[baited_genes$dpn_id==id,]$inttype,baited_genes[baited_genes$dpn_id==id,]$genetype,sep=';')

  k=121
  ##chr16:10387972−10587972
  xlim=c(10387972,10587972)
  zoom=100000
  ylim=c(0,200)
  plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name=name, col=col,hind=hind,xlim=xlim)
  
  ##chr16:9987972−10987972
  xlim=c(9987972,10987972)
  k=351
  zoom=500000
  ylim=c(0,50)
  plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name=name, col=col,hind=hind,xlim=xlim)
  
}

dev.off()

```


##5.4.2023: Interaction with the ATAC peak of interest:

```{r}

con=as.data.frame(L[['con']])
ifn=as.data.frame(L[['IFNg']])

##################################################################
##Now identify ATAC peaks within 20kb of different TSS-Dpnfrags:
##################################################################

##resize baited genes
# b=bed2gr(baited_genes[,1:3])
b=baited_genes
b20=baited_genes
# b20=bed2gr(b20[,1:3])
b20=resize(b20,width=40000,fix='center')
# values(b20)=baited_genes
b20.bed=gr2bed(b20)
names(b20.bed)=paste(names(b20.bed),'bait',sep='.')

##overlap with atac peaks
ov=as.matrix(findOverlaps(b20,grl[['ATAC_all']]))
a=grl[['ATAC_all']][ov[,2]]
values(a)=cbind(b20.bed[ov[,1],],values(b20[ov[,1]]))
##exclude all ATAC peaks overlapping TSSfrags
ov=as.matrix(findOverlaps(b,a))
a=a[-ov[,2]]

############################################
##Now check how many reads each promoter has
##Reads within peaks will be normalized
##to these reads
############################################

##in control sample
vc=vector()
for (id in baited_genes$dpn_id){
  vc=c(vc,sum(L[[1]][L[[1]]$baitID==id]$N))
}

##in IFNg sample
vi=vector()
for (id in baited_genes$dpn_id){
  vi=c(vi,sum(L[[8]][L[[8]]$baitID==id]$N))
}

baited_genes$cov_con=vc
baited_genes$cov_ifng=vi

#######################################################################################
##For these atac peaks find overlapping DpnII frags and two flanking frags on each side
##and quantify 
#######################################################################################

gr=GRanges()
for (i in 1:length(a)){
  ov=as.matrix(findOverlaps(a[i],hind))
  id=hind[ov[,2]]$id
  ids=(id[1]-2):(id[length(id)]+2)
  h=hind[ids]
  
  g=a[i]
  
  ##read counts in the control sample (norm to promoter read count and to id length):
  g$reads_con=sum(con[con$baitID==a[i]$dpn_id & con$otherEndID %in% ids,]$N)
  g$reads_con_normIds=g$reads_con/length(ids)
  g$reads_con_normIds_normCov=g$reads_con_normIds/baited_genes[baited_genes$dpn_id==a[i]$dpn_id,]$cov_con
  
  ##read counts in the IFNg sample (norm to promoter read count and to id length):
  g$reads_ifn=sum(ifn[ifn$baitID==a[i]$dpn_id & ifn$otherEndID %in% ids,]$N)
  g$reads_ifn_normIds=g$reads_ifn/length(ids)
  g$reads_ifn_normIds_normCov=g$reads_ifn_normIds/baited_genes[baited_genes$dpn_id==a[i]$dpn_id,]$cov_ifng
  
  ##score counts in the control sample (norm to promoter read count and to id length):
  g$scores_con=sum(con[con$baitID==a[i]$dpn_id & con$otherEndID %in% ids,]$score)
  g$scores_con_normIds=g$scores_con/length(ids)
  g$scores_con_normIds_normCov=g$scores_con_normIds/baited_genes[baited_genes$dpn_id==a[i]$dpn_id,]$cov_con
  
  ##score counts in the IFNg sample (norm to promoter read count and to id length):
  g$scores_ifn=sum(ifn[ifn$baitID==a[i]$dpn_id & ifn$otherEndID %in% ids,]$score)
  g$scores_ifn_normIds=g$scores_ifn/length(ids)
  g$scores_ifn_normIds_normCov=g$scores_ifn_normIds/baited_genes[baited_genes$dpn_id==a[i]$dpn_id,]$cov_ifng
  
  gr=c(gr,g)
}

##Mark old ATAC peak
gr$oldATACPeak=FALSE
ov=as.matrix(findOverlaps(gr,oldAtacPeak))
gr[ov[,1]]$oldATACPeak=TRUE

##Add distance center peak to bait peak:
gr$distance_atacPeakCenter_to_baitCenter=abs((start(gr)+end(gr))/2 - (gr$start.bait+gr$end.bait)/2)

saveRDS(gr,'ATAC_peaks_within_20kb_of_TSSfrags.rds')

#####################################################################
##Plots
#####################################################################

pdf('Read_count_quantification_within_promoter_proximal_ATAC_peaks_20kb.pdf',onefile=TRUE)

par(mfcol=c(2,3))

at=gr[gr$oldATACPeak]

xlim=c(0,max(c(gr$reads_con,gr$reads_ifn)))
plot(gr$reads_con,gr$reads_ifn,main='reads',xlim=xlim,ylim=xlim)
points(at$reads_con,at$reads_ifn,pch=20,col='red')
text(at$reads_con-2,at$reads_ifn+5,labels = at$V5,col='red')
abline(0,1,lty=5)

xlim=c(0,max(c(gr$score_con,gr$scores_ifn)))
plot(gr$scores_con,gr$scores_ifn,main='scores',xlim=xlim,ylim=xlim)
points(at$scores_con,at$scores_ifn,pch=20,col='red')
text(at$scores_con-0.5,at$scores_ifn+0.5,labels = at$V5,col='red')
abline(0,1,lty=5)

xlim=c(0,max(c(gr$reads_con_normIds,gr$reads_ifn_normIds)))
plot(gr$reads_con_normIds,gr$reads_ifn_normIds,main='reads_normIDs',xlim=xlim,ylim=xlim)
points(at$reads_con_normIds,at$reads_ifn_normIds,pch=20,col='red')
text(at$reads_con_normIds-1,at$reads_ifn_normIds+1,labels = at$V5,col='red')
abline(0,1,lty=5)

xlim=c(0,max(c(gr$score_con_normIds,gr$scores_ifn_normIds)))
plot(gr$scores_con_normIds,gr$scores_ifn_normIds,main='scores_normIDs',xlim=xlim,ylim=xlim)
points(at$scores_con_normIds,at$scores_ifn_normIds,pch=20,col='red')
text(at$scores_con_normIds-0.1,at$scores_ifn_normIds+0.2,labels = at$V5,col='red')
abline(0,1,lty=5)

xlim=c(0,max(c(gr$reads_con_normIds_normCov,gr$reads_ifn_normIds_normCov)))
plot(gr$reads_con_normIds_normCov,gr$reads_ifn_normIds_normCov,main='reads_normIDs_normCov',xlim=xlim,ylim=xlim)
points(at$reads_con_normIds_normCov,at$reads_ifn_normIds_normCov,pch=20,col='red')
text(at$reads_con_normIds_normCov-0.0001,at$reads_ifn_normIds_normCov+0.0002,labels = at$V5,col='red')
abline(0,1,lty=5)

xlim=c(0,max(c(gr$score_con_normIds_normCov,gr$scores_ifn_normIds_normCov)))
plot(gr$scores_con_normIds_normCov,gr$scores_ifn_normIds_normCov,main='scores_normIDs_normCov',xlim=xlim,ylim=xlim)
points(at$scores_con_normIds_normCov,at$scores_ifn_normIds_normCov,pch=20,col='red')
text(at$scores_con_normIds_normCov-0.00005,at$scores_ifn_normIds_normCov+0.00005,labels = at$V5,col='red')
abline(0,1,lty=5)

##relative to the distance
par(mfrow=c(2,2))

plot(gr$distance_atacPeakCenter_to_baitCenter,gr$reads_con_normIds_normCov,main='reads_normIDs_normCov con')
points(at$distance_atacPeakCenter_to_baitCenter,at$reads_con_normIds_normCov,pch=20,col='red')
text(at$distance_atacPeakCenter_to_baitCenter-0.0001,at$reads_con_normIds_normCov+0.0005,labels = at$V5,col='red')

##relative to the distance
plot(gr$distance_atacPeakCenter_to_baitCenter,gr$reads_ifn_normIds_normCov,main='reads_normIDs_normCov IFN')
points(at$distance_atacPeakCenter_to_baitCenter,at$reads_ifn_normIds_normCov,pch=20,col='red')
text(at$distance_atacPeakCenter_to_baitCenter+1500,at$reads_ifn_normIds_normCov+0.0001,labels = at$V5,col='red')

plot(gr$distance_atacPeakCenter_to_baitCenter,gr$scores_con_normIds_normCov,main='scores_normIDs_normCov con')
points(at$distance_atacPeakCenter_to_baitCenter,at$scores_con_normIds_normCov,pch=20,col='red')
text(at$distance_atacPeakCenter_to_baitCenter-0.0001,at$scores_con_normIds_normCov+0.00005,labels = at$V5,col='red')

plot(gr$distance_atacPeakCenter_to_baitCenter,gr$scores_ifn_normIds_normCov,main='scores_normIDs_normCov IFN')
points(at$distance_atacPeakCenter_to_baitCenter,at$scores_ifn_normIds_normCov,pch=20,col='red')
text(at$distance_atacPeakCenter_to_baitCenter-0.0001,at$scores_ifn_normIds_normCov+0.00005,labels = at$V5,col='red')

dev.off()



##Nice scatterplot:

names(gr)=1:length(gr)
dd=as.data.frame(gr[gr$oldATACPeak])
d=as.data.frame(gr)

##relative to the distance
plot(gr$distance_atacPeakCenter_to_baitCenter,gr$reads_ifn_normIds_normCov,main='reads_normIDs_normCov IFN')
points(at$distance_atacPeakCenter_to_baitCenter,at$reads_ifn_normIds_normCov,pch=20,col='red')
text(at$distance_atacPeakCenter_to_baitCenter+1500,at$reads_ifn_normIds_normCov+0.0001,labels = at$V5,col='red')

plot(gr$reads_con,gr$reads_ifn,main='reads',xlim=xlim,ylim=xlim)
points(at$reads_con,at$reads_ifn,pch=20,col='red')
text(at$reads_con-2,at$reads_ifn+5,labels = at$V5,col='red')
abline(0,1,lty=5)


library(ggplot2)

p1=ggplot()+
  theme_bw()+ theme(axis.line.x = element_line(color="black"),  axis.line.y = element_line(color="black"))+
  theme(legend.title=element_blank(),legend.position=c(0.85,0.2),legend.key = element_blank(),axis.line=element_line(colour="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
  
  geom_point(data=d,aes(distance_atacPeakCenter_to_baitCenter , reads_con_normIds_normCov),colour="black",alpha=0.3,size=I(3))+
  geom_point(data=dd,aes(distance_atacPeakCenter_to_baitCenter , reads_con_normIds_normCov),colour="red",alpha=1,size=I(3))+
  
  xlab('Distance bait to peak center [bp]')+ylab('Normalized Capture-C interaction frequency Con')+
  ylim(c(0,0.005))+
  
  theme(legend.position="none",text = element_text(size=12))+
  
  ggtitle('Control Capture-C \n red: Ciita-interacting ATAC peak \n left-to-right: Tss1, TSS2, TSS3')


p2=ggplot()+
  theme_bw()+ theme(axis.line.x = element_line(color="black"),  axis.line.y = element_line(color="black"))+
  theme(legend.title=element_blank(),legend.position=c(0.85,0.2),legend.key = element_blank(),axis.line=element_line(colour="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
  
  geom_point(data=d,aes(distance_atacPeakCenter_to_baitCenter , reads_ifn_normIds_normCov),colour="black",alpha=0.3,size=I(3))+
  geom_point(data=dd,aes(distance_atacPeakCenter_to_baitCenter , reads_ifn_normIds_normCov),colour="red",alpha=1,size=I(3))+
  
  xlab('Distance bait to peak center [bp]')+ylab('Normalized Capture-C interaction frequency IFNg')+
  ylim(c(0,0.005))+
  
  theme(legend.position="none",text = element_text(size=12))+
  
  ggtitle('IFNg Capture-C \n red: Ciita-interacting ATAC peak \n left-to-right: Tss1, TSS2, TSS3')


pdf('ggplot_scatters_distanceVsATACpeakInteractions.pdf',width=9,height=4.5)
multiplot(p1,p2,cols=2)
dev.off()


#############################################################
##Nice scatterplot marking all Ciita_TSS1 closeby ATAC peaks:
#############################################################

names(gr)=1:length(gr)
dd=as.data.frame(gr[gr$genenames=='Ciita_1'])
d=as.data.frame(gr)

library(ggplot2)

p1=ggplot()+
  theme_bw()+ theme(axis.line.x = element_line(color="black"),  axis.line.y = element_line(color="black"))+
  theme(legend.title=element_blank(),legend.position=c(0.85,0.2),legend.key = element_blank(),axis.line=element_line(colour="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
  
  geom_point(data=d,aes(distance_atacPeakCenter_to_baitCenter , reads_con_normIds_normCov),colour="black",alpha=0.3,size=I(3))+
  geom_point(data=dd,aes(distance_atacPeakCenter_to_baitCenter , reads_con_normIds_normCov),colour="green",alpha=0.5,size=I(3))+
  
  xlab('Distance bait to peak center [bp]')+ylab('Normalized Capture-C interaction frequency Con')+
  ylim(c(0,0.008))+
  
  theme(legend.position="none",text = element_text(size=12))+
  
  ggtitle('Control Capture-C \n green: Ciita-TSS1 proximal ATAC peaks')


p2=ggplot()+
  theme_bw()+ theme(axis.line.x = element_line(color="black"),  axis.line.y = element_line(color="black"))+
  theme(legend.title=element_blank(),legend.position=c(0.85,0.2),legend.key = element_blank(),axis.line=element_line(colour="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
  
  geom_point(data=d,aes(distance_atacPeakCenter_to_baitCenter , reads_ifn_normIds_normCov),colour="black",alpha=0.3,size=I(3))+
  geom_point(data=dd,aes(distance_atacPeakCenter_to_baitCenter , reads_ifn_normIds_normCov),colour="green",alpha=0.5,size=I(3))+
  
  xlab('Distance bait to peak center [bp]')+ylab('Normalized Capture-C interaction frequency IFNg')+
  ylim(c(0,0.005))+
  
  theme(legend.position="none",text = element_text(size=12))+
  
  ggtitle('IFNg Capture-C \n green: Ciita-TSS1 proximal ATAC peaks ')


########################################################
##Nice scatterplots only looking at genes with a similar
##promoter coverage as Ciita
########################################################

gn=baited_genes[baited_genes$cov_con>= 0.8*baited_genes[baited_genes$genenames=='Ciita_1']$cov_con & baited_genes$cov_con <= 1.2*baited_genes[baited_genes$genenames=='Ciita_1']$cov_con ]$genenames

names(gr)=1:length(gr)
dd=as.data.frame(gr[gr$genenames=='Ciita_1'])
d=as.data.frame(gr[gr$genenames %in% gn,])

library(ggplot2)

p3=ggplot()+
  theme_bw()+ theme(axis.line.x = element_line(color="black"),  axis.line.y = element_line(color="black"))+
  theme(legend.title=element_blank(),legend.position=c(0.85,0.2),legend.key = element_blank(),axis.line=element_line(colour="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
  
  geom_point(data=d,aes(distance_atacPeakCenter_to_baitCenter , reads_con_normIds_normCov),colour="black",alpha=0.3,size=I(3))+
  geom_point(data=dd,aes(distance_atacPeakCenter_to_baitCenter , reads_con_normIds_normCov),colour="green",alpha=0.5,size=I(3))+
  
  xlab('Distance bait to peak center [bp]')+ylab('Normalized Capture-C interaction frequency Con')+
  ylim(c(0,0.008))+
  
  theme(legend.position="none",text = element_text(size=12))+
  
  ggtitle('Control Capture-C genes with  promCov Ciita_1+/-20% \n green: Ciita-TSS1 proximal ATAC peaks')


p4=ggplot()+
  theme_bw()+ theme(axis.line.x = element_line(color="black"),  axis.line.y = element_line(color="black"))+
  theme(legend.title=element_blank(),legend.position=c(0.85,0.2),legend.key = element_blank(),axis.line=element_line(colour="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
  
  geom_point(data=d,aes(distance_atacPeakCenter_to_baitCenter , reads_ifn_normIds_normCov),colour="black",alpha=0.3,size=I(3))+
  geom_point(data=dd,aes(distance_atacPeakCenter_to_baitCenter , reads_ifn_normIds_normCov),colour="green",alpha=0.5,size=I(3))+
  
  xlab('Distance bait to peak center [bp]')+ylab('Normalized Capture-C interaction frequency IFNg')+
  ylim(c(0,0.005))+
  
  theme(legend.position="none",text = element_text(size=12))+
  
  ggtitle('IFNg Capture-C genes with promCov Ciita_1+/-20% \n green: Ciita-TSS1 proximal ATAC peaks ')



pdf('ggplot_scatters_distanceVsATACpeakInteractions_Ciita_1_ATAC_peaks_marked.pdf',width=12,height=9)
multiplot(p1,p3,p2,p4,cols=2)
dev.off()




```


##14.4.2023: Interaction with the ATAC peak of interest for individual replicates:

```{r}

##################################################################
##Identify ATAC peaks within 20kb of different TSS-Dpnfrags:
##################################################################

##resize baited genes
b=baited_genes
b20=baited_genes
b20=resize(b20,width=40000,fix='center')
b20.bed=gr2bed(b20)
names(b20.bed)=paste(names(b20.bed),'bait',sep='.')

##overlap with atac peaks
ov=as.matrix(findOverlaps(b20,grl[['ATAC_all']]))
a=grl[['ATAC_all']][ov[,2]]
values(a)=cbind(b20.bed[ov[,1],],values(b20[ov[,1]]))
ov=as.matrix(findOverlaps(b,a))
a=a[-ov[,2]]


for (j in 1:7){

  con=as.data.frame(L[[j]])
  ifn=as.data.frame(L[[(j+7)]])
  
  ############################################
  ##Now check how many reads each promoter has
  ##Reads within peaks will be normalized
  ##to these reads
  ############################################
  
  ##in control sample
  vc=vector()
  for (id in baited_genes$dpn_id){
    vc=c(vc,sum(L[[j]][L[[j]]$baitID==id]$N))
  }
  
  ##in IFNg sample
  vi=vector()
  for (id in baited_genes$dpn_id){
    vi=c(vi,sum(L[[(j+7)]][L[[(j+7)]]$baitID==id]$N))
  }
  
  baited_genes$cov_con=vc
  baited_genes$cov_ifng=vi
  
  #######################################################################################
  ##For these atac peaks find overlapping DpnII frags and two flanking frags on each side
  ##and quantify 
  #######################################################################################
  
  gr=GRanges()
  for (i in 1:length(a)){
    ov=as.matrix(findOverlaps(a[i],hind))
    id=hind[ov[,2]]$id
    ids=(id[1]-2):(id[length(id)]+2)
    h=hind[ids]
    
    g=a[i]
    
    ##read counts in the control sample (norm to promoter read count and to id length):
    g$reads_con=sum(con[con$baitID==a[i]$dpn_id & con$otherEndID %in% ids,]$N)
    g$reads_con_normIds=g$reads_con/length(ids)
    g$reads_con_normIds_normCov=g$reads_con_normIds/baited_genes[baited_genes$dpn_id==a[i]$dpn_id,]$cov_con
    
    ##read counts in the IFNg sample (norm to promoter read count and to id length):
    g$reads_ifn=sum(ifn[ifn$baitID==a[i]$dpn_id & ifn$otherEndID %in% ids,]$N)
    g$reads_ifn_normIds=g$reads_ifn/length(ids)
    g$reads_ifn_normIds_normCov=g$reads_ifn_normIds/baited_genes[baited_genes$dpn_id==a[i]$dpn_id,]$cov_ifng
    
    ##score counts in the control sample (norm to promoter read count and to id length):
    g$scores_con=sum(con[con$baitID==a[i]$dpn_id & con$otherEndID %in% ids,]$score)
    g$scores_con_normIds=g$scores_con/length(ids)
    g$scores_con_normIds_normCov=g$scores_con_normIds/baited_genes[baited_genes$dpn_id==a[i]$dpn_id,]$cov_con
    
    ##score counts in the IFNg sample (norm to promoter read count and to id length):
    g$scores_ifn=sum(ifn[ifn$baitID==a[i]$dpn_id & ifn$otherEndID %in% ids,]$score)
    g$scores_ifn_normIds=g$scores_ifn/length(ids)
    g$scores_ifn_normIds_normCov=g$scores_ifn_normIds/baited_genes[baited_genes$dpn_id==a[i]$dpn_id,]$cov_ifng
    
    gr=c(gr,g)
  }
  
  ##Mark old ATAC peak
  gr$oldATACPeak=FALSE
  ov=as.matrix(findOverlaps(gr,oldAtacPeak))
  gr[ov[,1]]$oldATACPeak=TRUE
  
  ##Add distance center peak to bait peak:
  gr$distance_atacPeakCenter_to_baitCenter=abs((start(gr)+end(gr))/2 - (gr$start.bait+gr$end.bait)/2)
  
  saveRDS(gr,paste0('ATAC_peaks_within_20kb_of_TSSfrags.',names(L)[j],'.rds'))
  
  #####################################################################
  ##Plots
  #####################################################################
  
  # pdf(paste0('Read_count_quant_within_promprox_ATAC_peaks_20kb.',names(L)[j],'.pdf'),onefile=TRUE)
  # 
  # par(mfcol=c(2,3))
  # 
  # at=gr[gr$oldATACPeak]
  # 
  # xlim=c(0,max(c(gr$reads_con,gr$reads_ifn)))
  # plot(gr$reads_con,gr$reads_ifn,main='reads',xlim=xlim,ylim=xlim)
  # points(at$reads_con,at$reads_ifn,pch=20,col='red')
  # text(at$reads_con-2,at$reads_ifn+5,labels = at$V5,col='red')
  # abline(0,1,lty=5)
  # 
  # xlim=c(0,max(c(gr$score_con,gr$scores_ifn)))
  # plot(gr$scores_con,gr$scores_ifn,main='scores',xlim=xlim,ylim=xlim)
  # points(at$scores_con,at$scores_ifn,pch=20,col='red')
  # text(at$scores_con-0.5,at$scores_ifn+0.5,labels = at$V5,col='red')
  # abline(0,1,lty=5)
  # 
  # xlim=c(0,max(c(gr$reads_con_normIds,gr$reads_ifn_normIds)))
  # plot(gr$reads_con_normIds,gr$reads_ifn_normIds,main='reads_normIDs',xlim=xlim,ylim=xlim)
  # points(at$reads_con_normIds,at$reads_ifn_normIds,pch=20,col='red')
  # text(at$reads_con_normIds-1,at$reads_ifn_normIds+1,labels = at$V5,col='red')
  # abline(0,1,lty=5)
  # 
  # xlim=c(0,max(c(gr$score_con_normIds,gr$scores_ifn_normIds)))
  # plot(gr$scores_con_normIds,gr$scores_ifn_normIds,main='scores_normIDs',xlim=xlim,ylim=xlim)
  # points(at$scores_con_normIds,at$scores_ifn_normIds,pch=20,col='red')
  # text(at$scores_con_normIds-0.1,at$scores_ifn_normIds+0.2,labels = at$V5,col='red')
  # abline(0,1,lty=5)
  # 
  # xlim=c(0,max(c(gr$reads_con_normIds_normCov,gr$reads_ifn_normIds_normCov)))
  # plot(gr$reads_con_normIds_normCov,gr$reads_ifn_normIds_normCov,main='reads_normIDs_normCov',xlim=xlim,ylim=xlim)
  # points(at$reads_con_normIds_normCov,at$reads_ifn_normIds_normCov,pch=20,col='red')
  # text(at$reads_con_normIds_normCov-0.0001,at$reads_ifn_normIds_normCov+0.0002,labels = at$V5,col='red')
  # abline(0,1,lty=5)
  # 
  # xlim=c(0,max(c(gr$score_con_normIds_normCov,gr$scores_ifn_normIds_normCov)))
  # plot(gr$scores_con_normIds_normCov,gr$scores_ifn_normIds_normCov,main='scores_normIDs_normCov',xlim=xlim,ylim=xlim)
  # points(at$scores_con_normIds_normCov,at$scores_ifn_normIds_normCov,pch=20,col='red')
  # text(at$scores_con_normIds_normCov-0.00005,at$scores_ifn_normIds_normCov+0.00005,labels = at$V5,col='red')
  # abline(0,1,lty=5)
  # 
  # ##relative to the distance
  # par(mfrow=c(2,2))
  # 
  # plot(gr$distance_atacPeakCenter_to_baitCenter,gr$reads_con_normIds_normCov,main='reads_normIDs_normCov con')
  # points(at$distance_atacPeakCenter_to_baitCenter,at$reads_con_normIds_normCov,pch=20,col='red')
  # text(at$distance_atacPeakCenter_to_baitCenter-0.0001,at$reads_con_normIds_normCov+0.0005,labels = at$V5,col='red')
  # 
  # ##relative to the distance
  # plot(gr$distance_atacPeakCenter_to_baitCenter,gr$reads_ifn_normIds_normCov,main='reads_normIDs_normCov IFN')
  # points(at$distance_atacPeakCenter_to_baitCenter,at$reads_ifn_normIds_normCov,pch=20,col='red')
  # text(at$distance_atacPeakCenter_to_baitCenter+1500,at$reads_ifn_normIds_normCov+0.0001,labels = at$V5,col='red')
  # 
  # plot(gr$distance_atacPeakCenter_to_baitCenter,gr$scores_con_normIds_normCov,main='scores_normIDs_normCov con')
  # points(at$distance_atacPeakCenter_to_baitCenter,at$scores_con_normIds_normCov,pch=20,col='red')
  # text(at$distance_atacPeakCenter_to_baitCenter-0.0001,at$scores_con_normIds_normCov+0.00005,labels = at$V5,col='red')
  # 
  # plot(gr$distance_atacPeakCenter_to_baitCenter,gr$scores_ifn_normIds_normCov,main='scores_normIDs_normCov IFN')
  # points(at$distance_atacPeakCenter_to_baitCenter,at$scores_ifn_normIds_normCov,pch=20,col='red')
  # text(at$distance_atacPeakCenter_to_baitCenter-0.0001,at$scores_ifn_normIds_normCov+0.00005,labels = at$V5,col='red')
  # 
  # dev.off()
  
  
  
  ##Nice scatterplot:
  
  names(gr)=1:length(gr)
  dd=as.data.frame(gr[gr$oldATACPeak])
  d=as.data.frame(gr)
  
  ##relative to the distance
  plot(gr$distance_atacPeakCenter_to_baitCenter,gr$reads_ifn_normIds_normCov,main='reads_normIDs_normCov IFN')
  points(at$distance_atacPeakCenter_to_baitCenter,at$reads_ifn_normIds_normCov,pch=20,col='red')
  text(at$distance_atacPeakCenter_to_baitCenter+1500,at$reads_ifn_normIds_normCov+0.0001,labels = at$V5,col='red')
  
  plot(gr$reads_con,gr$reads_ifn,main='reads',xlim=xlim,ylim=xlim)
  points(at$reads_con,at$reads_ifn,pch=20,col='red')
  text(at$reads_con-2,at$reads_ifn+5,labels = at$V5,col='red')
  abline(0,1,lty=5)
  
  
  library(ggplot2)
  
  p1=ggplot()+
    theme_bw()+ theme(axis.line.x = element_line(color="black"),  axis.line.y = element_line(color="black"))+
    theme(legend.title=element_blank(),legend.position=c(0.85,0.2),legend.key = element_blank(),axis.line=element_line(colour="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
    
    geom_point(data=d,aes(distance_atacPeakCenter_to_baitCenter , reads_con_normIds_normCov),colour="black",alpha=0.3,size=I(3))+
    geom_point(data=dd,aes(distance_atacPeakCenter_to_baitCenter , reads_con_normIds_normCov),colour="red",alpha=1,size=I(3))+
    
    xlab('Distance bait to peak center [bp]')+ylab('Normalized Capture-C interaction frequency Con')+
    ylim(c(0,0.009))+
    
    theme(legend.position="none",text = element_text(size=12))+
    
    ggtitle(paste(names(L)[j],'Control Capture-C \n red: Ciita-interacting ATAC peak \n left-to-right: Tss1, TSS2, TSS3'))
  
  
  p2=ggplot()+
    theme_bw()+ theme(axis.line.x = element_line(color="black"),  axis.line.y = element_line(color="black"))+
    theme(legend.title=element_blank(),legend.position=c(0.85,0.2),legend.key = element_blank(),axis.line=element_line(colour="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
    
    geom_point(data=d,aes(distance_atacPeakCenter_to_baitCenter , reads_ifn_normIds_normCov),colour="black",alpha=0.3,size=I(3))+
    geom_point(data=dd,aes(distance_atacPeakCenter_to_baitCenter , reads_ifn_normIds_normCov),colour="red",alpha=1,size=I(3))+
    
    xlab('Distance bait to peak center [bp]')+ylab('Normalized Capture-C interaction frequency IFNg')+
    ylim(c(0,0.009))+
    
    theme(legend.position="none",text = element_text(size=12))+
    
    ggtitle(paste(names(L)[j],'IFNg Capture-C \n red: Ciita-interacting ATAC peak \n left-to-right: Tss1, TSS2, TSS3'))
  
  
  pdf(paste0('ggplot_scatters_distanceVsATACpeakInteractions.',names(L)[j],'.pdf'),width=9,height=4.5)
  multiplot(p1,p2,cols=2)
  dev.off()

}

##TO DO: Calculate at the level of individual Replicates



```



############################################################
############################################################
##3. Make Data summary tables
############################################################
############################################################

! table creation requires a lot of memory !

```{r}

setwd(workdir)

##recommended, esp. if many view points:
memory.limit(size=32000)

##create table with all interactions:
##annotated with scores and reads from both replicates and pooled samples
##normalized to captured library size

# ints=getIntsBigdata(L,baited_genes,repscores=TRUE,LL=LL,scorecut=scorecut,readcut=readcut,mode=mode)
ints=getIntsBigdata(L[c('con','IFNg')],baited_genes,repscores=FALSE,scorecut=scorecut,readcut=readcut,mode=mode)

head(ints)

##annotate interactions with all your intervals
ints=annotateInts(ints=ints,grl=grl)

##look at your table
head(ints)
names(ints)

##save:
write.table(ints,'ints_all.txt',sep='\t',eol='\n',quote=FALSE)

```


##Make QC plots:

```{r}
##correlate quantification between datasets: reads
m=ints[,grep('N[.]',names(ints))]
mraw=log2(m[,-grep('downsampled',names(m))]+1)
mnorm=log2(m[,grep('downsampled',names(m))]+1)
pheatmap(mraw,show_rownames = F,filename = 'pheatmap_rawReads_byInt.png',main='log2(raw_reads+1)')
pheatmap(mnorm,show_rownames = F,filename = 'pheatmap_normReads_byInt.png',main='log2(downsampled_reads+1)')

##correlate quantification between datasets: scores
m=ints[,grep('score',names(ints))]
msc=asinh(m[,grep('rep',names(m))])
pheatmap(msc,show_rownames = F,filename = 'pheatmap_Scores_byInt.png',main='asinh(scores)')

# ##illustrate interactions in the line plots (requires ints to have a column called 'clusters_refined
# par(mfrow=c(2,1))
# ylim=c(0,40)
# k=10
# zoom=1000000
# id=5402489
# name=paste(baited_genes[baited_genes$dpn_id==id,]$genenames)
# plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name=name,intervals=grl,hind=hind,d=ints)
# xlim=c(14260000,14300000)
# plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name=name,intervals=grl,hind=hind,xlim=xlim,d=ints)

```



#########################################################
##23.04.2022
#########################################################

#####################################################################
#####################################################################
##4. Integration with K27ac peak sets
#####################################################################
#####################################################################

##Annotate K27ac peaks as connected or not to any pluripotency factor

```{r}

# ##Do the following just once, then load the k27ac_all at the beginning of the script!
# 
# range=c('Chr','Start','End')
# 
# k27ac_all$connected=FALSE
# ov=as.matrix(findOverlaps(bed2gr(k27ac_all[,range]), bed2gr(ints[,c('seqnames_otherEnd','start_otherEnd','end_otherEnd')])))
# k27ac_all[ov[,1],]$connected=TRUE
# 
# nconnected=vector()
# lfc=list()
# nam=vector()
# 
# for (g in as.character(unique(baited_genes$genenames))){
#   
#   print(g)
#   n=paste0(g,'_connected')
#   k27ac_all[,n]=FALSE
#   ov=as.matrix(findOverlaps(bed2gr(k27ac_all[,range]), bed2gr(ints[ints$genename==g,c('seqnames_otherEnd','start_otherEnd','end_otherEnd')])))
# k27ac_all[ov[,1],n]=TRUE
# 
#  print(summary(k27ac_all[,n]))
#  
#  nconnected=c(nconnected,nrow(k27ac_all[k27ac_all[,n]==TRUE,]))
#  lfc=c(lfc,list(k27ac_all[k27ac_all[,n]==TRUE,]$LFC))
#  nam=c(nam,g)    
# }
# 
# ##save once:
# write.csv(k27ac_all,paste0(base,'analysis/MouseGenomes/chip/masked_SNP.SV/ResultTables/Bl6Castaneus_k27ac_consensusPeaks_2reps_q',q,'_fe',fe,'.Peaks.Castaneus_vs_Bl6_annotated_pluripotencyGeneConnection_annotated.csv'),row.names=FALSE)

pdf('Summary_geneExpression_peakConnection.pdf',width=14,height=7)
par(mfrow=c(2,3))
barplot(baited_genes[order(baited_genes$LFC)]$LFC,names=baited_genes[order(baited_genes$LFC)]$genenames,las=3,ylim=c(-1.5,0.5),main='strain-specific gene expression',ylab='transcription LFC Bl6/Cast')
abline(h=0)
barplot(nconnected[order(baited_genes$LFC)],names=nam[order(baited_genes$LFC)],las=3,ylab='# connected K27ac peaks per promoter',ylim=c(0,12),main='total connected K27ac peaks')
abline(h=median(nconnected),lty=5)
boxplot(lfc[order(baited_genes$LFC)],names=nam[order(baited_genes$LFC)],las=3,ylab='LFC Bl6/Cast, K27ac in connecetd peaks',main='K27ac change at connected peaks')
abline(h=0,lty=2)
dev.off()

```


########################################################################
########################################################################
##Analysis of TFBSs at interactions with K27ac peaks
##Overlap with SNPs?
##Any candidate factors/SNPs that may drive differences in K27ac and interactions?
########################################################################
########################################################################

```{r}

# ##Do the following just once: 
# 
# ##load gene expression data:
# d=read.csv(paste0(base,'analysis/MouseGenomes/rnaseq/masked_SNP.SV/ResultTables/refSeq.geneBodyUnion.ignore.strand.TRUE.inter.feature.FALSE.Peaks.Castaneus_vs_Bl6_annotated_q0.005_fe1.5.csv'),stringsAsFactors=FALSE)
# 
# d$Bl6.pooled_log2fpkm=log2(d$Bl6.pooled_fpkm+0.01)
# d$Castaneus.pooled_log2fpkm=log2(d$Castaneus.pooled_fpkm+0.01)
# 
# plot(density(d$Bl6.pooled_log2fpkm))
# lines(density(d$Castaneus.pooled_log2fpkm))
# cutoff=0
# abline(v=0,lty=2,col='red')
# 
# d$expressed_Bl6=FALSE
# d$expressed_Castaneus=FALSE
# d[d$Castaneus.pooled_log2fpkm>0,]$expressed_Castaneus=TRUE
# d[d$Bl6.pooled_log2fpkm>0,]$expressed_Bl6=TRUE
# 
# head(d[d$expressed_Bl6 ==FALSE  & d$expressed_Castaneus,])
# 
# write.csv(d, 'RNAseq_data_exprAnnotated.csv',row.names=FALSE,quote=FALSE)

##load gene expression data:
d=read.csv( 'RNAseq_data_exprAnnotated.csv')

##expressed genes:
expr=d[d$expressed_Bl6 | d$expressed_Castaneus,]$genenames

library(TFBSTools)
suppressMessages(library(JASPAR2018))
opts <- list()
opts[["species"]] <- 10090
opts[["species"]] <- 'Mus musculus'
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
PFMatrixList
#> PFMatrixList of length 1
#> names(1): MA0002.1

##access individual entries:
name(PFMatrixList[[100]]) ##shows gene name
name(PFMatrixList) ##shows all names

pfm <- getMatrixByName(JASPAR2018, name="Foxo1")
seqLogo(toICM(pfm))

##identify TFs in the matrix that are also expressed in Bl6 or Castaneus:
names=unique(expr[expr %in% name(PFMatrixList)])
motifs=PFMatrixList[name(PFMatrixList) %in% names]

# ##motif matching:
# library(motifmatchr)
# library(GenomicRanges)
# library(SummarizedExperiment)
# library(BSgenome)
# 
# # load some example motifs
# data(example_motifs, package = "motifmatchr") 
# 

##############################################
##Next step:
##Scan genomic sites for specific TFBSs
#############################################


library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10.masked)

castaneus=
replaceLetterAt(c19,1:4,c('A','C','T','G'))

motifgrl=GRangesList()

for (g in as.character(baited_genes$genenames)){
  print(g)
  
  n=paste0(g,'_connected')
  
  ##select all k27ac peaks:
  k=k27ac_all[k27ac_all[,n],]
  
  motifgr=GRanges()
  
  if (nrow(k)>0){
    
    
    for (i in 1:nrow(k)){
      print(i)
      
      ##select one K27ac peak
      a=getSeq(Mmusculus,bed2gr(k[i,range]))
      
      ##search this peak for TFBSs
      siteset=searchSeq(toPWM(motifs),a)
      as(siteset,'GRanges')->gr.siteset   ##converts to GenomicRanges
      # as(siteset,'DataFrame')->df.siteset ##converts to data.frame
      gr.siteset
      
      # ##check sequence (1-bound or 0-bound?) - 1-bound!
      # 
      # ##a) for + strand
      # gr=bed2gr(k27ac_all[i,18:20])
      # grseq=gr
      # end(grseq)=start(gr)+end(gr.siteset[1])-1
      # start(grseq)=start(gr)+start(gr.siteset[1])-1
      # getSeq(Mmusculus,  grseq)
      # grseq
      # 
      # ##b) for sequences on - strand:
      # gr=bed2gr(k27ac_all[i,18:20])
      # grseq=gr
      # end(grseq)=start(gr)+end(gr.siteset[length(gr.siteset)])-1
      # start(grseq)=start(gr)+start(gr.siteset[length(gr.siteset)])-1
      # getSeq(Mmusculus,  grseq)
      # grseq
      
      ##assign proper coordinates to gr.siteset:
      gr=bed2gr(k[i,range])
      grseq=gr.siteset
      end(grseq)=start(gr)+end(gr.siteset)-1
      start(grseq)=start(gr)+start(gr.siteset)-1
      seqlevels(grseq)=c(seqlevels(gr),seqlevels(grseq))
      seqlevels(gr)=seqlevels(grseq)
      seqnames(grseq)=seqnames(gr)
      getSeq(Mmusculus,  grseq)
      gr.siteset=grseq
      
      ##check which motifs have SNPs or SVs
      ov=as.matrix(findOverlaps(gr.siteset,snp))
      gr.siteset$SNP=FALSE
    
      if (length(gr.siteset[ov[,1]])>0){
        gr.siteset[ov[,1]]$SNP=TRUE
      }
      ov=as.matrix(findOverlaps(gr.siteset,sv))
      gr.siteset$SV=FALSE
      gr.siteset[ov[,1]]
      if (length(gr.siteset[ov[,1]])>0){
        gr.siteset[ov[,1]]$SV=TRUE
      }
      
      
      gr.siteset$K27ac_peakID=k[i,]$X
      gr.siteset$connected_gene=g
      
      seqlevels(motifgr)=unique(c(seqlevels(motifgr),seqlevels(gr.siteset)))
      seqlevels(gr.siteset)=seqlevels(motifgr)
      
      motifgr=c(motifgr,gr.siteset)
    }
    
  }
  
  motifgrl=c(motifgrl,GRangesList(motifgr))

}

names(motifgrl)=as.character(baited_genes$genenames)
saveRDS(motifgrl,'Motif_analysis_pluripotencyGenesConnected_K27ac_peaks.rds')


###########################################################
##plots:
###########################################################

motifgrl=readRDS('Motif_analysis_pluripotencyGenesConnected_K27ac_peaks.rds')

##Plot all motifs within all K27ac peaks for each gene

pdf('K27acPeaks_and_motifs..pdf',width=12,height=6,onefile=TRUE)

layout.matrix <- matrix(c(1:3), nrow = 1, ncol = 3)

layout(mat = layout.matrix,
       heights = 1, # Heights of the two rows
       widths = c(6, 1,1)) # Widths of the two columns


for (i in 1:length(motifgrl)){
  a=motifgrl[[i]]
  name=names(motifgrl)[i]
  chr=as.character(unique(seqnames(a)))
  
  for (e in unique(a$K27ac_peakID)){
    b=a[a$K27ac_peakID==e]
    
    n=paste0(e,';',k27ac_all[k27ac_all$X==e,]$Chr,':', k27ac_all[k27ac_all$X==e,]$Start, '-', k27ac_all[k27ac_all$X==e,]$End)
    
    ylim=c(0,length(unique(b$TF))+3)
    xlim=c(k27ac_all[k27ac_all$X==e,]$Start, k27ac_all[k27ac_all$X==e,]$End)
    
    plot(1,1,ylim=ylim,xlim=xlim,col='white',main=paste(name,'E:',n),ylab='',xlab=paste('position on',chr))
    rect(k27ac_all[k27ac_all$X==e,]$Start,ylim[2]-1,k27ac_all[k27ac_all$X==e,]$End,ylim[2])
    ##plot SNPs and SV
    rect(start(snp[seqnames(snp)==chr]),rep(ylim[2]-2,length(snp[seqnames(snp)==chr])),end(snp[seqnames(snp)==chr]),rep(ylim[2]-1,length(snp[seqnames(snp)==chr])),col='black')
    text(x=k27ac_all[k27ac_all$X==e,]$Start, y=(ylim[2]-1.5), labels='SNPs', cex=0.5,pos=2)
      rect(start(sv[seqnames(sv)==chr]),rep(ylim[2]-3,length(sv[seqnames(sv)==chr])),end(sv[seqnames(sv)==chr]),rep(ylim[2]-2,length(sv[seqnames(sv)==chr])))
      text(x=k27ac_all[k27ac_all$X==e,]$Start, y=(ylim[2]-2.5), labels='SVs', cex=0.5,pos=2)
    
      k=0
      
      col=c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(12,'Set3'),brewer.pal(9,'Pastel1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Pastel2'))
    
    for (tf in unique(b$TF)){
      
      rect(start(b[b$TF==tf]),rep(k,length(b[b$TF==tf])),end(b[b$TF==tf]),rep(k+1,length(b[b$TF==tf])),col=col[k+1],border=NA)
      text(x=k27ac_all[k27ac_all$X==e,]$Start, y=(k+k+1)/2, labels=tf, col=col[k+1],cex=0.4,pos=2)
      k=k+1
    
    }
      
      vb=vector()
      vc=vector()
      for (tf in unique(b$TF)){
        bl=mean(d[d$genenames==tf,]$Bl6.pooled_log2fpkm)
        ca=mean(d[d$genenames==tf,]$Castaneus.pooled_log2fpkm)
        vb=c(vb,bl)
        vc=c(vc,ca)
      }
      
      xlim=c(min(c(vb,vc,0)),max(c(vb,vc,0)))
      names=c(unique(b$TF),rep('',3))
      barplot(c(vb,0,0,0),horiz=TRUE,col=col[1:length(unique(b$TF))],main='RNASeq Bl6',xlim=xlim,border=NA, xlab='log2(FPKM)',names=names,las=1)
      abline(v=0)
      barplot(c(vc,0,0,0),horiz=TRUE,col=col[1:length(unique(b$TF))],main='RNASeq Cast',xlim=xlim,border=NA,xlab='log2(FPKM)')
      abline(v=0)
    
  }
  
}

dev.off()



# ##2) Plot all motifs with a SNP within all K27ac peaks for each gene
# 
# motifgrl=readRDS('Motif_analysis_pluripotencyGenesConnected_K27ac_peaks.rds')
# 
# pdf('K27acPeaks_and_motifsWithSNPs.pdf',width=12,height=5,onefile=TRUE)
# 
# layout.matrix <- matrix(c(1:3), nrow = 1, ncol = 3)
# 
# layout(mat = layout.matrix,
#        heights = 1, # Heights of the two rows
#        widths = c(6, 1,1)) # Widths of the two columns
# 
# # layout.show(3)
# 
# for (i in 1:length(motifgrl)){
#   print(names(motifgrl)[i])
#   a=motifgrl[[i]]
#   name=names(motifgrl)[i]
#   chr=as.character(unique(seqnames(a)))
#   
#   for (e in unique(a$K27ac_peakID)){
#     b=a[a$K27ac_peakID==e]
#     c=b[b$SNP | b$SV] ##select only motifs with SNP/SV
#     b=b[b$TF %in% c$TF] ##plot only TFs whose  motifs also have SNP/SV
#     
#     if (length(b)>0){
#     
#       n=paste0(e,';',k27ac_all[k27ac_all$X==e,]$Chr,':', k27ac_all[k27ac_all$X==e,]$Start, '-', k27ac_all[k27ac_all$X==e,]$End)
#       
#       ylim=c(0,length(unique(b$TF))+3)
#       xlim=c(k27ac_all[k27ac_all$X==e,]$Start, k27ac_all[k27ac_all$X==e,]$End)
#       
#       plot(1,1,ylim=ylim,xlim=xlim,col='white',main=paste(name,'E:',n),ylab='',xlab=paste('position on',chr))
#       rect(k27ac_all[k27ac_all$X==e,]$Start,ylim[2]-1,k27ac_all[k27ac_all$X==e,]$End,ylim[2])
#       ##plot SNPs and SV
#       rect(start(snp[seqnames(snp)==chr]),rep(ylim[2]-2,length(snp[seqnames(snp)==chr])),end(snp[seqnames(snp)==chr]),rep(ylim[2]-1,length(snp[seqnames(snp)==chr])),col='black',border=NA)
#       text(x=k27ac_all[k27ac_all$X==e,]$Start, y=(ylim[2]-1.5), labels='SNPs', cex=0.4,pos=2)
#         rect(start(sv[seqnames(sv)==chr]),rep(ylim[2]-3,length(sv[seqnames(sv)==chr])),end(sv[seqnames(sv)==chr]),rep(ylim[2]-2,length(sv[seqnames(sv)==chr])),border=NA,col='black')
#         text(x=k27ac_all[k27ac_all$X==e,]$Start, y=(ylim[2]-2.5), labels='SVs', cex=0.4,pos=2)
#       
#         k=0
#         
#         col=c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(12,'Set3'),brewer.pal(9,'Pastel1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Pastel2'))
#       
#       for (tf in unique(b$TF)){
#         
#         rect(start(b[b$TF==tf]),rep(k,length(b[b$TF==tf])),end(b[b$TF==tf]),rep(k+1,length(b[b$TF==tf])),col=col[k+1],border=NA)
#         text(x=k27ac_all[k27ac_all$X==e,]$Start, y=(k+k+1)/2, labels=tf, col=col[k+1],cex=0.6,pos=2)
#         k=k+1
#       
#       }
#         
#       vb=vector()
#       vc=vector()
#       for (tf in unique(b$TF)){
#         bl=mean(d[d$genenames==tf,]$Bl6.pooled_log2fpkm)
#         ca=mean(d[d$genenames==tf,]$Castaneus.pooled_log2fpkm)
#         vb=c(vb,bl)
#         vc=c(vc,ca)
#       }
#       
#       xlim=c(min(c(vb,vc)),max(c(vb,vc)))
#       barplot(c(vb,0,0),horiz=TRUE,col=col[1:length(unique(b$TF))],main='FPKM Bl6',xlim=xlim)
#       barplot(c(vc,0,0),horiz=TRUE,col=col[1:length(unique(b$TF))],main='FPKM Cast',xlim=xlim)
#   
#       }
#     
#   }
#   
# }
# 
# 
# dev.off()



```


#####################################################################
##Scan Castaneus genome for specific TFBSs
##Repeat the same!
#####################################################################

```{r}


##load gene expression data:
d=read.csv( 'RNAseq_data_exprAnnotated.csv')

##expressed genes:
expr=d[d$expressed_Bl6 | d$expressed_Castaneus,]$genenames

library(TFBSTools)
suppressMessages(library(JASPAR2018))
opts <- list()
opts[["species"]] <- 10090
opts[["species"]] <- 'Mus musculus'
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
PFMatrixList
#> PFMatrixList of length 1
#> names(1): MA0002.1

##access individual entries:
name(PFMatrixList[[100]]) ##shows gene name
name(PFMatrixList) ##shows all names

pfm <- getMatrixByName(JASPAR2018, name="Sox17")
seqLogo(toICM(pfm))

##identify TFs in the matrix that are also expressed in Bl6 or Castaneus:
names=unique(expr[expr %in% name(PFMatrixList)])
motifs=PFMatrixList[name(PFMatrixList) %in% names]

# ##motif matching:
# library(motifmatchr)
# library(GenomicRanges)
# library(SummarizedExperiment)
# library(BSgenome)
# 
# # load some example motifs
# data(example_motifs, package = "motifmatchr") 
# 

##############################################
##Next step:
##Scan genomic sites for specific TFBSs
#############################################


library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10.masked)

# ##create Cast genome
# chr=paste0('chr',c(1:19,'X','Y'))
# 
# for (ch in chr){
#   c19=Mmusculus[[ch]]
#   pos=start(snp[seqnames(snp)==ch])
#   s=as.character(snp[seqnames(snp)==ch]$type)
#   
#   i=1
#   ci=c19[1:(pos[i]-1)]
#   ctest=c19[pos[i]]   ##TO DO: take into account 3 letters replaced by one!
#   print(paste('correct letter for SNP is:',as.character(ctest)==unlist(strsplit(s[i],split='/'))[1]))
#   cadd=DNAString(unlist(strsplit(s[i],split='/'))[2])
#   final=c(ci,cadd)
#   
#   for (i in 2:length(pos)){
#     ci=c19[(pos[(i-1)]+1):(pos[i]-1)]
#     ctest=c19[pos[i]]
#     print(paste('correct letter for SNP is:',as.character(ctest)==unlist(strsplit(s[i],split='/'))[1]))
#     cadd=DNAString(unlist(strsplit(s[i],split='/'))[2])
#     final=c(final,ci,cadd)
#     
#   }
# 
# }
# 
# # pos=c(10000000,10000002,10000005)
# # s=c('N/AT','N/C','N/T')
# 
# 
# motifgrl=GRangesList()
# 
# for (g in as.character(baited_genes$genenames)){
#   print(g)
#   
#   n=paste0(g,'_connected')
#   
#   ##select all k27ac peaks:
#   k=k27ac_all[k27ac_all[,n],]
#   
#   motifgr=GRanges()
#   
#   if (nrow(k)>0){
#     
#     
#     for (i in 1:nrow(k)){
#       print(i)
#       
#       ##select one K27ac peak
#       a=getSeq(Mmusculus,bed2gr(k[i,range]))
#       
#       ##search this peak for TFBSs
#       siteset=searchSeq(toPWM(motifs),a)
#       as(siteset,'GRanges')->gr.siteset   ##converts to GenomicRanges
#       # as(siteset,'DataFrame')->df.siteset ##converts to data.frame
#       gr.siteset
#       
#       # ##check sequence (1-bound or 0-bound?) - 1-bound!
#       # 
#       # ##a) for + strand
#       # gr=bed2gr(k27ac_all[i,18:20])
#       # grseq=gr
#       # end(grseq)=start(gr)+end(gr.siteset[1])-1
#       # start(grseq)=start(gr)+start(gr.siteset[1])-1
#       # getSeq(Mmusculus,  grseq)
#       # grseq
#       # 
#       # ##b) for sequences on - strand:
#       # gr=bed2gr(k27ac_all[i,18:20])
#       # grseq=gr
#       # end(grseq)=start(gr)+end(gr.siteset[length(gr.siteset)])-1
#       # start(grseq)=start(gr)+start(gr.siteset[length(gr.siteset)])-1
#       # getSeq(Mmusculus,  grseq)
#       # grseq
#       
#       ##assign proper coordinates to gr.siteset:
#       gr=bed2gr(k[i,range])
#       grseq=gr.siteset
#       end(grseq)=start(gr)+end(gr.siteset)-1
#       start(grseq)=start(gr)+start(gr.siteset)-1
#       seqlevels(grseq)=c(seqlevels(gr),seqlevels(grseq))
#       seqlevels(gr)=seqlevels(grseq)
#       seqnames(grseq)=seqnames(gr)
#       getSeq(Mmusculus,  grseq)
#       gr.siteset=grseq
#       
#       ##check which motifs have SNPs or SVs
#       ov=as.matrix(findOverlaps(gr.siteset,snp))
#       gr.siteset$SNP=FALSE
#     
#       if (length(gr.siteset[ov[,1]])>0){
#         gr.siteset[ov[,1]]$SNP=TRUE
#       }
#       ov=as.matrix(findOverlaps(gr.siteset,sv))
#       gr.siteset$SV=FALSE
#       gr.siteset[ov[,1]]
#       if (length(gr.siteset[ov[,1]])>0){
#         gr.siteset[ov[,1]]$SV=TRUE
#       }
#       
#       
#       gr.siteset$K27ac_peakID=k[i,]$X
#       gr.siteset$connected_gene=g
#       
#       seqlevels(motifgr)=unique(c(seqlevels(motifgr),seqlevels(gr.siteset)))
#       seqlevels(gr.siteset)=seqlevels(motifgr)
#       
#       motifgr=c(motifgr,gr.siteset)
#     }
#     
#   }
#   
#   motifgrl=c(motifgrl,GRangesList(motifgr))
# 
# }
# 
# names(motifgrl)=as.character(baited_genes$genenames)
# saveRDS(motifgrl,'Motif_analysis_pluripotencyGenesConnected_K27ac_peaks.rds')


###########################################################
##plots:
###########################################################

##Do the following once:

motifgrl=readRDS('Motif_analysis_pluripotencyGenesConnected_K27ac_peaks.rds')

m=unlist(motifgrl)
ov=as.matrix(findOverlaps(m,snp))
m$snp_pos='.'
m[ov[,1]]$snp_pos=start(snp[ov[,2]])
m$snp_type='.'
m[ov[,1]]$snp_type=as.character(snp[ov[,2]]$type)

mnosnp=m[m$SNP==FALSE]
msnp=m[m$SNP==TRUE]

msnp$Strand=strand(msnp)
strand(msnp)='*'
seqs=getSeq(Mmusculus,msnp)

##change sequene to incorporate SNPs:

msnp$TF.Cast='none'
msnp$site='destroyed'

# for (i in 4971:4975){
  for (i in 1:length(seqs)){
  print(i)
  
  seq=DNAString(seqs[[i]])
  pos=as.numeric(msnp[i]$snp_pos)-as.numeric(start(msnp[i]))
  x=unlist(strsplit(msnp[i]$snp_type,split='_'))[2]
  x=unlist(strsplit(x,split=','))[1]  ##sometimes there are two replacement strings
  toreplace=unlist(strsplit(x,split='/'))[1]
  len=nchar(toreplace)
  
  ##Careful with '-' strand motifs! How to translate?
  
  ##Take care to add some more sequence to the sequence if larger chunk is removed
  if(len>1){
    seq=unlist(getSeq(Mmusculus,GRanges(seqnames(msnp[i]),IRanges(start(msnp[i]),(end(msnp[i])+len) ))))
  }
  
  replacewith=DNAString(unlist(strsplit(x,split='/'))[2])
  ci=seq[1:pos]
  ctest=seq[(pos+1):(pos+len)]   ##TO DO: take into account 3 letters replaced by one!
  print(paste('correct letter for SNP is:',as.character(ctest)==toreplace))
  ##accounts for the possibility that the very last letter is a SNP:
  if (pos+len < length(seq)){
    after=seq[(pos+len+1):length(seq)]
  } else {
    after=DNAString('')
    }
  final=c(ci,replacewith,after)
  # seqs_cast=DNAStringSet(final)
  
  siteset=searchSeq(toPWM(motifs),final)
  ##search this peak for TFBSs
  as(siteset,'GRanges')->gr.siteset   ##converts to GenomicRanges
  # as(siteset,'DataFrame')->df.siteset ##converts to data.frame
  gr.siteset
  if(length(gr.siteset)>0){
    if (length(gr.siteset)==1){
      msnp[i]$TF.Cast=gr.siteset$TF
    } else {
      msnp[i]$TF.Cast=as.character(concatenate(gr.siteset$TF))
    }
    if (as.logical( msnp[i]$TF  %in% gr.siteset$TF)){
      msnp[i]$site='remains'
    }
  } 

}
      
saveRDS(msnp,'Motif_analysis_pluripotencyGenesConnected_K27ac_peaks_SNPsites.rds')

  
##Plot all motifs within all K27ac peaks for each gene
msnp=readRDS('Motif_analysis_pluripotencyGenesConnected_K27ac_peaks_SNPsites.rds')
snpgrl=GRangesList()
for (n in unique(names(msnp))){
  snpgrl=c(snpgrl,GRangesList(msnp[names(msnp)==n]))
}
names(snpgrl)=unique(names(msnp))

# pdf('K27acPeaks_and_motifs_Bl6vsCast_OnlySNPOv.pdf',width=10,height=5,onefile=TRUE)
# 
# layout.matrix <- matrix(c(1:3), nrow = 1, ncol = 3)
# 
# layout(mat = layout.matrix,
#        heights = 1, # Heights of the two rows
#        widths = c(6, 1,1)) # Widths of the two columns
# 
# 
# for (i in 1:length(snpgrl)){
#   a=snpgrl[[i]]
#   name=names(snpgrl)[i]
#   chr=as.character(unique(seqnames(a)))
#   
#   ##overlap a with snps:
#   snp$id=1:length(snp)
#   ov=as.matrix(findOverlaps(a,snp))
#   a=a[ov[,1]]
#   a$snp_id=snp[ov[,2]]$id
#   
#   for (e in unique(a$K27ac_peakID)){
#     bb=a[a$K27ac_peakID==e]
#     
#     nsnps=length(unique(bb$snp_id))
# 
#     
#     snpids=unique(bb$snp_id)
#     snpids=snpids[order(snpids)]
#     
#     for (s in snpids){
#       
#       b=bb[bb$snp_id==s]
#       
#       plotsite=resize(snp[snp$id==s],35,fix='center')
#       
#       n=paste0(e,';',k27ac_all[k27ac_all$X==e,]$Chr,':', k27ac_all[k27ac_all$X==e,]$Start, '-', k27ac_all[k27ac_all$X==e,]$End,' SNP',start(snp[snp$id==s]))
#       
#       ylim=c(0,length(unique(b$TF))+3)
#       xlim=c(start(plotsite), end(plotsite))
#       
#       plot(1,1,ylim=ylim,xlim=xlim,col='white',main=paste(name,'E:',n),ylab='',xlab=paste('position on',chr))
#       rect(k27ac_all[k27ac_all$X==e,]$Start,ylim[2]-1,k27ac_all[k27ac_all$X==e,]$End,ylim[2])
#       ##plot SNPs and SV
#       rect(start(snp[seqnames(snp)==chr]),rep(ylim[2]-2,length(snp[seqnames(snp)==chr])),end(snp[seqnames(snp)==chr]),rep(ylim[2]-1,length(snp[seqnames(snp)==chr])),col='black',border='black')
#       text(x=xlim[1], y=(ylim[2]-1.5), labels='SNPs', cex=1,pos=NULL)
#         rect(start(sv[seqnames(sv)==chr]),rep(ylim[2]-3,length(sv[seqnames(sv)==chr])),end(sv[seqnames(sv)==chr]),rep(ylim[2]-2,length(sv[seqnames(sv)==chr])),border=NA)
#         text(x=xlim[1], y=(ylim[2]-2.5), labels='SVs', cex=1,pos=NULL)
#       
#         k=0
#         
#         col=c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(12,'Set3'),brewer.pal(9,'Pastel1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Pastel2'))
#       
#       for (tf in unique(b$TF)){
#         
#         rect(start(b[b$TF==tf]),rep(k,length(b[b$TF==tf])),end(b[b$TF==tf]),rep(k+1,length(b[b$TF==tf])),col=col[k+1],border=NA)
#         text(x=xlim[1], y=(k+k+1)/2, labels=tf, col=col[k+1],cex=1,pos=NULL)
# 
#         if ( length(b[b$TF==tf & b$site=='destroyed' & b$TF.Cast=='none',])>0){
#           text(x=start(snp[snp$id==s]), y=(k+k+1)/2, labels='X',cex=3,pos=NULL)
#         } else if (length(b[b$TF==tf & b$site=='destroyed' & b$TF.Cast!='none',])>0){
#           text(x=start(snp[snp$id==s]), y=(k+k+1)/2, labels=b[b$TF==tf & b$site=='destroyed' & b$TF.Cast!='none',]$TF.Cast, cex=1,pos=NULL)
#         }
#         
#         
#         
#         k=k+1
#       
#       }
#         
#       vb=vector()
#       vc=vector()
#       for (tf in unique(b$TF)){
#         bl=mean(d[d$genenames==tf,]$Bl6.pooled_log2fpkm)
#         ca=mean(d[d$genenames==tf,]$Castaneus.pooled_log2fpkm)
#         vb=c(vb,bl)
#         vc=c(vc,ca)
#       }
#       
#       xlim=c(min(c(vb,vc)),max(c(vb,vc)))
#       barplot(c(vb,0,0),horiz=TRUE,col=col[1:length(unique(b$TF))],main='FPKM Bl6',xlim=xlim)
#       barplot(c(vc,0,0),horiz=TRUE,col=col[1:length(unique(b$TF))],main='FPKM Cast',xlim=xlim)
#     
#     
#   }
#   
# }
#   
# }
# 
# dev.off()



pdf('K27acPeaks_and_motifs_Bl6vsCast_OnlySNPOv.pdf',width=7,height=3.5,onefile=TRUE,pointsize=10)

layout.matrix <- matrix(c(1:3), nrow = 1, ncol = 3)

layout(mat = layout.matrix,
       heights = 1, # Heights of the two rows
       widths = c(5, 1.5,1.5)) # Widths of the two columns


for (i in 1:length(snpgrl)){
  a=snpgrl[[i]]
  name=names(snpgrl)[i]
  chr=as.character(unique(seqnames(a)))
  
  ##overlap a with snps:
  snp$id=1:length(snp)
  ov=as.matrix(findOverlaps(a,snp))
  a=a[ov[,1]]
  a$snp_id=snp[ov[,2]]$id
  
  for (e in unique(a$K27ac_peakID)){
    bb=a[a$K27ac_peakID==e]
    
    nsnps=length(unique(bb$snp_id))
    
    
    snpids=unique(bb$snp_id)
    snpids=snpids[order(snpids)]
    
    for (s in snpids){
      
      b=bb[bb$snp_id==s]
      
      plotsite=resize(snp[snp$id==s],35,fix='center')
      
      n=paste0(e,';',k27ac_all[k27ac_all$X==e,]$Chr,':', k27ac_all[k27ac_all$X==e,]$Start, '-', k27ac_all[k27ac_all$X==e,]$End,' SNP',start(snp[snp$id==s]))
      
      ylim=c(0,length(unique(b$TF))+2)
      xlim=c(start(plotsite), end(plotsite))
      
      plot(1,1,ylim=ylim,xlim=xlim,col='white',main=paste(name,'E:',n),ylab='',xlab=paste('motif position on',chr))
      ##plot SNPs and SV
      rect(start(snp[seqnames(snp)==chr]),rep(ylim[2]-1,length(snp[seqnames(snp)==chr])),end(snp[seqnames(snp)==chr]),rep(ylim[2],length(snp[seqnames(snp)==chr])),col='black',border='black')
      text(x=xlim[1], y=(ylim[2]-0.5), labels='SNPs', cex=0.75,pos=4)
      rect(start(sv[seqnames(sv)==chr]),rep(ylim[2]-2,length(sv[seqnames(sv)==chr])),end(sv[seqnames(sv)==chr]),rep(ylim[2]-1,length(sv[seqnames(sv)==chr])),border=NA)
      text(x=xlim[1], y=(ylim[2]-1.5), labels='SVs', cex=0.75,pos=4)
      
      k=0
      
      col=c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(12,'Set3'),brewer.pal(9,'Pastel1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Pastel2'))
      
      for (tf in unique(b$TF)){
        
        rect(start(b[b$TF==tf]),rep(k,length(b[b$TF==tf])),end(b[b$TF==tf]),rep(k+1,length(b[b$TF==tf])),col=col[k+1],border=NA)
        text(x=xlim[1], y=(k+k+1)/2, labels=tf, col=col[k+1],cex=0.75,pos=4)
        
        if ( length(b[b$TF==tf & b$site=='destroyed' & b$TF.Cast=='none',])>0){
          text(x=start(snp[snp$id==s]), y=(k+k+1)/2, labels='X',cex=3,pos=NULL)
        } else if (length(b[b$TF==tf & b$site=='destroyed' & b$TF.Cast!='none',])>0){
          text(x=start(snp[snp$id==s]), y=(k+k+1)/2, labels=b[b$TF==tf & b$site=='destroyed' & b$TF.Cast!='none',]$TF.Cast, cex=1,pos=NULL)
        }
        
        legend('topright',bty='n',legend='X=motif destroyed')
        
        
        k=k+1
        
      }
      
      vb=vector()
      vc=vector()
      for (tf in unique(b$TF)){
        bl=mean(d[d$genenames==tf,]$Bl6.pooled_log2fpkm)
        ca=mean(d[d$genenames==tf,]$Castaneus.pooled_log2fpkm)
        vb=c(vb,bl)
        vc=c(vc,ca)
      }
      
      xlim=c(floor(min(c(vb,vc,0))),ceiling(max(c(vb,vc,0))))
      # barplot(c(vb,0,0),horiz=TRUE,col=col[1:length(unique(b$TF))],xlab='FPKM',xlim=xlim,main='RNAseq Bl6',border=NA)
      # abline(v=0)
      # barplot(c(vc,0,0),horiz=TRUE,col=col[1:length(unique(b$TF))],xlab='FPKM',main='RNAseq Cast',xlim=xlim,border=NA)
      # abline(v=0)
      names=c(unique(b$TF),rep('',2))
      barplot(c(vb,0,0),horiz=TRUE,col=col[1:length(unique(b$TF))],main='RNASeq Bl6',xlim=xlim,border=NA, xlab='log2(FPKM)',names=names,las=1)
      abline(v=0)
      barplot(c(vc,0,0),horiz=TRUE,col=col[1:length(unique(b$TF))],main='RNASeq Cast',xlim=xlim,border=NA,xlab='log2(FPKM)')
      abline(v=0)
      
    }
    
  }
  
}

dev.off()


```









########################################################################
########################################################################
##4. Interval analysis (integration of ChIP-Seq peaks)
########################################################################
########################################################################

```{r}

names(grl)
num=2:8
num=2:4 ##k27ac only!

############################################
##1) get OneGeneOnePeak tables:
############################################

getIntervalInteractionPeaks(grl[num],hind,ints)

##look at your tables:
fls=list.files(pattern='oneGeneOnePeak')
f=read.delim(fls[1])
head(f)

# ##illustrate oneGeneOnePeak interactions:
# par(mfrow=c(3,1))
# ylim=c(0,40)
# k=5
# zoom=600000
# id=5402489
# name=paste(baited_genes[baited_genes$dpn_id==id,]$genenames)
# xlim=c(14260000,14300000)
# plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name=name,intervals=grl,hind=hind,d=ints,xlim=xlim)
# d=read.delim('oneGeneOnePeak_Ring1b_peaks_interacting.txt',stringsAsFactors=FALSE)
# plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name='Ring1b_peaks',intervals=grl,hind=hind,d=d,xlim=xlim)
# d=read.delim('oneGeneOnePeak_Cdk8_peaks_interacting.txt',stringsAsFactors=FALSE)
# plotInteractions(L,id,k,zoom,ylim=ylim,show.legend = TRUE,name='Cdk8_peaks',intervals=grl,hind=hind,d=d,xlim=xlim)

###############################################
##2) Make matrices containing each interaction
###############################################

memory.limit(size=32000)
zoom=100

for (i in num){

  print(names(grl)[i]) 
  dis=names(grl)[i]
  d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting.txt'),stringsAsFactors=FALSE)
  d=d[d$seqnames_bait==d$seqnames_otherEnd,]
  
  ##average data
  
  resfolder='matrices_newk27acPeaks/'
  dir.create(resfolder)
  
  ylim=c(0,250)
  m1=getMatrix(L,zoom,d,resfolder=resfolder,type='reads',norm=FALSE,name=names(grl)[i],ylim=ylim,xlim=c(-40,40),hind=hind)
  # getMatrix(L,zoom,d,resfolder=resfolder,type='reads',norm=FALSE,name=names(grl)[i],ylim=ylim,xlim=c(-40,40),hind=hind,troubleshooting=TRUE)
  ylim=c(0,5)
  m2=getMatrix(L,zoom,d,resfolder=resfolder,type='reads',norm=TRUE,name=names(grl)[i],normsam=2,xlim=c(-40,40),hind=hind,ylim=ylim)
  
  ylim=c(0,10)
  m3=getMatrix(L,zoom,d,resfolder=resfolder,type='scores',norm=FALSE,name=names(grl)[i],ylim=ylim,xlim=c(-40,40),hind=hind)
  
  ##replicate data
  
  resfolder='matrices_newk27acPeaks_reps/'
  dir.create(resfolder)
  
  ylim=c(0,250)
  m4=getMatrix(LL,zoom,d,resfolder=resfolder,type='reads',norm=FALSE,name=names(grl)[i],ylim=ylim,xlim=c(-40,40),hind=hind)
  m5=getMatrix(LL,zoom,d,resfolder=resfolder,type='reads',norm=TRUE,name=names(grl)[i],normsam=4,xlim=c(-40,40),hind=hind)
  ylim=c(0,10)
  m6=getMatrix(LL,zoom,d,resfolder=resfolder,type='scores',norm=FALSE,name=names(grl)[i],ylim=ylim,xlim=c(-40,40),hind=hind)
  
}

##annotate d with quantified interaction reads and scores

##average data:
resfolder='matrices_newk27acPeaks/'
getDataInPeaks(grl[num],zoom,resfolder=resfolder,overwrite=TRUE)

##replicates:
getDataInPeaks(grl[num],zoom,resfolder='matrices_newk27acPeaks_reps/')

#########
##plots:
#########

pdf('BoxplotsCC.pdf',onefile=TRUE)
boxplotsCC(grl,num=num,resfolder=resfolder)
dev.off()

##heatmaps to look at replicate correlations:
heatmapsCC(grl,num=num,resfolder='matrices_newk27acPeaks_reps/')

```


##Interaction strength with connected K27ac peaks by different genes:
##Barplots, no figure (but informative!)

```{r}


par(mfrow=c(2,2))

i=2
dr=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_reps_interacting_readsScoresInPeak.txt'),stringsAsFactors=FALSE)
d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_interacting_readsScoresInPeak.txt'),stringsAsFactors=FALSE)

##assign regulation of an interaction:
d$reg.int='stable'
d[d$Bl6.sum_bednorm_reads>d$Castaneus.sum_bednorm_reads & d$Bl6.sum_bednorm_scores > d$Castaneus.sum_bednorm_scores &
    dr$Bl6_rep3.sum_bednorm_reads>dr$Castaneus_rep2.sum_bednorm_reads & 
    dr$Bl6_rep4.sum_bednorm_reads>dr$Castaneus_rep3.sum_bednorm_reads &
    dr$Bl6_rep5.sum_bednorm_reads>dr$Castaneus_rep4.sum_bednorm_reads &
    dr$Bl6_rep6.sum_bednorm_reads>dr$Castaneus_rep5.sum_bednorm_reads &
    dr$Bl6_rep3.sum_bednorm_scores>dr$Castaneus_rep2.sum_bednorm_scores &
    dr$Bl6_rep4.sum_bednorm_scores>dr$Castaneus_rep3.sum_bednorm_scores &
    dr$Bl6_rep5.sum_bednorm_scores>dr$Castaneus_rep4.sum_bednorm_scores &
    dr$Bl6_rep6.sum_bednorm_scores>dr$Castaneus_rep5.sum_bednorm_scores ,]$reg.int='Bl6up'

d[d$Bl6.sum_bednorm_reads<d$Castaneus.sum_bednorm_reads & d$Bl6.sum_bednorm_scores < d$Castaneus.sum_bednorm_scores &
    dr$Bl6_rep3.sum_bednorm_reads<dr$Castaneus_rep2.sum_bednorm_reads & 
    dr$Bl6_rep4.sum_bednorm_reads<dr$Castaneus_rep3.sum_bednorm_reads &
    dr$Bl6_rep5.sum_bednorm_reads<dr$Castaneus_rep4.sum_bednorm_reads &
    dr$Bl6_rep6.sum_bednorm_reads<dr$Castaneus_rep5.sum_bednorm_reads &
    dr$Bl6_rep3.sum_bednorm_scores<dr$Castaneus_rep2.sum_bednorm_scores &
    dr$Bl6_rep4.sum_bednorm_scores<dr$Castaneus_rep3.sum_bednorm_scores &
    dr$Bl6_rep5.sum_bednorm_scores<dr$Castaneus_rep4.sum_bednorm_scores &
    dr$Bl6_rep6.sum_bednorm_scores<dr$Castaneus_rep5.sum_bednorm_scores ,]$reg.int='Castup'


bUpGenes=d[d$genetype=='down_strong',]
bUpGenes.bUpPeaks=d[d$genetype=='down_strong' & d$LFC_k27acPeak_oE <= -log2(1.5),]
cUpGenes=d[d$genetype=='up_strong',]
cUpGenes.cUpPeaks=d[d$genetype=='up_strong' & d$LFC_k27acPeak_oE >= log2(1.5),]

p=8
plot(log10(d$Bl6.sum_bednorm_reads+p),log10(d$Castaneus.sum_bednorm_reads+p))
abline(0,1,col='red')
points(log10(bUpGenes$Bl6.sum_bednorm_reads+p),log10(bUpGenes$Castaneus.sum_bednorm_reads+p),col='cornflowerblue',pch=20)
points(log10(bUpGenes.bUpPeaks$Bl6.sum_bednorm_reads+p),log10(bUpGenes.bUpPeaks$Castaneus.sum_bednorm_reads+p),col='darkblue',pch=20)

p=8
plot(log10(d$Bl6.sum_bednorm_reads+p),log10(d$Castaneus.sum_bednorm_reads+p))
abline(0,1,col='red')
points(log10(cUpGenes$Bl6.sum_bednorm_reads+p),log10(cUpGenes$Castaneus.sum_bednorm_reads+p),col='orange',pch=20)
points(log10(cUpGenes.cUpPeaks$Bl6.sum_bednorm_reads+p),log10(cUpGenes.cUpPeaks$Castaneus.sum_bednorm_reads+p),col='darkred',pch=20)


p=1
plot(log10(d$Bl6.sum_bednorm_scores+p),log10(d$Castaneus.sum_bednorm_scores+p))
abline(0,1,col='red')
points(log10(bUpGenes$Bl6.sum_bednorm_scores+p),log10(bUpGenes$Castaneus.sum_bednorm_scores+p),col='cornflowerblue',pch=20)
points(log10(bUpGenes.bUpPeaks$Bl6.sum_bednorm_scores+p),log10(bUpGenes.bUpPeaks$Castaneus.sum_bednorm_scores+p),col='darkblue',pch=20)

p=1
plot(log10(d$Bl6.sum_bednorm_scores+p),log10(d$Castaneus.sum_bednorm_scores+p))
abline(0,1,col='red')
points(log10(cUpGenes$Bl6.sum_bednorm_scores+p),log10(cUpGenes$Castaneus.sum_bednorm_scores+p),col='orange',pch=20)
points(log10(cUpGenes.cUpPeaks$Bl6.sum_bednorm_scores+p),log10(cUpGenes.cUpPeaks$Castaneus.sum_bednorm_scores+p),col='darkred',pch=20)

##Barplots
col=c('cornflowerblue','orange','grey')

##regulation of K27ac peak interactions (binary): 
tab=table(d$reg.int,d$peaktype_k27acPeak_oE)
plotNormBarplot(tab,ylab='fraction interacting peaks',legend.title = 'regulation of interactions',xlab='K27ac peaktype',col=col)

##Zoom in specifically onto Bl6-peaks with Bl6-up ints (and similar for Cast):
b=d[d$peaktype_k27acPeak_oE=='Bl6' & d$reg.int=='Bl6up',]
c=d[d$peaktype_k27acPeak_oE=='Castaneus' & d$reg.int=='Castup',]
tab=data.frame(table(b$genetype),table(c$genetype))
row.names(tab)=tab[,1]
tab=tab[,c(2,4)]
colnames(tab)=c('Bl6Peaks','CastPeaks')
plotNormBarplot(as.matrix(tab),ylab='fraction interacting peaks',legend.title = 'interacting gene type',xlab='K27ac peaktype ChIP+CapC',col=col)

##Zoom in specifically onto all Bl6-peaks (and similar for Cast):
b=d[d$peaktype_k27acPeak_oE=='Bl6' ,]
c=d[d$peaktype_k27acPeak_oE=='Castaneus',]
tab=data.frame(table(b$genetype),table(c$genetype))
row.names(tab)=tab[,1]
tab=tab[,c(2,4)]
colnames(tab)=c('Bl6Peaks','CastPeaks')
plotNormBarplot(as.matrix(tab),ylab='fraction interacting peaks',legend.title = 'interacting gene type',xlab='K27ac peaktype ChIP',col=col)

##Zoom in specifically onto Bl6-higher-interacting peaks (and similar for Cast):
b=d[d$reg.int=='Bl6up' ,]
c=d[d$reg.int=='Castup',]
tab=data.frame(table(b$genetype),table(c$genetype))
row.names(tab)=tab[,1]
tab=tab[,c(2,4)]
colnames(tab)=c('Bl6Peaks','CastPeaks')
plotNormBarplot(as.matrix(tab),ylab='fraction interacting peaks',legend.title = 'interacting gene type',xlab='K27ac peaktype CapC',col=col)

##Show all peaks:
tab=data.frame(table(d$genetype),table(baited_genes$genetype))
row.names(tab)=tab[,1]
tab=tab[,c(2,4)]
colnames(tab)=c('Connected_peaks','Captured_genes')
plotNormBarplot(as.matrix(tab),ylab='fraction',legend.title = 'gene type',col=col)


##############################################
##plot individual interactions from matrices!
##############################################

##Sort d to show similar types of interactions together:
d=d[order(d$genetype),]
d=d[order(d$LFC_k27acPeak_oE),]
d=d[order(d$reg.int),]

m=readRDS('matrices_newk27acPeaks/matrices_K27ac_2i_all_normTo_Bl6_noRep4_reads.rds')
m=m[c(2,3,5,6)]
col=c('cornflowerblue','darkorange','cornflowerblue','darkorange')
lty=c(1,1,2,2)

k=5
xlim=c(-20,20)
xlab='Distance from central DpnII Frag [#Frags]'
ylab='enrichment vs. enrichment center Bl6' 

frags=readRDS('K27ac_2i_all_interactingInData_hindFrag_annotated.rds') ##positioning of K27ac peaks

  
pdf('Example_ints_normTo.Bl6noRep4.pdf',width=18,height=10,onefile=TRUE)
par(mfrow=c(3,5))

for (i in 1:nrow(d)){
ylim=c(0,2)
  n=d[d$intID==row.names(m[[1]])[i],]
  main=paste(n$genename,n$genetype,n$intID,'int:',n$reg.int,'LFC:',round(n$LFC_k27acPeak_oE,2),sep='|')
  
  peak=n$otherEndID
  f=frags[frags$id==peak]
  if(n$baitID<n$otherEndID){
    max= -f$min
    min= -f$max
  } else {
    min=f$min
    max=f$max
  }
  
  colrect='grey'
  if(n$reg.int=='Castup'){colrect='yellow'}
  if(n$reg.int=='Bl6up'){colrect='green'}
  
  j=1
  plot( (-(ncol(m[[j]])-1)/2) : ((ncol(m[[j]])-1)/2),runmean(Rle(t(m[[j]][i,])), k=k,endrule = 'constant'),type='l',col=col[j],lty=lty[j],main=main,ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab)
  abline(v=0,lty=2)
  if (abs(n$LFC_k27acPeak_oE)<log2(1.5)){
    if(colrect=='yellow'){
      rect(min,-0.1,max,2.1,col=rgb(1,1,0,alpha=0.2),border='NA')
    } else if(colrect=='green'){
      rect(min,-0.1,max,2.1,col=rgb(0,1,0,alpha=0.2),border='NA')
    } else{rect(min,-0.1,max,2.1,col=rgb(0,0,0,alpha=0.2),border='NA')}
  } else if (n$LFC_k27acPeak_oE <= -log2(1.5)){
    if(colrect=='yellow'){
      rect(min,-0.1,max,2.1,col=rgb(1,1,0,alpha=0.2),border='blue')
    } else if(colrect=='green'){
      rect(min,-0.1,max,2.1,col=rgb(0,1,0,alpha=0.2),border='blue')
    } else{rect(min,-0.1,max,2.1,col=rgb(0,0,0,alpha=0.2),border='blue')}
  } else {
    if(colrect=='yellow'){
      rect(min,-0.1,max,2.1,col=rgb(1,1,0,alpha=0.2),border='orange')
    } else if(colrect=='green'){
      rect(min,-0.1,max,2.1,col=rgb(0,1,0,alpha=0.2),border='orange')
    } else{rect(min,-0.1,max,2.1,col=rgb(0,0,0,alpha=0.2),border='orange')}
  }
  for (j in 2:length(m)){
    lines( (-(ncol(m[[j]])-1)/2) : ((ncol(m[[j]])-1)/2),runmean(Rle(t(m[[j]][i,])), k=k,endrule = 'constant'),type='l',col=col[j],lty=lty[j])
  }
  
  legend('topright',legend=c('Bl6','Cast'),lwd=2,col=c('cornflowerblue','darkorange'),bty='n',cex=0.75)
  legend('topleft',fill=rep(c('green','yellow'),each=2),border = rep(c('blue','orange'),2),legend=c('intBl6up.peakBl6up','intBl6up.peakCastup','intCastup.peakBl6up','intCastup.peakCastup'),bty='n',cex=0.6)
  
}


dev.off()


```

##Boxplots of interaction change at all connected peaks by genes
##Boxplots, no figure, but informative!

```{r}

i=2
dr=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_reps_interacting_readsScoresInPeak.txt'),stringsAsFactors=FALSE)
d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_interacting_readsScoresInPeak.txt'),stringsAsFactors=FALSE)

##assign regulation of an interaction:
d$reg.int='stable'
d[d$Bl6.sum_bednorm_reads>d$Castaneus.sum_bednorm_reads & d$Bl6.sum_bednorm_scores > d$Castaneus.sum_bednorm_scores &
    dr$Bl6_rep3.sum_bednorm_reads>dr$Castaneus_rep2.sum_bednorm_reads & 
    dr$Bl6_rep4.sum_bednorm_reads>dr$Castaneus_rep3.sum_bednorm_reads &
    dr$Bl6_rep5.sum_bednorm_reads>dr$Castaneus_rep4.sum_bednorm_reads &
    dr$Bl6_rep6.sum_bednorm_reads>dr$Castaneus_rep5.sum_bednorm_reads &
    dr$Bl6_rep3.sum_bednorm_scores>dr$Castaneus_rep2.sum_bednorm_scores &
    dr$Bl6_rep4.sum_bednorm_scores>dr$Castaneus_rep3.sum_bednorm_scores &
    dr$Bl6_rep5.sum_bednorm_scores>dr$Castaneus_rep4.sum_bednorm_scores &
    dr$Bl6_rep6.sum_bednorm_scores>dr$Castaneus_rep5.sum_bednorm_scores ,]$reg.int='Bl6up'

d[d$Bl6.sum_bednorm_reads<d$Castaneus.sum_bednorm_reads & d$Bl6.sum_bednorm_scores < d$Castaneus.sum_bednorm_scores &
    dr$Bl6_rep3.sum_bednorm_reads<dr$Castaneus_rep2.sum_bednorm_reads & 
    dr$Bl6_rep4.sum_bednorm_reads<dr$Castaneus_rep3.sum_bednorm_reads &
    dr$Bl6_rep5.sum_bednorm_reads<dr$Castaneus_rep4.sum_bednorm_reads &
    dr$Bl6_rep6.sum_bednorm_reads<dr$Castaneus_rep5.sum_bednorm_reads &
    dr$Bl6_rep3.sum_bednorm_scores<dr$Castaneus_rep2.sum_bednorm_scores &
    dr$Bl6_rep4.sum_bednorm_scores<dr$Castaneus_rep3.sum_bednorm_scores &
    dr$Bl6_rep5.sum_bednorm_scores<dr$Castaneus_rep4.sum_bednorm_scores &
    dr$Bl6_rep6.sum_bednorm_scores<dr$Castaneus_rep5.sum_bednorm_scores ,]$reg.int='Castup'

##Calculate log2FC of interaction:
p=8
d$log2FC_CapC_Cast.vs.Bl6_sum_bednorm_reads=log2((d$Castaneus.sum_bednorm_reads+p)/(d$Bl6.sum_bednorm_reads+p))
p=1
d$log2FC_CapC_Cast.vs.Bl6_sum_bednorm_scores=log2((d$Castaneus.sum_bednorm_scores+p)/(d$Bl6.sum_bednorm_scores+p))

par(mfcol=c(2,4))
col=c( "cornflowerblue", "orange"  ,       "grey" )
ylim=c(-1,1)
boxplot(log2FC_CapC_Cast.vs.Bl6_sum_bednorm_reads~genetype,data=d,main='all connected K27ac peaks',outline=FALSE,notch=TRUE,col=col,xlab='genetype Cast.vs.Bl6',ylab='log2FC reads',ylim=ylim)
abline(h=0,lty=2)
legend('topleft',legend=c(paste0('nints=',nrow(d)),paste0('nproms=',length(unique(d$baitID)))  ),bty='n',cex=0.8 )
ylim=c(-2.5,2.5)
boxplot(log2FC_CapC_Cast.vs.Bl6_sum_bednorm_scores~genetype,data=d,main='all connected K27ac peaks',outline=FALSE,notch=TRUE,col=col,xlab='genetype Cast.vs.Bl6',ylab='log2FC scores',ylim=ylim)
abline(h=0,lty=2)
legend('topleft',legend=c(paste0('nints=',nrow(d)),paste0('nproms=',length(unique(d$baitID)))  ),bty='n',cex=0.8 )

for (cl in unique(d$peaktype_k27acPeak_oE)){
  ylim=c(-1,1)
  boxplot(log2FC_CapC_Cast.vs.Bl6_sum_bednorm_reads~genetype,data=d[d$peaktype_k27acPeak_oE==cl,],main=paste(cl,' connected K27ac peaks'),outline=FALSE,notch=TRUE,col=col,xlab='genetype Cast.vs.Bl6',ylab='log2FC reads',ylim=ylim)
  abline(h=0,lty=2)
  legend('topleft',legend=c(paste0('nints=',nrow(d[d$peaktype_k27acPeak_oE==cl,])),paste0('nproms=',length(unique(d[d$peaktype_k27acPeak_oE==cl,]$baitID)))  ),bty='n',cex=0.8 )
  ylim=c(-2.5,2.5)
   boxplot(log2FC_CapC_Cast.vs.Bl6_sum_bednorm_scores~genetype,data=d[d$peaktype_k27acPeak_oE==cl,],main=paste(cl,' connected K27ac peaks'),outline=FALSE,notch=TRUE,col=col,xlab='genetype Cast.vs.Bl6',ylab='log2FC scores',ylim=ylim)
  abline(h=0,lty=2)
  legend('topleft',legend=c(paste0('nints=',nrow(d[d$peaktype_k27acPeak_oE==cl,])),paste0('nproms=',length(unique(d[d$peaktype_k27acPeak_oE==cl,]$baitID)))  ),bty='n',cex=0.8 )
}

```


##Boxplots of interaction change for divergent peaks by different cutoffs for divergency

```{r}

pdf('Interactions_by_peakIdentity.pdf',width=16,height=10,onefile=TRUE)

i=2
d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_interacting_readsScoresInPeak.txt'),stringsAsFactors=FALSE)
dr=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_reps_interacting_readsScoresInPeak.txt'),stringsAsFactors=FALSE)

p=8
d$log2FC.sum_bednorm_reads.Cast.vs.Bl6 = log2((d$Castaneus.sum_bednorm_reads+p)/(d$Bl6_noRep4.sum_bednorm_reads+p))

p=1
d$log2FC.sum_bednorm_scores.Cast.vs.Bl6 = log2((d$Castaneus.sum_bednorm_scores+p)/(d$Bl6_noRep4.sum_bednorm_scores+p))

dd=d

ylimr=c(-1,1)
ylims=c(-3,3)

##All peaks regardless of interacting genes - are their interactions reduced?

par(mfcol=c(2,2))

##Bl6

l=list(d$log2FC.sum_bednorm_reads.Cast.vs.Bl6)
l=c(l,list(d[d$LFC_k27acPeak_oE < log2(1.15) & d$LFC_k27acPeak_oE > -log2(1.15),]$log2FC.sum_bednorm_reads.Cast.vs.Bl6))
ts=c(1,1.3,1.5,2,2.5,3)
for (tsh in ts){
  l=c(l,list(d[d$LFC_k27acPeak_oE <= -log2(tsh),]$log2FC.sum_bednorm_reads.Cast.vs.Bl6))
}

v=sapply(l,FUN=function(x) length(x))
names(l)=paste0(c('all','stab',ts),',n=',v)

boxplot(l,notch=TRUE,outline=FALSE,ylab='log2FC CapC Cast/Bl6 reads', xlab='FC ChIP-cutoff K27ac peak Bl6',ylim=ylimr)
abline(h=0,lty=2,col='red')


l=list(d$log2FC.sum_bednorm_scores.Cast.vs.Bl6)
l=c(l,list(d[d$LFC_k27acPeak_oE < log2(1.15) & d$LFC_k27acPeak_oE > -log2(1.15),]$log2FC.sum_bednorm_scores.Cast.vs.Bl6))
ts=c(1,1.3,1.5,2,2.5,3)
for (tsh in ts){
  l=c(l,list(d[d$LFC_k27acPeak_oE <= -log2(tsh),]$log2FC.sum_bednorm_scores.Cast.vs.Bl6))
}

v=sapply(l,FUN=function(x) length(x))
names(l)=paste0(c('all','stab',ts),',n=',v)

boxplot(l,notch=TRUE,outline=FALSE,ylab='log2FC CapC Cast/Bl6 scores', xlab='FC ChIP-cutoff K27ac peak Bl6',ylim=ylims)
abline(h=0,lty=2,col='red')


##Cast-specific peaks:

l=list(d$log2FC.sum_bednorm_reads.Cast.vs.Bl6)
l=c(l,list(d[d$LFC_k27acPeak_oE < log2(1.15) & d$LFC_k27acPeak_oE > -log2(1.15),]$log2FC.sum_bednorm_reads.Cast.vs.Bl6))
ts=c(1,1.3,1.5,2,2.5,3)
for (tsh in ts){
  l=c(l,list(d[d$LFC_k27acPeak_oE >= log2(tsh),]$log2FC.sum_bednorm_reads.Cast.vs.Bl6))
}

v=sapply(l,FUN=function(x) length(x))
names(l)=paste0(c('all','stab',ts),',n=',v)

boxplot(l, notch=TRUE,outline=FALSE, ylab='log2FC CapC Cast/Bl6 reads', xlab='FC ChIP-cutoff K27ac peak Cast',ylim=ylimr)
abline(h=0,lty=2,col='red')



l=list(d$log2FC.sum_bednorm_scores.Cast.vs.Bl6)
l=c(l,list(d[d$LFC_k27acPeak_oE < log2(1.15) & d$LFC_k27acPeak_oE > -log2(1.15),]$log2FC.sum_bednorm_scores.Cast.vs.Bl6))
ts=c(1,1.3,1.5,2,2.5,3)
for (tsh in ts){
  l=c(l,list(d[d$LFC_k27acPeak_oE >= log2(tsh),]$log2FC.sum_bednorm_scores.Cast.vs.Bl6))
}


v=sapply(l,FUN=function(x) length(x))
names(l)=paste0(c('all','stab',ts),',n=',v)


boxplot(l, notch=TRUE,outline=FALSE, ylab='log2FC CapC Cast/Bl6 scores', xlab='FC ChIP-cutoff K27ac peak Cast',ylim=ylims)
abline(h=0,lty=2,col='red')



##Look at peak interactions dependent on the genes these peaks interact:

for (gt in unique(dd$genetype)){

  main=paste(gt,'genes')
  
  d=dd[dd$genetype==gt,]
  
  par(mfcol=c(2,2))
  
  ##Bl6
  
  l=list(d$log2FC.sum_bednorm_reads.Cast.vs.Bl6)
  l=c(l,list(d[d$LFC_k27acPeak_oE < log2(1.15) & d$LFC_k27acPeak_oE > -log2(1.15),]$log2FC.sum_bednorm_reads.Cast.vs.Bl6))
  ts=c(1,1.3,1.5,2,2.5,3)
  for (tsh in ts){
    l=c(l,list(d[d$LFC_k27acPeak_oE <= -log2(tsh),]$log2FC.sum_bednorm_reads.Cast.vs.Bl6))
  }
  
  v=sapply(l,FUN=function(x) length(x))
  names(l)=paste0(c('all','stab',ts),',n=',v)

  
  boxplot(l,notch=TRUE,outline=FALSE,ylab='log2FC CapC Cast/Bl6 reads', xlab='FC ChIP-cutoff K27ac peak Bl6', main=main,ylim=ylimr)
  abline(h=0,lty=2,col='red')
  
  
  l=list(d$log2FC.sum_bednorm_scores.Cast.vs.Bl6)
  l=c(l,list(d[d$LFC_k27acPeak_oE < log2(1.15) & d$LFC_k27acPeak_oE > -log2(1.15),]$log2FC.sum_bednorm_scores.Cast.vs.Bl6))
  ts=c(1,1.3,1.5,2,2.5,3)
  for (tsh in ts){
    l=c(l,list(d[d$LFC_k27acPeak_oE <= -log2(tsh),]$log2FC.sum_bednorm_scores.Cast.vs.Bl6))
  }
  
  v=sapply(l,FUN=function(x) length(x))
  names(l)=paste0(c('all','stab',ts),',n=',v)

  
  boxplot(l,notch=TRUE,outline=FALSE,ylab='log2FC CapC Cast/Bl6 scores', xlab='FC ChIP-cutoff K27ac peak Bl6', main=main,ylim=ylims)
  abline(h=0,lty=2,col='red')
  
  
  ##Cast-specific peaks:
  
  l=list(d$log2FC.sum_bednorm_reads.Cast.vs.Bl6)
  l=c(l,list(d[d$LFC_k27acPeak_oE < log2(1.15) & d$LFC_k27acPeak_oE > -log2(1.15),]$log2FC.sum_bednorm_reads.Cast.vs.Bl6))
  ts=c(1,1.3,1.5,2,2.5,3)
  for (tsh in ts){
    l=c(l,list(d[d$LFC_k27acPeak_oE >= log2(tsh),]$log2FC.sum_bednorm_reads.Cast.vs.Bl6))
  }
  
  v=sapply(l,FUN=function(x) length(x))
  names(l)=paste0(c('all','stab',ts),',n=',v)

  
  boxplot(l, notch=TRUE,outline=FALSE, ylab='log2FC CapC Cast/Bl6 reads', xlab='FC ChIP-cutoff K27ac peak Cast', main=main,ylim=ylimr)
  abline(h=0,lty=2,col='red')
  
  
  
  l=list(d$log2FC.sum_bednorm_scores.Cast.vs.Bl6)
  l=c(l,list(d[d$LFC_k27acPeak_oE < log2(1.15) & d$LFC_k27acPeak_oE > -log2(1.15),]$log2FC.sum_bednorm_scores.Cast.vs.Bl6))
  ts=c(1,1.3,1.5,2,2.5,3)
  for (tsh in ts){
    l=c(l,list(d[d$LFC_k27acPeak_oE >= log2(tsh),]$log2FC.sum_bednorm_scores.Cast.vs.Bl6))
  }
  
  v=sapply(l,FUN=function(x) length(x))
  names(l)=paste0(c('all','stab',ts),',n=',v)

  
  boxplot(l, notch=TRUE,outline=FALSE, ylab='log2FC CapC Cast/Bl6 scores', xlab='FC ChIP-cutoff K27ac peak Cast', main=main,ylim=ylims)
  abline(h=0,lty=2,col='red')

}

dev.off()

```


#########################################################
##19.11.2021
#########################################################

########################################################
##TFBS enrichments in the differential peaks
##and in differential peaks connected to different genes
########################################################

```{r}

# ##load gene expression data:
# d=read.csv(paste0(base,'analysis/MouseGenomes/rnaseq/masked_SNP.SV/ResultTables/refSeq.geneBodyUnion.ignore.strand.TRUE.inter.feature.FALSE.Peaks.Castaneus_vs_Bl6_annotated_q0.005_fe1.5.csv'),stringsAsFactors=FALSE)

d=genes

d$Bl6.pooled_log2fpkm=log2(d$Bl6.pooled_fpkm+0.01)
d$Castaneus.pooled_log2fpkm=log2(d$Castaneus.pooled_fpkm+0.01)

plot(density(d$Bl6.pooled_log2fpkm))
lines(density(d$Castaneus.pooled_log2fpkm))
cutoff=0
abline(v=0,lty=2,col='red')

d$expressed_Bl6=FALSE
d$expressed_Castaneus=FALSE
d[d$Castaneus.pooled_log2fpkm>0,]$expressed_Castaneus=TRUE
d[d$Bl6.pooled_log2fpkm>0,]$expressed_Bl6=TRUE

head(d[d$expressed_Bl6 ==FALSE  & d$expressed_Castaneus,])

# write.csv(d, 'RNAseq_q0.005_fe1.5_forTFBSanalysis.csv',row.names=FALSE,quote=FALSE)

##expressed genes:
expr=d[d$expressed_Bl6 | d$expressed_Castaneus,]$genenames
expr=d[d$genenames %in% baited_genes$genenames]$genenames

library(TFBSTools)
suppressMessages(library(JASPAR2018))
opts <- list()
opts[["species"]] <- 10090
opts[["species"]] <- 'Mus musculus'
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
PFMatrixList
#> PFMatrixList of length 1
#> names(1): MA0002.1

##access individual entries:
name(PFMatrixList[[100]]) ##shows gene name
name(PFMatrixList) ##shows all names

pfm <- getMatrixByName(JASPAR2018, name="Foxo1")
seqLogo(toICM(pfm))

##identify TFs in the matrix that are also expressed in Bl6 or Castaneus:
names=unique(expr[expr %in% name(PFMatrixList),])
motifs=PFMatrixList[name(PFMatrixList) %in% names]

##motif matching:
library(motifmatchr)
library(GenomicRanges)
library(SummarizedExperiment)
library(BSgenome)

# load some example motifs
data(example_motifs, package = "motifmatchr") 

# # Load a set of peaks (oe)
# ints=read.delim('../../CaptureC/MyData/CapC11/postchicago/RAtimecourse/ints_all.txt',stringsAsFactors=FALSE)
# peaks=bed2gr(ints[,c('seqnames_otherEnd','start_otherEnd','end_otherEnd')])

# Define sets of peaks from K27ac data

range=c('Chr','Start','End')

k27ac_all$connected=FALSE
ov=as.matrix(findOverlaps(bed2gr(k27ac_all[,range]), bed2gr(ints[,c('seqnames_otherEnd','start_otherEnd','end_otherEnd')])))
k27ac_all[ov[,1],]$connected=TRUE


k27ac_all$connected_Nanog=FALSE
ov=as.matrix(findOverlaps(bed2gr(k27ac_all[,range]), bed2gr(ints[ints$genename=='Nanog',c('seqnames_otherEnd','start_otherEnd','end_otherEnd')])))
k27ac_all[ov[,1],]$connected_Nanog=TRUE


# k27ac_all$connected_stableGenes=FALSE
# ov=as.matrix(findOverlaps(bed2gr(k27ac_all[,range]), bed2gr(ints[ints$genetype=='very_stable',c('seqnames_otherEnd','start_otherEnd','end_otherEnd')])))
# k27ac_all[ov[,1],]$connected_stableGenes=TRUE
# 
# k27ac_all$connected_Bl6Genes=FALSE
# ov=as.matrix(findOverlaps(bed2gr(k27ac_all[,range]), bed2gr(ints[ints$genetype=='down_strong',c('seqnames_otherEnd','start_otherEnd','end_otherEnd')])))
# k27ac_all[ov[,1],]$connected_Bl6Genes=TRUE
# 
# k27ac_all$connected_CastaneusGenes=FALSE
# ov=as.matrix(findOverlaps(bed2gr(k27ac_all[,range]), bed2gr(ints[ints$genetype=='up_strong',c('seqnames_otherEnd','start_otherEnd','end_otherEnd')])))
# k27ac_all[ov[,1],]$connected_CastaneusGenes=TRUE

##make a list of peaksets
l=list(
  # peaks=bed2gr(k27ac_all[,range]),
  # peaks.stable=bed2gr(k27ac_all[k27ac_all$LFC> -log2(1.2) & k27ac_all$LFC < log2(1.2),range]),
  # peaks.Bl6.1=bed2gr(k27ac_all[k27ac_all$LFC< 0,range]),
  # peaks.Bl6.2=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(1.5),range]),
  # peaks.Bl6.3=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(2),range]),
  # peaks.Bl6.4=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(3),range]),
  # peaks.Bl6.5=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(4),range]),
  # peaks.Castaneus.1=bed2gr(k27ac_all[k27ac_all$LFC> 0,range]),
  # peaks.Castaneus.2=bed2gr(k27ac_all[k27ac_all$LFC>= log2(1.5),range]),
  # peaks.Castaneus.3=bed2gr(k27ac_all[k27ac_all$LFC>= log2(2),range]),
  # peaks.Castaneus.4=bed2gr(k27ac_all[k27ac_all$LFC>= log2(3),range]),
  # peaks.Castaneus.5=bed2gr(k27ac_all[k27ac_all$LFC>= log2(4),range]),
  # peaks.shared=bed2gr(k27ac_all[k27ac_all$peaktype=='Bl6+Castaneus',range]),
  # peaks.Bl6.specific=bed2gr(k27ac_all[k27ac_all$peaktype=='Bl6',range]),
  # peaks.Castaneus.specific=bed2gr(k27ac_all[k27ac_all$peaktype=='Castaneus',range]),
  # 
  # peaks.connected=bed2gr(k27ac_all[k27ac_all$connected,range]),
  # peaks.connectedToStableGenes=bed2gr(k27ac_all[k27ac_all$connected_stableGenes,range]),
  # peaks.connectedToBl6Genes=bed2gr(k27ac_all[k27ac_all$connected_Bl6Genes,range]),
  # peaks.connectedToCastaneusGenes=bed2gr(k27ac_all[k27ac_all$connected_CastaneusGenes,range]),
  # 
  # peaks.connectedToStableGenes_Bl6.3=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(2) & k27ac_all$connected_stableGenes,range]),
  # peaks.connectedToBl6Genes_Bl6.3=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(2) & k27ac_all$connected_Bl6Genes,range]),
  # peaks.connectedToCastaneusGenes_Bl6.3=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(2) & k27ac_all$connected_CastaneusGenes,range]),
  # 
  # peaks.connectedToStableGenes_Cast.3=bed2gr(k27ac_all[k27ac_all$LFC>= log2(2) & k27ac_all$connected_stableGenes,range]),
  # peaks.connectedToBl6Genes_Cast.3=bed2gr(k27ac_all[k27ac_all$LFC>= log2(2) & k27ac_all$connected_Bl6Genes,range]),
  # peaks.connectedToCastaneusGenes_Cast.3=bed2gr(k27ac_all[k27ac_all$LFC>= log2(2) & k27ac_all$connected_CastaneusGenes,range])
  # 
  
  peaks.connected=bed2gr(k27ac_all[k27ac_all$connected,range]),
  peaks.connected_Nanog=bed2gr(k27ac_all[k27ac_all$connected_Nanog,range])
  
)

names(l)=paste0(names(l),'_n=',sapply(l,FUN=function(x) length(x)))

##calculate motif enrichment for each TFBS expressed in ESCs:
i=1

peaks=l[[i]]
  
# Get motif matches for example motifs in peaks 
motif_ix <- matchMotifs(motifs, peaks, genome = "mm10") 
res=motifMatches(motif_ix) # Extract matches matrix from result
res[1,1]
name(motifs)

##convert the matrix in enrichments of different TFBSs among all K27ac peaks:
r=as.matrix(res)
summary(names(name(motifs))==colnames(r))  ##checks that the names of the matrix columns are the same as motif names
colnames(r)=name(motifs)  ##assigns TF names to matrix r
f=function(x){
  return(table(factor(x,levels=c(TRUE,FALSE)))[1]/length(x))
  }
fractions=apply(r,2,FUN=f )

m=matrix(fractions,nrow=1)
colnames(m)=names(fractions)

for (i in 2:length(l)){

  peaks=l[[i]]
    
  # Get motif matches for example motifs in peaks 
  motif_ix <- matchMotifs(motifs, peaks, genome = "mm10") 
  res=motifMatches(motif_ix) # Extract matches matrix from result
  res[1,1]
  name(motifs)
  
  ##convert the matrix in enrichments of different TFBSs among all K27ac peaks:
  r=as.matrix(res)
  summary(names(name(motifs))==colnames(r))  ##checks that the names of the matrix columns are the same as motif names
  colnames(r)=name(motifs)  ##assigns TF names to matrix r
  f=function(x){
    return(table(factor(x,levels=c(TRUE,FALSE)))[1]/length(x))
    }
  fractions=apply(r,2,FUN=f )
  
  mat=matrix(fractions,nrow=1)
  m=rbind(m,mat)
}
rownames(m)=names(l)

pheatmap(m,cluster_rows = FALSE,file='pheatmap_pluripotency_TFBS_K27ac_all.pdf')


##Now sort m by LFC of the TF expression:
expr=vector()
for (tf in colnames(m)){
  print(paste(tf,d[d$genenames==tf,]$LFC,sep=':'))
  expr=c(expr,mean(d[d$genenames==tf,]$LFC))
}
names(expr)=colnames(m)
order(expr)

par(mfrow=c(1,2))
pheatmap(m[,order(expr)],cluster_rows = FALSE,cluster_cols=FALSE,file='pheatmap_pluripotency_TFBS_K27ac_all_sortedIncreasingTFexpr.pdf')
pdf('Barplot_increasing_TF_expr_sortingUsedForHeatmap.pdf')
par(mfrow=c(2,1))
barplot(expr[order(expr)],las=3,ylab='LFC Cast.vs.Bl6')
dev.off()

###############################################################
##smaller list of peaksets (only Bl6 and Cast-specific peaks):
###############################################################

##make a list of peaksets
l=list(
  peaks.stable=bed2gr(k27ac_all[k27ac_all$LFC> -log2(1.2) & k27ac_all$LFC < log2(1.2),range]),
  peaks.Bl6.1=bed2gr(k27ac_all[k27ac_all$LFC< 0,range]),
  peaks.Bl6.2=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(1.5),range]),
  peaks.Bl6.3=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(2),range]),
  peaks.Bl6.4=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(3),range]),
  peaks.Bl6.5=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(4),range]),
  peaks.Castaneus.1=bed2gr(k27ac_all[k27ac_all$LFC> 0,range]),
  peaks.Castaneus.2=bed2gr(k27ac_all[k27ac_all$LFC>= log2(1.5),range]),
  peaks.Castaneus.3=bed2gr(k27ac_all[k27ac_all$LFC>= log2(2),range]),
  peaks.Castaneus.4=bed2gr(k27ac_all[k27ac_all$LFC>= log2(3),range]),
  peaks.Castaneus.5=bed2gr(k27ac_all[k27ac_all$LFC>= log2(4),range])
)

names(l)=paste0(names(l),'_n=',sapply(l,FUN=function(x) length(x)))

##calculate motif enrichment for each TFBS expressed in ESCs:
i=1

peaks=l[[i]]
  
# Get motif matches for example motifs in peaks 
motif_ix <- matchMotifs(motifs, peaks, genome = "mm10") 
res=motifMatches(motif_ix) # Extract matches matrix from result
res[1,1]
name(motifs)

##convert the matrix in enrichments of different TFBSs among all K27ac peaks:
r=as.matrix(res)
summary(names(name(motifs))==colnames(r))  ##checks that the names of the matrix columns are the same as motif names
colnames(r)=name(motifs)  ##assigns TF names to matrix r
f=function(x){
  return(table(factor(x,levels=c(TRUE,FALSE)))[1]/length(x))
  }
fractions=apply(r,2,FUN=f )

m=matrix(fractions,nrow=1)
colnames(m)=names(fractions)

for (i in 2:length(l)){

  peaks=l[[i]]
    
  # Get motif matches for example motifs in peaks 
  motif_ix <- matchMotifs(motifs, peaks, genome = "mm10") 
  res=motifMatches(motif_ix) # Extract matches matrix from result
  res[1,1]
  name(motifs)
  
  ##convert the matrix in enrichments of different TFBSs among all K27ac peaks:
  r=as.matrix(res)
  summary(names(name(motifs))==colnames(r))  ##checks that the names of the matrix columns are the same as motif names
  colnames(r)=name(motifs)  ##assigns TF names to matrix r
  f=function(x){
    return(table(factor(x,levels=c(TRUE,FALSE)))[1]/length(x))
    }
  fractions=apply(r,2,FUN=f )
  
  mat=matrix(fractions,nrow=1)
  m=rbind(m,mat)
}
rownames(m)=names(l)

pheatmap(m,cluster_rows = FALSE,file='pheatmap_pluripotency_TFBS_K27ac_all_Bl6orCastK27acPeaksOnly.pdf',cellheight = 7,gaps_row=c(1,6))


# ##Now sort m by LFC of the TF expression:
# expr=vector()
# for (tf in colnames(m)){
#   print(paste(tf,d[d$genenames==tf,]$LFC,sep=':'))
#   expr=c(expr,mean(d[d$genenames==tf,]$LFC))
# }
# names(expr)=colnames(m)
# order(expr)
# 
# par(mfrow=c(1,2))
# pheatmap(m[,order(expr)],cluster_rows = FALSE,cluster_cols=FALSE,file='pheatmap_pluripotency_TFBS_K27ac_all_sortedIncreasingTFexpr.pdf')
# pdf('Barplot_increasing_TF_expr_sortingUsedForHeatmap.pdf')
# par(mfrow=c(2,1))
# barplot(expr[order(expr)],las=3,ylab='LFC Cast.vs.Bl6')
# dev.off()



################################################################
##List with only Bl6/Cast-connected peaks:
################################################################

##make a list of peaksets
l=list(
 
  peaks.connected=bed2gr(k27ac_all[k27ac_all$connected,range]),
  
  peaks.connectedToStableGenes_Bl6.3=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(2) & k27ac_all$connected_stableGenes,range]),
  peaks.connectedToBl6Genes_Bl6.3=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(2) & k27ac_all$connected_Bl6Genes,range]),
  peaks.connectedToCastaneusGenes_Bl6.3=bed2gr(k27ac_all[k27ac_all$LFC<= -log2(2) & k27ac_all$connected_CastaneusGenes,range]),
  
  peaks.connectedToStableGenes_Cast.3=bed2gr(k27ac_all[k27ac_all$LFC>= log2(2) & k27ac_all$connected_stableGenes,range]),
  peaks.connectedToBl6Genes_Cast.3=bed2gr(k27ac_all[k27ac_all$LFC>= log2(2) & k27ac_all$connected_Bl6Genes,range]),
  peaks.connectedToCastaneusGenes_Cast.3=bed2gr(k27ac_all[k27ac_all$LFC>= log2(2) & k27ac_all$connected_CastaneusGenes,range])
)

names(l)=paste0(names(l),'_n=',sapply(l,FUN=function(x) length(x)))

##calculate motif enrichment for each TFBS expressed in ESCs:
i=1

peaks=l[[i]]
  
# Get motif matches for example motifs in peaks 
motif_ix <- matchMotifs(motifs, peaks, genome = "mm10") 
res=motifMatches(motif_ix) # Extract matches matrix from result
res[1,1]
name(motifs)

##convert the matrix in enrichments of different TFBSs among all K27ac peaks:
r=as.matrix(res)
summary(names(name(motifs))==colnames(r))  ##checks that the names of the matrix columns are the same as motif names
colnames(r)=name(motifs)  ##assigns TF names to matrix r
f=function(x){
  return(table(factor(x,levels=c(TRUE,FALSE)))[1]/length(x))
  }
fractions=apply(r,2,FUN=f )

m=matrix(fractions,nrow=1)
colnames(m)=names(fractions)

for (i in 2:length(l)){

  peaks=l[[i]]
    
  # Get motif matches for example motifs in peaks 
  motif_ix <- matchMotifs(motifs, peaks, genome = "mm10") 
  res=motifMatches(motif_ix) # Extract matches matrix from result
  res[1,1]
  name(motifs)
  
  ##convert the matrix in enrichments of different TFBSs among all K27ac peaks:
  r=as.matrix(res)
  summary(names(name(motifs))==colnames(r))  ##checks that the names of the matrix columns are the same as motif names
  colnames(r)=name(motifs)  ##assigns TF names to matrix r
  f=function(x){
    return(table(factor(x,levels=c(TRUE,FALSE)))[1]/length(x))
    }
  fractions=apply(r,2,FUN=f )
  
  mat=matrix(fractions,nrow=1)
  m=rbind(m,mat)
}
rownames(m)=names(l)

pheatmap(m,cluster_rows = FALSE,file='pheatmap_pluripotency_TFBS_K27ac_all_Bl6orCastConnectedK27acPeaksOnly.pdf',cellheight = 7,gaps_row=c(1,4),cellwidth=7)


##Now sort m by LFC of the TF expression:
expr=vector()
for (tf in colnames(m)){
  print(paste(tf,d[d$genenames==tf,]$LFC,sep=':'))
  expr=c(expr,mean(d[d$genenames==tf,]$LFC))
}
names(expr)=colnames(m)
order(expr)

par(mfrow=c(1,2))
pheatmap(m[,order(expr)],cluster_rows = FALSE,cluster_cols=FALSE,file='pheatmap_pluripotency_TFBS_K27ac_all_sortedIncreasingTFexpr.pdf')
pdf('Barplot_increasing_TF_expr_sortingUsedForHeatmap.pdf')
par(mfrow=c(2,1))
barplot(expr[order(expr)],las=3,ylab='LFC Cast.vs.Bl6')
dev.off()




###############################################################
##Next step: Look at actual (ChIP-Seq) binding of TFs
##paper with peaksets:
##Galonska,...Meissner (Cell Stem Cell 2015),
##Ground State Conditions Induce Rapid
##Reorganization of Core Pluripotency Factor Binding
##before Global Epigenetic Reprogramming
##downloaded mm9 peak sets (3passages Sox2 and Nanog, 24h Oct4 on 2i)
##
##paper with KLF4 peaks:
##Campigli Di Giammartino,...,Apostolou (2019, Nature Cell Biology):
##KLF4 is involved in the organization and regulation of pluripotency-associated three-dimensional enhancer networks
##already in mm10
###############################################################

#####################
##Meissner peaks:
#####################

##convert Meissner peaksets into bed files that can be lifted over
dir='peaksets_Meissner/'
for (f in list.files(path=dir)){
  saveBed(read.delim(paste0(dir,f),header=FALSE,stringsAsFactors=FALSE)[,1:4],paste0(dir,f,'.bed'))
}

##manually lift Meissner peaksets over from mm9 to mm10 on UCSC LiftOver
peaksets=GRangesList()
for (f in list.files(path=dir,pattern='liftOverMm9ToMm10.bed$')){
  peaksets=c(peaksets,GRangesList(bed2gr((read.delim(paste0(dir,f),header=FALSE,stringsAsFactors=FALSE)))))
}
names(peaksets)=c('Nanog','Pou5f1','Sox2')

#######################
##Apostolou peaks:
#######################

##Now add KLF4 peaksets. We have two replicates - use peaks detected by both!
dir='peaksets_Apostolou/'
lp=GRangesList()
for (f in list.files(path=dir)){
  lp=c(lp,GRangesList(bed2gr((read.delim(paste0(dir,f),header=FALSE,stringsAsFactors=FALSE)))))
}
merged=unique(reduce(c(lp[[1]],lp[[2]])))
merged=merged[as.matrix(findOverlaps(merged,lp[[1]]))[,1]]
merged=unique(merged[as.matrix(findOverlaps(merged,lp[[2]]))[,1]])
peaksets=c(peaksets,GRangesList(merged))
names(peaksets)[4]='Klf4'

###########################################################################################
##look at all defined sets of k27ac peaks, how many of them bind the 4 pluripotency factors
###########################################################################################

##load gene expression data:
genes=read.csv(paste0(base,'analysis/MouseGenomes/rnaseq/masked_SNP.SV/ResultTables/refSeq.geneBodyUnion.ignore.strand.TRUE.inter.feature.FALSE.Peaks.Castaneus_vs_Bl6_annotated_q0.005_fe1.5.csv'),stringsAsFactors=FALSE)
tss=bed2gr(genes[,c('Chr','Start','End')])
strand(tss)=genes$Strand
tss=resize(getTSSgr(tss),2001,fix='center')

l.tss=lapply(l, FUN=function(x) unique(x[as.matrix(findOverlaps(x,tss))[,1]]))
l.notss=lapply(l, FUN=function(x) unique(x[-as.matrix(findOverlaps(x,tss))[,1]]))
l.all=l

# l=l.all
# l=l.tss
l=l.notss
names(l)=paste0(names(l),'_n=',sapply(l,FUN=function(x) length(x)))


for (i in 1:length(l)){
  print(i)
  gr=l[[i]]
  gr$Nanog=FALSE
  gr$Pou5f1=FALSE
  gr$Sox2=FALSE
  gr$Klf4=FALSE
  ov=as.matrix(findOverlaps(gr,peaksets['Nanog']))
  gr[ov[,1]]$Nanog=TRUE
  ov=as.matrix(findOverlaps(gr,peaksets['Pou5f1']))
  gr[ov[,1]]$Pou5f1=TRUE
  ov=as.matrix(findOverlaps(gr,peaksets['Sox2']))
  gr[ov[,1]]$Sox2=TRUE
  ov=as.matrix(findOverlaps(gr,peaksets['Klf4']))
  gr[ov[,1]]$Klf4=TRUE
  l[[i]]=gr
}

ff=function(x){
    return(table(factor(x,levels=c(TRUE,FALSE)))[1]/length(x))
    }
f=function(x){
  df=as.matrix(values(x))
  return(apply(df,2,FUN=ff ))
}
m=do.call('rbind',lapply(l,FUN=f))


##plots:

# pheatmap(m,cluster_rows = FALSE,file='pheatmap_Nanog.Oct4.Sox2.Klf4.ChIP.peaks_K27ac_all.pdf',main = 'All K27ac peaks')
# pheatmap(m,cluster_rows = FALSE,file='pheatmap_Nanog.Oct4.Sox2.Klf4.ChIP.peaks_K27ac_all_TSS1kb.pdf',main='K27ac peaks within 1kb of a TSS')
# pheatmap(m,cluster_rows = FALSE,file='pheatmap_Nanog.Oct4.Sox2.Klf4.ChIP.peaks_K27ac_all_noTSS1kb.pdf',main='K27ac peaks > 1kb from a TSS')
# 

##Now sort m by LFC of the TF expression:
expr=vector()
d=genes
for (tf in colnames(m)){
  print(paste(tf,d[d$genenames==tf,]$LFC,sep=':'))
  expr=c(expr,mean(d[d$genenames==tf,]$LFC))
}
names(expr)=colnames(m)
order(expr)


# pheatmap(m[,order(expr)],cluster_rows = FALSE,cluster_cols=FALSE,file='pheatmap_Nanog.Oct4.Sox2.Klf4.ChIP.peaks_K27ac_all_sortedIncreasingTFexpr.pdf',main='All K27ac peaks')
# pheatmap(m[,order(expr)],cluster_rows = FALSE,cluster_cols=FALSE,file='pheatmap_Nanog.Oct4.Sox2.Klf4.ChIP.peaks_K27ac_all_TSS1kb_sortedIncreasingTFexpr.pdf',main='K27ac peaks within 1kb of a TSS')
pheatmap(m[,order(expr)],cluster_rows = FALSE,cluster_cols=FALSE,file='pheatmap_Nanog.Oct4.Sox2.Klf4.ChIP.peaks_K27ac_all_noTSS1kb_sortedIncreasingTFexpr.pdf',main='K27ac peaks > 1kb from a TSS')

pdf('Barplot_increasing_TF_expr_Nanog.Oct4.Sox2.Klf4.ChIP.peaks_sortingUsedForHeatmap.pdf')
barplot(expr[order(expr)],ylab='LFC Cast.vs.Bl6')
dev.off()


##Conclusion: It looks like Nanog expression is reduced an this correlates with a loss of K27ac.



##############################################################################################
##Now: See how interactions with k27ac peaks occupied or not occupied by Nanog are regulated
##############################################################################################

l=l.all
pdf('Boxplots_intLFC.all_plminNanog.pdf',onefile=TRUE)

l=l.tss
pdf('Boxplots_intLFC.1kbTSS_plminNanog.pdf',onefile=TRUE)

l=l.notss
pdf('Boxplots_intLFC.no1kbTSS_plminNanog.pdf',onefile=TRUE)

for (i in 1:length(l)){
  gr=l[[i]]
  gr$Nanog=FALSE
  gr$Pou5f1=FALSE
  gr$Sox2=FALSE
  gr$Klf4=FALSE
  ov=as.matrix(findOverlaps(gr,peaksets['Nanog']))
  gr[ov[,1]]$Nanog=TRUE
  ov=as.matrix(findOverlaps(gr,peaksets['Pou5f1']))
  gr[ov[,1]]$Pou5f1=TRUE
  ov=as.matrix(findOverlaps(gr,peaksets['Sox2']))
  gr[ov[,1]]$Sox2=TRUE
  ov=as.matrix(findOverlaps(gr,peaksets['Klf4']))
  gr[ov[,1]]$Klf4=TRUE
  l[[i]]=gr
}

nanog=l[[1]][l[[1]]$Nanog]
nonanog=l[[1]][l[[1]]$Nanog==FALSE]

i=2
dd=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_interacting_readsScoresInPeak.txt'),stringsAsFactors=FALSE)
dd$type_k27ac_otherEnd='any'
dd[dd$LFC_k27acPeak_oE <= -log2(2.5),]$type_k27ac_otherEnd='K27ac_Bl6'
dd[dd$LFC_k27acPeak_oE >= log2(2.5),]$type_k27ac_otherEnd='K27ac_Castaneus'
dd[dd$LFC_k27acPeak_oE <= log2(1.15) & dd$LFC_k27acPeak_oE >= -log2(1.15),]$type_k27ac_otherEnd='stable'

for (peaktype in unique(dd$type_k27ac_otherEnd)){
  
  d=dd[dd$type_k27ac_otherEnd==peaktype,]
  oe=bed2gr(d[,c('seqnames_otherEnd','start_otherEnd','end_otherEnd')])
  
  d$Nanog=FALSE
  ov=as.matrix(findOverlaps(oe,nanog))
  if (nrow(ov)>0){
    d[ov[,1],]$Nanog=TRUE
  }
  
  p=8
  d$LFC_sum_bednorm_reads=log2((d$Castaneus.sum_bednorm_reads+p)/(d$Bl6.sum_bednorm_reads+p))
  p=1
  d$LFC_sum_bednorm_scores=log2((d$Castaneus.sum_bednorm_scores+p)/(d$Bl6.sum_bednorm_scores+p))
  
  par(mfrow=c(2,2))
  
  ylim=c(-0.8,0.8)
  xlab='Nanog binding at oE'
  ylab='LFC sum bednorm reads Cast.vs.Bl6'
  boxplot(d$LFC_sum_bednorm_reads~d$Nanog,outline=FALSE,main=paste('all genes; peaktype:',peaktype),ylim=ylim,xlab=xlab,ylab=ylab)
  abline(h=0,lty=2,col='red')
  for (cl in unique(d$genetype)){
   boxplot(d[d$genetype==cl,]$LFC_sum_bednorm_reads~d[d$genetype==cl,]$Nanog,outline=FALSE,main=paste(cl,'genes; peaktype:', peaktype),ylim=ylim,xlab=xlab,ylab=ylab) 
    abline(h=0,lty=2,col='red')
  }
  
  ylim=c(-2,2)
  ylab='LFC sum bednorm scores Cast.vs.Bl6'
  boxplot(d$LFC_sum_bednorm_scores~d$Nanog,outline=FALSE,main=paste('all genes; peaktype:',peaktype),ylim=ylim,xlab=xlab,ylab=ylab)
  abline(h=0,lty=2,col='red')
  for (cl in unique(d$genetype)){
   boxplot(d[d$genetype==cl,]$LFC_sum_bednorm_scores~d[d$genetype==cl,]$Nanog,outline=FALSE,main=paste(cl,'genes; peaktype:', peaktype),ylim=ylim,xlab=xlab,ylab=ylab) 
    abline(h=0,lty=2,col='red')
  }

}

dev.off()



##############################################################################################
##Now: See how interactions with k27ac peaks occupied or not occupied by KLF4 are regulated
##############################################################################################

l=l.all
pdf('Boxplots_intLFC.all_plminKlf4.pdf',onefile=TRUE)

l=l.tss
pdf('Boxplots_intLFC.1kbTSS_plminKlf4.pdf',onefile=TRUE)

l=l.notss
pdf('Boxplots_intLFC.no1kbTSS_plminKlf4.pdf',onefile=TRUE)

for (i in 1:length(l)){
  gr=l[[i]]
  gr$Nanog=FALSE
  gr$Pou5f1=FALSE
  gr$Sox2=FALSE
  gr$Klf4=FALSE
  ov=as.matrix(findOverlaps(gr,peaksets['Nanog']))
  gr[ov[,1]]$Nanog=TRUE
  ov=as.matrix(findOverlaps(gr,peaksets['Pou5f1']))
  gr[ov[,1]]$Pou5f1=TRUE
  ov=as.matrix(findOverlaps(gr,peaksets['Sox2']))
  gr[ov[,1]]$Sox2=TRUE
  ov=as.matrix(findOverlaps(gr,peaksets['Klf4']))
  gr[ov[,1]]$Klf4=TRUE
  l[[i]]=gr
}

nanog=l[[1]][l[[1]]$Klf4]
nonanog=l[[1]][l[[1]]$Klf4==FALSE]

i=2
dd=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_interacting_readsScoresInPeak.txt'),stringsAsFactors=FALSE)
dd$type_k27ac_otherEnd='any'
dd[dd$LFC_k27acPeak_oE <= -log2(2.5),]$type_k27ac_otherEnd='K27ac_Bl6'
dd[dd$LFC_k27acPeak_oE >= log2(2.5),]$type_k27ac_otherEnd='K27ac_Castaneus'
dd[dd$LFC_k27acPeak_oE <= log2(1.15) & dd$LFC_k27acPeak_oE >= -log2(1.15),]$type_k27ac_otherEnd='stable'

for (peaktype in unique(dd$type_k27ac_otherEnd)){
  
  d=dd[dd$type_k27ac_otherEnd==peaktype,]
  oe=bed2gr(d[,c('seqnames_otherEnd','start_otherEnd','end_otherEnd')])
  
  d$Nanog=FALSE
  ov=as.matrix(findOverlaps(oe,nanog))
  if (nrow(ov)>0){
    d[ov[,1],]$Nanog=TRUE
  }
  
  p=8
  d$LFC_sum_bednorm_reads=log2((d$Castaneus.sum_bednorm_reads+p)/(d$Bl6.sum_bednorm_reads+p))
  p=1
  d$LFC_sum_bednorm_scores=log2((d$Castaneus.sum_bednorm_scores+p)/(d$Bl6.sum_bednorm_scores+p))
  
  par(mfrow=c(2,2))
  
  ylim=c(-0.8,0.8)
  xlab='Klf4 binding at oE'
  ylab='LFC sum bednorm reads Cast.vs.Bl6'
  boxplot(d$LFC_sum_bednorm_reads~d$Nanog,outline=FALSE,main=paste('all genes; peaktype:',peaktype),ylim=ylim,xlab=xlab,ylab=ylab)
  abline(h=0,lty=2,col='red')
  for (cl in unique(d$genetype)){
   boxplot(d[d$genetype==cl,]$LFC_sum_bednorm_reads~d[d$genetype==cl,]$Nanog,outline=FALSE,main=paste(cl,'genes; peaktype:', peaktype),ylim=ylim,xlab=xlab,ylab=ylab) 
    abline(h=0,lty=2,col='red')
  }
  
  ylim=c(-2,2)
  ylab='LFC sum bednorm scores Cast.vs.Bl6'
  boxplot(d$LFC_sum_bednorm_scores~d$Nanog,outline=FALSE,main=paste('all genes; peaktype:',peaktype),ylim=ylim,xlab=xlab,ylab=ylab)
  abline(h=0,lty=2,col='red')
  for (cl in unique(d$genetype)){
   boxplot(d[d$genetype==cl,]$LFC_sum_bednorm_scores~d[d$genetype==cl,]$Nanog,outline=FALSE,main=paste(cl,'genes; peaktype:', peaktype),ylim=ylim,xlab=xlab,ylab=ylab) 
    abline(h=0,lty=2,col='red')
  }

}

dev.off()



###################################################################
##General association between Oct4, Sox2, Klf4 or Nanog at otherEnd
##and regulation of interactions
###################################################################

l=l.all

for (i in 1:length(l)){
  gr=l[[i]]
  gr$Nanog=FALSE
  gr$Pou5f1=FALSE
  gr$Sox2=FALSE
  gr$Klf4=FALSE
  ov=as.matrix(findOverlaps(gr,peaksets['Nanog']))
  gr[ov[,1]]$Nanog=TRUE
  ov=as.matrix(findOverlaps(gr,peaksets['Pou5f1']))
  gr[ov[,1]]$Pou5f1=TRUE
  ov=as.matrix(findOverlaps(gr,peaksets['Sox2']))
  gr[ov[,1]]$Sox2=TRUE
  ov=as.matrix(findOverlaps(gr,peaksets['Klf4']))
  gr[ov[,1]]$Klf4=TRUE
  l[[i]]=gr
}

i=2
dd=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_interacting_readsScoresInPeak.txt'),stringsAsFactors=FALSE)
dd=annotateInts(dd,peaksets)

d=dd
p=8
d$LFC_sum_bednorm_reads=log2((d$Castaneus.sum_bednorm_reads+p)/(d$Bl6.sum_bednorm_reads+p))
p=1
d$LFC_sum_bednorm_scores=log2((d$Castaneus.sum_bednorm_scores+p)/(d$Bl6.sum_bednorm_scores+p))

ylimr=c(-1.3,1.3)
ylims=c(-4.5,4.5)

par(mfcol=c(2,4))
for (n in names(peaksets)){
  xlab=paste0(n,'_bait.',n,'_oE')
  boxplot(d$LFC_sum_bednorm_reads~d[,paste0(n,'_bait')]*d[,paste0(n,'_otherEnd')],data=d,outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_reads Cast.vs.Bl6',notch=TRUE,main=paste(n),las=3,ylim=ylimr)
  abline(h=0,lty=2,col='red')
  boxplot(d$LFC_sum_bednorm_scores~d[,paste0(n,'_bait')]*d[,paste0(n,'_otherEnd')],data=d,outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_scores Cast.vs.Bl6',notch=TRUE,main=paste(n),las=3,ylim=ylims)
  abline(h=0,lty=2,col='red')
}

##for Bl6-specific k27ac peaks oE
par(mfcol=c(2,4))
for (n in names(peaksets)){
  xlab=paste0(n,'_bait.',n,'_oE')
  a=d[d$LFC_k27acPeak_oE <= -log2(2),]
  boxplot(a$LFC_sum_bednorm_reads~a[,paste0(n,'_bait')]*a[,paste0(n,'_otherEnd')],outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_reads Cast.vs.Bl6',notch=TRUE,main=paste(n,'Bl6_k27ac_peak_oE'),las=3,ylim=ylimr)
  abline(h=0,lty=2,col='red')
  boxplot(a$LFC_sum_bednorm_scores~a[,paste0(n,'_bait')]*a[,paste0(n,'_otherEnd')],outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_scores Cast.vs.Bl6',notch=TRUE,main=paste(n,'Bl6_k27ac_peak_oE'),las=3,ylim=ylims)
  abline(h=0,lty=2,col='red')
}

##for Cast-specific k27ac peaks oE
par(mfcol=c(2,4))
for (n in names(peaksets)){
  xlab=paste0(n,'_bait.',n,'_oE')
  a=d[d$LFC_k27acPeak_oE >= log2(2),]
  boxplot(a$LFC_sum_bednorm_reads~a[,paste0(n,'_bait')]*a[,paste0(n,'_otherEnd')],outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_reads Cast.vs.Bl6',notch=TRUE,main=paste(n,'Cast_k27ac_peak_oE'),las=3,ylim=ylimr)
  abline(h=0,lty=2,col='red')
  boxplot(a$LFC_sum_bednorm_scores~a[,paste0(n,'_bait')]*a[,paste0(n,'_otherEnd')],outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_scores Cast.vs.Bl6',notch=TRUE,main=paste(n,'Cast_k27ac_peak_oE'),las=3,ylim=ylims)
  abline(h=0,lty=2,col='red')
}

##for stable k27ac peaks oE
par(mfcol=c(2,4))
for (n in names(peaksets)){
  xlab=paste0(n,'_bait.',n,'_oE')
  a=d[d$LFC_k27acPeak_oE >= -log2(1.15) & d$LFC_k27acPeak_oE <= log2(1.15),]
  boxplot(a$LFC_sum_bednorm_reads~a[,paste0(n,'_bait')]*a[,paste0(n,'_otherEnd')],outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_reads Cast.vs.Bl6',notch=TRUE,main=paste(n,'Stable_k27ac_peak_oE'),las=3,ylim=ylimr)
  abline(h=0,lty=2,col='red')
  boxplot(a$LFC_sum_bednorm_scores~a[,paste0(n,'_bait')]*a[,paste0(n,'_otherEnd')],outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_scores Cast.vs.Bl6',notch=TRUE,main=paste(n,'Stable_k27ac_peak_oE'),las=3,ylim=ylims)
  abline(h=0,lty=2,col='red')
}

##for Bl6-specific k27ac peaks oE interacting with down_strong genes
par(mfcol=c(2,4))
for (n in names(peaksets)){
  xlab=paste0(n,'_bait.',n,'_oE')
  a=d[d$LFC_k27acPeak_oE <= -log2(2) & d$genetype=='down_strong',]
  boxplot(a$LFC_sum_bednorm_reads~a[,paste0(n,'_bait')]*a[,paste0(n,'_otherEnd')],outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_reads Cast.vs.Bl6',notch=TRUE,main=paste(n,'Bl6gene;Bl6_k27ac_peak_oE'),las=3,ylim=ylimr)
  abline(h=0,lty=2,col='red')
  boxplot(a$LFC_sum_bednorm_scores~a[,paste0(n,'_bait')]*a[,paste0(n,'_otherEnd')],outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_scores Cast.vs.Bl6',notch=TRUE,main=paste(n,'Bl6gene;Bl6_k27ac_peak_oE'),las=3,ylim=ylims)
  abline(h=0,lty=2,col='red')
}

##for Bl6-specific k27ac peaks oE interacting with up_strong genes
par(mfcol=c(2,4))
for (n in names(peaksets)){
  xlab=paste0(n,'_bait.',n,'_oE')
  a=d[d$LFC_k27acPeak_oE <= -log2(2) & d$genetype=='up_strong',]
  boxplot(a$LFC_sum_bednorm_reads~a[,paste0(n,'_bait')]*a[,paste0(n,'_otherEnd')],outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_reads Cast.vs.Bl6',notch=TRUE,main=paste(n,'Castgene;Bl6_k27ac_peak_oE'),las=3,ylim=ylimr)
  abline(h=0,lty=2,col='red')
  boxplot(a$LFC_sum_bednorm_scores~a[,paste0(n,'_bait')]*a[,paste0(n,'_otherEnd')],outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_scores Cast.vs.Bl6',notch=TRUE,main=paste(n,'Castgene;Bl6_k27ac_peak_oE'),las=3,ylim=ylims)
  abline(h=0,lty=2,col='red')
}

##for Bl6-specific k27ac peaks oE interacting with stable genes
par(mfcol=c(2,4))
for (n in names(peaksets)){
  xlab=paste0(n,'_bait.',n,'_oE')
  a=d[d$LFC_k27acPeak_oE <= -log2(2) & d$genetype=='very_stable',]
  boxplot(a$LFC_sum_bednorm_reads~a[,paste0(n,'_bait')]*a[,paste0(n,'_otherEnd')],outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_reads Cast.vs.Bl6',notch=TRUE,main=paste(n,'Stabgene;Bl6_k27ac_peak_oE'),las=3,ylim=ylimr)
  abline(h=0,lty=2,col='red')
  boxplot(a$LFC_sum_bednorm_scores~a[,paste0(n,'_bait')]*a[,paste0(n,'_otherEnd')],outline=FALSE,xlab=xlab,ylab='LFC_sum_bednorm_scores Cast.vs.Bl6',notch=TRUE,main=paste(n,'Stabgene;Bl6_k27ac_peak_oE'),las=3,ylim=ylims)
  abline(h=0,lty=2,col='red')
}

 par(mfrow=c(2,2))
for (n in names(peaksets)){
  d[,paste(n,'bait.oE',sep='_')]=paste(d[,paste0(n,'_bait')],d[,paste0(n,'_otherEnd')],sep='.')
  tab=table(d[,paste(n,'bait.oE',sep='_')],d$genetype)
  plotNormBarplot(tab,ylab=paste('% ints w', n,'bait.oE'),main = n,ncol.legend=2,legend.title = paste(n,'bait.oE'),horiz=FALSE,xlab='genetype',col=c('grey','yellow','blue','green'),cex.legend = 0.75)

}

  par(mfrow=c(2,2))
for (n in names(peaksets)){
  peaktype='Bl6_specific_oE'
  d[,paste(n,'bait.oE',sep='_')]=paste(d[,paste0(n,'_bait')],d[,paste0(n,'_otherEnd')],sep='.')
  tab=table(d[d$LFC_k27acPeak_oE<= -log2(2),paste(n,'bait.oE',sep='_')],d[d$LFC_k27acPeak_oE<= -log2(2),]$genetype)
  plotNormBarplot(tab,ylab=paste('% ints w', n,'bait.oE'),main = paste(n,'peaktype:',peaktype),ncol.legend=2,legend.title = paste(n,'bait.oE'),horiz=FALSE,xlab='genetype',col=c('grey','yellow','blue','green'),cex.legend = 0.75)

}

 par(mfrow=c(2,2))
for (n in names(peaksets)){
  peaktype='Cast_specific_oE'
  d[,paste(n,'bait.oE',sep='_')]=paste(d[,paste0(n,'_bait')],d[,paste0(n,'_otherEnd')],sep='.')
  tab=table(d[d$LFC_k27acPeak_oE >= log2(2),paste(n,'bait.oE',sep='_')],d[d$LFC_k27acPeak_oE >= log2(2),]$genetype)
  plotNormBarplot(tab,ylab=paste('% ints w', n,'bait.oE'),main = paste(n,'peaktype:',peaktype),ncol.legend=2,legend.title = paste(n,'bait.oE'),horiz=FALSE,xlab='genetype',col=c('grey','yellow','blue','green'),cex.legend = 0.75)

}

 
 
 
 

```


####################
##04.01.2022
####################

#######################################################################################
##Testing the hypothesis that overall effect of distal enhancers
##varies for different types of promoters
##
##This means:
##(a) Bl6-specific genes either engage with Bl6-specific enhancers or  
##fail to engage with active enhancers in Castaneus cells and vice versa
##(b) Stable genes engage with active enhancers in every cell line to a similar degree
#######################################################################################

Here we calculate 

1. The sum of K27ac, interaction-reads and interaction-scores deltas for each K27ac-peak a gene interacts with (lfc.k27, lfc.reads,lfc.scores)
2. The sum of normalized (to the top min and max LFC) K27ac, interaction-reads and interaction-scores deltas for each K27ac-peak a gene interacts with (lfc.k27.norm, lfc.reads.norm,lfc.scores.norm)
4. SD of the sums, to test whether for each lost interaction stable genes gain another interaction
3. according to the ABC model of enahncer influence: Combined sum (pure addition) of lfc.k27.norm and lfc.reads.norm/lfc.scores.norm - this should test whether for stable genes this amounts to 0, suggesting that they are compensating for loss of activity or engagement of an enhancer


```{r}

pdf('Regulation_of_K27ac_and_interactions_perBait.pdf',onefile=TRUE)

# ##only once:
# i=2
# d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_interacting_readsScoresInPeak.txt'),stringsAsFactors=FALSE)
# 
# baited_genes$lfc.k27=0
# baited_genes$lfc.reads=0
# baited_genes$lfc.scores=0
# 
# for (n in baited_genes$genenames){
#   print(n)
# 
#   dd=d[d$genenames==n,]
#   baited_genes[baited_genes$genenames==n]$lfc.k27=sum(dd$LFC_k27acPeak_oE)
#   p=8
#   baited_genes[baited_genes$genenames==n]$lfc.reads=sum(log2((dd$Castaneus.sum_bednorm_reads+p)/(dd$Bl6.sum_bednorm_reads+p)))
#   p=1
#   baited_genes[baited_genes$genenames==n]$lfc.scores=sum(log2((dd$Castaneus.sum_bednorm_scores+p)/(dd$Bl6.sum_bednorm_scores+p)))
# }

##Average regulation of connected peaks per promoter:
# ##only once:
# i=2
# d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_interacting_readsScoresInPeak.txt'),stringsAsFactors=FALSE)
# 
# baited_genes$lfc.k27_mean=0
# baited_genes$lfc.reads_mean=0
# baited_genes$lfc.scores_mean=0
# 
# for (n in baited_genes$genenames){
#   print(n)
# 
#   dd=d[d$genenames==n,]
#   baited_genes[baited_genes$genenames==n]$lfc.k27_mean=mean(dd$LFC_k27acPeak_oE)
#   p=8
#   baited_genes[baited_genes$genenames==n]$lfc.reads_mean=mean(log2((dd$Castaneus.sum_bednorm_reads+p)/(dd$Bl6.sum_bednorm_reads+p)))
#   p=1
#   baited_genes[baited_genes$genenames==n]$lfc.scores_mean=mean(log2((dd$Castaneus.sum_bednorm_scores+p)/(dd$Bl6.sum_bednorm_scores+p)))
# }


##Average regulation of peaks surrounding promoters per promoter:
# ##only once:
# 
# baited_genes$lfc.k27_mean_peaks_within_100kb=0
# baited_genes$lfc.k27_mean_peaks_within_250kb=0
# baited_genes$lfc.k27_mean_peaks_within_500kb=0
# baited_genes$lfc.k27_mean_peaks_within_1Mb=0
# 
# gr=bed2gr(k27ac_all[,c('Chr','Start','End')])
# 
# for (n in baited_genes$genenames){
#   print(n)
#   
#   a=resize(baited_genes[baited_genes$genenames==n],width=200001,fix='center')
#   k=k27ac_all[as.matrix(findOverlaps(gr,a))[,1],]
#   dim(k)
#   baited_genes[baited_genes$genenames==n]$lfc.k27_mean_peaks_within_100kb=mean(k$LFC)
#   
#   a=resize(baited_genes[baited_genes$genenames==n],width=500001,fix='center')
#   k=k27ac_all[as.matrix(findOverlaps(gr,a))[,1],]
#   dim(k)
#   baited_genes[baited_genes$genenames==n]$lfc.k27_mean_peaks_within_250kb=mean(k$LFC)
#   
#   a=resize(baited_genes[baited_genes$genenames==n],width=1000001,fix='center')
#   k=k27ac_all[as.matrix(findOverlaps(gr,a))[,1],]
#   dim(k)
#   baited_genes[baited_genes$genenames==n]$lfc.k27_mean_peaks_within_500kb=mean(k$LFC)
#   
#   a=resize(baited_genes[baited_genes$genenames==n],width=2000001,fix='center')
#   k=k27ac_all[as.matrix(findOverlaps(gr,a))[,1],]
#   dim(k)
#   baited_genes[baited_genes$genenames==n]$lfc.k27_mean_peaks_within_1Mb=mean(k$LFC)
#   
# }

genetypes=unique(baited_genes$genetype)
cols=c('cornflowerblue','orange','green','red')

par(mfcol=c(2,4))

for (i in 1:length(genetypes)){

  gt=genetypes[i]
  col=cols[i]
  
  plot(baited_genes$lfc.k27,baited_genes$lfc.scores,main=gt)
  points(baited_genes[baited_genes$genetype==gt]$lfc.k27,baited_genes[baited_genes$genetype==gt]$lfc.scores,col=col,pch=20)
  abline(h=0,lty=2,col='grey')
  abline(v=0,lty=2,col='grey')
  
  plot(baited_genes$lfc.k27,baited_genes$lfc.reads,main=gt)
  points(baited_genes[baited_genes$genetype==gt]$lfc.k27,baited_genes[baited_genes$genetype==gt]$lfc.reads,col=col,pch=20)
  abline(h=0,lty=2,col='grey')
  abline(v=0,lty=2,col='grey')
}


##normalize log2 by adjusting the dynamic range of k27ac and interactions

# ##only once:
# d$LFC_k27acPeak_oE.norm=0
# p=8
# d$LFC.sum_bednorm_reads=log2((d$Castaneus.sum_bednorm_reads+p)/(d$Bl6.sum_bednorm_reads+p))
# p=1
# d$LFC.sum_bednorm_scores=log2((d$Castaneus.sum_bednorm_scores+p)/(d$Bl6.sum_bednorm_scores+p))
# d$LFC.sum_bednorm_reads.norm=0
# d$LFC.sum_bednorm_scores.norm=0
# 
# min=min(d$LFC_k27acPeak_oE)
# max=max(d$LFC_k27acPeak_oE)
# d[d$LFC_k27acPeak_oE<0,]$LFC_k27acPeak_oE.norm=-d[d$LFC_k27acPeak_oE<0,]$LFC_k27acPeak_oE/min
# d[d$LFC_k27acPeak_oE>0,]$LFC_k27acPeak_oE.norm=d[d$LFC_k27acPeak_oE>0,]$LFC_k27acPeak_oE/max
# 
# min=min(d$LFC.sum_bednorm_reads)
# max=max(d$LFC.sum_bednorm_reads)
# d[d$LFC.sum_bednorm_reads<0,]$LFC.sum_bednorm_reads.norm=-d[d$LFC.sum_bednorm_reads<0,]$LFC.sum_bednorm_reads/min
# d[d$LFC.sum_bednorm_reads>0,]$LFC.sum_bednorm_reads.norm=d[d$LFC.sum_bednorm_reads>0,]$LFC.sum_bednorm_reads/max
# 
# min=min(d$LFC.sum_bednorm_scores)
# max=max(d$LFC.sum_bednorm_scores)
# d[d$LFC.sum_bednorm_scores<0,]$LFC.sum_bednorm_scores.norm=-d[d$LFC.sum_bednorm_scores<0,]$LFC.sum_bednorm_scores/min
# d[d$LFC.sum_bednorm_scores>0,]$LFC.sum_bednorm_scores.norm=d[d$LFC.sum_bednorm_scores>0,]$LFC.sum_bednorm_scores/max
# 
# 
# baited_genes$lfc.k27.norm=0
# baited_genes$lfc.reads.norm=0
# baited_genes$lfc.scores.norm=0
# 
# for (n in baited_genes$genenames){
#   print(n)
# 
#   dd=d[d$genenames==n,]
#   baited_genes[baited_genes$genenames==n]$lfc.k27.norm=sum(dd$LFC_k27acPeak_oE.norm)
#   baited_genes[baited_genes$genenames==n]$lfc.reads.norm=sum(dd$LFC.sum_bednorm_reads.norm)
#   baited_genes[baited_genes$genenames==n]$lfc.scores.norm=sum(dd$LFC.sum_bednorm_scores.norm)
# }
# 
# ##SD as a measure of variability in the inputs (i.e.: do stable genes switch their input signals?):
# ##Conclusion: Yes, they seem to have more variability in how these genes engage enhancers!
# 
# baited_genes$sd.k27.norm=0
# baited_genes$sd.reads.norm=0
# baited_genes$sd.scores.norm=0
# 
# for (n in baited_genes$genenames){
#   print(n)
# 
#   dd=d[d$genenames==n,]
#   baited_genes[baited_genes$genenames==n]$sd.k27.norm=sd(dd$LFC_k27acPeak_oE.norm)
#   baited_genes[baited_genes$genenames==n]$sd.reads.norm=sd(dd$LFC.sum_bednorm_reads.norm)
#   baited_genes[baited_genes$genenames==n]$sd.scores.norm=sd(dd$LFC.sum_bednorm_scores.norm)
# }


genetypes=unique(baited_genes$genetype)
cols=c('cornflowerblue','orange','green','red')

par(mfcol=c(2,4))

for (i in 1:length(genetypes)){

  gt=genetypes[i]
  col=cols[i]
  
  plot(baited_genes$lfc.k27.norm,baited_genes$lfc.scores.norm,main=gt)
  points(baited_genes[baited_genes$genetype==gt]$lfc.k27.norm,baited_genes[baited_genes$genetype==gt]$lfc.scores.norm,col=col,pch=20)
  abline(h=0,lty=2,col='grey')
  abline(v=0,lty=2,col='grey')
  
  plot(baited_genes$lfc.k27.norm,baited_genes$lfc.reads.norm,main=gt)
  points(baited_genes[baited_genes$genetype==gt]$lfc.k27.norm,baited_genes[baited_genes$genetype==gt]$lfc.reads.norm,col=col,pch=20)
  abline(h=0,lty=2,col='grey')
  abline(v=0,lty=2,col='grey')
}


genetypes=unique(baited_genes$genetype)
cols=c('cornflowerblue','orange','green','red')

par(mfcol=c(2,4))

for (i in 1:length(genetypes)){

  gt=genetypes[i]
  col=cols[i]
  
  plot(baited_genes$sd.k27.norm,baited_genes$sd.scores.norm,main=gt)
  points(baited_genes[baited_genes$genetype==gt]$sd.k27.norm,baited_genes[baited_genes$genetype==gt]$sd.scores.norm,col=col,pch=20)
  abline(h=0,lty=2,col='grey')
  abline(v=0,lty=2,col='grey')
  
  plot(baited_genes$sd.k27.norm,baited_genes$sd.reads.norm,main=gt)
  points(baited_genes[baited_genes$genetype==gt]$sd.k27.norm,baited_genes[baited_genes$genetype==gt]$sd.reads.norm,col=col,pch=20)
  abline(h=0,lty=2,col='grey')
  abline(v=0,lty=2,col='grey')
}

par(mfcol=c(2,3))
boxplot(sd.k27.norm~genetype,data=as.data.frame(baited_genes),col=cols,main='sd.k27.norm',las=3,notch=TRUE)
boxplot(lfc.k27.norm~genetype,data=as.data.frame(baited_genes),col=cols,main='lfc.k27.norm',outline=FALSE,notch=TRUE)
abline(h=0,lty=2)
boxplot(sd.reads.norm~genetype,data=as.data.frame(baited_genes),col=cols,main='sd.reads.norm(ints)',las=3,notch=TRUE)
boxplot(lfc.reads.norm~genetype,data=as.data.frame(baited_genes),col=cols,main='lfc.reads.norm(ints)',outline=FALSE,notch=TRUE)
abline(h=0,lty=2)
boxplot(sd.scores.norm~genetype,data=as.data.frame(baited_genes),col=cols,main='sd.scores.norm(ints)',las=3,notch=TRUE)
boxplot(lfc.scores.norm~genetype,data=as.data.frame(baited_genes),col=cols,main='lfc.scores.norm(ints)',outline=FALSE,notch=TRUE)
abline(h=0,lty=2)


##combine loss of interaction with the loss of K27ac signal 

baited_genes$lfc.k27.intReads.norm=baited_genes$lfc.k27.norm+baited_genes$lfc.reads.norm
baited_genes$lfc.k27.intScores.norm=baited_genes$lfc.k27.norm+baited_genes$lfc.scores.norm

genetypes=unique(baited_genes$genetype)
cols=c('cornflowerblue','orange','green','red')

par(mfcol=c(1,4))

for (i in 1:length(genetypes)){

  gt=genetypes[i]
  col=cols[i]
  
  plot(baited_genes$lfc.k27.intReads.norm,baited_genes$lfc.k27.intScores.norm,main=gt)
  points(baited_genes[baited_genes$genetype==gt]$lfc.k27.intReads.norm,baited_genes[baited_genes$genetype==gt]$lfc.k27.intScores.norm,col=col,pch=20)
  abline(h=0,lty=2,col='grey')
  abline(v=0,lty=2,col='grey')
  
}


dev.off()


# ##only once:
# # write.table(d,'test.txt',sep='\t',eol='\n',quote=FALSE,row.names=FALSE)
# # head(read.delim('test.txt'))
# write.table(d,paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_interacting_readsScoresInPeak.txt'),sep='\t',eol='\n',quote=FALSE,row.names=FALSE)
# 
# saveRDS(baited_genes, '../designDir/baited_genes_annotated.rds')


##Simple Pdfs:

pdf('Regulation_of_K27ac_and_interactions_perBait_simple.pdf',onefile=TRUE, width=10,height=7)
par(mfrow=c(1,3))

main='Sum K27ac LFC at connected peaks'
boxplot(baited_genes$lfc.k27~baited_genes$genetype_simple,outline=FALSE,notch=TRUE, col=c('cornflowerblue','grey','orange'),ylab='cumulative LFC interacting K27ac (Cast/Bl6)',xlab='genetype',main=main)
abline(h=0,lty=2)

main='Sum interaction reads LFC at connected peaks'
boxplot(baited_genes$lfc.reads~baited_genes$genetype_simple,outline=FALSE,notch=TRUE, col=c('cornflowerblue','grey','orange'),ylab='cumulative LFC CapC reads (Cast/Bl6)',xlab='genetype',main=main)
abline(h=0,lty=2)

main='Sum interaction scores LFC at connected peaks'
boxplot(baited_genes$lfc.scores~baited_genes$genetype_simple,outline=FALSE,notch=TRUE, col=c('cornflowerblue','grey','orange'),ylab='cumulative LFC CapC scores (Cast/Bl6)',xlab='genetype',main=main)
abline(h=0,lty=2)

dev.off()


pdf('Regulation_of_K27ac_and_interactions_perBait_mean_simple.pdf',onefile=TRUE, width=10,height=7)
par(mfrow=c(1,3))

main='Mean K27ac LFC at connected peaks'
boxplot(baited_genes$lfc.k27_mean~baited_genes$genetype_simple,outline=FALSE,notch=TRUE, col=c('cornflowerblue','grey','orange'),ylab='mean LFC interacting K27ac per prom (Cast/Bl6)',xlab='genetype',main=main)
abline(h=0,lty=2)

main='Mean interaction reads LFC at connected peaks'
boxplot(baited_genes$lfc.reads_mean~baited_genes$genetype_simple,outline=FALSE,notch=TRUE, col=c('cornflowerblue','grey','orange'),ylab='mean LFC CapC reads per prom (Cast/Bl6)',xlab='genetype',main=main)
abline(h=0,lty=2)

main='Mean interaction scores LFC at connected peaks'
boxplot(baited_genes$lfc.scores_mean~baited_genes$genetype_simple,outline=FALSE,notch=TRUE, col=c('cornflowerblue','grey','orange'),ylab='mean LFC CapC scores per prom (Cast/Bl6)',xlab='genetype',main=main)
abline(h=0,lty=2)

dev.off()


pdf('Regulation_of_K27ac_surrounding_promoters_mean_simple.pdf',onefile=TRUE, width=13,height=7)
par(mfrow=c(1,4))

ylim=c(-4,4)

main='Mean K27ac LFC at peaks within 100kb of TSS'
boxplot(baited_genes$lfc.k27_mean_peaks_within_100kb~baited_genes$genetype_simple,outline=FALSE,notch=TRUE, col=c('cornflowerblue','grey','orange'),ylab='mean LFC interacting K27ac per prom (Cast/Bl6)',xlab='genetype',main=main,ylim=ylim)
abline(h=0,lty=2)

main='Mean K27ac LFC at peaks within 250kb of TSS'
boxplot(baited_genes$lfc.k27_mean_peaks_within_250kb~baited_genes$genetype_simple,outline=FALSE,notch=TRUE, col=c('cornflowerblue','grey','orange'),ylab='mean LFC interacting K27ac per prom (Cast/Bl6)',xlab='genetype',main=main,ylim=ylim)
abline(h=0,lty=2)

main='Mean K27ac LFC at peaks within 500kb of TSS'
boxplot(baited_genes$lfc.k27_mean_peaks_within_500kb~baited_genes$genetype_simple,outline=FALSE,notch=TRUE, col=c('cornflowerblue','grey','orange'),ylab='mean LFC interacting K27ac per prom (Cast/Bl6)',xlab='genetype',main=main,ylim=ylim)
abline(h=0,lty=2)

main='Mean K27ac LFC at peaks within 1Mb of TSS'
boxplot(baited_genes$lfc.k27_mean_peaks_within_1Mb~baited_genes$genetype_simple,outline=FALSE,notch=TRUE, col=c('cornflowerblue','grey','orange'),ylab='mean LFC interacting K27ac per prom (Cast/Bl6)',xlab='genetype',main=main,ylim=ylim)
abline(h=0,lty=2)

dev.off()

##Conclusion: Stable genes seem to have more variability in how these genes engage enhancers!




```


Since stable genes show more variable enhancer engagement (i.e. interactions):
Test whether each gene that loses an interaction also gains another interaction

```{r}


i=2
d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_interacting_readsScoresInPeak.txt'),stringsAsFactors=FALSE)

baited_genes$max.lfc.k27='NA'
baited_genes$min.lfc.k27='NA'
baited_genes$max.lfc.reads='NA'
baited_genes$min.lfc.reads='NA'
baited_genes$max.lfc.scores='NA'
baited_genes$min.lfc.scores='NA'

for (n in baited_genes$genenames){
  print(n)

  dd=d[d$genenames==n,]
  
  ##k27ac:
  ma=dd[dd$LFC_k27acPeak_oE==max(dd$LFC_k27acPeak_oE),]
  max=paste0(ma$seqnames_otherEnd,':',ma$start_otherEnd,'-',ma$end_otherEnd,';',ma$LFC_k27acPeak_oE)
  if (nrow(ma)>1){max=concatenate(max)}
  ma=dd[dd$LFC_k27acPeak_oE==min(dd$LFC_k27acPeak_oE),]
  min=paste0(ma$seqnames_otherEnd,':',ma$start_otherEnd,'-',ma$end_otherEnd,';',ma$LFC_k27acPeak_oE)
  if (nrow(ma)>1){min=concatenate(min)}
  baited_genes[baited_genes$genenames==n]$max.lfc.k27=max
  baited_genes[baited_genes$genenames==n]$min.lfc.k27=min
  
  ##reads:
  ma=dd[dd$LFC.sum_bednorm_reads==max(dd$LFC.sum_bednorm_reads),]
  max=paste0(ma$seqnames_otherEnd,':',ma$start_otherEnd,'-',ma$end_otherEnd,';',ma$LFC.sum_bednorm_reads)
  if (nrow(ma)>1){max=concatenate(max)}
  ma=dd[dd$LFC.sum_bednorm_reads==min(dd$LFC.sum_bednorm_reads),]
  min=paste0(ma$seqnames_otherEnd,':',ma$start_otherEnd,'-',ma$end_otherEnd,';',ma$LFC.sum_bednorm_reads)
  if (nrow(ma)>1){min=concatenate(min)}
  baited_genes[baited_genes$genenames==n]$max.lfc.reads=max
  baited_genes[baited_genes$genenames==n]$min.lfc.reads=min
  
  ##scores:
  ma=dd[dd$LFC.sum_bednorm_scores==max(dd$LFC.sum_bednorm_scores),]
  max=paste0(ma$seqnames_otherEnd,':',ma$start_otherEnd,'-',ma$end_otherEnd,';',ma$LFC.sum_bednorm_scores)
  if (nrow(ma)>1){max=concatenate(max)}
  ma=dd[dd$LFC.sum_bednorm_scores==min(dd$LFC.sum_bednorm_scores),]
  min=paste0(ma$seqnames_otherEnd,':',ma$start_otherEnd,'-',ma$end_otherEnd,';',ma$LFC.sum_bednorm_scores)
  if (nrow(ma)>1){min=concatenate(min)}
  baited_genes[baited_genes$genenames==n]$max.lfc.scores=max
  baited_genes[baited_genes$genenames==n]$min.lfc.scores=min
 
}


##function to test whether max and min are opposite (down vs up):
##test it by multiplying teh two extreme numbers


##1) K27ac

par(mfrow=c(2,3))
m=data.frame(baited_genes$max.lfc.k27, baited_genes$min.lfc.k27 )

##For Scatterplots: Convert m values to numeric values:

f <- function(x){
  
  x=as.character(x)
  
  if(grepl(',',x)){
    a=unlist(strsplit(x,split=','))[1]
  } else{
    a=x
  }
  a=as.numeric(unlist(strsplit(a,split=';'))[2])
  
  return(a)
}

m[,1]=as.character(m[,1])
m[,2]=as.character(m[,2])
m[,1]=sapply(m[,1],FUN=f)
m[,2]=sapply(m[,2],FUN=f)

genetypes=unique(baited_genes$genetype)

##For Barplots: check if the values are the opposite of each other:
baited_genes$lfc.k27.opposite=apply(m,1,FUN=function(x) x[1]*x[2]<0 )
baited_genes$lfc.k27.maxmin.bothover1.3=apply(m,1,FUN=function(x) abs(x[1])>log2(1.3) & abs(x[2]) >log2(1.3) )
baited_genes$lfc.k27.opposite.bothover1.3=FALSE
baited_genes[baited_genes$lfc.k27.opposite & baited_genes$lfc.k27.maxmin.bothover1.3 & is.na(baited_genes$lfc.k27.opposite)==FALSE]$lfc.k27.opposite.bothover1.3=TRUE
tab=table(baited_genes$lfc.k27.opposite.bothover1.3,baited_genes$genetype)
tab=tab[c(2,1),]
plotNormBarplot(tab,main='LFC k27ac oE bidirectional?',las=3)


##2) intsReads:

m=data.frame(baited_genes$max.lfc.reads, baited_genes$min.lfc.reads )

m[,1]=as.character(m[,1])
m[,2]=as.character(m[,2])
m[,1]=sapply(m[,1],FUN=f)
m[,2]=sapply(m[,2],FUN=f)

##For Barplots: check if the values are the opposite of each other:
baited_genes$lfc.intsReads.opposite=apply(m,1,FUN=function(x) x[1]*x[2]<0 )

baited_genes$lfc.intsReads.maxmin.bothover1.3=apply(m,1,FUN=function(x) abs(x[1])>log2(1.3) & abs(x[2]) >log2(1.3) )
baited_genes$lfc.intsReads.opposite.bothover1.3=FALSE
baited_genes[baited_genes$lfc.intsReads.opposite & baited_genes$lfc.intsReads.maxmin.bothover1.3 & is.na(baited_genes$lfc.intsReads.opposite)==FALSE]$lfc.intsReads.opposite.bothover1.3=TRUE

tab=table(baited_genes$lfc.intsReads.opposite.bothover1.3,baited_genes$genetype)
tab=tab[c(2,1),]
plotNormBarplot(tab,main='LFC intsReads oE bidirectional?',las=3)



##3) intsScores:

m=data.frame(baited_genes$max.lfc.scores, baited_genes$min.lfc.scores )

m[,1]=as.character(m[,1])
m[,2]=as.character(m[,2])
m[,1]=sapply(m[,1],FUN=f)
m[,2]=sapply(m[,2],FUN=f)

##For Barplots: check if the values are the opposite of each other:
baited_genes$lfc.intsScores.opposite=apply(m,1,FUN=function(x) x[1]*x[2]<0 )

baited_genes$lfc.intsScores.maxmin.bothover1.3=apply(m,1,FUN=function(x) abs(x[1])>log2(1.3) & abs(x[2]) >log2(1.3) )
baited_genes$lfc.intsScores.opposite.bothover1.3=FALSE
baited_genes[baited_genes$lfc.intsScores.opposite & baited_genes$lfc.intsScores.maxmin.bothover1.3 & is.na(baited_genes$lfc.intsScores.opposite)==FALSE]$lfc.intsScores.opposite.bothover1.3=TRUE

tab=table(baited_genes$lfc.intsScores.opposite.bothover1.3,baited_genes$genetype)
tab=tab[c(2,1),]
plotNormBarplot(tab,main='LFC intsScores oE bidirectional?',las=3)



pdf('Barplots_bidirectional_input.pdf')
par(mfrow=c(2,3))

##for empirical pvalues:
testcase=data.frame(values(baited_genes))
reg='genetype'
ntest=10000

char='lfc.k27.opposite'
tab=table(as.data.frame(baited_genes)[,char],baited_genes$genetype)
tab=tab[c(2,1),]
plotNormBarplot(tab,main='LFC K27ac bidirectional?',las=3,include.pvals=TRUE,df=testcase,reg=reg,char=char,ntest=ntest)

char='lfc.intsReads.opposite'
tab=table(as.data.frame(baited_genes)[,char],baited_genes$genetype)
tab=tab[c(2,1),]
plotNormBarplot(tab,main='LFC intsReads bidirectional?',las=3,include.pvals=TRUE,df=testcase,reg=reg,char=char,ntest=ntest)

char='lfc.intsScores.opposite'
tab=table(as.data.frame(baited_genes)[,char],baited_genes$genetype)
tab=tab[c(2,1),]
plotNormBarplot(tab,main='LFC intsScores bidirectional?',las=3,include.pvals=TRUE,df=testcase,reg=reg,char=char,ntest=ntest)

##bidirectional AND both directions >1.3 fold change

char='lfc.k27.opposite.bothover1.3'
tab=table(as.data.frame(baited_genes)[,char],baited_genes$genetype)
tab=tab[c(2,1),]
plotNormBarplot(tab,main='LFC K27ac >1.3 bidirectional?',las=3,include.pvals=TRUE,df=testcase,reg=reg,char=char,ntest=ntest)

char='lfc.intsReads.opposite.bothover1.3'
tab=table(as.data.frame(baited_genes)[,char],baited_genes$genetype)
tab=tab[c(2,1),]
plotNormBarplot(tab,main='LFC intsReads >1.3 bidirectional?',las=3,include.pvals=TRUE,df=testcase,reg=reg,char=char,ntest=ntest)

char='lfc.intsScores.opposite.bothover1.3'
tab=table(as.data.frame(baited_genes)[,char],baited_genes$genetype)
tab=tab[c(2,1),]
plotNormBarplot(tab,main='LFC intsScores >1.3 bidirectional?',las=3,include.pvals=TRUE,df=testcase,reg=reg,char=char,ntest=ntest)

dev.off()

```


##Hypothesis:




##Annotate baited_genes with n k27ac interactions, n k27ac[Bl6] interactions and n k27ac[Cast] interactions
##to see if there is a bias for genes that are Bl6 or Castaneus-specific

```{r}

i=2
d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_newk27acPeaks_interacting_readsScoresInPeak.txt'),stringsAsFactors=FALSE)

baited_genes$n_k27acPeaks_ints_cutoff3scores=0
baited_genes$n_k27ac_loss_cutoff3scores=0
baited_genes$n_k27ac_gain_cutoff3scores=0

cutoff=1.5
lost=d[d$LFC_k27acPeak_oE <= -log2(cutoff) & d$padj_k27acPeak_oE < 0.05,]
gained=d[d$LFC_k27acPeak_oE >= log2(cutoff) & d$padj_k27acPeak_oE < 0.05,]

for (j in 1:length(baited_genes)){
  print(j)
  
  dd=d[d$baitID==baited_genes[j]$dpn_id,]
  if ( nrow(dd)>0){
    baited_genes[j]$n_k27acPeaks_ints_cutoff3scores = nrow(dd)
  }
  
  dd=lost[lost$baitID==baited_genes[j]$dpn_id,]
  if ( nrow(dd)>0){
    baited_genes[j]$n_k27ac_loss_cutoff3scores = nrow(dd)
  }
  
  dd=gained[gained$baitID==baited_genes[j]$dpn_id,]
  if ( nrow(dd)>0){
    baited_genes[j]$n_k27ac_gain_cutoff3scores = nrow(dd)
  }
  
}

# ##just once:
# saveRDS(baited_genes,'baited_genes_annotated.rds')

plotwhat='n_k27ac_loss_cutoff3scores'
plotcumulative <- function(plotwhat){

  genetypes=unique(baited_genes$genetype)
  dd=data.frame(n_loss=0:max(as.data.frame(baited_genes)[,plotwhat]))
  
  for (g in genetypes){
    v=vector()
    for (j in dd[,1]){
      v=c(v,length(baited_genes[baited_genes$genetype==g & as.data.frame(baited_genes)[,plotwhat]>=j]))
    }
    dd[,g]=v
  }
  
  dd.norm=dd
  
  for (j in 2:ncol(dd)){
    dd.norm[,j]=dd[,j]/dd[1,j]
  }
  
  xlab=plotwhat
  ylab=paste('fraction genes with', plotwhat)
  cols=c('cornflowerblue','darkorange','green','red')
  plot(dd.norm[,1],dd.norm[,2],type='l',col=cols[1], main=plotwhat,xlab=xlab,ylab=ylab)
  for (j in 3:ncol(dd)){
    lines(dd[,1],dd.norm[,j],col=cols[j-1])
  }
  legend('topright',legend=genetypes,col=cols,lwd=2,bty='n')
  
} 

pdf('cumulative_lineplots_n_k27ac_differential_ints.pdf')
par(mfrow=c(1,2))
plotcumulative('n_k27ac_loss_cutoff3scores')
plotcumulative('n_k27ac_gain_cutoff3scores')
dev.off()

```



###########################################
##14.01.2022
###########################################

##A number of genes looks phenotypically like Scc1-degron
##Are cohesin genes misregulated? 

##Many changes occur at TAD-borders
##Is CTCF misregulated/motif mutated?

```{r}

##load gene expression data:
d=read.csv(paste0(base,'analysis/MouseGenomes/rnaseq/masked_SNP.SV/ResultTables/refSeq.geneBodyUnion.ignore.strand.TRUE.inter.feature.FALSE.Peaks.Castaneus_vs_Bl6_annotated_q0.005_fe1.5.csv'),stringsAsFactors=FALSE)

cohesin_genes=c('Rad21','Scc1','Scc3','Smc1','Smc3', 'Stag3', 'Wapl','Ctcf','Smc1a','Smc1b','Stag1','Stag2','Rec8')
cg=d[d$genenames %in% cohesin_genes,]

p=0.05
d$log2.Bl6.pooled_fpkm=log2(d$Bl6.pooled_fpkm+p)

pdf('Cohesin_domains_expression.pdf',width=10,height=7)

plot(d$log2.Bl6.pooled_fpkm,d$LFC,col=rgb(0,0,0,alpha=0.1))
abline(h=0,lty=2,lwd=2,col='darkgreen')
abline(h=log2(1.5),lty=2,lwd=2,col='green')
abline(h=-log2(1.5),lty=2,lwd=2,col='green')

cohesin_genes=c('Rad21','Scc1','Scc3','Smc1','Smc3', 'Stag3', 'Wapl','Ctcf','Smc1a','Smc1b','Stag1','Stag2','Rec8')
cg=d[d$genenames %in% cohesin_genes,]

points(cg$log2.Bl6.pooled_fpkm,cg$LFC,pch=20,col=rgb(1,0,0,alpha=0.8))
text(cg$log2.Bl6.pooled_fpkm,cg$LFC+0.3,labels=cg$genenames,col='red')

dev.off()


##genes that look like they have a cohesin phenotype:
cohesin_genes=c('Mtmr7','Adgre5','Mmp14', 'Mras', 'Capns','Ncoa7','Sumf1','Dusp4','Fbrsl1','Syne4','Agtrap','Trpc1','Rps24','Vapa','Stxbp2','Syt5','Ccnl1','Tmem17','Gars','Fam187b')

cg=d[d$genenames %in% cohesin_genes,]

p=0.05
d$log2.Bl6.pooled_fpkm=log2(d$Bl6.pooled_fpkm+p)

pdf('Cohesin_phenotype_gene_expression.pdf',width=10,height=7)

plot(d$log2.Bl6.pooled_fpkm,d$LFC,col=rgb(0,0,0,alpha=0.1))
abline(h=0,lty=2,lwd=2,col='darkgreen')
abline(h=log2(1.5),lty=2,lwd=2,col='green')
abline(h=-log2(1.5),lty=2,lwd=2,col='green')

cohesin_genes=c('Mtmr7','Adgre5','Mmp14', 'Mras', 'Capns','Ncoa7','Sumf1','Dusp4','Fbrsl1','Syne4','Agtrap','Trpc1','Rps24','Vapa','Stxbp2','Syt5','Ccnl1','Tmem17','Gars','Fam187b')
cg=d[d$genenames %in% cohesin_genes,]

points(cg$log2.Bl6.pooled_fpkm,cg$LFC,pch=20,col=rgb(1,0,0,alpha=0.8))
text(cg$log2.Bl6.pooled_fpkm,cg$LFC+0.3,labels=cg$genenames,col='red')

dev.off()


##Are there more SVs/SNPs at TSSs of genes that are affected cohesin-like?

d$cohesin_type='not_known'
d[d$genenames %in% baited_genes$genenames,]$cohesin_type='not_affected'
d[d$genenames %in% cg$genenames,]$cohesin_type='affected'

for (what in names(d)[29:35]){
tab=table(d[,what],d$cohesin_type)
  plotNormBarplot_sel(tab,sel='TRUE',main=paste('Fraction genes with',what),xlab='genes captureC cohesin phenotype')
}

baited_genes$cohesin_phenotype='not_affected'
baited_genes[baited_genes$genenames %in% cohesin_genes]$cohesin_phenotype='affected'
tab=table(baited_genes$cohesin_phenotype,baited_genes$genetype)
plotNormBarplot_sel(tab,sel='affected',main=paste('Fraction genes with cohesin phenotype'),xlab='genetype')

```



#################################################################################################

#####################################################
##Plot replicate line plots of random selected baits:
#####################################################

```{r}

ids=unique(L[[1]]$baitID)
sam=sample(ids,20)

pdf(paste('Lineplots_sample_reps.pdf',sep=''), width=20, height=14,onefile = TRUE,pointsize=12)

for (id in sam){
  # for (id in ids[92:length(ids)]){
  
  par(mfcol=c(2,3))
  print(id)
  name=paste(baited_genes[baited_genes$dpn_id==id,]$genenames,baited_genes[baited_genes$dpn_id==id,]$inttype,baited_genes[baited_genes$dpn_id==id,]$genetype,sep=';')
  
  colintervals=c('lightgrey','purple','cornflowerblue','orange', 'red', 'blue', 'green', 'pink', 'black')
  
  col=brewer.pal(9,'Blues')[c(3,5,7,9)]
    
  k=10
  zoom=200000
  ylim=c(0,0.1)
  plotInteractions(LL,id,k,zoom,unt=1,tam=2,unt2=3, tam2=4, ylim=ylim,show.legend = TRUE,name=name,intervals=grl,colintervals = colintervals, col=col)
  
  ylim=c(0,0.025)
  zoom=2000000
  k=40
  plotInteractions(LL,id,k,zoom,unt=1,tam=2,unt2=3, tam2=4, ylim=ylim,show.legend = TRUE,name=name,intervals=grl,colintervals = colintervals, col=col)
  
  col=brewer.pal(9,'Reds')[c(3,5,7,9)]
    
  k=10
  zoom=200000
  ylim=c(0,0.1)
  plotInteractions(LL,id,k,zoom,unt=5,tam=6,unt2=7, tam2=8, ylim=ylim,show.legend = TRUE,name=name,intervals=grl,colintervals = colintervals, col=col)
  
  ylim=c(0,0.025)
  zoom=2000000
  k=40
  plotInteractions(LL,id,k,zoom,unt=5,tam=6,unt2=7, tam2=8, ylim=ylim,show.legend = TRUE,name=name,intervals=grl,colintervals = colintervals, col=col)
  
  col=c('cornflowerblue','orange')
  
   k=10
  zoom=200000
  ylim=c(0,0.1)
  plotInteractions(L,id,k,zoom,unt=1,tam=2, ylim=ylim,show.legend = TRUE,name=name,intervals=grl,colintervals = colintervals, col=col)
  
  ylim=c(0,0.025)
  zoom=2000000
  k=40
  plotInteractions(L,id,k,zoom,unt=1,tam=2, ylim=ylim,show.legend = TRUE,name=name,intervals=grl,colintervals = colintervals, col=col)
  
}
dev.off()



ids=unique(L[[1]]$baitID)
sam=sample(ids,20)

pdf(paste('Lineplots_sample_reps_qualitative.pdf',sep=''), width=20, height=14,onefile = TRUE,pointsize=12)

for (id in sam){
  # for (id in ids[92:length(ids)]){
  
  par(mfrow=c(2,3))
  print(id)
  name=paste(baited_genes[baited_genes$dpn_id==id,]$genenames,baited_genes[baited_genes$dpn_id==id,]$inttype,baited_genes[baited_genes$dpn_id==id,]$genetype,sep=';')
  
  colintervals=c('lightgrey','purple','cornflowerblue','orange', 'red', 'blue', 'green', 'pink', 'black')
  col=c('cornflowerblue','orange')
  
  ylim=c(0,0.025)
  zoom=2000000
  k=40
  plotInteractions(L,id,k,zoom,unt=1,tam=2, ylim=ylim,show.legend = TRUE,name=name,intervals=grl,colintervals = colintervals, col=col) 
  
  ylim=c(0,0.05)
  for (i in 1:4){
    plotInteractions(LL,id,k,zoom,unt=i,tam=(i+4), ylim=ylim,show.legend = TRUE,name=name,intervals=grl,colintervals = colintervals, col=col)
  }
}
dev.off()

```



##Annotate ints with genetypes from the baited_genes file:

```{r}

i=GRanges(1,IRanges(ints$baitID,ints$baitID))
b=GRanges(1,IRanges(baited_genes$dpn_id,baited_genes$dpn_id))

ov=as.matrix(findOverlaps(i,b))
ints=ints[ov[,1],]
ints=cbind(ints,values(baited_genes[ov[,2]])[,c(26:33,42,43,46,48,50:63,65,67,69:72,77,79,89)])

##just once:
write.table(ints,'ints_all.txt',sep='\t',eol='\n',quote=FALSE)

```


##Annotate ints with k27ac data:

```{r}

##1) otherEnd

i=bed2gr(ints[,c("seqnames_otherEnd" ,        "start_otherEnd",            "end_otherEnd" )])
b=bed2gr(k27ac_all[,c("Chr"                  ,"Start"  ,               "End")])

ov=as.matrix(findOverlaps(i,b))
intsnotov=ints[-ov[,1],] ##interactions without K27ac peak at otherEnd
ints=ints[ov[,1],] ##interactions with K27ac peak at other End

k=k27ac_all[ov[,2],8:ncol(k27ac_all)]
names(k)=paste(names(k),'k27acPeak_oE',sep='_')
ints=cbind(ints,k)

##add empty values to all interactions without K27ac peak at otherEnd:
kk=as.data.frame(matrix(ncol=ncol(k),nrow=nrow(intsnotov)))
names(kk)=names(k)
intsnotov=cbind(intsnotov,kk)

ints=rbind(ints,intsnotov)


##2) bait

i=bed2gr(ints[,c("seqnames_bait" ,        "start_bait",            "end_bait" )])
b=bed2gr(k27ac_all[,c("Chr"                  ,"Start"  ,               "End")])

ov=as.matrix(findOverlaps(i,b))
intsnotov=ints[-ov[,1],]
ints=ints[ov[,1],]

k=k27ac_all[ov[,2],8:ncol(k27ac_all)]
names(k)=paste(names(k),'k27acPeak_bait',sep='_')
ints=cbind(ints,k)

##add empty values to all interactions without K27ac peak at otherEnd:
kk=as.data.frame(matrix(ncol=ncol(k),nrow=nrow(intsnotov)))
names(kk)=names(k)
intsnotov=cbind(intsnotov,kk)

ints=rbind(ints,intsnotov)


##just once:
write.table(ints,'ints_all.txt',sep='\t',eol='\n',quote=FALSE)

```



##Define clusters based on interactions:
##(not done!)

```{r}

ints$clustertype='stable'

ints[ints$Bl6_score > ints$Castaneus_score & ints$Bl6_N_downsampled>ints$Castaneus_N_downsampled & 
       ints$Bl6_N.1_downsampled > ints$Castaneus_N.1_downsampled & 
       ints$Bl6_N.2_downsampled > ints$Castaneus_N.2_downsampled &
       ints$Bl6_N.3_downsampled > ints$Castaneus_N.3_downsampled &
       ints$Bl6_N.4_downsampled > ints$Castaneus_N.4_downsampled &
       ints$Bl6_rep3_score > ints$Castaneus_rep2_score &
       ints$Bl6_rep4_score > ints$Castaneus_rep3_score &
       ints$Bl6_rep5_score > ints$Castaneus_rep4_score &
       ints$Bl6_rep6_score > ints$Castaneus_rep5_score ,]$clustertype='Bl6_up'

ints[ints$Bl6_score < ints$Castaneus_score & ints$Bl6_N_downsampled<ints$Castaneus_N_downsampled & 
       ints$Bl6_N.1_downsampled < ints$Castaneus_N.1_downsampled & 
       ints$Bl6_N.2_downsampled < ints$Castaneus_N.2_downsampled &
       ints$Bl6_N.3_downsampled < ints$Castaneus_N.3_downsampled &
       ints$Bl6_N.4_downsampled < ints$Castaneus_N.4_downsampled &
       ints$Bl6_rep3_score < ints$Castaneus_rep2_score &
       ints$Bl6_rep4_score < ints$Castaneus_rep3_score &
       ints$Bl6_rep5_score < ints$Castaneus_rep4_score &
       ints$Bl6_rep6_score < ints$Castaneus_rep5_score ,]$clustertype='Castaneus_up'

ints[ints$Bl6_score >=5  & ints$Castaneus_score < 4 & 
        log2((ints$Bl6_N_downsampled+1)/(ints$Castaneus_N_downsampled+1)) >= log2(1.5) &
       ints$clustertype == 'Bl6_up',]$clustertype='Bl6_specific'

ints[ints$Bl6_score < 4  & ints$Castaneus_score >= 5 & 
        log2((ints$Bl6_N_downsampled+1)/(ints$Castaneus_N_downsampled+1)) <= -log2(1.5) &
       ints$clustertype == 'Castaneus_up',]$clustertype='Castaneus_specific'


ints$clustertype=factor(ints$clustertype, levels=c('stable','Bl6_specific', 'Bl6_up','Castaneus_specific',       'Castaneus_up'))


##just once:
write.table(ints,'ints_all.txt',sep='\t',eol='\n',quote=FALSE)

```


##Enrichments of different clustertypes within genetypes

```{r}

# defaultparameters:
#   par(mar=c(5, 4, 4, 2) + 0.1)

set.seed(100)

##empirical pvalues:
reg='genetype'
char='clustertype'
empir(ints,reg,char,1000)

ints$clustertype=factor(ints$clustertype, levels=c('stable','Bl6_specific', 'Bl6_up','Castaneus_specific',       'Castaneus_up'))


pdf('Barplots_clustertype_CapC_genetype.pdf',onefile=TRUE,width=14)

par(mar=c(13,4,4,2)+0.1)

tab=table(ints[,reg],ints[,char])
tab=cbind(tab,data.frame(table(ints[,reg]))[,2])
colnames(tab)[ncol(tab)]='all'
tab=tab[rowSums(tab)>0,]
col=brewer.pal(nrow(tab),'Set1')

par(mfrow=c(1,1))
plotNormBarplot(tab,col=col,ylab='fraction ints with genetype',ylim=c(0,1.3),ncol.legend = 3,cex.legend = 0.9,horiz=FALSE,las=3,xlab='clustertype interactions',legend.title = 'genetype')

for (i in 1:nrow(tab)){
  test=paste('empirPvalues',rownames(tab)[i],'genes',char,'1000_tests.txt',sep='_')
  test=read.delim(test)
  text=apply(test[,c('pvalue_enriched', 'pvalue_depleted')],1,FUN=function(x) min(x)) ##select minimum pvalue
  plotNormBarplot_sel(tab,main=paste('enrichment',rownames(tab)[i], 'genes'),col=col[i],las=3,sel=rownames(tab)[i],ylab='fraction ints with clustertype',pval=text)

}


tab=t(tab[,-ncol(tab)])
tab=cbind(tab,data.frame(table(ints[,char]))[,2])
colnames(tab)[ncol(tab)]='all'
tab=tab[rowSums(tab)>0,]
col=c(brewer.pal(6,'Set1'))

par(mfrow=c(1,1))
plotNormBarplot(tab,col=col,ylab='fraction ints with clustertype CapC',ylim=c(0,1.4),ncol.legend = 3,cex.legend = 0.7,horiz=FALSE,las=3,xlab='clustertype genes',legend.title = 'genetype')

for (i in 1:nrow(tab)){
  plotNormBarplot_sel(tab,main=paste('enrichment ',rownames(tab)[i], 'ints'),col=col[i],las=3,sel=rownames(tab)[i],ylab='fraction ints with clustertype CapC',xlab='clustertype genes')
}

dev.off()

## At end of plotting, reset to default settings:
par(mar=c(5, 4, 4, 2) + 0.1)

```


################################################
#################################
##Regulation of interactions
##in relationship to regulation of 
##k27ac peaks
#################################
################################################


##Single interaction space:
##(not done yet)

```{r}


# par(mar=c(13,4,4,2)+0.1)

##A) regulation by peak types:

d=ints[ints$k27ac_2i_otherEnd,]
d$k27ac_reg_oE='stable'
d[d$k27ac_Bl6_otherEnd,]$k27ac_reg_oE='Bl6_up'
d[d$k27ac_Cast_otherEnd,]$k27ac_reg_oE='Castaneus_up'
dd=data.frame( score=c(d$Bl6_score,d$Castaneus_score), celline=rep(c('Bl6','Castaneus'),each=nrow(d)), peaktype=rep(d$peaktype_k27acPeak_oE,2), reads_log2=c( log2(d$Bl6_N_downsampled+1),log2(d$Castaneus_N_downsampled+1) ), k27ac_reg_oE=c(d$k27ac_reg_oE,d$k27ac_reg_oE) )

col=c('cornflowerblue','orange')
boxplot(  score~celline*peaktype,data=dd,outline=FALSE,las=3,col=col,notch=TRUE)
abline(h=5,lty=2)

p <- ggplot(dd,aes(x=peaktype, y=score, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,20))
p1 <- p +  theme_classic() + scale_fill_manual(values=col) + ggtitle('regulation of interactions by peaktypes')

p <- ggplot(dd,aes(x=peaktype, y=reads_log2, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,8))
p2 <- p +  theme_classic() + scale_fill_manual(values=col) + ggtitle('regulation of interactions by peaktypes')

p <- ggplot(dd,aes(x=k27ac_reg_oE, y=score, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,20))
p3 <- p +  theme_classic() + scale_fill_manual(values=col) + ggtitle('regulation of interactions by peaktypes')

p <- ggplot(dd,aes(x=k27ac_reg_oE, y=reads_log2, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,8))
p4 <- p +  theme_classic() + scale_fill_manual(values=col) + ggtitle('regulation of interactions by peaktypes')

multiplot(p1,p2,p3,p4,cols=2)


##B) regulation by peak types for specific gene types:

for (typ in unique(ints$genetype)){
  d=ints[ints$k27ac_2i_otherEnd & ints$genetype==typ,]
  d$k27ac_reg_oE='stable'
  d[d$k27ac_Bl6_otherEnd,]$k27ac_reg_oE='Bl6_up'
  d[d$k27ac_Cast_otherEnd,]$k27ac_reg_oE='Castaneus_up'
  dd=data.frame( score=c(d$Bl6_score,d$Castaneus_score), celline=rep(c('Bl6','Castaneus'),each=nrow(d)), peaktype=rep(d$peaktype_k27acPeak_oE,2), reads_log2=c( log2(d$Bl6_N_downsampled+1),log2(d$Castaneus_N_downsampled+1) ), k27ac_reg_oE=c(d$k27ac_reg_oE,d$k27ac_reg_oE) )
  
  col=c('cornflowerblue','orange')
  
  p <- ggplot(dd,aes(x=peaktype, y=score, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,20))
  p1 <- p +  theme_classic() + scale_fill_manual(values=col) + ggtitle(paste(typ,'regulation of interactions by peaktypes'))
  
  p <- ggplot(dd,aes(x=peaktype, y=reads_log2, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,8))
  p2 <- p +  theme_classic() + scale_fill_manual(values=col) + ggtitle(paste(typ,'regulation of interactions by peaktypes'))
  
  p <- ggplot(dd,aes(x=k27ac_reg_oE, y=score, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,20))
  p3 <- p +  theme_classic() + scale_fill_manual(values=col) + ggtitle(paste(typ,'regulation of interactions by peaktypes'))
  
  p <- ggplot(dd,aes(x=k27ac_reg_oE, y=reads_log2, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,8))
  p4 <- p +  theme_classic() + scale_fill_manual(values=col) + ggtitle(paste(typ,'regulation of interactions by peaktypes'))
  
  multiplot(p1,p2,p3,p4,cols=2)

}

```

###############################################
##How are typical K27ac interactions regulated?
###############################################

##1) Determine all K27ac peaks that interact with a gene

##2) Make interaction matrices for each OneGeneOnePeak interaction

##3) Calculate scores and reads within these peaks


```{r}

############################################
##1) get OneGeneOnePeak tables:
############################################

getIntervalInteractionPeaks(grl[2:length(grl)],hind,ints)

############################################
##2)create metaprofile matrices and plot them:
############################################

##Make non-normalized matrices

##Plot aggregate profiles norm to Ring1bfl_UNT_G2

memory.limit(size=32000)

zoom=100
resfolder='matrices/'
dir.create(resfolder)
normsam='Bl6'

for (i in 2:(length(grl)-1)){
  # for (i in 8){
 print(names(grl)[i]) 
  dis=names(grl)[i]

  pdf(paste('Metaprofiles',names(grl)[i],'normalizedTo',normsam,'.pdf',sep='_'),onefile=TRUE,width=10,height=10)
  par(mfrow=c(2,2))
  ylim=c(0,2)
  
    
    
    d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting.txt'),stringsAsFactors=FALSE)
    
    ##only take peaks that are on the same chromosome:
    d=d[d$seqnames_bait==d$seqnames_otherEnd,]
    dim(d)
    
    
    test=getMatrix(L,zoom,d,resfolder=resfolder,type='reads',norm=FALSE,name=names(grl)[i])
    
    ##normalize 
    testn=test
    norm=test[[normsam]][,(zoom+1)]
    pseudo=1
    for (j in 1:length(testn)){
      a=cbind(testn[[j]],norm)
      b=t(apply(a,1,FUN=function(x) (x+pseudo)/(x[length(x)]+pseudo)))
      testn[[j]]=b[,1:(ncol(b)-1)]
    }
    
    col=c(brewer.pal(9,'Set1'))
    plotAggregatePeaks(testn,col=rep(col,2),lty=rep(c(1,2),each=length(testn)/2),mainprefix=paste(names(grl)[i], ' '),ylim=ylim)
    
    
    
    test=getMatrix(L,zoom,d,resfolder=resfolder,type='scores',norm=FALSE,name=names(grl)[i])
    
    ##normalize 
    testn=test
    norm=test[[normsam]][,(zoom+1)]
    pseudo=1
    for (j in 1:length(testn)){
      a=cbind(testn[[j]],norm)
      b=t(apply(a,1,FUN=function(x) (x+pseudo)/(x[length(x)]+pseudo)))
      testn[[j]]=b[,1:(ncol(b)-1)]
    }
    
    col=c(brewer.pal(9,'Set1'))
    plotAggregatePeaks(testn,col=rep(col,2),lty=rep(c(1,2),each=length(testn)/2),mainprefix=paste('scores,', names(grl)[i], ' '),ylim=ylim)
    
    
    dev.off()

}


############################################
##3) Calculate: reads and scores in features
############################################

getDataInPeaks(grl[2:(length(grl)-1)],zoom)

##4) Plot boxplots:

pdf('Boxplots_readsScoresInFeatures.pdf',onefile=TRUE,width=7, height=10)

for (i in 2:(length(grl)-1)){
  d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),stringsAsFactors = FALSE)
  # names(d)[grep('reads.1',names(d))]=paste(strsplit_string(names(d)[grep('reads.1',names(d))],s2='reads.1'),'scores',sep='')
  # write.table(d, paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),sep='\t',eol='\n',quote=FALSE,row.names=FALSE)
  # 
  par(mfcol=c(3,2))
  for (type in c('reads','scores')){
    # for (type in c('reads')){
    for (cl in c('sum','mean','summit')){
      y=paste(cl,'bednorm',type,sep='_')
      df=d[,grep(y,names(d))]
      names(df)=unlist(lapply(strsplit(names(df),split='[.]'),FUN=function(x) x[1]))
      col=rep(c('cornflowerblue','orange'),2)
      boxplot(df,outline=FALSE,las=3,notch=TRUE,ylab=y,main=paste(names(grl)[i]),col=col,border=rep(c('black','darkgrey'),each=2))
      legend('topright',bty='n',legend=c(paste0('n=',nrow(df)),paste0('nprom=',length(unique(d$genenames))))  )
    }
  }
}
dev.off()


for (g in unique(baited_genes$genetype)){
  
  pdf(paste('Boxplots_readsScoresInFeatures_',g,'.genes.pdf',sep=''),onefile=TRUE,width=7, height=10)
  
  for (i in 2:(length(grl)-1)){
    d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),stringsAsFactors = FALSE)
    d=d[d$genetype==g,]
    # names(d)[grep('reads.1',names(d))]=paste(strsplit_string(names(d)[grep('reads.1',names(d))],s2='reads.1'),'scores',sep='')
    # write.table(d, paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),sep='\t',eol='\n',quote=FALSE,row.names=FALSE)
    # 
    par(mfcol=c(3,2))
    for (type in c('reads','scores')){
      # for (type in c('reads')){
      for (cl in c('sum','mean','summit')){
        y=paste(cl,'bednorm',type,sep='_')
        df=d[,grep(y,names(d))]
        names(df)=unlist(lapply(strsplit(names(df),split='[.]'),FUN=function(x) x[1]))
        col=rep(c('cornflowerblue','orange'),2)
        boxplot(df,outline=FALSE,las=3,notch=TRUE,ylab=y,main=paste(g, names(grl)[i],sep=':'),col=col,border=rep(c('black','darkgrey'),each=2))
        legend('topright',bty='n',legend=c(paste0('n=',nrow(df)),paste0('nprom=',length(unique(d$genenames))))  )
      }
    }
  }
  dev.off()

}


##5) Plot Scatterplots:

pdf('Scatterplots_readsScoresInK27ac_peaks.pdf',onefile=TRUE,width=10, height=10)

plotScatter <- function (what){
  
  if(what=='reads'){
    xlim=c(0,10)
  } else {xlim=c(0,5)}
  ylim=xlim
  
  i=2
  d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),stringsAsFactors = FALSE)
  
  plot(log2(d[,paste0('Bl6.mean_bednorm_',what)]+1),log2(d[,paste0('Castaneus.mean_bednorm_',what)]+1), main=paste(names(grl)[i], what,sep=':' ), xlab=paste('log2 Bl6', what), ylab= paste('log2 Castaneus', what),xlim=xlim,ylim=ylim)
  abline(0,1,col='red')
  
  i=3
  d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),stringsAsFactors = FALSE)
  points(log2(d[,paste0('Bl6.mean_bednorm_',what)]+1),log2(d[,paste0('Castaneus.mean_bednorm_',what)]+1),col='cornflowerblue',pch=20)
  
  i=4
  d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),stringsAsFactors = FALSE)
  points(log2(d[,paste0('Bl6.mean_bednorm_',what)]+1),log2(d[,paste0('Castaneus.mean_bednorm_',what)]+1),col='orange',pch=20)
  
  legend('topleft',bty='n',pch=20,col=c('cornflowerblue','orange'),legend=names(grl)[3:4],cex=0.8)  
}


par(mfcol=c(1,2))
plotScatter('reads')
plotScatter('scores')

dev.off()

```


##################################################
##Annotate baited_genes with n interaction peaks
##n Bl6-specific/Cast-specific interaction peaks
##################################################

```{r}

##Annotate baited genes with peak numbers they interact with!

baited_genes=baited_genes[order(baited_genes$genenames)]

i=grep("K27ac_2i_all",names(grl))
d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),stringsAsFactors = FALSE)
baited_genes$n_k27acPeaks_ints=sapply(baited_genes$genenames, FUN=function(x) nrow(d[d$genename==x,]))

i=grep("k27ac_Bl6",names(grl))
d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),stringsAsFactors = FALSE)
baited_genes$n_Bl6_k27acPeaks_ints=sapply(baited_genes$genenames, FUN=function(x) nrow(d[d$genename==x,]))

i=grep("k27ac_Cast",names(grl))
d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),stringsAsFactors = FALSE)
baited_genes$n_Castaneus_k27acPeaks_ints=sapply(baited_genes$genenames, FUN=function(x) nrow(d[d$genename==x,]))


baited_genes$n_stable_k27acPeaks_ints=apply(as.data.frame(baited_genes)[,(ncol(as.data.frame(baited_genes))-2):ncol(as.data.frame(baited_genes))],1,FUN=function(x) x[1]-x[2]-x[3])

# ##only once:
# saveRDS(baited_genes,'../designDir/baited_genes_annotated.rds')


w=1.4
h=2.1
pdf('Boxplots_numberK27acpeaks_comparison_connnected_vs_closeby.pdf',width=w*7,height=h*7,pointsize=12)
##Compare how peaks in the vicinity of a gene compare to connected peaks info
p=0.001

par(mfrow=c(3,2))
baited_genes$genetype=factor(baited_genes$genetype,levels=c('down_strong','up_strong','very_stable'))
col=c('cornflowerblue','orange','grey')

main='Bl6 K27ac peaks within 1Mb'
boxplot((n_Bl6_k27acPeaks_plmin1Mb+p)/(n_k27acPeaks_plmin1Mb+p)~genetype,data=as.data.frame(baited_genes),outline=FALSE,notch=TRUE,col=col,main=main, ylab='fraction Bl-specific peaks')
legend('topright',bty='n',fill=col,legend=levels(baited_genes$genetype),cex=0.75)

main='Cast K27ac peaks within 1Mb'
boxplot((n_Castaneus_k27acPeaks_plmin1Mb+p)/(n_k27acPeaks_plmin1Mb+p)~genetype,data=as.data.frame(baited_genes),outline=FALSE,notch=TRUE,col=col,main=main, ylab='fraction Castaneus-specific peaks')
legend('topright',bty='n',fill=col,legend=levels(baited_genes$genetype),cex=0.75)


main='Bl6 K27ac peaks within 0.5Mb'
boxplot((n_Bl6_k27acPeaks_plmin0.5Mb+p)/(n_k27acPeaks_plmin0.5Mb+p)~genetype,data=as.data.frame(baited_genes),outline=FALSE,notch=TRUE,col=col,main=main, ylab='fraction Bl-specific peaks')
legend('topright',bty='n',fill=col,legend=levels(baited_genes$genetype),cex=0.75)

main='Cast K27ac peaks within 0.5Mb'
boxplot((n_Castaneus_k27acPeaks_plmin0.5Mb+p)/(n_k27acPeaks_plmin0.5Mb+p)~genetype,data=as.data.frame(baited_genes),outline=FALSE,notch=TRUE,col=col,main=main, ylab='fraction Castaneus-specific peaks')
legend('topright',bty='n',fill=col,legend=levels(baited_genes$genetype),cex=0.75)


main='Bl6 K27ac peaks within 0.25Mb'
boxplot((n_Bl6_k27acPeaks_plmin0.25Mb+p)/(n_k27acPeaks_plmin0.25Mb+p)~genetype,data=as.data.frame(baited_genes),outline=FALSE,notch=TRUE,col=col,main=main, ylab='fraction Bl-specific peaks')
legend('topright',bty='n',fill=col,legend=levels(baited_genes$genetype),cex=0.75)

main='Cast K27ac peaks within 0.25Mb'
boxplot((n_Castaneus_k27acPeaks_plmin0.25Mb+p)/(n_k27acPeaks_plmin0.25Mb+p)~genetype,data=as.data.frame(baited_genes),outline=FALSE,notch=TRUE,col=col,main=main, ylab='fraction Castaneus-specific peaks')
legend('topright',bty='n',fill=col,legend=levels(baited_genes$genetype),cex=0.75)



main='connected Bl6 K27ac peaks'
boxplot((n_Bl6_k27acPeaks_ints+p)/(n_k27acPeaks_ints+p)~genetype,data=as.data.frame(baited_genes),outline=FALSE,notch=TRUE,col=col,main=main, ylab='fraction Bl-specific peaks')
legend('topright',bty='n',fill=col,legend=levels(baited_genes$genetype),cex=0.75)

main='connected Cast K27ac peaks'
boxplot((n_Castaneus_k27acPeaks_ints+p)/(n_k27acPeaks_ints+p)~genetype,data=as.data.frame(baited_genes),outline=FALSE,notch=TRUE,col=col,main=main, ylab='fraction Castaneus-specific peaks')
legend('topright',bty='n',fill=col,legend=levels(baited_genes$genetype),cex=0.75)


main='genes w connected K27ac peaks'
boxplot((n_Bl6_k27acPeaks_ints+p)/(n_k27acPeaks_ints+p)~genetype,data=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0,]),outline=FALSE,notch=TRUE,col=col,main=main, ylab='fraction Bl-specific peaks')
legend('topright',bty='n',fill=col,legend=levels(baited_genes$genetype),cex=0.75)
legend('top',bty='n',legend=paste('n=',length(baited_genes[baited_genes$n_k27acPeaks_ints>0,])))

main='genes w connected K27ac peaks'
boxplot((n_Castaneus_k27acPeaks_ints+p)/(n_k27acPeaks_ints+p)~genetype,data=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0,]),outline=FALSE,notch=TRUE,col=col,main=main, ylab='fraction Castaneus-specific peaks')
legend('topright',bty='n',fill=col,legend=levels(baited_genes$genetype),cex=0.75)
legend('topleft',bty='n',legend=paste('n=',length(baited_genes[baited_genes$n_k27acPeaks_ints>0,])))


main='connected genes K27ac peaks within 0.25Mb'
boxplot((n_Bl6_k27acPeaks_plmin0.25Mb+p)/(n_k27acPeaks_plmin0.25Mb+p)~genetype,data=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0,]),outline=FALSE,notch=TRUE,col=col,main=main, ylab='fraction Bl-specific peaks')
legend('topright',bty='n',fill=col,legend=levels(baited_genes$genetype),cex=0.75)

main='connected genes K27ac peaks within 0.25Mb'
boxplot((n_Castaneus_k27acPeaks_plmin0.25Mb+p)/(n_k27acPeaks_plmin0.25Mb+p)~genetype,data=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0,]),outline=FALSE,notch=TRUE,col=col,main=main, ylab='fraction Castaneus-specific peaks')
legend('topright',bty='n',fill=col,legend=levels(baited_genes$genetype),cex=0.75)

dev.off()


w=2.5
h=1
pdf('Scatters_comparison_fractionBl6_Castaneus_peaks.pdf',width=w*7,height=h*7,pointsize=12)

p=0.001

##Scatter comparing fraction Bl6-specific ints with fraction Cast-specific ints
par(mfrow=c(1,3))
xlim=c(0,1)
ylim=xlim
xlab='Bl-6 specific K27ac peaks/all K27ac peaks'
ylab='Castaneus specific K27ac peaks/all K27ac peaks'
main='fraction strain-specific interacting peaks'
df=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0])
x=jitter((df$n_Bl6_k27acPeaks_ints+p)/(df$n_k27acPeaks_ints+p),amount=0.01)
y=jitter((df$n_Castaneus_k27acPeaks_ints+p)/(df$n_k27acPeaks_ints+p),amount=0.01)
plot(x,  y,xlim=xlim,ylim=ylim , col=rgb(0,0,0,alpha=0.15),cex=2,xlab=xlab,ylab=ylab,main=main)
xx=x[df$genetype=='down_strong']
yy=y[df$genetype=='down_strong']
points(xx,  yy,col=  rgb(0.2,0.2,1,alpha=0.15) ,pch=20,cex=2)
xx=x[df$genetype=='up_strong']
yy=y[df$genetype=='up_strong']
points(xx,  yy,col= rgb(1,0.3,0,alpha=0.15),pch=20,cex=2)
abline(0,1,lty=2,col='darkgreen')
legend('topright',bty='n',pch=20,col=col,legend=levels(baited_genes$genetype),title='genetype')

##Scatter comparing fraction Bl6-specific peaks in vicinity with fraction Cast-specific peaks in vicinity
xlim=c(0,1)
ylim=xlim
xlab='Bl-6 specific K27ac peaks/all K27ac peaks'
ylab='Castaneus specific K27ac peaks/all K27ac peaks'
main='fraction strain-specific peaks within 0.25Mb'
p=0.001
df=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0,])
x=jitter((df$n_Bl6_k27acPeaks_plmin0.25Mb+p)/(df$n_k27acPeaks_plmin0.25Mb+p),amount=0.01)
y=jitter((df$n_Castaneus_k27acPeaks_plmin0.25Mb+p)/(df$n_k27acPeaks_plmin0.25Mb+p),amount=0.01)
plot(x,  y,xlim=xlim,ylim=ylim , col=rgb(0,0,0,alpha=0.15),cex=2,xlab=xlab,ylab=ylab,main=main)
xx=x[df$genetype=='down_strong']
yy=y[df$genetype=='down_strong']
points(xx,  yy,col=  rgb(0.2,0.2,1,alpha=0.15) ,pch=20,cex=2)
xx=x[df$genetype=='up_strong']
yy=y[df$genetype=='up_strong']
points(xx,  yy,col= rgb(1,0.3,0,alpha=0.15),pch=20,cex=2)
abline(0,1,lty=2,col='darkgreen')
legend('topright',bty='n',pch=20,col=col,legend=levels(baited_genes$genetype),title='genetype')

##Scatter comparing fraction Bl6-specific peaks in vicinity with fraction Cast-specific peaks in vicinity (1Mb)
xlim=c(0,1)
ylim=xlim
xlab='Bl-6 specific K27ac peaks/all K27ac peaks'
ylab='Castaneus specific K27ac peaks/all K27ac peaks'
main='fraction strain-specific peaks within 1Mb'
p=0.001
df=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0,])
x=jitter((df$n_Bl6_k27acPeaks_plmin1Mb+p)/(df$n_k27acPeaks_plmin1Mb+p),amount=0.01)
y=jitter((df$n_Castaneus_k27acPeaks_plmin1Mb+p)/(df$n_k27acPeaks_plmin1Mb+p),amount=0.01)
plot(x,  y,xlim=xlim,ylim=ylim , col=rgb(0,0,0,alpha=0.15),cex=2,xlab=xlab,ylab=ylab,main=main)
xx=x[df$genetype=='down_strong']
yy=y[df$genetype=='down_strong']
points(xx,  yy,col=  rgb(0.2,0.2,1,alpha=0.15) ,pch=20,cex=2)
xx=x[df$genetype=='up_strong']
yy=y[df$genetype=='up_strong']
points(xx,  yy,col= rgb(1,0.3,0,alpha=0.15),pch=20,cex=2)
abline(0,1,lty=2,col='darkgreen')
legend('topright',bty='n',pch=20,col=col,legend=levels(baited_genes$genetype),title='genetype')


##Scatter comparing fraction stable ints with fraction stable peaks in vicinity
par(mfrow=c(1,3))
xlim=c(0,1)
ylim=xlim
xlab='stable K27ac peaks/all K27ac peaks interactions'
ylab='stable K27ac peaks/all K27ac peaks within 0.25Mb'
main='fraction stable peaks'
p=0.001
df=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0,])
x=jitter((df$n_stable_k27acPeaks_ints+p)/(df$n_k27acPeaks_ints+p),amount=0.01)
y=jitter((df$n_stable_k27acPeaks_plmin0.25Mb+p)/(df$n_k27acPeaks_plmin0.25Mb+p),amount=0.01)
plot(x,  y,xlim=xlim,ylim=ylim , col=rgb(0,0,0,alpha=0.15),cex=2,xlab=xlab,ylab=ylab,main=main)
xx=x[df$genetype=='down_strong']
yy=y[df$genetype=='down_strong']
points(xx,  yy,col=  rgb(0.2,0.2,1,alpha=0.15) ,pch=20,cex=2)
xxx=x[df$genetype=='up_strong']
yyy=y[df$genetype=='up_strong']
points(xxx,  yyy,col= rgb(1,0.3,0,alpha=0.15),pch=20,cex=2)
abline(0,1,lty=2,col='darkgreen')
legend('topright',bty='n',pch=20,col=col,legend=levels(baited_genes$genetype),title='genetype')

##boxplots looking at stable interactions:
xxxx=x[df$genetype=='very_stable']
yyyy=y[df$genetype=='very_stable']
boxplot(xx,xxx,xxxx,notch=TRUE,main='fraction ints w stable peaks',col=col)
boxplot(yy,yyy,yyyy,notch=TRUE,main='fraction stable peaks within 0.25Mb',col=col)


##Scatter comparing fraction stable ints with fraction stable peaks in vicinity (1Mb)
par(mfrow=c(1,3))
xlim=c(0,1)
ylim=xlim
xlab='stable K27ac peaks/all K27ac peaks interactions'
ylab='stable K27ac peaks/all K27ac peaks within 1Mb'
main='fraction stable peaks'
p=0.001
df=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0,])
x=jitter((df$n_stable_k27acPeaks_ints+p)/(df$n_k27acPeaks_ints+p),amount=0.01)
y=jitter((df$n_stable_k27acPeaks_plmin1Mb+p)/(df$n_k27acPeaks_plmin1Mb+p),amount=0.01)
plot(x,  y,xlim=xlim,ylim=ylim , col=rgb(0,0,0,alpha=0.15),cex=2,xlab=xlab,ylab=ylab,main=main)
xx=x[df$genetype=='down_strong']
yy=y[df$genetype=='down_strong']
points(xx,  yy,col=  rgb(0.2,0.2,1,alpha=0.15) ,pch=20,cex=2)
xxx=x[df$genetype=='up_strong']
yyy=y[df$genetype=='up_strong']
points(xxx,  yyy,col= rgb(1,0.3,0,alpha=0.15),pch=20,cex=2)
abline(0,1,lty=2,col='darkgreen')
legend('topright',bty='n',pch=20,col=col,legend=levels(baited_genes$genetype),title='genetype')

##boxplots looking at stable interactions:
xxxx=x[df$genetype=='very_stable']
yyyy=y[df$genetype=='very_stable']
boxplot(xx,xxx,xxxx,notch=TRUE,main='fraction ints w stable peaks',col=col)
boxplot(yy,yyy,yyyy,notch=TRUE,main='fraction stable peaks within 1Mb',col=col)

dev.off()
```



##calculate the effect each interacting peak should have on gene expression

##Use the formula from Jesse Engreitz (K27ac RPM * Cap-C signal) 

##To avoid one of the measures skewing the results, always normalize to the mean signal of both, Cap-C and K27ac

```{r}

# ##Normalize signal of K27ac enrichment under K27ac peaks to mean enrichment
# k27ac$Bl6_pooled_normToMean=k27ac$Bl6_pooled/mean(k27ac$Bl6_pooled)
# k27ac$Castaneus_pooled_normToMean=k27ac$Castaneus_pooled/mean(k27ac$Castaneus_pooled)
# 
# ##Normalize sum of interaction for k27ac peaks (reads and scores) to mean sum of interaction
# 
# i=grep("K27ac_2i_all",names(grl))
# d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),stringsAsFactors = FALSE)
# d=d[order(d$genetype),]
# 
# ##reads
# d$Bl6.sum_bednorm_reads_normToMean=d$Bl6.sum_bednorm_reads/mean(d$Bl6.sum_bednorm_reads)
# d$Castaneus.sum_bednorm_reads_normToMean=d$Castaneus.sum_bednorm_reads/mean(d$Castaneus.sum_bednorm_reads)
# 
# ##scores
# d$Bl6.sum_bednorm_scores_normToMean=d$Bl6.sum_bednorm_scores/mean(d$Bl6.sum_bednorm_scores)
# d$Castaneus.sum_bednorm_scores_normToMean=d$Castaneus.sum_bednorm_scores/mean(d$Castaneus.sum_bednorm_scores)


##Normalize signal of K27ac enrichment under K27ac peaks to mean K27ac enrichment of both strains
k27ac=read.csv(paste0(base, 'analysis/MouseGenomes/chip/masked_SNP.SV/ResultTables/Bl6Castaneus_k27ac.Peaks.Castaneus_vs_Bl6_annotated.csv'))
k27ac$Bl6_pooled_normToMean=k27ac$Bl6_pooled/mean(c(k27ac$Bl6_pooled,k27ac$Castaneus_pooled))
k27ac$Castaneus_pooled_normToMean=k27ac$Castaneus_pooled/mean(c(k27ac$Bl6_pooled,k27ac$Castaneus_pooled))

##Normalize sum of interaction for k27ac peaks (reads and scores) to mean sum of interaction of both strains

i=grep("K27ac_2i_all",names(grl))
d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),stringsAsFactors = FALSE)
d=d[order(d$genetype),]

##reads
d$Bl6.sum_bednorm_reads_normToMean=d$Bl6.sum_bednorm_reads/mean(c(d$Bl6.sum_bednorm_reads, d$Castaneus.sum_bednorm_reads))
d$Castaneus.sum_bednorm_reads_normToMean=d$Castaneus.sum_bednorm_reads/mean(c(d$Bl6.sum_bednorm_reads, d$Castaneus.sum_bednorm_reads))

##scores
d$Bl6.sum_bednorm_scores_normToMean=d$Bl6.sum_bednorm_scores/mean(c(d$Bl6.sum_bednorm_scores,d$Castaneus.sum_bednorm_scores))
d$Castaneus.sum_bednorm_scores_normToMean=d$Castaneus.sum_bednorm_scores/mean(c(d$Bl6.sum_bednorm_scores,d$Castaneus.sum_bednorm_scores))

##annotate d with normalized k27ac enrichment counts within the peaks
gr=bed2gr(d[,c("chr_K27ac_2i_all","start_K27ac_2i_all","end_K27ac_2i_all" )])
k=bed2gr(k27ac[,c("Chr"     ,                    "Start"         ,              "End" )])
ov=as.matrix(findOverlaps(gr,k))
d=d[ov[,1],]
dd=k27ac[ov[,2],grep('normToMean',names(k27ac))]
names(dd)=paste('k27ac_withinPeak_oE',names(dd),sep='.')
d=cbind(d,dd)


##Calculate the contribution of peak to gene activity/activation
d$k27ac_peak_contribution_reads_Bl6=d$k27ac_withinPeak_oE.Bl6_pooled_normToMean*d$Bl6.sum_bednorm_reads_normToMean
d$k27ac_peak_contribution_scores_Bl6=d$k27ac_withinPeak_oE.Bl6_pooled_normToMean*d$Bl6.sum_bednorm_scores_normToMean
d$k27ac_peak_contribution_reads_Castaneus=d$k27ac_withinPeak_oE.Castaneus_pooled_normToMean*d$Castaneus.sum_bednorm_reads_normToMean
d$k27ac_peak_contribution_scores_Castaneus=d$k27ac_withinPeak_oE.Castaneus_pooled_normToMean*d$Castaneus.sum_bednorm_scores_normToMean

##Just once:
write.table(d,paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),sep='\t',eol='\n',quote=FALSE,row.names=FALSE)

pdf('Cor_peakContribution_reads_vs_scores.pdf')
par(mfrow=c(1,2))
p=0.001
plot(log2(d$k27ac_peak_contribution_reads_Bl6+p),log2(d$k27ac_peak_contribution_scores_Bl6+p),main='K27ac peak contribution Bl6')
sp=cor(log2(d$k27ac_peak_contribution_reads_Bl6+p),log2(d$k27ac_peak_contribution_scores_Bl6+p),method='spearman')
legend('topleft',legend=paste('corSp=',round(sp,3),sep=''),bty='n')
plot(log2(d$k27ac_peak_contribution_reads_Castaneus+p),log2(d$k27ac_peak_contribution_scores_Castaneus+p),main='K27ac peak contribution Castaneus')
sp=cor(log2(d$k27ac_peak_contribution_reads_Castaneus+p),log2(d$k27ac_peak_contribution_scores_Castaneus+p),method='spearman')
legend('topleft',legend=paste('corSp=',round(sp,3),sep=''),bty='n')
dev.off()


##Annotate baited genes with the sum of all contributions for each gene!
baited_genes$totalInteracting_k27ac_peak_contribution_reads_Bl6=sapply(baited_genes$genenames, FUN=function(x) sum(d[d$genename==x,'k27ac_peak_contribution_reads_Bl6']))
baited_genes$totalInteracting_k27ac_peak_contribution_scores_Bl6=sapply(baited_genes$genenames, FUN=function(x) sum(d[d$genename==x,'k27ac_peak_contribution_scores_Bl6']))
baited_genes$totalInteracting_k27ac_peak_contribution_reads_Castaneus=sapply(baited_genes$genenames, FUN=function(x) sum(d[d$genename==x,'k27ac_peak_contribution_reads_Castaneus']))
baited_genes$totalInteracting_k27ac_peak_contribution_scores_Castaneus=sapply(baited_genes$genenames, FUN=function(x) sum(d[d$genename==x,'k27ac_peak_contribution_scores_Castaneus']))

##Plot contribution to gene activity for each gene type

pdf('Boxplots_deltaPeakContribution_byGeneType.pdf')
par(mfrow=c(2,2))

p=0.001
plot(log2(baited_genes$totalInteracting_k27ac_peak_contribution_reads_Bl6+p),log2(baited_genes$totalInteracting_k27ac_peak_contribution_reads_Castaneus+p),col=rgb(0,0,0,alpha=0.15),main='reads K27ac peak contribution',xlab='Bl6',ylab='Castaneus')
points(log2(baited_genes[baited_genes$genetype=='down_strong']$totalInteracting_k27ac_peak_contribution_reads_Bl6+p),log2(baited_genes[baited_genes$genetype=='down_strong']$totalInteracting_k27ac_peak_contribution_reads_Castaneus+p),col=rgb(0.2,0.2,1,alpha=0.2),pch=20)
points(log2(baited_genes[baited_genes$genetype=='up_strong']$totalInteracting_k27ac_peak_contribution_reads_Bl6+p),log2(baited_genes[baited_genes$genetype=='up_strong']$totalInteracting_k27ac_peak_contribution_reads_Castaneus+p),col=rgb(1,0.2,0,alpha=0.2),pch=20)
abline(0,1,lty=2,col='darkgreen')

p=0.001
plot(log2(baited_genes$totalInteracting_k27ac_peak_contribution_scores_Bl6+p),log2(baited_genes$totalInteracting_k27ac_peak_contribution_scores_Castaneus+p),col=rgb(0,0,0,alpha=0.15),main='scores K27ac peak contribution',xlab='Bl6',ylab='Castaneus')
points(log2(baited_genes[baited_genes$genetype=='down_strong']$totalInteracting_k27ac_peak_contribution_scores_Bl6+p),log2(baited_genes[baited_genes$genetype=='down_strong']$totalInteracting_k27ac_peak_contribution_scores_Castaneus+p),col=rgb(0.2,0.2,1,alpha=0.2),pch=20)
points(log2(baited_genes[baited_genes$genetype=='up_strong']$totalInteracting_k27ac_peak_contribution_scores_Bl6+p),log2(baited_genes[baited_genes$genetype=='up_strong']$totalInteracting_k27ac_peak_contribution_scores_Castaneus+p),col=rgb(1,0.2,0,alpha=0.2),pch=20)
abline(0,1,lty=2,col='darkgreen')

col=c('cornflowerblue', 'orange','grey')
boxplot( log2((totalInteracting_k27ac_peak_contribution_reads_Castaneus+p)/(totalInteracting_k27ac_peak_contribution_reads_Bl6+p))~genetype,data=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0]),outline=FALSE,col=col,las=3,main='k27ac peak contribution',ylab='log2 reads contribution Castaneus/Bl6',notch=TRUE)
abline(h=0,lty=2,col='darkgreen')
boxplot( log2((totalInteracting_k27ac_peak_contribution_scores_Castaneus+p)/(totalInteracting_k27ac_peak_contribution_scores_Bl6+p))~genetype,data=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0]),outline=FALSE,col=col,las=3,main='k27ac peak contribution',ylab='log2 scores contribution Castaneus/Bl6',notch=TRUE)
abline(h=0,lty=2,col='darkgreen')

##plotting genes that have no SNP in the bait:
boxplot( log2((totalInteracting_k27ac_peak_contribution_reads_Castaneus+p)/(totalInteracting_k27ac_peak_contribution_reads_Bl6+p))~genetype,data=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0 & baited_genes$snp_TSSplmin1kb_bait==FALSE]),outline=FALSE,col=col,las=3,main='k27ac peak contribution noSNPtss',ylab='log2 reads contribution Castaneus/Bl6',notch=TRUE)
abline(h=0,lty=2,col='darkgreen')
boxplot( log2((totalInteracting_k27ac_peak_contribution_scores_Castaneus+p)/(totalInteracting_k27ac_peak_contribution_scores_Bl6+p))~genetype,data=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0 & baited_genes$snp_TSSplmin1kb_bait==FALSE]),outline=FALSE,col=col,las=3,main='k27ac peak contribution noSNPtss',ylab='log2 scores contribution Castaneus/Bl6',notch=TRUE)
abline(h=0,lty=2,col='darkgreen')


##plotting genes that have no SNP in the bait:
boxplot( log2((totalInteracting_k27ac_peak_contribution_reads_Castaneus+p)/(totalInteracting_k27ac_peak_contribution_reads_Bl6+p))~genetype,data=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0 & baited_genes$snp_TSSplmin1kb_bait]),outline=FALSE,col=col,las=3,main='k27ac peak contribution SNPtss',ylab='log2 reads contribution Castaneus/Bl6',notch=TRUE)
abline(h=0,lty=2,col='darkgreen')
boxplot( log2((totalInteracting_k27ac_peak_contribution_scores_Castaneus+p)/(totalInteracting_k27ac_peak_contribution_scores_Bl6+p))~genetype,data=as.data.frame(baited_genes[baited_genes$n_k27acPeaks_ints>0 & baited_genes$snp_TSSplmin1kb_bait]),outline=FALSE,col=col,las=3,main='k27ac peak contribution SNPtss',ylab='log2 scores contribution Castaneus/Bl6',notch=TRUE)
abline(h=0,lty=2,col='darkgreen')

dev.off()

```



##For stable genes it is clear that they do not on average compensate for a connection to a strain-specific peak by connecting to another strain-specific peak (for instance: Bl6-peak connection in Bl6 --> Castaneus-peak connection increases in Castaneus) (see aggregate peaks and scatters)

##Question: Do they increase their connection to any other active peak to compensate?

##Strategy: Scatters for each gene/peak interaction looking at reads Castaneus vs Bl6 for all stable genes and all other genes

```{r}

i=grep("K27ac_2i_all",names(grl))
d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),stringsAsFactors = FALSE)
d=d[order(d$genetype),]
gr=bed2gr(d[,c("chr_K27ac_2i_all","start_K27ac_2i_all","end_K27ac_2i_all" )])

i=grep("k27ac_Bl6",names(grl))
d[,"k27ac_Bl6_otherEnd"]=FALSE
ov=as.matrix(findOverlaps(gr,grl[[i]]))
d[ov[,1],"k27ac_Bl6_otherEnd"]=TRUE

i=grep("k27ac_Cast",names(grl))
d[,"k27ac_Castaneus_otherEnd"]=FALSE
ov=as.matrix(findOverlaps(gr,grl[[i]]))
d[ov[,1],"k27ac_Castaneus_otherEnd"]=TRUE

##Just once:
i=grep("K27ac_2i_all",names(grl))
write.table(d,paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),sep='\t',eol='\n',quote=FALSE,row.names=FALSE)


w=1
h=1.5
pdf('Scatters_Interactions_with_K27acPeaks_perGene_numbersForPeakContribution.pdf',width=w*7,height=h*7)

par(mfrow=c(3,2))

for (id in unique(d$baitID)){
  # for (id in unique(d$baitID)[1:3]){
  print(id)
 a=d[d$baitID==id,] 
 
 p=1
  x=log2(d$Bl6.mean_bednorm_reads+p)
  y=log2(d$Castaneus.mean_bednorm_reads+p)
  xlim=c(min(c(x,y)), max(c(x,y))*1.1)
  ylim=xlim
 x=log2(a$Bl6.mean_bednorm_reads+p)
 y=log2(a$Castaneus.mean_bednorm_reads+p)
 xlab=paste0('K27ac peak int log2 reads Bl6+',p)
 ylab=paste0('K27ac peak int log2 reads Castaneus+',p)
 main=paste(baited_genes[baited_genes$dpn_id==id]$genenames,baited_genes[baited_genes$dpn_id==id]$genetype )
 plot(x ,y , ylim=ylim,xlim=xlim,pch=21,col=rgb(0,0,0,alpha=0.3),cex=3,xlab=xlab,ylab=ylab,main=main )
 x=log2(a[a$k27ac_Bl6_otherEnd,]$Bl6.mean_bednorm_reads+p)
 y=log2(a[a$k27ac_Bl6_otherEnd,]$Castaneus.mean_bednorm_reads+p)
 points(x ,y , ylim=ylim,xlim=xlim,pch=20,col=rgb(0.2,0.2,1,alpha=0.3),cex=3,xlab=xlab,ylab=ylab,main=main )
 x=log2(a[a$k27ac_Castaneus_otherEnd,]$Bl6.mean_bednorm_reads+p)
 y=log2(a[a$k27ac_Castaneus_otherEnd,]$Castaneus.mean_bednorm_reads+p)
 points(x ,y , ylim=ylim,xlim=xlim,pch=20,col=rgb(1,0.3,0,alpha=0.3),cex=3,xlab=xlab,ylab=ylab,main=main )
 abline(0,1,lty=2,col='darkgreen')
 legend('topleft',col=c(rgb(0.2,0.2,1,alpha=0.3),rgb(1,0.3,0,alpha=0.3)),pch=20,legend=c('Bl6-specific','Castaneus-specific'),title='K27ac peak type',bty='n')
 x=log2(a$Bl6.mean_bednorm_reads+p)
 y=log2(a$Castaneus.mean_bednorm_reads+p)
 text(x-0.5,y+0.5, labels=round(a$k27ac_peak_contribution_reads_Bl6,2),col='blue')
 text(x+0.5,y-0.5, labels=round(a$k27ac_peak_contribution_reads_Castaneus,2),col='red')
 ##mark the peak with top contribution to gene activity
 maxbl6=max(a$k27ac_peak_contribution_reads_Bl6)
 maxc=max(a$k27ac_peak_contribution_reads_Castaneus)
 x=log2(a[a$k27ac_peak_contribution_reads_Bl6==maxbl6,]$Bl6.mean_bednorm_reads+p)
 y=log2(a[a$k27ac_peak_contribution_reads_Bl6==maxbl6,]$Castaneus.mean_bednorm_reads+p)
 text(x,y,labels='X',col=rgb(0,0,1,alpha=0.6),cex=2)
 x=log2(a[a$k27ac_peak_contribution_reads_Castaneus==maxc,]$Bl6.mean_bednorm_reads+p)
 y=log2(a[a$k27ac_peak_contribution_reads_Castaneus==maxc,]$Castaneus.mean_bednorm_reads+p)
 text(x,y,labels='X',col=rgb(1,0,0,alpha=0.6),cex=2)
 
  x=asinh(d$Bl6.mean_bednorm_scores)
  y=asinh(d$Castaneus.mean_bednorm_scores)
  xlim=c(min(c(x,y)), max(c(x,y))*1.1)
  ylim=xlim
 x=asinh(a$Bl6.mean_bednorm_scores)
 y=asinh(a$Castaneus.mean_bednorm_scores)
 xlab='K27ac peak int asinh scores Bl6'
 ylab='K27ac peak int asinh scores Castaneus'
 main=paste(baited_genes[baited_genes$dpn_id==id]$genenames,baited_genes[baited_genes$dpn_id==id]$genetype )
 plot(x ,y , ylim=ylim,xlim=xlim,pch=21,col=rgb(0,0,0,alpha=0.3),cex=3,xlab=xlab,ylab=ylab,main=main )
 x=asinh(a[a$k27ac_Bl6_otherEnd,]$Bl6.mean_bednorm_scores)
 y=asinh(a[a$k27ac_Bl6_otherEnd,]$Castaneus.mean_bednorm_scores)
 points(x ,y , ylim=ylim,xlim=xlim,pch=20,col=rgb(0.2,0.2,1,alpha=0.3),cex=3,xlab=xlab,ylab=ylab,main=main )
 x=asinh(a[a$k27ac_Castaneus_otherEnd,]$Bl6.mean_bednorm_scores)
 y=asinh(a[a$k27ac_Castaneus_otherEnd,]$Castaneus.mean_bednorm_scores)
 points(x ,y , ylim=ylim,xlim=xlim,pch=20,col=rgb(1,0.3,0,alpha=0.3),cex=3,xlab=xlab,ylab=ylab,main=main )
 abline(0,1,lty=2,col='darkgreen')
 legend('topleft',col=c(rgb(0.2,0.2,1,alpha=0.3),rgb(1,0.3,0,alpha=0.3)),pch=20,legend=c('Bl6-specific','Castaneus-specific'),title='K27ac peak type',bty='n')
 x=asinh(a$Bl6.mean_bednorm_scores)
 y=asinh(a$Castaneus.mean_bednorm_scores)
 text(x-0.2,y+0.2, labels=round(a$k27ac_peak_contribution_scores_Bl6,2),col='blue')
 text(x+0.2,y-0.2, labels=round(a$k27ac_peak_contribution_scores_Castaneus,2),col='red')
 ##mark the peak with top contribution to gene activity
 maxbl6=max(a$k27ac_peak_contribution_scores_Bl6)
 maxc=max(a$k27ac_peak_contribution_scores_Castaneus)
 x=asinh(a[a$k27ac_peak_contribution_scores_Bl6==maxbl6,]$Bl6.mean_bednorm_scores)
 y=asinh(a[a$k27ac_peak_contribution_scores_Bl6==maxbl6,]$Castaneus.mean_bednorm_scores)
 text(x,y,labels='X',col=rgb(0,0,1,alpha=0.6),cex=2)
 x=asinh(a[a$k27ac_peak_contribution_scores_Castaneus==maxc,]$Bl6.mean_bednorm_scores)
 y=asinh(a[a$k27ac_peak_contribution_scores_Castaneus==maxc,]$Castaneus.mean_bednorm_scores)
 text(x,y,labels='X',col=rgb(1,0,0,alpha=0.6),cex=2)
 
 
}

##############################################
##All gene interactions:
##############################################

a=d
main='all gene interactions'

p=1
  x=log2(d$Bl6.mean_bednorm_reads+p)
  y=log2(d$Castaneus.mean_bednorm_reads+p)
  xlim=c(min(c(x,y)), max(c(x,y)))
  ylim=xlim
 x=log2(a$Bl6.mean_bednorm_reads+p)
 y=log2(a$Castaneus.mean_bednorm_reads+p)
 xlab=paste0('K27ac peak int log2 reads Bl6+',p)
 ylab=paste0('K27ac peak int log2 reads Castaneus+',p)
 plot(x ,y , ylim=ylim,xlim=xlim,pch=21,col=rgb(0,0,0,alpha=0.3),cex=1.5,xlab=xlab,ylab=ylab,main=main )
 x=log2(a[a$k27ac_Bl6_otherEnd,]$Bl6.mean_bednorm_reads+p)
 y=log2(a[a$k27ac_Bl6_otherEnd,]$Castaneus.mean_bednorm_reads+p)
 points(x ,y , ylim=ylim,xlim=xlim,pch=20,col=rgb(0.2,0.2,1,alpha=0.3),cex=1.5,xlab=xlab,ylab=ylab,main=main )
 x=log2(a[a$k27ac_Castaneus_otherEnd,]$Bl6.mean_bednorm_reads+p)
 y=log2(a[a$k27ac_Castaneus_otherEnd,]$Castaneus.mean_bednorm_reads+p)
 points(x ,y , ylim=ylim,xlim=xlim,pch=20,col=rgb(1,0.3,0,alpha=0.3),cex=1.5,xlab=xlab,ylab=ylab,main=main )
 abline(0,1,lty=2,col='darkgreen')
 legend('topleft',col=c(rgb(0.2,0.2,1,alpha=0.3),rgb(1,0.3,0,alpha=0.3)),pch=20,legend=c('Bl6-specific','Castaneus-specific'),title='K27ac peak type',bty='n')
 
 x=asinh(d$Bl6.mean_bednorm_scores)
  y=asinh(d$Castaneus.mean_bednorm_scores)
  xlim=c(min(c(x,y)), max(c(x,y)))
  ylim=xlim
 x=asinh(a$Bl6.mean_bednorm_scores)
 y=asinh(a$Castaneus.mean_bednorm_scores)
 xlab='K27ac peak int asinh scores Bl6'
 ylab='K27ac peak int asinh scores Castaneus'
 plot(x ,y , ylim=ylim,xlim=xlim,pch=21,col=rgb(0,0,0,alpha=0.3),cex=1.5,xlab=xlab,ylab=ylab,main=main )
 x=asinh(a[a$k27ac_Bl6_otherEnd,]$Bl6.mean_bednorm_scores)
 y=asinh(a[a$k27ac_Bl6_otherEnd,]$Castaneus.mean_bednorm_scores)
 points(x ,y , ylim=ylim,xlim=xlim,pch=20,col=rgb(0.2,0.2,1,alpha=0.3),cex=1.5,xlab=xlab,ylab=ylab,main=main )
 x=asinh(a[a$k27ac_Castaneus_otherEnd,]$Bl6.mean_bednorm_scores)
 y=asinh(a[a$k27ac_Castaneus_otherEnd,]$Castaneus.mean_bednorm_scores)
 points(x ,y , ylim=ylim,xlim=xlim,pch=20,col=rgb(1,0.3,0,alpha=0.3),cex=1.5,xlab=xlab,ylab=ylab,main=main )
 abline(0,1,lty=2,col='darkgreen')
 legend('topleft',col=c(rgb(0.2,0.2,1,alpha=0.3),rgb(1,0.3,0,alpha=0.3)),pch=20,legend=c('Bl6-specific','Castaneus-specific'),title='K27ac peak type',bty='n')
 

 ####################################
 ##Different gene type interactions
 ####################################
 
 for (typ in unique(d$genetype)){
   a=d[d$genetype==typ,]
  main=paste(typ, 'gene interactions')

p=1
  x=log2(d$Bl6.mean_bednorm_reads+p)
  y=log2(d$Castaneus.mean_bednorm_reads+p)
  xlim=c(min(c(x,y)), max(c(x,y)))
  ylim=xlim
 x=log2(a$Bl6.mean_bednorm_reads+p)
 y=log2(a$Castaneus.mean_bednorm_reads+p)
 xlab=paste0('K27ac peak int log2 reads Bl6+',p)
 ylab=paste0('K27ac peak int log2 reads Castaneus+',p)
 plot(x ,y , ylim=ylim,xlim=xlim,pch=21,col=rgb(0,0,0,alpha=0.3),cex=1.5,xlab=xlab,ylab=ylab,main=main )
 x=log2(a[a$k27ac_Bl6_otherEnd,]$Bl6.mean_bednorm_reads+p)
 y=log2(a[a$k27ac_Bl6_otherEnd,]$Castaneus.mean_bednorm_reads+p)
 points(x ,y , ylim=ylim,xlim=xlim,pch=20,col=rgb(0.2,0.2,1,alpha=0.3),cex=1.5,xlab=xlab,ylab=ylab,main=main )
 x=log2(a[a$k27ac_Castaneus_otherEnd,]$Bl6.mean_bednorm_reads+p)
 y=log2(a[a$k27ac_Castaneus_otherEnd,]$Castaneus.mean_bednorm_reads+p)
 points(x ,y , ylim=ylim,xlim=xlim,pch=20,col=rgb(1,0.3,0,alpha=0.3),cex=1.5,xlab=xlab,ylab=ylab,main=main )
 abline(0,1,lty=2,col='darkgreen')
 legend('topleft',col=c(rgb(0.2,0.2,1,alpha=0.3),rgb(1,0.3,0,alpha=0.3)),pch=20,legend=c('Bl6-specific','Castaneus-specific'),title='K27ac peak type',bty='n')
 
 x=asinh(d$Bl6.mean_bednorm_scores)
  y=asinh(d$Castaneus.mean_bednorm_scores)
  xlim=c(min(c(x,y)), max(c(x,y)))
  ylim=xlim
 x=asinh(a$Bl6.mean_bednorm_scores)
 y=asinh(a$Castaneus.mean_bednorm_scores)
 xlab='K27ac peak int asinh scores Bl6'
 ylab='K27ac peak int asinh scores Castaneus'
 plot(x ,y , ylim=ylim,xlim=xlim,pch=21,col=rgb(0,0,0,alpha=0.3),cex=1.5,xlab=xlab,ylab=ylab,main=main )
 x=asinh(a[a$k27ac_Bl6_otherEnd,]$Bl6.mean_bednorm_scores)
 y=asinh(a[a$k27ac_Bl6_otherEnd,]$Castaneus.mean_bednorm_scores)
 points(x ,y , ylim=ylim,xlim=xlim,pch=20,col=rgb(0.2,0.2,1,alpha=0.3),cex=1.5,xlab=xlab,ylab=ylab,main=main )
 x=asinh(a[a$k27ac_Castaneus_otherEnd,]$Bl6.mean_bednorm_scores)
 y=asinh(a[a$k27ac_Castaneus_otherEnd,]$Castaneus.mean_bednorm_scores)
 points(x ,y , ylim=ylim,xlim=xlim,pch=20,col=rgb(1,0.3,0,alpha=0.3),cex=1.5,xlab=xlab,ylab=ylab,main=main )
 abline(0,1,lty=2,col='darkgreen')
 legend('topleft',col=c(rgb(0.2,0.2,1,alpha=0.3),rgb(1,0.3,0,alpha=0.3)),pch=20,legend=c('Bl6-specific','Castaneus-specific'),title='K27ac peak type',bty='n')
 
 }
 
 
dev.off()



#######################################################
##In a second analysis:
##Compare whether a gene takes its primary
##enhancer or its primary enhancer starts contributing
##more or less to gene activity
#######################################################

##CAlculate the change in contribution of primary enhancer in Bl6:

i=grep("K27ac_2i_all",names(grl))
d=read.delim(paste0('oneGeneOnePeak_',names(grl)[i],'_interacting_readsScoresInPeak.txt'),stringsAsFactors = FALSE)
d=d[order(d$genetype),]
p=0.001

f=function(x,celline='Bl6',what='reads',returnwhich='Bl6'){
 a=d[d$genename==x,]
 if (nrow(a)>0){
   ci=paste('k27ac_peak_contribution',what,celline,sep='_')
   m=max(a[,ci])
   y=a[a[,ci] ==m ,paste('k27ac_peak_contribution',what,returnwhich,sep='_')]
 } else {y=0}
 return(y)
}


##Get contributions of the primary enhancers!

baited_genes$Bl6.primaryE_peakContribution_reads_Bl6=sapply(baited_genes$genenames, FUN=f,celline='Bl6',what='reads',returnwhich='Bl6'   )

baited_genes$Bl6.primaryE_peakContribution_reads_Castaneus=sapply(baited_genes$genenames, FUN=f,celline='Bl6',what='reads',returnwhich='Castaneus'   )

baited_genes$Castaneus.primaryE_peakContribution_reads_Bl6=sapply(baited_genes$genenames, FUN=f,celline='Castaneus',what='reads',returnwhich='Bl6'   )

baited_genes$Castaneus.primaryE_peakContribution_reads_Castaneus=sapply(baited_genes$genenames, FUN=f,celline='Castaneus',what='reads',returnwhich='Castaneus'   )


baited_genes$Bl6.primaryE_peakContribution_scores_Bl6=sapply(baited_genes$genenames, FUN=f,celline='Bl6',what='scores',returnwhich='Bl6'   )

baited_genes$Bl6.primaryE_peakContribution_scores_Castaneus=sapply(baited_genes$genenames, FUN=f,celline='Bl6',what='scores',returnwhich='Castaneus'   )

baited_genes$Castaneus.primaryE_peakContribution_scores_Bl6=sapply(baited_genes$genenames, FUN=f,celline='Castaneus',what='scores',returnwhich='Bl6'   )

baited_genes$Castaneus.primaryE_peakContribution_scores_Castaneus=sapply(baited_genes$genenames, FUN=f,celline='Castaneus',what='scores',returnwhich='Castaneus'   )

##Now get information about the primary enhancer switch

f=function(x,what='reads'){
 a=d[d$genename==x,]
 if (nrow(a)>0){
   
   ci=paste('k27ac_peak_contribution',what,'Bl6',sep='_')
   m=max(a[,ci])
   y=a[a[,ci] ==m ,]
   
   ci=paste('k27ac_peak_contribution',what,'Castaneus',sep='_')
   m=max(a[,ci])
   z=a[a[,ci] ==m ,]
   
   return(y$intID!=z$intID)
   
 } else {return(FALSE)}
 
}

baited_genes$primaryE_switch_reads=sapply(baited_genes$genenames, FUN=f,what='reads' )
baited_genes$primaryE_switch_scores=sapply(baited_genes$genenames, FUN=f,what='scores' )

##Get information about the identity of the primary enhancer:

f=function(x,celline='Bl6',what='reads'){
 a=d[d$genename==x,]
 if (nrow(a)>0){
   ci=paste('k27ac_peak_contribution',what,celline,sep='_')
   m=max(a[,ci])
   y='stable'
   if(a[a[,ci] ==m ,"k27ac_Bl6_otherEnd"]){
     y='Bl6_specific'
   } else if (a[a[,ci] ==m ,"k27ac_Castaneus_otherEnd"]){
     y='Castaneus_specific' 
   }
 } else {y='none'}
 return(y)
}

baited_genes$primaryEtype_Bl6_reads=sapply(baited_genes$genenames, FUN=f,what='reads',celline='Bl6' )
baited_genes$primaryEtype_Castaneus_reads=sapply(baited_genes$genenames, FUN=f,what='reads',celline='Castaneus' )
baited_genes$primaryEtype_Bl6_scores=sapply(baited_genes$genenames, FUN=f,what='scores',celline='Bl6' )
baited_genes$primaryEtype_Castaneus_scores=sapply(baited_genes$genenames, FUN=f,what='scores',celline='Castaneus' )


##Get information about the interaction strengths of the primary enhancer:

f=function(x,celline='Bl6',what='reads',whichlinetoreport='Bl6'){
 a=d[d$genename==x,]
 if (nrow(a)>0){
   ci=paste('k27ac_peak_contribution',what,celline,sep='_')
   m=max(a[,ci])
   y=a[a[,ci] ==m ,paste0(whichlinetoreport,'.sum_bednorm_',what)]
 } else {y=0}
 return(y)
}

baited_genes$primaryE_Bl6_int_readsBl6=sapply(baited_genes$genenames, FUN=f,what='reads',celline='Bl6' ,whichlinetoreport='Bl6')
baited_genes$primaryE_Bl6_int_readsCastaneus=sapply(baited_genes$genenames, FUN=f,what='reads',celline='Bl6' ,whichlinetoreport='Castaneus')

baited_genes$primaryE_Castaneus_int_readsBl6=sapply(baited_genes$genenames, FUN=f,what='reads',celline='Castaneus' ,whichlinetoreport='Bl6')
baited_genes$primaryE_Castaneus_int_readsCastaneus=sapply(baited_genes$genenames, FUN=f,what='reads',celline='Castaneus' ,whichlinetoreport='Castaneus')

baited_genes$primaryE_Bl6_int_scoresBl6=sapply(baited_genes$genenames, FUN=f,what='scores',celline='Bl6' ,whichlinetoreport='Bl6')
baited_genes$primaryE_Bl6_int_scoresCastaneus=sapply(baited_genes$genenames, FUN=f,what='scores',celline='Bl6' ,whichlinetoreport='Castaneus')

baited_genes$primaryE_Castaneus_int_scoresBl6=sapply(baited_genes$genenames, FUN=f,what='scores',celline='Castaneus',whichlinetoreport='Bl6' )
baited_genes$primaryE_Castaneus_int_scoresCastaneus=sapply(baited_genes$genenames, FUN=f,what='scores',celline='Castaneus',whichlinetoreport='Castaneus' )

###########################################
##Plots:
###########################################


scatters <- function(celline,what){
  p=0.001
  xlab='Bl6'
  ylab='Castaneus'

  x=log2(as.data.frame(baited_genes)[,paste0(celline,'.primaryE_peakContribution_',what,'_Bl6')] +p)
  y=log2(as.data.frame(baited_genes)[,paste0(celline,'.primaryE_peakContribution_',what,'_Castaneus')] +p) 
  
  baited_genes$change=y-x
  
  xlim=c(min(x,y),max(x,y))
  ylim=xlim
  main=paste0(celline, '.primaryE contribution, ',what)
  plot( x,y , cex=2, col=rgb(0,0,0,alpha=0.2),ylim=ylim,xlim=xlim,main=main,xlab=xlab,ylab=ylab  )
  x=log2(as.data.frame(baited_genes)[as.data.frame(baited_genes)$genetype=='down_strong',paste0(celline,'.primaryE_peakContribution_',what,'_Bl6')] +p)
  y=log2(as.data.frame(baited_genes)[as.data.frame(baited_genes)$genetype=='down_strong',paste0(celline,'.primaryE_peakContribution_',what,'_Castaneus')] +p) 
  points( x,y , cex=2, pch=20, col=rgb(0.2,0.2,1,alpha=0.2))
  x=log2(as.data.frame(baited_genes)[as.data.frame(baited_genes)$genetype=='up_strong',paste0(celline,'.primaryE_peakContribution_',what,'_Bl6')] +p)
  y=log2(as.data.frame(baited_genes)[as.data.frame(baited_genes)$genetype=='up_strong',paste0(celline,'.primaryE_peakContribution_',what,'_Castaneus')] +p) 
  points( x,y , cex=2, pch=20, col=rgb(1,0.2,0,alpha=0.2))
  abline(0,1,lty=2,col='darkgreen')
  legend('topleft',col=c(rgb(0.2,0.2,1,alpha=0.3),rgb(1,0.3,0,alpha=0.3)),pch=20,legend=c('Bl6-up','Castaneus-up'),title='gene type',bty='n')
  
  boxplot(change~genetype,data=as.data.frame(baited_genes)[baited_genes$n_k27acPeaks_ints>0,], col=c('cornflowerblue','orange','grey'),outline=FALSE,notch=TRUE,main=main,ylab='Delta primaryE contribution log2(Cast/Bl6)')
  abline(h=0,lty=2,col='darkgreen')
}


h=1
w=1
pdf('PrimaryE_regulation.pdf',onefile=TRUE,width=w*7,height=h*7)

par(mfrow=c(2,2))
scatters('Bl6','reads')
scatters('Castaneus','reads')
scatters('Bl6','scores')
scatters('Castaneus','scores')

##Plot fraction of genes with a switch in enhancer!
col=c('cornflowerblue','orange','grey','white')

tab=table(baited_genes[baited_genes$n_k27acPeaks_ints>0]$primaryE_switch_reads, baited_genes[baited_genes$n_k27acPeaks_ints>0]$genetype)
tab=cbind(tab,all=rowSums(tab))
plotNormBarplot_sel(tab,main='primaryE switch (contribution: reads)',sel='TRUE',ylab='fraction genes with primaryE switch',col=col)
tab=table(baited_genes[baited_genes$n_k27acPeaks_ints>0]$primaryE_switch_scores, baited_genes[baited_genes$n_k27acPeaks_ints>0]$genetype)
tab=cbind(tab,all=rowSums(tab))
plotNormBarplot_sel(tab,main='primaryE switch (contribution: scores)',sel='TRUE',ylab='fraction genes with primaryE switch',col=col)

##Plot genes with strongest enhancer being K27ac-specific

par(mfrow=c(2,2))

col=c('cornflowerblue','orange','grey','white')

tab=table(baited_genes[baited_genes$n_k27acPeaks_ints>0]$primaryEtype_Bl6_reads, baited_genes[baited_genes$n_k27acPeaks_ints>0]$genetype)
tab=cbind(tab,all=rowSums(tab))
# sel='Bl6_specific'
# plotNormBarplot_sel(tab,main='primaryE in Bl6 (contribution: reads)',sel=sel,ylab=paste('fraction genes with',sel,'primaryE'),col=col)
plotNormBarplot(tab,main='primaryE identity in Bl6 (contribution: reads)',ylab=paste('fraction genes'),col=col)

tab=table(baited_genes[baited_genes$n_k27acPeaks_ints>0]$primaryEtype_Bl6_scores, baited_genes[baited_genes$n_k27acPeaks_ints>0]$genetype)
tab=cbind(tab,all=rowSums(tab))
plotNormBarplot(tab,main='primaryE identity in Bl6 (contribution: scores)',ylab=paste('fraction genes'),col=col)

tab=table(baited_genes[baited_genes$n_k27acPeaks_ints>0]$primaryEtype_Castaneus_reads, baited_genes[baited_genes$n_k27acPeaks_ints>0]$genetype)
tab=cbind(tab,all=rowSums(tab))
plotNormBarplot(tab,main='primaryE identity in Castaneus (contribution: reads)',ylab=paste('fraction genes'),col=col)

tab=table(baited_genes[baited_genes$n_k27acPeaks_ints>0]$primaryEtype_Castaneus_scores, baited_genes[baited_genes$n_k27acPeaks_ints>0]$genetype)
tab=cbind(tab,all=rowSums(tab))
plotNormBarplot(tab,main='primaryE identity in Castaneus (contribution: scores)',ylab=paste('fraction genes'),col=col)

##Plot change in interaction strengths (sum reads)
par(mfrow=c(2,2))

p=1
b=baited_genes[baited_genes$n_k27acPeaks_ints>0]

x=log2(b$primaryE_Bl6_int_readsBl6+p)
y=log2(b$primaryE_Bl6_int_readsCastaneus+p)
main='primaryE.Bl6 interaction reads.sum'
xlab='Bl6'
ylab='Castaneus'
plot(x,y, main=main,xlab=xlab,ylab=ylab,cex=2,col=rgb(0,0,0,alpha=0.2))
abline(0,1,col='darkgreen')
points(x[b$genetype=='down_strong'], y[b$genetype=='down_strong'],col=rgb(0.2,0.2,1,alpha=0.2),pch=20,cex=2)
points(x[b$genetype=='up_strong'], y[b$genetype=='up_strong'],col=rgb(1,0.2,0,alpha=0.2),pch=20,cex=2)

boxplot(y-x~b$genetype,col=col,outline=FALSE,notch=TRUE,xlab='genetype',main=main,ylab='change interaction Cast/Bl6')
abline(h=0,lty=2,col='darkgreen')

x=log2(b$primaryE_Bl6_int_scoresBl6+p)
y=log2(b$primaryE_Bl6_int_scoresCastaneus+p)
main='primaryE.Bl6 interaction scores.sum'
xlab='Bl6'
ylab='Castaneus'
plot(x,y, main=main,xlab=xlab,ylab=ylab,cex=2,col=rgb(0,0,0,alpha=0.2))
abline(0,1,col='darkgreen')
points(x[b$genetype=='down_strong'], y[b$genetype=='down_strong'],col=rgb(0.2,0.2,1,alpha=0.2),pch=20,cex=2)
points(x[b$genetype=='up_strong'], y[b$genetype=='up_strong'],col=rgb(1,0.2,0,alpha=0.2),pch=20,cex=2)

boxplot(y-x~b$genetype,col=col,outline=FALSE,notch=TRUE,xlab='genetype',main=main,ylab='change interaction Cast/Bl6')
abline(h=0,lty=2,col='darkgreen')

x=log2(b$primaryE_Castaneus_int_readsBl6+p)
y=log2(b$primaryE_Castaneus_int_readsCastaneus+p)
main='primaryE.Castaneus interaction reads.sum'
xlab='Bl6'
ylab='Castaneus'
plot(x,y, main=main,xlab=xlab,ylab=ylab,cex=2,col=rgb(0,0,0,alpha=0.2))
abline(0,1,col='darkgreen')
points(x[b$genetype=='down_strong'], y[b$genetype=='down_strong'],col=rgb(0.2,0.2,1,alpha=0.2),pch=20,cex=2)
points(x[b$genetype=='up_strong'], y[b$genetype=='up_strong'],col=rgb(1,0.2,0,alpha=0.2),pch=20,cex=2)

boxplot(y-x~b$genetype,col=col,outline=FALSE,notch=TRUE,xlab='genetype',main=main,ylab='change interaction Cast/Bl6')
abline(h=0,lty=2,col='darkgreen')

x=log2(b$primaryE_Castaneus_int_scoresBl6+p)
y=log2(b$primaryE_Castaneus_int_scoresCastaneus+p)
main='primaryE.Castaneus interaction scores.sum'
xlab='Bl6'
ylab='Castaneus'
plot(x,y, main=main,xlab=xlab,ylab=ylab,cex=2,col=rgb(0,0,0,alpha=0.2))
abline(0,1,col='darkgreen')
points(x[b$genetype=='down_strong'], y[b$genetype=='down_strong'],col=rgb(0.2,0.2,1,alpha=0.2),pch=20,cex=2)
points(x[b$genetype=='up_strong'], y[b$genetype=='up_strong'],col=rgb(1,0.2,0,alpha=0.2),pch=20,cex=2)

boxplot(y-x~b$genetype,col=col,outline=FALSE,notch=TRUE,xlab='genetype',main=main,ylab='change interaction Cast/Bl6')
abline(h=0,lty=2,col='darkgreen')

dev.off()       

```




















```{r}

##1) Determine hind frags that overlap the center of k27ac peaks

##for each k27ac peak: determine the DpnII fragment its center overlaps:
b=bed2gr(k27ac_all[,c("Chr"                  ,"Start"  ,               "End")])
b=resize(b,1,fix='center')

ov=as.matrix(findOverlaps(b,hind))
k27ac_all=k27ac_all[ov[,1],]
k27ac_all$dpn_id=hind[ov[,2]]$id

##for each interaction determine the name of k27ac peak the gene interacts
b=bed2gr(k27ac_all[,c("Chr"                  ,"Start"  ,               "End")])
i=bed2gr(ints[,c("seqnames_otherEnd" ,        "start_otherEnd",            "end_otherEnd" )])

ov=as.matrix(findOverlaps(i,b))
ints.k=ints[ov[,1],]
ints.k$dpn_id_k27acPeak=k27ac_all[ov[,2],]$dpn_id
ints.k$name_k27acPeak=k27ac_all[ov[,2],]$X
ints.k$LFC_k27acPeak=k27ac_all[ov[,2],]$LFC
ints.k$padj_k27acPeak=k27ac_all[ov[,2],]$padj
ints.k$Bl6_pooled_fpkm_k27acPeak=k27ac_all[ov[,2],]$Bl6_pooled_fpkm
ints.k$Castaneus_pooled_fpkm_k27acPeak=k27ac_all[ov[,2],]$Castaneus_pooled_fpkm
ints.k= cbind(ints.k,k27ac_all[ov[,2],c('Chr','Start','End')])
names(ints.k)[(ncol(ints.k)-2):ncol(ints.k)]=paste(names(ints.k)[(ncol(ints.k)-2):ncol(ints.k)],'k27acPeak',sep='_')

write.table(ints.k,'ints_k27ac_interacting.txt',sep='\t',eol='\n',quote=FALSE,row.names=FALSE)

```

############################################
##create metaprofile matrices and plot them:
############################################

Use the wrapping function! 

define intID as baitID;dpn_id_k27acPeak (unique interactions only!)

define otherEndID as dpn_id_k27acPeak

First: reads

```{r}

peaks=read.delim('ints_k27ac_interacting.txt',stringsAsFactors=FALSE)

##make data frame d with clusters_refined and other end data encompassing the peak summit!
d=peaks
d$otherEndID=d$dpn_id_k27acPeak
d[,c('seqnames_otherEnd','start_otherEnd','end_otherEnd')]=d[,c("Chr_k27acPeak"  ,  "Start_k27acPeak"  ,"End_k27acPeak" )]
dim(d)
d=d[d$seqnames_bait==d$seqnames_otherEnd,]
dim(d)
d$intID=paste(d$baitID,d$dpn_id_k27acPeak,sep=';')

id=unique(d$intID)[1]
a=d[d$intID==id,]
df=a[1,]
for (id in unique(d$intID)[2:length(unique(d$intID))] ){
  df=rbind(df,d[d$intID==id,][1,])
}

dim(df)
d=df
row.names(d)=d$intID
 
ids=unique(peaks$baitID)
 
zoom=40

dis='k27ac_peak_center'

Lm=matrices_stranded_InteractionsDensityGeneral(L=L,ids=ids,zoom=zoom, unt=1,tam=2,d=d)

col=c('cornflowerblue','orange')

pdf('Metaprofile_H3K27ac_peak_interactions_colMedians.pdf',onefile=TRUE)

up=row.names(d[d$LFC_k27acPeak >= 1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE,])
down=row.names(d[d$LFC_k27acPeak <=  -1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE,])

upgenes=row.names(d[d$genetype=='up_strong',])
downgenes=row.names(d[d$genetype=='down_strong',])
stablegenes=row.names(d[d$genetype=='very_stable',])

ylim=c(0,1.2)
ylab='norm read counts'
xlab='dist to K27ac peak center [DpnII frag]'

dist=(ncol(Lm[[1]])-1)/2
x=-dist:dist

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]]),type='l',main=paste('all',nrow(Lm[[1]]), 'H3K27ac peaks interactions'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]]),col=col[2])
lines(x,colMedians(Lm[[3]]),col=col[1],lty=3)
lines(x,colMedians(Lm[[4]]),col=col[2],lty=3)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% up,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% up,]), 'Cast-up H3K27ac peaks'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% up,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% up,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% up,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions')

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% down,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% down,]), 'Bl6-up H3K27ac peaks'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% down,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% down,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% down,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% c(up,down)==FALSE,]), 'stable H3K27ac peaks'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions')


##


plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% upgenes,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% upgenes,]), 'Cast-up genes'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% upgenes,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% upgenes,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% upgenes,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% downgenes,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% downgenes,]), 'Bl6-up genes'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% downgenes,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% downgenes,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% downgenes,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% stablegenes,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% stablegenes,]), 'stable genes'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% stablegenes,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% stablegenes,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% stablegenes,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )

################################################################################
##How are interactions with K27ac peaks regulated at 
##stable genes with or without changing interactions?
################################################################################

##stable genes with downregulated peak interactions:

stablegenes=unique(d[d$genetype=='very_stable' ,]$baitID)

withdown=vector()
withup=vector()
without=vector()
for (id in stablegenes){
  a=d[d$baitID==id,]
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ){
    withdown=c(withdown,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak >= 1) ){
    withup=c(withup,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ==FALSE ){
    if ( TRUE %in% (a$LFC_k27acPeak >= 1) ==FALSE){
      without=c(without,id)
    }
  }  
}

##plot all interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown,])
ints.withup=row.names(d[d$baitID %in% withup,])
ints.without=row.names(d[d$baitID %in% without,])

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'stable with Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'stable with Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]), 'stable with stable.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.without,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.without,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.without,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )


##plot only differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak <= -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak >= 1,])

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'stable only Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'stable only Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )


##plot only non-differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak > -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak < 1,])

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'stable with Bl6.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'stable with Cast.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )




################################################################################
##How are interactions with K27ac peaks regulated at 
##down genes with or without changing interactions?
################################################################################

##down genes with peak interactions:

stablegenes=unique(d[d$genetype=='down_strong' ,]$baitID)

withdown=vector()
withup=vector()
without=vector()
for (id in stablegenes){
  a=d[d$baitID==id,]
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ){
    withdown=c(withdown,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak >= 1) ){
    withup=c(withup,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ==FALSE ){
    if ( TRUE %in% (a$LFC_k27acPeak >= 1) ==FALSE){
      without=c(without,id)
    }
  }  
}

##plot all interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown,])
ints.withup=row.names(d[d$baitID %in% withup,])
ints.without=row.names(d[d$baitID %in% without,])

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'down with Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'down with Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]), 'down with stable.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.without,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.without,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.without,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )


##plot only differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak <= -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak >= 1,])

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'down only Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'down only Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )


##plot only non-differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak > -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak < 1,])

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'down with Bl6.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'down with Cast.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )



################################################################################
##How are interactions with K27ac peaks regulated at 
##up genes with or without changing interactions?
################################################################################

##up genes with peak interactions:

stablegenes=unique(d[d$genetype=='up_strong' ,]$baitID)

withdown=vector()
withup=vector()
without=vector()
for (id in stablegenes){
  a=d[d$baitID==id,]
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ){
    withdown=c(withdown,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak >= 1) ){
    withup=c(withup,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ==FALSE ){
    if ( TRUE %in% (a$LFC_k27acPeak >= 1) ==FALSE){
      without=c(without,id)
    }
  }  
}

##plot all interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown,])
ints.withup=row.names(d[d$baitID %in% withup,])
ints.without=row.names(d[d$baitID %in% without,])

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'up with Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'up with Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]), 'up with stable.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.without,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.without,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.without,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )


##plot only differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak <= -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak >= 1,])

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'up only Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'up only Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )


##plot only non-differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak > -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak < 1,])

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'up with Bl6.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'up with Cast.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

dev.off()




##########################################
##smoothed lines:
##########################################


k=3

pdf(paste0('Metaprofile_H3K27ac_peak_interactions_colMedians_k',k,'.pdf'),onefile=TRUE)

up=row.names(d[d$LFC_k27acPeak >= 1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE,])
down=row.names(d[d$LFC_k27acPeak <=  -1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE,])

upgenes=row.names(d[d$genetype=='up_strong',])
downgenes=row.names(d[d$genetype=='down_strong',])
stablegenes=row.names(d[d$genetype=='very_stable',])

ylim=c(0,1)
ylab='norm read counts'
xlab='dist to K27ac peak center [DpnII frag]'

dist=(ncol(Lm[[1]])-1)/2
x=-dist:dist

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]])),k=k,endrule='constant'),type='l',main=paste('all',nrow(Lm[[1]]), 'H3K27ac peaks interactions'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]])),k=k,endrule='constant'),col=col[1],lty=3)
lines(x,runmean(Rle(colMedians(Lm[[4]])),k=k,endrule='constant'),col=col[2],lty=3)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% up,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% up,]), 'Cast-up H3K27ac peaks'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% up,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% up,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% up,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions')

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% down,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% down,]), 'Bl6-up H3K27ac peaks'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% down,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% down,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% down,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% c(up,down)==FALSE,]), 'stable H3K27ac peaks'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions')


##


plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% upgenes,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% upgenes,]), 'Cast-up genes'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% upgenes,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% upgenes,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% upgenes,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% downgenes,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% downgenes,]), 'Bl6-up genes'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% downgenes,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% downgenes,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% downgenes,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% stablegenes,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% stablegenes,]), 'stable genes'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% stablegenes,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% stablegenes,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% stablegenes,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )

################################################################################
##How are interactions with K27ac peaks regulated at 
##stable genes with or without changing interactions?
################################################################################

##stable genes with downregulated peak interactions:

stablegenes=unique(d[d$genetype=='very_stable' ,]$baitID)

withdown=vector()
withup=vector()
without=vector()
for (id in stablegenes){
  a=d[d$baitID==id,]
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ){
    withdown=c(withdown,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak >= 1) ){
    withup=c(withup,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ==FALSE ){
    if ( TRUE %in% (a$LFC_k27acPeak >= 1) ==FALSE){
      without=c(without,id)
    }
  }  
}

##plot all interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown,])
ints.withup=row.names(d[d$baitID %in% withup,])
ints.without=row.names(d[d$baitID %in% without,])

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'stable with Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'stable with Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]), 'stable with stable.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )


##plot only differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak <= -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak >= 1,])

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'stable only Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'stable only Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )


##plot only non-differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak > -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak < 1,])

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'stable with Bl6.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'stable with Cast.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )




################################################################################
##How are interactions with K27ac peaks regulated at 
##down genes with or without changing interactions?
################################################################################

##down genes with peak interactions:

stablegenes=unique(d[d$genetype=='down_strong' ,]$baitID)

withdown=vector()
withup=vector()
without=vector()
for (id in stablegenes){
  a=d[d$baitID==id,]
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ){
    withdown=c(withdown,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak >= 1) ){
    withup=c(withup,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ==FALSE ){
    if ( TRUE %in% (a$LFC_k27acPeak >= 1) ==FALSE){
      without=c(without,id)
    }
  }  
}

##plot all interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown,])
ints.withup=row.names(d[d$baitID %in% withup,])
ints.without=row.names(d[d$baitID %in% without,])

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'down with Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'down with Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]), 'down with stable.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )


##plot only differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak <= -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak >= 1,])

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'down only Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'down only Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )


##plot only non-differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak > -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak < 1,])

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'down with Bl6.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'down with Cast.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )



################################################################################
##How are interactions with K27ac peaks regulated at 
##up genes with or without changing interactions?
################################################################################

##up genes with peak interactions:

stablegenes=unique(d[d$genetype=='up_strong' ,]$baitID)

withdown=vector()
withup=vector()
without=vector()
for (id in stablegenes){
  a=d[d$baitID==id,]
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ){
    withdown=c(withdown,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak >= 1) ){
    withup=c(withup,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ==FALSE ){
    if ( TRUE %in% (a$LFC_k27acPeak >= 1) ==FALSE){
      without=c(without,id)
    }
  }  
}

##plot all interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown,])
ints.withup=row.names(d[d$baitID %in% withup,])
ints.without=row.names(d[d$baitID %in% without,])

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'up with Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'up with Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]), 'up with stable.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )


##plot only differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak <= -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak >= 1,])

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'up only Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'up only Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )


##plot only non-differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak > -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak < 1,])

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'up with Bl6.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'up with Cast.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

dev.off()




##load matrices:

Lm=loadMatrices(L,1)

b=Lm[[1]]
c=Lm[[2]]
bc=Lm[[3]]
cc=Lm[[4]]

```


Then: scores

```{r}

peaks=read.delim('ints_k27ac_interacting.txt',stringsAsFactors=FALSE)

##make data frame d with clusters_refined and other end data encompassing the peak summit!
d=peaks
d$otherEndID=d$dpn_id_k27acPeak
d[,c('seqnames_otherEnd','start_otherEnd','end_otherEnd')]=d[,c("Chr_k27acPeak"  ,  "Start_k27acPeak"  ,"End_k27acPeak" )]
dim(d)
d=d[d$seqnames_bait==d$seqnames_otherEnd,]
dim(d)
d$intID=paste(d$baitID,d$dpn_id_k27acPeak,sep=';')

id=unique(d$intID)[1]
a=d[d$intID==id,]
df=a[1,]
for (id in unique(d$intID)[2:length(unique(d$intID))] ){
  df=rbind(df,d[d$intID==id,][1,])
}

dim(df)
d=df
row.names(d)=d$intID
 
ids=unique(peaks$baitID)
 
zoom=40

dis='k27ac_peak_center'

Lm=matrices_stranded_InteractionsDensityGeneral(L=L,ids=ids,zoom=zoom, unt=1,tam=2,d=d,type='scores')

col=c('cornflowerblue','orange')

pdf('Metaprofile_H3K27ac_peak_interactions_colMedians_scores.pdf',onefile=TRUE)

up=row.names(d[d$LFC_k27acPeak >= 1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE,])
down=row.names(d[d$LFC_k27acPeak <=  -1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE,])

upgenes=row.names(d[d$genetype=='up_strong',])
downgenes=row.names(d[d$genetype=='down_strong',])
stablegenes=row.names(d[d$genetype=='very_stable',])

ylim=c(0,1.2)
ylab='norm read counts'
xlab='dist to K27ac peak center [DpnII frag]'

dist=(ncol(Lm[[1]])-1)/2
x=-dist:dist

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]]),type='l',main=paste('all',nrow(Lm[[1]]), 'H3K27ac peaks interactions'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]]),col=col[2])
lines(x,colMedians(Lm[[3]]),col=col[1],lty=3)
lines(x,colMedians(Lm[[4]]),col=col[2],lty=3)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% up,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% up,]), 'Cast-up H3K27ac peaks'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% up,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% up,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% up,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions')

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% down,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% down,]), 'Bl6-up H3K27ac peaks'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% down,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% down,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% down,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% c(up,down)==FALSE,]), 'stable H3K27ac peaks'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions')


##


plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% upgenes,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% upgenes,]), 'Cast-up genes'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% upgenes,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% upgenes,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% upgenes,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% downgenes,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% downgenes,]), 'Bl6-up genes'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% downgenes,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% downgenes,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% downgenes,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% stablegenes,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% stablegenes,]), 'stable genes'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% stablegenes,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% stablegenes,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% stablegenes,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )

################################################################################
##How are interactions with K27ac peaks regulated at 
##stable genes with or without changing interactions?
################################################################################

##stable genes with downregulated peak interactions:

stablegenes=unique(d[d$genetype=='very_stable' ,]$baitID)

withdown=vector()
withup=vector()
without=vector()
for (id in stablegenes){
  a=d[d$baitID==id,]
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ){
    withdown=c(withdown,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak >= 1) ){
    withup=c(withup,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ==FALSE ){
    if ( TRUE %in% (a$LFC_k27acPeak >= 1) ==FALSE){
      without=c(without,id)
    }
  }  
}

##plot all interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown,])
ints.withup=row.names(d[d$baitID %in% withup,])
ints.without=row.names(d[d$baitID %in% without,])

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'stable with Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'stable with Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]), 'stable with stable.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.without,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.without,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.without,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )


##plot only differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak <= -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak >= 1,])

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'stable only Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'stable only Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )


##plot only non-differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak > -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak < 1,])

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'stable with Bl6.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'stable with Cast.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )




################################################################################
##How are interactions with K27ac peaks regulated at 
##down genes with or without changing interactions?
################################################################################

##down genes with peak interactions:

stablegenes=unique(d[d$genetype=='down_strong' ,]$baitID)

withdown=vector()
withup=vector()
without=vector()
for (id in stablegenes){
  a=d[d$baitID==id,]
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ){
    withdown=c(withdown,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak >= 1) ){
    withup=c(withup,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ==FALSE ){
    if ( TRUE %in% (a$LFC_k27acPeak >= 1) ==FALSE){
      without=c(without,id)
    }
  }  
}

##plot all interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown,])
ints.withup=row.names(d[d$baitID %in% withup,])
ints.without=row.names(d[d$baitID %in% without,])

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'down with Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'down with Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]), 'down with stable.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.without,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.without,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.without,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )


##plot only differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak <= -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak >= 1,])

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'down only Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'down only Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )


##plot only non-differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak > -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak < 1,])

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'down with Bl6.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'down with Cast.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )



################################################################################
##How are interactions with K27ac peaks regulated at 
##up genes with or without changing interactions?
################################################################################

##up genes with peak interactions:

stablegenes=unique(d[d$genetype=='up_strong' ,]$baitID)

withdown=vector()
withup=vector()
without=vector()
for (id in stablegenes){
  a=d[d$baitID==id,]
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ){
    withdown=c(withdown,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak >= 1) ){
    withup=c(withup,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ==FALSE ){
    if ( TRUE %in% (a$LFC_k27acPeak >= 1) ==FALSE){
      without=c(without,id)
    }
  }  
}

##plot all interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown,])
ints.withup=row.names(d[d$baitID %in% withup,])
ints.without=row.names(d[d$baitID %in% without,])

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'up with Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'up with Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]), 'up with stable.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.without,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.without,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.without,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )


##plot only differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak <= -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak >= 1,])

par(mfrow=c(2,2))
plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'up only Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'up only Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )


##plot only non-differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak > -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak < 1,])

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'up with Bl6.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'up with Cast.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2])
lines(x,colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[1],lty=5)
lines(x,colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,]),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

dev.off()




##########################################
##smoothed lines:
##########################################


k=3

pdf(paste0('Metaprofile_H3K27ac_peak_interactions_colMedians_scores_k',k,'.pdf'),onefile=TRUE)

up=row.names(d[d$LFC_k27acPeak >= 1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE,])
down=row.names(d[d$LFC_k27acPeak <=  -1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE,])

upgenes=row.names(d[d$genetype=='up_strong',])
downgenes=row.names(d[d$genetype=='down_strong',])
stablegenes=row.names(d[d$genetype=='very_stable',])

ylim=c(0,1)
ylab='norm read counts'
xlab='dist to K27ac peak center [DpnII frag]'

dist=(ncol(Lm[[1]])-1)/2
x=-dist:dist

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]])),k=k,endrule='constant'),type='l',main=paste('all',nrow(Lm[[1]]), 'H3K27ac peaks interactions'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]])),k=k,endrule='constant'),col=col[1],lty=3)
lines(x,runmean(Rle(colMedians(Lm[[4]])),k=k,endrule='constant'),col=col[2],lty=3)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% up,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% up,]), 'Cast-up H3K27ac peaks'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% up,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% up,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% up,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions')

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% down,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% down,]), 'Bl6-up H3K27ac peaks'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% down,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% down,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% down,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% c(up,down)==FALSE,]), 'stable H3K27ac peaks'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% c(up,down) == FALSE,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions')


##


plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% upgenes,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% upgenes,]), 'Cast-up genes'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% upgenes,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% upgenes,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% upgenes,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% downgenes,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% downgenes,]), 'Bl6-up genes'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% downgenes,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% downgenes,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% downgenes,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% stablegenes,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% stablegenes,]), 'stable genes'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% stablegenes,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% stablegenes,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% stablegenes,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )

################################################################################
##How are interactions with K27ac peaks regulated at 
##stable genes with or without changing interactions?
################################################################################

##stable genes with downregulated peak interactions:

stablegenes=unique(d[d$genetype=='very_stable' ,]$baitID)

withdown=vector()
withup=vector()
without=vector()
for (id in stablegenes){
  a=d[d$baitID==id,]
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ){
    withdown=c(withdown,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak >= 1) ){
    withup=c(withup,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ==FALSE ){
    if ( TRUE %in% (a$LFC_k27acPeak >= 1) ==FALSE){
      without=c(without,id)
    }
  }  
}

##plot all interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown,])
ints.withup=row.names(d[d$baitID %in% withup,])
ints.without=row.names(d[d$baitID %in% without,])

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'stable with Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'stable with Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]), 'stable with stable.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )


##plot only differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak <= -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak >= 1,])

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'stable only Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'stable only Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )


##plot only non-differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak > -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak < 1,])

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'stable with Bl6.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'stable with Cast.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )




################################################################################
##How are interactions with K27ac peaks regulated at 
##down genes with or without changing interactions?
################################################################################

##down genes with peak interactions:

stablegenes=unique(d[d$genetype=='down_strong' ,]$baitID)

withdown=vector()
withup=vector()
without=vector()
for (id in stablegenes){
  a=d[d$baitID==id,]
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ){
    withdown=c(withdown,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak >= 1) ){
    withup=c(withup,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ==FALSE ){
    if ( TRUE %in% (a$LFC_k27acPeak >= 1) ==FALSE){
      without=c(without,id)
    }
  }  
}

##plot all interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown,])
ints.withup=row.names(d[d$baitID %in% withup,])
ints.without=row.names(d[d$baitID %in% without,])

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'down with Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'down with Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]), 'down with stable.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )


##plot only differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak <= -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak >= 1,])

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'down only Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'down only Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )


##plot only non-differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak > -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak < 1,])

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'down with Bl6.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'down with Cast.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )



################################################################################
##How are interactions with K27ac peaks regulated at 
##up genes with or without changing interactions?
################################################################################

##up genes with peak interactions:

stablegenes=unique(d[d$genetype=='up_strong' ,]$baitID)

withdown=vector()
withup=vector()
without=vector()
for (id in stablegenes){
  a=d[d$baitID==id,]
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ){
    withdown=c(withdown,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak >= 1) ){
    withup=c(withup,id)
  } 
  if ( TRUE %in% (a$LFC_k27acPeak <= -1) ==FALSE ){
    if ( TRUE %in% (a$LFC_k27acPeak >= 1) ==FALSE){
      without=c(without,id)
    }
  }  
}

##plot all interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown,])
ints.withup=row.names(d[d$baitID %in% withup,])
ints.without=row.names(d[d$baitID %in% without,])

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'up with Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'up with Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.without,]), 'up with stable.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.without,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions' )


##plot only differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak <= -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak >= 1,])

par(mfrow=c(2,2))
plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'up only Bl6.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'up only Cast.K27ac-ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )


##plot only non-differential interactions of such genes:

ints.withdown=row.names(d[d$baitID %in% withdown & d$LFC_k27acPeak > -1,])
ints.withup=row.names(d[d$baitID %in% withup & d$LFC_k27acPeak < 1,])

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withdown,]), 'up with Bl6.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withdown,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6 ,title='interactions' )

plot(x,runmean(Rle(colMedians(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),type='l',main=paste(nrow(Lm[[1]][row.names(Lm[[1]]) %in% ints.withup,]), 'up with Cast.K27ac-ints other ints'),ylim=ylim,col=col[1], ylab=ylab,xlab=xlab)
lines(x,runmean(Rle(colMedians(Lm[[2]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2])
lines(x,runmean(Rle(colMedians(Lm[[3]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[1],lty=5)
lines(x,runmean(Rle(colMedians(Lm[[4]][row.names(Lm[[1]]) %in% ints.withup,])),k=k,endrule='constant'),col=col[2],lty=5)
legend('topleft',lty=rep(c(1,3),each=2), legend= c('Bl6','Cast','Bl6_con','Cast_con'), col=rep(col,2), bty='n',cex=0.6,title='interactions'  )

dev.off()




##load matrices:

Lm=loadMatrices(L,1,type='scores')

b=Lm[[1]]
c=Lm[[2]]
bc=Lm[[3]]
cc=Lm[[4]]

```

##!!!Comment: heck out that control matrices for scores are done for scores, not for reads!
##!!! Do scores matrices need to be normalized to peak wt?

##Use these matrices to quantify normalized reads or scores of interaction with K27ac peaks!

```{r}

peaks=read.delim('ints_k27ac_interacting.txt',stringsAsFactors=FALSE)

##make data frame d with clusters_refined and other end data encompassing the peak summit!
d=peaks
d$otherEndID=d$dpn_id_k27acPeak
d[,c('seqnames_otherEnd','start_otherEnd','end_otherEnd')]=d[,c("Chr_k27acPeak"  ,  "Start_k27acPeak"  ,"End_k27acPeak" )]
dim(d)
d=d[d$seqnames_bait==d$seqnames_otherEnd,]
dim(d)
d$intID=paste(d$baitID,d$dpn_id_k27acPeak,sep=';')

id=unique(d$intID)[1]
a=d[d$intID==id,]
df=a[1,]
for (id in unique(d$intID)[2:length(unique(d$intID))] ){
  df=rbind(df,d[d$intID==id,][1,])
}

dim(df)
d=df
row.names(d)=d$intID
 
ids=unique(peaks$baitID)
 
zoom=40

dis='k27ac_peak_center'


##load normalized matrices:
Lms=loadMatrices(L,1,type='scores')
Lmr=loadMatrices(L,1,type='reads')

##for each interacting K27ac peak identify all fragments covered and their position relative to (zoom+1) position of the peak center:
gr=bed2gr(d[,c("Chr_k27acPeak" ,"Start_k27acPeak" , "End_k27acPeak")])
n=vector()

d[,paste0('reads_within_k27acPeak_normToPeakCenterBl6_',names(Lmr)[1])]=0
d[,paste0('reads_within_k27acPeak_normToPeakCenterBl6_',names(Lmr)[2])]=0
d[,paste0('score_within_k27acPeak_normToPeakCenterBl6_',names(Lms)[1])]=0
d[,paste0('score_within_k27acPeak_normToPeakCenterBl6_',names(Lms)[2])]=0

# for (i in 1:nrow(d)){
  for (i in 45:nrow(d)){
  print(i)
  id=d[i,]$dpn_id_k27acPeak
  low=(zoom+1)-(id-min(hind[as.matrix(findOverlaps(gr[i],hind))[,2]]$id))
  hi=(zoom+1)+max(hind[as.matrix(findOverlaps(gr[i],hind))[,2]]$id)-id
  n=c(n,length(low:hi))
  d[i,paste0('reads_within_k27acPeak_normToPeakCenterBl6_',names(Lmr)[1])]=sum(Lmr[[1]][i,][low:hi])
  d[i,paste0('reads_within_k27acPeak_normToPeakCenterBl6_',names(Lmr)[2])]=sum(Lmr[[2]][i,][low:hi])
  d[i,paste0('score_within_k27acPeak_normToPeakCenterBl6_',names(Lms)[1])]=sum(Lms[[1]][i,][low:hi])
  d[i,paste0('score_within_k27acPeak_normToPeakCenterBl6_',names(Lms)[2])]=sum(Lms[[2]][i,][low:hi])
}


d[,paste0('reads_center_k27acPeak_normToPeakCenterBl6_',names(Lmr)[1])]=0
d[,paste0('reads_center_k27acPeak_normToPeakCenterBl6_',names(Lmr)[2])]=0
d[,paste0('score_center_k27acPeak_normToPeakCenterBl6_',names(Lms)[1])]=0
d[,paste0('score_center_k27acPeak_normToPeakCenterBl6_',names(Lms)[2])]=0

for (i in 1:nrow(d)){
  # for (i in 45:nrow(d)){
  print(i)
  low=41
  hi=41
  d[i,paste0('reads_center_k27acPeak_normToPeakCenterBl6_',names(Lmr)[1])]=sum(Lmr[[1]][i,][low:hi])
  d[i,paste0('reads_center_k27acPeak_normToPeakCenterBl6_',names(Lmr)[2])]=sum(Lmr[[2]][i,][low:hi])
  d[i,paste0('score_center_k27acPeak_normToPeakCenterBl6_',names(Lms)[1])]=sum(Lms[[1]][i,][low:hi])
  d[i,paste0('score_center_k27acPeak_normToPeakCenterBl6_',names(Lms)[2])]=sum(Lms[[2]][i,][low:hi])
}



##Boxplots:

##downregulated, upregulated and stable peaks:
Ld=list(down=row.names(d[d$LFC_k27acPeak <= -1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE,]),
  up=row.names(d[d$LFC_k27acPeak >= 1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE,]),
stable=row.names(d[ abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE,]),

down.genes_down=row.names(d[d$LFC_k27acPeak <= -1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='down_strong',]),
down.genes_up=row.names(d[d$LFC_k27acPeak >= 1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='down_strong',]),
down.genes_stable=row.names(d[ abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='down_strong',]),

up.genes_down=row.names(d[d$LFC_k27acPeak <= -1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='up_strong',]),
up.genes_up=row.names(d[d$LFC_k27acPeak >= 1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='up_strong',]),
up.genes_stable=row.names(d[ abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='up_strong',]),

stable.genes_down=row.names(d[d$LFC_k27acPeak <= -1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='very_stable',]),
stable.genes_up=row.names(d[d$LFC_k27acPeak >= 1 & d$padj_k27acPeak < 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='very_stable',]),
stable.genes_stable=row.names(d[ abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='very_stable',])

)

pdf('Boxplots_peakCenterNormalised_interactions_with_K27ac_peaks.pdf',onefile=TRUE,pointsize=12,width=7,height=10)
par(mfrow=c(2,2))
names=c('Bl6','Cast')
col=c('cornflowerblue','orange')
for (i in 1:length(Ld)){
  boxplot(d[row.names(d) %in% Ld[[i]],169:170]/d[row.names(d) %in% Ld[[i]],]$n_dpn_id_within_k27acPeak,outline=FALSE,notch=TRUE,main=paste('reads,',names(Ld)[i], 'peaks'), ylab='reads in peak norm',col=col,names=names)
  boxplot(d[row.names(d) %in% Ld[[i]],171:172]/d[row.names(d) %in% Ld[[i]],]$n_dpn_id_within_k27acPeak,outline=FALSE,notch=TRUE,main=paste('scores,',names(Ld)[i], 'peaks'), ylab='scores in peak norm',col=col,names=names)
  
  boxplot(d[row.names(d) %in% Ld[[i]],174:175],outline=FALSE,notch=TRUE,main=paste('reads,',names(Ld)[i], 'peaks'), ylab='reads in peak center norm',col=col,names=names)
  boxplot(d[row.names(d) %in% Ld[[i]],174:175],outline=FALSE,notch=TRUE,main=paste('scores,',names(Ld)[i], 'peaks'), ylab='scores in peak center norm',col=col,names=names)
}
dev.off()

write.table(d,'ints_k27ac_interacting_annotated.txt',sep='\t',eol='\n',quote=FALSE,row.names=FALSE)

```


##Look at individual genes:

##If they have an interaction with a Bl6-specific peak, but are tx not affected - do they also have an interaction with another peak that increases? With K27ac stable/Cast-specific or with an ATAC peak?

```{r}

###################################################
##Annotation of baited_genes with interactions
###################################################

peaks=read.delim('ints_k27ac_interacting_annotated.txt',stringsAsFactors = FALSE)

names=c('CC12_K27acPeak_interaction',
        'CC12_Bl6_K27acPeak_interaction',
        'CC12_Cast_K27acPeak_interaction',
        'CC12_stable_K27acPeak_interaction',
        
        'CC12_Bl6_K27acPeak_Bl6_interaction',
        'CC12_Cast_K27acPeak_Bl6_interaction',
        'CC12_stable_K27acPeak_Bl6_interaction',
        
        'CC12_Bl6_K27acPeak_Castaneus_interaction',
        'CC12_Cast_K27acPeak_Castaneus_interaction',
        'CC12_stable_K27acPeak_Castaneus_interaction',
        
        'CC12_Bl6_K27acPeak_stable_interaction',
        'CC12_Cast_K27acPeak_stable_interaction',
        'CC12_stable_K27acPeak_stable_interaction'
        
        )
        
        
        
df=matrix(FALSE,ncol=length(names),nrow=length(baited_genes))
df=as.data.frame(df)
names(df)=names



for (i in 1:length(baited_genes))
{
  print(i)
  id=baited_genes[i]$dpn_id
  a=peaks[peaks$baitID==id,]
  
  if (nrow(a)>0)
  {
    df[i,'CC12_K27acPeak_interaction']=TRUE
    
    b=ints[ints$baitID==id,]
    
    ##to assign peak regulation to interacting peaks
    gr=bed2gr(a[,c('Chr_k27acPeak', 'Start_k27acPeak','End_k27acPeak')])
    gr2=bed2gr(b[,c('seqnames_otherEnd','start_otherEnd','end_otherEnd')])
    ov=as.matrix(findOverlaps(gr,gr2))
    a=a[ov[,1],]
    a$clustertype=b[ov[,2],]$clustertype
  
    
    types=c('Bl6','Cast')
    clust=c('Bl6','Castaneus')
    
    for (n in types)
    { 
      
      if (TRUE %in% a[,paste0('k27ac_',n,'_otherEnd')])
      {
        df[i,paste0('CC12_',n,'_K27acPeak_interaction')]=TRUE
        
        l=list()
        
        for (typ in clust)
        {
          names=names(l)
          l=c(l,list( unique(a[ a[,paste0('k27ac_',n,'_otherEnd')] & a[,paste0('reads_within_k27acPeak_normToPeakCenterBl6_',typ)] > a[,paste0('reads_within_k27acPeak_normToPeakCenterBl6_', clust[-grep(typ,clust)] )]   & 
                        a[,paste0('score_within_k27acPeak_normToPeakCenterBl6_',typ)] > a[,paste0('score_within_k27acPeak_normToPeakCenterBl6_',clust[-grep(typ,clust)])]   & 
                          a$clustertype %in% unique(ints$clustertype)[grep(typ,unique(ints$clustertype))],]$name_k27acPeak)))
          
          names(l)=c(names,typ)
          
          if ( length(l[[typ]]) >0 )
          {
            df[i,paste0('CC12_',n,'_K27acPeak_',typ,'_interaction')]=TRUE
          }
          
        }
        
        if ( length(unique(a[a[,paste0('k27ac_',n,'_otherEnd')],]$name_k27acPeak)) > sum(sapply(l,FUN=function(x) length(x))) ) 
        {
          typ='stable'
          df[i,paste0('CC12_',n,'_K27acPeak_',typ,'_interaction')]=TRUE
        }
        
    }
  }
    
    n='stable'
    st=a[ a[,paste0('k27ac_Bl6_otherEnd')] ==FALSE & a[,paste0('k27ac_Cast_otherEnd')] ==FALSE,]
    if ( nrow(st)>0 )
      {
        df[i,paste0('CC12_',n,'_K27acPeak_interaction')]=TRUE
        
        l=list()
        
        for (typ in clust)
        {
          names=names(l)
          l=c(l,list( unique(st[ st[,paste0('reads_within_k27acPeak_normToPeakCenterBl6_',typ)] > st[,paste0('reads_within_k27acPeak_normToPeakCenterBl6_', clust[-grep(typ,clust)] )]   & 
                        st[,paste0('score_within_k27acPeak_normToPeakCenterBl6_',typ)] > st[,paste0('score_within_k27acPeak_normToPeakCenterBl6_',clust[-grep(typ,clust)])]   & 
                          st$clustertype %in% unique(ints$clustertype)[grep(typ,unique(ints$clustertype))],]$name_k27acPeak)))
          
          names(l)=c(names,typ)
          
          if ( length(l[[typ]]) >0 )
          {
            df[i,paste0('CC12_',n,'_K27acPeak_',typ,'_interaction')]=TRUE
          }
          
        }
        
        if ( length(unique(st$name_k27acPeak)) > sum(sapply(l,FUN=function(x) length(x))) ) 
        {
          typ='stable'
          df[i,paste0('CC12_',n,'_K27acPeak_',typ,'_interaction')]=TRUE
        }
        
      }
    
  }
}

summary(df)

##To DO:
# ##Annotate specifically with atac interactions (our FCS data)
# names=c(
#         'CC12_ATACPeak_interaction',
#         'CC12_ATACPeak_Bl6_interaction',
#         'CC12_ATACPeak_Castaneus_interaction',
#         'CC12_ATACPeak_stable_interaction')
#         
#         
#         
# df=matrix(FALSE,ncol=length(names),nrow=length(baited_genes))
# df=as.data.frame(df)
# names(df)=names

# for (i in 1:length(baited_genes))
# {
#   print(i)
#   id=baited_genes[i]$dpn_id
#   a=ints[ints$baitID==id,]
#   if (TRUE %in% a$atac_otherEnd)
#   {
#     df[i,'CC12_ATACPeak_interaction']=TRUE
#     
#     clust=c('Bl6','Cast','stable')
#       
#         for (typ in clust)
#         {
#           if ( nrow(a[ a[,paste0('atac_otherEnd')] & a$clustertype %in% unique(ints$clustertype)[grep(typ,unique(ints$clustertype))],])>0)
#           {
#             df[i,paste0('CC12_ATACPeak_',typ,'_interaction')]=TRUE
#           }
#           
#         }
#       
#     
#   }
# }

values(baited_genes)=cbind(values(baited_genes),df)

#############
##plots:
#############

pdf('Barplots_genes_interaction_with_K27ac_peaks.pdf',onefile=TRUE,height=14,width=10,pointsize=12)

##Which genes are more likely to interact with which type of peaks:
f <- function(cols,rows){
  col=c('cornflowerblue','orange','grey','white')
  tab=table(values(baited_genes)[,rows],values(baited_genes)[,cols])
  tab=cbind(tab,data.frame(table(values(baited_genes)[,rows]))[,2])
  colnames(tab)[ncol(tab)]='all'
  test=empir(as.data.frame(baited_genes),reg=cols,char=rows,n=1000)
  test=do.call('rbind',test)
  test=test[grep('TRUE',row.names(test)),]
  text=apply(test[,c('pvalue_enriched', 'pvalue_depleted')],1,FUN=function(x) min(x)) ##select minimum pvalue
  text=text[match(colnames(tab), sapply(strsplit(names(text),split='_genes'), FUN=function(x) x[1])  )]
  plotNormBarplot_sel(tab,xlab=cols,ylab='Fraction genes',main=rows,sel='TRUE',las=3,col=col,pval=text)
}

par(mfrow=c(2,2))
cols='genetype'
rows='CC12_K27acPeak_interaction'
f(cols,rows)
rows='CC12_stable_K27acPeak_interaction'
f(cols,rows)
rows='CC12_Bl6_K27acPeak_interaction'
f(cols,rows)
rows='CC12_Cast_K27acPeak_interaction'
f(cols,rows)

##Conclusion: Interactions are ususally with cis-regulated peaks: up genes interact with upregulated peaks etc.
##Question: Are genes more likely to have cis-peaks in their vicinity? 

##A)Make a simple calculation by looking at K27ac peaks within 1Mb from the gene promoter

b=resize(baited_genes,2000000,fix='center')
baited_genes$n_k27acPeaks_plmin1Mb=countOverlaps(b,k27ac)
baited_genes$n_Bl6_k27acPeaks_plmin1Mb=countOverlaps(b,k27ac_Bl6)
baited_genes$n_Castaneus_k27acPeaks_plmin1Mb=countOverlaps(b,k27ac_Cast)

##substract K27ac overlapping TSS
baited_genes[baited_genes$k27ac_bait]$n_k27acPeaks_plmin1Mb= baited_genes[baited_genes$k27ac_bait]$n_k27acPeaks_plmin1Mb-1
baited_genes[baited_genes$k27ac_loss_bait]$n_Bl6_k27acPeaks_plmin1Mb= baited_genes[baited_genes$k27ac_loss_bait]$n_Bl6_k27acPeaks_plmin1Mb-1
baited_genes[baited_genes$k27ac_gain_bait]$n_Castaneus_k27acPeaks_plmin1Mb= baited_genes[baited_genes$k27ac_gain_bait]$n_Castaneus_k27acPeaks_plmin1Mb-1
baited_genes$n_stable_k27acPeaks_plmin1Mb=baited_genes$n_k27acPeaks_plmin1Mb-baited_genes$n_Bl6_k27acPeaks_plmin1Mb-baited_genes$n_Castaneus_k27acPeaks_plmin1Mb

dd=as.data.frame(values(baited_genes))
col=c('cornflowerblue','orange','grey')
by='genetype'
what='n_k27acPeaks_plmin1Mb'
boxplot(dd[,what]~dd[,by],outline=FALSE,las=3,ylab=what,xlab=by,notch=TRUE,col=col)
what='n_stable_k27acPeaks_plmin1Mb'
boxplot(dd[,what]~dd[,by],outline=FALSE,las=3,ylab=what,xlab=by,notch=TRUE,col=col)
what='n_Bl6_k27acPeaks_plmin1Mb'
boxplot(dd[,what]~dd[,by],outline=FALSE,las=3,ylab=what,xlab=by,notch=TRUE,col=col)
what='n_Castaneus_k27acPeaks_plmin1Mb'
boxplot(dd[,what]~dd[,by],outline=FALSE,las=3,ylab=what,xlab=by,notch=TRUE,col=col)


##B)Make a simple calculation by looking at K27ac peaks within 500kb from the gene promoter

b=resize(baited_genes,1000000,fix='center')
baited_genes$n_k27acPeaks_plmin0.5Mb=countOverlaps(b,k27ac)
baited_genes$n_Bl6_k27acPeaks_plmin0.5Mb=countOverlaps(b,k27ac_Bl6)
baited_genes$n_Castaneus_k27acPeaks_plmin0.5Mb=countOverlaps(b,k27ac_Cast)

##substract K27ac overlapping TSS
baited_genes[baited_genes$k27ac_bait]$n_k27acPeaks_plmin0.5Mb= baited_genes[baited_genes$k27ac_bait]$n_k27acPeaks_plmin0.5Mb-1
baited_genes[baited_genes$k27ac_loss_bait]$n_Bl6_k27acPeaks_plmin0.5Mb= baited_genes[baited_genes$k27ac_loss_bait]$n_Bl6_k27acPeaks_plmin0.5Mb-1
baited_genes[baited_genes$k27ac_gain_bait]$n_Castaneus_k27acPeaks_plmin0.5Mb= baited_genes[baited_genes$k27ac_gain_bait]$n_Castaneus_k27acPeaks_plmin0.5Mb-1
baited_genes$n_stable_k27acPeaks_plmin0.5Mb=baited_genes$n_k27acPeaks_plmin0.5Mb-baited_genes$n_Bl6_k27acPeaks_plmin0.5Mb-baited_genes$n_Castaneus_k27acPeaks_plmin0.5Mb

dd=as.data.frame(values(baited_genes))
col=c('cornflowerblue','orange','grey')
by='genetype'
what='n_k27acPeaks_plmin0.5Mb'
boxplot(dd[,what]~dd[,by],outline=FALSE,las=3,ylab=what,xlab=by,notch=TRUE,col=col)
what='n_stable_k27acPeaks_plmin0.5Mb'
boxplot(dd[,what]~dd[,by],outline=FALSE,las=3,ylab=what,xlab=by,notch=TRUE,col=col)
what='n_Bl6_k27acPeaks_plmin0.5Mb'
boxplot(dd[,what]~dd[,by],outline=FALSE,las=3,ylab=what,xlab=by,notch=TRUE,col=col)
what='n_Castaneus_k27acPeaks_plmin0.5Mb'
boxplot(dd[,what]~dd[,by],outline=FALSE,las=3,ylab=what,xlab=by,notch=TRUE,col=col)

##C)Make a simple calculation by looking at K27ac peaks within 250kb from the gene promoter

b=resize(baited_genes,500000,fix='center')
baited_genes$n_k27acPeaks_plmin0.25Mb=countOverlaps(b,k27ac)
baited_genes$n_Bl6_k27acPeaks_plmin0.25Mb=countOverlaps(b,k27ac_Bl6)
baited_genes$n_Castaneus_k27acPeaks_plmin0.25Mb=countOverlaps(b,k27ac_Cast)

##substract K27ac overlapping TSS
baited_genes[baited_genes$k27ac_bait]$n_k27acPeaks_plmin0.25Mb= baited_genes[baited_genes$k27ac_bait]$n_k27acPeaks_plmin0.25Mb-1
baited_genes[baited_genes$k27ac_loss_bait]$n_Bl6_k27acPeaks_plmin0.25Mb= baited_genes[baited_genes$k27ac_loss_bait]$n_Bl6_k27acPeaks_plmin0.25Mb-1
baited_genes[baited_genes$k27ac_gain_bait]$n_Castaneus_k27acPeaks_plmin0.25Mb= baited_genes[baited_genes$k27ac_gain_bait]$n_Castaneus_k27acPeaks_plmin0.25Mb-1
baited_genes$n_stable_k27acPeaks_plmin0.25Mb=baited_genes$n_k27acPeaks_plmin0.25Mb-baited_genes$n_Bl6_k27acPeaks_plmin0.25Mb-baited_genes$n_Castaneus_k27acPeaks_plmin0.25Mb

dd=as.data.frame(values(baited_genes))
col=c('cornflowerblue','orange','grey')
by='genetype'
what='n_k27acPeaks_plmin0.25Mb'
boxplot(dd[,what]~dd[,by],outline=FALSE,las=3,ylab=what,xlab=by,notch=TRUE,col=col)
what='n_stable_k27acPeaks_plmin0.25Mb'
boxplot(dd[,what]~dd[,by],outline=FALSE,las=3,ylab=what,xlab=by,notch=TRUE,col=col)
what='n_Bl6_k27acPeaks_plmin0.25Mb'
boxplot(dd[,what]~dd[,by],outline=FALSE,las=3,ylab=what,xlab=by,notch=TRUE,col=col)
what='n_Castaneus_k27acPeaks_plmin0.25Mb'
boxplot(dd[,what]~dd[,by],outline=FALSE,las=3,ylab=what,xlab=by,notch=TRUE,col=col)



#################################################
##Include how the peaks that genes interact with
##are regulated
#################################################

par(mfrow=c(4,3))

##plotting function:
f <- function(cols,rows,ylim=c(0,1)){
  col=c('blue','green','grey')
  values(baited_genes)[,rows]=factor(values(baited_genes)[,rows],levels=c(FALSE,TRUE))
  values(s)[,rows]=factor(values(s)[,rows],levels=c(FALSE,TRUE))
  values(u)[,rows]=factor(values(u)[,rows],levels=c(FALSE,TRUE))
  values(d)[,rows]=factor(values(d)[,rows],levels=c(FALSE,TRUE))
  
  tab=table(values(s)[,rows],values(s)[,cols])
  tab=cbind(tab,data.frame(table(values(baited_genes)[,rows]))[,2])
  colnames(tab)[ncol(tab)]='all'
  test=empir(as.data.frame(baited_genes),reg=cols,char=rows,ntest=1000)
  test=do.call('rbind',test)
  test=test[grep('TRUE$',row.names(test)),]
  text=apply(test[,c('pvalue_enriched', 'pvalue_depleted')],1,FUN=function(x) min(x)) ##select minimum pvalue
  text=text[match( colnames(tab),sapply(strsplit(names(text),split='_genes'), FUN=function(x) x[1]) )]
  plotNormBarplot_sel(tab,xlab=cols,ylab=paste('genes with',rows),main='stable genes',col=col,sel='TRUE',ylim=ylim,pval=text)
  
  tab=table(values(d)[,rows],values(d)[,cols])
  tab=cbind(tab,data.frame(table(values(baited_genes)[,rows]))[,2])
  colnames(tab)[ncol(tab)]='all'
  test=empir(as.data.frame(baited_genes),reg=cols,char=rows,ntest=1000)
  test=do.call('rbind',test)
  test=test[grep('TRUE$',row.names(test)),]
  text=apply(test[,c('pvalue_enriched', 'pvalue_depleted')],1,FUN=function(x) min(x)) ##select minimum pvalue
  text=text[match( colnames(tab),sapply(strsplit(names(text),split='_genes'), FUN=function(x) x[1]) )]
  plotNormBarplot_sel(tab,xlab=cols,ylab=paste('genes with',rows),main='Bl6 up genes',col=col,sel='TRUE',ylim=ylim,pval=text)
  
  tab=table(values(u)[,rows],values(u)[,cols])
  tab=cbind(tab,data.frame(table(values(baited_genes)[,rows]))[,2])
  colnames(tab)[ncol(tab)]='all'
  test=empir(as.data.frame(baited_genes),reg=cols,char=rows,ntest=1000)
  test=do.call('rbind',test)
  test=test[grep('TRUE$',row.names(test)),]
  text=apply(test[,c('pvalue_enriched', 'pvalue_depleted')],1,FUN=function(x) min(x)) ##select minimum pvalue
  text=text[match( colnames(tab),sapply(strsplit(names(text),split='_genes'), FUN=function(x) x[1]) )]
  plotNormBarplot_sel(tab,xlab=cols,ylab=paste('genes with',rows),main='Cast up genes',col=col,sel='TRUE',ylim=ylim,pval=text)
  
}

s=baited_genes[baited_genes$genetype=='very_stable']
d=baited_genes[baited_genes$genetype=='down_strong']
u=baited_genes[baited_genes$genetype=='up_strong']

##1a) Interaction with stable peaks for genes that interact with Bl6-specific K27ac peaks:
par(mfrow=c(4,3))
cols='CC12_Bl6_K27acPeak_interaction'

rows='CC12_stable_K27acPeak_interaction'
f(cols,rows)
rows='CC12_stable_K27acPeak_stable_interaction'
f(cols,rows)
rows='CC12_stable_K27acPeak_Castaneus_interaction'
f(cols,rows)
rows='CC12_stable_K27acPeak_Bl6_interaction'
f(cols,rows)

##1b) Interaction with Bl6 peaks for genes that interact with Bl6-specific K27ac peaks:
par(mfrow=c(4,3))
cols='CC12_Bl6_K27acPeak_interaction'

rows='CC12_Bl6_K27acPeak_interaction'
f(cols,rows)
rows='CC12_Bl6_K27acPeak_stable_interaction'
f(cols,rows)
rows='CC12_Bl6_K27acPeak_Castaneus_interaction'
f(cols,rows)
rows='CC12_Bl6_K27acPeak_Bl6_interaction'
f(cols,rows)

##1c) Interaction with Cast peaks for genes that interact with Bl6-specific K27ac peaks:
par(mfrow=c(4,3))
cols='CC12_Bl6_K27acPeak_interaction'

rows='CC12_Cast_K27acPeak_interaction'
f(cols,rows)
rows='CC12_Cast_K27acPeak_stable_interaction'
f(cols,rows)
rows='CC12_Cast_K27acPeak_Castaneus_interaction'
f(cols,rows)
rows='CC12_Cast_K27acPeak_Bl6_interaction'
f(cols,rows)

##2a) Interaction with stable peaks for genes that interact with Cast-specific K27ac peaks:
par(mfrow=c(4,3))
cols='CC12_Cast_K27acPeak_interaction'

rows='CC12_stable_K27acPeak_interaction'
f(cols,rows)
rows='CC12_stable_K27acPeak_stable_interaction'
f(cols,rows)
rows='CC12_stable_K27acPeak_Castaneus_interaction'
f(cols,rows)
rows='CC12_stable_K27acPeak_Bl6_interaction'
f(cols,rows)

##2b) Interaction with Bl6 peaks for genes that interact with Cast-specific K27ac peaks:
par(mfrow=c(4,3))
cols='CC12_Cast_K27acPeak_interaction'

rows='CC12_Bl6_K27acPeak_interaction'
f(cols,rows)
rows='CC12_Bl6_K27acPeak_stable_interaction'
f(cols,rows)
rows='CC12_Bl6_K27acPeak_Castaneus_interaction'
f(cols,rows)
rows='CC12_Bl6_K27acPeak_Bl6_interaction'
f(cols,rows)

##2c) Interaction with Cast peaks for genes that interact with Cast-specific K27ac peaks:
par(mfrow=c(4,3))
cols='CC12_Cast_K27acPeak_interaction'

rows='CC12_Cast_K27acPeak_interaction'
f(cols,rows)
rows='CC12_Cast_K27acPeak_stable_interaction'
f(cols,rows)
rows='CC12_Cast_K27acPeak_Castaneus_interaction'
f(cols,rows)
rows='CC12_Cast_K27acPeak_Bl6_interaction'
f(cols,rows)


####################################################
##Boxplots with normalized read and score counts for
##genes with different types of peak interactions
####################################################

d=peaks
row.names(d)=d$intID

##downregulated, upregulated and stable peaks:
Ld=list(
  stable.genes_withBl6Peaks_stablePeak=row.names(d[d$baitID %in%  baited_genes[baited_genes$CC12_Bl6_K27acPeak_interaction]$dpn_id & abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='very_stable',]),
  stable.genes_withoutBl6Peaks_stablePeak=row.names(d[d$baitID %in%  baited_genes[baited_genes$CC12_Bl6_K27acPeak_interaction==FALSE]$dpn_id & abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='very_stable',]),
  
  Bl6up.genes_withBl6Peaks_stablePeak=row.names(d[d$baitID %in%  baited_genes[baited_genes$CC12_Bl6_K27acPeak_interaction]$dpn_id & abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='down_strong',]),
  Bl6up.genes_withoutBl6Peaks_stablePeak=row.names(d[d$baitID %in%  baited_genes[baited_genes$CC12_Bl6_K27acPeak_interaction==FALSE]$dpn_id & abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='down_strong',]),
  
  Castup.genes_withBl6Peaks_stablePeak=row.names(d[d$baitID %in%  baited_genes[baited_genes$CC12_Bl6_K27acPeak_interaction]$dpn_id & abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='up_strong',]),
  Castup.genes_withoutBl6Peaks_stablePeak=row.names(d[d$baitID %in%  baited_genes[baited_genes$CC12_Bl6_K27acPeak_interaction==FALSE]$dpn_id & abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='up_strong',]),
  
  
  stable.genes_withCastPeaks_stablePeak=row.names(d[d$baitID %in%  baited_genes[baited_genes$CC12_Cast_K27acPeak_interaction]$dpn_id & abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='very_stable',]),
  stable.genes_withoutCastPeaks_stablePeak=row.names(d[d$baitID %in%  baited_genes[baited_genes$CC12_Cast_K27acPeak_interaction==FALSE]$dpn_id & abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='very_stable',]),
  
  Bl6up.genes_withCastPeaks_stablePeak=row.names(d[d$baitID %in%  baited_genes[baited_genes$CC12_Cast_K27acPeak_interaction]$dpn_id & abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='down_strong',]),
  Bl6up.genes_withoutCastPeaks_stablePeak=row.names(d[d$baitID %in%  baited_genes[baited_genes$CC12_Casr_K27acPeak_interaction==FALSE]$dpn_id & abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='down_strong',]),
  
  Castup.genes_withCastPeaks_stablePeak=row.names(d[d$baitID %in%  baited_genes[baited_genes$CC12_Cast_K27acPeak_interaction]$dpn_id & abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='up_strong',]),
  Castup.genes_withoutCastPeaks_stablePeak=row.names(d[d$baitID %in%  baited_genes[baited_genes$CC12_Cast_K27acPeak_interaction==FALSE]$dpn_id & abs(d$LFC_k27acPeak) < 1 & d$padj_k27acPeak >= 0.05 & is.na(d$padj_k27acPeak)==FALSE & d$genetype=='up_strong',]),
  
)

par(mfrow=c(2,2))
names=c('Bl6','Cast')
col=c('cornflowerblue','orange')
for (i in 1:length(Ld)){
  if (length(Ld[[i]])>0){
  boxplot(d[row.names(d) %in% Ld[[i]],169:170]/d[row.names(d) %in% Ld[[i]],]$n_dpn_id_within_k27acPeak,outline=FALSE,notch=TRUE,main=paste('reads,',names(Ld)[i]), ylab='reads in stable peak norm',col=col,names=names)
  boxplot(d[row.names(d) %in% Ld[[i]],171:172]/d[row.names(d) %in% Ld[[i]],]$n_dpn_id_within_k27acPeak,outline=FALSE,notch=TRUE,main=paste('scores,',names(Ld)[i]), ylab='scores in stable peak norm',col=col,names=names)
  
  boxplot(d[row.names(d) %in% Ld[[i]],174:175],outline=FALSE,notch=TRUE,main=paste('reads,',names(Ld)[i]), ylab='reads in stable peak center norm',col=col,names=names)
  boxplot(d[row.names(d) %in% Ld[[i]],174:175],outline=FALSE,notch=TRUE,main=paste('scores,',names(Ld)[i]), ylab='scores in stable peak center norm',col=col,names=names)
  }
}



dev.off()


# ##Just once:
# saveRDS(baited_genes,'../designDir/baited_genes_annotated.rds')

```


##plot interactions of genes to have a look at whether all peak interactions are being captured:

```{r}

peaks=read.delim('ints_k27ac_interacting_annotated.txt',stringsAsFactors = FALSE)

intervals=GRangesList(stripGR(tad),stripGR(k27ac),k27ac_Bl6, k27ac_Cast,stripGR(refseqtss))
names(intervals)=c('TADs_Bonev','K27ac_2i_all','k27ac_Bl6','k27ac_Cast','RefseqTSS')


##define different genes:

l=list(
  stableGenes=baited_genes[baited_genes$genetype=='very_stable'],
  Bl6UpGenes=baited_genes[baited_genes$genetype=='down_strong'],
  CastaneusUpGenes=baited_genes[baited_genes$genetype=='up_strong']
)

for ( cl in 1:length(l)){
  print(cl)
  
  sam=sample(l[[cl]]$dpn_id,20)

  pdf(paste0('Lineplots_sample_',names(l)[cl],'2.pdf',sep=''), width=20, height=14,onefile = TRUE,pointsize=12)
   par(mfrow=c(2,2))
   
  for (id in sam){
      
      loc=c("Chr_k27acPeak" , "Start_k27acPeak" ,"End_k27acPeak")
      
      a=peaks[peaks$baitID==id ,]
      b6=bed2gr(a[ grep('Bl6',a$clustertype), loc])
      ca=bed2gr(a[ grep('Cast',a$clustertype), loc])
      st=bed2gr(a[ grep('stable',a$clustertype), loc])
      iv=c(intervals, GRangesList(st),GRangesList(b6),GRangesList(ca))
      names(iv)=c(names(intervals),'stable_k27acInts','Bl6_up_k27acInts','Cast_up_k27acInts')
      
      print(id)
      name=paste(baited_genes[baited_genes$dpn_id==id,]$genenames,baited_genes[baited_genes$dpn_id==id,]$inttype,baited_genes[baited_genes$dpn_id==id,]$genetype,sep=';')
      
      colintervals=c('lightgrey','purple','cornflowerblue','orange', 'red', 'blue', 'green', 'pink', 'black')
      
      col=c('cornflowerblue','orange')
       
      k=10
      zoom=200000
      ylim=c(0,100)
      plotInteractions(L,id,k,zoom,unt=1,tam=2, ylim=ylim,show.legend = TRUE,name=name,intervals=iv,colintervals = colintervals, col=col)
      
      ylim=c(0,20)
      zoom=1000000
      k=39
      plotInteractions(L,id,k,zoom,unt=1,tam=2, ylim=ylim,show.legend = TRUE,name=name,intervals=iv,colintervals = colintervals, col=col)
      
      
    }
    dev.off()

}


```

##TO DO:
##Quantify interactions with all K27ac peaks within a certain distance from the gene that have >=3 reads and a score of >=4 at at least one position



##Bigwig URLs for all samples:

```{r}

# setwd('~/')
dir.create('bed')
dir.create('bed/norm')
setwd('./bed/norm')

# for (id in baited_genes$dpn_id){
for (id in baited_genes$dpn_id){
  print(id)
  getBedGraphsNormGetbednorm(id,L,hind)
  getBedGraphsChicagoScore(id,L,hind)
}


##later: rename all 1:189 bg so that instead of the refseq ID they are called by the genename!!

###################################
##now make UCSC compatible tracks:
###################################

cd /omics/groups/OE0624/internal/angelika/analysis/CaptureC/MyData/CapC16/bed/norm

mkdir UCSC_compatible
chromSizes=/omics/groups/OE0624/internal/databank/chrom.sizes/mm10.chrom.sizes.txt

# module load bedtools/2.29.2
# module load ensembl_tools/76
module load ucsc.genome-userapps/v392

for f in *bg.bed
do 
sort -k1,1 -k2,2 -n "$f" > "$f"_sorted.bed
bedGraphToBigWig "$f"_sorted.bed ${chromSizes}  UCSC_compatible/"$f"_bg_sorted.bw
done

# bedGraphToBigWig 9230102O04Rik_id3011673_Bl6_bg.bed_sorted.bed ${chromSizes}  UCSC_compatible/"9230102O04Rik_id3011673_Bl6_bg.bed_sorted.bed_bg_sorted.bw"

for f in *viewpoint.bed
do 
sort -k1,1 -k2,2 -n "$f" > "$f"_sorted.bed
bedToBigBed "$f"_sorted.bed ${chromSizes}  UCSC_compatible/"$f"_bg_sorted.bb
done


########################################################################################################################
##transfer the files to the galaxy server and write a hub.txt based on their location on the server, but programmed in R  
#########################################################################################################################


ids=baited_genes$dpn_id

cols=c('cornflowerblue','orange')
names(cols)=c('Bl6','Cast') 
df=data.frame(cols=cols,
              rgb=c('67,162,202', '217,95,14') 
              )
row.names(df)=names(cols)

######################
##SuperTrack
######################

l=list()

for (id in ids){

  name=paste(names(baited_genes[baited_genes$dpn_id==id]),id,sep='_id')
  shortname=names(baited_genes[baited_genes$dpn_id==id])
  
  v=c(paste('track ', shortname),
      'superTrack on show',
      'group regulation',
      paste('shortLabel', shortname),
      paste('longLabel',shortname),
      'windowingFunction mean',
      'allButtonPair on',
      'Configurable on',
      'autoScale on'
      )
  
  
  l=c(l,list(v))
  
  a=c(
      paste('track composite_',name,sep=''),
      paste('parent ',shortname,sep=''),
      'type bigWig',
      'compositeTrack on show'
      )
  
  l=c(l,list(a))
  
  for (sample in row.names(df)){
    
    col=cols[names(cols)==sample]
    ru=as.character(df[df$cols==col,'rgb'])
    
    a=c(
      paste('track composite_',name,'_',sample,sep=''),
      paste('bigDataUrl http://nasmythserver1.bioch.ox.ac.uk/UCSC/angelika.feldmann@bioch.ox.ac.uk/CapC/CapC16/bigwig/',name,'_',sample,'_bg.bed_bg_sorted.bw',sep=''),
      paste('shortLabel ',sample,sep=''),
      paste('longLabel ',shortname,'_',sample,sep=''),
      'type bigWig',
      'visibility full',
      paste('parent composite_',name, ' on',sep=''),
      paste('color',ru ))
    
    l=c(l,list(a))
  
    a=c(
      paste('track composite_',name,'_',sample,'_scores',sep=''),
      paste('bigDataUrl http://nasmythserver1.bioch.ox.ac.uk/UCSC/angelika.feldmann@bioch.ox.ac.uk/CapC/CapC16/bigwig/',name,'_',sample,'_scores_bg.bed_bg_sorted.bw',sep=''),
      paste('shortLabel ',sample,'_scores',sep=''),
      paste('longLabel ',shortname,'_',sample,'_scores',sep=''),
      'type bigWig',
      'visibility full',
      paste('parent composite_',name, ' on',sep=''),
      paste('color',ru ))
    
    l=c(l,list(a))
  
  }
    
  
  a=c(paste('track multiWig_reads_',name,sep=''),
      paste('parent ',shortname, ' on',sep=''),
    'type bigWig',
    'container multiWig',
    'aggregate transparentOverlay',
    'showSubtrackColorOnUi on',
    'maxHeightPixels 500:100:8'
  )
  
  
  l=c(l,list(a))
  
  for (sample in row.names(df)){
    
    col=cols[names(cols)==sample]
    ru=as.character(df[df$cols==col,'rgb'])
    
    a=c(
      paste('track multiWig_reads_',name,'_',sample,sep=''),
      paste('bigDataUrl http://nasmythserver1.bioch.ox.ac.uk/UCSC/angelika.feldmann@bioch.ox.ac.uk/CapC/CapC16/bigwig/',name,'_',sample,'_bg.bed_bg_sorted.bw',sep=''),
      paste('shortLabel reads_',sample,sep=''),
      paste('longLabel reads_',shortname,'_',sample,sep=''),
      'type bigWig',
      'visibility full',
      paste('parent multiWig_reads_',name, ' on',sep=''),
      paste('color',ru ))
    
    l=c(l,list(a))
  
  }
  
  a=c(paste('track multiWig_scores_',name,sep=''),
      paste('parent ',shortname, ' on',sep=''),
    'type bigWig',
    'container multiWig',
    'aggregate transparentOverlay',
    'showSubtrackColorOnUi on',
    'maxHeightPixels 500:100:8'
  )
  
  
  l=c(l,list(a))
  
  for (sample in row.names(df)){
    
    col=cols[names(cols)==sample]
    ru=as.character(df[df$cols==col,'rgb'])
    
    a=c(
      paste('track multiWig_scores_',name,'_',sample,sep=''),
      paste('bigDataUrl http://nasmythserver1.bioch.ox.ac.uk/UCSC/angelika.feldmann@bioch.ox.ac.uk/CapC/CapC16/bigwig/',name,'_',sample,'_scores_bg.bed_bg_sorted.bw',sep=''),
      paste('shortLabel scores_',sample,sep=''),
      paste('longLabel scores_',shortname,'_',sample,sep=''),
      'type bigWig',
      'visibility full',
      paste('parent multiWig_scores_',name, ' on',sep=''),
      paste('color',ru ))
    
    l=c(l,list(a))
  
  }
  
  
  a=c(
      paste('track viewpoint_',name,sep=''),
      paste('bigDataUrl http://nasmythserver1.bioch.ox.ac.uk/UCSC/angelika.feldmann@bioch.ox.ac.uk/CapC/CapC16/bigwig/',name,'_viewpoint.bed_bg_sorted.bb',sep=''),
      paste('shortLabel viewpoint_',shortname,sep=''),
      paste('longLabel viewpoint_',name,sep=''),
      'type bigBed',
      'visibility full',
      paste('parent ', shortname, ' on',sep=''),
      paste('color 0,200,0'))
    

  l=c(l,list(a))
  
  
  
}


# dir.create('mm10')
lapply(l, function(x) cat(c(x, " "), file = "mm10/tracks.txt", sep = "\n", append = TRUE))

##transfer the tracks.txt file to the server!  




######################
##SuperTrack by chromosome
######################

l=list()
path='bigDataUrl http://nasmythserver1.bioch.ox.ac.uk/UCSC/angelika.feldmann@bioch.ox.ac.uk/CapC/CapC12/bigwig/'
limitsreads='0:100'
limitsscores='0:10'


for (ch in unique(seqnames(baited_genes))){
  
  print(ch)
  
  
  v=c(paste('track ', ch),
      'superTrack on show',
      'group regulation',
      paste('shortLabel', ch),
      paste('longLabel',ch),
      'windowingFunction mean',
      'allButtonPair on',
      'Configurable on',
      'autoScale off',
      'maxHeightPixels 0:50:128'
      )
  
  
  l=c(l,list(v))
  
  ids=baited_genes[seqnames(baited_genes)==ch]$dpn_id
  
  for (id in ids){
  
  name=paste(names(baited_genes[baited_genes$dpn_id==id]),id,sep='_id')
  shortname=names(baited_genes[baited_genes$dpn_id==id])
  
  a=c(
      paste('track composite_',name,sep=''),
      paste('parent ',ch,sep=''),
      'type bigWig',
      'compositeTrack on show'
      )
  
  l=c(l,list(a))
  
  for (sample in row.names(df)){
    
    col=cols[names(cols)==sample]
    ru=as.character(df[df$cols==col,'rgb'])
    
    a=c(
      paste('track composite_',name,'_',sample,sep=''),
      paste(path,name,'_',sample,'_bg.bed_bg_sorted.bw',sep=''),
      paste('shortLabel ',sample,sep=''),
      paste('longLabel ',shortname,'_',sample,sep=''),
      'type bigWig',
      'visibility full',
      paste('parent composite_',name, ' on',sep=''),
      paste('color',ru ),
      paste('viewLimits', limitsreads ))
    
    l=c(l,list(a))
  
    a=c(
      paste('track composite_',name,'_',sample,'_scores',sep=''),
      paste(path,name,'_',sample,'_scores_bg.bed_bg_sorted.bw',sep=''),
      paste('shortLabel ',sample,'_scores',sep=''),
      paste('longLabel ',shortname,'_',sample,'_scores',sep=''),
      'type bigWig',
      'visibility full',
      paste('parent composite_',name, ' on',sep=''),
      paste('color',ru ),
      paste('viewLimits', limitsscores ))
    
    l=c(l,list(a))
  
  }
    
  
  a=c(paste('track multiWig_reads_',name,sep=''),
      paste('parent ',ch, ' on',sep=''),
    'type bigWig',
    'container multiWig',
    'aggregate transparentOverlay',
    'showSubtrackColorOnUi on',
    paste('viewLimits', limitsreads )
  )
  
  
  l=c(l,list(a))
  
  for (sample in row.names(df)){
    
    col=cols[names(cols)==sample]
    ru=as.character(df[df$cols==col,'rgb'])
    
    a=c(
      paste('track multiWig_reads_',name,'_',sample,sep=''),
      paste(path,name,'_',sample,'_bg.bed_bg_sorted.bw',sep=''),
      paste('shortLabel reads_',sample,sep=''),
      paste('longLabel reads_',shortname,'_',sample,sep=''),
      'type bigWig',
      'visibility full',
      paste('parent multiWig_reads_',name, ' on',sep=''),
      paste('color',ru ))
    
    l=c(l,list(a))
  
  }
  
  a=c(paste('track multiWig_scores_',name,sep=''),
      paste('parent ',ch, ' on',sep=''),
    'type bigWig',
    'container multiWig',
    'aggregate transparentOverlay',
    'showSubtrackColorOnUi on',
    paste('viewLimits', limitsscores ))
    
  
  
  l=c(l,list(a))
  
  for (sample in row.names(df)){
    
    col=cols[names(cols)==sample]
    ru=as.character(df[df$cols==col,'rgb'])
    
    a=c(
      paste('track multiWig_scores_',name,'_',sample,sep=''),
      paste(path,name,'_',sample,'_scores_bg.bed_bg_sorted.bw',sep=''),
      paste('shortLabel scores_',sample,sep=''),
      paste('longLabel scores_',shortname,'_',sample,sep=''),
      'type bigWig',
      'visibility full',
      paste('parent multiWig_scores_',name, ' on',sep=''),
      paste('color',ru ))
    
    l=c(l,list(a))
  
  }
  
  
  a=c(
      paste('track viewpoint_',name,sep=''),
      paste(path,name,'_viewpoint.bed_bg_sorted.bb',sep=''),
      paste('shortLabel viewpoint_',shortname,sep=''),
      paste('longLabel viewpoint_',name,sep=''),
      'type bigBed',
      'visibility full',
      paste('parent ', ch, ' on',sep=''),
      paste('color 0,200,0'))
    

  l=c(l,list(a))
  
  
  
  }
}


# dir.create('mm10_byChr')
lapply(l, function(x) cat(c(x, " "), file = "mm10_byChr/tracks.txt", sep = "\n", append = TRUE))

##transfer the tracks.txt file to the server!  



##Bedfile with Nanog interactions:
saveBed(bed2gr(ints[ints$genename=='Nanog',c('seqnames_otherEnd','start_otherEnd','end_otherEnd')]),'Nanog_interactome.bed')

```










```{r}

dis=2
peaks=read.delim(paste('ints_all_aggregated_distance',dis,'_summits_annotated.txt',sep=''),stringsAsFactors=FALSE)

##make data frame d with clusters_refined and other end data encompassing the peak summit!
d=peaks
d$otherEndID=d$summit
d[,c('seqnames_otherEnd','start_otherEnd','end_otherEnd')]=d[,c('seqnames_summit','start_summit','end_summit')]
dim(d)
d=d[d$seqnames_bait==d$seqnames_otherEnd,]
dim(d)


##load matrices:

Lmed=loadMatrices(L[1:2],1)
Lpc=loadMatrices(L[3:4],1)
untm=Lmed[[1]]
tamm=Lmed[[2]]
untmc=Lmed[[3]]
tammc=Lmed[[4]]

untp=Lpc[[1]]
tamp=Lpc[[2]]
untpc=Lpc[[3]]
tampc=Lpc[[4]]

##remove rows with all 0 values in UNT:


pdf('Metaprofiles_oE_by_IntervalTypes_normTo_UNT.pdf')
           
##make a List (LL) with different types of intervals:
LL=list(
  Ring1b_peaks_oE=row.names(d[d$Ring1b_peaks_oE,]),
  Cdk8_peaks_oE=row.names(d[d$Cdk8_peaks_oE,]),
  Ring1b.Cdk8_peaks_oE=row.names(d[d$Ring1b_peaks_oE ==TRUE & d$Cdk8_peaks_oE ==TRUE ,]),
  noRing1b.Cdk8_peaks_oE=row.names(d[d$Ring1b_peaks_oE ==FALSE & d$Cdk8_peaks_oE,]),
  Ring1b.Pcgf2_peaks_oE=row.names(d[d$Pcgf2_peaks_oE & d$Ring1b_peaks_oE,]),
  Ring1b.Pcgf2.Suz12_peaks_oE=row.names(d[d$Pcgf2_peaks_oE & d$Ring1b_peaks_oE & d$Suz12_peaks_oE,]),
  Pcgf2_peaks_oE=row.names(d[d$Pcgf2_peaks_oE,]),
  Cdk8.Pcgf2_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==TRUE & d$Cdk8_peaks_oE ==TRUE ,]),
  
  Cdk8.noPcgf2_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==FALSE & d$Cdk8_peaks_oE,]),
  
  Atac_peaks_oE=row.names(d[d$atac_oE,]),
  
  Cdk8.Pcgf2.Ring1b_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==TRUE & d$Cdk8_peaks_oE ==TRUE & d$Ring1b_peaks_oE,]),
  
  Cdk8.Pcgf2.Ring1b.Suz12_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==TRUE & d$Cdk8_peaks_oE ==TRUE & d$Ring1b_peaks_oE & d$Suz12_peaks_oE,]),
  
  Cdk8.noPcgf2.noRing1b.noSuz12_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==FALSE & d$Cdk8_peaks_oE ==TRUE & d$Ring1b_peaks_oE ==FALSE & d$Suz12_peaks_oE ==FALSE ,]),
  
  ATAC.noPcgf2.noRing1b.noSuz12_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==FALSE & d$Suz12_peaks_oE ==FALSE & d$Ring1b_peaks_oE==FALSE & d$atac_otherEnd,]),
  
  ATAC.noPcgf2.noRing1b.noSuz12.noCdk8_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==FALSE & d$atac_oE ==TRUE & d$Ring1b_peaks_oE==FALSE & d$Cdk8_peaks_oE==FALSE & d$Suz12_peaks_oE ==FALSE,]),
  
  TSS_oE=row.names(d[d$RefseqTSS_oE,]),
  
  ATAC.noTSS_oE=row.names(d[d$RefseqTSS_oE==FALSE & d$atac_oE,])
)

u=3
t=4

for (i in 1:length(LL)){

  par(mfrow=c(2,2))
  xlab='distanceFromSummit [DpnII frags]'
  ylab='enrichment vs UNT summit'
  ylim=c(0,1.1)
  
  c=LL[[i]]
  
  main=names(LL)[i]
  n=length(c)
  
  plot(1:ncol(untm),base::colMedians(untm[row.names(untm)%in% c,],na.rm=TRUE),type='l',col='blue',xlab=xlab,ylab=ylab,ylim=ylim,main=main)
  lines(1:ncol(untm),base::colMeans(tamm[row.names(tamm)%in% c,],na.rm=TRUE),type='l',col='red')
  legend('topleft',col=c('white','blue','red'),lty=1,cex=0.7,legend=c(paste0('n=',n),'Med13fl_UNT','Med13fl_TAM'),bty='n')
  
  plot(1:ncol(untp),base::colMeans(untp[row.names(untp)%in% c,],na.rm=TRUE),type='l',col='black',xlab=xlab,ylab=ylab,ylim=ylim,main=main)
  lines(1:ncol(untp),base::colMeans(tamp[row.names(tamp)%in% c,],na.rm=TRUE),type='l',col='orange')
  legend('topleft',col=c('white','black','orange'),lty=1,cex=0.7,legend=c(paste0('n=',n),names(L)[u:t]),bty='n')
  
  plot(1:ncol(untm),base::colMeans(untmc[row.names(untmc)%in% c,],na.rm=TRUE),type='l',col='blue',xlab=xlab,ylab=ylab,ylim=ylim,main=paste(main,'rand Con'))
  lines(1:ncol(untm),base::colMeans(tammc[row.names(tammc)%in% c,],na.rm=TRUE),type='l',col='red')
  legend('topleft',col=c('white','blue','red'),lty=1,cex=0.7,legend=c(paste0('n=',n),'Med13fl_UNT','Med13fl_TAM'),bty='n')
  
  plot(1:ncol(untp),base::colMeans(untpc[row.names(untpc)%in% c,],na.rm=TRUE),type='l',col='black',xlab=xlab,ylab=ylab,ylim=ylim,main=paste(main,'rand Con'))
  lines(1:ncol(untp),base::colMeans(tampc[row.names(tampc)%in% c,],na.rm=TRUE),type='l',col='orange')
  legend('topleft',col=c('white','black','orange'),lty=1,cex=0.7,legend=c(paste0('n=',n),names(L)[u:t]),bty='n')

}


dev.off()

```



##Quantify within peaks scores and norm reads (as for CapC8, Neil's paper)
##plot scatterplot

```{r}


##Do the following just once:

dis=2
peaks=read.delim(paste('ints_all_aggregated_distance',dis,'_summits_annotated.txt',sep=''),stringsAsFactors=FALSE)

d=peaks[peaks$seqnames_bait==peaks$seqnames_otherEnd,]


##1) normalize all L objects:

sk=c(getsk(L))


L_norm=L
for (i in 1:length(L)){
  a=L[[i]]
  a$N=a$N*sk[i]
  a$intID=paste(a$baitID,a$otherEndID,sep=';')
  L_norm[[i]]=a

}

L=L_norm
rm(L_norm)
gc()


##2) annotate peaks with reads and scores in summit and within peak for all samples

for (name in names(L)){
  d[,paste('N_summit',name,sep='_')]=0
  d[,paste('N_tot',name,sep='_')]=0
  d[,paste('score_summit',name,sep='_')]=0
  d[,paste('score_tot',name,sep='_')]=0
}



d$hindids=0

dd=d

for (i in 1:nrow(d)){


  print(i)

  summit=d[i,]$summit
  bait=d[i,]$baitID
  bsID=paste(bait,summit,sep=';')
  st=d[i,]$start_otherEnd
  en=d[i,]$end_otherEnd
  hindids=hind[start(hind)>=st & end(hind)<=en & seqnames(hind)==seqnames(hind[hind$id==bait])]$id
  d[i,]$hindids=length(hindids)

  for (name in names(L)){
    if (nrow(L[[name]][ L[[name]]$intID ==bsID ])>0){
      d[i,][,paste('N_summit',name,sep='_')]=L[[name]][ L[[name]]$intID ==bsID ]$N
      d[i,][,paste('score_summit',name,sep='_')]=L[[name]][ L[[name]]$intID ==bsID ]$score
    }

    if (length(hindids)>1){
      d[i,][,paste('N_tot',name,sep='_')]= sum(L[[name]][ L[[name]]$baitID ==bait & L[[name]]$otherEndID %in% hindids ]$N)
      d[i,][,paste('score_tot',name,sep='_')]= sum(L[[name]][ L[[name]]$baitID ==bait & L[[name]]$otherEndID %in% hindids ]$score)
    }  else {
      d[i,][,paste('N_tot',name,sep='_')]= d[i,][,paste('N_summit',name,sep='_')]
      d[i,][,paste('score_tot',name,sep='_')]= d[i,][,paste('score_summit',name,sep='_')]
    }

  }

}


write.table(d,paste('Peaks_annotated_with_N_and_scores_dis',dis,'.txt',sep=''),sep='\t',eol='\n',quote=FALSE)

d=read.delim(paste('Peaks_annotated_with_N_and_scores_dis',dis,'.txt',sep=''),stringsAsFactors = FALSE)


####################################
##make a table with control peaks:
####################################

##define distance-matched controls:
d=d[order(d$baitID),]
ids=unique(d$baitID)
d$id=1:nrow(d)
d$summit_con=NA
for (id in ids){
  
  df=data.frame(d[d$baitID==id,]$summit,d[d$baitID==id,]$id)
  controls=df[,1]
  con=controls
  con[controls-id > 0]=id-(controls[controls-id > 0]-id)
  con[controls-id < 0]=id+(id-controls[controls-id < 0])
  controls=con
  ##
  chr=seqnames(hind[id])
  cids=hind[seqnames(hind)==chr]$id
  
  if (length(controls[controls %in% cids ==FALSE])>0){
    
  	  message('Controls go outside of the chromosome. Assign chromosome ends to controls...')
      con=controls
  	  con[controls %in% cids ==FALSE & controls>cids[length(cids)] ] = cids[length(cids)]
  	  con[controls %in% cids ==FALSE & controls<cids[1] ] = cids[1]
  	  controls=con
  
    }
  
  d[df[,2],]$summit_con=controls
  
  }

##resize controls to match the size of peaks:

##a) define min and max fragID for each chromosome
chr=unique(d$seqnames_bait)
st=vector()
en=vector()
for (c in chr){
  st=c(st,min(hind[as.character(seqnames(hind))==c]$id))
  en=c(en,max(hind[as.character(seqnames(hind))==c]$id))
}
minmaxids=data.frame(chr=chr,min=st,max=en)

##b) resize the controls to appropriate peak size
gr=GRanges(d$seqnames_bait,IRanges(d$summit_con,d$summit_con))  ##granges with summit_con ID as start/end
gr=resize(gr,d$hindids,fix='center')  ##resize gr based on the number of hindids in the peak (fix center)
gr$min=0
gr$max=0
for (i in 1:nrow(minmaxids)){
  gr[seqnames(gr)==minmaxids[i,]$chr]$min=minmaxids[i,]$min
  gr[seqnames(gr)==minmaxids[i,]$chr]$max=minmaxids[i,]$max
}

##if gr goes beyond the ends of the chromosome, shift the peak to start/end at the end of the chromosome!
gr[start(gr)<gr$min] =shift(gr[start(gr)<gr$min],  gr[start(gr)<gr$min]$min-start(gr[start(gr)<gr$min])     )
gr[end(gr)>gr$max] =shift(gr[end(gr)>gr$max],  gr[end(gr)>gr$max]$max-end(gr[end(gr)>gr$max])     ) 			

##control start fragID and control end fragID are start and end of the gr:
starts=sapply(gr,FUN=function(x) start(hind[start(x)]))
ends=sapply(gr,FUN=function(x) end(hind[end(x)]))
d$seqnames_con=d$seqnames_bait
d$start_con=starts
d$end_con=ends
d$conID=start(gr)

d[width(gr)>1,]$conID=sapply(gr[width(gr)>1],FUN=function(x) concatenate(start(x):end(x)))


##TO DO:

######################################
##Annotate controls with scores and N:
#######################################

# L_norm=L
# for (i in 1:length(L)){
#   a=L[[i]]
#   a$N=a$N*sk[i]
#   a$intID=paste(a$baitID,a$otherEndID,sep=';')
#   L_norm[[i]]=a
# 
# }
# 
# L=L_norm
# rm(L_norm)
# gc()


##2) annotate peaks with reads and scores in summit and within peak for all samples

for (name in names(L)){
  d[,paste('N_summit_con',name,sep='_')]=0
  d[,paste('N_tot_con',name,sep='_')]=0
  d[,paste('score_summit_con',name,sep='_')]=0
  d[,paste('score_tot_con',name,sep='_')]=0
}




dd=d

for (i in 1:nrow(d)){


  print(i)

  summit=d[i,]$summit_con
  bait=d[i,]$baitID
  bsID=paste(bait,summit,sep=';')
  st=d[i,]$start_con
  en=d[i,]$end_con
  hindids=hind[start(hind)>=st & end(hind)<=en & seqnames(hind)==seqnames(hind[hind$id==bait])]$id
  # d[i,]$hindids=length(hindids)

  for (name in names(L)){
    if (nrow(L[[name]][ L[[name]]$intID ==bsID ])>0){
      d[i,][,paste('N_summit_con',name,sep='_')]=L[[name]][ L[[name]]$intID ==bsID ]$N
      d[i,][,paste('score_summit_con',name,sep='_')]=L[[name]][ L[[name]]$intID ==bsID ]$score
    } 
    if (length(hindids)>1){
      d[i,][,paste('N_tot_con',name,sep='_')]= sum(L[[name]][ L[[name]]$baitID ==bait & L[[name]]$otherEndID %in% hindids ]$N)
      d[i,][,paste('score_tot_con',name,sep='_')]= sum(L[[name]][ L[[name]]$baitID ==bait & L[[name]]$otherEndID %in% hindids ]$score)
    }  else {
      d[i,][,paste('N_tot_con',name,sep='_')]= d[i,][,paste('N_summit_con',name,sep='_')]
      d[i,][,paste('score_tot_con',name,sep='_')]= d[i,][,paste('score_summit_con',name,sep='_')]
    }

  }

}


write.table(d,paste('Peaks_annotated_with_N_and_scores_con_dis',dis,'.txt',sep=''),sep='\t',eol='\n',quote=FALSE)
```

##Plots with reads and scores quantifications (ints within specific intervals): 

```{r}

dis=2
d=read.delim(paste('Peaks_annotated_with_N_and_scores_con_dis',dis,'.txt',sep=''),stringsAsFactors = FALSE)

##function to make a list with different types of intervals:

##make a List (LL) with different types of intervals:
List.intIDs <- function(d) {
  
  LL=list(
  Ring1b_peaks_oE=row.names(d[d$Ring1b_peaks_oE,]),
  
  Cdk8_peaks_oE=row.names(d[d$Cdk8_peaks_oE,]),
  
  Ring1b.Cdk8_peaks_oE=row.names(d[d$Ring1b_peaks_oE ==TRUE & d$Cdk8_peaks_oE ==TRUE ,]),
  
  noRing1b.Cdk8_peaks_oE=row.names(d[d$Ring1b_peaks_oE ==FALSE & d$Cdk8_peaks_oE,]),
  
  Ring1b.noCdk8_peaks_oE=row.names(d[d$Ring1b_peaks_oE ==TRUE & d$Cdk8_peaks_oE ==FALSE ,]),
  
  Ring1b.Pcgf2_peaks_oE=row.names(d[d$Pcgf2_peaks_oE & d$Ring1b_peaks_oE,]),
  
  Ring1b.Pcgf2.Suz12_peaks_oE=row.names(d[d$Pcgf2_peaks_oE & d$Ring1b_peaks_oE & d$Suz12_peaks_oE,]),
  
  Pcgf2_peaks_oE=row.names(d[d$Pcgf2_peaks_oE,]),
  
  Cdk8.Pcgf2_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==TRUE & d$Cdk8_peaks_oE ==TRUE ,]),
  
  Cdk8.noPcgf2_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==FALSE & d$Cdk8_peaks_oE,]),
  
  noCdk8.Pcgf2_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==TRUE & d$Cdk8_peaks_oE==FALSE,]),
  
  
  Cdk8.Pcgf2.Ring1b_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==TRUE & d$Cdk8_peaks_oE ==TRUE & d$Ring1b_peaks_oE,]),
  
  Cdk8.Pcgf2.Ring1b.Suz12_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==TRUE & d$Cdk8_peaks_oE ==TRUE & d$Ring1b_peaks_oE & d$Suz12_peaks_oE,]),
  
  Cdk8.noPcgf2.noRing1b.noSuz12_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==FALSE & d$Cdk8_peaks_oE ==TRUE & d$Ring1b_peaks_oE ==FALSE & d$Suz12_peaks_oE ==FALSE ,]),
  
  Atac_peaks_oE=row.names(d[d$atac_oE,]),
  
  noATAC_peaks_oE=row.names(d[d$atac_oE == FALSE,]),
  
  noATAC.noPcG.noCdk8_noTSS_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==FALSE & d$atac_oE ==FALSE & d$Ring1b_peaks_oE==FALSE & d$Cdk8_peaks_oE==FALSE & d$Suz12_peaks_oE ==FALSE & d$RefseqTSS_oE==FALSE,]),
  
  ATAC.noPcgf2.noRing1b.noSuz12_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==FALSE & d$Suz12_peaks_oE ==FALSE & d$Ring1b_peaks_oE==FALSE & d$atac_otherEnd,]),
  
  ATAC.noPcgf2.noRing1b.noSuz12.noCdk8_peaks_oE=row.names(d[d$Pcgf2_peaks_oE ==FALSE & d$atac_oE ==TRUE & d$Ring1b_peaks_oE==FALSE & d$Cdk8_peaks_oE==FALSE & d$Suz12_peaks_oE ==FALSE,]),
  
  TSS_oE=row.names(d[d$RefseqTSS_oE,]),
  
  ATAC.noTSS_oE=row.names(d[d$RefseqTSS_oE==FALSE & d$atac_oE,])
  
)
  
  return(LL)
  
}


###################################
##1) All genes scatterplots
###################################


pdf('Scatterplots_scoresInPeak_byIntervals.pdf',onefile=TRUE)

LL=List.intIDs(d)

for (i in 1:length(LL)){
  c=LL[[i]]
  
  par(mfrow=c(2,2))
  xlab='UNT (asinh)'
  ylab='TAM (asinh)'
  
  dd=d[d$intID %in% c,]
  dd=asinh(dd[,grep('score_tot',names(dd))])
  
  ylims=c(0,max(dd))
  xlims=ylims
  
  main=names(LL)[i]
  n=length(c)
  
 
  smoothScatter( dd$score_tot_Med13fl_UNT, dd$score_tot_Med13fl_TAM,  main=paste('Score',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch='*')
  abline(0,1,lty=2,col='red')
  legend('topleft',bty='n',legend=paste0('n=',n),cex=0.8)
  smoothScatter( dd$score_tot_Pcgf2fl_UNT, dd$score_tot_Pcgf2fl_TAM,  main=paste('score',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch='*')
  abline(0,1,lty=2,col='red')
  
  smoothScatter( dd$score_tot_con_Med13fl_UNT, dd$score_tot_con_Med13fl_TAM,  main=paste('score',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch='*')
  abline(0,1,lty=2,col='red')
  legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
  smoothScatter( dd$score_tot_con_Pcgf2fl_UNT, dd$score_tot_con_Pcgf2fl_TAM,  main=paste('score',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch='*')
  abline(0,1,lty=2,col='red')
  legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
}

dev.off()



pdf('Scatterplots_readsInPeak_byIntervals.pdf',onefile=TRUE)

LL=List.intIDs(d)

for (i in 1:length(LL)){
  c=LL[[i]]
  
  par(mfrow=c(2,2))
  xlab='UNT (log2)'
  ylab='TAM (log2)'
  
  dd=d[d$intID %in% c,]
  p=0.001616
  dd=log2(dd[,grep('N_tot',names(dd))]+p)
  
  ylims=c(min(dd),max(dd))
  xlims=ylims
  
  main=names(LL)[i]
  n=length(c)
  
 
  smoothScatter( dd$N_tot_Med13fl_UNT, dd$N_tot_Med13fl_TAM,  main=paste('Reads',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch='*')
  abline(0,1,lty=2,col='red')
  legend('topleft',bty='n',legend=paste0('n=',n),cex=0.8)
  smoothScatter( dd$N_tot_Pcgf2fl_UNT, dd$N_tot_Pcgf2fl_TAM,  main=paste('Reads',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch='*')
  abline(0,1,lty=2,col='red')
  
  smoothScatter( dd$N_tot_con_Med13fl_UNT, dd$N_tot_con_Med13fl_TAM,  main=paste('Reads',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch='*')
  abline(0,1,lty=2,col='red')
  legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
  smoothScatter( dd$N_tot_con_Pcgf2fl_UNT, dd$N_tot_con_Pcgf2fl_TAM,  main=paste('Reads',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch='*')
  abline(0,1,lty=2,col='red')
  legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
}

dev.off()



########################################
##2) Different genes marked scatterplots
########################################

for (gt in unique(d$genetype)){

pdf(paste0('Scatterplots_scoresInPeak_byIntervals.',gt,'_marked.pdf'),onefile=TRUE)

  LL=List.intIDs(d)
  LLL=List.intIDs(d[d$genetype==gt,])
  
  for (i in 1:length(LL)){
    c=LL[[i]]
    cc=LLL[[i]]
    
    par(mfrow=c(2,2))
    xlab='UNT (asinh)'
    ylab='TAM (asinh)'
    
    dd=d[d$intID %in% c,]
    dd=asinh(dd[,grep('score_tot',names(dd))])
    
    ##only ints of gt:
    ddd=d[d$intID %in% cc,]
    ddd=asinh(ddd[,grep('score_tot',names(ddd))])
    
    ylims=c(0,max(dd))
    xlims=ylims
    
    main=names(LL)[i]
    n=length(c)
    nn=length(cc)
    
   
    smoothScatter( dd$score_tot_Med13fl_UNT, dd$score_tot_Med13fl_TAM,  main=paste('Score',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    points( ddd$score_tot_Med13fl_UNT, ddd$score_tot_Med13fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend=paste0('n=',n),cex=0.8)
    smoothScatter( dd$score_tot_Pcgf2fl_UNT, dd$score_tot_Pcgf2fl_TAM,  main=paste('score',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    points( ddd$score_tot_Pcgf2fl_UNT, ddd$score_tot_Pcgf2fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend=paste0(gt, ', n=',nn),cex=0.8,col='red',pch='*')
    
    smoothScatter( dd$score_tot_con_Med13fl_UNT, dd$score_tot_con_Med13fl_TAM,  main=paste('score',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    points( ddd$score_tot_con_Med13fl_UNT, ddd$score_tot_con_Med13fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
    smoothScatter( dd$score_tot_con_Pcgf2fl_UNT, dd$score_tot_con_Pcgf2fl_TAM,  main=paste('score',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    points( ddd$score_tot_con_Pcgf2fl_UNT, ddd$score_tot_con_Pcgf2fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
  }
  
  dev.off()
  
} 
  


for (gt in unique(d$genetype)){
  
  pdf(paste0('Scatterplots_readsInPeak_byIntervals.',gt,'_marked.pdf'),onefile=TRUE)
  
  LL=List.intIDs(d)
  LLL=List.intIDs(d[d$genetype==gt,])
  
  
  for (i in 1:length(LL)){
    c=LL[[i]]
    cc=LLL[[i]]
    
    par(mfrow=c(2,2))
    xlab='UNT (log2)'
    ylab='TAM (log2)'
    
    dd=d[d$intID %in% c,]
    p=0.001616
    dd=log2(dd[,grep('N_tot',names(dd))]+p)
    
    ##only ints of gt:
    ddd=d[d$intID %in% cc,]
    ddd=log2(ddd[,grep('N_tot',names(ddd))]+p)
    
    ylims=c(min(dd),max(dd))
    xlims=ylims
    
    main=names(LL)[i]
    n=length(c)
    nn=length(cc)
    
   
    smoothScatter( dd$N_tot_Med13fl_UNT, dd$N_tot_Med13fl_TAM,  main=paste('Reads',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    points( ddd$N_tot_Med13fl_UNT, ddd$N_tot_Med13fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend=paste0('n=',n),cex=0.8)
    smoothScatter( dd$N_tot_Pcgf2fl_UNT, dd$N_tot_Pcgf2fl_TAM,  main=paste('Reads',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    points( ddd$N_tot_Pcgf2fl_UNT, ddd$N_tot_Pcgf2fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend=paste0(gt, ', n=',nn),cex=0.8,col='red',pch='*')
    
    
    smoothScatter( dd$N_tot_con_Med13fl_UNT, dd$N_tot_con_Med13fl_TAM,  main=paste('Reads',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    points( ddd$N_tot_con_Med13fl_UNT, ddd$N_tot_con_Med13fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
    smoothScatter( dd$N_tot_con_Pcgf2fl_UNT, dd$N_tot_con_Pcgf2fl_TAM,  main=paste('Reads',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    points( ddd$N_tot_con_Pcgf2fl_UNT, ddd$N_tot_con_Pcgf2fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend=paste0(gt, ', n=',nn),cex=0.8,col='red',pch='*')
    
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
  }
  
  dev.off()


}


########################################
##3) Different genes only scatterplots
########################################

for (gt in unique(d$genetype)){

pdf(paste0('Scatterplots_scoresInPeak_byIntervals.',gt,'.pdf'),onefile=TRUE)

  LL=List.intIDs(d)
  LLL=List.intIDs(d[d$genetype==gt,])
  
  for (i in 1:length(LL)){
    c=LL[[i]]
    cc=LLL[[i]]
    
    par(mfrow=c(2,2))
    xlab='UNT (asinh)'
    ylab='TAM (asinh)'
    
    dd=d[d$intID %in% c,]
    dd=asinh(dd[,grep('score_tot',names(dd))])
    
    ##only ints of gt:
    ddd=d[d$intID %in% cc,]
    ddd=asinh(ddd[,grep('score_tot',names(ddd))])
    
    ylims=c(0,max(dd))
    xlims=ylims
    
    main=names(LL)[i]
    n=length(cc)
    
    
    # 
    # smoothScatter( ddd$score_tot_Med13fl_UNT, ddd$score_tot_Med13fl_TAM,  main=paste('Score',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    # # points( ddd$score_tot_Med13fl_UNT, ddd$score_tot_Med13fl_TAM,pch='*',col='red')
    # abline(0,1,lty=2,col='red')
    # legend('topleft',bty='n',legend=paste0('n=',n),cex=0.8)
    # smoothScatter( ddd$score_tot_Pcgf2fl_UNT, ddd$score_tot_Pcgf2fl_TAM,  main=paste('score',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    # # points( ddd$score_tot_Pcgf2fl_UNT, ddd$score_tot_Pcgf2fl_TAM,pch='*',col='red')
    # abline(0,1,lty=2,col='red')
    # 
    # smoothScatter( ddd$score_tot_con_Med13fl_UNT, ddd$score_tot_con_Med13fl_TAM,  main=paste('score',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    # #points( ddd$score_tot_con_Med13fl_UNT, ddd$score_tot_con_Med13fl_TAM,pch='*',col='red')
    # abline(0,1,lty=2,col='red')
    # legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
    # smoothScatter( ddd$score_tot_con_Pcgf2fl_UNT, ddd$score_tot_con_Pcgf2fl_TAM,  main=paste('score',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    # #points( ddd$score_tot_con_Pcgf2fl_UNT, ddd$score_tot_con_Pcgf2fl_TAM,pch='*',col='red')
    # abline(0,1,lty=2,col='red')
    # legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
    
    plot( ddd$score_tot_Med13fl_UNT, ddd$score_tot_Med13fl_TAM,  main=paste('Score',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch=21, col=rgb(0,0,0,0.3))
    # points( ddd$score_tot_Med13fl_UNT, ddd$score_tot_Med13fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend=paste0('n=',n),cex=0.8)
    plot( ddd$score_tot_Pcgf2fl_UNT, ddd$score_tot_Pcgf2fl_TAM,  main=paste('score',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch=21, col=rgb(0,0,0,0.3))
    # points( ddd$score_tot_Pcgf2fl_UNT, ddd$score_tot_Pcgf2fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    
    plot( ddd$score_tot_con_Med13fl_UNT, ddd$score_tot_con_Med13fl_TAM,  main=paste('score',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch=21, col=rgb(0,0,0,0.3))
    #points( ddd$score_tot_con_Med13fl_UNT, ddd$score_tot_con_Med13fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
    plot( ddd$score_tot_con_Pcgf2fl_UNT, ddd$score_tot_con_Pcgf2fl_TAM,  main=paste('score',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch=21, col=rgb(0,0,0,0.3))
    #points( ddd$score_tot_con_Pcgf2fl_UNT, ddd$score_tot_con_Pcgf2fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
  }
  
  dev.off()
  
} 
  


for (gt in unique(d$genetype)){
  
  pdf(paste0('Scatterplots_readsInPeak_byIntervals.',gt,'.pdf'),onefile=TRUE)
  
  LL=List.intIDs(d)
  LLL=List.intIDs(d[d$genetype==gt,])
  
  
  for (i in 1:length(LL)){
    c=LL[[i]]
    cc=LLL[[i]]
    
    par(mfrow=c(2,2))
    xlab='UNT (log2)'
    ylab='TAM (log2)'
    
    dd=d[d$intID %in% c,]
    p=0.001616
    dd=log2(dd[,grep('N_tot',names(dd))]+p)
    
    ##only ints of gt:
    ddd=d[d$intID %in% cc,]
    ddd=log2(ddd[,grep('N_tot',names(ddd))]+p)
    
    ylims=c(min(dd),max(dd))
    xlims=ylims
    
    main=names(LL)[i]
    n=length(cc)
    # nn=length(cc)
    
   
    # smoothScatter( ddd$N_tot_Med13fl_UNT, ddd$N_tot_Med13fl_TAM,  main=paste('Reads',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    # # points( ddd$N_tot_Med13fl_UNT, ddd$N_tot_Med13fl_TAM,pch='*',col='red')
    # abline(0,1,lty=2,col='red')
    # legend('topleft',bty='n',legend=paste0('n=',n),cex=0.8)
    # smoothScatter( ddd$N_tot_Pcgf2fl_UNT, ddd$N_tot_Pcgf2fl_TAM,  main=paste('Reads',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    # # points( ddd$N_tot_Pcgf2fl_UNT, ddd$N_tot_Pcgf2fl_TAM,pch='*',col='red')
    # abline(0,1,lty=2,col='red')
    # 
    # 
    # smoothScatter( ddd$N_tot_con_Med13fl_UNT, ddd$N_tot_con_Med13fl_TAM,  main=paste('Reads',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    # # points( ddd$N_tot_con_Med13fl_UNT, ddd$N_tot_con_Med13fl_TAM,pch='*',col='red')
    # abline(0,1,lty=2,col='red')
    # legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
    # smoothScatter( ddd$N_tot_con_Pcgf2fl_UNT, ddd$N_tot_con_Pcgf2fl_TAM,  main=paste('Reads',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch='*')
    # # points( ddd$N_tot_con_Pcgf2fl_UNT, ddd$N_tot_con_Pcgf2fl_TAM,pch='*',col='red')
    # abline(0,1,lty=2,col='red')
    # legend('topleft',bty='n',legend=paste0(gt, ', n=',nn),cex=0.8,col='red',pch='*')
    
     plot( ddd$N_tot_Med13fl_UNT, ddd$N_tot_Med13fl_TAM,  main=paste('Reads',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch=21,col=rgb(0,0,0,0.3))
    # points( ddd$N_tot_Med13fl_UNT, ddd$N_tot_Med13fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend=paste0('n=',n),cex=0.8)
    plot( ddd$N_tot_Pcgf2fl_UNT, ddd$N_tot_Pcgf2fl_TAM,  main=paste('Reads',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch=21,col=rgb(0,0,0,0.3))
    # points( ddd$N_tot_Pcgf2fl_UNT, ddd$N_tot_Pcgf2fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    
    
    plot( ddd$N_tot_con_Med13fl_UNT, ddd$N_tot_con_Med13fl_TAM,  main=paste('Reads',main), xlab=paste('Med13fl', xlab), ylab=paste('Med13fl', ylab),ylim=ylims,xlim=xlims,pch=21,col=rgb(0,0,0,0.3))
    # points( ddd$N_tot_con_Med13fl_UNT, ddd$N_tot_con_Med13fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    legend('topleft',bty='n',legend='dist_matched con',cex=0.8)
    plot( ddd$N_tot_con_Pcgf2fl_UNT, ddd$N_tot_con_Pcgf2fl_TAM,  main=paste('Reads',main), xlab=paste('Pcgf2fl', xlab), ylab=paste('Pcgf2fl', ylab),ylim=ylims,xlim=xlims,pch=21,col=rgb(0,0,0,0.3))
    # points( ddd$N_tot_con_Pcgf2fl_UNT, ddd$N_tot_con_Pcgf2fl_TAM,pch='*',col='red')
    abline(0,1,lty=2,col='red')
    
  }
  
  dev.off()

}
```


#####################################
##ggboxplots reads and scores
#####################################

```{r}

dis=2
d=read.delim(paste('Peaks_annotated_with_N_and_scores_con_dis',dis,'.txt',sep=''),stringsAsFactors = FALSE)

cols=c('red','orange','darkblue','blue')

##make list with intIDs in intervals of interest!
LL=List.intIDs(d)
Lg=list()

##all interactions:

name='all_interactions'
a=d
##make sure to only look at interactions with a score of >0 in at least one treatment in the resp cell line:
am=a[(a$score_tot_Med13fl_UNT+a$score_tot_Med13fl_TAM)>0,]
ap=a[(a$score_tot_Pcgf2fl_UNT+a$score_tot_Pcgf2fl_TAM)>0,]

dd=data.frame(score_per_frag=c(am$score_tot_Med13fl_UNT/am$hindids,am$score_tot_Med13fl_TAM/am$hindids, ap$score_tot_Pcgf2fl_UNT/ap$hindids, ap$score_tot_Pcgf2fl_TAM/ap$hindids),
                # celline=rep(c('Med13fl_UNT','Med13fl_TAM','Pcgf2fl_UNT','Pcgf2fl_TAM'),each=nrow(a)),
                celline=c(rep(c('Med13fl_UNT','Med13fl_TAM'),each=nrow(am)),rep(c('Pcgf2fl_UNT','Pcgf2fl_TAM'),each=nrow(ap))),
                norm_reads_per_frag=c(am$N_tot_Med13fl_UNT/am$hindids,am$N_tot_Med13fl_TAM/am$hindids, ap$N_tot_Pcgf2fl_UNT/ap$hindids, ap$N_tot_Pcgf2fl_TAM/ap$hindids),
                genetype=c(rep(am$genetype,2),rep(ap$genetype,2)),
                oE_type=c(rep(am$type_oE,2),rep(ap$type_oE,2)))
  dd$celline=factor(dd$celline,levels=c('Med13fl_UNT','Med13fl_TAM','Pcgf2fl_UNT','Pcgf2fl_TAM'))
  dd$oE_type=factor(dd$oE_type,levels=c('atac', 'Cdk8_noPcG',     'Cdk8_andPcG',  'Cdk8_and.cPRC1',     'polycomb',        'cPRC1', 'none'))
  dd$genetype=factor(dd$genetype,levels=c('inactive','constitutively_active','pluripotent','RA72h_repressed','induced_Med13ESderepressed','induced_Med13independent','induced_Med13dependent'))

p <- ggplot(dd,aes(x=oE_type, y=score_per_frag, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,12.5))
p <- p +  theme_minimal() + scale_fill_manual(values=cols) + ggtitle(paste(name, 'scores'))
Lg=c(Lg,list(p))

p <- ggplot(dd,aes(x=oE_type, y=norm_reads_per_frag, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,0.13))
p <- p +  theme_minimal() + scale_fill_manual(values=cols) + ggtitle(paste(name, 'reads'))
Lg=c(Lg,list(p))

  
p <- ggplot(dd,aes(x=genetype, y=score_per_frag, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,12.5))
p <- p +  theme_minimal() + scale_fill_manual(values=cols) + ggtitle(paste(name, 'scores'))
Lg=c(Lg,list(p))

p <- ggplot(dd,aes(x=genetype, y=norm_reads_per_frag, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,0.13))
p <- p +  theme_minimal() + scale_fill_manual(values=cols) + ggtitle(paste(name, 'reads'))
Lg=c(Lg,list(p))




##special interactions:

for (i in 1:length(LL)){

  name=names(LL)[i]
  a=d[d$intID %in% LL[[i]],]
  ##make sure to only look at interactions with a score of >0 in at least one treatment in the resp cell line:
  am=a[(a$score_tot_Med13fl_UNT+a$score_tot_Med13fl_TAM)>0,]
  ap=a[(a$score_tot_Pcgf2fl_UNT+a$score_tot_Pcgf2fl_TAM)>0,]
  
  dd=data.frame(score_per_frag=c(am$score_tot_Med13fl_UNT/am$hindids,am$score_tot_Med13fl_TAM/am$hindids, ap$score_tot_Pcgf2fl_UNT/ap$hindids, ap$score_tot_Pcgf2fl_TAM/ap$hindids),
                # celline=rep(c('Med13fl_UNT','Med13fl_TAM','Pcgf2fl_UNT','Pcgf2fl_TAM'),each=nrow(a)),
                celline=c(rep(c('Med13fl_UNT','Med13fl_TAM'),each=nrow(am)),rep(c('Pcgf2fl_UNT','Pcgf2fl_TAM'),each=nrow(ap))),
                norm_reads_per_frag=c(am$N_tot_Med13fl_UNT/am$hindids,am$N_tot_Med13fl_TAM/am$hindids, ap$N_tot_Pcgf2fl_UNT/ap$hindids, ap$N_tot_Pcgf2fl_TAM/ap$hindids),
                genetype=c(rep(am$genetype,2),rep(ap$genetype,2)))
  dd$celline=factor(dd$celline,levels=c('Med13fl_UNT','Med13fl_TAM','Pcgf2fl_UNT','Pcgf2fl_TAM'))
  dd$genetype=factor(dd$genetype,levels=c('inactive','constitutively_active','pluripotent','RA72h_repressed','induced_Med13ESderepressed','induced_Med13independent','induced_Med13dependent'))

  p <- ggplot(dd,aes(x=genetype, y=score_per_frag, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,12.5))
  p <- p +  theme_minimal() + scale_fill_manual(values=cols) + ggtitle(paste(name, 'scores'))
  Lg=c(Lg,list(p))
  
  p <- ggplot(dd,aes(x=genetype, y=norm_reads_per_frag, fill=celline)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + coord_cartesian(ylim=c(0,0.13))
  p <- p +  theme_minimal() + scale_fill_manual(values=cols) + ggtitle(paste(name, 'reads'))
  Lg=c(Lg,list(p))

}

pdf('ggBoxplots_byIntervals_byGenetypes.pdf',onefile=TRUE, width=14,height=10)
for (i in (1:length(Lg))[(1:length(Lg))%%2==1]){
  multiplot(Lg[[i]],Lg[[(i+1)]],cols=1)
}
dev.off()



```









#################################
##interesting params
#################################

##gene types:
         induced    increased_expr         repressed         active_ES     stable_active   stable_inactive 
"TRUE :27       " "TRUE :64       " "TRUE :45       " "TRUE :179      " "TRUE :97       " "TRUE :14       "

```{r}



samples=c('0h', '2h', '4h','6h','24h','48h','72h')    

##just once:
# L_interesting_params=list()
# for (s in samples){
#   interesting_params=getInterestingParams(s,hind,DesignDir)
#   L_interesting_params=c(L_interesting_params,list(interesting_params))
# }

##after everything has been saved:

L_interesting_params=list()
for (s in samples){
  L_interesting_params=c(L_interesting_params,list(read.delim(paste0('results/interesting_params_',s,'.txt'))))
}

names(L_interesting_params)=samples

plotInterestingParams(L_interesting_params)


##general plotting parameters:
col=brewer.pal(7,'PiYG')
timepoints=c(0,2,4,6,24,48,72)

##Function to plot barplots with numbers:
barplot_numbers=function(v,names,col,main=NULL,xlab=NULL,ylab=NULL){
  bp=barplot(v,names=timepoints,col=col, main=main,xlab=xlab,ylab=ylab,ylim=c(0,max(v)+0.05*(max(v))))
  text(bp,max(v)+0.04*(max(v)),labels=v)
}

##Function to plot the number and distances of significant interactions as barplots:

bp_interestingparams <- function(L_interesting_params){

  v=vector()
  for (i in 1:length(L_interesting_params)){
    a=L_interesting_params[[i]]
    v=c(v,sum(a$Number_signif_interactions))
  }
  
  barplot_numbers(v=v,names=timepoints,col=col,main='total sig ints', xlab='hours RA treatment', ylab='number')
  
  v=vector()
  for (i in 1:length(L_interesting_params)){
    a=L_interesting_params[[i]]
    v=c(v,median(a$Number_signif_interactions))
  }
  
  barplot_numbers(v=v,names=timepoints,col=col,main='median sig ints per prom', xlab='hours RA treatment', ylab='number')
  
  v=vector()
  for (i in 1:length(L_interesting_params)){
    a=L_interesting_params[[i]]
    v=c(v,round(mean(a$Number_signif_interactions),2))
  }
  
  barplot_numbers(v=v,names=timepoints,col=col,main='mean sig ints per prom', xlab='hours RA treatment', ylab='number')
  
  
  v=vector()
  for (i in 1:length(L_interesting_params)){
    a=L_interesting_params[[i]]
    v=c(v,round(median(a$Distance_signif_median,na.rm=TRUE)/1000,1))
  }
  
  barplot_numbers(v=v,names=timepoints,col=col,main='median median Distance sig ints per prom [kb]', xlab='hours RA treatment', ylab='number')
  
  v=vector()
  for (i in 1:length(L_interesting_params)){
    a=L_interesting_params[[i]]
    v=c(v,round(mean(a$Distance_signif_median,na.rm=TRUE)/1000,1))
  }
  barplot_numbers(v=v,names=timepoints,col=col,main='mean median Distance sig ints per prom [kb]', xlab='hours RA treatment', ylab='number')
  
  
  v=vector()
  for (i in 1:length(L_interesting_params)){
    a=L_interesting_params[[i]]
    v=c(v,round(median(a$Distance_signif_mean,na.rm=TRUE)/1000,1))
  }
  
  barplot_numbers(v=v,names=timepoints,col=col,main='median mean Distance sig ints per prom [kb]', xlab='hours RA treatment', ylab='number')
  
  v=vector()
  for (i in 1:length(L_interesting_params)){
    a=L_interesting_params[[i]]
    v=c(v,round(mean(a$Distance_signif_mean,na.rm=TRUE)/1000,1))
  }
  barplot_numbers(v=v,names=timepoints,col=col,main='mean mean Distance sig ints per prom [kb]', xlab='hours RA treatment', ylab='number')
  
}


##now make interesting plots for different types of genes (interesting genes)
h=3
w=3

##1) all genes:
png('plots/interesting_params_all_baited_genes.png',width=350*w, height=400*h)
par(mfrow=c(h,w))
bp_interestingparams(L_interesting_params)
dev.off()

##2) different sets of genes:

interesting_genes=c("induced" ,"increased_expr",  "repressed" ,      "active_ES" ,      "stable_active" ,  "stable_inactive")

for (ig in interesting_genes){
  
  LL=L_interesting_params
  
  for (i in 1:length(LL)){
    LL[[i]]=LL[[i]][LL[[i]]$baitID %in% as.data.frame(values(baited_genes))[as.data.frame(values(baited_genes))[,ig]==TRUE,]$dpn_id ,]
  
  }
  
png(paste0('plots/interesting_params_',ig,'_baited_genes.png'),width=350*w, height=400*h)
par(mfrow=c(h,w))
bp_interestingparams(LL)
dev.off()

}

######################################
##Function to plot the number and distances of significant interactions as lineplots
######################################

Lp_interestingparams <- function(L_interesting_params,param,main,ylim,ylab,type='median'){
  ##param: character, name of a column in data farmes from L_interesting_param: which parameter should be plotted?
  ##main: character, plot title
  ##type: character, median (default) or mean

    interesting_genes=c("induced" ,"increased_expr",  "repressed" ,      "active_ES" ,      "stable_active" ,  "stable_inactive")
    cols=c('green','darkgreen','red','blue','purple','grey')
    df=data.frame(interesting_genes=interesting_genes,col=cols)
  
    timepoint=c(0,2,4,6,24,48,72)
    
    v=vector()
    vv=vector()
    for (i in 1:length(L_interesting_params)){
      a=L_interesting_params[[i]]
      if (type=='median'){
      v=c(v,median(a[,param],na.rm=TRUE))
      } else {v=c(v,mean(a[,param],na.rm=TRUE))}
      # vv=c(vv,sd(a$Number_signif_interactions))
    }
    
    xlab='timepoint'
    # ylab='n sigInts per Prom (median)'
    plot( timepoint , v, col='black', ylim=ylim , xlab=xlab,ylab=ylab, main=main,pch=20)
    lines( timepoint , v,col='black')
    legend('topright',legend=c('all',interesting_genes),col=c('black',cols),lty=1,bty='n', cex=0.7)
    
    
  for (ig in interesting_genes){
    LL=L_interesting_params
    col=as.character(df[df$interesting_genes==ig,]$col)
    
    v=vector()
    vv=vector()
    for (i in 1:length(LL)){
        a=LL[[i]][LL[[i]]$baitID %in% as.data.frame(values(baited_genes))[as.data.frame(values(baited_genes))[,ig]==TRUE,]$dpn_id ,]
        if (type=='median'){
      v=c(v,median(a[,param],na.rm=TRUE))
      } else {v=c(v,mean(a[,param],na.rm=TRUE))}
        # vv=c(vv,sd(a$Number_signif_interactions))
    }
    points( timepoint,v, col=col,pch=20 )
    lines( timepoint,v, col=col )
    
  }
  
}

type='median'
main=paste(type, 'signif Ints per promoter')
param='Number_signif_interactions'
ylab=paste0('n sigInts per Prom (',type, ')')
ylim=c(0,11)
Lp_interestingparams(L_interesting_params,param,main,ylim,ylab,type)

type='mean'
main=paste(type, 'signif Ints per promoter')
param='Number_signif_interactions'
ylab=paste0('n sigInts per Prom (',type, ')')
ylim=c(8,23)
Lp_interestingparams(L_interesting_params,param,main,ylim,ylab,type)

type='median'
main=paste(type, 'meanDistance signif Ints per promoter (prom with Ints)')
param='Distance_signif_mean'
ylab=paste0('meanDistance sigInt per Prom with sigInt (',type, ')')
ylim=c(0,700000)
Lp_interestingparams(L_interesting_params,param,main,ylim,ylab,type)

type='mean'
main=paste(type, 'meanDistance signif Ints per promoter (prom with Ints)')
param='Distance_signif_mean'
ylab=paste0('meanDistance sigInt per Prom with sigInt (',type, ')')
ylim=c(0,2700000)
Lp_interestingparams(L_interesting_params,param,main,ylim,ylab,type)


################################################
##Function to plot density plots for all samples
##at all timepoints
###############################################

dp_interestingparams <- function(L_interesting_params,param,main,ylim,xlim,xlab){
  ##param: character, name of a column in data farmes from L_interesting_param: which parameter should be plotted?
  ##main: character, plot title
  ##type: character, median (default) or mean

    interesting_genes=c("induced" ,"increased_expr",  "repressed" ,      "active_ES" ,      "stable_active" ,  "stable_inactive")
    cols=c('green','darkgreen','red','blue','purple','grey')
    df=data.frame(interesting_genes=interesting_genes,col=cols)
  
    timepoint=c(0,2,4,6,24,48,72)
  
    for (i in 1:length(L_interesting_params)){  
      
      plot(density(L_interesting_params[[i]][,param]),main=paste(main,names(L_interesting_params)[i]),xlab=xlab,ylab=ylab,ylim=ylim,xlim=xlim)
      legend('topright',legend=c('all',interesting_genes),col=c('black',cols),lty=1,bty='n')
    
      for (ig in interesting_genes){ 
        col=as.character(df[df$interesting_genes==ig,]$col)
        lines(density(L_interesting_params[[i]][L_interesting_params[[i]]$baitID %in% as.data.frame(values(baited_genes))[as.data.frame(values(baited_genes))[,ig]==TRUE,]$dpn_id , param]),col=col)
      }
    }
  
}

xlim=c(-20,140)
ylim=c(0,0.09)
main='signif Ints per promoter'
param='Number_signif_interactions'
xlab=paste('n', main)
pdf(paste0('plots/',param,'_perProm_each_timePoint.pdf'),width=14,height=14)
par(mfrow=c(2,4))
dp_interestingparams(L_interesting_params, ylim=ylim,xlab=xlab,xlim=xlim,main=main,param=param)
dev.off()

################################################
##Function to plot density plots for each sample
##over many timepoints
###############################################

dp_interestingparams <- function(L_interesting_params,param,main,xlim,col,ylim){
  ##param: character, name of a column in data farmes from L_interesting_param: which parameter should be plotted?
  ##main: character, plot title
  ##type: character, median (default) or mean

    interesting_genes=c("induced" ,"increased_expr",  "repressed" ,      "active_ES" ,      "stable_active" ,  "stable_inactive")
    cols=c('green','darkgreen','red','blue','purple','grey')
    df=data.frame(interesting_genes=interesting_genes,col=cols)
    
  for (ig in interesting_genes){ 
    
    i=1
    plot(density(L_interesting_params[[i]][L_interesting_params[[i]]$baitID %in% as.data.frame(values(baited_genes))[as.data.frame(values(baited_genes))[,ig]==TRUE,]$dpn_id , param]),col=col[i],main=paste(main,ig),xlim=xlim,ylim=ylim)
    legend('topright',bty='n',col=col,lwd=2,legend=names(L_interesting_params))
    
    for (i in 2:length(L_interesting_params)){  
      
        lines(density(L_interesting_params[[i]][L_interesting_params[[i]]$baitID %in% as.data.frame(values(baited_genes))[as.data.frame(values(baited_genes))[,ig]==TRUE,]$dpn_id , param]),col=col[i])
      }
    }
  
}

##general plotting parameters:
col=brewer.pal(7,'PiYG')
timepoints=c(0,2,4,6,24,48,72)

xlim=c(-10,130)
ylim=c(0,0.09)
main='signif Ints per promoter'
param='Number_signif_interactions'
xlab=paste('n', main)
pdf(paste0('plots/',param,'_perProm_each_geneType.pdf'),width=14,height=14)
par(mfrow=c(2,3))
dp_interestingparams(L_interesting_params, ylim=ylim,xlim=xlim,main=main,param=param,col=col)
dev.off()


#################################
##Function to plot density plots
##for distance of interactions
##each gene type
################################

dp_interestingparams <- function(L_interesting_params,param,main,xlim,col,ylim){
  ##param: character, name of a column in data farmes from L_interesting_param: which parameter should be plotted?
  ##main: character, plot title
  ##type: character, median (default) or mean

    interesting_genes=c("induced" ,"increased_expr",  "repressed" ,      "active_ES" ,      "stable_active" ,  "stable_inactive")
    cols=c('green','darkgreen','red','blue','purple','grey')
    df=data.frame(interesting_genes=interesting_genes,col=cols)
    
  for (ig in interesting_genes){ 
    
    i=1
    plot(density(log10(L_interesting_params[[i]][L_interesting_params[[i]]$baitID %in% as.data.frame(values(baited_genes))[as.data.frame(values(baited_genes))[,ig]==TRUE,]$dpn_id , param]/1000),na.rm=TRUE),col=col[i],main=paste(main,ig),xlim=xlim,ylim=ylim)
    legend('topright',bty='n',col=col,lwd=2,legend=names(L_interesting_params))
    
    for (i in 2:length(L_interesting_params)){  
      
        lines(density(log10(L_interesting_params[[i]][L_interesting_params[[i]]$baitID %in% as.data.frame(values(baited_genes))[as.data.frame(values(baited_genes))[,ig]==TRUE,]$dpn_id , param]/1000),na.rm=TRUE),col=col[i])
      }
    }
  
}

##general plotting parameters:
col=brewer.pal(7,'PiYG')
timepoints=c(0,2,4,6,24,48,72)

xlim=c(1,5)
ylim=c(0,2)
main='mean Distance per Prom (log10 kbp)'
param='Distance_signif_mean'
xlab=paste('n', main)
pdf(paste0('plots/',param,'_perProm_each_geneType.pdf'),width=14,height=14)
par(mfrow=c(2,3))
dp_interestingparams(L_interesting_params, ylim=ylim,xlim=xlim,main=main,param=param,col=col)
dev.off()

##########################################################
##number of genes with 0 interactions at each time point:
##########################################################

pdf(paste0('plots/Fraction_genes_with_0Ints_each_geneType.pdf'),width=14,height=14)
par(mfrow=c(2,3))

 interesting_genes=c("induced" ,"increased_expr",  "repressed" ,      "active_ES" ,      "stable_active" ,  "stable_inactive")
    cols=c('green','darkgreen','red','blue','purple','grey')
    df=data.frame(interesting_genes=interesting_genes,col=cols)
    
  for (ig in interesting_genes){ 
    
    v=vector()
    for (i in 1:length(L_interesting_params)){
      a=L_interesting_params[[i]][L_interesting_params[[i]]$baitID %in% as.data.frame(values(baited_genes))[as.data.frame(values(baited_genes))[,ig]==TRUE,]$dpn_id , ]
      v=c(v, round(nrow(a[a[,param]==0,])/nrow(a),2))
    }
    barplot_numbers(v,names=names(L_interesting_params),col=col,main=paste('Fraction genes with 0 interactions',ig),ylab='Fraction of all genes',xlab='hours RA treatment')
      
  }

dev.off()

    
    
```


##Clustering of interactions:

```{r}

mat=as.matrix(ints[,2:8])
row.names(mat)=row.names(ints)

pheatmap(asinh(mat))
pheatmap(asinh(mat),kmeans_k = 5)

set.seed(1984)

##hierarchical clustering:
my_hclust_gene <- hclust(dist(asinh(mat)), method = "complete")
 
# # install if necessary
# install.packages("dendextend")
#  
# load package
library(dendextend)
 
png('plots/Cluster_analysis_hclust_tree.png')

as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)

dev.off()

##cut the tree into 3 clusters and assign a cluster to each interaction:
clmax=10
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = clmax)
head(my_gene_col)

my_gene_col=data.frame(name=names(my_gene_col),cluster=paste0('cluster',my_gene_col))
row.names(my_gene_col)=my_gene_col[,1]

##add random "clusters"
my_random <- as.factor(sample(x = 1:clmax, size = nrow(my_gene_col), replace = TRUE))
my_gene_col$random <- my_random

##add annotation with gene type:
b=unlist(strsplit(as.character(my_gene_col$name),split=';'))[1:length(unlist(strsplit(as.character(my_gene_col$name),split=';')))%%2==1]
my_gene_col$genetype='NA'
my_gene_col[b %in% c(baited_genes[baited_genes$induced]$dpn_id,baited_genes[baited_genes$increased_expr]$dpn_id) ,]$genetype='induced'
my_gene_col[b %in% c(baited_genes[baited_genes$repressed]$dpn_id) ,]$genetype='repressed'
my_gene_col[b %in% c(baited_genes[baited_genes$stable_active]$dpn_id) ,]$genetype='stably_active'
my_gene_col[b %in% c(baited_genes[baited_genes$stable_inactive]$dpn_id) ,]$genetype='stably_inactive'
my_gene_col=my_gene_col[,-1]


#plot heatmaps:
pheatmap(asinh(mat), annotation_row = my_gene_col,cluster_cols = FALSE,show_rownames = FALSE,filename='plots/pheatmap_allScores.png')
pheatmap(asinh(mat), kmeans_k=clmax,cluster_cols=FALSE,filename='plots/pheatmap_allScores_kmeans.png')

for (cl in 1:clmax){
  pheatmap(asinh(mat[row.names(mat) %in% row.names(my_gene_col[my_gene_col$cluster==paste0('cluster',cl),]),]), annotation_row = my_gene_col[my_gene_col$cluster==paste0('cluster',cl),],cluster_cols = FALSE,show_rownames = FALSE,cluster_rows=FALSE,main = paste('cluster',cl),filename=paste0('plots/pheatmap_allScores_cluster',cl,'.png'))
}

##plot average cluster profiles:

##mean:
cl=1
d=data.frame(cluster1=colMeans(asinh(mat[row.names(mat) %in% row.names(my_gene_col[my_gene_col$cluster==paste0('cluster',cl),]),])))
for (cl in 2:clmax){
  d=cbind(d,colMeans(asinh(mat[row.names(mat) %in% row.names(my_gene_col[my_gene_col$cluster==paste0('cluster',cl),]),])))
}
names(d)=paste0('cluster',1:clmax)
row.names(d)=colnames(mat)
dd=t(d)
m=as.matrix(dd)
pheatmap(m,cluster_rows = FALSE,cluster_cols = FALSE,display_numbers = TRUE,main='colmean each cluster sigInt: >2.31',filename='plots/pheatmap_allScores_hclust_meanPerCluster.png')

##median:
cl=1
d=data.frame(cluster1=colMedians(asinh(mat[row.names(mat) %in% row.names(my_gene_col[my_gene_col$cluster==paste0('cluster',cl),]),])))
for (cl in 2:clmax){
  d=cbind(d,colMedians(asinh(mat[row.names(mat) %in% row.names(my_gene_col[my_gene_col$cluster==paste0('cluster',cl),]),])))
}
names(d)=paste0('cluster',1:clmax)
row.names(d)=colnames(mat)
dd=t(d)
m=as.matrix(dd)
pheatmap(m,cluster_rows = FALSE,cluster_cols = FALSE,display_numbers = TRUE,main='colmedian each cluster sigInt: >2.31',filename='plots/pheatmap_allScores_hclust_medianPerCluster.png')

##fraction each gene type per profile:
my_gene_col$genetype=factor(my_gene_col$genetype,levels=c('induced','repressed','stably_active','stably_inactive'))
tab=table(my_gene_col$genetype,my_gene_col$cluster)

png('plots/hclust_allScores_enrichment_of_genetypes.png')
plotNormBarplot(tab,main='hclust scores',las=3,ylab='Fraction genes',col=c('green','red','purple','grey'))
dev.off()


```


##Identify all significant interactions between all samples:

```{r}

LL=L

ids=vector()
for (i in 1:length(L)){
  a=L[[i]]
  a$intID=paste(a$baitID,a$otherEndID,sep=';')
  LL[[i]]=a
  a=a[a$score>=5,]
  ids=c(ids,a$intID)
}

L=LL

ids=unique(ids)

ints=data.frame(intID=ids[order(ids)])

for (i in 1:length(L)){
  a=L[[i]]
  a=a[order(a$intID),]
  a=a[a$intID %in% ids,]
  ints[,paste(names(L)[i],'score',sep='_')]=0
  ints[,paste(names(L)[i],'N',sep='_')]=0
  
  print(summary(ints[ints$intID %in% a$intID,]$intID==a$intID))
  
  ints[ints$intID %in% a$intID,paste(names(L)[i],'score',sep='_')]=a$score
  ints[ints$intID %in% a$intID,paste(names(L)[i],'N',sep='_')]=a$N
  
  ##get pooled and replicate read numbers
  
  ##how many replicates?
  x=names(a)
  x=x[grep('N',x)]
  x=x[grep('[.]',x)]
  n=as.numeric(unlist(strsplit(x,split='[.]'))[length(unlist(strsplit(x,split='[.]')))])
  for (j in 1:n){
    ints[,paste0(names(L)[i],'_N.',j)]=0
    ints[ints$intID %in% a$intID,paste0(names(L)[i],'_N.',j)]=as.data.frame(a)[,paste0('N.',j)]
  }
  
}

ints$intID=as.character(ints$intID)
ints$baitID=as.numeric(unlist(strsplit(ints$intID,split=';'))[1:length(unlist(strsplit(ints$intID,split=';')))%%2==1])
ints$otherEndID=as.numeric(unlist(strsplit(ints$intID,split=';'))[1:length(unlist(strsplit(ints$intID,split=';')))%%2==0])

oe=gr2bed(hind[findOverlapsDf(ints,hind,'otherEndID','id')[,2]])
b=gr2bed(hind[findOverlapsDf(ints,hind,'baitID','id')[,2]])
names(oe)=paste(names(oe),'otherEnd',sep='_')
names(b)=paste(names(b),'bait',sep='_')
names(oe)[1]='seqnames_otherEnd'
names(b)[1]='seqnames_bait'
ints=cbind(ints,oe,b)
row.names(ints)=ints$intID

ints=ints[order(ints$otherEndID),]
ints=ints[order(ints$baitID),]

##normalize read counts by the read numbers in the respective libraries (written to vecor v):

m=ints[,grep('[.]',names(ints))]

v=vector()
for (i in 1:length(L)){
  a=L[[i]]
  
  ##how many replicates?
  x=names(a)
  x=x[grep('N',x)]
  x=x[grep('[.]',x)]
  n=as.numeric(unlist(strsplit(x,split='[.]'))[length(unlist(strsplit(x,split='[.]')))])
  for (j in 1:n){
    v=c(v,sum(as.data.frame(a)[,paste0('N.',j)]))
    names(v)[length(v)]=paste0(names(L)[i],'_N.',j)
  }
}
libsize=v

##downsample to smallest:
names(m)==names(libsize)
mat=m
for (i in 1:ncol(m)){
  m[,i]=m[,i]/libsize[i]*min(libsize)
}

names(m)=paste(names(m),'downsampled',sep='_')
ints=cbind(ints,m)

##annotate ints with gene:
ints$genename='NA'
ints$genetype='NA'
for (b in unique(ints$baitID)){
  print(b)
  n=as.character(baited_genes[baited_genes$dpn_id==b,]$genename)
  typ=as.character(baited_genes[baited_genes$dpn_id==b,]$type)
  ints[ints$baitID==b,]$genename=n
  ints[ints$baitID==b,]$genetype=typ
}
head(ints)

##annotate with peaks of PRC1, CDK8 etc
b=bed2gr(ints[,c('seqnames_bait','start_bait','end_bait')])

b$cdk8=FALSE
b$ring1b=FALSE
b$pcgf2=FALSE
b$suz12=FALSE
b$atac=FALSE
b$tss=FALSE
b[as.matrix(findOverlaps(b,cdk))[,1]]$cdk8=TRUE
b[as.matrix(findOverlaps(b,ring))[,1]]$ring1b=TRUE
b[as.matrix(findOverlaps(b,pcgf))[,1]]$pcgf2=TRUE
b[as.matrix(findOverlaps(b,suz))[,1]]$suz12=TRUE
b[as.matrix(findOverlaps(b,atac))[,1]]$atac=TRUE
b[as.matrix(findOverlaps(b,refseqtss))[,1]]$tss=TRUE

names(values(b))=paste0(names(values(b)),'_bait')
ints=cbind(ints,as.data.frame(values(b)))

b=bed2gr(ints[,c('seqnames_otherEnd','start_otherEnd','end_otherEnd')])

b$cdk8=FALSE
b$ring1b=FALSE
b$pcgf2=FALSE
b$suz12=FALSE
b$atac=FALSE
b$tss=FALSE
b[as.matrix(findOverlaps(b,cdk))[,1]]$cdk8=TRUE
b[as.matrix(findOverlaps(b,ring))[,1]]$ring1b=TRUE
b[as.matrix(findOverlaps(b,pcgf))[,1]]$pcgf2=TRUE
b[as.matrix(findOverlaps(b,suz))[,1]]$suz12=TRUE
b[as.matrix(findOverlaps(b,atac))[,1]]$atac=TRUE
b[as.matrix(findOverlaps(b,refseqtss))[,1]]$tss=TRUE

names(values(b))=paste0(names(values(b)),'_otherEnd')
ints=cbind(ints,as.data.frame(values(b)))


##normalize total read counts by the read numbers in the respective libraries (written to vecor v):

m=ints[,grep('_N$',names(ints))]

v=vector()
for (i in 1:length(L)){
  a=L[[i]]
  
    v=c(v,sum(as.data.frame(a)[,'N']))
    names(v)[length(v)]=paste0(names(L)[i],'_N')
  
}
libsize=v

##downsample to smallest:
names(m)==names(libsize)
mat=m
for (i in 1:ncol(m)){
  m[,i]=m[,i]/libsize[i]*min(libsize)
}

names(m)=paste(names(m),'downsampled',sep='_')
ints=cbind(ints,m)

##just once:
# write.table(ints,'ints_all.txt',sep='\t',eol='\n',quote=FALSE)

```


####################################################


###################################################
########################################
##Identify and plot "interaction peaks"
########################################
###################################################


Strategy:

1. Take all "significant interactions"

2. Identify "summits": If several adjacent fragments are significant take the one with the highest read number as "summit"

3. Plot summit +/- n fragments around the summit



```{r}

LL=L

ids=vector()
for (i in 1:length(L)){
  a=L[[i]]
  a$intID=paste(a$baitID,a$otherEndID,sep=';')
  LL[[i]]=a
  a=a[a$score>=5,]
  ids=c(ids,a$intID)
}

L=LL

ids=unique(ids)

ints=data.frame(intID=ids[order(ids)])

for (i in 1:length(L)){
  a=L[[i]]
  a=a[order(a$intID),]
  a=a[a$intID %in% ids,]
  ints[,names(L)[i]]=0
  
  print(summary(ints[ints$intID %in% a$intID,]$intID==a$intID))
  
  ints[ints$intID %in% a$intID,names(L)[i]]=a$score
}

ints$intID=as.character(ints$intID)
ints$baitID=as.numeric(unlist(strsplit(ints$intID,split=';'))[1:length(unlist(strsplit(ints$intID,split=';')))%%2==1])
ints$otherEndID=as.numeric(unlist(strsplit(ints$intID,split=';'))[1:length(unlist(strsplit(ints$intID,split=';')))%%2==0])

oe=gr2bed(hind[findOverlapsDf(ints,hind,'otherEndID','id')[,2]])
b=gr2bed(hind[findOverlapsDf(ints,hind,'baitID','id')[,2]])
names(oe)=paste(names(oe),'otherEnd',sep='_')
names(b)=paste(names(b),'bait',sep='_')
names(oe)[1]='seqnames_otherEnd'
names(b)[1]='seqnames_bait'
ints=cbind(ints,oe,b)
row.names(ints)=ints$intID

ints=ints[order(ints$otherEndID),]
ints=ints[order(ints$baitID),]

##summary how many ints are common among all treatments
v=vector()
m=as.matrix(ints[,2:8])
shared=apply(m,1,FUN=function(x) length(x[x>=5]) )

png('plots/shared_interactions.png')
barplot(table(shared),xlab='shared between n treatments',ylab='# prom-interacting DpnII frags',main='summary shared interactions')
dev.off()

##replace the numbers by RA[number]
for (i in 2:8){names(ints)[i]=paste0('RA',names(ints)[i])}

##define clusters of interaction type

ints$clusters='na'
ints[ints$RA0h>=5 & ints$RA72h<4 ,]$clusters='ES_specific'
ints[ints$RA0h<4 & ints$RA72h>=5 ,]$clusters='RA_specific'
ints[ints$RA0h<4 & ints$RA72h<4 ,]$clusters='transient'
ints[ints$RA0h>=5 & ints$RA72h>=5 ,]$clusters='stable'

ints$clusters=factor(ints$clusters,levels = c('na','ES_specific','RA_specific','stable','transient'))

png('plots/Clusters_ints_frags.png')
tab=table(factor(ints$clusters,levels = c('na','ES_specific','RA_specific','stable','transient')))
barplot(tab,xlab='type interaction', ylab='# interactions', main='clusters timecourse DpnII frags')
dev.off()

png('plots/shared_interactions_byCluster.png')
par(mfrow=c(2,3))

for (cl in levels(ints$clusters)){
  v=vector()
  m=as.matrix(ints[ints$clusters==cl,2:8])
  shared=apply(m,1,FUN=function(x) length(x[x>=5]) )
  barplot(table(shared),xlab='shared between n treatments',ylab='# prom-interacting DpnII frags',main=paste('cluster:', cl))
}

dev.off()



png('plots/Boxplots_scores_all_significant_interactions.png',width=600,height=600)
par(mfrow=c(2,3))

boxplot(ints[,2:8],outline=FALSE, ylab='chiscores',col=brewer.pal(7,'PiYG'),ylim=c(0,14),main='all')
legend('topright',fill=brewer.pal(7,'PiYG'),legend=names(L),ncol=4,bty='n',title='RA treatment')
abline(h=5,lty=2)

for (cl in levels(ints$clusters)){
  boxplot(ints[ints$clusters==cl,2:8],outline=FALSE, ylab='chiscores',col=brewer.pal(7,'PiYG'),ylim=c(0,16),main=cl)
  legend('topright',fill=brewer.pal(7,'PiYG'),legend=names(L),ncol=4,bty='n',title='RA treatment')
  abline(h=5,lty=2)
}

dev.off()


##annotate with genename 
ints$genename='na'

for (id in unique(ints$baitID)){
  ints[ints$baitID==id,]$genename=names(baited_genes[baited_genes$dpn_id==id])
}

##annotate with induction status
ints$induced=FALSE
ints[ints$baitID %in% baited_genes[baited_genes$induced]$dpn_id,]$induced=TRUE

##make bedfiles for all possible clusters
for (cl in levels(ints$clusters)){
  bed=ints[ints$clusters==cl,11:13]
  bed$name=paste(ints[ints$clusters==cl,]$genename,1:nrow(bed),sep='_')
  saveBed(bed,paste0('ints.',cl,'.bed'))
}

write.table(ints,'ints_all.txt',sep='\t',eol='\n',quote=FALSE)


```


##Merge significant closeby interactions and assign the most prevalent cluster to them

```{r}

ints=read.delim('ints_all.txt',stringsAsFactors=FALSE)
# trans=ints[ints$seqnames_bait!=ints$seqnames_otherEnd,]
# ints=ints[ints$seqnames_bait==ints$seqnames_otherEnd,]

ids=unique(ints$baitID)

##distance fragments over which to merge (1 means two sig ints fragments with one unsig ints fragment inbetween will be merged)

aggregatePeaks <- function(ints, dis,plot=FALSE){
  ##aggregates significant interactions to peaks based on distance between peaks (1=no unsignificant fragment inbetween) and gives them a specific cluster
  ##Arguments:
  ##ints: ints table of interactions with columns: intID, seqnames_otherEnd, start_otherEnd, end_otherEnd, clusters
  ##dis: distance between significantly interacting fragments that are still considered one peak
  ##plot: Should a lineplot be made for each id?
  
  d=ints
  L_peaks=list()
  
  
  
  for (id in ids){
    print(id)
    
  
    a=d[d$baitID==id,]
    gr=GRanges(a$seqnames_otherEnd, IRanges(a$otherEndID,a$otherEndID))
    values(gr)=a
    gr=reduce(gr,min.gapwidth=dis)
    
    L_ints=list()
    for (i in 1:length(gr)){
      oe=start(gr[i]):end(gr[i]) 
      L_ints[[i]]=a[a$otherEndID %in% oe,]
    }
    
    ##aggregate to peaks, averaging everything
   
    
    for (i in 1:length(L_ints)){
      
      b=L_ints[[i]]
      c=b[1,]
      
      if (nrow(b)>1){
        oe=b[1,]$otherEndID
        for (j in 2:nrow(b)){oe=paste(oe,b[j,]$otherEndID,sep=',')}
        c$otherEndID=oe
      
        oe=b[1,]$intID
        for (j in 2:nrow(b)){oe=paste(oe,b[j,]$intID,sep=',')}
        c$intID=oe
      }
      
      c$start_otherEnd=min(b$start_otherEnd)
      c$end_otherEnd=max(b$end_otherEnd)
      
      for (j in 2:(length(types)+1)){
        c[,j]=mean(b[,j])
      }
      
      cl=names(table(b$clusters)[table(b$clusters)==max(table(b$clusters))])
      if (length(cl)>1){
        if (length(cl)==nrow(b)){
          c=b
        } else {
          for (j in 2:length(cl)){c=rbind(c,c)}
        }
      }
      c$clusters=cl
      
      ##Resolve ambiguous clusters:
  
      
      L_peaks=c(L_peaks,list(c))
      
    }
    
    
  }
  
  df=do.call('rbind',L_peaks)
  d=df
  # d=d[d$clusters %in% c('other')==FALSE,]
  d$clusters_refined=d$clusters
  
  if (plot==TRUE){
    png(paste0('Aggregated_peaks_dis',dis,'.png'),width=1000, height=800)
    par(mfrow=c(2,1))
    id=834429
    plotInteractions(L,id,k=5,zoom=1000000,unt=1,tam=2,show.legend = TRUE,ylim=c(0,0.06),d=d,xlim=c(146500000,147000000),intervals=grl)
    plotInteractions(L,id,k=5,zoom=1000000,unt=1,tam=2,show.legend = TRUE,ylim=c(0,0.06),d=d,xlim=c(148000000,148200000),intervals=grl)
    dev.off()
  }
  
  return(d)
}


for (dis in 1:10){
  print(dis)
  d=aggregatePeaks(ints,dis,plot=TRUE)
  write.table(d, paste('ints_all_aggregated_distance',dis,'.txt',sep=''),sep='\t',eol='\n',quote=FALSE)
}

```


##From merged peaks identify summits as locally maximum read count in the wt (if peaks present there) or the ko in which teh peaks are present):

(selected dis: 5)

```{r}

dis=5

dd=read.delim(paste('ints_all_aggregated_distance',dis,'.txt',sep=''),stringsAsFactors=FALSE)
dd$summit=0

size=2  ##peak+/- how many DpnII frags should be considered for summit?

d=dd[dd$clusters %in% c('Med13ko_specific')==FALSE,]

##For peaks that are present in UNT:
for (i in 1:nrow(d)){
  print(i)
  a=d[i,]
  oe=as.numeric(unlist(strsplit(a$otherEndID,split=',')))
  if (length(oe)>1){
    b=L[[1]]
    bb=b[b$baitID==a$baitID,]
    bb=bb[bb$otherEndID %in% (min(oe)-size):(max(oe)+size),]
    summit=bb[bb$N==max(bb$N),]$otherEndID
    
    if (length(summit)>1){
      
      if (length(summit[summit %in% oe])==1){
          summit=summit[summit %in% oe]
        } else if (length(summit[summit %in% oe])>0){
          summit=summit[summit %in% oe]
          summit=sample(summit,1)
        } else {
          summit=summit[summit %in% min(oe):max(oe)]
          if (length(summit)>1){
            summit=sample(summit,1)
          }
        }
    } 
    
    d[i,]$summit=summit
    
  } else {
    d[i,]$summit=oe
  }
}

peaks=d



##For Med13ko specific:
d=dd[dd$clusters %in% c('Med13ko_specific')==TRUE,]

##For peaks that are present in Tir1:
for (i in 1:nrow(d)){
  print(i)
  a=d[i,]
  oe=as.numeric(unlist(strsplit(a$otherEndID,split=',')))
  if (length(oe)>1){
    b=L[[2]]
    bb=b[b$baitID==a$baitID,]
    bb=bb[bb$otherEndID %in% (min(oe)-size):(max(oe)+size),]
    summit=bb[bb$N==max(bb$N),]$otherEndID
    
    if (length(summit)>1){
      
      if (length(summit[summit %in% oe])==1){
          summit=summit[summit %in% oe]
        } else if (length(summit[summit %in% oe])>0){
          summit=summit[summit %in% oe]
          summit=sample(summit,1)
        } else {
          summit=summit[summit %in% min(oe):max(oe)]
          if (length(summit)>1){
            summit=sample(summit,1)
          }
        }
    } 
    
    d[i,]$summit=summit
    
  } else {
    d[i,]$summit=oe
  }
}

peaks=rbind(peaks,d)




ov=findOverlapsDf(peaks,hind,'summit','id')
bed.summit=gr2bed(hind[ov[,2]])
names(bed.summit)=paste(names(bed.summit),'summit',sep='_')
names(bed.summit)[1]='seqnames_summit'
peaks=cbind(peaks,bed.summit)

grl2=GRangesList(bed2gr(bed.summit))
names(grl)='summits'
d=peaks[peaks$clusters != 'stable',]
d$clusters_refined=d$clusters
id=834429
plotInteractions(L,id,k=20,zoom=200000,unt=1,tam=2,show.legend = TRUE,ylim=c(0,0.1),d=d
                 
```
