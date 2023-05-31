
params0<-read.delim('aa.txt',header = T,row.names = 1,stringsAsFactors = F)
params0<-params0[-nrow(params0),]

params1<-read.delim('growth.txt',header = T,row.names = 1,stringsAsFactors = F)
otu_table<-NULL
taxon_table<-NULL
tree<-NULL
export_path<-getwd()
export_path1<-paste(export_path,"/zhj",sep = "")
library(microchat)
library(tidyverse)
###训练样本
#params<-addMicrochatSample(params,maxnum=1.05,minnum=0.95,
#                           corrected=TRUE,
#                           robust=0.05,
#                           decim=2,
#                           addTo.sample.size=6,
#                           add.sample.size=6)
###生成microchat对象
mchat1<-setParamchat(otu_table,taxon_table,tree=tree,params=params1)

###筛选样本
tidymchat1<-tidyMicrochat(mchat1,
                          group.keep=NULL)
submchat1<-subsampleMicrochat(tidymchat1,sample.size=NULL,export_path=export_path)

###1.计算
microchatParamobj<-calcMicrochatParam(submchat1,
                                      export_path=export_path1)


dd.dat<-microchatParamobj$statistics
dd.dat<-subset(dd.dat,select = c(group,mean,index))

dd.dat.d<-spread(dd.dat, group, mean)%>%column_to_rownames(var = "index")%>%t()%>%data.frame()
colnames(params0)<-rownames(dd.dat.d)
ws<-cbind(t(params0),dd.dat.d)


mes1<-t(dd.dat.d)%>%t()
mat<-t(params0)
dim(mes1)
dim(mat)

corr_me_env<-psych::corr.test(mat,mes1,adjust="fdr",method = 'pearson')
resr<-corr_me_env$r%>%data.frame()
resp<-corr_me_env$p%>%data.frame()

resr1<-resr
resr1$group<-rownames(resr1)
cxcr<-reshape2::melt(resr1)
colnames(cxcr)[3]<-"r"
resp1<-resp
resp1$group<-rownames(resp1)
cxcp<-reshape2::melt(resp1)
colnames(cxcp)[3]<-"p.value"
com<-cbind(cxcr,cxcp)
com<-subset(com,select = c(2,1,3,6))
#if (length(unique(com$group))==1) com$group<-paste0(id.sel$trait,id.sel$module)[1]
com$p.value<-round(com$p.value,2)
com$sig<-ifelse(com$p.value<0.001,"***",
                ifelse(com$p.value<=0.01,"**",
                       ifelse(com$p.value<0.05,"*","")))

com$remark<-ifelse(com$r>0 & com$sig!="","up",
                   ifelse(com$r<0 & com$sig!="","down","not"))

sigcolor<-c("grey","white","red")
ggplot(com, mapping = aes(x=variable,y=group)) +
  geom_tile(linetype="dashed",
            fill="white",aes(color=r))+
  scale_color_gradient2(
    limits=c(-1,1),
    midpoint = 0,
    low = sigcolor[1],
    mid = "white",
    high = sigcolor[3],
    space = "Lab" ,n.breaks=6)+
  geom_tile(aes(height=r,width=r,fill = r))+
  scale_fill_gradient2(
    limits=c(-1,1),
    midpoint = 0,
    low = sigcolor[1],
    mid = "white",
    high = sigcolor[3],
    space = "Lab" ,n.breaks=6)+
  geom_text(aes(label=sig),vjust=0.75,color="black",angle=45,size=5)+
  theme(text = element_text(family = "serif"),
        title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        #plot.margin = margin(1, 1, 1, 1, "cm"),
        #panel.background = element_rect(fill = "white",color="white"),
        panel.grid = element_line(color=NA),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle=270,vjust = 0.5),
        aspect.ratio = length(unique(com$group))/length(unique(com$variable)))


