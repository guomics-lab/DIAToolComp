fileRT<-list.files(path = '../../Matrix/',pattern = "RT",full.names = T)
fileRTS<-fileRT[grepl("_S_",fileRT)]
library(corrplot)
library(RColorBrewer)

MatToList<-function(MatIn){
  ListOut<- unlist(MatIn)
  names(ListOut)<-as.vector(sapply(1:ncol(MatIn),paste,row.names(MatIn),sep="_"))
  ListOut
}

corRTMats<-list()
for(condName in LETTERS[1:6]){
  fileRTSSel<-fileRTS[grepl(condName,fileRTS)]
  
  corRTlist<-list()
  for(m in 1:length(fileRTSSel)){
    readIn1<-read.table(fileRTSSel[m],stringsAsFactors = F,sep="\t",header = T,row.names = 1)
    corRTlist[[m]]<-MatToList(readIn1)
  }
  
  corRTMat<-matrix(NA,nrow=length(corRTlist),ncol=length(corRTlist))
  fileNames<-sapply(strsplit(fileRTSSel, "\\/"), function(x) x[5], simplify=T)
  names(corRTlist)<-row.names(corRTMat)<-colnames(corRTMat)<-paste(sapply(strsplit(fileNames, "_")
                                                                          , function(x) x[c(2)], simplify=T),
                                                    sapply(strsplit(fileNames, "_"), function(x) x[c(4)], simplify=T),sep="_")
  # RT cor mat
  for(i in 1:nrow(corRTMat)){
    for(j in 1:ncol(corRTMat)){
      list1<-corRTlist[[row.names(corRTMat)[i]]]
      # list1<-gsub("#N/A","NA",list1)
      list2<-corRTlist[[colnames(corRTMat)[j]]]
      # list2<-gsub("#N/A","NA",list2)
      interList<-intersect(names(list1),names(list2))
      
      list1I<-list1[interList]
      list2I<-list2[interList]
      
      corRTMat[j,i]<-corRTMat[i,j]<-cor(list1I,list2I,use="na.or.complete")
    }
  }
  
  corRTMats[[which(condName==condNames)]]<-corRTMat
  print(condName)
}

# pdf("Figure4A_corRTMat_FZ20210729.pdf",width=10,height=30)
# par(mfrow=c(3,2),mar=c(4,4,4,4))
# for(m in 1:6){
#   corRTMat<-corRTMats[[m]]
#   corRTMat[is.na(corRTMat)]<-0
#   #RT ABC
#   corrplot(corRTMat, order = "original", type = "lower", method = "ellipse", diag = FALSE,
#            tl.srt = 60,tl.col = "black",main=paste(condNames[m],"RT",sep="\t")
#            , cl.pos = "r", cl.lim = c(0,1), cl.length = 21
#            , bg = "grey", addCoef.col = "red", col=brewer.pal(n=10, name="RdYlBu"),na.label = "*")
#   
# }
# dev.off()

#RT difference
diffRTMats<-list()
for(condName in LETTERS[c(1:6)]){
  fileRTSSel<-fileRTS[grepl(condName,fileRTS)]
  
  diffRTlist<-list()
  for(m in 1:length(fileRTSSel)){
    readIn1<-read.table(fileRTSSel[m],stringsAsFactors = F,sep="\t",header = T,row.names = 1)
    diffRTlist[[m]]<-MatToList(readIn1)
  }
  
  diffRTMat<-matrix(NA,nrow=length(diffRTlist),ncol=length(diffRTlist))
  fileNames<-sapply(strsplit(fileRTSSel, "\\/"), function(x) x[5], simplify=T)
  names(diffRTlist)<-row.names(diffRTMat)<-colnames(diffRTMat)<-paste(sapply(strsplit(fileNames, "_")
                                                                          , function(x) x[c(2)], simplify=T),
                                                                   sapply(strsplit(fileNames, "_"), function(x) x[c(4)], simplify=T),sep="_")
  # RT cor mat
  pdf(paste0("Figure4_RTcordiff_",condName,"_FZ20210729.pdf"),width = 15,height = 15)
  par(mfrow=c((length(diffRTlist)),(length(diffRTlist))))
  for(i in 1:nrow(diffRTMat)){
    for(j in 1:i){
      if(i==j){
        plot.new()
      }else{
        list1<-diffRTlist[[row.names(diffRTMat)[i]]]
        # list1<-gsub("#N/A","NA",list1)
        list2<-diffRTlist[[colnames(diffRTMat)[j]]]
        # list2<-gsub("#N/A","NA",list2)
        interList<-intersect(names(list1),names(list2))
        
        list1I<-list1[interList]
        list2I<-list2[interList]
        r <- cor(list1I,list2I,use="na.or.complete")
        if(length(list1I)>0){
          smoothScatter(list1I, list2I, nrpoints = 100,cex = 2,
                        colramp = colorRampPalette(c(blues9,"yellow", "red")),
                        main = paste("r =",round(r,4)),xlab=row.names(diffRTMat)[i]
                        ,ylab=colnames(diffRTMat)[j],cex.main=2,cex.lab=2)
        }else{
          plot.new()
        }
      }
    }
    if(i<nrow(diffRTMat)){
      for(j in (i+1):ncol(diffRTMat)){
        list1<-diffRTlist[[row.names(diffRTMat)[i]]]
        # list1<-gsub("#N/A","NA",list1)
        list2<-diffRTlist[[colnames(diffRTMat)[j]]]
        # list2<-gsub("#N/A","NA",list2)
        interList<-intersect(names(list1),names(list2))
        
        list1I<-list1[interList]
        list2I<-list2[interList]
        boxplot(list1I-list2I,outline = F,main=paste(row.names(diffRTMat)[i],
                                                     colnames(diffRTMat)[j],sep="-")
                ,cex.main=2,cex.axis=2,cex.lab=2
                )
      }
    }
  }
  dev.off()
}
 