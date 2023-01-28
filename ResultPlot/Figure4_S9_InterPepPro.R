library(pheatmap)
pchLabel<-c(1,20)
names(pchLabel)<-c("free","lib")
colorLabel<-c("#2697c1","#eb9396","#a5c48d","#2ab58d","#d4c41d")
names(colorLabel)<-c("1eny","2dnn","3osw","4sky","5spn")
HInd<-c(1.8,-4.5,-3.5,-3.5,2.5,-3.5,-3.5,-0.4,-3.2,4.5,3.8,-3.9,1.9,2.8,-1.6,-0.8,-0.7,-0.9,-1.3,4.2)
names(HInd)<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
IntersectList <- function(listListSel){
  intAll0<-listListSel[[1]]
  if(length(listListSel)>1){
    for(i in 2:(length(listListSel))){
      intAll0<-intersect(intAll0,listListSel[[i]])
    }
  }
  intAll0
}

GravyCal<-function(seqin){
  seq_sp<-HInd[unlist(strsplit(seqin,""))]
  round(sum(seq_sp)/length(seq_sp),2)
}

PepLength<-function(seqin){
  seq_sp<-unlist(strsplit(seqin,""))
  length(seq_sp)
}

MisCleave<-function(seqin){
  seq_sp<-unlist(strsplit(seqin,""))
  sum(seq_sp%in%c("K","R"))
}

VennOd90<-function(listList, ylab=""){
  nList<-length(listList)
  df_vennAll<-data.frame(matrix(NA,nrow=nList,ncol=2^nList-1))
  m<-1
  for(i in nList:1){
    comMat<-combn(nList,i)
    for(j in 1:ncol(comMat)){
      df_vennAll[comMat[,j],m]<-1
      m<-m+1
    }
  }
  row.names(df_vennAll)<-names(listList)
  AllElement<-unique(unlist(listList))
  if(nList>0){
    peptAll<-list()
    for(j in 1:ncol(df_vennAll)){
      indSel<-which(!is.na(df_vennAll[,j]))
      indUnSel<-which(is.na(df_vennAll[,j]))
      if(j>1){
        EleSel<-AllElement[AllElement%in%(IntersectList(listList[indSel]))&(!AllElement%in%(unlist((listList[indUnSel]))))]
        peptAll[[j]]<-EleSel
      }else {
        peptAll[[j]]<-IntersectList(listList[indSel])
      }
    }
  }
  numAll<-sapply(peptAll,length)
  ind0<-which(numAll==0)
  od<-order(numAll,decreasing = T)
  od<-od[!od%in%ind0]
  ind90<-cumsum(numAll[od])/max(cumsum(numAll[od]))<0.90
  peptAllod90<-peptAll[od][ind90]
  
  if(nList>1){
    df_vennAllOd<-df_vennAll[,od][,ind90]
  }
  df_vennAllOd[is.na(df_vennAllOd)]<-0
  names(peptAllod90)<-paste0("X",apply(df_vennAllOd,2,paste0,collapse=""))
  return(peptAllod90)
}

vennMap<-function(xnames,rowVenn){
  mapNames<-gsub("X","",xnames)
  strsplit(mapNames,"")
  dfPlot<-as.data.frame(strsplit(mapNames,""))
  dfPlot2<-mapply(dfPlot, FUN=as.numeric)
  row.names(dfPlot2)<-rowVenn
  pheatmap(dfPlot2,show_rownames=T,show_colnames=F,cluster_rows=F,cluster_cols = F
           ,color=c("#5254AD","#FA6216"),border_color = "white",legend=F
           )
}

# pdf("PepProInt90BD_20210919.pdf",width = 4,height = 4)
  
# par(mfrow=c(2,3))
#pepList
tsvAllPep<-list.files("../../Matrix/",pattern = "PepInt.tsv")
tsvAllPro<-list.files("../../Matrix/",pattern = "ProInt.tsv")

# pdf(paste0("PepIntCount.pdf"),width =12,height = 6)
# pdf(paste0("PepIntMisClv.pdf"),width =12,height = 6)
# par(mfrow=c(2,3),mar=c(2,4,4,4))
for(DataID in LETTERS[1]){#LETTERS[1:6]
  # pdf(paste0("PepIntMap",DataID,".pdf"),width = 6,height = 6)
  tsvSel<-tsvAllPep[grepl(DataID,tsvAllPep)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  tsvNames<-paste(sapply(strsplit(tsvSel,"_"),"[[",2),sapply(strsplit(tsvSel,"_"),"[[",4)
         ,sep="_")
  ratSel<-gsub("Int","RatioT",tsvSel)
  
  pepList<-list()#pepListT<-
  for(i in 1:length(tsvSel)){
    pepIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F)
    pepOut<-as.data.frame(pepIn[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn)),]
                          ,StringAsFactors=F)
    # ratIn<-read.table(paste0("../../Matrix/",ratSel[i]),stringsAsFactors = F,header = T,row.names = NULL)
    colnames(pepOut)<-colnames(pepIn)
    row.names(pepOut)<-row.names(pepIn)[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn))]
    pepList[[i]]<-row.names(pepOut)
    # pepListT[[i]]<-row.names(pepOut)[row.names(pepOut)%in%ratIn$row.names[which(ratIn$TF)]]
  }
  names(pepList)<-tsvNames
  pepVenn90<-VennOd90(pepList,paste0("pep",DataID))
  # meanClv<-mean(sapply(pepVenn90[[1]],MisCleave))-1
  # barplot(sapply(pepVenn90,length),horiz = T,col = 0,las=2,main=DataID)
  # vennMap(names(pepVenn90),tsvNames)
  # dev.off()
  meanClv<-c()
  for(i in 1:length(pepVenn90)){
    meanClv[i]<-mean(sapply(pepVenn90[[i]],MisCleave))-1
  }
  pdf(paste0("PepIntMis2",DataID,".pdf"),width = 6,height = 6)
  pheatmap(as.matrix(meanClv),cluster_rows=F,cluster_cols=F,color =0
           ,legend = T, display_numbers = TRUE,fontsize = 12)
  dev.off()
  
  #Gravy####
  # df0<-data.frame(vennName=paste0((LETTERS)[1],names(pepVenn90)[1]),GravyScore=sapply(pepVenn90[[1]],GravyCal))
  # for(i in 2:length(pepVenn90)){
  #   df1<-data.frame(vennName=paste0((LETTERS)[i],names(pepVenn90)[i]),GravyScore=sapply(pepVenn90[[i]],GravyCal))
  #   df0<-rbind(df0,df1)
  # }
  # 
  # library(ggridges)
  # pdf(paste0("PepGravy",DataID,".pdf"),width = 6,height = 12)
  # ggplot(df0, aes(x = GravyScore, y = vennName,fill=vennName)) +
  #   geom_density_ridges() +
  #   theme_ridges() + 
  #   theme(legend.position = "none")
  # dev.off()
  
  #PepLength####
  # df0<-data.frame(vennName=paste0((LETTERS)[1],names(pepVenn90)[1]),PepLen=sapply(pepVenn90[[1]],PepLength))
  # for(i in 2:length(pepVenn90)){
  #   df1<-data.frame(vennName=paste0((LETTERS)[i],names(pepVenn90)[i]),PepLen=sapply(pepVenn90[[i]],PepLength))
  #   df0<-rbind(df0,df1)
  # }
  # ggplot(df0, aes(x = PepLen, y = vennName,fill=vennName)) +
  #   geom_density_ridges() +
  #   theme_ridges() +
  #   theme(legend.position = "none")
  # ggsave(paste0("PepLength",DataID,".pdf"))
}
# dev.off()

#retention time
for(DataID in LETTERS[5:6]){#LETTERS[1:6]
  tsvSel<-tsvAllPep[grepl(DataID,tsvAllPep)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  tsvNames<-paste(sapply(strsplit(tsvSel,"_"),"[[",2),sapply(strsplit(tsvSel,"_"),"[[",4)
                  ,sep="_")
  ratSel<-gsub("Int","RatioT",tsvSel)
  libIn<-read.delim(list.files(path=paste0("/Volumes/ToolComp/ToolComp/Library/lib",DataID,"/5spn"),pattern = ".xls$",full.names = T)[1],stringsAsFactors = F)
libRT<-unique(libIn[,c("StrippedPeptide","iRT")])
  libRTVal<-libRT$iRT
  libRTVal<-as.numeric(gsub(",",".",libRTVal))
  names(libRTVal)<-libRT$StrippedPeptide
#F
  libRTVal2<-as.numeric(gsub(",",".",libRTVal))
  names(libRTVal2)<-libRT$StrippedPeptide
  libRTVal<-libRTVal2
  
  pepList<-list()
  for(i in 1:length(tsvSel)){
    pepIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F)
    pepOut<-as.data.frame(pepIn[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn)),]
                          ,StringAsFactors=F)
    colnames(pepOut)<-colnames(pepIn)
    row.names(pepOut)<-row.names(pepIn)[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn))]
    pepList[[i]]<-row.names(pepOut)
  }
  names(pepList)<-tsvNames
  pepVenn90<-VennOd90(pepList,paste0("pep",DataID))
  #RT####
  #1-5
  # df0<-data.frame(vennName=paste0((LETTERS)[1],names(pepVenn90)[1]),RT=libRTVal[pepVenn90[[1]]])
  # for(i in 2:length(pepVenn90)){
  #   df1<-data.frame(vennName=paste0((LETTERS)[i],names(pepVenn90)[i]),RT=libRTVal[pepVenn90[[i]]])
  #   df0<-rbind(df0,df1)
  # }
  # df0<-df0[!is.na(df0$RT),]
  #6
  df0<-data.frame(vennName=paste0((LETTERS)[1],names(pepVenn90)[1]),RT=libRTVal[pepVenn90[[1]]])
  for(i in 2:length(pepVenn90)){
    df1<-data.frame(vennName=paste0((LETTERS)[i],names(pepVenn90)[i]),RT=libRTVal[pepVenn90[[i]]])
    df0<-rbind(df0,df1)
  }
  df0<-df0[!is.na(df0$RT),]
  
  
  ggplot(df0, aes(x = RT, y = vennName,fill=vennName)) +
    geom_density_ridges() +
    theme_ridges() +
    theme(legend.position = "none")
  ggsave(paste0("PepIntRT",DataID,".pdf"))
}

# Intensity
for(DataID in LETTERS[6]){#LETTERS[1:6]
  tsvSel<-tsvAllPep[grepl(DataID,tsvAllPep)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  tsvNames<-paste(sapply(strsplit(tsvSel,"_"),"[[",2),sapply(strsplit(tsvSel,"_"),"[[",4)
                  ,sep="_")
  pepList<-list()
  pepInList<-list()
  for(i in 1:length(tsvSel)){
    pepIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F)
    pepOut<-as.data.frame(pepIn[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn)),]
                          ,StringAsFactors=F)
    colnames(pepOut)<-colnames(pepIn)
    row.names(pepOut)<-row.names(pepIn)[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn))]
    pepList[[i]]<-row.names(pepOut)
    pepInList[[i]]<-pepIn
  }
  names(pepList)<-tsvNames
  pepVenn90<-VennOd90(pepList,paste0("pep",DataID))
  
  for(j in 1:length(tsvNames)){
    df0<-data.frame(vennName=paste0((LETTERS)[1],names(pepVenn90)[1])
                    ,Intensity=log10(pepInList[[j]][pepVenn90[[1]],]))
    for(i in 2:length(pepVenn90)){
      df1<-data.frame(vennName=paste0((LETTERS)[i],names(pepVenn90)[i])
                      ,Intensity=log10(pepInList[[j]][pepVenn90[[i]],]))
      df0<-rbind(df0,df1)
    }
    ggplot(df0, aes(x = Intensity, y = vennName,fill=vennName)) +
      geom_density_ridges() +
      theme_ridges() +
      theme(legend.position = "none")
    # ggsave(paste0("PepIntensity",DataID,"-",tsvNames[j],".pdf"))
  }
}

for(DataID in LETTERS[1:5]){#LETTERS[1:6]
  tsvSel<-tsvAllPep[grepl(DataID,tsvAllPep)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  tsvNames<-paste(sapply(strsplit(tsvSel,"_"),"[[",2),sapply(strsplit(tsvSel,"_"),"[[",4)
                  ,sep="_")
  pepList<-list()
  pepInList<-list()
  for(i in 1:length(tsvSel)){
    pepIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F)
    pepOut<-as.data.frame(pepIn[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn)),]
                          ,StringAsFactors=F)
    colnames(pepOut)<-colnames(pepIn)
    row.names(pepOut)<-row.names(pepIn)[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn))]
    pepList[[i]]<-row.names(pepOut)
    pepInList[[i]]<-pepIn
  }
  names(pepList)<-tsvNames
  pepVenn90<-VennOd90(pepList,paste0("pep",DataID))
  
  for(j in 1:length(tsvNames)){
    df0<-data.frame(vennName=paste0((LETTERS)[1],names(pepVenn90)[1])
                    ,Intensity=log10(apply(pepInList[[j]][pepVenn90[[1]],],1,mean,na.rm=T)))
    for(i in 2:length(pepVenn90)){
      df1<-data.frame(vennName=paste0((LETTERS)[i],names(pepVenn90)[i])
                      ,Intensity=log10(apply(pepInList[[j]][pepVenn90[[i]],],1,mean,na.rm=T)))
      df0<-rbind(df0,df1)
    }
    df0<-df0[df0$Intensity!="NaN",]
    df0<-df0[df0$Intensity!="-Inf",]
    
    ggplot(df0, aes(x = Intensity, y = vennName,fill=vennName)) +
      geom_density_ridges() +
      theme_ridges() +
      theme(legend.position = "none")
    ggsave(paste0("PepIntensity",DataID,"-",tsvNames[j],".pdf"))
  }
  
}

# pie chart####
for(DataID in LETTERS[1:6]){#LETTERS[1:6]
  tsvSel<-tsvAllPep[grepl(DataID,tsvAllPep)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  tsvNames<-paste(sapply(strsplit(tsvSel,"_"),"[[",2),sapply(strsplit(tsvSel,"_"),"[[",4)
                  ,sep="_")
  pepList<-list()
  for(i in 1:length(tsvSel)){
    pepIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F)
    pepOut<-as.data.frame(pepIn[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn)),]
                          ,StringAsFactors=F)
    colnames(pepOut)<-colnames(pepIn)
    row.names(pepOut)<-row.names(pepIn)[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn))]
    pepList[[i]]<-row.names(pepOut)
  }
  names(pepList)<-tsvNames
  pepVenn90<-VennOd90(pepList,paste0("pep",DataID))
  
  pepFull<-unique(unlist(pepList))
  cumPies<-cumsum(sapply(pepVenn90,length))/length(pepFull)
  
  pdf(paste0("PepPie",DataID,".pdf"),width = 6,height = 12)
  par(mfrow=c(length(cumPies),1),mar=c(1,1,1,1))
  for(i in 1:length(cumPies)){
    pie(c(1-cumPies[i],cumPies[i]),labels=NA)
  }
  dev.off()
}

# peptide per protein####
#lib based
tsvAllPep<-list.files("../../Matrix/",pattern = "PepInt.tsv")

for(DataID in LETTERS[1:6]){
  tsvSel<-tsvAllPep[grepl(DataID,tsvAllPep)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  tsvSel<-tsvSel[grepl("lib",tsvSel)]
  # 
  tsvNames<-paste(sapply(strsplit(tsvSel,"_"),"[[",2),sapply(strsplit(tsvSel,"_"),"[[",4)
                  ,sep="_")
  
  for(i in 1:length(tsvSel)){
    pepIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F)
    pepOut<-as.data.frame(pepIn[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn)),]
                          ,StringAsFactors=F)
    colnames(pepOut)<-colnames(pepIn)
    row.names(pepOut)<-row.names(pepIn)[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn))]
    
    proIn<-read.table(paste0("../../Matrix/",gsub("Pep","Pro",tsvSel[i])),stringsAsFactors = F,sep="\t",header = T,row.names = 1)
    proOut<-as.data.frame(proIn[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn))),],StringAsFactors=F)
    colnames(proOut)<-colnames(proIn)
    row.names(proOut)<-row.names(proIn)[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn)))]
    
    libIn<-read.delim(list.files(path=paste0("/Volumes/ToolComp/ToolComp/Library/lib",DataID,"/5spn"),pattern = ".xls$",full.names = T)[1],stringsAsFactors = F)
    libPepPro<-unique(libIn[,c("StrippedPeptide","UniProtIds")])
    
    pepNameLib<-libPepPro$UniProtIds
    names(pepNameLib)<-libPepPro$StrippedPeptide
    
    pepPro<-pepNameLib[row.names(pepOut)]
    pepProL<-pepPro[!is.na(pepPro)]
    
    proPepName<-unique(pepProL)
    proPep<-c()
    for(m in 1:length(proPepName)){
      proPep[m]<-paste(names(pepProL)[which(pepProL==proPepName[m])],collapse = ";")
    }
    names(proPep)<-proPepName
   
    proPepMat<-proPep[row.names(proOut)]
    proPepMat2<-proPepMat[!is.na(proPepMat)]
    
    #
    tsvOut<-gsub("PepInt","PepPro",tsvSel[i])
    print(tsvOut)
    write.csv(proPepMat2,tsvOut)
  }
}

#lib-free
# library(seqinr)
# tsvAllPep<-list.files("../../Matrix/",pattern = "PepInt.tsv")
# for(DataID in LETTERS[1:6]){
#   tsvSel<-tsvAllPep[grepl(DataID,tsvAllPep)]
#   tsvSel<-tsvSel[grepl("S",tsvSel)]
#   tsvSel<-tsvSel[grepl("free",tsvSel)]
#   tsvNames<-paste(sapply(strsplit(tsvSel,"_"),"[[",2),sapply(strsplit(tsvSel,"_"),"[[",4)
#                   ,sep="_")
#   fastaIn<-list.files(path=paste0("/Volumes/ToolComp/ToolComp/Library/lib",DataID)
#                       ,pattern = ".fasta",full.names = T)
#   libFasta<-read.fasta(fastaIn[1],seqtype = "AA",as.string = T)
#   fastaAll<-unlist(libFasta)
# 
#   for(i in 1:length(tsvSel)){
#     pepIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F)
#     pepOut<-as.data.frame(pepIn[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn)),]
#                           ,StringAsFactors=F)
#     colnames(pepOut)<-colnames(pepIn)
#     row.names(pepOut)<-row.names(pepIn)[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn))]
# 
#     pepPro<-c()
#     for(j in 1:length(row.names(pepOut))){
#       indSel<-which(grepl(row.names(pepOut)[j],fastaAll))
#       if(length(indSel)>0){
#         pepPro[j]<-paste(names(fastaAll)[indSel],collapse = ";")
#       }else{
#         pepPro[j]<-NA
#       }
#       if(j%%1000==0) print(j)
#     }
#     print(tsvSel[i])
# 
#     names(pepPro)<-row.names(pepOut)
#     pepProL<-pepPro[!is.na(pepPro)]
# 
#     proPepName<-unique(pepProL)
#     proPep<-c()
#     for(m in 1:length(proPepName)){
#       proPep[m]<-paste(names(pepProL)[which(pepProL==proPepName[m])],collapse = ";")
#     }
#     names(proPep)<-proPepName
# 
#     tsvOut<-gsub("PepInt","PepPro",tsvSel[i])
#     write.csv(proPep,tsvOut)
#   }
# }

#ProIntensity
library(ggridges)
library(ggplot2)
tsvAllPro<-list.files("../../Matrix/",pattern = "ProInt.tsv")

for(DataID in LETTERS[1:6]){#
  tsvSel<-tsvAllPro[grepl(DataID,tsvAllPro)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  tsvNames<-paste(sapply(strsplit(tsvSel,"_"),"[[",2),sapply(strsplit(tsvSel,"_"),"[[",4)
                  ,sep="_")
  proList<-list()
  proInList<-list()
  for(i in 1:length(tsvSel)){
    proIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F,sep="\t",header = T,row.names = 1)
    proOut<-as.data.frame(proIn[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn))),]
                          ,StringAsFactors=F)
    colnames(proOut)<-colnames(proIn)
    row.names(proOut)<-row.names(proIn)[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn)))]
    proList[[i]]<-row.names(proOut)
    proInList[[i]]<-proIn
  }
  names(proInList)<-names(proList)<-tsvNames
  proVenn90<-VennOd90(proList,paste0("pro",DataID))

  for(j in 1:length(tsvNames)){
    df0<-data.frame(vennName=paste0((LETTERS)[1],names(proVenn90)[1])
                    # ,Intensity=log10(apply(proInList[[j]][proVenn90[[1]],],1,mean,na.rm=T)))
                    ,Intensity=log10((proInList[[j]][proVenn90[[1]],])))
    for(i in 2:length(proVenn90)){
      df1<-data.frame(vennName=paste0((LETTERS)[i],names(proVenn90)[i])
                      # ,Intensity=log10(apply(proInList[[j]][proVenn90[[i]],],1,mean,na.rm=T)))
                      ,Intensity=log10((proInList[[j]][proVenn90[[i]],])))
      df0<-rbind(df0,df1)
    }
    ggplot(df0, aes(x = Intensity, y = vennName,fill=vennName)) +
      geom_density_ridges() +
      theme_ridges() +
      theme(legend.position = "none")
    ggsave(paste0("ProIntensity3",DataID,"-",tsvNames[j],".pdf"))
  }
}

# peptide per protein venn####
pepProFiles<-list.files(pattern = "_PepPro.tsv")
for(DataID in LETTERS[1:6]){
  tsvSel<-tsvAllPro[grepl(DataID,tsvAllPro)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  tsvNames<-paste(sapply(strsplit(tsvSel,"_"),"[[",2),sapply(strsplit(tsvSel,"_"),"[[",4)
                  ,sep="_")
  pepProList<-proList<-list()
  for(i in 1:length(tsvSel)){
    proIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F,sep="\t",header = T,row.names = 1)
    if(grepl("lib",tsvSel[i])){
      pepProIn<-read.csv(gsub("ProInt","PepPro",tsvSel[i]),stringsAsFactors = F,row.names = 1)
    }else{
      pepProIn<-read.table(gsub("ProInt","PepPro",tsvSel[i]),stringsAsFactors = F,header = T)
    }
    
    proOut<-as.data.frame(proIn[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn))),]
                          ,StringAsFactors=F)
    colnames(proOut)<-colnames(proIn)
    row.names(proOut)<-row.names(proIn)[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn)))]
    proList[[i]]<-row.names(proOut)
    pepProList[[i]]<-pepProIn
  }
  names(pepProList)<-names(proList)<-tsvNames
  proVenn90<-VennOd90(proList,paste0("pro",DataID))
  
  for(j in 1:length(tsvNames)){
    pepProCount<-sapply(pepProList[[j]][,1],countPep)
    names(pepProCount)<-row.names(pepProList[[j]])
    
    pepProCountSel<-pepProCount[proVenn90[[1]]]
    pepProCountSel<-pepProCountSel[!is.na(pepProCountSel)]
    df0<-data.frame(vennName=paste0((LETTERS)[1],names(proVenn90)[1])
                    ,PepPerPro=pepProCountSel)
    
    for(i in 2:length(proVenn90)){
      pepProCountSel<-pepProCount[proVenn90[[i]]]
      pepProCountSel<-pepProCountSel[!is.na(pepProCountSel)]
      if(length(pepProCountSel)>0){
        df1<-data.frame(vennName=paste0((LETTERS)[i],names(proVenn90)[i])
                        ,PepPerPro=pepProCountSel)
      }else{
        df1<-data.frame(vennName=paste0((LETTERS)[i],names(proVenn90)[i])
                        ,PepPerPro=0)
        print(i)
      }
      df0<-rbind(df0,df1)
    }
    ggplot(df0, aes(x = PepPerPro, y = vennName,fill=vennName)) +
      geom_density_ridges() +
      theme_ridges() +  xlim(0, 50)+
      theme(legend.position = "none")
    ggsave(paste0("PepPerPro2",DataID,"-",tsvNames[j],".pdf"))
  }
}

#MW count####
MWs<-c(71.03711,156.10111,114.04293,115.02694,103.00919,129.04259,128.05858,57.02146,137.05891,113.08406,113.08406,128.09496,131.04049,147.06841,97.05276,87.03203,101.04768,186.07931,163.06333,99.06841,150.953636,237.147727)
names(MWs)<-c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","U","O") 
MWcal<-function(seqin){
  pep_sp1<-unlist(strsplit(seqin," "))
  if(length(pep_sp1)>1){
    pep_sp2<-unlist(strsplit(pep_sp1[length(pep_sp1)],""))
    pep_sp<-pep_sp2[which(pep_sp2=="M")[1]:length(pep_sp2)]
  }else{
    pep_sp<-unlist(strsplit(pep_sp1[length(pep_sp1)],""))
  }
  sum(MWs[pep_sp])+18.01056
}
#fastaLib MW cal
for(DataID in LETTERS[6]){#F only
    fastaIn<-list.files(path=paste0("/Volumes/ToolComp/ToolComp/Library/lib",DataID)
                        ,pattern = ".fasta$",full.names = T)
    libFasta<-read.fasta(fastaIn[1],seqtype = "AA",as.string = T,whole.header = T)
    fastaAll<-unlist(libFasta)
    mwAll<-sapply(fastaAll,MWcal)
    protNameAll<-sapply(sapply(names(fastaAll),strsplit,"\\|"),"[",2)
    protNameAll2<-sapply(sapply(protNameAll,strsplit," "),"[",1)
    protNameAll3<-protNameAll2
    for(i in 1:length(protNameAll2)){
      if(grepl(":",protNameAll2[i])){
        protNameAll3[i]<-unlist(strsplit(protNameAll2[i],":"))[2]
      }
    }
    names(mwAll)<-protNameAll3
    mwOut<-mwAll[!is.na(names(mwAll))]
    mwOut2<-mwOut[which(!duplicated(names(mwOut)))]
    write.csv(mwOut2,paste0(DataID,"_MW.csv"))
}

library(ggplot2)
library(ggridges)
tsvAllPro<-list.files("../../Matrix/",pattern = "ProInt.tsv")

for(DataID in LETTERS[1:6]){
  tsvSel<-tsvAllPro[grepl(DataID,tsvAllPro)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  tsvNames<-paste(sapply(strsplit(tsvSel,"_"),"[[",2),sapply(strsplit(tsvSel,"_"),"[[",4)
                  ,sep="_")
  proList<-list()
  MWin<-read.csv(paste0(DataID,"_MW.csv"),stringsAsFactors = F)
  MWlib<-MWin[,2]
  names(MWlib)<-MWin[,1]
  MWlib<-MWlib[!is.na(names(MWlib))]
  
  for(i in 1:length(tsvSel)){
    proIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F,sep="\t",header = T,row.names = 1)
    proOut<-as.data.frame(proIn[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn))),]
                          ,StringAsFactors=F)
    colnames(proOut)<-colnames(proIn)
    row.names(proOut)<-row.names(proIn)[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn)))]
    proList[[i]]<-row.names(proOut)
  }
  proVenn90<-VennOd90(proList,paste0("pro",DataID))
  
  df0<-data.frame(vennName=paste0((LETTERS)[1],names(proVenn90)[1]),MW=MWlib[proVenn90[[1]]])
  for(i in 2:length(proVenn90)){
    newMW<-MWlib[proVenn90[[i]]]
    newMW<-newMW[!is.na(newMW)]
    df1<-data.frame(vennName=paste0((LETTERS)[i],names(proVenn90)[i]),MW=newMW)
    df0<-rbind(df0,df1)
  }
  df0<-df0[!is.na(df0$MW),]
  
  ggplot(df0, aes(x = MW, y = vennName,fill=vennName)) +
    geom_density_ridges() +
    theme_ridges() +  xlim(0, 500000)+
    theme(legend.position = "none")
  ggsave(paste0("mwPro_",DataID,".pdf"))
}


countPep<-function(pepsin){
  sum(unlist(strsplit(pepsin,""))==";")+1
}

for(DataID in LETTERS[1:6]){#LETTERS[1:6]
  tsvSel<-tsvAllPro[grepl(DataID,tsvAllPro)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  tsvNames<-paste(sapply(strsplit(tsvSel,"_"),"[[",2),sapply(strsplit(tsvSel,"_"),"[[",4)
                  ,sep="_")
  proList<-list()
  for(i in 1:length(tsvSel)){
    proIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F,sep="\t",header = T,row.names = 1)
    proOut<-as.data.frame(proIn[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn))),]
                          ,StringAsFactors=F)
    colnames(proOut)<-colnames(proIn)
    row.names(proOut)<-row.names(proIn)[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn)))]
    proList[[i]]<-row.names(proOut)
  }
  names(proList)<-tsvNames
  proVenn90<-VennOd90(proList,paste0("pro",DataID))
  
  proFull<-unique(unlist(proList))
  cumPies<-cumsum(sapply(proVenn90,length))/length(proFull)
  
  pdf(paste0("ProPie2",DataID,".pdf"),width = 6,height = 12)
  par(mfrow=c(length(cumPies),1),mar=c(1,1,1,1))
  for(i in 1:length(cumPies)){
    pie(c(1-cumPies[i],cumPies[i]),labels=NA)
  }
  dev.off()
}

# pdf(paste0("ProIntCount.pdf"),width =12,height = 6)
# par(mfrow=c(2,3),mar=c(2,4,4,4))
for(DataID in LETTERS[1:6]){#LETTERS[1:6]
  tsvSel<-tsvAllPro[grepl(DataID,tsvAllPro)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  tsvNames<-paste(sapply(strsplit(tsvSel,"_"),"[[",2),sapply(strsplit(tsvSel,"_"),"[[",4)
                  ,sep="_")
  ratSel<-gsub("Int","RatioT",tsvSel)
  
  proList<-list()#proListT<-
  for(i in 1:length(tsvSel)){
    proIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F,sep="\t",header = T,row.names = 1)
    proOut<-as.data.frame(proIn[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn))),]
                          ,StringAsFactors=F)
    colnames(proOut)<-colnames(proIn)
    row.names(proOut)<-row.names(proIn)[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn)))]
    proList[[i]]<-row.names(proOut)
  }
  names(proList)<-tsvNames
  proVenn90<-VennOd90(proList,paste0("pro",DataID))
  # barplot(sapply(proVenn90,length),horiz = T,col = 0,las=2,main=DataID)
  
  pdf(paste0("ProIntMap2",DataID,".pdf"))
  vennMap(names(proVenn90),tsvNames)
  dev.off()
}
# dev.off()

#MW vs Intensity Protein
for(DataID in LETTERS[1:5]){
  tsvSel<-tsvAllPro[grepl(DataID,tsvAllPro)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  tsvNames<-paste(sapply(strsplit(tsvSel,"_"),"[[",2),sapply(strsplit(tsvSel,"_"),"[[",4)
                  ,sep="_")
  MWin<-read.csv(paste0(DataID,"_MW.csv"),stringsAsFactors = F)
  MWlib<-MWin[,2]
  names(MWlib)<-MWin[,1]
  MWlib<-MWlib[!is.na(names(MWlib))]
  # proIntList<-list()
  for(i in 1:length(tsvSel)){
    proIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F,sep="\t",header = T,row.names = 1)
    proOut<-as.data.frame(proIn[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn))),]
                          ,StringAsFactors=F)
    colnames(proOut)<-colnames(proIn)
    row.names(proOut)<-row.names(proIn)[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn)))]
    
    proInt<-(apply(proOut,1,mean,na.rm=T))
    MWall<-MWlib[names(proInt)]
    plot((proInt),(MWall))
    print(cor((proInt),(MWall), use = "na.or.complete"))
  }
}
