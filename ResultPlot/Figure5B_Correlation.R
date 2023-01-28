library(corrplot)
library(RColorBrewer)
tsvAllPep<-list.files(pattern = "PepInt.tsv" ,path="../../Matrix/")
tsvAllPro<-list.files(pattern = "ProInt.tsv",path="../../Matrix/")

pep_files <-tsvAllPep[grepl("S",tsvAllPep)]
pro_files <-tsvAllPro[grepl("S",tsvAllPro)]

MatToList<-function(MatIn){
  ListOut<- unlist(MatIn)
  names(ListOut)<-as.vector(sapply(paste0("X",1:ncol(MatIn)),paste,row.names(MatIn),sep="_"))
  ListOut
}

condNames<-LETTERS[1:6]
corPepProMats<-list()
for(condName in condNames[6]){
  filePep<-pep_files[grepl(condName,pep_files)]
  filePro<-pro_files[grepl(condName,pro_files)]
  
  corPeplist<-corProlist<-list()
  for(m in 1:length(filePep)){
    readIn1 <-read.delim(paste0("../../Matrix/",filePep[m]),header =T,sep="\t",stringsAsFactors = F,row.names = 1)
    readIn1T<-read.delim(paste0("../../Matrix/",gsub("Int","RatioT",filePep[m])),header =T,sep="\t",stringsAsFactors = F,row.names = NULL)
    readIn1Sel<-readIn1[row.names(readIn1)%in%readIn1T$row.names[readIn1T$TF],]

    readIn2<-read.delim(paste0("../../Matrix/",filePro[m]),header =T,sep="\t",stringsAsFactors = F,row.names = 1)
    readIn2T<-read.delim(paste0("../../Matrix/",gsub("Int","RatioT",filePro[m])),header =T,sep="\t",stringsAsFactors = F,row.names = NULL)
    readIn2Sel<-readIn2[row.names(readIn2)%in%readIn2T$row.names[readIn2T$TF],]
    if(condName=="F"){#one column
      names(readIn1Sel)<-row.names(readIn1)[row.names(readIn1)%in%readIn1T$row.names[readIn1T$TF]]
      names(readIn2Sel)<-row.names(readIn2)[row.names(readIn2)%in%readIn2T$row.names[readIn2T$TF]]
      
      corPeplist[[m]]<-readIn1Sel
      corProlist[[m]]<-readIn2Sel
    }else{
      corPeplist[[m]]<-MatToList(readIn1Sel)
      corProlist[[m]]<-MatToList(readIn2Sel)
    }
  }
  
  corPepProMat<-matrix(NA,nrow=length(corProlist),ncol=length(corProlist))
  names(corPeplist)<-names(corProlist)<-row.names(corPepProMat)<-colnames(corPepProMat)<-sapply(strsplit(filePro, "_"), function(x)  paste(x[2],x[4],sep = "_"), simplify=T)
  
  #pep cormat
  for(i in 2:nrow(corPepProMat)){ 
    for(j in 1:(i-1)){
      list1<-corPeplist[[row.names(corPepProMat)[i]]]
      list2<-corPeplist[[colnames(corPepProMat)[j]]]
      interList<-intersect(names(list1),names(list2))
      
      list1I<-list1[interList]
      # list1I[list1I=="#N/A"]=NA
      # list1I <- as.numeric(list1I)
      list2I<-list2[interList]
      # list2I[list2I=="#N/A"]=NA
      # list2I <- as.numeric(list2I)
      corPepProMat[i,j]<-cor(list1I,list2I,use="na.or.complete")
    }
  }
  
  #pro cormat
  for(i in 1:(nrow(corPepProMat)-1)){
    for(j in (i+1):ncol(corPepProMat)){
      list1<-corProlist[[row.names(corPepProMat)[i]]]
      list2<-corProlist[[colnames(corPepProMat)[j]]]
      
      interList<-intersect(names(list1),names(list2))
      
      list1I<-list1[interList]
      list2I<-list2[interList]
      
      corPepProMat[i,j]<-cor(list1I,list2I,use="na.or.complete")
    }
  }
  
  corPepProMats[[which(condName==condNames)]]<-corPepProMat
  print(condName)
}

pdf("Figure3CorrelationPlot_FZ20210729.pdf",width=10,height=20)
par(mfrow=c(3,2),mar=c(2,2,6,2))
for(m in 1:length(corPepProMats)){
  corPepProMat<-corPepProMats[[m]]
  corPepProMat[is.na(corPepProMat)]<-0
  library(corrplot)
  #peptide
  corrplot(corPepProMat, order = "original", type = "lower", method = "ellipse", diag = FALSE,
           tl.srt = 60,tl.col = "black",main=paste(condNames[m],"peptide")
           , cl.pos = "r", col.lim = c(-1,1), cl.length = 21
           , bg = "grey", addCoef.col = "white", col=brewer.pal(n=10, name="RdYlBu"),na.label = "*",number.cex=2,tl.cex	=2)
  
  #protein
  corrplot(corPepProMat, order = "original", type = "upper", method = "ellipse", diag = FALSE,
           tl.srt = 60,tl.col = "black",main=paste(condNames[m],"protein")
           , cl.pos = "r", col.lim = c(-1,1), cl.length = 21
           , bg = "grey", addCoef.col = "white", col=rev(brewer.pal(n=10, name="RdYlBu")),na.label = "*",number.cex=2,tl.cex	=2)
}
dev.off()

