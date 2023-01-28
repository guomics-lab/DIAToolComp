prot_files<-list.files(pattern = "_ProInt",path="./Matrix")

sampleAList<-list(c(1,3,5),c(1,2,3),c(1,2,3))
sampleBList<-list(c(2,4,6),c(4,5,6),c(4,5,6))
names(sampleAList)<-names(sampleBList)<-LETTERS[1:3]
fastaAll<-read.csv("Library/FastaHYEEM.csv",stringsAsFactors = F,header = F,row.names = 1)
for(condName in LETTERS[1:3]){
  fileAll<-prot_files[grepl(condName,prot_files)]
  
  library(seqinr)
  fastaLib<-list.files(pattern = ".fasta$",paste0("Library/lib",condName)
                       ,full.names = T)  
  fastaIn<-read.fasta(file = fastaLib,seqtype = "AA",as.string = T)
  fasta_seqAll<-sapply(fastaIn,"[[",1)
  fastapro<-names(fasta_seqAll)
  fastaspec<-sapply(strsplit(fastapro, "\\_"), function(x) x[2], simplify=T)
  fastaID<-sapply(strsplit(fastapro, "\\|"), function(x) x[2], simplify=T)
  names(fastaspec)<-fastaID      
  
  for(m in 1:length(fileAll)){
    fileout<-gsub("ProInt","ProRatio",fileAll[m])
    
    if(!file.exists(paste0('Matrix/',fileout))){
      readIn<-read.table(paste0('Matrix/',fileAll[m]),stringsAsFactors = F,sep='\t',header = T,row.names=1)
      #remove protein groups
      readIn2<-readIn[!grepl(";",row.names(readIn)),]
      
      specAll<-fastaspec[row.names(readIn2)]
      specAll<-gsub("YEAS8","YEAST",specAll)
      
      Amean<-apply(readIn2[,sampleAList[[condName]]],1,mean,na.rm=T)
      Bmean<-apply(readIn2[,sampleBList[[condName]]],1,mean,na.rm=T)
      
      ABratio<-Amean/Bmean
      RatMat<-cbind(Amean,Bmean,ABratio,specAll)
      RatMat[RatMat=="NaN"]<-NA
      row.names(RatMat)<-row.names(readIn2)
      
      RatMat<-RatMat[!is.na(RatMat[,3]),]
      RatMat<-RatMat[!is.na(row.names(RatMat)),]
      write.table(RatMat,paste0('Matrix/',fileout),sep="\t")
      print(fileout)
    }
  }
}

prot_files<-list.files(pattern = "_ProInt",path="./Matrix")
for(condName in LETTERS[4:6]){
  fileAll<-prot_files[grepl(condName,prot_files)]
  fileAllS<-fileAll[grepl("S",fileAll)]
  fileAllL<-fileAll[grepl("L",fileAll)]
  
  for(m in 1:length(fileAllS)){
    fileout<-gsub("ProInt","ProRatio",fileAllS[m])
    # if(!file.exists(paste0('../../Matrix/',fileout))){
    readInS0<-read.table(paste0('./Matrix/',fileAllS[m]),stringsAsFactors = F,sep='\t'
                         ,header = T,row.names = 1)
    readInL0<-read.table(paste0('./Matrix/',fileAllL[m]),stringsAsFactors = F,sep='\t'
                         ,header = T,row.names = 1)
    #remove protein groups
    readInS<-as.data.frame(readInS0[which(apply(is.na(readInS0),1,sum)!=ncol(readInS0)&!grepl(";",row.names(readInS0))),]
                           ,StringAsFactors=F)
    readInL<-as.data.frame(readInL0[which(apply(is.na(readInL0),1,sum)!=ncol(readInL0)&!grepl(";",row.names(readInL0))),]
                           ,StringAsFactors=F)
    row.names(readInS)<-row.names(readInS0)[which(apply(is.na(readInS0),1,sum)!=ncol(readInS0)&!grepl(";",row.names(readInS0)))]
    row.names(readInL)<-row.names(readInL0)[which(apply(is.na(readInL0),1,sum)!=ncol(readInL0)&!grepl(";",row.names(readInL0)))]
    colnames(readInS)<-colnames(readInS0)
    colnames(readInL)<-colnames(readInL0)
    
    proAll<-unique(c(row.names(readInS),row.names(readInL)))
    
    RatMat<-data.frame(matrix(NA,nrow=length(proAll),ncol = 3),stringsAsFactors = F)
    colnames(RatMat)<-c("Amean","Bmean","ABratio")
    row.names(RatMat)<-proAll
    
    Amean<-apply(readInS,1,mean,na.rm=T)
    Bmean<-apply(readInL,1,mean,na.rm=T)
    RatMat[names(Amean),"Amean"]<-Amean
    RatMat[names(Bmean),"Bmean"]<-Bmean
    RatMat$ABratio<-(RatMat$Amean)/(RatMat$Bmean)
    
    RatMat[RatMat=="NaN"]<-NA
    RatMat<-RatMat[!is.na(RatMat$ABratio),]
    RatMat<-RatMat[!is.na(row.names(RatMat)),]
    
    print(fileout)
    write.table(RatMat,paste0('./Matrix/',fileout),sep="\t")
    # }
  }
}
