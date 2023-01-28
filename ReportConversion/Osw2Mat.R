#OpenSWATH####
#functions
Osw2PreInt<-function(readPepList){
  PepIntList<-list()
  for(m in 1:length(readPepList)){
    readIn<-readPepList[[m]]
    readIn<-readIn[readIn$decoy==0&readIn$m_score<0.01&readIn$peak_group_rank==1,]
    readIn<-unique(readIn)
    
    readPep<-unique(readIn[,grepl("FullPeptideName|Charge|Intensity" ,colnames(readIn))])
    readPep2<-unique(cbind(paste0(readPep$FullPeptideName,readPep$Charge),readPep$Intensity))
    readPepInt<-as.numeric(readPep2[,2])
    names(readPepInt)<-readPep2[,1]
    PepIntList[[m]]<-readPepInt
  }
  
  preAll<-unique(unlist(sapply(PepIntList,names)))
  dfOut<-data.frame(matrix(NA,nrow=length(preAll),ncol=length(PepIntList)),stringsAsFactors = F)
  row.names(dfOut)<-preAll
  
  for(m in 1:length(PepIntList)){
    dfOut[names(PepIntList[[m]]),m]<-PepIntList[[m]]
  }
  dfOut
}
StripeCharge<-function(preIn){
  preIn_sp<-unlist(strsplit(preIn,""))
  paste0(preIn_sp[1:(length(preIn_sp)-1)],collapse = "")
}

Pre2Pep<-function(PreIn){
  pepPre<-as.vector(sapply(row.names(PreIn),StripeCharge))
  pepAll<-unique(pepPre)
  
  pepOut<-data.frame(matrix(NA,nrow=length(pepAll),ncol=ncol(PreIn)))
  colnames(pepOut)<-colnames(PreIn)
  row.names(pepOut)<-pepAll
  
  for(i in 1:length(pepAll)){
    indPep<-which(pepPre==pepAll[i])
    if(length(indPep)>1){
      df_preInSel<-as.data.frame(PreIn[indPep,],StringAsFactors=F)
      meanAll<-c()
      for(j in 1:ncol(df_preInSel)){
        meanAll[j]<-mean(as.numeric(df_preInSel[,j]),na.rm=T)
      }
      pepOut[i,]<-meanAll
    }else{
      pepOut[i,]<-PreIn[indPep,]
    }
    # print(i)
  }
  pepOut[pepOut=="NaN"]<-NA
  pepOut
}

Osw2PE<-function(readPepList){
  PepIntList<-list()
  for(m in 1:length(readPepList)){
    readIn<-readPepList[[m]]
    #remove decoy m score <0.01 p rank=1
    readIn<-readIn[readIn$decoy==0&readIn$m_score<0.01&readIn$peak_group_rank==1,]
    readIn<-unique(readIn)
    
    readPep<-unique(readIn[,grepl("FullPeptideName|Charge|Intensity|Protein" ,colnames(readIn))])
    readPep2<-unique(cbind(paste(paste(readPep$FullPeptideName,readPep$Charge,sep="_")
                                 ,readPep$ProteinName,sep="-")
                           ,readPep$Intensity))
    readPepInt<-as.numeric(readPep2[,2])
    names(readPepInt)<-readPep2[,1]
    PepIntList[[m]]<-readPepInt
  }
  
  preproAll<-unique(unlist(sapply(PepIntList,names)))
  peOut<-data.frame(matrix(NA,nrow=length(preproAll),ncol=length(PepIntList)),stringsAsFactors = F)
  row.names(peOut)<-preproAll
  
  for(m in 1:length(PepIntList)){
    peOut[names(PepIntList[[m]]),m]<-PepIntList[[m]]
  }
  
  preAll<-sapply(strsplit(preproAll,"-"),"[[",1)
  proAll<-sapply(strsplit(preproAll,"-"),"[[",2)
  
  PEOut<-cbind(proAll,peOut)
  row.names(PEOut)<-preAll
  uniAll<-sapply(sapply(PEOut$proAll,strsplit,"\\|"),"[[",1)
  PEOut$proAll<-uniAll
  PEOut
}

Osw2RT<-function(readPepList){
  PepRTList<-list()
  for(m in 1:length(readPepList)){
    readIn<-readPepList[[m]]
    #remove decoy m score <0.01 p rank=1
    readIn<-readIn[readIn$decoy==0&readIn$m_score<0.01&readIn$peak_group_rank==1,]
    readIn<-unique(readIn)
    
    readPep<-unique(readIn[,grepl("FullPeptideName|Charge|^RT" ,colnames(readIn))])
    readPep2<-unique(cbind(paste0(readPep$FullPeptideName,readPep$Charge),readPep$RT))
    readPepRT<-as.numeric(readPep2[,2])
    names(readPepRT)<-readPep2[,1]
    PepRTList[[m]]<-readPepRT
  }
  
  preAll<-unique(unlist(sapply(PepRTList,names)))
  rtOut<-data.frame(matrix(NA,nrow=length(preAll),ncol=length(PepRTList)),stringsAsFactors = F)
  row.names(rtOut)<-preAll
  MSnames<-sapply(strsplit(names(readPepList),"\\/"),"[[",5)
  MSnames<-sapply(strsplit(MSnames,"\\."),"[[",1)
  colnames(rtOut)<-MSnames
  
  for(m in 1:length(PepRTList)){
    rtOut[names(PepRTList[[m]]),m]<-PepRTList[[m]]/60
  }
  rtOut
}

# running for samples
DataIDs<-c("A_S","B_S","D_S","D_L","E_S","E_L")#""D_L","E_S","E_L",
for(DataID in DataIDs){
  DataID_sp<-unlist(strsplit(DataID,"_"))
  dirAll<-paste0("Result/",DataID_sp[1],"/3osw/",DataID_sp[2])
  resultFiles<-list.files(dirAll,pattern = ".tsv",recursive = F,full.names = T)
  readPepList<-list()
  for(resultFile in resultFiles){
    readIn<-read.delim(resultFile,stringsAsFactors = F,sep='\t')
    readPepList[[which(resultFiles==resultFile)]]<-readIn
  }
  names(readPepList)<-resultFiles
  
  nameOut<-paste0(DataID_sp[1],"_lib_",DataID_sp[2],"_3osw")
  FnamePre<-paste0("Matrix20210724/",nameOut,"_PreInt.tsv")
  FnamePep<-paste0("Matrix20210724/",nameOut,"_PepInt.tsv")
  FnamePE<-paste0("Matrix20210724/",nameOut,"_ToPE.csv")
  FnameRT<-paste0("Matrix20210724/",nameOut,"_PreRT.tsv")
  
  # if(!file.exists(FnamePre)){
  #   PreIn<-Osw2PreInt(readPepList)
  #   write.table(PreIn,FnamePre,sep="\t")
  #   print(FnamePre)
  # }
  
  if(!file.exists(FnamePep)){
    PepOut<-Pre2Pep(Osw2PreInt(readPepList))
    write.table(PepOut,FnamePep,sep="\t")
    print(FnamePep)
  }
  
  # if(!file.exists(FnamePE)){
  #   proOut<-Osw2PE(readPepList)
  #   write.csv(proOut,FnamePE)
  #   print(FnamePE)
  # }
  
  # if(!file.exists(FnameRT)){
  #   PreRT<-Osw2RT(readPepList)
  #   write.table(PreRT,FnameRT,sep="\t")
  #   print(FnameRT)
  # }
}

