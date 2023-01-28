#DIANN####
#functions
Sky2PreInt<-function(MatIn){
  MatOut<-MatIn[MatIn$Protein!="Decoys",]
  MatQvalue<-MatOut[,grepl("Q.Value",colnames(MatOut))]
  MatArea<-MatOut[,grepl("Total.Area",colnames(MatOut))]
  MatArea[MatQvalue>0.01]<-NA
  MatArea[MatArea=="#N/A"]<-NA
  MatArea[MatArea==0]<-NA
  
  preIDAll<-paste0(MatOut$Modified.Sequence,MatOut$Precursor.Charge)
  preIDAllUniq<-preIDAll[!duplicated(preIDAll)]
  MatAreaOut<-MatArea[!duplicated(preIDAll),]
  
  preIDAllUniq<-gsub("\\[\\+57\\]","(UniMod:4)",preIDAllUniq)
  preIDAllUniq<-gsub("\\[\\+16\\]","(UniMod:35)",preIDAllUniq)
  
  row.names(MatAreaOut)<-preIDAllUniq
  MatAreaOut<-MatAreaOut[apply(is.na(MatAreaOut),1,sum)!=ncol(MatAreaOut),]
  
  MatAreaOut
}

Pre2Pep<-function(PreIn){
  pepPre<-as.vector(sapply(row.names(PreIn),StripeCharge))
  pepAll<-unique(pepPre)
  
  pepOut<-data.frame(matrix(NA,nrow=length(pepAll),ncol=ncol(PreIn)))
  colnames(pepOut)<-colnames(PreIn)
  row.names(pepOut)<-pepAll
  
  for(i in 1:length(pepAll)){#
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
    if(i%%1000==0) print(i/length(pepAll)*100)
  }
  pepOut[pepOut=="NaN"]<-NA
  pepOut<-pepOut[apply(is.na(pepOut),1,sum)!=ncol(pepOut),]
  
  pepOut
}

Sky2PE<-function(MatIn){
  #remove decoys
  MatOut<-MatIn[MatIn$Protein!="Decoys",]
  MatOut$Protein<-gsub("sp\\|","",MatOut$Protein)
  MatOut$Protein<-gsub("tr\\|","",MatOut$Protein)
  
  MatOut<-MatOut[grepl("\\|",MatOut$Protein),]
  
  MatOut$Protein<- as.vector(sapply(strsplit(MatOut$Protein,"\\|"),"[[",1))
  
  MatQvalue<-MatOut[,grepl("Q.Value",colnames(MatOut))]
  MatArea<-MatOut[,grepl("Total.Area",colnames(MatOut))]
  MatArea[MatQvalue>0.01]<-NA
  MatArea[MatArea=="#N/A"]<-NA
  MatArea[MatArea==0]<-NA
  
  preIDAll<-paste(MatOut$Modified.Sequence,MatOut$Precursor.Charge,sep="_")
  preIDAll<-gsub("\\[\\+57\\]","(UniMod:4)",preIDAll)
  preIDAll<-gsub("\\[\\+16\\]","(UniMod:35)",preIDAll)
  
  preIDAllUniq<-preIDAll[!duplicated(preIDAll)]
  MatAreaOut<-MatArea[!duplicated(preIDAll),]
  
  proUniq<-c()
  for(i in 1:length(preIDAllUniq)){#
    indSel<-which(preIDAll==preIDAllUniq[i])
    proUniq[i]<-paste(MatOut$Protein[indSel],collapse = ";")
    if(i%%1000==0) print(paste(round(i/length(preIDAllUniq)*100,2),"%"))
  }
  prePE<-cbind(proUniq,MatAreaOut)
  row.names(prePE)<-preIDAllUniq
  prePEOut<-prePE[apply(is.na(prePE),1,sum)!=(ncol(prePE)-1),]
  
  prePEOut
}

Sky2RT<-function(MatIn){
  MatOut<-MatIn[MatIn$Protein!="Decoys",]
  MatQvalue<-MatOut[,grepl("Q.Value",colnames(MatOut))]
  MatRT<-MatOut[,grepl("Retention.Time",colnames(MatOut))]
  MatRT[MatQvalue>0.01]<-NA
  MatRT[MatRT=="#N/A"]<-NA
  
  preIDAll<-paste0(MatOut$Modified.Sequence,MatOut$Precursor.Charge)
  preIDAll_sp<-as.vector(unlist(sapply(preIDAll,StripePep)))
  preIDAllUniq<-preIDAll[!duplicated(preIDAll)]
  
  preIDAllUniq<-gsub("\\[\\+57\\]","(UniMod:4)",preIDAllUniq)
  preIDAllUniq<-gsub("\\[\\+16\\]","(UniMod:35)",preIDAllUniq)
  
  MatRTOut<-MatRT[!duplicated(preIDAll),]
  
  row.names(MatRTOut)<-preIDAllUniq
  MatRTOut<-MatRTOut[apply(is.na(MatRTOut),1,sum)!=ncol(MatRTOut),]
  
  MatRTOut
}

StripePep<-function(pepIn){
  pepIn_sp<-unlist(strsplit(pepIn,""))
  leftBrk<-which(grepl("\\[",pepIn_sp))
  rightBrk<-which(grepl("\\]",pepIn_sp))
  
  if(length(leftBrk)>0){
    indDel<-c()
    for(j in 1:length(leftBrk)){
      indDel<-c(indDel,leftBrk[j]:rightBrk[j])
    }
    pepOut<- paste0(pepIn_sp[-indDel],collapse = "") 
  }else{
    pepOut<-pepIn
  }
  pepOut
}
StripeCharge<-function(preIn){
  preIn_sp<-unlist(strsplit(preIn,""))
  paste0(preIn_sp[1:(length(preIn_sp)-1)],collapse = "")
}

# running for samples
dirAll<-paste0("Result/",LETTERS[c(1)],"/4sky/")
prcsvAll<- list.files(dirAll,pattern="_PreAreaRTQvalue20210512.csv", recursive = T,full.names = T)
prcsv<-list.files(dirAll,pattern="_PreAreaRTQvalue20210512.csv", recursive = T)
nameOutAll<-paste0(gsub("_PreAreaRTQvalue20210512.csv","",prcsv),"_4sky")

for(m in 1:1){
  MatIn<-read.delim(prcsvAll[m],header = T
                    ,stringsAsFactors = F,sep = ",")
  
  FnamePre<-paste0("Matrix20210724/",nameOutAll[m],"_PreInt.tsv")
  FnamePep<-paste0("Matrix20210724/",nameOutAll[m],"_PepInt.tsv")
  FnamePE<-paste0("Matrix20210724/",nameOutAll[m],"_ToPE.csv")
  FnameRT<-paste0("Matrix20210724/",nameOutAll[m],"_PreRT.tsv")
  
  if(!file.exists(FnamePre)){
    PreIn<-Sky2PreInt(MatIn)
    write.table(PreIn,FnamePre,sep="\t")
    print(FnamePre)
  }

  if(!file.exists(FnamePep)){
    PepOut<-Pre2Pep(Sky2PreInt(MatIn))
    write.table(PepOut,FnamePep,sep="\t")
    print(FnamePep)
  }
  
  # if(!file.exists(FnamePE)){
    proOut<-Sky2PE(MatIn)
    proOut[proOut==0]<-NA
    
    write.csv(proOut,FnamePE)
    print(FnamePE)
  # }
  
  # if(!file.exists(FnameRT)){
  #   PreRT<-Sky2RT(MatIn)
  #   write.table(PreRT,FnameRT,sep="\t")
  #   print(FnameRT)
  # }
}
