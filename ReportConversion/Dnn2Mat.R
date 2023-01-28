#DIANN####
#functions
Dnn2PreInt<-function(prIn){
  indPre<-which(colnames(prIn)=="Precursor.Id")
  preOut<-as.data.frame(prIn[,(1+indPre):ncol(prIn)],StringAsFactors=F)
  row.names(preOut)<-prIn$Precursor.Id
  
  MSnames_sp<-sapply(colnames(prIn)[(1+indPre):ncol(prIn)],strsplit,"\\.")
  MSnames<-as.vector(sapply(MSnames_sp, "[[", 7))
  colnames(preOut)<-MSnames
  preOut
}

Dnn2PepInt<-function(prIn){
  indPre<-which(colnames(prIn)=="Precursor.Id")
  preOut<-as.data.frame(prIn[,(1+indPre):ncol(prIn)],StringAsFactors=F)
  
  pepAll<-unique(prIn$Stripped.Sequence)
  pepOut<-data.frame(matrix(NA,nrow=length(pepAll),ncol=ncol(preOut)))
  row.names(pepOut)<-pepAll
  
  for(i in 1:length(pepAll)){
    indPep<-which(prIn$Stripped.Sequence==pepAll[i])
    if(length(indPep)>1){
      df_preInSel<-as.data.frame(preOut[indPep,],StringAsFactors=F)
      meanAll<-c()
      for(j in 1:ncol(df_preInSel)){
        meanAll[j]<-mean(as.numeric(df_preInSel[,j]),na.rm=T)
      }
      pepOut[i,]<-meanAll
    }else{
      pepOut[i,]<-preOut[indPep,]
    }
  }
  pepOut
}

Dnn2ProInt<-function(pgIn){
  proOut<-as.data.frame(pgIn[,6:ncol(pgIn)],StringAsFactors=F)
  MSnames_sp<-sapply(colnames(pgIn)[6:ncol(pgIn)],strsplit,"\\.")
  MSnames<-as.vector(sapply(MSnames_sp, "[[", 7))
  colnames(proOut)<-MSnames
  row.names(proOut)<-pgIn$Protein.Group
  proOut
}

Dnn2RT<-function(MatIn){
  MatRT<-MatIn[,c("File.Name","Precursor.Id","RT")]
  preAll<-unique(MatIn$Precursor.Id)
  fileAll<-unique(MatIn$File.Name)
  rtOut<-data.frame(matrix(NA,nrow=length(preAll),ncol=length(fileAll)))
  row.names(rtOut)<-preAll
  colnames(rtOut)<-fileAll
  
  for(i in 1:nrow(MatRT)){
    rtOut[MatRT$Precursor.Id[i],MatRT$File.Name[i]]<-MatRT$RT[i]
  }
  
  MSnames_sp<-sapply(colnames(rtOut),strsplit,"\\\\")
  MSnames<-as.vector(sapply(MSnames_sp, "[[", 6))
  MSnames_sp2<-sapply(MSnames,strsplit,"-")
  MSnames2<-as.vector(sapply(MSnames_sp2, "[[", 1))
  
  colnames(rtOut)<-MSnames2
  rtOut
}

# running for samples
dirAll<-paste0("Result/",LETTERS[1:6],"/2dnn/ForUse")
prtsvAll<- list.files(dirAll,pattern="pr_matrix.tsv", recursive = T,full.names = T)
pgtsvAll<-list.files(dirAll,pattern="pg_matrix.tsv", recursive = T,full.names = T)
mattsvAll<-list.files(dirAll,pattern=".tsv", recursive = T,full.names = T)
mattsvAll<-mattsvAll[!grepl("matrix",mattsvAll)]
mattsv<-list.files(dirAll,pattern=".tsv", recursive = T)
mattsv<-mattsv[!grepl("matrix",mattsv)]
nameOutAll<-paste0(gsub(".tsv","",mattsv),"_2dnn")

for(m in 1:length(mattsv)){
  # pgIn<-read.delim(pgtsvAll[m],sep="\t",header = T,stringsAsFactors = F)
  prIn<-read.delim(prtsvAll[m],sep="\t",header = T,stringsAsFactors = F)
  # MatIn<-read.delim(mattsvAll[m],sep="\t",header = T,stringsAsFactors = F)

  FnamePre<-paste0("Matrix20210724/",nameOutAll[m],"_PreInt.tsv")
  FnamePep<-paste0("Matrix20210724/",nameOutAll[m],"_PepInt.tsv")
  FnamePro<-paste0("Matrix20210724/",nameOutAll[m],"_ProInt.tsv")
  FnameRT<-paste0("Matrix20210724/",nameOutAll[m],"_PreRT.tsv")
  
  # if(!file.exists(FnamePre)){
  #   prOut<-Dnn2PreInt(prIn)
  #   write.table(prOut,FnamePre,sep="\t")
  #   print(FnamePre)
  # }

  if(!file.exists(FnamePep)){
    pepOut<-Dnn2PepInt(prIn)
    write.table(pepOut,FnamePep,sep="\t")
    print(FnamePep)
  }
  
  # if(!file.exists(FnamePro)){
  #   proOut<-Dnn2ProInt(pgIn)
  #   write.table(proOut,FnamePro,sep="\t")
  #   print(FnamePro)
  # }
  
  # if(!file.exists(FnameRT)){
  #   rtOut<-Dnn2RT(MatIn)
  #   write.table(rtOut,FnameRT,sep="\t")
  #   print(FnameRT)
  # }
}
