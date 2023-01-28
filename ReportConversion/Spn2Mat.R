#Spectronaut####
#functions
Spn2PreInt<-function(MatIn){
  #remove decoy proteins
  MatOut<-MatIn[!grepl("Decoy",MatIn$PG.ProteinGroups),]
  #remove proteins with [Acetyl (Protein N-term)]
  MatOut2<-MatOut[!grepl("\\[Acetyl \\(Protein N-term\\)\\]",MatOut$EG.PrecursorId),]
  
  MatPre<-as.data.frame(MatOut2[,grepl("EG.TotalQuantity",colnames(MatOut2))],stringsAsFactors=F)
  colnames(MatPre)<-colnames(MatOut2)[grepl("EG.TotalQuantity",colnames(MatOut2))]
  PreNames<-MatOut2$EG.PrecursorId
  
  #change precursor naming format
  PreNames<-gsub("_","",PreNames)
  PreNames<-gsub("\\.","",PreNames)
  
  #change mod names
  PreNames<-gsub("\\[Carbamidomethyl \\(C\\)\\]","(UniMod:4)",PreNames)
  PreNames<-gsub("\\[Oxidation \\(M\\)\\]","(UniMod:35)",PreNames)
  
  row.names(MatPre)<-PreNames
  
  #change column names removed suffix
  
  MatPre[MatPre=="NaN"]<-NA
  MatPre[MatPre=="Filtered"]<-NA
  
  MatPre
}

Spn2PepInt<-function(MatIn){
  #remove decoy proteins
  MatOut<-MatIn[!grepl("Decoy",MatIn$PG.ProteinGroups),]
  #remove proteins with [Acetyl (Protein N-term)]
  MatOut2<-MatOut[!grepl("\\[Acetyl \\(Protein N-term\\)\\]",MatOut$EG.PrecursorId),]
  
  PepNames<-MatOut2$EG.PrecursorId
  PepNames<-gsub("_","",PepNames)
  PepNames<-as.vector(sapply(sapply(PepNames,strsplit,"\\."),"[[",1))
  PepNames<-as.vector(sapply(PepNames,StripePep))
  #change mod names
  PepNames<-gsub("\\[Carbamidomethyl \\(C\\)\\]","(UniMod:4)",PepNames)
  PepNames<-gsub("\\[Oxidation \\(M\\)\\]","(UniMod:35)",PepNames)
  PepNamesNoDup<-PepNames[!duplicated(PepNames)]
  
  MatPep<-as.data.frame(MatOut2[which(!duplicated(PepNames)),grepl("PEP",colnames(MatOut2))]
                        ,stringsAsFactors=F)
  colnames(MatPep)<-colnames(MatOut2)[grepl("PEP",colnames(MatOut2))]
  MatPep[MatPep=="NaN"]<-NA
  MatPep2<-as.data.frame(MatPep[which(apply(is.na(MatPep),1,sum)!=ncol(MatPep)),] ,stringsAsFactors=F)
  PepNamesNoDup2<-PepNamesNoDup[which(apply(is.na(MatPep),1,sum)!=ncol(MatPep))]
  
  #change column names removed suffix
  # MSnames_sp<-sapply(colnames(MatPep),strsplit,"\\.")
  MSnames<-gsub(".PEP.Quantity","",colnames(MatPep))
  
  row.names(MatPep2)<-PepNamesNoDup2
  colnames(MatPep2)<-MSnames
  
  MatPep2
}

Spn2ProInt<-function(MatIn){
  #remove decoy proteins
  MatOut<-MatIn[!grepl("Decoy",MatIn$PG.ProteinGroups),]
  #remove proteins with [Acetyl (Protein N-term)]
  MatOut2<-MatOut[!grepl("\\[Acetyl \\(Protein N-term\\)\\]",MatOut$EG.PrecursorId),]
  
  ProNames<-MatOut2$PG.ProteinGroups
  MatPro<-as.data.frame(MatOut2[which(!duplicated(ProNames)),grepl("PG.Quantity",colnames(MatOut2))],stringsAsFactors=F)
  colnames(MatPro)<-colnames(MatOut2)[grepl("PG.Quantity",colnames(MatOut2))]
  
  ProNames<-ProNames[!duplicated(ProNames)]
  
  #change column names removed suffix
  # MSnames_sp<-sapply(colnames(MatPro),strsplit,"\\.")
  MSnames<-gsub(".PG.Quantity","",colnames(MatPro))
  MatPro[MatPro=="NaN"]<-NA
  ProNames<-ProNames[which(apply(is.na(MatPro),1,sum)!=ncol(MatPro))]
  MatPro<-as.data.frame(MatPro[which(apply(is.na(MatPro),1,sum)!=ncol(MatPro)),],stringsAsFactors=F)
  
  row.names(MatPro)<-ProNames 
  colnames(MatPro)<-MSnames

  MatPro
}

Spn2RT<-function(MatIn){
  #remove decoy proteins
  MatOut<-MatIn[!grepl("Decoy",MatIn$PG.ProteinGroups),]
  #remove proteins with [Acetyl (Protein N-term)]
  MatOut2<-MatOut[!grepl("\\[Acetyl \\(Protein N-term\\)\\]",MatOut$EG.PrecursorId),]
  
  MatRT<-as.data.frame(MatOut2[,grepl("EG.ApexRT",colnames(MatOut2))],stringsAsFactors=F)
  colnames(MatRT)<-colnames(MatOut2)[grepl("EG.ApexRT",colnames(MatOut2))]
  
  PreNames<-MatOut2$EG.PrecursorId
  
  #change precursor naming format
  PreNames<-gsub("_","",PreNames)
  PreNames<-gsub("\\.","",PreNames)
  
  #change mod names
  # PreNames[grepl("\\(",PreNames)]
  PreNames<-gsub("\\[Carbamidomethyl \\(C\\)\\]","(UniMod:4)",PreNames)
  PreNames<-gsub("\\[Oxidation \\(M\\)\\]","(UniMod:35)",PreNames)
  
  row.names(MatRT)<-PreNames
  
  #change column names removed suffix
  MSnames_sp<-sapply(colnames(MatRT),strsplit,"\\.")
  MSnames<-as.vector(sapply(MSnames_sp, "[[", 4))
  
  colnames(MatRT)<-MSnames
  
  MatRT[MatRT=="NaN"]<-NA
  MatRT[MatRT=="Filtered"]<-NA
  MatRT
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

# running for samples
dirAll<-paste0("Result/",LETTERS[6],"/5spn/Export")
xlsAll<-list.files(dirAll,pattern="xls", recursive = T)
xlsFullAll<-list.files(dirAll,pattern="xls",full.names = T, recursive = T)
nameOutAll<-paste(gsub(".xls","",xlsAll),"5spn",sep="_")

for(m in 3:length(xlsFullAll)){
  MatIn<-read.delim(xlsFullAll[m],stringsAsFactors = F,header = T)

  FnamePre<-paste0("Matrix20210724/",nameOutAll[m],"_PreInt.tsv")
  FnamePep<-paste0("Matrix20210724/",nameOutAll[m],"_PepInt.tsv")
  FnamePro<-paste0("Matrix20210724/",nameOutAll[m],"_ProInt.tsv")
  FnameRT<-paste0("Matrix20210724/",nameOutAll[m],"_PreRT.tsv")

  if(!file.exists(FnamePre)){
    MatPre<-Spn2PreInt(MatIn)
    write.table(MatPre,FnamePre,sep="\t")
    print(FnamePre)
  }

  if(!file.exists(FnamePep)){
    MatPep<-Spn2PepInt(MatIn)
    write.table(MatPep,FnamePep,sep="\t")
    print(FnamePep)
  }
  
  if(!file.exists(FnamePro)){
    MatPro<-Spn2ProInt(MatIn)
    write.table(MatPro,FnamePro,sep="\t")
    print(FnamePro)
  }
  # 
  if(!file.exists(FnameRT)){
    MatRT<-Spn2RT(MatIn)
    write.table(MatRT,FnameRT,sep="\t")
    print(FnameRT)
  }
}
