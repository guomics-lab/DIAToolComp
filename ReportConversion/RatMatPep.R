StripePep2<-function(pepIn){
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
StripePep<-function(pepIn){
  pepIn_sp<-unlist(strsplit(pepIn,""))
  leftBrk<-which(grepl("\\(",pepIn_sp))
  rightBrk<-which(grepl("\\)",pepIn_sp))
  
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

pep_files = list.files(pattern = "PepInt",path = "./Matrix/")
sampleAList<-list(c(1,3,5),c(1,2,3),c(1,2,3))
sampleBList<-list(c(2,4,6),c(4,5,6),c(4,5,6))
names(sampleAList)<-names(sampleBList)<-LETTERS[1:3]

#free ABC_2dnn####
for(condName in LETTERS[1:3]){
  #peptide
  fileAll<-pep_files[grepl(condName,pep_files)]
  fileFree<-fileAll[grepl("free",fileAll,fixed = T)]
  
  fastaLib<-list.files(pattern = ".fasta$",paste0("./Library/lib",condName)
                       ,full.names = T)
  fastaIn<-read.fasta(file = fastaLib,seqtype = "AA",as.string = T)
  fasta_seqAll<-sapply(fastaIn,"[[",1)
  fastapro<-names(fasta_seqAll)
  fastaspec<-sapply(strsplit(fastapro, "\\_"), function(x) x[2], simplify=T)
  fastaID<-sapply(strsplit(fastapro, "\\|"), function(x) x[2], simplify=T)
  names(fastaspec)<-fastaID
  
  for(m in which(grepl("2dnn",fileFree))){
    fileout<-gsub("PepInt","PepRatio",fileFree[m])
    if(!file.exists(paste0('./Matrix/',fileout))){
      readIn<-read.table(paste0('./Matrix/',fileFree[m]),stringsAsFactors = F,sep='\t')
      pepAll<-sapply(row.names(readIn),StripePep)#without modifications
      #for library free
      #read Dnn's result
      dnnExport<-paste0("./Result/",condName,"/2dnn/ForUse/",condName,"_free_S.pr_matrix.tsv")
      readDnn<-read.delim(dnnExport,stringsAsFactors = F
                          ,header = T)    
      readDnn2<-unique(cbind(readDnn$Protein.Ids,readDnn$Stripped.Sequence))
      row.names(readDnn2)<-readDnn2[,2]
      pgAll<-readDnn2[pepAll,1] 
      
      pgAll1<-sapply(strsplit(pgAll, ";"), function(x) x[1], simplify=T) # first protein in protein group
      specAll<-fastaspec[pgAll1]
      specAll<-gsub("YEAS8","YEAST",specAll)
      
      Amean<-apply(readIn[,sampleAList[[condName]]],1,mean,na.rm=T)
      Bmean<-apply(readIn[,sampleBList[[condName]]],1,mean,na.rm=T)
      ABratio<-Amean/Bmean
      RatMat<-cbind(Amean,Bmean,ABratio,specAll)
      RatMat[RatMat=="NaN"]<-NA
      row.names(RatMat)<-pepAll
     
      write.table(RatMat,paste0('./Matrix/',fileout),sep="\t")
      print(fileout)
    }
  }
}

# free_ABC_5spn#### 
for(condName in LETTERS[1:3]){
  #peptide
  fileAll<-pep_files[grepl(condName,pep_files)]
  fileFree<-fileAll[grepl("free",fileAll,fixed = T)]
  
  fastaLib<-list.files(pattern = ".fasta$",paste0("./Library/lib",condName)
                       ,full.names = T)
  fastaIn<-read.fasta(file = fastaLib,seqtype = "AA",as.string = T)
  fasta_seqAll<-sapply(fastaIn,"[[",1)
  fastapro<-names(fasta_seqAll)
  fastaspec<-sapply(strsplit(fastapro, "\\_"), function(x) x[2], simplify=T)
  fastaID<-sapply(strsplit(fastapro, "\\|"), function(x) x[2], simplify=T)
  names(fastaspec)<-fastaID
  
  for(m in which(grepl("5spn",fileFree))){
    fileout<-gsub("PepInt","PepRatio",fileFree[m])
    if(!file.exists(paste0('./Matrix/',fileout))){
      readIn<-read.table(paste0('./Matrix/',fileFree[m]),stringsAsFactors = F,sep='\t')
      pepAll<-sapply(row.names(readIn),StripePep)#without modifications
      #for library free
      #read Spn's result
      spnExport<-paste0("./Result/",condName,"/5spn/Export/",condName,"_free_S.xls")
      readSpn<-read.delim(spnExport,stringsAsFactors = F
                          ,header = T)    
      readSpn2<-unique(cbind(readSpn$PG.ProteinGroups,readSpn$EG.PrecursorId))
      readSpn2[,2]<-gsub("_","",readSpn2[,2])
      readSpn2[,2]<-sapply(strsplit(readSpn2[,2], "\\."), function(x) x[1], simplify=T)
      readSpn2[,2]<-sapply(readSpn2[,2], StripePep2)
      readSpn3<-unique(readSpn2)
      row.names(readSpn3)<-readSpn3[,2]
      pgAll<-readSpn3[pepAll,1]
      
      pgAll1<-sapply(strsplit(pgAll, ";"), function(x) x[1], simplify=T) # first protein in protein group
      specAll<-fastaspec[pgAll1]
      specAll<-gsub("YEAS8","YEAST",specAll)
      
      Amean<-apply(readIn[,sampleAList[[condName]]],1,mean,na.rm=T)
      Bmean<-apply(readIn[,sampleBList[[condName]]],1,mean,na.rm=T)
      ABratio<-Amean/Bmean
      RatMat<-cbind(Amean,Bmean,ABratio,specAll)
      RatMat[RatMat=="NaN"]<-NA
      row.names(RatMat)<-pepAll
      
      print(fileout)
      write.table(RatMat,paste0('./Matrix/',fileout),sep="\t")
    }
  }
}

# lib_ABC ####
# Ratio Matrix
for(condName in LETTERS[2]){
  #for library based
  fileAll<-pep_files[grepl(condName,pep_files)]
  fileAll<-fileAll[grepl("lib",fileAll)]
  
  libSpn<-list.files(pattern = ".xls$",paste0("./Library/lib",condName,"/5spn/")
                     ,full.names = T)
  libIn<-read.delim(libSpn,stringsAsFactors = F)
  firstName<-sapply(strsplit(libIn$Protein.Name, ","), function(x) x[1], simplify=T)
  firstName<-sapply(strsplit(firstName, ";"), function(x) x[1], simplify=T)
  libIn$Species<-sapply(strsplit(firstName, "_"), function(x) x[2], simplify=T)
  libUse<-unique(libIn[,c("StrippedPeptide","Species")])
  libSeqAll<-libUse[,2]
  names(libSeqAll)<-libUse[,1]
  
  for(m in 1:length(fileAll)){
    fileout<-gsub("PepInt","PepRatio",fileAll[m])
    # if(!file.exists(paste0('./Matrix/',fileout))){
      readIn<-read.table(paste0('./Matrix/',fileAll[m]),stringsAsFactors = F,sep='\t',header = T,row.names = 1)
      pepAll<-sapply(row.names(readIn),StripePep)
      specAll<-libSeqAll[pepAll]
      specAll<-gsub("YEAS8","YEAST",specAll)
      
      readA<-readIn[,sampleAList[[condName]]]
      readB<-readIn[,sampleBList[[condName]]]
      
      Amean<-apply(readA,1,mean,na.rm=T)
      Bmean<-apply(readB,1,mean,na.rm=T)
      ABratio<-Amean/Bmean
      RatMat<-cbind(Amean,Bmean,ABratio,specAll)
      RatMat[RatMat=="NaN"]<-NA
      row.names(RatMat)<-pepAll
      
      print(fileout)
      write.table(RatMat,paste0('./Matrix/',fileout),sep="\t")
    }
  # }
}

# freeOrLib_DEF#### 
for(condName in LETTERS[4:6]){
  fileAll<-pep_files[grepl(condName,pep_files)]
  fileAllS<-fileAll[grepl("S",fileAll)]
  fileAllL<-fileAll[grepl("L",fileAll)]
  
  for(m in 1:length(fileAllS)){
    fileout<-gsub("PepInt","PepRatio",fileAllS[m])
    if(!file.exists(paste0('./Matrix/',fileout))){
      readInS<-read.table(paste0('./Matrix/',fileAllS[m]),stringsAsFactors = F,sep='\t')
      readInL<-read.table(paste0('./Matrix/',fileAllL[m]),stringsAsFactors = F,sep='\t')
      
      pepAllS<-sapply(row.names(readInS),StripePep)
      pepAllL<-sapply(row.names(readInL),StripePep)
      pepAll<-unique(c(pepAllS,pepAllL))
       
      RatMat<-data.frame(matrix(NA,nrow=length(pepAll),ncol = 3),stringsAsFactors = F)
      colnames(RatMat)<-c("Amean","Bmean","ABratio")
      row.names(RatMat)<-pepAll
      
      Amean<-apply(readInS,1,mean,na.rm=T)
      Bmean<-apply(readInL,1,mean,na.rm=T)
      RatMat[names(Amean),"Amean"]<-Amean
      RatMat[names(Bmean),"Bmean"]<-Bmean
      RatMat$ABratio<-RatMat$Amean/RatMat$Bmean
      RatMat[RatMat=="NaN"]<-NA
      
      write.table(RatMat,paste0('./Matrix/',fileout),sep="\t")
      print(fileout)
    }
  }
}
