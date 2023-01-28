#DIANN####
#functions
pepIn<-read.delim("Result/A/1eny/Export/A_lib_S.elib.peptides.txt",stringsAsFactors = F)
proIn<-read.delim("Result/A/1eny/Export/A_lib_S.elib.proteins.txt",stringsAsFactors = F)

rtAll<-list.files(path = "Result/A/1eny/Export/",pattern = "rt_fit", recursive = T,full.names = T)

Eny2PepInt<-function(pepIn){
  pepOut<-as.data.frame(pepIn[,4:ncol(pepIn)],StringAsFactors=F)
  
  pepAll<-pepIn$Peptide
  pepAll<-gsub("C\\[\\+57\\.021464\\]","C(UniMod:4)",pepAll)
  pepAll<-gsub("M\\[\\+15\\.994915\\]","M(UniMod:35)",pepAll)
  
  MSnames_sp<-sapply(colnames(pepOut),strsplit,"\\.")
  MSnames<-as.vector(sapply(MSnames_sp, "[[", 1))
  colnames(pepOut)<-MSnames
  row.names(pepOut)<-pepAll
  
  pepOut
}

Eny2ProInt<-function(proIn){
  proOut<-as.data.frame(proIn[,4:ncol(proIn)],StringAsFactors=F)
  MSnames_sp<-sapply(colnames(proIn)[4:ncol(proIn)],strsplit,"\\.")
  MSnames<-as.vector(sapply(MSnames_sp, "[[", 1))
  colnames(proOut)<-MSnames
  uniprotAll<-sapply(proIn$Protein,strsplit,"\\|")
  uniprotAll<-as.vector(sapply(uniprotAll, "[[", 2))
  row.names(proOut)<-uniprotAll
  
  proOut
}

Eny2RT<-function(rtAll){
  pepSeq<-c()
  rtIn<-list()
  for(i in 1:length(rtAll)){
    MatRT<-read.delim(rtAll[i],stringsAsFactors = F)
    rtIn[[i]]<-MatRT
    pepSeq<-c(pepSeq, rtIn[[i]]$sequence)
  }
  pepSeqAll<-unique(pepSeq)
  
  rtOut<-data.frame(matrix(NA,nrow=length(pepSeqAll),ncol=length(rtAll)))
  for(m in 1:length(rtIn)){
    MatRT<-rtIn[[m]]
    indSel<-which(pepSeq%in%MatRT$sequence)
    length(indSel)
    length(MatRT$actual)
    nrow(rtOut[indSel,])
    MatRT$actual
  }
 
  MatRT<-MatIn[,c("File.Name","Precursor.Id","RT")]
  preAll<-unique(MatIn$Precursor.Id)
  fileAll<-unique(MatIn$File.Name)
  row.names(rtOut)<-preAll
  colnames(rtOut)<-fileAll
  
  MSnames_sp<-sapply(colnames(rtOut),strsplit,"\\\\")
  MSnames<-as.vector(sapply(MSnames_sp, "[[", 6))
  MSnames_sp2<-sapply(MSnames,strsplit,"-")
  MSnames2<-as.vector(sapply(MSnames_sp2, "[[", 1))
  
  colnames(rtOut)<-MSnames2
  rtOut
}

# running for samples
dirAll<-paste0("Result/",LETTERS[c(1,2,4,5)],"/1eny/Export")
proteinAll<-list.files(dirAll,pattern="protein", recursive = T,full.names = T)
peptideAll<-list.files(dirAll,pattern="peptide", recursive = T,full.names = T)

mattsv<-list.files(dirAll,pattern="protein", recursive = T)
nameOutAll<-paste0(gsub(".elib.proteins.txt","",mattsv),"_1eny")

for(m in 1:length(mattsv)){
  proIn<-read.delim(proteinAll[m],sep="\t",header = T,stringsAsFactors = F)
  pepIn<-read.delim(peptideAll[m],sep="\t",header = T,stringsAsFactors = F)
  
  FnamePro<-paste0("Matrix/",nameOutAll[m],"_ProInt.tsv")
  FnamePep<-paste0("Matrix/",nameOutAll[m],"_PepInt.tsv")
  
  if(!file.exists(FnamePep)){
    pepOut<-Eny2PepInt(pepIn)
    write.table(pepOut,FnamePep,sep="\t")
    print(FnamePep)
  }
  
  if(!file.exists(FnamePro)){
    proOut<-Eny2ProInt(proIn)
    write.table(proOut,FnamePro,sep="\t")
    print(FnamePro)
  }
  
}
