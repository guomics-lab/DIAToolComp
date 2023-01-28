libNames<-c("01A","01B","01C","02A","02B","02C")
fileAll<-list.files(path = '../spn_lib6/',pattern = "xls$")
libAll<-list()
for(m in 1:length(libNames)){
  readin<-read.delim(paste0("../spn_lib6/",fileAll[m]),stringsAsFactors = F,sep='\t')[1:10,]
  libAll[[m]]<-readin
}

# m=2
# write.csv(fileAll[m],paste0("spnCSV",fileAll[m],".csv"))

library(readr)
library(plyr)
library(readxl)
library(stringr)
library(magrittr)

for(m in 6:6){
  df0 <- read.delim(paste0("../spn_lib6/",fileAll[m]),stringsAsFactors = F,sep='\t')
  df<-df0[!as.logical(df0$ExcludeFromAssay),]
  df1 <- df[1]
  
  df1$PrecursorMz <- df$PrecursorMz
  df1$ProductMz <- df$FragmentMz
  df1$PrecursorCharge <- df$PrecursorCharge
  df1$ProductChange <- df$FragmentCharge
  df1$LibraryIntensity <- df$RelativeIntensity
  df1$NormalizedRetentionTime <- df$iRT
  df1$PeptideSequence <- df$StrippedPeptide
  
  if(grepl("_",df$ModifiedPeptide[1])){
    df1$ModifiedPeptideSequence <- sapply(df$ModifiedPeptide,function(x){str_split(x,"_")[[1]][2]})
  }  else{
    df1$ModifiedPeptideSequence <- df$ModifiedPeptide 
  }

  
  df1$ModifiedPeptideSequence <- gsub("\\[Carbamidomethyl \\(C\\)\\]", "(UniMod:4)",df1$ModifiedPeptideSequence)   
  df1$ModifiedPeptideSequence <- gsub(pattern="\\[Oxidation \\(M\\)\\]","(UniMod:35)",df1$ModifiedPeptideSequence)
  df1$ModifiedPeptideSequence <- gsub(pattern="\\[Acetyl \\(Protein N-term\\)\\]","(UniMod:1)",df1$ModifiedPeptideSequence)
  
  df1$PeptideGroupLabel <- paste0(df1$ModifiedPeptideSequence,df$PrecursorCharge)
  df1$ProteinId <- paste0(df$ProteinGroups,"|",df$Protein.Name)
  df1$UniprotId <- paste0(df$UniProtIds,"|",df$Protein.Name)
  df1$FragmentType <- df$FragmentType
  df1$FragmentSeriesNumber <- df$FragmentNumber
  df1$CollisionEnergy <- rep(-1,nrow(df))
  df1$PrecursorIonMobility <- rep(-1,nrow(df))
  if(sum(grepl("IonMobility",colnames(df)))>0){
    if(df$IonMobility[1]!='NaN'){
      df1$PrecursorIonMobility <- df$IonMobility
    } 
  }
  df1$TransitionGroupId <- df1$PeptideGroupLabel
  df1$Decoy <- rep(0,nrow(df))
  d <- grep("(UniMod:1)",df1$ModifiedPeptideSequence)
  if(length(d)>0){
    df1 <- df1[-d,]
  }
  df10<-df1
  # write_tsv(df1,"Testis_Library_oldos_geTOP6.tsv")
  
  df <- df10
  nm <- sapply(df$PeptideSequence, function(x){str_length(x)})
  df1 <- df[nm>7,]
  df1$PeptideGroupLabel <- sapply(df1$PeptideGroupLabel,function(x){paste(str_sub(x,1,-2),str_sub(x,-1,-1),sep  ="_")})
  df1$TransitionGroupId <- df1$PeptideGroupLabel
  scannm <- c(rep(1:length(table(df1$PeptideGroupLabel)),as.numeric(table(df1$PeptideGroupLabel))))
  df2 <- df1[order(df1$PeptideGroupLabel),]
  df2$PeptideGroupLabel <- paste0(scannm,"_",df2$PeptideGroupLabel)
  df2$TransitionGroupId <- paste0(scannm,"_",df2$TransitionGroupId)
  write_tsv(df2,paste0("oldos_delect7_geTOP6",fileAll[m],".tsv"))
  
  df2$PrecursorMz<-gsub(",",".",df2$PrecursorMz)
  df2$ProductMz<-gsub(",",".",df2$ProductMz)
  
}

