PEoutAllF<-list.files(pattern = "_PEout",path="./Matrix20210724",full.names = T)
for(i in 1:length(PEoutAllF)){
  FnameOutF<-gsub("PEout","ProInt",PEoutAllF[i])
  if(!file.exists(FnameOutF)){
    PEIn<-read.table(PEoutAllF[i],header = T,row.names = 1
                     ,stringsAsFactors = F,sep = "\t")
    FileOut<-round(2^PEIn,0)
    write.table(FileOut,FnameOutF,sep="\t")
    print(FnameOutF)
  }
}
