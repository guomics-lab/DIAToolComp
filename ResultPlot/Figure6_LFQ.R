fileRat<-list.files(path = "../../Matrix/",pattern = "Ratio.tsv")
fileRat<-fileRat[grepl("A|B|C",fileRat)]
fileCond<-c()
for(i in 1:length(fileRat)){
  fileCond[i]<-paste(unlist(strsplit(fileRat[i],"_"))[1:4],collapse =  "_")
}
fileCondAll<-unique(fileCond)

colorAll<-c("#2697c1","#eb9396","#a5c48d","#d4c41d","#eb9396")
names(colorAll)<-c("HUMAN","YEAS8","ECOLI","CAEEL","YEAST")
LFQPlot<-function(MatIn,tit="LFQ"){
  
  specNames<-unique(MatIn$specAll)
  specNames<-specNames[!is.na(specNames)]
  
  
  #log10B vs log2(A/B)
  plot(log10(MatIn$Bmean),log2(MatIn$ABratio)
       ,col=colorAll[MatIn$specAll]
       ,xlab="log10 B",ylab = "log2 A/B",cex.lab=2
       ,main=tit,cex.axis=2,cex.main=2)

}

pdf("LFQ20210802.pdf",height = 6,width=6)
par(mar=c(4,6,3,1))
for(i in 1:length(fileCondAll)){
  fileToRead<-fileRat[grepl(fileCondAll[i],fileCond)]
  
  readIn<-read.table(paste0("../../Matrix/",fileToRead[1]),stringsAsFactors = F,row.names = NULL,header = T)
  LFQPlot(readIn,paste0(fileCondAll[i],"_Pep"))
  readIn<-read.table(paste0("../../Matrix/",fileToRead[2]),stringsAsFactors = F,row.names = NULL,header = T)
  LFQPlot(readIn,paste0(fileCondAll[i],"_Pro"))
}
dev.off()
