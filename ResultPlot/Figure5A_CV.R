library(vioplot)
colorLabel<-c("#2697c1","#eb9396","#a5c48d","#2ab58d","#d4c41d")
names(colorLabel)<-c("1eny","2dnn","3osw","4sky","5spn")

CoV<-function(vectorIn){
  vectorIn<-as.numeric(vectorIn)
  sd(vectorIn,na.rm=T)/mean(vectorIn,na.rm=T)
}

for(fileCond in c("PepInt","ProInt","PreInt")){#,"PepInt","ProInt"
  filePre<-list.files(path = '../../Matrix/',pattern = fileCond,full.names = T)
  filePre<-filePre[grepl("_S_",filePre)]
  
  pdf(paste0("Figure3C_QunatCV_",fileCond,"_FZ20210729.pdf"),height = 20,width = 10)
  par(mfrow=c(5,1))
  for(condName in LETTERS[1:5]){
    filePreSel<-filePre[grepl(condName,filePre)]
    fileNames<-sapply(strsplit(filePreSel, "\\/"), function(x) x[5], simplify=T)
    selNames<-paste(sapply(strsplit(fileNames, "_")
                           , function(x) x[c(2)], simplify=T),
                    sapply(strsplit(fileNames, "_"), function(x) x[c(4)], simplify=T),sep="_")
    toolNames<-sapply(strsplit(selNames,"_"),"[[",2)
    libOrfree<-sapply(strsplit(selNames,"_"),"[[",1)
    
    colorOut<-colorIn<-colorLabel[toolNames]
    colorIn[libOrfree=="free"]<-0
    colorOut[libOrfree!="free"]<-0
    
    CVlist<-list()
    if(condName=="A"){
      for(m in 1:length(filePreSel)){
        readIn<-read.table(filePreSel[m],stringsAsFactors = F,sep="\t",header = T,row.names = 1)
        CVall<-c(apply(readIn[,c(1,3,5)],1,CoV),
                 apply(readIn[,c(2,4,6)],1,CoV))
        CVlist[[m]]<-CVall[!is.na(CVall)]
        print(paste(condName,m,fileCond))
        # print(paste(colnames(readIn),row.names(readIn)[1]))
      }
    }else if(condName%in%c("B","C")){
      for(m in 1:length(filePreSel)){
        readIn<-read.table(filePreSel[m],stringsAsFactors = F,sep="\t",header = T,row.names = 1)
        CVall<-c(apply(readIn[,c(1:3)],1,CoV),
                       apply(readIn[,c(4:6)],1,CoV))
        CVlist[[m]]<-CVall[!is.na(CVall)]
        print(paste(condName,m,fileCond))
        # print(paste(colnames(readIn),row.names(readIn)[1]))
      }
    }else{
      for(m in 1:length(filePreSel)){
        readIn<-read.table(filePreSel[m],stringsAsFactors = F,sep="\t",header = T,row.names = 1)
        CVall<-c(apply(readIn,1,CoV))
        CVlist[[m]]<-CVall[!is.na(CVall)]
        print(paste(condName,m,fileCond))
        # print(paste(colnames(readIn),row.names(readIn)[1]))
      }
    }
    names(CVlist)<-selNames
    CVlist2<-CVlist[sapply(CVlist,length)>0]
    vioplot(CVlist2,main=paste(fileCond),las=2, plotCentre = "line",col=colorIn
                   ,border=colorOut,areaEqual=F,drawRect=T)
  }
  dev.off()
}

