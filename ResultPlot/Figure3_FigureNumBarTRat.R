#TrueList processing
#ABC
ratPepAll<-list.files("../../Matrix/",pattern = "PepRatio.tsv")
ratProAll<-list.files("../../Matrix/",pattern = "ProRatio.tsv")

ATratio<-c(1,2,0.25);names(ATratio)<-c("HUMAN","YEAST","ECOLI")
BTratio<-c(1,1/1.3,2);names(BTratio)<-c("HUMAN","CAEEL","YEAST")
CTratio<-c(1,1/3);names(CTratio)<-c("HUMAN","YEAST")

Tratio<-list(ATratio,BTratio,CTratio)
names(Tratio)<-LETTERS[1:3]

for(condName in LETTERS[3]){
  ratPepSel<-ratPepAll[grepl(condName,ratPepAll)]
  ratProSel<-ratProAll[grepl(condName,ratProAll)]
  
  for(m in 1:length(ratPepSel)){
    if(!file.exists(paste0('../../Matrix/',gsub("Ratio","RatioT",ratPepSel[m])))){
      ratPep<-read.table(paste0("../../Matrix/",ratPepSel[m]),stringsAsFactors = F,row.names =NULL)
      ratPro<-read.table(paste0("../../Matrix/",ratProSel[m]),stringsAsFactors = F,row.names =NULL)
      #peptide
      TraioSel<-Tratio[[condName]]
      TraioSelUp<-TraioSel*1.3
      TraioSelDw<-TraioSel*0.7
      
      ratPep$TF<-NA
      for(i in 1:nrow(ratPep)){
        if(!is.na(ratPep[i,5])&!(is.na(ratPep$ABratio[i]))&ratPep[i,5]%in%names(TraioSel)){
          if(ratPep$ABratio[i]<TraioSelUp[ratPep[i,5]]&ratPep$ABratio[i]>TraioSelDw[ratPep[i,5]]){
            ratPep$TF[i]<-T
          }else{
            ratPep$TF[i]<-F
          }
        }
      }
      
      #protein
      ratPro$TF<-NA
      for(i in 1:nrow(ratPro)){
        if(!is.na(ratPro$specAll[i])&!(is.na(ratPro$ABratio[i]))&ratPro[i,5]%in%names(TraioSel)){
          if(ratPro$ABratio[i]<TraioSelUp[ratPro[i,5]]&ratPro$ABratio[i]>TraioSelDw[ratPro[i,5]]){
            ratPro$TF[i]<-T
          }else{
            ratPro$TF[i]<-F
          }
        }
      }
      write.table(ratPep,paste0('../../Matrix/',gsub("Ratio","RatioT",ratPepSel[m])),sep="\t",row.names = F)
      write.table(ratPro,paste0('../../Matrix/',gsub("Ratio","RatioT",ratProSel[m])),sep="\t",row.names = F)
      print(paste0('../../Matrix/',gsub("Ratio","RatioT",ratPepSel[m])))
      print(paste0('../../Matrix/',gsub("Ratio","RatioT",ratProSel[m])))
    }
  }
}

#DEF
ratPepAll<-list.files("../../Matrix/",pattern = "PepRatio.tsv")
ratProAll<-list.files("../../Matrix/",pattern = "ProRatio.tsv")

for(condName in LETTERS[4:6]){
  ratPepSel<-ratPepAll[grepl(condName,ratPepAll)]
  ratProSel<-ratProAll[grepl(condName,ratProAll)]
  
  for(m in 1:length(ratPepSel)){
    # if(!file.exists(paste0('../../Matrix/',gsub("Ratio","RatioT",ratPepSel[m])))){
      ratPep<-read.table(paste0("../../Matrix/",ratPepSel[m]),stringsAsFactors = F,row.names =NULL)
      ratPro<-read.table(paste0("../../Matrix/",ratProSel[m]),stringsAsFactors = F,row.names =NULL)
      #peptide
      ratPep$TF<-NA
      for(i in 1:nrow(ratPep)){
        if(!is.na(ratPep$ABratio[i])&!is.na(ratPep$Amean[i])){
          ratPep$TF[i]<-T
        }
      }
      
      #protein
      ratPro$TF<-NA
      for(i in 1:nrow(ratPro)){
        if(!is.na(ratPro$ABratio[i])&!is.na(ratPro$Amean[i])){
          ratPro$TF[i]<-T
        }
      }
      write.table(ratPep,paste0('../../Matrix/',gsub("Ratio","RatioT",ratPepSel[m])),sep="\t",row.names = F)
      write.table(ratPro,paste0('../../Matrix/',gsub("Ratio","RatioT",ratProSel[m])),sep="\t",row.names = F)
    # }
  }
}


colorLabel<-c("#2697c1","#eb9396","#a5c48d","#2ab58d","#d4c41d")
names(colorLabel)<-c("1eny","2dnn","3osw","4sky","5spn")

pdf("PepProNum20210726.pdf",width = 16,height = 12)
# pdf("PepNum20210726.pdf",width = 8,height = 12)
par(mfcol=c(3,4),lwd = 4)
#peptide
tsvAll<-list.files("../../Matrix/",pattern = "PepInt.tsv")
for(DataID in LETTERS[1:6]){
  tsvSel<-tsvAll[grepl(DataID,tsvAll)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  ratSel<-gsub("Int","RatioT",tsvSel)
  
  pepTList<-pepList<-list()
  for(i in 1:length(tsvSel)){
    pepIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F,header = T,row.names = 1)
    ratIn<-read.table(paste0("../../Matrix/",ratSel[i]),stringsAsFactors = F,header = T,row.names = NULL)
    
    pepOut<-as.data.frame(pepIn[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn)),]
                          ,StringAsFactors=F)
    colnames(pepOut)<-colnames(pepIn)
    row.names(pepOut)<-row.names(pepIn)[which(apply(is.na(pepIn),1,sum)!=ncol(pepIn))]
    pepList[[i]]<-row.names(pepOut)
    pepTList[[i]]<-row.names(pepOut)[row.names(pepOut)%in%ratIn$row.names[which(ratIn$TF)]]
  }
 
  
  #color told by names
  toolNames<-sapply(strsplit(tsvSel,"_"),"[[",4)
  colorOut<-colorIn<-colorLabel[toolNames]
  libOrfree<-sapply(strsplit(tsvSel,"_"),"[[",2)
  colorIn[libOrfree=="free"]<-NA
  colorOut[libOrfree!="free"]<-NA
  
  # number to plot
  numToPlot<-sapply(pepList,length)
  numToPlotT<-sapply(pepTList,length)
  
  Od<-order(numToPlot)
  
  # to plot
  b1<-barplot(numToPlot[Od]
          ,col=colorIn[Od]
          ,border = colorOut[Od],space=0.2
          ,cex.axis = 2,cex.lab=2,cex.main=2
          ,horiz = T
          ,xlab = "# peptides"
          ,main = DataID
          )
  colorShadow<-colorOut[Od]
  colorShadow[is.na(colorShadow)]<-"#FCFAFA"
  borderShadow<-colorOut[Od]
  borderShadow[is.na(borderShadow)]<-"#FCFAFA"
  
  barplot(numToPlotT[Od],axes=F
          ,col=colorShadow,space=0.2,density = 5
          ,border =borderShadow ,horiz = T
          ,add=T)
  text(numToPlot[Od]/2,as.vector(b1),numToPlot[Od],cex=2)
}
# dev.off()
#protein
# pdf("ProNum20210726.pdf",width = 8,height = 12)
# par(mfcol=c(3,2),lwd = 4)
tsvAll<-list.files("../../Matrix/",pattern = "ProInt.tsv")
for(DataID in LETTERS[1:6]){
  tsvSel<-tsvAll[grepl(DataID,tsvAll)]
  tsvSel<-tsvSel[grepl("S",tsvSel)]
  ratSel<-gsub("Int","RatioT",tsvSel)
  
  proTList<-proList<-list()
  for(i in 1:length(tsvSel)){
    proIn<-read.table(paste0("../../Matrix/",tsvSel[i]),stringsAsFactors = F,sep="\t",header = T,row.names = 1)
    ratIn<-read.table(paste0("../../Matrix/",ratSel[i]),stringsAsFactors = F,header = T,row.names = NULL)
    
    proOut<-as.data.frame(proIn[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn))),]#remove protein groups
                          ,StringAsFactors=F)
    colnames(proOut)<-colnames(proIn)
    row.names(proOut)<-row.names(proIn)[which(apply(is.na(proIn),1,sum)!=ncol(proIn)&!grepl(";",row.names(proIn)))]
    proList[[i]]<-row.names(proOut)
    proTList[[i]]<-row.names(proOut)[row.names(proOut)%in%ratIn$row.names[which(ratIn$TF)]]
  }

  #color told by names
  toolNames<-sapply(strsplit(tsvSel,"_"),"[[",4)
  colorOut<-colorIn<-colorLabel[toolNames]
  libOrfree<-sapply(strsplit(tsvSel,"_"),"[[",2)
  colorIn[libOrfree=="free"]<-NA
  colorOut[libOrfree!="free"]<-NA
  
  # number to plot
  numToPlot<-sapply(proList,length)
  numToPlotT<-sapply(proTList,length)
  
  Od<-order(numToPlot)
  
  # to plot
  b1<-barplot(numToPlot[Od]
              ,col=colorIn[Od]
              ,border = colorOut[Od],space=0.2
              ,cex.axis = 2,cex.lab=2,cex.main=2
              ,horiz = T
              ,xlab = "# proteins"
              ,main = DataID
  )
  colorShadow<-colorOut[Od]
  colorShadow[is.na(colorShadow)]<-"#FCFAFA"
  borderShadow<-colorOut[Od]
  borderShadow[is.na(borderShadow)]<-"#FCFAFA"
  barplot(numToPlotT[Od],axes=F
          ,col=colorShadow,space=0.2,density = 6
          ,border =borderShadow,horiz = T
          ,add=T)
  text(numToPlot[Od]/2,as.vector(b1),numToPlot[Od],cex=2)
}
dev.off()

# pdf("Legend.pdf",width = 15,height = 15)
# plot.new()
# colorOut<-colorIn<-colorLabel
# legendIn<-c("#2697c1","#eb9396","#a5c48d","#2ab58d","#d4c41d",NA,NA,NA,NA,NA)
# legendOut<-c(NA,NA,NA,NA,NA,"#eb9396","#d4c41d",NA,NA,NA)
# legend("center",legend=c("EncyclopeDIA","DIA-NN","OpenSWATH","Skyline","Spectronaut"
#                          ,"DIA-NN lib-free","Spectronaut lib-free","","","")
#        ,fill=legendIn,border=legendOut
#        ,cex=2,bty="n",ncol=2 )
# dev.off()