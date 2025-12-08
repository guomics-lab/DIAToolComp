packages <- c("shiny", "shinythemes", "dplyr")

# install.packages("vctrs")

for(pkg in packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
  }
}
library(shiny)
library(shinythemes)
library(dplyr)
options(shiny.maxRequestSize=60*1024^2)
load("template.RData")
ui <- bootstrapPage(
  navbarPage(theme = shinytheme("cosmo"), collapsible = TRUE,
             "DIA Tool Comparison", id="nav"
             ,tabPanel("Conversion",
                       sidebarLayout(
                         sidebarPanel(
                           tags$hr(),
                           radioButtons("tool", "Select DIA tool",
                                        choices = c(EncyclopeDIA = "eny",
                                                    DIANN = "dnn",
                                                    Skyline = "sky",
                                                    OpenSWATH = "osw",
                                                    Spectronaut = "spn"
                                                    ),
                                        selected = "dnn"),
                           tags$hr(),
                           uiOutput("Upload_file"),
                         ),
                         mainPanel(
                           fluidRow(
                             column(4, offset=1, ( downloadButton("downloadPre", "Download precursor matrix"))),
                           ),
                           tags$hr(),
                           fluidRow(
                             column(4, offset=1, ( downloadButton("downloadPep", "Download peptide matrix"))),
                           ),
                           tags$hr(),
                           fluidRow(
                             column(4, offset=1, ( downloadButton("downloadPro", "Download protein matrix"))),
                           ),
                         ),
                       )
             ),
             # tabPanel("Comparison",
             #   
             # ),
             tabPanel("Template",
                      fluidRow(
                        column(6, offset=1, ( actionButton("do", "Load sample data"))),
                      ),
                      tags$hr(),
                      fluidRow(
                        column(4, offset=1, ( downloadButton("downloadEnyPep", "Download EncyclopeDIA's peptide matrix sample"))),
                        column(4, offset=1, ( downloadButton("downloadEnyPro", "Download EncyclopeDIA's protein matrix sample"))),
                      ),
                      tags$hr(),
                      
                      fluidRow(
                        column(4, offset=1, ( downloadButton("downloadDnnPep", "Download DIANN's precursor matrix sample"))),
                        column(4, offset=1, ( downloadButton("downloadDnnPro", "Download DIANN's protein group matrix sample"))),
                      ),
                      tags$hr(),
                      
                      fluidRow(
                        column(4, offset=1, ( downloadButton("downloadOsw", "Download OpenSWATH matrix sample"))),
                      ),
                      tags$hr(),
                      
                      fluidRow(
                        column(4, offset=1, ( downloadButton("downloadSky", "Download Skyline matrix sample"))),
                      ),
                      tags$hr(),
                      
                      fluidRow(
                        column(4, offset=1, ( downloadButton("downloadSpn", "Download Spectronaut matrix sample"))),
                      ),
             )
	     )
)
server <- function(input, output) {
 
  Eny2PepInt<-function(pepIn){
    pepOut<-as.data.frame(pepIn[,4:ncol(pepIn)],stringsAsFactors=F)

    pepAll<-pepIn$Peptide
    pepAll<-gsub("C\\[\\+57\\.021464\\]","C(UniMod:4)",pepAll)
    pepAll<-gsub("M\\[\\+15\\.994915\\]","M(UniMod:35)",pepAll)

    # MSnames_sp<-sapply(colnames(pepOut),strsplit,"\\.")
    # MSnames<-as.vector(sapply(MSnames_sp, "[[", 1))
    # colnames(pepOut)<-MSnames
    row.names(pepOut)<-pepAll
    pepOut
  }
  
  Eny2ProInt<-function(proIn){
    proOut<-as.data.frame(proIn[,4:ncol(proIn)],stringsAsFactors=F)
    # MSnames_sp<-sapply(colnames(proIn)[4:ncol(proIn)],strsplit,"\\.")
    # MSnames<-as.vector(sapply(MSnames_sp, "[[", 1))
    # colnames(proOut)<-MSnames
    colnames(proOut)<-colnames(proIn)[4:ncol(proIn)]
    uniprotAll<-sapply(proIn$Protein,strsplit,"\\|")
    uniprotAll<-as.vector(sapply(uniprotAll, "[[", 2))
    row.names(proOut)<-uniprotAll
    proOut
  }
  
  Dnn2PreInt<-function(prIn){
    indPre<-which(colnames(prIn)=="Precursor.Id")
    ind_prot_row<-which(colnames(prIn)=="Protein.Group")
   
    preOut<-as.data.frame(prIn[,c(ind_prot_row,(indPre):ncol(prIn))],stringsAsFactors=F)#1+
    # row.names(preOut)<-prIn$Precursor.Id
    
    # MSnames_sp<-sapply(colnames(prIn)[(1+indPre):ncol(prIn)],strsplit,"\\.")
    # MSnames<-as.vector(sapply(MSnames_sp, "[[", 7))
    # colnames(preOut)<-MSnames
    colnames(preOut)<-colnames(prIn)[(indPre):ncol(prIn)]
    preOut
  }
  Dnn2PepInt<-function(prIn){
    # prIn<-read.csv("data-DnnPre2025-12-08.csv",header = T,stringsAsFactors = F, check.names=FALSE, sep=",",nrow=5000)
    # dim(prIn)
    indPre<-which(colnames(prIn)=="Precursor.Id")
    indseq<-which(colnames(prIn)=="Stripped.Sequence")
    # preOut<-as.data.frame(prIn[,(indPre):ncol(prIn)],stringsAsFactors=F)
    ind_mean<-prIn[,c(indseq,(indPre+1):ncol(prIn))]
    # 直接使用 group_by 和 summarise
    pepOut <- ind_mean %>%
      group_by(Stripped.Sequence) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE) ) %>%
      as.data.frame()
    pepOut$Protein.Group<-prIn$Protein.Group[match(pepOut$Stripped.Sequence,prIn$Stripped.Sequence)]
    
    pepOut<-pepOut[,c(ncol(pepOut),1:(ncol(pepOut)-1))]
    
    pepOut
  }
  
  Dnn2ProInt<-function(pgIn){
    
    # pgIn<-read.csv("data-DnnPro2025-12-08.csv",header = T,stringsAsFactors = F, check.names=FALSE, sep=",")
    indPro<-which(colnames(pgIn)=="Protein.Group")
    proOut<-as.data.frame(pgIn[,c(indPro,6:ncol(pgIn))],stringsAsFactors=F)
    colnames(proOut)<-colnames(pgIn)[c(indPro,6:ncol(pgIn))]
    proOut
  }
  
  Osw2PreInt<-function(readPepList){
    PepIntList<-list()
    for(m in 1:length(readPepList)){
      readIn<-readPepList[[m]]
      readIn<-readIn[readIn$decoy==0&readIn$m_score<0.01&readIn$peak_group_rank==1,]
      readIn<-unique(readIn)
      
      readPep<-unique(readIn[,grepl("FullPeptideName|Charge|Intensity" ,colnames(readIn))])
      readPep2<-unique(cbind(paste0(readPep$FullPeptideName,readPep$Charge),readPep$Intensity))
      readPepInt<-as.numeric(readPep2[,2])
      names(readPepInt)<-readPep2[,1]
      PepIntList[[m]]<-readPepInt
    }
    
    preAll<-unique(unlist(sapply(PepIntList,names)))
    dfOut<-data.frame(matrix(NA,nrow=length(preAll),ncol=length(PepIntList)),stringsAsFactors = F)
    row.names(dfOut)<-preAll
    
    for(m in 1:length(PepIntList)){
      dfOut[names(PepIntList[[m]]),m]<-PepIntList[[m]]
    }
    dfOut
  }
  StripeChargeOsw<-function(preIn){
    preIn_sp<-unlist(strsplit(preIn,""))
    paste0(preIn_sp[1:(length(preIn_sp)-1)],collapse = "")
  }
  
  Osw2PE<-function(readPepList){
    PepIntList<-list()
    for(m in 1:length(readPepList)){
      readIn<-readPepList[[m]]
      #remove decoy m score <0.01 p rank=1
      readIn<-readIn[readIn$decoy==0&readIn$m_score<0.01&readIn$peak_group_rank==1,]
      readIn<-unique(readIn)
      
      readPep<-unique(readIn[,grepl("FullPeptideName|Charge|Intensity|Protein" ,colnames(readIn))])
      readPep2<-unique(cbind(paste(paste(readPep$FullPeptideName,readPep$Charge,sep="_")
                                   ,readPep$ProteinName,sep="-")
                             ,readPep$Intensity))
      readPepInt<-as.numeric(readPep2[,2])
      names(readPepInt)<-readPep2[,1]
      PepIntList[[m]]<-readPepInt
    }
    
    preproAll<-unique(unlist(sapply(PepIntList,names)))
    peOut<-data.frame(matrix(NA,nrow=length(preproAll),ncol=length(PepIntList)),stringsAsFactors = F)
    row.names(peOut)<-preproAll
    
    for(m in 1:length(PepIntList)){
      peOut[names(PepIntList[[m]]),m]<-PepIntList[[m]]
    }
    
    preAll<-sapply(strsplit(preproAll,"-"),"[[",1)
    proAll<-sapply(strsplit(preproAll,"-"),"[[",2)
    
    PEOut<-cbind(proAll,peOut)
    row.names(PEOut)<-preAll
    uniAll<-sapply(sapply(PEOut$proAll,strsplit,"\\|"),"[[",1)
    PEOut$proAll<-uniAll
    PEOut
  }
  
  Pre2PepOsw<-function(PreIn){
    pepPre<-as.vector(sapply(row.names(PreIn),StripeChargeOsw))
    pepAll<-unique(pepPre)
    
    pepOut<-data.frame(matrix(NA,nrow=length(pepAll),ncol=ncol(PreIn)))
    colnames(pepOut)<-colnames(PreIn)
    row.names(pepOut)<-pepAll
    
    for(i in 1:length(pepAll)){
      indPep<-which(pepPre==pepAll[i])
      if(length(indPep)>1){
        df_preInSel<-as.data.frame(PreIn[indPep,],stringsAsFactors=F)
        meanAll<-c()
        for(j in 1:ncol(df_preInSel)){
          meanAll[j]<-mean(as.numeric(df_preInSel[,j]),na.rm=T)
        }
        pepOut[i,]<-meanAll
      }else{
        pepOut[i,]<-PreIn[indPep,]
      }
      # print(i)
    }
    pepOut[pepOut=="NaN"]<-NA
    pepOut
  }
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
  
  Pre2PepSky<-function(PreIn){
    pepPre<-as.vector(sapply(row.names(PreIn),StripeChargeSky))
    pepAll<-unique(pepPre)
    
    pepOut<-data.frame(matrix(NA,nrow=length(pepAll),ncol=ncol(PreIn)))
    colnames(pepOut)<-colnames(PreIn)
    row.names(pepOut)<-pepAll
    
    for(i in 1:length(pepAll)){#
      indPep<-which(pepPre==pepAll[i])
      if(length(indPep)>1){
        df_preInSel<-as.data.frame(PreIn[indPep,],stringsAsFactors=F)
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
  StripeChargeSky<-function(preIn){
    preIn_sp<-unlist(strsplit(preIn,""))
    paste0(preIn_sp[1:(length(preIn_sp)-1)],collapse = "")
  }
  
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
    PepNames<-as.vector(sapply(PepNames,StripePepSpn))
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
  StripePepSpn<-function(pepIn){
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
  
  # read
  datasetInput <- reactive({
    if(input$tool == "dnn") {
      prIn<-read.delim(input$file1$datapath, header = T, check.names=FALSE, sep=",")
      withProgress(message = 'Calculation in precurosr',
                   detail = 'This may take a while...', value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                       Sys.sleep(0.25)
                     }
                   })
      prOut<-Dnn2PreInt(prIn)
      withProgress(message = 'Calculation in peptide',
                   detail = 'This may take a while...', value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                       Sys.sleep(0.25)
                     }
                   })
      pepOut<-Dnn2PepInt(prIn)
      
      pgIn<-read.delim(input$file2$datapath,header = T,stringsAsFactors = F, check.names=FALSE, sep=",")
      
      proOut<-Dnn2ProInt(pgIn)
      
      list(prOut,pepOut,proOut)
    }else if(input$tool=="eny"){
      pepIn<-read.delim(input$file1$datapath,header = T,stringsAsFactors = F, check.names=FALSE, sep=",")
      proIn<-read.delim(input$file2$datapath,header = T,stringsAsFactors = F, check.names=FALSE, sep=",")
      
      pepOut<-Eny2PepInt(pepIn)
      proOut<-Eny2ProInt(proIn)
      
      prOut<-"No precursor matrix for this tool!"
      
      list(prOut,pepOut,proOut)
    }else if(input$tool=="osw"){
      resultFiles<-input$file1$datapath
      readPepList<-list()
      for(resultFile in resultFiles){
        readIn<-read.delim(resultFile,stringsAsFactors = F, check.names=FALSE, sep=",")
        readPepList[[which(resultFiles==resultFile)]]<-readIn
      }
      names(readPepList)<-resultFiles
      
      PreIn<-Osw2PreInt(readPepList)
      PepOut<-Pre2PepOsw(PreIn)
      proOut<-Osw2PE(readPepList)
      
      list(PreIn,PepOut,proOut)
    }else if(input$tool=="sky"){
      MatIn<-read.delim(input$file1$datapath,header = T
                        ,stringsAsFactors = F, check.names=FALSE, sep=",")
      PreIn<-Sky2PreInt(MatIn)
      PepOut<-Pre2PepSky(Sky2PreInt(MatIn))
      proOut<-Sky2PE(MatIn)
      proOut[proOut==0]<-NA
      
      list(PreIn,PepOut,proOut)
    }else if(input$tool=="spn"){
      MatIn<-read.delim(input$file1$datapath,header = T
                        ,stringsAsFactors = F, check.names=FALSE, sep=",")
      MatPre<-Spn2PreInt(MatIn)
      MatPep<-Spn2PepInt(MatIn)
      MatPro<-Spn2ProInt(MatIn)
      list(MatPre,MatPep,MatPro)
    }
  })
  
  
  output$Upload_file = renderUI({
    if(input$tool=="dnn"){
      tagList(
        fileInput("file1",label="precursor file",
                  multiple = F,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")
        ),
        fileInput("file2",label="protein group file",
                  multiple = F,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv"))
      )
    }else if(input$tool=="eny"){
      tagList(
        fileInput("file1",label="peptide file",
                  multiple = F,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")
        ),
        fileInput("file2",label="protein file",
                  multiple = F,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv"))
      )
    }else if(input$tool=="osw"){
      fileInput("file1",label="tsv files",
                multiple = T,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")
      )
    }else if(input$tool=="sky"){
      fileInput("file1",label="csv file",
                multiple = F,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")
      )
    }else if(input$tool=="spn"){
      fileInput("file1",label="xls file",
                multiple = F,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")
      )
    }
  })
  
  # download
  output$downloadPre <- downloadHandler(
    filename = function() {
      paste(input$tool, "_PreInt.tsv", sep = "")
    },
    content = function(file) {
      shiny::withProgress(
        message = paste0("Downloading", input$tool, "_PreInt.tsv Data"),
        value = 0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(5/10)
          
          dds = datasetInput()
          write.table(dds[[1]], file,sep="\t",row.names = F)
        }
      )
    }
  )
  
  output$downloadPep <- downloadHandler(
    filename = function() {
      paste(input$tool, "_PepInt.tsv", sep = "")
    },
    content = function(file) {
      shiny::withProgress(
        message = paste0("Downloading", input$tool, "_PepInt.tsv Data"),
        value = 0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(5/10)
          write.table(datasetInput()[[2]], file,sep="\t",row.names = F)
        }
      )
    }
  )
  
  output$downloadPro <- downloadHandler(
    filename = function() {
      paste(input$tool, "_ProInt.tsv", sep = "")
    },
    content = function(file) {
      shiny::withProgress(
        message = paste0("Downloading", input$tool, "_ProInt.tsv Data"),
        value = 0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(5/10)
          write.table(datasetInput()[[3]], file,sep="\t",row.names = F)
        }
      )
    }
  )
  
  output$downloadSky <- downloadHandler(
    filename = function() {
      paste("data-Sky", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(skyin, file,row.names = F)
    }
  )
  output$downloadSpn <- downloadHandler(
    filename = function() {
      paste("data-Spn", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(spnin, file,row.names = F)
    }
  )
  output$downloadOsw <- downloadHandler(
    filename = function() {
      paste("data-Osw", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(oswin, file,row.names = F)
    }
  )
  output$downloadEnyPep <- downloadHandler(
    filename = function() {
      paste("data-EnyPep", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(enypep, file,row.names = F)
    }
  )
  output$downloadEnyPro <- downloadHandler(
    filename = function() {
      paste("data-EnyPro", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(enypro, file,row.names = F)
    }
  )
  output$downloadDnnPep<- downloadHandler(
    filename = function() {
      paste("data-DnnPre", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(dnnpr, file,row.names = F)
    }
  )
  output$downloadDnnPro<- downloadHandler(
    filename = function() {
      paste("data-DnnPro", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(dnnpg, file,row.names = F)
    }
  )
  observeEvent(input$do, {
    shiny::withProgress(
      message = paste0("Loading sample datasest"),
      value = 0,
      {
        shiny::incProgress(1/10)
        shiny::incProgress(5/10)
        
      }
    )
  }
  )
}
shinyApp(ui,server)
