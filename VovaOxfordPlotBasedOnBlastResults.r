rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

#install.packages("seqinr")
library(seqinr)
library(plotrix)

#read contig names and lengths from fasta file
contigNamesAndLengths_fromFastaFile <- function(sFileFasta){
  if(FALSE){
    con = file(sFileFasta, "r")
    while ( TRUE ) {
      line = readLines(con, n = 1)
      if ( length(line) == 0 ) {
        break
      }
      print(line)
      break
    }
    close(con)
  }
  
  #sFileFasta="C:/Frenkel/LTCPython/VovaPy/20211004/Cnig_gn3.1/Cnig_gn3.1.fasta"
  myData <- read.fasta(sFileFasta) 
  nContigs=length(myData)
  datalist = list()
  i=1
  coorStart=0
  for (contig in myData){
    name=attr(contig, "name")
    n=length(contig)
    #g=list(name,n)
    dat <- data.frame(name = name, len = n, coorStart=coorStart)
    dat$i <- i  # maybe you want to keep track of which iteration produced it?
    datalist[[i]] <- dat
    i<-i+1
    coorStart<-coorStart+n
  }
  df = do.call(rbind, datalist)
  return(df)
}


#graph

startMyPlotOxford<-function(contigsX,sx="x",sy="",colorOfPoint="black",ycoef=2,xcoef=1,xcoefMin=0,ycoefMin=0){
  fontSize=2#0
  nCtg1=nrow(contigsX)
  xmax=contigsX$coorStart[nCtg1]+contigsX$len[nCtg1]
  xlim=c(xmax*xcoef*xcoefMin,xmax*xcoef)
  ylim=c(xmax*xcoef*ycoef*ycoefMin,xmax*xcoef*ycoef)#1.5), 
  plot(
    #no need this point (need only to start graph)
    c(2*xmax),#vx 
    c(2*xmax), #vy
    #
    # log="x", #axis x in logarithmic scale 
    xlab=sx , #caption of axis x
    ylab=sy, #caption of axis y
    xaxt="n", #no caption foir axis x (will be added manually)
    col=colorOfPoint,#"yellow", #color of points
    xlim=xlim, #xMin,xMax
    ylim=ylim, #xMin,xMax
    #cex=2,#scaling of font size
    #cex.lab=fontSize/10,#scaling of font size (lable of axes)
    #cex.axis=fontSize/10,#scaling of font size (numbers)
    #mar=c(5.1,8.1,4.1,2.1),#margines:bootom, left,top,right
  )
  return(list(xlim,ylim))
}

addGridBasedOnContigs<-function(contigs,sColor = "darkgreen",bHorizontal=TRUE,lenMin=1000000){
  nCtg=nrow(contigs)
  for (iCtg in 1:nCtg){
    x<-contigs[iCtg,"coorStart"]
    l<-contigs[iCtg,"len"]
    if(l>=lenMin){
      if(bHorizontal){
        #abline(h = x, col = sColor)
        segments(0, x, x1 = 1000000000, col = sColor)#, y1 = y0)
      }else{
        #abline(v = x, col = sColor)
        segments(x, 0, y1 = 1000000000, col = sColor)
      }
    }else{
      break
    }
  }
}

addPointsFromdf<-function(df,sx,sy,sColor="black"){
  vx<-df[,sx]
  vy<-df[,sy]
  #points(x=vx, y=vy, type = "p", pch=".",col=sColor)
  points(x=vx, y=vy, type = "p", pch=1,cex=0.01,col=sColor)
}
addPointsFromMarkers<-function(dfMarkers,sMap_x,sMap_y,sColor="black"){
  sx<-paste("x_",sMap_x, sep = "")
  sy<-paste("x_",sMap_y, sep = "")
  addPointsFromdf(dfMarkers,sx,sy,sColor)
}






addPointsByBlast<-function(contigsX,contigsY,sFileBlast){
  
}

addCoorColumnToMarkers<-function(dfMarkers,sMap="Cnig_gn3.1",contigs,bType=TRUE){
  #c<-dfMarkers["Cnig_gn3.1_pos"]
  #sMap="Cnig_gn3.1"
  #contigs<-contigs1
  sShapka_part<-paste(sMap,"_part", sep = "")
  sShapka_pos<-paste(sMap,"_pos", sep = "")
  sShapka_n<-paste(sMap,"_iType", sep = "")
  sx<-paste("x_",sMap, sep = "")
  print(sx)
  print(sMap)
  #print()
  
  nm_mm<-nrow(dfMarkers)
  dfMarkers[sx]<-NA
  #vx=list()
  #vy=list()
  for (im in 1:nm_mm){
    i<-0
    if (bType){
      i<-dfMarkers[im,sShapka_n]
    }
    if(!(is.na(i))){
      if(i==0){
        sCtg<-dfMarkers[im,sShapka_part]
        iCtg<-which(contigs$name == sCtg)
        x0<-contigs[iCtg,"coorStart"]
        x<-dfMarkers[im,sShapka_pos]
        if(FALSE){
          print(paste("im=",im, sep = ""))
          print(paste("i=",i, sep = ""))
          print(paste("sCtg=",sCtg, sep = ""))
          print(paste("iCtg=",iCtg, sep = ""))
          print(paste("x0=",x0, sep = ""))
          print(paste("x=",x, sep = ""))
        }
        dfMarkers[im,sx]=x+x0
      }
    }
  }
  print("privet")
  return(dfMarkers)
}

PrivmanLab_makeOxfordPlots<-function(){
  
}

oxfordPlot_compareGenomesByMarkers<-function(){
  sFileFasta1_Cnig_gn3.1="C:/Frenkel/LTCPython/VovaPy/20211004/Cnig_gn3.1/Cnig_gn3.1.fasta"
  sFileFasta2_Cnig_gn3a.1="C:/Frenkel/LTCPython/VovaPy/20211004/Cnig_gn3a.1/Cnig_gn3a.1.fasta"
  sFileFasta3_E="C:/Frenkel/Privman/LabMeetingPaperDiscussion/Darras_et_al_2022_CthaglipsisGenome/JAJUXE01.1.fsa_nt"
  sFileFasta4_C="C:/Frenkel/Privman/LabMeetingPaperDiscussion/Darras_et_al_2022_CthaglipsisGenome/JAJUXC01.1.fsa_nt"
  sFileFasta5_Formica="C:/Frenkel/Privman/Cnig_gn1/formica_selysi/genome_assemblies_genome_fasta/ncbi-genomes-2020-05-04/GCA_009859135.1_ASM985913v1_genomic.fna"
  sFileFasta6_SolenopsisInvicta_ref="C:/Frenkel/Privman/NanoporeSequencingSupergene/TrainingData/referenceGenome/JAEQMK01.1.fsa_nt/JAEQMK01.1.fsa_nt"
  sFileFasta6_SolenopsisInvicta_b="C:/Frenkel/Privman/NanoporeSequencingSupergene/TrainingData/VDGK01.1.fsa_nt"
  sFileFasta6_SolenopsisInvicta_bb="C:/Frenkel/Privman/NanoporeSequencingSupergene/TrainingData/VDGJ01.1.fsa_nt"
  
  
  contigs1_Cnig_gn3.1 <- contigNamesAndLengths_fromFastaFile(sFileFasta1_Cnig_gn3.1)
  #
  contigs3_JAJUXE <- contigNamesAndLengths_fromFastaFile(sFileFasta3_E)
  contigs4_JAJUXC <- contigNamesAndLengths_fromFastaFile(sFileFasta4_C)
  contigs5_Formica <- contigNamesAndLengths_fromFastaFile(sFileFasta5_Formica)
  contigs6_SolenopsisInvicta <- contigNamesAndLengths_fromFastaFile(sFileFasta6_SolenopsisInvicta_ref)
  contigs6_b <- contigNamesAndLengths_fromFastaFile(sFileFasta6_SolenopsisInvicta_b)
  contigs6_bb <- contigNamesAndLengths_fromFastaFile(sFileFasta6_SolenopsisInvicta_bb)
  #n=nchar(c)
  #n=length(b)
  #name=attr(b, "name")
  
  sFileMarkers="C:/Frenkel/LTCPython/VovaPy/MarkersWithMultiplePos_JAJUXC_JAJUXE.txt"
  dfMarkers<-read.csv(sFileMarkers, header = TRUE, sep = "\t", quote = "\"",
                      dec = ".", fill = TRUE, comment.char = "#")
  dfMarkers<-addCoorColumnToMarkers(dfMarkers,sMap="Cnig_gn3.1",contigs1_Cnig_gn3.1)
  dfMarkers<-addCoorColumnToMarkers(dfMarkers,sMap="JAJUXE",contigs3_JAJUXE)
  dfMarkers<-addCoorColumnToMarkers(dfMarkers,sMap="JAJUXC",contigs4_JAJUXC) 
  dfMarkers<-addCoorColumnToMarkers(dfMarkers,sMap="Formica",contigs5_Formica)
  dfMarkers<-addCoorColumnToMarkers(dfMarkers,sMap="Solenopsis_invicta_ref",contigs6_SolenopsisInvicta)
  
  
  #"Solenopsis_invicta_ref" vs Formica
  if(FALSE){
    #"Solenopsis_invicta_ref" vs Formica
    #oxford_Formica_vs_Solenopsis_ref.pdf
    startMyPlotOxford(contigs6_SolenopsisInvicta,sx="Solenopsis_invicta_ref",sy="Formica",ycoef=0.9)
    addGridBasedOnContigs(contigs6_SolenopsisInvicta,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
    addGridBasedOnContigs(contigs5_Formica,sColor = "darkgreen",bHorizontal=TRUE,lenMin=2000000)
    addPointsFromMarkers(dfMarkers,"Solenopsis_invicta_ref","Formica",sColor="black")
  }
  
  #Cnig_gn3.1 vs Formica
  if(FALSE){
    #Cnig_gn3.1 vs Formica
    #oxford_Formica_vs_Cnig_gn3.1.pdf
    startMyPlotOxford(contigs1_Cnig_gn3.1,sx="Cnig_gn3.1",sy="Formica",ycoef=1.15)
    addGridBasedOnContigs(contigs1_Cnig_gn3.1,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
    addGridBasedOnContigs(contigs5_Formica,sColor = "darkgreen",bHorizontal=TRUE,lenMin=2000000)
    addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","Formica",sColor="black")
  }
  
  #JAJUXE vs JAJUXC (Cathagliphis 2 by Darras)
  if(FALSE){
    #JAJUXE vs JAJUXC (Cathagliphis 2 by Darras)
    #oxford_JAJUXC_vs_JAJUXE____twoCathagliphisIspanica_Darras.pdf
    startMyPlotOxford(contigs3_JAJUXE,sx="JAJUXE",sy="JAJUXC",ycoef=1)
    addGridBasedOnContigs(contigs3_JAJUXE,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
    addGridBasedOnContigs(contigs4_JAJUXC,sColor = "darkgreen",bHorizontal=TRUE,lenMin=2000000)
    addPointsFromMarkers(dfMarkers,"JAJUXE","JAJUXC",sColor="black")
    
  }
  
  #Cnig_gn3.1 vs JAJUXC,JAJUXE
  if(FALSE){
    #Cnig_gn3.1 vs JAJUXC,JAJUXE
    #oxford_twoCathagliphisIspanica_Darras_JAJUXC_and_JAJUXE_vs_Cnig_gn3.1.pdf
    ll=startMyPlotOxford(contigs1_Cnig_gn3.1,sx="Cnig_gn3.1",sy="Cathagliphis Ispanica",ycoef=1)
    xlim=ll[[1]]
    ylim=ll[[2]]
    addGridBasedOnContigs(contigs1_Cnig_gn3.1,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
    #addGridBasedOnContigs(contigs5_Formica,sColor = "darkgreen",bHorizontal=TRUE,lenMin=2000000)
    #addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","Formica",sColor="black")
    addGridBasedOnContigs(contigs4_JAJUXC,sColor = "black",bHorizontal=TRUE,lenMin=2000000)
    addGridBasedOnContigs(contigs3_JAJUXE,sColor = "red",bHorizontal=TRUE,lenMin=2000000)
    #
    addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","JAJUXC",sColor="black")
    addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","JAJUXE",sColor="red")
    
    vs=c("JAJUXC","JAJUXE")
    vtc=c("black","red")
    

    px=0.7
    xLegend=xlim[1]*(1-px)+xlim[2]*px
    py=1.03
    yLegend=ylim[1]*(1-py)+ylim[2]*py
    legend(xLegend,yLegend,#coordinates of left top corner 
           legend=vs,#captions
           box.lty=0,
           bg="transparent",
           text.col=vtc,#"red",
           #col=vColor, #colors
           #lty=vLineType,#types of line 
           #cex.lab=20,#cex=2,#fontSize,
           #box.lty=1,#(0 - no border of legend)
           #lwd=vwidthOfLines
    )

  }
  
  #Cnig_gn3.1 vs JAJUXC,JAJUXE, Formica,Solinopsis
  if(FALSE){
    #Cnig_gn3.1 vs JAJUXC,JAJUXE, Formica,Solinopsis
    #oxford_twoCathagliphisIspanica_Darras_JAJUXC_and_JAJUXE_and_Formica_and_Solinopsis_vs_Cnig_gn3.1.pdf
    ll=startMyPlotOxford(contigs1_Cnig_gn3.1,sx="Cnig_gn3.1",sy="",ycoef=2)
    xlim=ll[[1]]
    ylim=ll[[2]]
    addGridBasedOnContigs(contigs1_Cnig_gn3.1,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
    #addGridBasedOnContigs(contigs5_Formica,sColor = "darkgreen",bHorizontal=TRUE,lenMin=2000000)
    #addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","Formica",sColor="black")
    addGridBasedOnContigs(contigs4_JAJUXC,sColor = "cyan",bHorizontal=TRUE,lenMin=2000000)
    addGridBasedOnContigs(contigs3_JAJUXE,sColor = "green",bHorizontal=TRUE,lenMin=2000000)
    addGridBasedOnContigs(contigs5_Formica,sColor = "red",bHorizontal=TRUE,lenMin=2000000)
    addGridBasedOnContigs(contigs6_SolenopsisInvicta,sColor = "black",bHorizontal=TRUE,lenMin=1000000)
    #
    addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","JAJUXC",sColor="cyan")
    addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","JAJUXE",sColor="green")
    addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","Formica",sColor="red")
    addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","Solenopsis_invicta_ref",sColor="black")
    
    vs=c("JAJUXC","JAJUXE","Formica","Solenopsis")
    vtc=c("cyan","green","red","black")
    #
    px=0.7
    xLegend=xlim[1]*(1-px)+xlim[2]*px
    py=1.03
    yLegend=ylim[1]*(1-py)+ylim[2]*py
    legend(xLegend,yLegend,#coordinates of left top corner 
           legend=vs,#captions
           box.lty=0,
           bg="transparent",
           text.col=vtc,#"red",
           #col=vColor, #colors
           #lty=vLineType,#types of line 
           #cex.lab=20,#cex=2,#fontSize,
           #box.lty=1,#(0 - no border of legend)
           #lwd=vwidthOfLines
    )
    
  }
  
  if(FALSE){
    #uzhe ne nado
    if(FALSE){
      #Cnig_gn3.1 vs JAJUXC,JAJUXE,Formica
      startMyPlotOxford(contigs1_Cnig_gn3.1,sx="Cnig_gn3.1")
      lwd=0.1
      addGridBasedOnContigs(contigs1_Cnig_gn3.1,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
      addGridBasedOnContigs(contigs5_Formica,sColor = "green",bHorizontal=TRUE,lenMin=2000000)
      addGridBasedOnContigs(contigs4_JAJUXC,sColor = "blue",bHorizontal=TRUE,lenMin=2000000)
      addGridBasedOnContigs(contigs3_JAJUXE,sColor = "red",bHorizontal=TRUE,lenMin=2000000)
      addGridBasedOnContigs(contigs6_SolenopsisInvicta,sColor = "darkgray",bHorizontal=TRUE,lenMin=2000000)
      #
      addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","Formica",sColor="green")
      addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","JAJUXC",sColor="blue")
      addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","JAJUXE",sColor="red")
      addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","Solenopsis_invicta_ref",sColor="darkgray")
      vs=c("Formica","JAJUXC","JAJUXE","Solenopsis")
      vtc=c("green","blue","red","darkgray")
      
      nCtg1=nrow(contigs1_Cnig_gn3.1)
      xmax=contigs1_Cnig_gn3.1$coorStart[nCtg1]+contigs1_Cnig_gn3.1$len[nCtg1]
      xLegend=xmax*0.7
      yLegend=xmax*0.6#1.45
      legend(xLegend,yLegend,#coordinates of left top corner 
             legend=vs,#captions
             box.lty=0,
             bg="transparent",
             text.col=vtc,#"red",
             #col=vColor, #colors
             #lty=vLineType,#types of line 
             #cex.lab=20,#cex=2,#fontSize,
             #box.lty=1,#(0 - no border of legend)
             #lwd=vwidthOfLines
      )
      
      
      
      
      #JAJUXE vs JAJUXC
      startMyPlotOxford(contigs4_JAJUXC,sx="JAJUXC",ycoef=1.1)
      addGridBasedOnContigs(contigs4_JAJUXC,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
      addGridBasedOnContigs(contigs3_JAJUXE,sColor = "darkgreen",bHorizontal=TRUE,lenMin=2000000)
      addPointsFromMarkers(dfMarkers,"JAJUXC","JAJUXE",sColor="black")
      
      
      #JAJUXE,"Formica" vs x= JAJUXC
      startMyPlotOxford(contigs4_JAJUXC,sx="JAJUXC",ycoef=1.1)
      addGridBasedOnContigs(contigs4_JAJUXC,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
      addGridBasedOnContigs(contigs3_JAJUXE,sColor = "darkgreen",bHorizontal=TRUE,lenMin=2000000)
      addGridBasedOnContigs(contigs5_Formica,sColor = "green",bHorizontal=TRUE,lenMin=2000000)
      addPointsFromMarkers(dfMarkers,"JAJUXC","JAJUXE",sColor="black")
      addPointsFromMarkers(dfMarkers,"JAJUXC","Formica",sColor="green")
      
      #JAJUXE,"Formica" vs x= JAJUXC
      startMyPlotOxford(contigs4_JAJUXC,sx="JAJUXC",ycoef=1.1)
      addGridBasedOnContigs(contigs4_JAJUXC,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
      addGridBasedOnContigs(contigs3_JAJUXE,sColor = "darkgreen",bHorizontal=TRUE,lenMin=2000000)
      addGridBasedOnContigs(contigs5_Formica,sColor = "green",bHorizontal=TRUE,lenMin=2000000)
      addGridBasedOnContigs(contigs1_Cnig_gn3.1,sColor = "red",bHorizontal=TRUE,lenMin=1000000)
      addPointsFromMarkers(dfMarkers,"JAJUXC","JAJUXE",sColor="black")
      addPointsFromMarkers(dfMarkers,"JAJUXC","Formica",sColor="green")
      addPointsFromMarkers(dfMarkers,"JAJUXC","Cnig_gn3.1",sColor="red")
      
      
      
      
      
    
    }
  }
  
  #b and bb vs ref
  if(FALSE){
    #b and bb vs ref
    #Solinopsis_b_and_Solinopsis_bb_and_Cnig_gn3.1_and_Formica_vs_Solinopsis_ref.pdf
    sFileMarkers6="C:/Frenkel/LTCPython/VovaPy/MarkersWithMultiplePos_Solenopsis_invicta.txt"
    dfMarkers6<-read.csv(sFileMarkers6, header = TRUE, sep = "\t", quote = "\"",
                         dec = ".", fill = TRUE, comment.char = "#")
    dfMarkers6<-addCoorColumnToMarkers(dfMarkers6,sMap="Solenopsis_invicta_ref",contigs6_SolenopsisInvicta) 
    dfMarkers6<-addCoorColumnToMarkers(dfMarkers6,sMap="Solenopsis_invicta_b",contigs6_b)
    dfMarkers6<-addCoorColumnToMarkers(dfMarkers6,sMap="Solenopsis_invicta_bb",contigs6_bb)
    
    ll=startMyPlotOxford(contigs6_SolenopsisInvicta,sx="Solenopsis_invicta_ref",ycoef=1.1)
    xlim=ll[[1]]
    ylim=ll[[2]]
    addGridBasedOnContigs(contigs6_SolenopsisInvicta,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
    addGridBasedOnContigs(contigs6_b,sColor = "cyan",bHorizontal=TRUE,lenMin=2000000)
    addGridBasedOnContigs(contigs6_bb,sColor = "green",bHorizontal=TRUE,lenMin=2000000)
    addGridBasedOnContigs(contigs1_Cnig_gn3.1,sColor = "red",bHorizontal=TRUE,lenMin=1000000)
    
    addPointsFromMarkers(dfMarkers6,"Solenopsis_invicta_ref","Solenopsis_invicta_b",sColor="cyan")
    addPointsFromMarkers(dfMarkers6,"Solenopsis_invicta_ref","Solenopsis_invicta_bb",sColor="green")
    addPointsFromMarkers(dfMarkers,"Solenopsis_invicta_ref","Cnig_gn3.1",sColor="red")
    
    vs=c("b","bb","Cnig_gn3.1")
    vtc=c("cyan","green","red")
    #
    px=0.7
    xLegend=xlim[1]*(1-px)+xlim[2]*px
    py=1.03
    yLegend=ylim[1]*(1-py)+ylim[2]*py
    legend(xLegend,yLegend,#coordinates of left top corner 
           legend=vs,#captions
           box.lty=0,
           bg="transparent",
           text.col=vtc,#"red",
           #col=vColor, #colors
           #lty=vLineType,#types of line 
           #cex.lab=20,#cex=2,#fontSize,
           #box.lty=1,#(0 - no border of legend)
           #lwd=vwidthOfLines
    )
  }
  
  #LD
  if(FALSE){
    #LD_plot_Cnig_gn3.1_RADmarkers_669diploids.pdf
    #LD_plot_Cnig_gn3.1_RADmarkers_669diploids_and_twoCathagliphisIspanica_and_Formica.pdf
    bNotOnlyLD=FALSE
    #bNotOnlyLD=TRUE
    
    
    sFileLD="C:/Frenkel/Privman/NanoporeSequencingSupergene/TrainingData/LD_allSignif.txt"
    dfMarkersLD<-read.csv(sFileLD, header = TRUE, sep = "\t", quote = "\"",
                          dec = ".", fill = TRUE, comment.char = "#")
    dfMarkersLD<-addCoorColumnToMarkers(dfMarkersLD,sMap="m1Cnig_gn3.1",contigs1_Cnig_gn3.1,bType=FALSE) 
    dfMarkersLD<-addCoorColumnToMarkers(dfMarkersLD,sMap="m2Cnig_gn3.1",contigs1_Cnig_gn3.1,bType=FALSE)
    #
    ll=startMyPlotOxford(contigs1_Cnig_gn3.1,sx="Cnig_gn3.1",xcoef=1,ycoef=1.2)
    xlim=ll[[1]]
    ylim=ll[[2]]
    
    lwd=0.01
    addGridBasedOnContigs(contigs1_Cnig_gn3.1,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
    addGridBasedOnContigs(contigs1_Cnig_gn3.1,sColor = "darkgreen",bHorizontal=TRUE,lenMin=1000000)
    if(bNotOnlyLD){
      addGridBasedOnContigs(contigs4_JAJUXC,sColor = "cyan",bHorizontal=TRUE,lenMin=2000000)
      addGridBasedOnContigs(contigs3_JAJUXE,sColor = "green",bHorizontal=TRUE,lenMin=2000000)
      addGridBasedOnContigs(contigs5_Formica,sColor = "black",bHorizontal=TRUE,lenMin=2000000)
    }
    #
    addPointsFromMarkers(dfMarkersLD,"m1Cnig_gn3.1","m2Cnig_gn3.1",sColor="red")
    if(bNotOnlyLD){
      addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","JAJUXC",sColor="cyan")
      addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","JAJUXE",sColor="green")
      addPointsFromMarkers(dfMarkers,"Cnig_gn3.1","Formica",sColor="black")
      
      vs=c("Cnig_gn3.1","JAJUXC","JAJUXE","Formica")
      vtc=c("red","cyan","green","black")
      #
      px=0.75
      xLegend=xlim[1]*(1-px)+xlim[2]*px
      py=1.03
      yLegend=ylim[1]*(1-py)+ylim[2]*py
      legend(xLegend,yLegend,#coordinates of left top corner 
             legend=vs,#captions
             box.lty=0,
             bg="transparent",
             text.col=vtc,#"red",
             #col=vColor, #colors
             #lty=vLineType,#types of line 
             #cex.lab=20,#cex=2,#fontSize,
             #box.lty=1,#(0 - no border of legend)
             #lwd=vwidthOfLines
      )
      
    }
  }
  
  #no need
  if(FALSE){
    #geneticMap testing (not Good)
    sFileGM="C:/Frenkel/Privman/NanoporeSequencingSupergene/TrainingData/rec.txt"
    dfMarkersGM<-read.csv(sFileGM, header = TRUE, sep = "\t", quote = "\"",
                          dec = ".", fill = TRUE, comment.char = "#")
    nm_mm<-nrow(dfMarkersGM)
    dfMarkersGM["x"]<-NA
    dfMarkersGM["y"]<-NA
    #vx=list()
    #vy=list()
    for (im in 1:nm_mm){
      #im=1
      #mName1	X2	p	rec	d_cM
      df=dfMarkers[dfMarkers$marker == dfMarkersGM$mName[im], ]
      x<-df$x_Cnig_gn3.1[1]
      df=dfMarkers[dfMarkers$marker == dfMarkersGM$mName1[im], ]
      y<-df$x_Cnig_gn3.1[1]
      dfMarkersGM[im,"x"]=x
      dfMarkersGM[im,"y"]=y
      print(im)
    }
    startMyPlotOxford(contigs1_Cnig_gn3.1,sx="Cnig_gn3.1",xcoef=0.2,ycoef=1.2)
    lwd=0.01
    addGridBasedOnContigs(contigs1_Cnig_gn3.1,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
    addGridBasedOnContigs(contigs1_Cnig_gn3.1,sColor = "darkgreen",bHorizontal=TRUE,lenMin=1000000)
    #vx<-dfMarkersGM[,"x"]
    #vy<-dfMarkersGM[,"y"]
    #
    #Create a function to generate a continuous color palette
    rbPal <- colorRampPalette(c('red','blue',"green","yellow"))
    #
    #This adds a column of color values
    # based on the y values
    dfMarkersGM$Col <- rbPal(100)[as.numeric(cut(dfMarkersGM$d_cM,breaks = 100))]
    #
    #plot(dat$x,dat$y,pch = 20,col = dat$Col)
    #
    #points(x=vx, y=vy, type = "p", pch=".",col=sColor)
    points(dfMarkersGM$x, dfMarkersGM$y, type = "p", pch=1,cex=0.01,col=dfMarkersGM$Col)
  }
  
  #geneticMap testing
  if(FALSE){
    sFileGM="C:/Frenkel/Privman/NanoporeSequencingSupergene/TrainingData/recPP2.txt"
    dfMarkersGMpp<-read.csv(sFileGM, header = TRUE, sep = "\t", quote = "\"",
                            dec = ".", fill = TRUE, comment.char = "#")
    rbPal <- colorRampPalette(c('red','blue',"green","yellow"))
    dfMarkersGMpp$Col <- rbPal(100)[as.numeric(cut(dfMarkersGMpp$d_cM,breaks = 100))]
    #dfMarkersGMpp$Col<-"black"
    #
    #plot(dat$x,dat$y,pch = 20,col = dat$Col)
    #
    #points(x=vx, y=vy, type = "p", pch=".",col=sColor)
    
    #chr2
    #startMyPlotOxford(contigs1_Cnig_gn3.1,sx="Cnig_gn3.1, chr02",xcoef=0.12,ycoef=1,xcoefMin=0.55,ycoefMin=0.55)
    
    #all
    startMyPlotOxford(contigs1_Cnig_gn3.1,sx="Cnig_gn3.1, all",xcoef=0.8,ycoef=1)
    
    #chr3
    #startMyPlotOxford(contigs1_Cnig_gn3.1,sx="Cnig_gn3.1, chr03",xcoef=0.172,ycoef=1,xcoefMin=0.7,ycoefMin=0.7)
    
    #chr9
    #startMyPlotOxford(contigs1_Cnig_gn3.1,sx="Cnig_gn3.1, chr09",xcoef=0.44,ycoef=1,xcoefMin=0.86,ycoefMin=0.86)
    
    
    
    lwd=0.01
    addGridBasedOnContigs(contigs1_Cnig_gn3.1,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
    addGridBasedOnContigs(contigs1_Cnig_gn3.1,sColor = "darkgreen",bHorizontal=TRUE,lenMin=1000000)
    points(dfMarkersGMpp$coorSS1, dfMarkersGMpp$coorSS2, type = "p", pch=1,cex=0.01,col=dfMarkersGMpp$Col)
    points(dfMarkersGMpp$coorSS2, dfMarkersGMpp$coorSS1, type = "p", pch=1,cex=0.01,col=dfMarkersGMpp$Col)
  }
  
  
}

{
  #plot dist control
  
  sFile="C:/Frenkel/Privman/Alexandra/GeneticMaps/camponotus_202205/camponotus_MP.txt.distControl.txt"
  myData<-read.csv(sFile, header = TRUE, sep = "\t", quote = "\"",
             dec = ".", fill = TRUE, comment.char = "#")
  plot(myData$dPhys,myData$dGenetic)#marker1	marker2	ctg	coor1	coor2	dPhys	r	dGenetic
}







addLineToPlot<-function(x0,y0,x1,y1,col="black",w=1,sLineType="solid"){
  #plot(0, 0, col = "white", xlab = "", ylab = "")
  # Draw one line
  # Draw one line as in Example 1
  segments(x0 = x0, y0 = y0, x1 =x1, y1 = y1,
           
           # Color of line
           col = col,#"darkgreen",      
           
           # Thickness of line
           lwd = w,#5,                
           
           # Line type
           lty = sLineType)#"dotted")
}
addTextToPlot<-function(x,y,sText,color="black",pTextSizeRelativelyToCurrentFontZize=0.65,iTypeFromTheObject=0){
  #https://stat.ethz.ch/R-manual/R-devel/library/graphics/html/text.html
  text(x, y,  sText,
       cex=0.65, pos=3,col="red")
  #pos:
  #1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the specified (x,y)
}
drewMaps<-function(){
  plot(
    #no need this point (need only to start graph)
    0,#vx 
    0, #vy
    col="white",#"yellow", #color of points
    #
    # log="x", #axis x in logarithmic scale 
    #xlab=sx , #caption of axis x
    #ylab=sy, #caption of axis y
    xaxt="n", #no caption foir axis x (will be added manually)
    yaxt="n", #no caption foir axis y (will be added manually)
    
    xlim=c(-10,100+10), #xMin,xMax
    ylim=c(-10,100+10), #xMin,xMax
    #cex=2,#scaling of font size
    #cex.lab=fontSize/10,#scaling of font size (lable of axes)
    #cex.axis=fontSize/10,#scaling of font size (numbers)
    #mar=c(5.1,8.1,4.1,2.1),#margines:bootom, left,top,right
  )
  addLineToPlot(10,0,10,30,col="blue",w=5)
  text(10, 0,  "ctg", cex=0.65, pos=1,col="black")
}


drewObjectsVovaFromTxtFile<-function(sFileName,bHeader = FALSE){
  bHeader = FALSE
  #sFileName="C:/Frenkel/LTCPython/vovaGraphs/graphDataExample.txt"
  #sFileName="C:/Frenkel/LTCPython/VovaPy/graphProba___1.txt"
  dfMyTable<-read.csv(sFileName, header = bHeader, sep = "\t", quote = "\"",
                      dec = ".", fill = TRUE, comment.char = "#")
  nRows<-nrow(dfMyTable)
  nCols<-ncol(dfMyTable)
  
  xlim<-c(0,10)
  ylim<-c(0,10)
  for (iRow in 1:nRows){
    #if (bType){
    #  i<-dfMarkers[im,sShapka_n]
    s<-dfMyTable[iRow,1]
    
    iCol<-2
    b<-iCol<=nCols
    if(b){
      b<-dfMyTable[iRow,iCol]!=""
    }
    n<-1
    while(b){
      #print(paste("im=",im, sep = ""))
      
      #print(dfMyTable[iRow,iCol])
      iCol<-iCol+1
      if(b){
        b<-iCol<=nCols
        if(b){b<-!is.na(dfMyTable[iRow,iCol])}
        if(b){
          b<-dfMyTable[iRow,iCol]!=""
        }
        if(b){
          n<-iCol
        }
      }
    }
    #print(n)
    
    if(s=="xlim"){
      dx<-as.numeric(dfMyTable[iRow,4])
      xlim[1]<-as.numeric(dfMyTable[iRow,2])-dx
      xlim[2]<-as.numeric(dfMyTable[iRow,3])+dx
    }
    if(s=="ylim"){
      dy<-as.numeric(dfMyTable[iRow,4])
      ylim[1]<-as.numeric(dfMyTable[iRow,2])-dy
      ylim[2]<-as.numeric(dfMyTable[iRow,3])+dy
      plot(
        #no need this point (need only to start graph)
        0,#vx 
        0, #vy
        col="white",#"yellow", #color of points
        xaxt="n", #no caption foir axis x (will be added manually)
        yaxt="n", #no caption foir axis y (will be added manually)
        
        xlim=xlim, #xMin,xMax
        ylim=ylim, #xMin,xMax
      )
    }
    if(s=="text"){
      #text	10	0	"ctg"	cex	0.65	pos	1	col	"black"
      x<-as.numeric(dfMyTable[iRow,2])
      y<-as.numeric(dfMyTable[iRow,3])
      ss<-dfMyTable[iRow,4]
      cex=0.65
      pos=1
      col="green"
      for (iCol in 5:n){
        if(dfMyTable[iRow,iCol]=="cex"){cex<-as.numeric(dfMyTable[iRow,iCol+1])}
        if(dfMyTable[iRow,iCol]=="pos"){pos<-as.numeric(dfMyTable[iRow,iCol+1])}
        if(dfMyTable[iRow,iCol]=="col"){col<-dfMyTable[iRow,iCol+1]}
      }
      text(x, y,  ss, cex=cex, pos=pos,col=col)
    }
    if(s=="segments"){
      #segments	10	0	10	30	col	"darkgreen"	lwd	5	lty	sLineType	#"dotted"
      x0<-as.numeric(dfMyTable[iRow,2])
      y0<-as.numeric(dfMyTable[iRow,3])
      x1<-as.numeric(dfMyTable[iRow,4])
      y1<-as.numeric(dfMyTable[iRow,5])
      lwd=1
      col="green"
      lty="solid"#"dotted"
      for (iCol in 6:n){
        if(dfMyTable[iRow,iCol]=="lwd"){lwd<-as.numeric(dfMyTable[iRow,iCol+1])}
        if(dfMyTable[iRow,iCol]=="lty"){lty<-dfMyTable[iRow,iCol+1]}
        if(dfMyTable[iRow,iCol]=="col"){col<-dfMyTable[iRow,iCol+1]}
        
        #h_hexNumber -> #hexNumber
        if(substr(col,1,2)=="h_"){col=paste("#", substr(col,3,nchar(col)), sep = "")}  #{col="#"+col}
      }
      print(col)
      segments(x0 = x0, y0 = y0, x1 =x1, y1 = y1,col=col,lwd=lwd,lty=lty)
    }
    if(s=="arrows"){
      #arrows	10	0	10	30	col	"darkgreen"	lwd	5	lty	sLineType	#"dotted"
      x0<-as.numeric(dfMyTable[iRow,2])
      y0<-as.numeric(dfMyTable[iRow,3])
      x1<-as.numeric(dfMyTable[iRow,4])
      y1<-as.numeric(dfMyTable[iRow,5])
      length<-0.25
      lwd=1
      col="green"
      lty="solid"#"dotted"
      code=2#0 -no, 1 - from 1 to 0, 2 -from 0 to 1, 3 - both
      for (iCol in 6:n){
        if(dfMyTable[iRow,iCol]=="length"){length<-as.numeric(dfMyTable[iRow,iCol+1])}
        if(dfMyTable[iRow,iCol]=="lwd"){lwd<-as.numeric(dfMyTable[iRow,iCol+1])}
        if(dfMyTable[iRow,iCol]=="lty"){lty<-dfMyTable[iRow,iCol+1]}
        if(dfMyTable[iRow,iCol]=="col"){col<-dfMyTable[iRow,iCol+1]}
        if(dfMyTable[iRow,iCol]=="code"){code<-as.numeric(dfMyTable[iRow,iCol+1])}
      }
      arrows(x0 = x0, y0 = y0, x1 =x1, y1 = y1, length=length,col=col,lwd=lwd,lty=lty,code=code)
    }
    #https://www.rdocumentation.org/packages/plotrix/versions/3.8-2/topics/draw.arc
    #draw.arc(x=1,y=NULL,radius=1,angle1=deg1*pi/180,angle2=deg2*pi/180,deg1=0,deg2=45,n=0.05,col=NA,lwd=NA,...)
    if(s=="draw.arc"){
      #draw.arc	50	50	20	30	90	col	"green"	lwd	"solid"
      x0<-as.numeric(dfMyTable[iRow,2])
      y0<-as.numeric(dfMyTable[iRow,3])
      radius<-as.numeric(dfMyTable[iRow,4])
      alphaGrad0<-as.numeric(dfMyTable[iRow,5])
      alphaGrad1<-as.numeric(dfMyTable[iRow,6])
      lwd=1
      col="green"
      lty="solid"#"dotted"
      for (iCol in 7:n){
        if(dfMyTable[iRow,iCol]=="lwd"){lwd<-as.numeric(dfMyTable[iRow,iCol+1])}
        if(dfMyTable[iRow,iCol]=="lty"){lty<-dfMyTable[iRow,iCol+1]}
        if(dfMyTable[iRow,iCol]=="col"){col<-dfMyTable[iRow,iCol+1]}
      }
      draw.arc(x=x0,y=y0,radius=radius,deg1=alphaGrad0,deg2=alphaGrad1,col=col,lwd=lwd,lty=lty)
    }
    #points(x=vx, y=vy, type = "p", pch=1,cex=0.01,col=sColor)
    if(s=="points"){
      x0<-as.numeric(dfMyTable[iRow,2])
      y0<-as.numeric(dfMyTable[iRow,3])
      pch=1
      cex=0.01
      col="green"
      type="p"
      for (iCol in 4:n){
        if(dfMyTable[iRow,iCol]=="pch"){pch<-as.numeric(dfMyTable[iRow,iCol+1])}
        if(dfMyTable[iRow,iCol]=="cex"){cex<-as.numeric(dfMyTable[iRow,iCol+1])}
        if(dfMyTable[iRow,iCol]=="col"){col<-dfMyTable[iRow,iCol+1]}
        if(dfMyTable[iRow,iCol]=="type"){type<-dfMyTable[iRow,iCol+1]}
      }
      points(x=x0,y=y0,type = type, pch=pch,cex=cex,col=col)
    }
  }
}

drewObjectsVovaFromTxtFile("C:/Frenkel/LTCPython/VovaPy/graphProba___1.txt")




arrows(0, 0, 30, 30, length = 0.25, angle = 30,code = 3)














b<-a

nCtg1=nrow(a)
nCtg2=nrow(b)
xmax=a$coorStart[nCtg1]+a$len[nCtg1]
ymax=a$coorStart[nCtg2]+a$len[nCtg2]
vx=c(0,xmax)
vy=c(ymax,0)
plot(vx,vy)
#grid on x
addGridBasedOnContigs(a,sColor = "darkgreen",bHorizontal=FALSE,lenMin=1000000)
addGridBasedOnContigs(b,sColor = "darkgreen",bHorizontal=TRUE,lenMin=1000000)


sFileBlast="C:/Frenkel/Privman/Cnig_gn1/res_Cnig_gn1_vs_FormicaSelysi_e70.out"
dfblast<-read.csv(sFileBlast, header = FALSE, sep = "\t", quote = "\"",
                  dec = ".", fill = TRUE, comment.char = "#")







#c<-dfMarkers["Cnig_gn3.1_pos"]






#dfMarkers$x_temp[1]=10000000
#dfMarkers$y_temp[1]=10000000
vx<-dfMarkers$x_temp
vy<-dfMarkers$y_temp
points(x=vx, y=vy, type = "p", pch=".",col="black")




#sFile="C:/Frenkel/Privman/GWAS/ind/LD_nOfStateMin_bestIndsFromColonies_25.txtLTCgraph.txt"
#myTable=read.table(file = sFile, sep = '\t', header = TRUE)
#myTable=read.table(file = sFile, sep = '\t', header = FALSE)
