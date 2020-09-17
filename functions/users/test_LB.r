
if(TRUE){#select subsets of lb data
  markerInCohort<-function(comDD, #the combined data gain in getData_dy.r
                           subjid, #subject id (multiple selection)
                           armcd,  #cohort information (multiple selection)
                           lbspec, #specific (single selection)
                           lbtest,  #lab test marker (single selection)
                           selColumn,#select the other columns (multiple selection)
                           sex     #gender
  ){
    subd<-comDD
    subd$SUBJID<-as.character(subd$SUBJID)
    #print('markerInCohort_begin');print(sort(unique(subd$SUBJID)));print('markerInCohort_end');
    subd<-subd[subd$SUBJID%in%subjid,]
    subd<-subd[subd$ARMCD%in%armcd,]
    subd<-subd[subd$LBSPEC==lbspec&subd$LBTEST==lbtest,]
    subd<-subd[subd$SEX%in%sex,]
    subd<-subd[,c(keepColumn,selColumn)]
    return(subd)
  }  
  
}

if(TRUE){#draw profile plot for the selected lb markers
  
  #----Univariate Y-axis plot----#
  #----Univariate Y-axis plot----#
  #----Univariate Y-axis plot----#
  profilePlot_1Y<-function(comDD, #the combined data gain in getData_dy.r
                           subjid, #subject id (multiple selection)
                           armcd,  #cohort information (multiple selection)
                           lbspec, #specific (single selection)
                           lbtest, #lab test marker (single selection)
                           sex,    #gender
                           Xaxis='numeric',
                           strip_cex=0.7, pt_cex=0.7, islog=FALSE, yrg=NA, xrg=NA,
                           showkey=TRUE,
                           dataonly=FALSE,
                           xlog=FALSE,
                           changeFromBL=FALSE,
                           changeFromBL.rf2="screening",  #or "day1", "both"
                           addBoxplot=TRUE,
                           onlyBoxplot=TRUE,
                           panelByPart=NULL,
                           returnData=FALSE,
                           addPval = TRUE,   #add pvalues
                           drawNumSamples=FALSE, #create a plot showing numbers of samples in Y-axi
                           refDay=1,         #the reference day as the baseline
                           exposure=NULL,       #add the exposure data
                           addExposure=TRUE,
                           visitAnno=NULL,    #a table of user-defined visit.anno
                           avgDup=FALSE       #using the averge value over duplicates                           
  ){
    
    subd<-comDD
    
    #clean up user-defined visit annotation
    if(!is.null(visitAnno)){
      #colnames: LBDY.L  LBDY.H  visit.anno
    }
    
    #clean up exposure data
    if(!is.null(exposure) && addExposure){
      myex <- unique(exposure[,c('USUBJID','EXDOSE', 'EXDOSU', 'EXTRT', 'EXSTDY')])
      myex <- reshape(data=myex, v.names = "EXDOSE", timevar='EXTRT',
                      idvar = c("USUBJID", 'EXSTDY'), direction='wide')
      myex$LBDY <- myex$EXSTDY
      subd <- merge(subd, myex, by=c('USUBJID','LBDY'), all=TRUE)
      library(latticeExtra)
    }
    
    library(lattice)
    library(latticeExtra)
    
    if(is.null(panelByPart)) panelByPart<- FALSE
    if(is.null(onlyBoxplot)) onlyBoxplot<- FALSE
    if(is.null(addBoxplot)) addBoxplot<- FALSE
    
    
    
    subd$LBTESTCD[is.na(subd$LBTESTCD)]<-''
    subd$SUBJID<-as.character(subd$SUBJID)
    subd<-subd[subd$SUBJID%in%subjid,]
    subd$arm<-subd$ARMCD
    subd<-subd[subd$ARMCD%in%armcd,]
    subd<-subd[subd$LBSPEC%in%lbspec&subd$LBTEST%in%lbtest,]
    subd<-subd[subd$SEX%in%sex,]
    
    #clean up duplicate values using LBSEQ and LBREFID
    #a problem identifed from JYCA data transfer
    #only keep the duplicat with the lowest LBSEQ
    if(all(c('LBSEQ','LBREFID') %in% colnames(subd))){
      if(!all(subd$LBREFID == "")){
        orig.nrow <- nrow(subd)
        smp.dup   <- unique(subd[,c('LBSEQ','LBREFID')])
        smp.seq   <- tapply(smp.dup$LBSEQ, smp.dup$LBREFID, 
                            function(x){sort(x)[1]})
        smp.dup   <- data.frame(LBSEQ=smp.seq, LBREFID=names(smp.seq))
        subd      <- merge(subd, smp.dup, by=c('LBSEQ','LBREFID'))
        print(paste('Remove',
                    orig.nrow-nrow(subd), 
                    'duplicate sample(s) for',
                    lbtest))
      }
    }
    
    subd<-subd[order(subd$LBDY),]
    un1<-unique(subd$LBSTRESU)[1]
    subd$SUBJID<-as.factor(subd$SUBJID)
    subd$LBSTRESN<-as.numeric(subd$LBSTRESN)
    subd<-subd[!is.na(subd$LBSTRESN)&subd$LBSTRESN!="",]
    subd$arm<-gsub('COHORT ','C', subd$arm, fixed=TRUE)
    if( tolower(substring(unique(subd$ARM), 1, 4))=='part' && panelByPart){
      subd$arm <- substring(subd$ARM, 1, 5)
      if(all(grepl(":", (subd$ARM)))){
        tma1 <- sapply(strsplit(dm$ARM, split=":", fixed=TRUE), function(x){x[1]})
        subd$arm <- paste(subd$arm, gsub(" ", "", tma1))
      }
    }else if (all(grepl(":", (subd$ARM))) && panelByPart){
      tma1 <- sapply(strsplit(subd$ARM, split=":", fixed=TRUE), function(x){x[1]})
      subd$arm <- gsub(" ", "", tma1)
    }
    
    if(drawNumSamples){
      subd$SUBJID <- subd$LBTEST
      subd$LBSTRESN <- 1
      un1 <- 'count'
    }
    
    if(changeFromBL && !drawNumSamples){
      if(changeFromBL.rf2=='screening'){
        s1<-subd$LBDY<1
        if(sum(s1)>0)
          tm1<-aggregate(LBSTRESN~SUBJID, data=subd[s1,c('SUBJID','LBSTRESN')], FUN=mean, na.rm=TRUE)
        else
          tm1<-data.frame(matrix( ncol=2, dimnames=list(NULL, c('SUBJID','LBSTRESN'))))
        s2<-subd$LBDY==1&!(subd$SUBJID%in%tm1$SUBJID)
        if(sum(s2)>0)
          tm2<-aggregate(LBSTRESN~SUBJID, data=subd[s2,c('SUBJID','LBSTRESN')], FUN=mean, na.rm=TRUE)
        else
          tm2<-data.frame(matrix( ncol=2, dimnames=list(NULL, c('SUBJID','LBSTRESN'))))
        if(nrow(tm1)>0 & nrow(tm2)==0) 
          tm12<-tm1
        else if (nrow(tm1)==0 & nrow(tm2)>0) 
          tm12<-tm2
        else
          tm12<-rbind(tm1, tm2)
      } else if(changeFromBL.rf2=="day1" ){
        if(!is.null(subd$LBBLFL) & sum(subd$LBBLFL=="Y")>0){ #use baseline flags
          tm12 <- unique(subd[subd$LBDY==refDay & subd$LBBLFL=='Y', c('SUBJID', 'LBSTRESN')])
        }else{
          tm12 <- unique(subd[subd$LBDY==refDay, c('SUBJID', 'LBSTRESN')])
        }
        try( tm12 <- aggregate(LBSTRESN~SUBJID, data=tm12, FUN=mean, na.rm=TRUE))
      } else {
        #or both
        tm12 <- unique(subd[subd$LBDY<=1, c('SUBJID', 'LBSTRESN')])
        tm12 <- aggregate(LBSTRESN~SUBJID, data=tm12, FUN=mean, na.rm=TRUE)
      }
      
      names(tm12)[names(tm12)=='LBSTRESN']<-'LBSTRESN.start'
      tm12$LBSTRESN.start <- tm12$LBSTRESN.start + .05
      subd<-merge(subd, tm12, by='SUBJID')
      #percent change from baseline
      subd$LBSTRESN<- 100*(subd$LBSTRESN-subd$LBSTRESN.start)/subd$LBSTRESN.start
      subd <- subd[order(subd[,c('LBDY')]),]
      subd$LBTESTCD[subd$LBTESTCD==""]<-subd$LBTEST[subd$LBTESTCD==""]
      subd$LBTESTCD <- paste("%chgFromBL for", subd$LBTESTCD)
      
      if(changeFromBL.rf2=="day1" ){
        subd <- subd[subd$LBDY> -1,]
      }
      
      if(avgDup){
        cc.na <- NULL
        for(cc in 1:ncol(subd)) {
          if(nrow(subd)>nrow(unique(subd[,-cc]))) cc.na <- c(cc.na, colnames(subd[cc]))
        }
        subd <- subd[,~colnames(subd)%in%c('LBBLFL')]
        try( subd <- aggregate(LBSTRESN~., data=subd, FUN=mean, na.rm=TRUE) )
      }
      
    }
    
    if(length(xrg)==0||any(is.na(xrg)|toupper(xrg)=='NA'|tolower(xrg)=='~')){
      xrg<-range(subd$LBDY,na.rm=TRUE,finite=TRUE)
    }else{
      xrg<-paste(xrg, collapse=',')
      xrg<-gsub("(", "", xrg, fixed=TRUE)
      xrg<-gsub(")", "", xrg, fixed=TRUE)
      xrg<-as.numeric(strsplit(xrg,split=',', fixed=TRUE)[[1]])
      xrg<-xrg[!is.na(xrg)]
      if(length(xrg)<2){
        xrg<-range(subd$LBDY,na.rm=TRUE,finite=TRUE)
      }else{
        subd <- subd[subd$LBDY>=xrg[1]&subd$LBDY<=xrg[2],]
      }
    }
    
    lattice.theme <- trellis.par.get()
    #col <- lattice.theme$superpose.symbol$col[1:nS]
    if(islog && !drawNumSamples){
      subd$LBSTRESN<-log(subd$LBSTRESN)
      yL<-paste0("log(",subd$LBTESTCD[1],"): ", subd$LBTEST[1])
      if(length(yrg)==0||any(is.na(yrg)|toupper(yrg)=='NA'|tolower(yrg)=='~')){
        yrg<-range(subd$LBSTRESN,na.rm=TRUE,finite=TRUE)
      }else{
        yrg<-paste(yrg, collapse=',')
        yrg<-gsub("(", "", yrg, fixed=TRUE)
        yrg<-gsub(")", "", yrg, fixed=TRUE)
        yrg<-as.numeric(strsplit(yrg,split=',', fixed=TRUE)[[1]])
        yrg<-yrg[!is.na(yrg)]
        if(length(yrg)<2){yrg<-range(subd$LBSTRESN,na.rm=TRUE,finite=TRUE)}
      }
      if(length(xrg)==0||any(is.na(xrg)|toupper(xrg)=='NA'|tolower(xrg)=='~')){
        xrg<-range(subd$LBDY,na.rm=TRUE,finite=TRUE)
      }else{
        xrg<-paste(xrg, collapse=',')
        xrg<-gsub("(", "", xrg, fixed=TRUE)
        xrg<-gsub(")", "", xrg, fixed=TRUE)
        xrg<-as.numeric(strsplit(xrg,split=',', fixed=TRUE)[[1]])
        xrg<-xrg[!is.na(xrg)]
        if(length(xrg)<2){xrg<-range(subd$LBDY,na.rm=TRUE,finite=TRUE)}
      }
      if(showkey){
        p1<-xyplot(LBSTRESN~LBDY|arm, data=subd, groups=SUBJID,
                   ylab=yL, xlab="LBDY", type=c('o','g'), cex=pt_cex, xlim=xrg,ylim=yrg,
                   scales=list(x=list(cex=1, rot=90, log=xlog), y=list(cex=1)),
                   strip=strip.custom(par.strip.text=list(cex=strip_cex)),
                   par.settings = list(superpose.symbol = list(pch = 0:20, alhpa=0.75),
                                       superpose.line=list(lty=1:8)),
                   auto.key = list(space = "right", lines=TRUE, cex=strip_cex)  );
      }else{
        p1<-xyplot(LBSTRESN~LBDY|arm, data=subd, groups=SUBJID,
                   ylab=yL, xlab="LBDY", type=c('o','g'), cex=pt_cex, xlim=xrg,ylim=yrg,
                   scales=list(x=list(cex=1, rot=90, log=xlog), y=list(cex=1)),
                   strip=strip.custom(par.strip.text=list(cex=strip_cex)));
      }
    }else{
      if(drawNumSamples){
        yL <- 'Samples'
      }else{
        yL<-paste0(subd$LBTESTCD[1],": ", subd$LBTEST[1])
      }
      if(length(yrg)==0||any(is.na(yrg)|toupper(yrg)=='NA'|tolower(yrg)=='~')){
        yrg<-range(subd$LBSTRESN,na.rm=TRUE,finite=TRUE)
      }else{
        yrg<-paste(yrg, collapse=',')
        yrg<-gsub("(", "", yrg, fixed=TRUE)
        yrg<-gsub(")", "", yrg, fixed=TRUE)
        yrg<-as.numeric(strsplit(yrg,split=',', fixed=TRUE)[[1]])
        yrg<-yrg[!is.na(yrg)]
        if(length(yrg)<2){yrg<-range(subd$LBSTRESN,na.rm=TRUE,finite=TRUE)}
      }
      
      subd$LBDY2 <- subd$LBDY
      subd$LBDY2[subd$LBDY<0] <- -14
      
      #add test P-value
      if(addPval){
        subd$LBDY3 <- subd$LBDY-mean(subd$LBDY, na.rm=T)
        if(length(unique(subd$arm))>1){
          mod1 <- aov(LBSTRESN~LBDY3 * arm, data=subd)
          mod2 <- aov(LBSTRESN~LBDY3 + arm, data=subd)
          rs1  <- anova(mod1, mod2) 
          pv1  <- rs1[2, "Pr(>F)"]
          xL   <- paste0("Days after 1st dose. F_Pvalue(day*arm)=", round(pv1,4))
        }else{
          require(lme4)
          mod1 <- lme4::lmer(LBSTRESN~LBDY3+(1|SUBJID), data=subd)
          mod2 <- lme4::lmer(LBSTRESN~(1|SUBJID), data=subd)
          rs1  <- anova(mod1, mod2, refit=FALSE) #use REML
          #rs1  <- anova(mod1, mod2, refit=TRUE) #use ML
          pv1  <- rs1[2, "Pr(>Chisq)"]
          xL   <- paste0("Days after 1st dose. Chisq_Pvalue(day)=", round(pv1,4))
        }
        
      } else {xL   <- paste0("Days after 1st dose")}
      
      
      if(Xaxis!='numeric'){
        subd<-subd[order(subd$LBDY),]
        if(!is.null(visitAnno)){
          subd$LBDY2 <- ''
          my.anno.L <- NULL
          for(va.i in 1:nrow(visitAnno)){
            if(is.na(visitAnno$LBDY.L[va.i] )) next
            my.anno <- as.character(visitAnno$visit.anno[va.i])
            my.lw   <- visitAnno$LBDY.L[va.i]
            my.hg   <- visitAnno$LBDY.H[va.i]
            #enfore the anno to be screening if VISIT defined
            if(grepl('scr', tolower(my.anno))&'VISIT'%in%colnames(subd)){
              subd$LBDY2[subd$VISIT=='SCREENING']<-my.anno
              subd$LBDY[subd$VISIT=='SCREENING'] <- my.lw
            }
            if('VISIT'%in%colnames(subd)){
              subd$LBDY2[subd$LBDY>=my.lw & subd$LBDY<=my.hg & 
                           subd$VISIT!='SCREENING']<-my.anno
              subd$LBDY[subd$LBDY>=my.lw & subd$LBDY<=my.hg & 
                          subd$VISIT!='SCREENING']<-my.lw
            }else{
              subd$LBDY2[subd$LBDY>=visitAnno$LBDY.L[va.i] & 
                           subd$LBDY<=visitAnno$LBDY.H[va.i]]<-my.anno
              subd$LBDY[subd$LBDY>=visitAnno$LBDY.L[va.i] & 
                          subd$LBDY<=visitAnno$LBDY.H[va.i]]<-my.lw
            }#..
            my.anno.L <- c(my.anno.L, my.anno)
          }
          subd <- subd[subd$LBDY2!='', ]
          subd<-subd[subd$LBDY>=xrg[1]&subd$LBDY<=xrg[2],]
          my.anno.L <- my.anno.L[my.anno.L%in%unique(subd$LBDY2)]
          subd$LBDY2 <- factor(subd$LBDY2, level=my.anno.L)
        }else if('VISIT'%in%colnames(subd)){
          sel.day1 <- !is.na(subd$LBDY) && subd$LBDY==1
          sel.dy <- subd$LBDY<0 | subd$VISIT=='SCREENING'
          subd$LBDY2[!sel.dy]<-paste0("D",subd$LBDY2[!sel.dy])
          #subd$LBDY2[!sel.dy]<-paste0("W",ceiling(subd$LBDY2[!sel.dy]/7))
          visit <- gsub('SCREENING', 'SCRN', toupper(subd$VISIT))
          visit[grepl('follow', tolower(visit))] <- 'FU'
          visit <- gsub('CYCLE', 'C', visit)
          #visit[sel.day1] <- 'C1'
          sel.997 <- grepl('997', visit)
          visit[sel.997]<-paste0("C",ceiling(subd$LBDY[sel.997]/28))
          #visit <- gsub('997', '1', visit)
          head( lbd <- strptime(subd$LBDTC,format="%Y-%m-%dT%H:%M") )
          head( std <- strptime(subd$RFSTDTC,format="%Y-%m-%dT%H:%M") )
          prpo <- ifelse(std<=lbd, "", ""); #'pre')
          visit[sel.day1] <- paste0(visit[sel.day1], prpo[sel.day1])
          subd$LBDY2<-paste0(subd$LBDY2, visit)
          #if(addBoxplot){
          subd$LBDY2[sel.dy]<-'SCRN'
          #} else {
          #  subd$LBDY2[sel.dy]<-paste0("D", subd$LBDY[sel.dy], 'SCRN')
          #}
        }
        
        #order the X-axis variable
        subd <- subd[order(subd$LBDY),]
        print( ord1 <- unique(subd[,c("LBDY2", "LBDY")]) )
        subd$LBDY2 <- factor(subd$LBDY2, level=ord1$LBDY2[order(ord1$LBDY)])
        subd$LBDY  <-subd$LBDY2
        
        if(drawNumSamples){
          subd<-aggregate(LBSTRESN~arm+ARMCD+LBDY+SUBJID, data=subd, FUN=sum, na.rm=T)
        }
        
        if(dataonly){
          return(list(subd=subd, yL=yL, pt_cex=pt_cex,
                      yrg=yrg, xrg=xrg, strip_cex=strip_cex))
        }
        if(onlyBoxplot){
          subd2<-subd
          subd <- aggregate(LBSTRESN~arm+ARMCD+LBDY, data=subd, FUN = mean, na.rm=T)
          subd$SUBJID <- subd$ARMCD
        }
        if(showkey){
          p1<-xyplot(LBSTRESN~factor(LBDY)|arm, data=subd, groups=SUBJID,
                     ylab=yL, xlab=xL, type=c('o','g'), cex=pt_cex, ylim=yrg, 
                     scales=list(x=list(cex=0.65, rot=30, log=xlog), y=list(cex=1)),
                     strip=strip.custom(par.strip.text=list(cex=strip_cex)),
                     par.settings = list(superpose.symbol = list(pch = 0:20, alhpa=0.75),
                                         superpose.line=list(lty=1:8)),
                     auto.key = list(space = "right", lines=TRUE, cex=strip_cex)  );
        }else{
          p1<-xyplot(LBSTRESN~factor(LBDY)|arm, data=subd, groups=SUBJID,
                     ylab=yL, xlab=xL, type=c('o','g'), cex=pt_cex, ylim=yrg, 
                     scales=list(x=list(cex=0.65, rot=30, log=xlog), y=list(cex=1)),
                     strip=strip.custom(par.strip.text=list(cex=strip_cex)));
        }
        
        if(addBoxplot){
          if(onlyBoxplot){
            p1 <- p1+bwplot(LBSTRESN~ factor(LBDY)|arm, data=subd2, cex=0.9, pch='*', box.ratio=0.25,
                            fill='transparent', 
                            lattice.options=lattice.options(key=list(border='black')) )
          }else{
            p1 <- p1+ bwplot(LBSTRESN~ factor(LBDY)|arm, data=subd, cex=0.9, pch='*', box.ratio=0.25,
                             fill='transparent', 
                             lattice.options=lattice.options(key=list(border='black')) )
          }
        }
      }else{
        if(dataonly){
          return(list(subd=subd, yL=yL, pt_cex=pt_cex,
                      yrg=yrg, xrg=xrg, strip_cex=strip_cex))}
        
        if(drawNumSamples){
          subd<-aggregate(LBSTRESN~arm+ARMCD+LBDY+SUBJID, data=subd, FUN=sum, na.rm=T)
        }
        
        if(showkey){
          p1<-xyplot(LBSTRESN~LBDY|arm, data=subd, groups=SUBJID,
                     ylab=yL, xlab=xL, type=c('o','g'), cex=pt_cex, ylim=yrg, xlim=xrg,
                     scales=list(x=list(cex=1, rot=90, log=xlog), y=list(cex=1)),
                     strip=strip.custom(par.strip.text=list(cex=strip_cex)),
                     par.settings = list(superpose.symbol = list(pch = 0:20, alhpa=0.75),
                                         superpose.line=list(lty=1:8)),
                     auto.key = list(space = "right", lines=TRUE, cex=strip_cex)  );
        } else {
          p1<-xyplot(LBSTRESN~LBDY|arm, data=subd, groups=SUBJID,
                     ylab=yL, xlab=xL, type=c('o','g'), cex=pt_cex, ylim=yrg, xlim=xrg,
                     scales=list(x=list(cex=1, rot=90, log=xlog), y=list(cex=1)),
                     strip=strip.custom(par.strip.text=list(cex=strip_cex)));
        }
      }
    }
    
    if(!is.null(exposure) && addExposure && Xaxis=='numeric'){
      subd <- subd[!is.na(subd$EXSTDY), colSums(!is.na(subd))>1]
      trt.var <- colnames(subd)[grepl('EXDOSE.',colnames(subd), fixed=TRUE)]
      col.ref <- c('darkblue','red', 'darkgreen',  'orange')
      subd$col.ref1 <- 'white'
      subd$trtn     <- 0
      for(i in 1:length(trt.var)){
        subd$col.ref1[!is.na(subd[,trt.var[i]])] <- col.ref[i]
        subd$trtn[!is.na(subd[,trt.var[i]])] <- subd$trtn[!is.na(subd[,trt.var[i]])] +1
      }
      subd$col.ref1[subd$trtn>1] <- col.ref[i+1]
      
      p2 <- xyplot(LBSTRESN~EXSTDY|arm, 
                   data=subd,
                   #groups=EXTRT, 
                   fill.color = as.character(subd$col.ref1),
                   panel=function(x,y, fill.color,...,subscripts){
                     fill = fill.color[subscripts] 
                     panel.abline(v=x, col=fill, lwd=0.3)
                   })
      p1 <- p1+as.layer(p2)
      
      
    }
    
    if(returnData){
      return(list(p1=p1, dat1=subd))
    }else{
      return(p1)
    }
  }
  
  
  
  #----Summary Univariate Y-axis with error bar plot----#
  profilePlot_1Yer<-function(comDD, #the combined data gain in getData_dy.r
                             subjid, #subject id (multiple selection)
                             armcd,  #cohort information (multiple selection)
                             lbspec, #specific (single selection)
                             lbtest,  #lab test marker (single selection)
                             strip_cex=0.8, pt_cex=0.7,
                             Xaxis="numeric", #only "numeric" or "factor"
                             islog=FALSE, yrg=NA, xrg=NA
  ){
    subd<-comDD[comDD$SUBJID%in%subjid,]
    subd<-subd[subd$ARMCD%in%armcd,]
    subd<-subd[subd$LBSPEC==lbspec&subd$LBTEST==lbtest,]
    subd<-subd[order(subd$LBDY),]
    un1<-unique(subd$LBSTRESU)[1]
    subd$SUBJID<-as.factor(subd$SUBJID)
    subd$LBSTRESN<-as.numeric(subd$LBSTRESN)
    subd<-subd[!is.na(subd$LBSTRESN)&subd$LBSTRESN!="",]
    subd<-subd[subd$LBDY>- 21,]
    subd$week<-floor(subd$LBDY/7)+1
    #subd$week<-floor(subd$LBDY)
    if(islog){
      subd$LBSTRESN<-log(subd$LBSTRESN)
      yL<-paste0("log(",gsub(" ", '', subd$LBTEST[1]), ") ", 
                 " +/- 2*sdErr")
    }else{
      yL<-paste(gsub(" ", "", subd$LBTEST[1]), " +/- 2*sdErr")
    }
    
    #get averages
    subd<-unique(subd[,c("SUBJID","ARM","week","LBSTRESN","LBSTRESU")])
    tm1<-aggregate(x=subd[,"LBSTRESN"], by=list(subd$ARM, subd$week),
                   FUN=mean, na.rm=TRUE)
    colnames(tm1)<-c("ARM", "week", "avg_LBSTRESN")
    tm2<-aggregate(x=subd[,"LBSTRESN"], by=list(subd$ARM, subd$week),
                   FUN=sd, na.rm=TRUE)
    colnames(tm2)<-c("ARM", "week", "sd_LBSTRESN")
    tm3<-aggregate(x=subd[,"LBSTRESN"], by=list(subd$ARM, subd$week),
                   FUN=function(x){sum(!is.na(x))})
    colnames(tm3)<-c("ARM", "week", "freq")
    tm4<-aggregate(x=subd[,"SUBJID"], by=list(subd$ARM),
                   FUN=function(x){length(unique(x))})
    colnames(tm4)<-c("ARM", "ARM2")
    tm4$ARM2<-paste0(tm4$ARM,", N=", tm4$ARM2)
    
    tm12<-merge(tm1,tm2, by=c("ARM", "week"))
    tm12<-merge(tm12,tm3, by=c("ARM", "week"))
    tm12<-merge(tm12,tm4, by=c("ARM"))
    
    tm12<-tm12[!is.na(tm12$avg_LBSTRESN)&!is.na(tm12$sd_LBSTRESN),]
    tm12$low<-pmax(0, tm12$avg_LBSTRESN-2*tm12$sd_LBSTRESN)
    tm12$up<- tm12$avg_LBSTRESN + 2*tm12$sd_LBSTRESN
    tm12<-tm12[order(tm12$week),]
    
    y1rg<-range(c(tm12$low, tm12$up))
    if(TRUE){#change xlim and ylim
      if(length(yrg)==0||(is.na(yrg)|toupper(yrg)=='NA'|tolower(yrg)=='~')){
        yrg<-y1rg
      }else{
        yrg<-gsub("(", "", yrg, fixed=TRUE)
        yrg<-gsub(")", "", yrg, fixed=TRUE)
        yrg<-as.numeric(strsplit(yrg,split=',',fixed=TRUE)[[1]])
        if(any(is.na(yrg))){ yrg<-y1rg }
      }
      if(length(xrg)==0||(is.na(xrg)|toupper(xrg)=='NA'|tolower(xrg)=='~')){
        xrg<-range(tm12$week,na.rm=TRUE,finite=TRUE)
      }else{
        xrg<-gsub("(", "", xrg, fixed=TRUE)
        xrg<-gsub(")", "", xrg, fixed=TRUE)
        xrg<-as.numeric(strsplit(xrg,split=',',fixed=TRUE)[[1]])
        if(any(is.na(xrg))){ xrg<-range(tm12$week,na.rm=TRUE,finite=TRUE)}
      }
    }
    if(Xaxis=="factor"){
      p1<-xyplot(avg_LBSTRESN~factor(week)|ARM2, data=tm12, type=c('g',"b"), col="blue4",
                 ylab=yL, xlab="week",
                 scales=list(x=list(cex=1, rot=90),y=list(cex=1)),
                 ylim=yrg, groups=ARM,  cex=pt_cex,
                 strip=strip.custom(par.strip.text=list(cex=strip_cex)),
                 panel.groups=function(x, y, subscripts, col, ...){
                   low<-tm12$low[subscripts]
                   up<- tm12$up[subscripts]
                   label<-as.character(tm12$freq[subscripts])
                   yLabel<-rep(y1rg[2]*.95, length(label))
                   panel.xyplot(x, y, col=col, ...)
                   panel.arrows(x, low, x, up, angle=90, code=3, length=0.05, col=col) 
                   panel.text(x, yLabel, label=label, col=col, srt=90, cex=strip_cex)
                 }, panel='panel.superpose')
    }else{
      p1<-xyplot(avg_LBSTRESN~week|ARM2, data=tm12, type=c('g',"b"), col="blue4",
                 ylab=yL, xlab="week",
                 scales=list(x=list(cex=1, rot=90),y=list(cex=1)),
                 groups=ARM,  cex=pt_cex, xlim=xrg, ylim=yrg, 
                 strip=strip.custom(par.strip.text=list(cex=strip_cex)),
                 panel.groups=function(x, y, subscripts, col, ...){
                   low<-tm12$low[subscripts]
                   up<- tm12$up[subscripts]
                   label<-as.character(tm12$freq[subscripts])
                   yLabel<-rep(y1rg[2]*.95, length(label))
                   panel.xyplot(x, y, col=col, ...)
                   panel.arrows(x, low, x, up, angle=90, code=3, length=0.05, col=col) 
                   panel.text(x, yLabel, label=label, col=col, srt=90, cex=strip_cex)
                 }, panel='panel.superpose')
    } 
    return(p1)
  }
  
}

if(TRUE){#pairwise profile plot for two markers
  pkpd<-function(comDD, #the combined data gain in getData_dy.r
                 subjid, #subject id (multiple selection)
                 armcd,  #cohort information (multiple selection)
                 lbspec, #specific (single selection)
                 lbtest1,  #lab test marker (single selection)
                 lbtest2,  #lab test marker (single selection)
                 byVar='SUBJID'
  ){
    subd<-comDD[comDD$SUBJID%in%subjid,]
    subd<-subd[subd$ARMCD%in%armcd,]
    subd$LBTESTCD[is.na(subd$LBTESTCD)] <- ''
    subd<-subd[subd$LBSPEC==lbspec&subd$LBTEST==lbtest,]
    #get common subject
    s1<-unique(subd$SUBJID)
    s2<-unique(pk$SUBJID)
    s12<-s1[s1%in%s2]
    subd<-subd[subd$SUBJID%in%s12,]
    pk1<-pk[pk$SUBJID%in%s12,]
    
    if(!is.null(subd$LBELTM)){
      subd$LBELTM<-gsub("P", "", subd$LBELTM)
      subd$LBELTM<-gsub("T", "", subd$LBELTM)
      subd$LBELTM<-as.numeric(subd$LBELTM)/24
      subd$LBELTM[is.na(subd$LBELTM)]<-0
      subd$LBDY<-subd$LBDY+subd$LBELTM
      subd$LBELTM<-NULL
    }
    
    subd<-subd[order(subd$LBDY),]
    un1<-unique(subd$LBSTRESU)[1]
    subd$SUBJID<-as.factor(subd$SUBJID)
    subd$LBSTRESN<-as.numeric(subd$LBSTRESN)
    yL<-paste0(subd$LBTESTCD[1],": ", subd$LBTEST[1])
    p1<-xyplot(LBSTRESN~LBDY|subd[,byVar], data=subd, groups=SUBJID,
               ylab=yL, xlab="DY", type=c('o','g'), col="red",
               scales=list(x=list(cex=0.5, rot=90), y=list(cex=1)))
    
    if(!is.null(pk1$PCELTM)){
      pk1$PCELTM<-gsub("P", "", pk1$PCELTM)
      pk1$PCELTM<-gsub("T", "", pk1$PCELTM)
      pk1$PCELTM<-as.numeric(pk1$PCELTM)/24
    }
    pk1$PCELTM[is.na(pk1$PCELTM)]<-0
    pk1$PCDY<-pk1$PCDY+pk1$PCELTM
    pk1<-pk1[order(pk1$PCDY),]
    un2<-unique(pk1$PCSTRESU)[1]
    pk1$SUBJID<-as.factor(pk1$SUBJID)
    pk1$PCSTRESN<-as.numeric(pk1$PCSTRESN)
    yL2<-paste0(pk1$PCTESTCD[1],": ", pk1$PCTEST[1])
    p2<-xyplot(PCSTRESN~PCDY|pk1[,byVar], data=pk1, groups=SUBJID,
               ylab=yL2, xlab="Day of Study",type=c('o','g'), col="black", pch=16, cex=0.6,
               scales=list(x=list(cex=0.5, rot=90), y=list(cex=1)))
    
    p12<-doubleYScale(p2,p1, style1=0, style2=4,add.ylab2=TRUE, 
                      main=paste(yL2,"and",yL))
    
    return(p12)
  }
  
}

if(TRUE){#QC plot (histogram of selected lb markers)
  histQC<-function(comDD, #the combined data gain in getData_dy.r
                   subjid, #subject id (multiple selection)
                   armcd,  #cohort information (multiple selection)
                   lbspec, #specific (single selection)
                   lbtest  #lab test marker (single selection)
  ){
    subd<-comDD[comDD$SUBJID%in%subjid,]
    subd<-subd[subd$ARMCD%in%armcd,]
    subd<-subd[subd$LBSPEC==lbspec&subd$LBTEST==lbtest,]
    subd<-subd[order(subd$LBDY),]
    un1<-unique(subd$LBSTRESU)[1]
    subd$SUBJID<-as.factor(subd$SUBJID)
    subd$LBSTRESN<-as.numeric(subd$LBSTRESN)
    yL<-paste0(subd$LBTESTCD[1],": ", subd$LBTEST[1], " (", un1, " )")
    p1<-histogram(~LBSTRESN, data =subd, main=yL, col='red', type = "density",
                  panel=function(x, ...) {
                    panel.histogram(x, ...)
                    panel.grid()
                    panel.mathdensity(dmath = dnorm, col = "black", lwd=2,
                                      args = list(mean=mean(x),sd=sd(x)))
                  })
    return(p1)
  }
  
  histQC_arm<-function(comDD, #the combined data gain in getData_dy.r
                       subjid, #subject id (multiple selection)
                       armcd,  #cohort information (multiple selection)
                       lbspec, #specific (single selection)
                       lbtest  #lab test marker (single selection)
  ){
    subd<-comDD[comDD$SUBJID%in%subjid,]
    subd<-subd[subd$ARMCD%in%armcd,]
    subd<-subd[subd$LBSPEC==lbspec&subd$LBTEST==lbtest,]
    subd<-subd[order(subd$LBDY),]
    un1<-unique(subd$LBSTRESU)[1]
    subd$SUBJID<-as.factor(subd$SUBJID)
    subd$LBSTRESN<-as.numeric(subd$LBSTRESN)
    yL<-paste0(subd$LBTESTCD[1],": ", subd$LBTEST[1], " (", un1, " )")
    p1<-histogram(~LBSTRESN|factor(ARMCD), data =subd, main=yL, type = "density",
                  panel=function(x, ...) {
                    panel.histogram(x, ...)
                    panel.grid()
                    panel.mathdensity(dmath = dnorm, col = "red", lwd=2,
                                      args = list(mean=mean(x),sd=sd(x)))
                  })
    return(p1)
  }
}

if(TRUE){#Marker distribution with boxplots
  boxp_lb<-function(comDD, lbspec, lbtest, byVar){
    subd<-comDD[comDD$LBSPEC==lbspec&comDD$LBTEST==lbtest,]
    un1<-unique(subd$LBSTRESU)[1]
    subd$LBSTRESN<-as.numeric(subd$LBSTRESN)
    yL<-paste0(subd$LBTESTCD[1],": ", subd$LBTEST[1])   
    if(byVar%in%colnames(subd)){
      p1<-qplot(factor(subd[,byVar]),LBSTRESN,data=subd,geom="boxplot", 
                main=yL, xlab=byVar, ylab=yL)+
        theme(axis.text.x=element_text(angle=90))
    }else{ #draw overall
      subd$overall<-'overall'
      p1<-qplot(overall,LBSTRESN,data=subd,geom="boxplot", 
                main=yL, xlab=byVar, ylab=yL)+
        theme(axis.text.x=element_text(angle=90))
    }
    return(p1)
  }
}

if(TRUE){#PK/PD
  pkpd<-function(comDD, #the combined data gain in getData_dy.r
                 subjid, #subject id (multiple selection)
                 armcd,  #cohort information (multiple selection)
                 lbspec, #specific (single selection)
                 lbtest,  #lab test marker (single selection)
                 pk, byVar='SUBJID'
  ){
    subd<-comDD[comDD$SUBJID%in%subjid,]
    subd<-subd[subd$ARMCD%in%armcd,]
    subd<-subd[subd$LBSPEC==lbspec&subd$LBTEST==lbtest,]
    #get common subject
    s1<-unique(subd$SUBJID)
    s2<-unique(pk$SUBJID)
    s12<-s1[s1%in%s2]
    subd<-subd[subd$SUBJID%in%s12,]
    pk1<-pk[pk$SUBJID%in%s12,]
    
    if(!is.null(subd$LBELTM)){
      subd$LBELTM<-gsub("P", "", subd$LBELTM)
      subd$LBELTM<-gsub("T", "", subd$LBELTM)
      subd$LBELTM<-as.numeric(subd$LBELTM)/24
      subd$LBELTM[is.na(subd$LBELTM)]<-0
      subd$LBDY<-subd$LBDY+subd$LBELTM
      subd$LBELTM<-NULL
    }
    
    subd<-subd[order(subd$LBDY),]
    un1<-unique(subd$LBSTRESU)[1]
    subd$SUBJID<-as.factor(subd$SUBJID)
    subd$LBSTRESN<-as.numeric(subd$LBSTRESN)
    yL<-paste0(subd$LBTESTCD[1],": ", subd$LBTEST[1])
    p1<-xyplot(LBSTRESN~LBDY|subd[,byVar], data=subd, groups=SUBJID,
               ylab=yL, xlab="DY", type=c('o','g'), col="red",
               scales=list(x=list(cex=0.5, rot=90), y=list(cex=1)))
    
    if(!is.null(pk1$PCELTM)){
      pk1$PCELTM<-gsub("P", "", pk1$PCELTM)
      pk1$PCELTM<-gsub("T", "", pk1$PCELTM)
      pk1$PCELTM<-as.numeric(pk1$PCELTM)/24
    }
    pk1$PCELTM[is.na(pk1$PCELTM)]<-0
    pk1$PCDY<-pk1$PCDY+pk1$PCELTM
    pk1<-pk1[order(pk1$PCDY),]
    un2<-unique(pk1$PCSTRESU)[1]
    pk1$SUBJID<-as.factor(pk1$SUBJID)
    pk1$PCSTRESN<-as.numeric(pk1$PCSTRESN)
    yL2<-paste0(pk1$PCTESTCD[1],": ", pk1$PCTEST[1])
    p2<-xyplot(PCSTRESN~PCDY|pk1[,byVar], data=pk1, groups=SUBJID,
               ylab=yL2, xlab="Day of Study",type=c('o','g'), col="black", pch=16, cex=0.6,
               scales=list(x=list(cex=0.5, rot=90), y=list(cex=1)))
    
    p12<-doubleYScale(p2,p1, style1=0, style2=4,add.ylab2=TRUE, 
                      main=paste(yL2,"and",yL))
    
    return(p12)
  }
  
}

if(TRUE){#safety review
  # Create heatmap
  doPlot2<-function(mydata, usubjid1, ori='portrait'){
    newdata<-subset(mydata, usubjid==usubjid1)
    # Create a variable, for.labels, that will be used to match 
    # colors to labels in the legend
    if('cmdecod'%in%colnames(newdata)){
      newdata$cmdecod[is.na(newdata$cmdecod)]<-""
    } else {
      newdata$cmdecod <- ""
    }
    if('aedecod'%in%colnames(newdata)){
      newdata$aedecod[is.na(newdata$aedecod)]<-""
      newdata$aesev[is.na(newdata$aesev)]<-""
    }else{
      newdata$aedecod<-""
      newdata$aesev<-""
    }
    if(!'lbnrind'%in%colnames(newdata)){
      newdata$lbnrind <- ""
    }
    newdata<-mutate(newdata, cmdecod=CA(cmdecod),
                    aedecod=CA(aedecod),
                    for.labels=ifelse(domain=="AE" & aesev=="", as.character(aetoxgr),#"NA",
                                      ifelse(domain=="LB" & lbnrind=="", "NA",
                                             ifelse(domain=="CM", "CM",
                                                    ifelse(domain=="AE", as.character(aesev),
                                                           ifelse(domain=="LB", as.character(lbnrind),""))))))
    sVar<-c("days","daystr",'daysend',"cmstrdy","cmenddy")
    if(all(sVar%in%colnames(newdata)))
      newdata<-newdata[rowSums(newdata[,sVar]< -20, na.rm=TRUE)==0,]
    # Create factors to set the measurement order on the y-axis
    newdata<-newdata[order(newdata$domain, newdata$aedecod, 
                           newdata$cmdecod, newdata$param),]
    newdata$aedecod<- factor(newdata$aedecod, ordered=TRUE)
    newdata$cmdecod<- factor(newdata$cmdecod, ordered=TRUE)
    newdata$param<- factor(newdata$param, ordered=TRUE)
    
    have.lb<- ifelse("LB" %in% unique(newdata$domain), 1, 0)
    have.cm<- ifelse("CM" %in% unique(newdata$domain), 1, 0)
    have.ae<- ifelse("AE" %in% unique(newdata$domain), 1, 0)
    
    if(ori=="portrait"){
      bp<- ggplot()+
        theme(text=element_text(size=6),
              axis.text.y=element_text(angle=15, vjust=0.5))+
        facet_grid(domain ~ ., scales = "free", space = "free")+
        scale_y_discrete(expand=c(0.05, 0.6))+
        scale_color_manual("Color Coding", 
                           values=c("NA"="tan4", "0"="gray60", "1"="green", "2"="yellow", 
                                    "3"="orange", "4"="red", 
                                    "ABSENT"="gray40","MILD"="green","MODERATE"="yellow","SEVERE"="red",
                                    "LOW"="blue","NORMAL"="gray50",  "HIGH"="red3", "CM"="skyblue4"),
                           breaks=c("NA", "0", "1", "2", "3", "4", 
                                    "ABSENT","MILD","MODERATE","SEVERE",
                                    "LOW", "NORMAL", "HIGH", "CM"))+
        labs(title=paste("Patient Profile for ", usubjid1, sep=""),
             x="Study Day", y="Test Name", fill="Color Coding")+
        theme_bw()+
        theme(legend.position="bottom")+
        theme(legend.key=element_blank())+
        theme(legend.text=element_text(size = 10))+
        theme(legend.title=element_text(size = 10, face="italic"))+
        theme(panel.background=element_rect(fill='white', colour='black'))+
        theme(panel.border=element_blank(), axis.line=element_line(colour="black"))+
        theme(axis.text.x=element_text(size=rel(0.9)))+
        theme(axis.text.y=element_text(size=rel(0.9)))+
        theme(axis.title.x=element_text(size=11))+
        theme(axis.title.y=element_text(size=11))+
        theme(plot.title=element_text(size=12, face="bold"))
      
      if(have.ae==1){
        bp<- bp +
          geom_segment(data=subset(newdata,domain=="AE"), 
                       aes(x=daystr, xend=daysend, y=aedecod, yend=aedecod, color=for.labels), 
                       size=1)+
          geom_point(data=subset(newdata, domain=="AE"), 
                     mapping=aes(x=daystr, y=aedecod, color=for.labels), shape=24, size=2)+
          geom_point(data=subset(newdata, domain=="AE"), 
                     mapping=aes(x=daysend, y=aedecod, color=for.labels), 
                     shape=25, size=2, fill="grey")
      }
      
      if(have.cm==1){
        bp<- bp +
          geom_segment(data=subset(newdata,domain=="CM"), 
                       aes(x=cmstrdy, xend=cmenddy, y=cmdecod, yend=cmdecod, color=for.labels), 
                       size=1)+
          geom_point(data=subset(newdata, domain=="CM"), 
                     mapping=aes(x=cmstrdy, y=cmdecod, color=for.labels), 
                     shape=24, size=2)+
          geom_point(data=subset(newdata, domain=="CM"), 
                     mapping=aes(x=cmenddy, y=cmdecod, color=for.labels), 
                     shape=25, size=2, fill="skyblue4")
      }
      
      if(have.lb==1){
        bp<- bp +
          geom_text(data=subset(newdata,domain=="LB"), 
                    aes(x=days, y=param, label=lbstresn, color=for.labels, angle=35), 
                    size=2.5)
      }
      
    }else{ #landscape
      bp<- ggplot()+
        theme(text=element_text(size=6),
              axis.text.x=element_text(angle=0, vjust=0.5))+
        facet_grid(.~ domain , scales = "free", space = "free")+
        scale_x_discrete(expand=c(0.05, 0.6))+
        scale_color_manual("Color Coding", 
                           values=c("NA"="tan4", "0"="gray60", "1"="green", "2"="yellow", 
                                    "3"="orange", "4"="red", 
                                    "ABSENT"="gray40","MILD"="green","MODERATE"="yellow","SEVERE"="red",
                                    "LOW"="blue","NORMAL"="gray50",  "HIGH"="red3", "CM"="skyblue4"),
                           breaks=c("NA", "0", "1", "2", "3", "4", 
                                    "ABSENT","MILD","MODERATE","SEVERE",
                                    "LOW", "NORMAL", "HIGH", "CM"))+
        labs(title=paste("Patient Profile for ", usubjid1, sep=""),
             y="Study Day", x="Test Name", fill="Color Coding")+
        theme_bw()+
        theme(legend.position="bottom")+
        theme(legend.key=element_blank())+
        theme(legend.text=element_text(size = 10))+
        theme(legend.title=element_text(size = 10, face="italic"))+
        theme(panel.background=element_rect(fill='white', colour='black'))+
        theme(panel.border=element_blank(), axis.line=element_line(colour="black"))+
        theme(axis.text.x=element_text(size=rel(0.9)))+
        theme(axis.text.y=element_text(size=rel(0.9)))+
        theme(axis.title.x=element_text(size=11))+
        theme(axis.title.y=element_text(size=11))+
        theme(plot.title=element_text(size=12, face="bold"))
      
      if(have.ae==1){
        bp<- bp +
          geom_segment(data=subset(newdata,domain=="AE"), 
                       aes(y=daystr, yend=daysend, x=aedecod, xend=aedecod, color=for.labels, angle=90), 
                       size=1)+
          geom_point(data=subset(newdata, domain=="AE"), 
                     mapping=aes(y=daystr, x=aedecod, color=for.labels), shape=24, size=2)+
          geom_point(data=subset(newdata, domain=="AE"), 
                     mapping=aes(y=daysend, x=aedecod, color=for.labels), 
                     shape=25, size=2, fill="grey")
      }
      
      if(have.cm==1){
        bp<- bp +
          geom_segment(data=subset(newdata,domain=="CM"), 
                       aes(y=cmstrdy, yend=cmenddy, x=cmdecod, xend=cmdecod, color=for.labels, angle=90), 
                       size=1)+
          geom_point(data=subset(newdata, domain=="CM"), 
                     mapping=aes(y=cmstrdy, x=cmdecod, color=for.labels), 
                     shape=24, size=2)+
          geom_point(data=subset(newdata, domain=="CM"), 
                     mapping=aes(y=cmenddy, x=cmdecod, color=for.labels), 
                     shape=25, size=2, fill="skyblue4")
      }
      
      if(have.lb==1){
        bp<- bp +
          geom_text(data=subset(newdata,domain=="LB"), 
                    aes(y=days, x=param, label=lbstresn, color=for.labels, angle=90), 
                    size=2.5)
      }
      
    }
    
    
    return(bp)
  }
  
  
}
