
options(stringsAsFactors=FALSE)

DD<<-indataset
  try(dm<<-DD[["dm.csv"]]) #must have
  try(lb<<-DD[["lb.csv"]]) #must have
  try(ae<<-DD[["ae.csv"]])
  try(cm0<<-DD[["cm.csv"]])
  try(ex<<-DD[["ex.csv"]])
  try(pc<<-DD[["pc.csv"]])
  try(pf0<<-DD[["pf.csv"]] ) #biomarker: rtPCR datasets for PD effect
  try(pf_ref<<-DD[["pf_ref.csv"]])

  if(is.null(dm)) try(dm<<-DD[["dm.sas7bdat"]]) #must have
  if(is.null(lb)) try(lb<<-DD[["lb.sas7bdat"]]) #must have
  if(is.null(ae)) try(ae<<-DD[["ae.sas7bdat"]])
  if(is.null(ex)) try(ex<<-DD[["ex.sas7bdat"]])
  if(is.null(cm0)) try(cm0<<-DD[["cm.sas7bdat"]])
  if(is.null(pc)) try(pc<<-DD[["pc.sas7bdat"]])
  if(is.null(pf0)) try(pf0<<-DD[["pf.sas7bdat"]] ) #biomarker: rtPCR datasets for PD effect
  if(is.null(pf_ref)) try(pf_ref<<-DD[["pf_ref.sas7bdat"]])

  if(is.null(dm))  try(dm<<-DD[grepl('dm.',tolower(names(DD)))][[1]])
  if(is.null(lb))  try(lb<<-DD[grepl('lb.',tolower(names(DD)))][[1]])
  if(is.null(ae))  try(ae<<-DD[grepl('ae.',tolower(names(DD)))][[1]])
  if(is.null(cm0)) try(cm0<<-DD[grepl('cm.',tolower(names(DD)))][[1]])
  if(is.null(pc))  try(pc <<-DD[grepl('pc.',tolower(names(DD)))][[1]])

  inh.sel<-grepl("inhibition", names(DD))
  if(any(inh.sel)){
    try(inhib<<-DD[inh.sel][[1]])
  }else {inhib<<-NULL}

  bm1s1<-grepl("biomarker_cst", names(DD), fixed=TRUE)
  if(any(bm1s1)){bm1<<-DD[bm1s1][[1]]}else{bm1<<-NULL}
if(!is.null(bm1)){
  
  colnames(bm1)<-gsub(" ", "", colnames(bm1), fixed=TRUE)
  colnames(bm1) <-toupper(gsub(".", "", colnames(bm1), fixed=TRUE))
  selcol<-toupper(c("patientid", "colldate", "referencetimingcode",
    "testname", "result", "unitcode"))
  if(any(!selcol%in%colnames(bm1))){
      selcol<-toupper(c("PATIENTNUMBER", "DATECOLLECTED", "TIMECOLLECTED",
      "TESTabbrev", "LABRESULT", "RESULTUOM"))
  }
  if(any(!selcol%in%colnames(bm1))){
      selcol<-toupper(c("PATIENTID", "COLLDATE", "COLLTIME",
      "TESTNAME", "RESULT", 'LBSTRESU'))
      #bm1$LBSTRESU<-''
  }
  bm1<-bm1[,selcol]
  selColnm<-c('SUBJID', 'LBDTC', 'LBELTM','LBTESTCD','LBSTRESN', 'LBSTRESU')
  colnames(bm1)<-selColnm

  bm1<-bm1[!is.na(bm1$SUBJID),]
  bm1<-merge(unique(dm[,c("STUDYID","USUBJID",'SUBJID')]), bm1, by="SUBJID", all.y=TRUE)
  bm1<-bm1[,!colnames(bm1)%in%"SUBJID"]
  bm1$LBTEST<-bm1$LBTESTCD

  lbcol<-colnames(lb)
  n_lbcol<-lbcol[!lbcol%in%colnames(bm1)]
  tmdd<-matrix(NA, nrow=nrow(bm1), ncol=length(n_lbcol))
  colnames(tmdd)<-n_lbcol
  bm1<-data.frame(bm1, tmdd)
  bm1<-bm1[,lbcol]
  bm1$LBSPEC<-"biomarkers"
  bm1$LBSTRESN<-as.numeric(as.character(bm1$LBSTRESN), na.rm=TRUE)
  bm1$LBTEST<-gsub(" ", "", bm1$LBTEST)
  bm1$LBTESTCD<-gsub(" ", "", bm1$LBTESTCD)  
  bm1$LBDTC<-as.Date(bm1$LBDTC, format="%d-%b-%y")

  #get lbdy value
  #dat<-unique(lb[grepl("CYCLE1 DAY1",lb$LBTPTREF),c("USUBJID","LBDTC")])
  head(dat<-unique(lb[lb$LBDY==1,c("USUBJID","LBDTC")]))
  dat<-dat[!is.na(dat$LBDTC),]
  dat<-tapply(as.Date(dat$LBDTC,format="%Y-%m-%d"),   
    dat$USUBJID, function(x){as.character(min(x,na.rm=TRUE))})
  #colnames(dat)<-c("USUBJID","startdate")
  dat<-data.frame(USUBJID=names(dat), startdate=dat)
  dat$startdate<-as.Date(dat$startdate, format="%Y-%m-%d")
  dat<-unique(dat)
  #remove patients with different starting date of study
  frq_dat<-tapply(rep(1,nrow(dat)), dat$USUBJID, sum)
  dat<-dat[dat$USUBJID%in%names(frq_dat)[frq_dat==1],]
  lb$LBDTC<-as.Date(lb$LBDTC, format="%Y-%m-%d")

  lb<<-rbind(lb, bm1[,lbcol])

  lb<<-merge(lb, dat,  all.x=TRUE, by="USUBJID")
  lb$LBDTC<-as.Date(lb$LBDTC, format="%Y-%m-%d")
  lb$LBDY<- as.numeric(lb$LBDTC-lb$startdate+1)

  } #if biomarker_cst.csv is available

#for the biomarker rtPCR data
if(!is.null(pf0)){
  #use MEANCT
  pf0<<-pf0[pf0$PFTESTCD_qPCR=="MEANCT",]
  pf_ref$LBDAT<-paste(pf_ref$SPECIMEN_COLL_DATE, pf_ref$SPECIMEN_COLL_TIME,sep="-")
  pf1<<-merge(pf0, pf_ref[,c('SUBJID','LAB_REQ_NUMBER','LBDAT')],
    by=c('SUBJID','LAB_REQ_NUMBER'))
  dim(pf1<-pf1[pf1$LBDAT!='',])
  dim(hh<-unique(pf1[,c("SUBJID","LBDAT")]))
  hh$LBDAT<-strptime(hh$LBDAT,format='%d-%b-%y-%H:%M')
  hh<-hh[!is.na(hh$LBDAT),]#the useful date and time in gene expression

  #get date for LBDY==1
  dat<-unique(lb[!is.na(lb$LBDY)&lb$LBDY==1, c("USUBJID","LBDTC")])
  dat$LBDTC<-as.Date(dat$LBDTC,format=c("%Y-%m-%d"))
  dat<-unique(dat)
  nchar1<-nchar(dat$USUBJID)
  dat$SUBJID<-substr(dat$USUBJID, nchar1-3, nchar1)

  #add the LBSPEC column
  if(!'LBSPEC'%in%colnames(pf1)){
    pf1$LBSPEC<-'PF domain'
  }
  
  pf2<<-merge(pf1, dat, by='SUBJID')
  pf2$LBDAT<-strptime(pf2$LBDAT,format='%d-%b-%y-%H:%M')
  pf2<-pf2[!is.na(pf2$LBDAT),]
  pf2$LBDTC<-strptime(pf2$LBDTC,format=c("%Y-%m-%d"))
  pf2$LBDY<-as.numeric(pf2$LBDAT-pf2$LBDTC)/24+1
 
  colnames(pf2)<-gsub("PF","LB", colnames(pf2))
  colnames(pf2)<-gsub("_qPCR","", colnames(pf2))
  pf2$LBTEST<-pf2$LBTESTCD<-pf2$LBGENROI
  col1<-colnames(lb)
  pf3<<-pf2[,colnames(pf2)%in%col1]
  pf3$LBSTRESN<-pf3$LBORRES
  pf3$LBSTRESN[gsub(" ", "", pf3$LBSTRESN,fixed=TRUE)==""]<-40
  pf3$LBSTRESU<-pf3$LBORRESU
  col2<-col1[!col1%in%colnames(pf3)]
  pf4<-cbind(pf3,matrix(NA, nrow=nrow(pf3),ncol=length(col2)))
  colnames(pf4)<-c(colnames(pf3),col2)

  lb$LBDTC<-as.Date(lb$LBDTC, format="%Y-%m-%d")
  pf4$LBDTC<-as.Date(pf4$LBDTC, format="%Y-%m-%d")
  lb<<- rbind(lb, pf4)

}

#merge the data
noCol1 <- c("STUDYID",'SUBJID', 'DOMAIN')
comDD<<-merge(dm[,!colnames(dm)%in%noCol1], 
              lb[,!colnames(lb)%in%noCol1], 
              all.y=TRUE, by=c("USUBJID"))
comDD$USUBJID<-as.character(comDD$USUBJID)
comDD$SUBJID <- substring(comDD$USUBJID, nchar(comDD$USUBJID)-3)

print(head(comDD))

  if(is.null(comDD$startdate)){
    comDD$startdate<-as.Date(NA)
  }else{
    comDD$startdate<-as.character(comDD$startdate)
  }

  if(TRUE){#re-calculate lbdy
    na.sub<-unique(comDD$SUBJID[is.na(comDD$LBDY)])
    for(i in na.sub){
      sel1<-comDD$SUBJID==i
      sel2<-sel1&comDD$LBDY==1  #&comDD$LBELTM>0&comDD$LBELTM<15
      #try(comDD$startdate[sel1]<-min(as.Date(comDD$LBDTC[sel2]),na.rm=T))
      tm<-class(try(stt<-min(as.Date(comDD$LBDTC[sel2]),na.rm=T)))
      #try(comDD$LBDY[sel1]<-as.numeric(comDD$LBDTC[sel1]-comDD$startdate[sel1]+1))
      #if(tm!="try-error")
      comDD<-comDD[!is.na(comDD$LBDY),]
    }
    
    if(!is.null(comDD$LBELTM)){
      comDD$LBELTM<-gsub("P", "", comDD$LBELTM)
      comDD$LBELTM<-gsub("T", "", comDD$LBELTM)
      comDD$LBELTM<-as.numeric(comDD$LBELTM)/24
      comDD$LBELTM[is.na(comDD$LBELTM)]<- 0
      comDD$LBELTM[comDD$LBELTM>=1] <- 0
      comDD$LBDY<-comDD$LBDY+comDD$LBELTM
      comDD$LBELTM<-NULL
  }
  comDD$SUBJID<-paste(comDD$ARMCD, comDD$SUBJID)
  ss1<-is.na(comDD$LBSTRESN)&!is.na(comDD$LBSTRESC)
  comDD$LBSTRESN[ss1]<-as.numeric(gsub("<","", comDD$LBSTRESC[ss1]))

}

if(TRUE){#revise comDD data
  comDD$LBTEST<-paste0(as.character(comDD$LBTEST), ' (',
    gsub(' ', '', as.character(comDD$LBSTRESU), fixed=TRUE), ')')
   fill.na <- !is.na(comDD$LBSTRESN)&gsub(' ','',comDD$LBSTRESN)!=''
   comDD$LBSTRESN[!fill.na]<-comDD$LBORRES[!fill.na]
   comDD$LBSTRESU[!fill.na]<-comDD$LBORRESU[!fill.na]
  comDD<-comDD[!is.na(comDD$LBDY)&comDD$LBDY>-50,]
  comDD$LBSPEC[comDD$LBSPEC=='']<-'blank'
  comDD$LBTEST[comDD$LBTEST=='']<-'blank'
}

if(!is.null(pc)){#concatenate the data, get PK/LB combinatin long table
  pk<<-merge(dm, pc, by=c("USUBJID","STUDYID"))
  ss2<-is.na(pk$PCSTRESN)&!is.na(pk$PCSTRESC)
  pk$PCSTRESN[ss2]<-as.numeric(gsub("<","", pk$PCSTRESC[ss2]))
  pk$SUBJID<-paste(pk$ARMCD, pk$SUBJID)
  comDD2<<-comDD[comDD$SUBJID%in%unique(pk$SUBJID),]
}

if(TRUE){#define demographic variables of boxplot selection
  dmVar<-c("SITEID","SEX","RACE","ETHNIC","ARMCD","COUNTRY")
  dmVar<<-dmVar[dmVar%in%colnames(comDD)]
  uni_subjid<<-sort(unique(comDD$SUBJID))
  uni_armcd<<-unique(comDD$ARMCD)
  uni_lbspec<<-unique(comDD$LBSPEC)
  #uni_lbtest<-unique(comDD$LBTEST);#depend on uni_lbspec
  keepColumn<<-c("SUBJID","ARM","ARMCD","LBDY","LBSTRESN","LBSTRESU","LBTEST","LBTESTCD")
  sel_col<<-colnames(comDD)[!colnames(comDD)%in%c(keepColumn,"DOMAIN.x","DOMAIN.y")]

}



if(TRUE){#get data for safety review
  if(!is.null(ae) && !"SUBJID"%in%colnames(ae)){
    ae$SUBJID <- substring(ae$USUBJID, nchar(ae$USUBJID)-3)
  }
  
  
  lb2.col <- c("LBDY","USUBJID","SUBJID","LBSTRESN","LBNRIND", "LBTESTCD","LBTEST","LBSTRESU","LBSPEC")
  lb2.col <- lb2.col[lb2.col%in%colnames(comDD)]
  lb2<-comDD[,lb2.col]
  colnames(lb2)[1]<-c("days")
  lb2$DOMAIN<-'LB'
  lb2$param<-paste0(lb2$LBTESTCD,": ", lb2$LBTEST, " (", lb2$LBSTRESU, ")")
  lb2.col2 <- c("DOMAIN","USUBJID","LBSTRESN","LBNRIND","param", "days", "SUBJID","LBSPEC")
  lb2.col2 <- lb2.col2[lb2.col2%in%colnames(lb2)]
  lb2<-lb2[, lb2.col2]
  

  if(!is.null(ae)){
    ae2<-ae[,c("AESTDY", "AEENDY", "DOMAIN","USUBJID", "AEDECOD", "AESEV", "AETOXGR")]  
    colnames(ae2)[1:2]<-c("daystr","daysend")
    ae2<-merge(ae2, unique(lb2[,c("USUBJID","SUBJID")]), by="USUBJID")
  }else{ae2 <- NULL}
  
  if(!is.null(cm0)){
    cm2<-cm0[,c("CMSTDY","CMENDY","DOMAIN","USUBJID","CMDECOD")]
    colnames(cm2)[1:2]<-c("CMSTRDY","CMENDDY")
    cm2<-merge(cm2, unique(lb2[,c("USUBJID","SUBJID")]), by="USUBJID")
  }else{cm2 <- NULL}

  by123<-c("DOMAIN","USUBJID", "SUBJID")
  if(!is.null(ae)){
    mydata<<-merge(lb2, ae2, by=by123, all=TRUE, sort=FALSE)
    mydata<-mydata[!(mydata$DOMAIN=="AE"&is.na(mydata$daystr)),]
  }else{
    mydata<<-lb2;
  }
  
  if(!is.null(cm2)){
    mydata<<-merge(mydata, cm2, by=by123, all=TRUE, sort=FALSE)
    mydata<-mydata[!(mydata$DOMAIN=="CM"&is.na(mydata$CMSTRDY)),]
  }
  
  mydata$USUBJID<-mydata$SUBJID
  mydata<-mydata[!(mydata$DOMAIN=="LB"&is.na(mydata$days)),]
  #mydata<-mydata[!(mydata$DOMAIN=="AE"&is.na(mydata$daystr)),]
  
  colnames(mydata)<-tolower(colnames(mydata))

  subject.id <<- unique(mydata$usubjid)  
  mydata<<-mydata
  
  #baseline lb data
  day1befor<-strptime(comDD$RFSTDTC,format='%Y-%m-%dT%H:%M')>=strptime(comDD$LBDTC,format='%Y-%m-%dT%H:%M')
  dim(mylb.log<<-unique(comDD[,c('SUBJID', "LBTEST", "LBSPEC", "LBSTRESN", "LBDY")]))
  names(mylb.log)<-c('usubjid','param','lbspec','lbstresn', 'lbdy')
  dim(mylb<<-unique(comDD[comDD$LBDY<1|(!is.na(day1befor)&day1befor),c('SUBJID', "LBTEST", "LBSPEC", "LBSTRESN")]))
  names(mylb)<-c('usubjid','param','lbspec','lbstresn')
  #dim(mylb<<-mydata[mydata$domain=='LB',])
  
  #baseline lb data and ae data
  if(!is.null(ae)){
    dim(myae<<-unique(mydata[mydata$domain=='AE',c("usubjid", "aedecod", "aetoxgr")]))
    mylbae.log<<-merge(mylb.log, myae)
    dim(mylb<<-aggregate(lbstresn~usubjid+param+lbspec, data=mylb, FUN=mean, na.rm=TRUE))
    mylbae<<-merge(mylb, myae)
  }else{mylbae<-NULL}
  
  if(FALSE){
    # sort(unique(mylb.log$param))
    #head(tm<-mylb.log[grepl('CK18', mylb.log$param),])
    head(tm<-mylb.log[mylb.log$param=="Amyloid, Beta (pg/mL)",])
    # sort(unique(myae$aedecod))
    ae1term<-'VOMITING'
    ae.sub<-unique(myae$usubjid[myae$aedecod==ae1term])
    tm$subgrp<-'other'
    tm$subgrp[tm$usubjid%in%ae.sub]<-ae1term
    dim(tm2<-unique(tm[,c('lbdy','lbstresn','usubjid','subgrp', 'param')]))
    qplot(lbdy, lbstresn, data=tm2, group=usubjid, xlim=c(0,100))+
      geom_line()+geom_point()#+facet_grid(subgrp~param)+geom_hline(yintercept=c(500, 250), colour=c('red', 'blue'))
  }
  #ggplot(dynamicData, aes(x=subgroup, y=lbstresn)) +geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+guides(colour=FALSE)+geom_point()+theme(plot.margin=unit(c(1,1,1,3),'cm')) +coord_flip()
}


if(TRUE){#load relative libraries.
  library(beanplot)

  library(corrplot)  #create heatmap for correlation matrix#

  library(grid)
  library(gridExtra)
  library(ggplot2)

  library(Heatplus)  #creating customized heatmaps with covariate profile#
  library(Hmisc)

  library(lattice)
  library(latticeExtra)
  library(lmtest)

  library(multtest)

  library(plotrix)   #use 2D-heatmap#
  library(plyr)
  library(psych)

  library(reshape)
  library(reports)
  library(rtf)	 #creating rtf files#

  library(shiny)
  library(sp)
  library(survival)  #time to event analysis#

  library(vioplot)

  library(xtable)
}

  mydata <<- mydata

if(!is.null(inhib) & !is.null(dm)){
  inhib1<-inhib[inhib$blfl==0,]
  inhib2<<-merge(unique(inhib1[,c('SUBJID', 'inhib2', 'geneid')]), unique(dm[,c('SUBJID', 'ARMCD', 'ARM')]))
  uni_gene_inhib2<<-sort(unique(inhib2$geneid))
  uni_arm_inhib2<<-sort(unique(inhib2$ARMCD))
  #
#  dynamicData<<-inhib2[inhib2$geneid=='NRARP',]
#  qplot(ARMCD, inhib2, data=dynamicData, geom='boxplot')+guides(colour=FALSE)+theme(plot.margin=unit(c(1,2,1,1),'cm'))+ylab('%Inhibition')+xlab("ARM")+geom_abline(intercept=50, slope=0)+geom_point(colour='red')
}

