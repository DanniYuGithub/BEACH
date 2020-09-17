if(TRUE){#header
  #*****************************************************************************
  #Eli Lilly and Company - GLOBAL STATISTICAL SCIENCES - PROGRAM
  #CODE NAME (required)                : LightON1_meta.R
  #PROJECT NAME (required)             : geneAtlas
  #STUDY NUMBER                        : 
  #DESCRIPTION (required)              : update / organize / integrate data
  #SPECIFICATIONS(optional)            : 
  #VALIDATION TYPE (required)          : 
  #INDEPENDENT REPLICATION (optional)  : 
  #ORIGINAL CODE (required)            : 
  #EXTERNAL CODE FILES THAT ARE NOT 
  #RE-USEABLE CODE MODULES             : 
  #SOFTWARE/VERSION# (required)        : R Version 3.5.0
  #INFRASTRUCTURE                      : R studio server
  #DATA INPUT location                 : gene expression and mutation
  #SDD: 
  #                                    
  #OUTPUT location: 
  #
  #
  #R SCRIPT location
  #
  #
  #-----------------------------------------------------------------------------
  #     #Author                 Code History Description             Date Change
  #---- ------------       ---------------------------------------  ------------
  #1.1  Danni Yu            Code Creation                             02/01/2019
  #1.2  Danni Yu            Revise code to for feather data analysis  02/28/2019
  #*****************************************************************************
}

#dataset was integrated by LightON1_integrate.R 
#set data access permition as chmod a+r+x-w LightON1



#define paths and load clinical & mutation datasets
if(F){#tbd after run LightON1_integrate.R
  projP <- "TBD_by_user"
  progfile<-'LightON1_meta.R'
  
  #define data paths  
  dir(frompath<<-"data")
  smfile.nm<-'projects_summary.csv'
  gnameFile<-file.path(frompath, 'gene_ensemble.Rdata')
  
  #the selected list of datasets under project selections
  sprjL<<-list(clin=NULL, mut=NULL, gexp_array=NULL, gexps_rnase=NULL,
               mut_cna=NULL, pexp_antib=NULL)
  
  #innitiate global vectors that will be updated in select.prj function
  uni.mutg<<-uni.gexpA<<-uni.gexpR <<- uni.pexpA <<- uni.gene <<-NULL
  uni.gene.default<<-c('RB1', 'KMT2C', 'KMT2D', 'ARID1A')
  grpVar.clin<<-numVar.clin<<-sprgVar<<-NULL
  grpVar.clin0<<-c('proj_id','proj_name','sex','vital_status', 
                   'agec_at_enrollment','agec_at_diagnosis',
                   #'agec_at_last_followup',
                  'disease_status_last_followup', 'tumour_stage_at_diagnosis',
                  'tumour_stage_system', 'tumour_histological_type', 
                  'tumour_stage', 'specimen_type', 'specimen_treatment_type',
                  'specimen_processing', 'first_therapy_type',
                  'frist_therapy_therapeutic_intent',  'first_therapy_response',
                  'second_therapy_type', 'second_therapy_therapeutic_intent',
                  'second_therapy_response', 'other_therapy', 'dummy','db')
  numVar.clin0<<-c('age_at_enrollment', 'age_at_diagnosis',
                   'age_at_last_followup',
                   'relapse_interval', 'survival_time',
                  'specimen_interval', 'first_therapy_duration',
                  'first_therapy_start_interval', 'second_therapy_start_interval',
                  'second_therapy_duration')
  sgrpVar0<<-c('mutation_type','aa_mutation','cds_mutation','chromosome')
  
  #global mapping variable
  conVar1 <<- "submitted_id"
  #sample level mapping variable
  conVar2<<-"submitted_sample_id"
    
}

#select projects from the project summary table
if(F){#tbd
  require(survival)
  require(plyr)
  require(dplyr)
  require(gplots)
  library(data.table)
  library(feather)
  library(glmnet)
  
  #get gene match to chromosom
  load(gnameFile)
  g.chrom<<-distinct(data.table(Hugo_Symbol=g_list3$Hugo_Symbol,
                       chromosome=g_list3$chromosome_name))
  
  #read in the project summary table
  prjsm<<-read.csv(file.path(frompath,smfile.nm), 
                   h=T, stringsAsFactors=F)
  
  #function    sel.id="75,5"
  select.prj<<-function(sel.id, #input$text, the 'id' value in prjsm
                        prjTbl=prjsm #the project summary table
  ){
    #setup an internal function for read in data
    rbF<-function(kw='_clin.feather', 
                  fp=tp, #file path and file names
                  pj=NULL,   #project id and project name
                  subj=NULL){
      t1<-fp[[2]][grepl(kw, fp[[2]])] #get data name
      if(length(t1)==0){return(NULL)}
      tm<-read_feather(file.path(fp[[1]], t1))
      if(!is.null(pj)){
        tm$proj_id<-pj[1]
        tm$proj_name<-pj[2]
        tm$db<-pj[3]
        tm$db_id<-pj[4]
        tm$dummy<-'integrate'
      }
      tm<-tm[colSums(!is.na(tm))>0]#tm is a data.table object, not data.frame
      if(!is.null(subj)){
        tm<-tm[tm$submitted_sample_id%in%subj,]
        if(nrow(tm)==0) return(NULL)
      }
      return(tm)
    }
    #convert id1 into a numeric vector
    if(length(sel.id)>0 && grepl('all',tolower(sel.id))){sel.id<-prjTbl$id}
    if(is.null(sel.id)||gsub(' ','',sel.id)==''){sel.id<-1}
    if(length(sel.id)>0 && !is.numeric(sel.id)){
      sel.id<-as.numeric(strsplit(as.character(sel.id),split=',')[[1]])
      sel.id<-sel.id[!is.na(sel.id)]
    }
    if(length(sel.id)==0){  sel.id<-1 }
    Oclin<-Omut<-OgexpA<-OgexpR<-OmutCNV<-OpexpA<-list()
    for(i in 1:length(sel.id)){
      p0<-prjTbl[sel.id[i],]
      tmp<-file.path(frompath, paste(p0[,c('db','db_id')],collapse='/'))
      tp<-list(tmp, dir(tmp))
      pj1<-c(p0$proj_id[1], p0$proj_name[1], p0$db[1], p0$db_id[1])
      Oclin[[i]]<-rbF(kw='_clin.feather', pj=pj1) #read in clin data
      ds<-rbF(kw='_datasum.feather', fp=tp)  #read in data summary data
      Omut[[i]]<-rbF(kw='_mut.feather')#read in mut&cnv&fus
      # #pull out the expression datasets match to mutational datasets
      # sub1<-unique(mut1$submitted_sample_id)
      if(F){#this takes additional 2 minutes to load
        OgexpA[[i]]<-rbF(kw='_mut_gexp_array.feather')
        OgexpR[[i]]<-rbF(kw='_mut_gexp_rnaseq.feather')
        OpexpA[[i]]<-rbF(kw='_mut_pexp_antib.feather')
      }
    }
    #remove columns with missing data
    rbL<-function(dd){as.data.frame(rbindlist(dd, fill=TRUE))}
    mut<-rbL(Omut); s1.mut<<-unique(mut[,conVar2])#global variable for sampleID
    clin<-rbL(Oclin); clin<-clin[clin[,conVar2]%in%s1.mut,]
    outL<-list(clin=clin, mut=mut)
               # gexpA=rbL(OgexpA),gexpR=rbL(OgexpR), pexpA=rbL(OpexpA))
    #lapply(outL, dim)
    #update the global vectors for gene names
    if(!is.null(outL$mut)){uni.mutg<<-unique(outL$mut$Hugo_Symbol)}

    uni.gene<<-sort(uni.mutg,method="quick")
    uni.gene.default<<-uni.gene.default[uni.gene.default%in%uni.gene]
    if(length(uni.gene.default)==0){uni.gene.default<<-uni.gene[1]}

    #update the global vectors for selecting variable names    
    grpVar.clin<<-grpVar.clin0[grpVar.clin0%in%colnames(outL$clin)]
    numVar.clin<<-numVar.clin0[numVar.clin0%in%colnames(outL$clin)]
    sgrpVar<<-sgrpVar0[sgrpVar0%in%colnames(outL$mut)]
    
    #change clin data into numberic
    for(c1 in numVar.clin){
      outL$clin[,c1]<-as.numeric(outL$clin[,c1])
    }
      
    return(outL)
  }
  #example list
  #tt1<-Sys.time(); sprjL<<-select.prj('7,8,9,10,11'); print(Sys.time()-tt1)
  #tt1<-Sys.time(); sprjL<<-select.prj('All'); print(Sys.time()-tt1)
  #lapply(sprjL, dim)
  
  
  #function for data summary heatmap
  see.click.prj<<-function(chk.id=NULL, #input$text, the 'id' value in prjsm
                          prjTbl=prjsm, #the project summary table
                           showMutOnly=T
  ){
    ff1<-function(tt='no project is selected yet.'){
      plot(0~0, axes=F, ylab='',xlab='', col='white')
      mtext(text=tt)
    }
    if(is.null(chk.id)){ ff1(); return(NULL) }
    #choose the last selected id and the load the datasummary data
    chk.id<-chk.id[length(chk.id)]
    p0<-prjTbl[chk.id,]
    tt<-paste0(p0$db, p0$db_id, "_", p0$proj_id,'\nData Availability')
    tmp<-file.path(frompath, paste(p0[,c('db','db_id')],collapse='/'))
    fn<-dir(tmp); fn<-fn[grepl('_datasum.feather',fn)]
    if(length(fn)>0) ds0<-as.data.frame(read_feather(file.path(tmp, fn[1])))
    if(nrow(ds0)>5000){showMutOnly<-T}
    if(showMutOnly){
      if('mut'%in%colnames(ds0)){
        ds0<-ds0[ds0$mut==1,]
      }else{ff1('mut data is not available');return(NULL)}
    }
    #create the heatmap
    c1n<-c(conVar1, conVar2, 'submitted_specimen_id')
    ds1<-t(apply(ds0[, !colnames(ds0)%in%c1n], 2, as.numeric))
    rs<-rowSums(ds1)
    ds1<-ds1[rs>0,]
    colnames(ds1)<-ds0[,c1n[1]]
    rownames(ds1) <- paste(rownames(ds1) , rs[rs>0])
    heatmap.2(ds1, col=c('white','lightblue'), scale='none', labCol=FALSE,
            Rowv=T,  Colv=T, dendrogram='none', margins=c(3,10), offsetRow=0, 
            key=F,xlab='samples', main=tt, trace='none',
            lmat=rbind(c(5,3,4),c(2,1,4)), lwid=c(0.5,30,5), lhei=c(1,6.5))
            #rowsep=1:nrow(ds1), colsep=1:ncol(ds1), sepwidth=c(0.01, 0.005))
  }
  # see.click.prj(chk.id=c(3,75)) #only show the plot for the last click

  see.select.prj<<-function(chk.id=NULL, #input$text, the 'id' value in prjsm
                             prjTbl=prjsm, #the project summary table
                            showMutOnly=T
  ){
    if(length(chk.id)>0 && is.character(chk.id))
      try(chk.id<-as.numeric(strsplit(as.character(chk.id), split=',')[[1]]))
    if(is.null(chk.id) || length(chk.id)<1 || is.na(chk.id)){
      plot(0~0, axes=F, ylab='',xlab='', col='white')
      mtext(text='no project is selected yet.')
      return(NULL)
    }
    tt<-nn<-NULL; dsL<-list();
    c1n<-c("submitted_id", 'submitted_sample_id', 'submitted_specimen_id')
    #choose the last selected id and load the datasummary data
    for(i in 1:length(chk.id)){
      p0<-prjTbl[chk.id[i],]
      tt<-c(tt, paste0(p0$db, p0$db_id, "_", p0$proj_id))
      tmp<-file.path(frompath, paste(p0[,c('db','db_id')],collapse='/'))
      fn<-dir(tmp); fn<-fn[grepl('_datasum.feather',fn)]
      if(length(fn)>0){
        ds0<-as.data.frame(read_feather(file.path(tmp, fn[1])))
        if(showMutOnly){
          if('mut'%in%colnames(ds0)){
            ds0<-ds0[ds0$mut==1,]
          }else{next}
        }
        c2n<-colnames(ds0)[!colnames(ds0)%in%c1n]
        ds0[,c2n]<-apply(ds0[,c2n],2, as.numeric)
        ds0<-ds0[order(rowSums(ds0[,c2n])),]
        dsL[[i]]<-ds0
        nn<-c(nn, nrow(ds0))
      }
    }
    tt1<-paste0('Data Availability\n',paste(tt,collapse=', '))
    lab1<-rep('',sum(nn))
    nn2<-cumsum(nn[-length(nn)])
    lab1[round(nn/2)+c(0,nn2)]<-tt
    ds0<-do.call('rbind.fill', dsL)
    ds0[is.na(ds0)]<-0
    ds1<-t(apply(ds0[, !colnames(ds0)%in%c1n], 2, as.numeric))
    #create the heatmap
    rs<-rowSums(ds1)
    rownames(ds1) <- paste(rownames(ds1) , rs)
    ds1<-ds1[order(rowSums(ds1)),]
    heatmap.2(ds1, col=c('white','lightblue'), scale='none', 
              labCol=lab1, srtCol=0, cexCol=1,
              Rowv=F,  Colv=F, dendrogram='none', margins=c(4,10), offsetRow=0, 
              key=F,xlab='samples', main=tt1, trace='none',
              colsep=nn2, sepcolor='black',
              lmat=rbind(c(5,3,4),c(2,1,4)), lwid=c(0.5,30,3), lhei=c(1.5,4))
  }
  # see.select.prj(chk.id=c(3,75)) #show the integrated plot for all selection

}


#for realtime summary table, visualization and and download. 
if(TRUE){
  
  #function to grab matched gene names for mutation data,
  #and find and return hugo gene names 
  findHugo<<-function(x='CDK7, ARID1A, MYCN', HG.v=uni.gene, 
                      dum=NULL #input$userdo will be used to trigger the update
                      ){
    #x='AAA1,BBB'
    #HG.v: a vector of gene names
    if(length(x)==0)return(NULL)
    x<-strsplit(x, split=',', fixed=T)[[1]]
    x<-toupper(gsub(' ', '', x, fixed=T))
    x2<-lapply(x, function(s){HG.v[grepl(s,HG.v, fixed=T)]})
    y<-sapply(x2, function(s){paste(s,collapse='||')})
    y1<-paste('gene names matched to the search:\n', paste(y, collapse='\n'))
    y2<-x[x%in%HG.v]
    return(list(y1, y2, unlist(x2)))
  } 
  # findHugo("aka, csm,or2b, cdk7, arid1a")
  
  #function to pull out available levels from group variable
  #sprjL$clin is used here
  gLv1<<-function(varName=NA, limit=200, 
                  dum=NULL #input$userdo will be used to trigger the update
                  ){
    if(is.null(varName) || is.na(varName) || 
       gsub(' ', '', varName, fixed=T)==''){
      return(NA)}
    tmp<-sprjL$clin[,varName]
    if(is.numeric(tmp)){
      return(NA)
    }else{
      tmp<-sort(unique(as.character(tmp)),method="quick")
      return(tmp[1:min(length(tmp), limit)])
    }
  }
  #function to pull out available levels from sub group variable
  #sprjL$mut is used here
  gLv2<<-function(varName, geneName, limit=200,
                  dum=NULL #input$userdo will be used to trigger the update
                  ){
    if(is.null(varName) || is.na(varName)){return(NULL)}
    tmp<-sprjL$mut[sprjL$mut$Hugo_Symbol%in%geneName,varName]
    if(is.numeric(tmp)){
      return(NULL)
    }else{
      tmp<-sort(unique(as.character(tmp)),method="quick")
      return(tmp[1:min(length(tmp), limit)])
    }
  }
  #function to get the range of number variable.
  #....tbd
  
  #construct the levels in a listing object x into matrix
  lev2mat<<-function(x){
    if(!is.list(x)){return(data.frame(x='x is not a list in lev2mat'))}
    x<-lapply(x,function(x){x[!is.na(x)]})
    x<-lapply(x,function(x){x[!toupper(x)%in%c('NULL','NA')]})
    len.x<-sapply(x, length)
    x<-x[len.x>0]
    if(length(x)==1){
      o1<-data.frame(x[[1]])
      colnames(o1)<-names(x[1])
      return(o1)
    }
    if(length(x)<1){
      return(NULL)
    }
    x2<-list(); 
    len.x<-sapply(x, length)
    for(m in 1:length(x)){
      x2[[m]]<-rep(rep(x[[m]], each=prod(len.x[-(1:m)])),
                   times=prod(len.x[0:(m-1)]))
    }
    o1<-do.call("cbind",x2)
    colnames(o1)<-names(x)
    return(o1)
  }
  
  #pull mOS from survfit summary table
  survSum1<<-function(sel1){
    if(is.null(sel1)){return('NA(NA, NA)')}
    paste0(round(sel1['median'],2), '(',
           round(sel1['0.95LCL'],2), ', ',
           round(sel1['0.95UCL'],2), ')')
  }
  
  #Function for clinical summary table at subtypes nested in groups Bootstrap
  sgrpBSsum<<-function(
    d0=sdat2.1,   #data to be summarized with OS_mon and OS_cnr and bmk
    d0all=sMut1.sa, #all candidate genes and samples
    size1=n1, m1=b1,
    var1=sel_grpVar,
    mutFreqTh=0.01
  ){
    #group levels
    if(length(var1)==0){
      v1<-'NULL'
    }else if(length(var1)==1){
      v1<-paste(unique(d0[,var1]), collapse=' | ')
    }else{
      v1<-t(apply(d0[,var1],2, function(x){paste(unique(x), collapse=' | ')}))
    }
    #biomarker labels
    v2<-strsplit(d0$labelBmk[d0$labelBmk!=''], split='|', fixed=T)
    v2<-paste(unique(unlist(v2)), collapse='|')
    #get parameters
    N <-length(unique(d0[,conVar1]))
    n <-length(unique(d0[d0$testBmk==1, conVar1]))
    if(n==0){return(NULL)}
    #bootstrap estimate prevalence or mutaion rate
    pctv<-rep(NA, m1)
    a<-rep(0, N)
    a[sample(1:length(a),n)]<-1
    for(i in 1:m1){ #BS samples
      # d0s<-d0[sample(1:nrow(d0),size1, T),]
      # Ni <-length(unique(d0s[,conVar1]))
      # ni <-length(unique(d0s[d0$testBmk==1, conVar1]))
      # pctv[i]<-round(ni/Ni,3)
      b<-sample(a, size1, T)
      pctv[i]<-round(sum(b)/size1,3)
    }
    pctv<-pctv[!(is.na(pctv)|is.nan(pctv))]
    pct<-round(quantile(pctv, prob=c(0.5,0.025,0.975)),3)
    pct0<-round(n/N, 3)
    bs1<-paste0(m1, ' times, size=',size1)
    if(sum(d0$OS_mon, na.rm=T)==0){
      sr6mon_BMKpos<-sr6mon_BMKneg<-mOS1<-mOS0<-HR<-pval<-NA
      f0<-st.f0<-leg1<-col1<-lty1<-NULL
    }else{
      #6 month survival rate
      sr6mon_BMKpos<-round(mean(d0$OS_mon[d0$testBmk==1]>=6,na.rm=T),2);
      sr6mon_BMKneg<-round(mean(d0$OS_mon[d0$testBmk==0]>=6,na.rm=T),2);
      #get mOS values
      d0$bmk <- as.factor(d0$testBmk)
      d00<-distinct(d0[!is.na(d0$OS_mon)&!is.na(d0$OS_cnr),
                       c('bmk','OS_mon', 'OS_cnr', conVar1)])
      if(nrow(d00)<3){f0<-st.f0<-NULL}else{
        try(f0<-survfit(Surv(OS_mon, OS_cnr)~bmk, data=d00))
        try(st.f0 <- summary(f0)$table)
      }
      tb.bmk<-table(d00$bmk)
      if(all(tb.bmk>0)&length(tb.bmk)==2){
        mOS1<-survSum1(st.f0['bmk=1',])
        mOS0<-survSum1(st.f0['bmk=0',])
        col1<-c('blue','red')
        lty1<-c(2,1)
        leg1<-paste(c(tb.bmk['0'],tb.bmk['1']), c('bmk0 mOS=','bmk1 mOS='), 
                    c(mOS0, mOS1))
      }else if('1'%in%names(tb.bmk) && tb.bmk['1']>0){
        mOS1<-survSum1(st.f0); mOS0<-'no obs'
        col1<-c('red')
        lty1<-c(1)
        leg1<-paste(c(tb.bmk['0']), c('bmk0 mOS='), c(mOS0))
      }else if ('0'%in%names(tb.bmk) && tb.bmk['0']>0){
        mOS0<-survSum1(st.f0); mOS1<-'no obs'
        col1<-c('blue')
        lty1<-c(2)
        leg1<-paste(c(tb.bmk['0']), c('bmk0 mOS='), c(mOS0))
      }else{
        mOS0<-mOS1<-'no obs'
        f0<-col1<-lty1<-leg1<-NULL
      }
      #get HR
      if(sum(!is.na(d00$OS_mon))<5){
        HR<-pval<-NA
      }else{
        # hh<<-d0
        t1<-try(f1<-coxph(Surv(OS_mon, OS_cnr)~bmk, data=d00))
        if(class(t1)=='try-error') {
          HR<-pval<-NA
        }else if (any(is.na(f1$coefficients))){
          HR<-pval<-NA
        }else{
          sc.f1<-summary(f1)$coefficients
          HR<-round(sc.f1['bmk1','exp(coef)'],2)
          pval<-round(sc.f1['bmk1','Pr(>|z|)'],2)
        }
      }
    }
    #set up the summary in a row
    otb<-data.frame(v1, sgrp=v2, n=N, n_BMKpos=n, 
                    prevl_BMKpos_95CI=paste0(pct0,' (',pct[2],',',pct[3],')'), 
                    HR=HR, pval=pval, 
                    mOS_BMKpos=mOS1, mOS_BMKneg=mOS0, 
                    SurvivalRate6mon_BMKpos=sr6mon_BMKpos,
                    SurvivalRate6mon_BMKneg=sr6mon_BMKneg)
    outb.c1<-grepl('prevl_BMK_pos_95CI', colnames(otb))
    colnames(otb)[outb.c1]<-paste0('prevl_BMKpos_95CI_n',size1)
    #get gene freqency in the group
    #can be very slow if the mutation numbers are large. !!!!
    # gfreq1<-merge(d0all[,c(conVar2, 'Hugo_Symbol')], d0[,c(conVar1,conVar2)])
    gfreq1<-inner_join(d0all[,c(conVar2, 'Hugo_Symbol')], 
                       d0[,c(conVar1,conVar2)], by=conVar2)
    gfreq1<-gfreq1[,c(conVar1,'Hugo_Symbol')]
    gfreq1<-distinct(gfreq1)
    gfreq2<-table(gfreq1$Hugo_Symbo)/N
    gfreq2<-sort(gfreq2[gfreq2>mutFreqTh], decreasing=T, method="quick")
    if(length(gfreq2)==0){gfreq2<-NULL}
    gHR2<-NULL
    if(sum(d0$OS_mon, na.rm=T)>10){
      # ghr<-ghrp<-NULL
      # for(x in names(gfreq2[1:max(30,length(gfreq2))])){
      #   s1<-d0all[d0all[,'Hugo_Symbol']==x, conVar2]
      #   bm1<-rep(0, nrow(d0))
      #   bm1[d0[,conVar2]%in%s1]<-1
      #   t1<-try(f1<-coxph(Surv(OS_mon, OS_cnr)~bm1, data=d0))
      #   if(class(t1)=='try-error' | any(is.na(f1$coefficients))){
      #     ghr<-c(ghr,NA)
      #     ghrp<-c(ghrp,NA)
      #   }else{
      #     sc.f1<-summary(f1)$coefficients
      #     ghr<-c(ghr, round(sc.f1['bm1','exp(coef)'],2))
      #     ghrp<-c(ghrp,round(sc.f1['bm1','Pr(>|z|)'],2))
      #   }
      # } #needs to improve with fast fit
    }
    #get significant correlations with interested biomarker
    #get non-NA Hazard ratios with p-values
    
    #components for KM plot
    leg2<-paste0('HR=',HR, ' p=',pval)
    forplot<-list(f0=f0,leg1=leg1,leg2=leg2, col1=col1,lty1=lty1, tt1=v1,
                  BS_size_number=bs1, pctv=pctv)
    return(list(otb=otb, forplot=forplot, bmk1freq=gfreq2, bmk1HR=gHR2))
  }
  
  #function for get correlated genes
  #output a one-row matrix for multiple genes
  scorT<<-function(yy=distinct(sdat2.1[,c('submitted_sample_id','testBmk')]),
        xx=sMut2t[,c('submitted_sample_id','labelBmk')], 
        th=0.05){
    xx<-xx[xx$submitted_sample_id%in%yy$submitted_sample_id,]
    if(nrow(xx)<11){return(NULL)}
    xx2<-as.data.table(xx); xx2$dum<-1
    xx2w<-dcast(xx2, submitted_sample_id~labelBmk, value.var='dum')
    yy2<-yy[yy$submitted_sample_id%in%xx2w$submitted_sample_id,]
    o1<-order(yy2$submitted_sample_id)
    t1<-try(f1<-cv.glmnet(as.matrix(xx2w[,-1]), yy2[o1,2], 
                      family='binomial', alpha=.5))
    if(class(t1)=='try-error'){return(NULL)}
    cef1<-coef(f1, s=0)[-1,1]
    cef1<-cef1[abs(cef1)>th]
    if(length(cef1)>0){
      cef2<-data.frame(matrix(cef1,nrow=1));colnames(cef2)<-names(cef1);
      return(cef2)
    }else{return(NULL)}
  }
  
  #Function for group summary tables
  #global objects: sprjL$clin, sprjL$mut, conVar1 
  grpSum<<-function(
    geneText='rb1,mycn',  #for mutated gene names
    geneHS=NULL, #exact Hugo_Symbol
    n1=40, b1=500, #sample size and numbers of times for bootstrap
    sel_grpVar=input$dropdown, #select the group variable names (upto 3)
    sel_grpVar_lev=list(v1=gLv1(sel_grpVar), v2=NA,
                        v3=NA),#select levels under the grpVar 
    maxLev_grp=1000, #the maximum number of groups or sub groups
    sel_sgrpVar=input$dropdown4, #var names for biomarker selection
    sel_sgrpVar_lev=gLv2(sel_sgrpVar,c('RB1','MYCN')), #levs in bmk
    min_n_sgrp=5, #the maximum number of subgroups or sub groups
    runCorMut=FALSE, #whether run corMut analysis
    mutFreqTh1=0
  ){
    #just take upto 3 variables for grouping, upto 1 variables
    len.sel_grpVar<-length(sel_grpVar)
    if(len.sel_grpVar>3){sel_grpVar<-sel_grpVar[1:3]}
    if(length(sel_sgrpVar)>1){sel_sgrpVar<-sel_sgrpVar[1]}
    if(!is.list(sel_grpVar_lev)){
      sel_grpVar_lev<-list(v1=sel_grpVar_lev,v2=NA,v3=NA)
    }else{
      l0<-list(v1=NA,v2=NA,v3=NA)
      if(length(sel_grpVar)>0)l0$v1<-sel_grpVar_lev[[1]]
      if(length(sel_grpVar)>1)l0$v2<-sel_grpVar_lev[[2]]
      if(length(sel_grpVar)>2)l0$v3<-sel_grpVar_lev[[3]]
      sel_grpVar_lev<-l0
    }
    
    o1<-'fail to find any matched mutated genes for '
    o1<-data.frame(id=1, x1=paste0(o1,geneText))
    
    #find matched mutated genes
    print(geneHS)
    geneHS<-gsub(' ', '', geneHS, fixed=T)
    geneHS<-geneHS[geneHS!='']
    if(!is.null(geneHS) && length(geneHS)==0){
      return(list(stb=o1,kmplot=NULL, cor_mut=NULL, bmk1freq=NULL))
    }
    if(!is.null(geneHS)){
      gene1 <- geneHS
    }else{
      if(!is.null(geneText)){
        gene1<-findHugo(geneText, uni.mutg)[[2]]
      }else{
        gene1<-NULL
      }
    }
    if(length(gene1)==0){
      return(list(stb=o1,kmplot=NULL, cor_mut=NULL, bmk1freq=NULL))
    }
    #confirm the existence of the group variables
    if(any(!(has.sel_grpVar<-sel_grpVar%in%grpVar.clin))){
      return(list(stb=data.frame(id=1,
        x1=paste("failed to find group variables:",
                  paste(sel_grpVar[has.sel_grpVar],collapse=', '))),
        kmplot=NULL, cor_mut=NULL, bmk1freq=NULL))
    }
    #confirm the existence of sub group variable for biomarkers
    if(!sel_sgrpVar%in%sgrpVar){
      return( list(stb=data.frame(id=1,
        x1=paste('failed to find sub group variables:', sel_sgrpVar)), 
        kmplot=NULL, cor_mut=NULL, bmk1freq=NULL))
    }
    #~~~~stat to use sprjL
    #setup the mutation data
    sMut1<-sprjL$mut; #s1.mut<-unique(sMut1[,conVar2]);#get sample ids
    #setup the clinical data
    sdat1<-sprjL$clin; #sdat1<-sdat1[sdat1[,conVar2]%in%s1.mut,];
    #correct typo in colnames
    colnames(sdat1)[colnames(sdat1)=='cna']<-'cnv'
    if(!'cnv'%in%colnames(sdat1)){sdat1$cnv<-NA}
    if(nrow(sdat1)<3){
      return( list(stb=data.frame(id=1,
        x1=paste('failed to find 3+ clin/mut matched samples.')), 
        kmplot=NULL, cor_mut=NULL, bmk1freq=NULL))
    }
    #remove NA values in survival_time
    if(!"survival_time"%in%colnames(sprjL$clin)){
      sdat1$survival_time<-NA
      sdat1$vital_status<-NA
    }
    sdat1$survival_time<-as.numeric(as.character(sdat1$survival_time))
    sdat1$vital_status<-as.character(sdat1$vital_status)
    #pull out and clean up the levels in group variables
    for(m in 1:length(sel_grpVar)){
      u.lev<-as.character(unique(sdat1[,sel_grpVar[m]])) 
      if(length(u.lev)==0) next
      #limit number of group levels to show
      tmp<-sel_grpVar_lev[[m]][sel_grpVar_lev[[m]]%in%u.lev]
      if(length(tmp)>maxLev_grp){sel_grpVar_lev[[m]]<-tmp[1:maxLev_grp]}
      #update the subset of clinical data with group variables
      sdat1<-sdat1[sdat1[,sel_grpVar[m]]%in%sel_grpVar_lev[[m]],]
    }
    #update patient&sample ids after group selection
    sdat1[,conVar1]<-as.character(sdat1[,conVar1])
    #shrink down if only CNV marker
    if(all(grepl('CNV:', sel_sgrpVar_lev))){
      sdat1<-sdat1[!is.na(sdat1$cnv),]
    }else if(all(grepl('FUS:', sel_sgrpVar_lev))){
      #shrink down if only fusion marker
      sdat1<-sdat1[!is.na(sdat1$fus),]
    }else if(sum(grepl('CNV:|FUS:', sel_sgrpVar_lev))==0){
     #shrink down if only somatic mutation samples are available
      sdat1<-sdat1[!is.na(sdat1$mut),]
    }
    uni.ids<-unique(sdat1[,conVar1])
    #pull data for mutant frequency estimate across genes
    sa0<-sMut1[,sel_sgrpVar[1]]%in%sel_sgrpVar_lev
    sMut1.sa<-sMut1[sa0,c(conVar2, 'Hugo_Symbol', sel_sgrpVar[1])]
    #sMut1.sa$labelBmk<-paste(sMut1.sa$Hugo_Symbol, sMut1.sa[,sel_sgrpVar[1]])
    #let the mutation data have only the selected genes
    s1<- sMut1$Hugo_Symbol%in%gene1 
    s11<-sa0 & s1
    if(runCorMut){
      sMut1t<-sMut1
      sMut1t$labelBmk<-paste(sMut1t$Hugo_Symbol, sMut1t[,sel_sgrpVar[1]])
      sMut2t<-sMut1t[!s11,]#pull out data with other markers
      sMut1t$labelBmk[!s11]<-rep('', sum(!s11))
      #set numeric 0,1 value for biomarker positive samples
      sMut1t$testBmk<-0; sMut1t$testBmk[s11]<- 1;
      sMut1t<-as.data.table(sMut1t)
      sMut1s<-sMut1t[,1*(sum(testBmk)>0), by=submitted_sample_id]
      colnames(sMut1s)[2]<-'testBmk'
      sMut1s1<-sMut1t[labelBmk!='',paste(unique(labelBmk),collapse='|'), 
                      by=submitted_sample_id]
      colnames(sMut1s)[2]<-'labelBmk'
      sMut1s<-inner_join(sMut1s, sMut1s1);
    }else{
      sMut1t<-sMut1[s11,]
      sMut1t$labelBmk<-paste(sMut1t$Hugo_Symbol, sMut1t[,sel_sgrpVar])
      sMut1t<-as.data.table(sMut1t)
      sMut1s<-sMut1t[,paste(unique(labelBmk),collapse='|'), 
                     by=submitted_sample_id]
      colnames(sMut1s)[2]<-'labelBmk'
      sMut1s$testBmk<-1
    }
    # sdat1.1<-merge(sdat1, sMut1s, by=conVar2, all.x=T)
    sdat1.1<-left_join(sdat1, sMut1s, by=conVar2)
    if(nrow(sdat1.1)==0){
      return(list(stb=o1,kmplot=NULL, cor_mut=NULL, bmk1freq=NULL))
    }
    #derive tte in months for overal survival
    sdat1.1$OS_mon<-sdat1.1$survival_time/30
    sdat1.1$OS_cnr<-NA #censor variable: 0=alive, 1=deceased
    sdat1.1$OS_cnr[sdat1.1$vital_status=='alive']<-0
    sdat1.1$OS_cnr[sdat1.1$vital_status=='deceased']<-1
    #colnames in the data ready for sub-group summary
    sdat2.c <-c('OS_mon','OS_cnr','testBmk', 'labelBmk',sel_grpVar)
    #get the data ready for group summary
    col1.1<-c(conVar1, conVar2, sdat2.c, 'mut','cnv')
    #col1.1<-col1.1[!col1.1%in%colnames(sdat1.1)]
    sdat2<-distinct(sdat1.1[,col1.1])
    s22<-is.na(sdat2$testBmk)
    if(length(s22))
    sdat2$testBmk[s22]<-0
    sdat2$labelBmk[s22]<-''
    
    #construct unique levels across group variables
    grpVar_l2m<-lev2mat(sel_grpVar_lev)
    
    #~~start to get summary table within each group~~#
    oo<-pplot<-corM<-bmk1freq<-list() #the output object of summary table
    pplot.id<-1;#the object for KM plots
    for(m1 in 1:nrow(grpVar_l2m)){
      print(m1) #here is very slow
      sdat2.1<-sdat2[!(is.na(sdat2$mut)&is.na(sdat2$cnv)),]
      for(m2 in 1:ncol(grpVar_l2m)){ #only keep obs in sel_grpVar
        sdat2.1s<-sdat2.1[,sel_grpVar[m2]]%in%grpVar_l2m[m1,m2]
        sdat2.1<-sdat2.1[sdat2.1s, ]
      }
      n.sdat2.1<-nrow(sdat2.1)
      if(n.sdat2.1<5){oo[[m1]]<-NULL;next}
      # sdat2.1$bmk<-ifelse(sdat2.1[,sel_sgrpVar]%in%sel_sgrpVar_lev,1,0)
      tm.bmk<-table(sdat2.1$testBmk)
      if(length(tm.bmk[names(tm.bmk)==1])==0){next}
      oo1<-sgrpBSsum(d0=sdat2.1, var1=sel_grpVar, d0all=sMut1.sa, #total gene
                     size1=n1, m1=b1, mutFreqTh=mutFreqTh1)#get bootstrap result
      if(is.null(oo1)){next}
      #correlational study with other mutations
      if(runCorMut){
        scort1<-scorT(yy=distinct(sdat2.1[,c('submitted_sample_id','testBmk')]),
              xx=sMut2t[,c("submitted_sample_id",'labelBmk')], th=0.05)
        if(!is.null(scort1)){scort1<-data.frame(id=pplot.id,scort1)}
        corM[[pplot.id]]<-scort1 #1 row for 1 group
      }
      oo[[pplot.id]]<-oo1$otb
      pplot[[pplot.id]]<-oo1$forplot
      bmk1freq[[pplot.id]]<-oo1$bmk1freq
      pplot.id<-pplot.id+1
    }
    oo<-do.call('rbind',oo)
    colnames(oo)[1:(len.sel_grpVar+1)]<-c(sel_grpVar, sel_sgrpVar)
    oo<-data.frame(id=1:nrow(oo), oo)
    if(runCorMut){corMm<-rbindlist(corM, fill=TRUE)}else{corMm<-NULL}
    return(list(stb=oo, kmplot=pplot, cor_mut=corMm, bmk1freq=bmk1freq))
  }
  
  #draw KM plot sgrp_sum
  sgrp_km_plot<<-function(refL){ #refL=forplot
    #refL is a list entry from a list of KM reference collections
    col1<-refL$col1
    lty1<-refL$lty1
    tt1<-refL$tt1
    f0<-refL$f0
    leg1<-refL$leg1
    leg2<-refL$leg2
    if(is.null(col1)){
      plot(0~0,col='white',ylab='',xlab='',axes=F)
      mtext('no data')
    }else{
      plot(f0, xlab=paste0('months'), 
           ylab='Probability of Survival',
           main=tt1,
           col=col1, lty=lty1)
      legend('topright',legend=leg1,col=col1, lty=lty1, bty='n')
      legend('bottomleft',legend=leg2, bty='n')
    }
  }
  
  #draw hist for bootstrap prevalence
  sgrp_hist_bspct<<-function(refL){
    hist(refL$pctv,
         xlab=paste0('Bootstrap estimated prevalence','\n',
                     refL$BS_size_number),
         main=refL$tt1)
  }
  
  #get BS size and number
  gBS1<<-function(x='40,500'){
    if(is.null(x)){x='40, 100'}
    x1<-as.numeric(strsplit(x,split=',')[[1]])
    x1<-x1[!is.na(x1)]
    if(length(x1)==0){
      return(c(40,1000))
    }else if(length(x1)==1){
      return(c(x1,1000))
    }else{
      return(x1[1:2])
    }
  }
  
  #for BS, KM, bar, heatmap
  sgrpPlot<<-function(ids=input$outTbl_rows_selected,
                      freqTH=as.numeric(input$text3)){
    ids1<-ids[length(ids)]
    if(is.null(freqTH) || is.na(freqTH)){freqTH<-0.1}
    if(freqTH<=0 | freqTH>=1){freqTH<-0.1}
    forBarpFreq<-km_out$bmk1freq[[ids1]][km_out$bmk1freq[[ids1]]>freqTH]
    if(length(forBarpFreq)>0){
      nf<-layout(t(matrix(c(1,2,3,3),nrow=2)),widths=c(1,1),
                 heights=c(2,2))
      # layout.show(nf)
      par(mai=c(0.5,0.5,0.5,0.5))
      sgrp_hist_bspct(km_out$kmplot[[ids1]])
      sgrp_km_plot(km_out$kmplot[[ids1]])
      par(mai=c(1.3,0.7,0.5,0))
      barplot(forBarpFreq,las=2)
    }else{
      par(mfrow=c(1,2)) 
      sgrp_hist_bspct(km_out$kmplot[[ids1]])
      sgrp_km_plot(km_out$kmplot[[ids1]])
    }
  }
  
}


#pull out all variables as a reference for new database add-in
if(F){
  tt1<-Sys.time();try(sprjL<<-select.prj('all'));Sys.time()-tt1
  allVars<-rbind(
    data.frame(X1=colnames(sprjL$clin), X2=rep('clin',ncol(sprjL$clin))),
    data.frame(X1=colnames(sprjL$mut), X2=rep('mut', ncol(sprjL$mut))) )
  write.csv(allVars, file=file.path(frompath,'LightON_variables.csv'),
            row.names=F)
}


#test codes to check bugs
if(F){
  source(file.path(projP, progfile))
  
  input<-NULL
  #1st module
  #input$text<- '85, 11,10,9,8,7' #'all'
  tx1<-'85, 11,10,9,8,7'
  tt1<-Sys.time();try(sprjL<<-select.prj(tx1));Sys.time()-tt1
  # see.click.prj(c(78), showMutOnly = F)
  
  #2nd module
  input$text<-'RB1,ARID1A,KMT2C,KMT2D, MYCN, AURKA, AURKB'
  # paste(findHugo(input$text)[[1]])
  
  input$text2<-'40,1000' #parameters for bootstrap
  input$dropdown<-'proj_id'
  input$dropdown2<-gLv1(input$dropdown)
  input$dropdown3<-"mutation_type"; #sgrpVar
  input$dropdown7<-'RB1'  #c('RB1','AURKA','KMT2C','KMT2D','ARID1A')#, 'MYCN')
#  input$dropdown7<-c('NTRK1','NTRK2', 'NTRK3')
  input$dropdown4<-gLv2(input$dropdown3, input$dropdown7)
  
  n1b1<<-gBS1(input$text2)
  ttt1<-Sys.time()
  try(km_out<<-grpSum(geneText=NULL, 
                      geneHS=input$dropdown7,
                      n1=n1b1[1], b1=n1b1[2],
                      sel_grpVar=input$dropdown, 
                      sel_grpVar_lev=list(v1=input$dropdown2,
                                          v2=input$dropdown5,
                                          v3=input$dropdown6),
                      sel_sgrpVar=input$dropdown3, #varnames for sel bmk
                      sel_sgrpVar_lev=input$dropdown4,
                      mutFreqTh1=0))
  Sys.time()-ttt1
  print( km_out$stb )
  print( km_out$bmk1freq )
  par(mfrow=c(1,2))
  sgrp_hist_bspct(km_out$kmplot[[1]])
  sgrp_km_plot(km_out$kmplot[[1]])
  sgrpPlot(ids=1)
  
  if(F){
    geneText=NULL; 
    geneHS=input$dropdown7;
    n1=n1b1[1]; b1=n1b1[2];
    sel_grpVar=input$dropdown; 
    sel_grpVar_lev=list(v1=input$dropdown2,
                        v2=input$dropdown5,
                        v3=input$dropdown6);
    sel_sgrpVar=input$dropdown3; #varnames for sel bmk
    sel_sgrpVar_lev=input$dropdown4;
    maxLev_grp=500;
    min_n_sgrp=5;
    runCorMut=FALSE;
    mutFreqTh1=0.02
  }
  
}

if(F){ #for KRAS G12C
  source(file.path(projP, progfile))
  input<-NULL
  
  #1st module
  tx1<-'all'
  tt1<-Sys.time();try(sprjL<<-select.prj(tx1));Sys.time()-tt1
  # lapply(sprjL, dim)
  # head(sprjL$mut)
  
  #2nd module
  input$text<-'KRAS'
  print(paste(findHugo(input$text)[[1]]))
  input$text2<-'40,1000' #parameters for bootstrap
  input$dropdown<-c('db','proj_id','proj_name')
  input$dropdown2<-gLv1('db')
  input$dropdown5<-gLv1('proj_id')
  input$dropdown6<-gLv1('proj_name')
  input$dropdown3<-"aa_mutation"; #sgrpVar
  input$dropdown7<-'KRAS'
  mtype<-gLv2(input$dropdown3, input$dropdown7)
  print(input$dropdown4<-mtype[grep('g12c',tolower(mtype))])
  
  n1b1<<-gBS1(input$text2)
  ttt1<-Sys.time()
  try(km_out<<-grpSum(geneText=NULL, 
                      geneHS=input$dropdown7,
                      n1=n1b1[1], b1=n1b1[2],
                      sel_grpVar=input$dropdown, 
                      sel_grpVar_lev=list(v1=input$dropdown2,
                                          v2=input$dropdown5,
                                          v3=input$dropdown6),
                      sel_sgrpVar=input$dropdown3, #varnames for sel bmk
                      sel_sgrpVar_lev=input$dropdown4,
                      mutFreqTh1=0))
  Sys.time()-ttt1
  head( outtb1<-km_out$stb )
  #colnames(outtb1)[grepl('95CI', colnames(outtb1))]<-c('N','prevl_BMKpos_95CI_n40')

}


#set default project 
try(sprjL<<-select.prj("1"))

#function for update covariate names
grpVar.clin.fun<<-function(
  dum=NULL
){
  print(dum)
  return(grpVar.clin0[grpVar.clin0%in%colnames(sprjL$clin)])
}






