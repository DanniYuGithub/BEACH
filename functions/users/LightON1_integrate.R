if(TRUE){#header
  #*****************************************************************************
  #Eli Lilly and Company - GLOBAL STATISTICAL SCIENCES - PROGRAM
  #CODE NAME (required)                : LightON1_integrate.R
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
  #---- ------------       ---------------------------------------  ---------
  #1.1  Danni Yu            Code Creation                             01/10/2019
  #1.2  Danni Yu            Add ICGC projects                         01/29/2019
  #1.3  Danni Yu            Revision and debug                        03/12/2019
  #*****************************************************************************
}

#manage folders by project index 1, 2, .... in LightON1 folder
#Each project folder includes:
# X_datasum.feather for project summary
# X_clin.feather for patients and samples information (all with mut test)
# X_mut.feather for somatic mutation data
# X_gexp_array.feather, X_gexp_rnaseq.feather  for gene expression data
# X_pexp_antib.feather for protein expression named by antibody (optional)
# X_mut_cna.feather for copy numbers (optional)
# X_fus.feather for fusion (optional)
# where 'X' project names; 'xx' gene names

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#!!!!!!! survival_time in days !!!!!!!!!!#
# c('alive','deceased') in vital_status  #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

#define paths, files, program,...
if(TRUE){
  rm(list=ls());
  
  beginTime<-Sys.time()
  
  runDB1<-T #convert the 1st database ICGC to innitiate the project summary file
  runDB2<-FALSE #convert IMPACT database
  runDB3<-TRUE #convert GENIE database
  
  rerunClin<-F #enforce to update clin data
  rerunMut<-F #enforce to update mutation data
  rerunExp<-F #enforce to update gene expression data
  rerunPExp<-F #enforce to update protein expression data
  
  #define comment key words in output files
  com.clin<-'_clin.feather'
  com.gexp1<-'_gexp_array.feather'
  com.gexp1m<-'_mut_gexp_array.feather'
  com.gexp2<-'_gexp_rnaseq.feather'
  com.gexp2m <- '_mut_gexp_rnaseq.feather'
  com.mut1<-'_mut.feather'
  com.mut2<-'_mut_cnv.feather'
  com.cnv<-'_cnv.feather'
  com.pexp<-'_pexp_antib.feather'
  com.pexpm<-'_mut_pexp_antib.feather'
  com.sum<-'_datasum.feather'
  
  #global mapping variable
  conVar1 <- "submitted_id"
  #sample level mapping variable
  conVar2 <-"submitted_sample_id"
  
  options(stringsAsFactors=F)
  
  projP <- "local1" #to be defined by users with a local drive
  
  dir(progpath<-file.path(projP,'prog'))
  progfile<-'LightON1_integrate.R'
  
  dir(frompath<-"data")
  dir(topath<-'local2')#to be defined by users with a local drive
  projfile<-'projects_summary.csv'
  dsfile<-'dataspec.csv'
  gnameFile<-file.path(topath, 'gene_ensemble.Rdata')
  
  #icgc raw data
  db1<-'ICGC'
  dir(path1<-file.path(frompath, db1, 'release_27/Projects'))
  path1file<-c('icgc_fnm.Rdata',"icgc_project_code.csv",
               'icgc_protein_antibody.csv', 'illumina_probeID.csv')
  
  #MSK-impact raw data
  db2<-'MSK-IMPACT'
  path2<-file.path(topath, db2)
  
  #GINIE raw data
  db3<-'GENIE'
  path3<-file.path(topath, db3)
  path3f<-file.path(frompath, 'GENIE/data')
  
  #for a quick summary about platform
  cnv.plat<-mut.plat<-fus.plat<-NULL
  
  #grep all existing file names
  if(TRUE){
    f1s<-dir(file.path(topath,'ICGC')); f1nms<-list();
    for(i in 1:length(f1s)){
      f1nms[[i]]<-dir(file.path(topath,'ICGC', f1s[i]))
    }
    names(f1nms)<-f1s
    ##check file names for protein expression
    # unlist(sapply(f1nms, function(x){x[grep('_pexp',x)]}))
  }
  
  
}

# predefine functions
if(TRUE){
  library(data.table)
  library(feather)
  library(reshape)
  library(dplyr)
  
  #get or update datasum file
  getDSum<-function(din, kvar=NULL, v2=NULL, ds=datasum){
    if(is.null(din)||length(din)==0){return(ds)}
    if(is.null(kvar)){
      if(is.null(v2)){stop('define v2 in getDSum')}
      o<-data.frame(x1=unique(din), x2=1)
      colnames(o)<-v2
    }else{
      tmpv<-din[,!colnames(din)%in%kvar]
      tmpv[!is.na(tmpv)]<-1
      tmpv[is.na(tmpv)]<-0
      o<-data.frame(din[,kvar], tmpv)
      colnames(o)[1:length(kvar)]<-kvar
    }
    if(!is.null(ds)){
      ds<-merge(ds,o,all=T)
      ds[is.na(ds)]<-0
    }else{ds<-o}
    return(ds)
  }
  
  
  #update data spec file
  up.dataspec<-function(d1=dataspec, #
                        d2, #new data
                        v1){
    if(is.null(d2))return(d1)
    ds1<-data.frame(DataName=rep(v1,ncol(d2)),VarName=colnames(d2))
    d1<- unique(rbind(d1,ds1))
    return(d1)
  }
  
  #control gene names as Hugo_Symbol
  gnmIcgc<-function(d1=mut1, #or mut2 for expa
                    v1='gene_affected', #'gene_id'
                    gl=g_list12 ){
    d1v<-d1[,v1]
    d1v<-as.vector(as.character(d1v))
    d1v<-gsub(' ', '', d1v)
    d1<-d1[!is.na(d1v)&d1v!='',]
    d1v<-d1[,v1]
    d1v<-gsub('\\..*','', d1v)
    d1[,v1]<-d1v
    if(any(unique(d1v)%in%gl$refseq$Hugo_Symbol)){
      d1$Hugo_Symbol<-d1v
    }else if(any(unique(substring(d1v,1,3))%in%'NM_')){
      d1<-merge(d1, gl$refseq, by.x=v1, by.y='refseq')
    }else if(substring(d1v[1],1,4)=='ENSG'){
      d1<-merge(d1, gl$ensembl, by.x=v1, by.y='ensembl')
    }else if(any(unique(d1v)%in%gl$illprob$ILL_prob)){
      d1<-merge(d1, gl$illprob, by.x=v1, by.y='ILL_prob')
    }else if(any(unique(d1v)%in%gl$ensemblTr$ensemblTr)){
      d1<-merge(d1, gl$ensemblTr, by.x=v1, by.y='ensemblTr')
    }else{}
    return(d1[,!colnames(d1)%in%v1])
  }
  
  #calculate z score
  getZ<-function(d1=gexp, v1='submitted_sample_id', v2='Hugo_Symbol',  
                 v3='gexpZ'){
    d2<-d1[,c(v1, v2,v3)]
    for(q in unique(d2[,v1])){
      s<-d2[,v1]==q
      d2[s,v3] <- (d2[s,v3]-mean(d2[s,v3],na.rm=T))/sd(d2[s,v3],na.rm=T)
    }
    return(d2)
  }
  
  #map from chromosomal location to Hugo Symbol
  #exaustive mapping method: as long as having an overlap with the genes
  chrm2hs<-function(d1, vs=NULL, v2='mutation_type', vid='submitted_sample_id',
                    g1=g_list3){
    #internal function to filter out problomatic segments
    filt1<-function(x){
      x1s<-g2[x,]
      dis1<-max(x1s$end_position) - min(x1s$start_position)
      dis2<-max(x1s$end_position-x1s$start_position)
      x1<-x1s$Hugo_Symbol 
      #remove the pieces too long to be true
      if(dis1<=5*dis2){return(x1)}else{return(NULL)}
    }
    #d1=mut2; 
    vid<-vid[1]
    if(is.null(vs)){
      vs<-c('chromosome','chromosome_start','chromosome_end')
    }
    if(any(!vs%in%colnames(d1))){stop('chrm2hs, not matched vs var names')}
    #get LOH
    loh<-unique(d1[grepl('loh', tolower(d1[,v2])),c(vid,v2,vs)])
    loh3<-NULL
    if(nrow(loh)>0){
      loh$lab<-gsub(' ','_',loh[,v2])
      loh$Hugo_Symbol<-NA
      c1<-c2<-c3<-c4<-NULL
      for(chr1 in unique(loh[,vs[1]])){
        g2<-g1[g1$chromosome_name==chr1,]
        loh1<-loh[loh[,vs[1]]==chr1,]
        n1<-nrow(loh1); n2<-nrow(g2)
        if(any(n1==0|n2==0)) next
        mat1s<-matrix(rep(loh1[,vs[2]],n2), ncol=n2)
        mat1e<-matrix(rep(loh1[,vs[3]],n2), ncol=n2)
        mat2s<-t(matrix(rep(g2$start_position,n1), ncol=n1))
        mat2e<-t(matrix(rep(g2$end_position,n1), ncol=n1))
        s1<-mat1s<=mat2s&mat1e>=mat2s | mat1s<=mat2e&mat1e>=mat2e
        s2<-apply(s1, 1, filt1)
        if(is.null(s2) | length(s2)==0) next
        ls2<-sapply(s2, length);
        MM<-1:length(ls2);
        s2.ul<-unlist(s2)
        c1<-c(c1,unlist(lapply(MM,function(x){rep(loh1[x,vid],ls2[[x]])})))
        c2<-c(c2,s2.ul)
        c3<-c(c3,unlist(lapply(MM,function(x){rep(loh1[x,'lab'],ls2[[x]])})))
        c4<-c(c4,rep(chr1,length(s2.ul)))
      }
      d2<-data.frame(c1, c2,c3, c4)
      if(nrow(d2)>0){
        colnames(d2)<-c(vid,'Hugo_Symbol','cnv','chromosome')
        d2$mutation_type<-paste0('CNV:',d2$cnv)
        loh3<-d2
      }
    }
    
    #get gain or loss result
    d1[,v2]<-tolower(as.character(d1[,v2]))
    d1<-unique(d1[d1[,v2]%in%c('gain','loss'),c(vid,v2,vs)])
    tb1<-table(d1$mutation_type)
    print(tb1)
    if(all(is.na(tb1[c('gain','loss')]))) return(NULL)
    d1$dum<-ifelse(d1[,v2]=='gain', 1, -1) 
    d1$Hugo_Symbol<-NA
    c1<-c2<-c3<-c4<-d3<-NULL
    for(chr1 in unique(d1[,vs[1]])){#procece chrom by chrom
      g2<-g1[g1$chromosome_name==chr1,]
      dd1<-d1[d1[,vs[1]]==chr1,]
      n1<-nrow(dd1); n2<-nrow(g2)
      if(any(n1==0|n2==0)) next
      mat1s<-matrix(rep(dd1[,vs[2]],n2), ncol=n2)
      mat1e<-matrix(rep(dd1[,vs[3]],n2), ncol=n2)
      mat2s<-t(matrix(rep(g2$start_position,n1), ncol=n1))
      mat2e<-t(matrix(rep(g2$end_position,n1), ncol=n1))
      s1<-mat1s<=mat2s&mat1e>=mat2s | mat1s<=mat2e&mat1e>=mat2e
      s2<-apply(s1, 1, filt1)
      if(is.null(s2) | length(s2)==0) next
      ls2<-sapply(s2, length); 
      MM<-1:length(ls2);
      s2.ul<-unlist(s2)
      c1<-c(c1,unlist(lapply(MM,function(x){rep(dd1[x,vid],ls2[[x]])})))
      c2<-c(c2,s2.ul)
      c3<-c(c3,unlist(lapply(MM,function(x){rep(dd1[x,'dum'],ls2[[x]])})))
      c4<-c(c4,rep(chr1,length(s2.ul)))
    }
    d2<-data.frame(c1, c2,c3, c4)
    if(nrow(d2)>0){
      colnames(d2)<-c(vid,'Hugo_Symbol','cnv','chromosome')
      d3<-aggregate(cnv~., data=d2,FUN=sum)
      d3$mutation_type<-'CNV:eql'
      d3$mutation_type[d3$cnv>0]<-'CNV:amp'
      d3$mutation_type[d3$cnv<0]<-'CNV:del'
    }
    d4<-rbindlist(list(d3,loh3))
    return(d4)
  }
  
  #function convert numeric age to categories
  age2c<-function(ag1=clin$donor_age_at_enrollment){
    ag1C<-rep(NA,length(ag1))
    ag1C[ag1<=1]<-'0~1'
    ag1C[ag1>1 & ag1<=5]<-'2~5'
    ag1C[ag1>5 & ag1<=10]<-'6~10'
    ag1C[ag1>10 & ag1<=15]<-'11~15'
    ag1C[ag1>15 & ag1<=20]<-'16~20'
    ag1C[ag1>20 & ag1<=25]<-'21~25'
    ag1C[ag1>25 & ag1<=30]<-'26~30'
    ag1C[ag1>30 & ag1<=35]<-'31~35'
    ag1C[ag1>35 & ag1<=40]<-'36~40'
    ag1C[ag1>40 & ag1<=45]<-'41~45'
    ag1C[ag1>45 & ag1<=50]<-'46~50'
    ag1C[ag1>50 & ag1<=55]<-'51~55'
    ag1C[ag1>55 & ag1<=60]<-'56~60'
    ag1C[ag1>60 & ag1<=65]<-'61~65'
    ag1C[ag1>65 & ag1<=70]<-'66~70'
    ag1C[ag1>70 & ag1<=75]<-'71~75'
    ag1C[ag1>75 & ag1<=80]<-'76~80'
    ag1C[ag1>80 & ag1<=85]<-'81~85'
    ag1C[ag1>85 & ag1<=90]<-'86~90'
    ag1C[ag1>95]<-'96+'
    return(ag1C)
  }
  
  
  #---get ranges for age, os---#
  ifunc1<-function(x){
    x<-x[is.finite(x)&!is.na(x)]
    if(length(x)==0) return(NA) else return(range(x))
  }
  
}

#get mappings between genes symbols and ensembles
if(FALSE){
  load(gnameFile)#for gene.ensemble
  library('biomaRt')
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  # listAttributes(mart)
  G_list<-getBM(attributes=c("ensembl_gene_id", 'refseq_mrna', "hgnc_symbol",
                             'ensembl_transcript_id', 'chromosome_name',
                             'start_position', 'end_position'),
                filters ="ensembl_gene_id",
                values  =gene.ensemble, 
                mart    =mart)
  colnames(G_list)[1:4] <- c('ensembl','refseq', 'Hugo_Symbol', 'ensemblTr')
  G_list<-G_list[G_list$Hugo_Symbol!='',]
  g_list<-unique(G_list[,c('ensembl','Hugo_Symbol')])
  g_list2<-unique(G_list[G_list$refseq!='',c('refseq','Hugo_Symbol')])
  g_list5<-unique(G_list[G_list$ensemblTr!='',c('ensemblTr','Hugo_Symbol')])
  
  #for mapping between HG and chromosome location
  g_list3<-unique(G_list[,c('Hugo_Symbol','chromosome_name','start_position',
                            'end_position', 'band')])
  g_list3<-g_list3[!is.na(g_list3$Hugo_Symbol),]
  g_list3<-g_list3[g_list3$Hugo_Symbol!='',]
  st1<-aggregate(start_position~Hugo_Symbol+chromosome_name, 
                 data=g_list3, FUN=min)
  ed1<-aggregate(end_position~Hugo_Symbol+chromosome_name, 
                 data=g_list3, FUN=max)
  g_list3<-merge(st1,ed1)

  #get illumina probe ID
  g_list4<-read.csv(file.path(path1, path1file[4]), h=T, comment.char='#')
  g_list4<-unique(g_list4[,c('PROBE_ID','GENE')])
  colnames(g_list4)<-c('ILL_prob','Hugo_Symbol')
  
  #save the gene mapping file
  g_list12<-list(ensembl=g_list, refseq=g_list2, illprob=g_list4, ensemblTr=g_list5)
  save(gene.ensemble, g_list, g_list2, g_list12, g_list3, file=gnameFile)

  

}else{load(gnameFile)} #for the gene mapping: g_list

#build up the innitial data summary file
prj.sum0<-data.frame(id=NA, db=NA, db_id=NA, 
                     proj_id=NA, proj_name=NA, country=NA,
                     age_range=NA, #at enrollment age
                     OS_range=NA,  #survival_time
                     N=NA, #total number of subjects
                     N_mut=NA, N_mut_exp_array=NA, N_mut_exp_seq=NA, # Num pts
                     m_mutgene=NA, #number of unique genes
                     N_bmk=NA, #number of patients in bmk data
                     m_bmk=NA  #number of unique biomarkers
                     )

if(runDB1){
  prj.sum<-NULL
  dataspec<-data.frame(DataName=NULL, VarName=NULL)
}else{
  prj.sum<-read.csv(file.path(topath, projfile), h=T)
}


#---ICGC---#
#https://docs.icgc.org/submission/guide/donor-clinical-data-guidelines/
if(runDB1){
  icgc.stid<-0 #start id
  load(file.path(path1, path1file[1]))#fnm will be loaded here.
  prj_code<-read.csv(file.path(path1,path1file[2]), h=T)
  pexp_antb<-read.csv(file.path(path1,path1file[3]),h=T, comment.char='#')
  icgc.id<-1:length(fnm)
  if(!dir.exists(file.path(topath, db1))){dir.create(file.path(topath, db1))}
  #fnm1.cnt; fnm1.cnms; fnm1; fnm;
  prjid<-1
  for(i in icgc.id){
    #get the project name
    #i=78 #have donor therapy data
    pnm1<-names(fnm)[i]
    fin1<-file.path(path1,pnm1)#the input project folder
    ff1<-file.path(topath,db1,i)#output folder name
    print(paste(i, pnm1,'start.'))
    
    #define the output file names
    of.clin<-file.path(ff1, paste0(pnm1,com.clin))
    of.gexp1<-file.path(ff1, paste0(pnm1,com.gexp1)) #zscore only
    of.gexp1m<-file.path(ff1, paste0(pnm1,com.gexp1m)) 
    of.gexp2<-file.path(ff1, paste0(pnm1,com.gexp2)) #zscore only
    of.gexp2m<-file.path(ff1, paste0(pnm1,com.gexp2m))
    of.mut1<-file.path(ff1, paste0(pnm1,com.mut1))#including all the mutations
    of.mut2<-file.path(ff1, paste0(pnm1,com.mut2))
    of.cnv<-file.path(ff1, paste0(pnm1,com.cnv))
    of.pexp<-file.path(ff1, paste0(pnm1,com.pexp))
    of.pexpm<-file.path(ff1, paste0(pnm1,com.pexpm))
    of.sum<-file.path(ff1, paste0(pnm1,com.sum))#available for heatmap
    
    #define all the originial file names 
    clinf1<-file.path(fin1, paste0('donor.',pnm1,'.tsv'))
    clinf2<-file.path(fin1, paste0('donor_family.',pnm1,'.tsv'))
    clinf3<-file.path(fin1,paste0('donor_therapy.',pnm1,'.tsv'))
    clinf4<-file.path(fin1, paste0('sample.',pnm1,'.tsv'))
    clinf5<-file.path(fin1, paste0('specimen.',pnm1,'.tsv'))
    gexpf1<-file.path(fin1,paste0("exp_array.",pnm1, '.tsv'))
    gexpf2<-file.path(fin1,paste0("exp_seq.",pnm1, '.tsv'))
    pexpf1<-file.path(fin1,paste0("protein_expression.",pnm1, '.tsv'))
    mutf1<-file.path(fin1, paste0('simple_somatic_mutation.open.',pnm1,'.tsv'))
    mutf2<-file.path(fin1, paste0('copy_number_somatic_mutation.',pnm1,'.tsv'))
    
    #if no clin data, then skip the project
    if(!file.exists(clinf1)){next}
    if(!dir.exists(ff1)){dir.create(ff1)}

    #---setup original mutation files---#
    if(file.exists(of.mut1) & !rerunMut){
      t1<-try(mut<-as.data.frame(read_feather(of.mut1)))
      if(class(t1)=='try-error') file.remove(of.mut1)
    }else{mut<-NULL}
    subj1<-NULL
    if(file.exists(mutf1)){#simple somatic
      if(rerunMut|!file.exists(of.mut1)){
        mut1<-as.data.frame(fread(mutf1, h=T))
        mut1.nm<-c('submitted_sample_id', 'gene_affected', 'chromosome',
                   'mutation_type', 'aa_mutation','cds_mutation')
        mut<-gnmIcgc(unique(mut1[,mut1.nm]), 'gene_affected')
        write_feather(mut, of.mut1)#write down the mut data
      }
      u1s<-!grepl('CNV:|FUS:',mut$mutation_type)
      u1<-unique(mut$submitted_sample_id[u1s])
      if(length(u1)>0){
        subj1<-data.frame(submitted_sample_id=u1)
        subj1$mut<-1
      }
    }
    #setup mut_pts
    datasum<-NULL
    datasum<-getDSum(mut$submitted_sample_id, v2=c('submitted_sample_id','mut'))
    #---End original mutation data setup---#
    
    #---setup cnv files---#
    if(file.exists(of.mut2)){
      file.remove(of.mut2)
      # t1<-try(mut_cnv<-as.data.frame(read_feather(of.mut2)))
      # if(class(t1)=='try-error') file.remove(of.mut2)
    }
    subj2<-NULL; mut2<-mut_cnv<-NULL; 
    if(file.exists(mutf2)){#copy number somatic
      # if(rerunMut|!file.exists(of.mut2)){
      if(rerunMut&file.exists(of.mut1)){ #if mut.feather is ready
        mut2<-as.data.frame(fread(mutf2, h=T))
        print(tmp1<-paste('CNV: Proj',i,paste(unique(mut2$platform, collapse=', '))))
        cnv.plat<-c(cnv.plat,tmp1)
        print(paste(unique(mut2$mutation_type, collapse=', ')))
        mut_cnv<-chrm2hs(d1=mut2, g1=g_list3,
                         vs=c("chromosome", "chromosome_start","chromosome_end"), 
                         v2='mutation_type', vid='submitted_sample_id')
        if(!is.null(mut_cnv) && !is.null(nrow(mut_cnv)) && nrow(mut_cnv)>0){
          # write_feather(mut_cnv, of.mut2)
          cnv1c<-colnames(mut_cnv)[colnames(mut_cnv)%in%colnames(mut)]
          mut_cnv<-as.data.frame(mut_cnv)
          mut<-rbindlist(list(mut, mut_cnv[,cnv1c]),fill=T)
          write_feather(mut, of.mut1)#append to mut and overwrite.
        }
      }else if(!rerunMut){
        mut2<-mut_cnv<-mut[grepl('CNV:',mut$mutation_type,fixed=T),]
      }else{#when mut not exist but still require to update mut data
        mut2<-as.data.frame(fread(mutf2, h=T))
        print(tmp1<-paste('CNV: Proj',i,paste(unique(mut2$platform, collapse=', '))))
        cnv.plat<-c(cnv.plat,tmp1)
        mut_cnv<-chrm2hs(d1=mut2, g1=g_list3,
                         vs=c("chromosome", "chromosome_start","chromosome_end"), 
                         v2='mutation_type', vid='submitted_sample_id')
        write_feather(mut_cnv, of.cnv)
      }
      if(!is.null(mut_cnv) && nrow(mut_cnv)>0){
        u1<-unique(mut_cnv$submitted_sample_id)
        subj2<-data.frame(submitted_sample_id=u1)
        subj2$cnv<-1
      }
    }
    #setup mut_cnv_pts
    # datasum<-getDSum(mut2$submitted_sample_id, 
    #                  v2=c('submitted_sample_id','cnv'))
    #---End CNV data setup---#
    
    #---get sample id for mutation data---#
    subj.mut<-unique(mut$submitted_sample_id)
    has.subj.mut<-length(subj.mut)>0
    #---End sample id setup---#
    
    #---set gene expression array data---#
    if (file.exists(of.gexp1) & !rerunExp){
      try1<-try(gexpa<-as.data.frame(read_feather(of.gexp1)))
      if(class(try1)=='try-error')file.remove(of.gexp1)
    }else{gexpa<-NULL}
    subj3<-NULL
    if(file.exists(gexpf1)){
      if(rerunExp|!file.exists(of.gexp1)){
        gexp1<-as.data.frame(fread(gexpf1, h=T))
        gexp1.nm<-c('submitted_sample_id', 'gene_id',
                    'normalized_expression_value')
        if(!gexp1.nm[3]%in%colnames(gexp1)){stop('check gexp1.nm')}
        gexp1t<-gexp1[,gexp1.nm]
        colnames(gexp1t)[3]<-'gexpZ'
        gexpa<-gnmIcgc(d1=gexp1t, v1='gene_id')
        gexpa<-getZ(d1=gexpa)
        write_feather(gexpa, of.gexp1)
      }
      if(has.subj.mut){#samples with both mut and expa
        tmp<-gexpa[gexpa$submitted_sample_id%in%subj.mut,]
        if(nrow(tmp)>0){write_feather(tmp, of.gexp1m)}
      }
      if(!is.null(gexpa) && nrow(gexpa)>0){
        u1<-unique(gexpa$submitted_sample_id)
        subj3<-data.frame(submitted_sample_id=u1)
        subj3$gexpa<-1
      }
    }
    #setup gexpa_pts
    # datasum<-getDSum(gexpa$submitted_sample_id, 
    #                  v2=c('submitted_sample_id','gexpa'))
    #---End gene expression micro array data setup---#
    
    #---set gene expression rnaseq data---#
    if(file.exists(of.gexp2) & !rerunExp){
      try1<-try(gexps<-as.data.frame(read_feather(of.gexp2)))
      if(class(try1)=='try-error')file.remove(of.gexp2)
    }else{gexps<-NULL}
    subj4<-NULL
    if(file.exists(gexpf2)){
      if(rerunExp|!file.exists(of.gexp2)){
        gexp2<-as.data.frame(fread(gexpf2, h=T))
        gexp2.nm1<-c('normalized_expression_value', 'normalized_read_count')
        gexp2.nm1<-gexp2.nm1[gexp2.nm1%in%colnames(gexp2)]
        if(length(gexp2.nm1)==0){stop('check gexp2.nm1')}
        gexp2.nm<-c('submitted_sample_id', 'gene_id',gexp2.nm1[1])
        gexp2t1<-gexp2[,gexp2.nm]
        colnames(gexp2t1)[3]<-'gexpZ'
        gexps<-gnmIcgc(gexp2t1, 'gene_id')
        gexps<-getZ(d1=gexps)
        write_feather(gexps, of.gexp2)
      }
      if(has.subj.mut){
        tmp<-gexps[gexps$submitted_sample_id%in%subj.mut,]
        if(nrow(tmp)>0){write_feather(tmp, of.gexp2m)}
      }
      if(!is.null(gexps) && nrow(gexps)>0){
        u1<-unique(gexps$submitted_sample_id)
        subj4<-data.frame(submitted_sample_id=u1)
        subj4$gexps<-1
      }
    }
    #setup gexps_pts
    # datasum<-getDSum(gexps$submitted_sample_id, 
    #                  v2=c('submitted_sample_id','gexps'))
    #---End gene expression rnaseq data setup---#
    
    #---setup protein data---#
    if (file.exists(of.pexp) & !rerunPExp){
      try1<-try(pExp<-as.data.frame(read_feather(of.pexp)))
      if(class(try1)=='try-error')file.remove(of.pexp)
    }else{pExp<-NULL}
    subj5<-NULL
    if(file.exists(pexpf1)){
      if(rerunPExp|!file.exists(of.pexp)){
        pexp1<-as.data.frame(fread(pexpf1, h=T))
        pexp1.nm<-c('submitted_sample_id', 'gene_name', 
                    'normalized_expression_level')
        pexp1<-pexp1[,pexp1.nm]
        colnames(pexp1)[3]<-'pexpZ'
        pexp1.nm1<-c('gene_name','Hugo_Symbol')
        if(any(unique(pexp1$gene_name)[1:40]%in%g_list$Hugo_Symbol)){
          colnames(pexp1)[colnames(pexp1)=='gene_name']<-'Hugo_Symbol'
          pExp<-pexp1
        }else{
          pExp<-merge(pexp1,pexp_antb[,pexp1.nm1],by='gene_name')
        }
        write_feather(pExp, of.pexp)
      }
      if(has.subj.mut){
        tmp<-pExp[pExp$submitted_sample_id%in%subj.mut,]
        if(nrow(tmp)>0){write_feather(tmp, of.pexpm)}
      }
      if(!is.null(pExp) && nrow(pExp)>0){
        u1<-unique(pExp$submitted_sample_id)
        subj5<-data.frame(submitted_sample_id=u1)
        subj5$pExp<-1
      }
    }#ready to save
    #setup pexp_pts
    #---End protein data setup---#
    
    #---import and set up clin files---#
    if(file.exists(of.clin) & !rerunClin){
      t1<-try(clin<-as.data.frame(read_feather(of.clin)))
      if(class(t1)=='try-error') file.remove(of.clin)
    }
    if(rerunClin|!file.exists(of.clin)){
      clin1<-as.data.frame(fread(clinf1, h=T))
      clin1.nm<-c('submitted_donor_id', 'donor_sex', 
                  'donor_vital_status', 'disease_status_last_followup',
                  'donor_age_at_enrollment', 'donor_age_at_diagnosis',
                  'donor_age_at_last_followup', 'donor_relapse_interval',
                  'donor_tumour_stage_at_diagnosis', 'donor_survival_time')
      clin<-unique(clin1[,colnames(clin1)%in%clin1.nm])
      #categorize age level
      clin.nm.age<-c( 'donor_age_at_enrollment', 'donor_age_at_diagnosis',
                      'donor_age_at_last_followup')
      c1n<-colnames(clin)[colnames(clin)%in%clin.nm.age]
      if(length(c1n)>0){
        clin.agc<-apply(as.matrix(clin[,c1n],ncol=length(c1n)), 2, age2c)
        colnames(clin.agc)<-gsub('age','agec',c1n)
        clin<-data.frame(clin,clin.agc)
      }
      if(file.exists(clinf2)){#donor family
        clin2<-as.data.frame(fread(clinf2, h=T))
        clin2.nm<-c("submitted_donor_id",'donor_has_relative_with_cancer_history')
        clin2<-unique(clin2[, clin2.nm])
        clin2<-unique(clin2[clin2[,2]%in%c('yes','no'),])
        clin<-merge(clin, clin2, by='submitted_donor_id',all.x=T)
      }
      if(file.exists(clinf3)){#donor therapy
        clin3<-as.data.frame(fread(clinf3, h=T))
        clin3.nm<-c("first_therapy_type", "first_therapy_therapeutic_intent",
                    "first_therapy_start_interval", "first_therapy_duration",
                    "first_therapy_response", "second_therapy_type",
                    "second_therapy_therapeutic_intent",
                    "second_therapy_start_interval", "second_therapy_duration",
                    "second_therapy_response","other_therapy" )
        clin3.nm <- clin3.nm[clin3.nm%in%colnames(clin3)]
        clin3<-unique(clin3[,c("submitted_donor_id",clin3.nm)])
        clin3[clin3=='NA'|clin3==''|clin3=='unknown']<-NA
        clin<-merge(clin, clin3, by='submitted_donor_id',all.x=T)
      }
      if(file.exists(clinf4)){#samples
        clin4<-as.data.frame(fread(clinf4,h=T))
        clin4.nm<-c('submitted_donor_id','submitted_sample_id',
                    'submitted_specimen_id')
        clin4<-unique(clin4[,clin4.nm])
        clin<-merge(clin,clin4, by='submitted_donor_id',all.x=T)
      }
      if(file.exists(clinf5)){#specificity
        clin5<-as.data.frame(fread(clinf5,h=T))
        clin5.nm1<-c('submitted_donor_id', 'submitted_specimen_id')
        clin5.nm<-c('specimen_type', 'specimen_donor_treatment_type',
                    'specimen_interval', 'specimen_processing',
                    'tumour_stage_system','tumour_histological_type', 
                    'tumour_stage')
        clin5.nm<-clin5.nm[clin5.nm%in%colnames(clin5)]
        clin5<-unique(clin5[,c(clin5.nm1,clin5.nm)])
        clin<-merge(clin, clin5, by=clin5.nm1, all.x=T)
      }
      #clean up the final clin data
      colnames(clin)<-gsub('donor_','',colnames(clin))
      if(!is.null(subj1)){clin<-merge(clin,subj1, all.x=T)}
      if(!is.null(subj2)){clin<-merge(clin,subj2, all.x=T)}
      if(!is.null(subj3)){clin<-merge(clin,subj3, all.x=T)}
      if(!is.null(subj4)){clin<-merge(clin,subj4, all.x=T)}
      if(!is.null(subj5)){clin<-merge(clin,subj5, all.x=T)}
      write_feather(clin, of.clin)
    }
    #setup clin_pts
    kv1<-c('submitted_id', 'submitted_sample_id','submitted_specimen_id')
    datasum<-getDSum(clin, kvar=kv1, ds=NULL)
    #---End clin setup---#
    
    
    #---save the data summary file---#
    write_feather(datasum, of.sum)
    #---End: save the data summary file---#
    
    #---get ranges for age, os---#
    ifunc1<-function(x){
      x<-x[is.finite(x)&!is.na(x)]
      if(length(x)==0) return(NA) else return(range(x))
    }
    age.rg<-ifunc1(clin[,grepl('age_at_enrollment',colnames(clin))])
    os.rg<-ifunc1(clin[,grepl('survival_time', colnames(clin))])
    #---End get ranges for age, os---#
       
    #---update the summary table---#
    prj.sum<-rbind(prj.sum, prj.sum0)
    prj.sum$id[prjid] <- prjid + icgc.stid
    prj.sum$db[prjid] <- db1
    prj.sum$db_id[prjid] <- prjid
    prj.sum$proj_id[prjid]<-pnm1
    prj.sum$proj_name[prjid]<-prj_code$Project_Name[prj_code[,1]==pnm1]
    prj.sum$country[prjid]<-prj_code$Country[prj_code[,1]==pnm1]
    prj.sum$age_range[prjid]<-paste(age.rg,collapse='~')
    prj.sum$OS_range[prjid]<-paste(os.rg,collapse='~')
    prj.sum$N[prjid]<-length(unique(clin$submitted_id))
    prj.sum$N_mut[prjid]<-sum(datasum$mut==1)
    prj.sum$N_mut_exp_array[prjid]<-sum(datasum$mut==1&datasum$gexpa==1)
    prj.sum$N_mut_exp_seq[prjid]<-sum(datasum$mut==1&datasum$gexps==1)
    prj.sum$m_mutgene[prjid]<-length(unique(mut$Hugo_Symbol))
    print(prj.sum[prjid,])
    #---End: update the summary table---#
    
    prjid<-prjid+1
    print(paste(i, pnm1,'is done.'))
  }
  
  if(T){
    write.csv(prj.sum, 
              file=file.path(topath, gsub('.csv',paste0('_',db1,'.csv'),projfile)), 
              row.names=F)
  }
}
#---End: ICGC---#

#---MSK-IMPACT---#
if(runDB2){
  #load all impack data
  msk0<-get(load(file.path(frompath,db2, 'all_msk.Rdata')))
  lapply(msk0, dim)
  #create the folder for data depository
  if(!dir.exists(path2)){dir.create(path2)}
  
  #pull out mutation data
  mut<-distinct(msk0$mut[,c('Tumor_Sample_Barcode','Hugo_Symbol', 
                   'Variant_Classification', 'HGVSp_Short', 'cDNA_change')])
  colnames(mut)[1]<-conVar2
  colnames(mut)[3:5]<-c('mutation_type','aa_mutation','cds_mutation')
  #pull the samples had mut test
  mut.subj<-as.character(levels(mut[,conVar2]))
  
  #pull out samples had cna test
  cna.subj<-colnames(msk0$cna)[-1]
  cna.subj<-gsub('.', '-', cna.subj, fixed=T)
  #setup cna data
  cna<-melt(msk0$cna, id='Hugo_Symbol');  #proc.time()-t1
  cna<-cna[cna$value!=0,]
  colnames(cna)[2]<-conVar2
  cna$mutation_type<-'CNV:amp'
  cna$mutation_type[cna$value<0]<-'CNV:del'
  cna1<-cna[,c(conVar2, 'Hugo_Symbol','mutation_type')]
  cna1[,conVar2]<-gsub('.','-',cna1[,conVar2], fixed=T)
  
  #pull out samples had fusion test
  fus.subj<-as.character(levels(msk0$fus$Tumor_Sample_Barcode))
  #setup fusion data
  fus<-distinct(msk0$fus[,c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Fusion')])
  fus$mutation_type<-paste0('FUS:',fus$Fusion)
  colnames(fus)[1]<-conVar2
  fus1<-fus[,c(conVar2, 'Hugo_Symbol', 'mutation_type')]
  
  #add fus1 and cna1 into mut
  mut1<-data.frame(rbind_list(list(mut, cna1, fus1))) #a daa.table object
  
  #set up clinical data
  #get mapping between sample id and gene name
  clin.c1<-c('Tumor_Sample_Barcode', 'Hugo_Symbol')
  mut_ext<-distinct(msk0$mut_ext[,clin.c1])
  colnames(mut_ext)[1]<-conVar2
  #get tumor type variables
  clin.c2<-c("CANCER_TYPE", 'CANCER_TYPE_DETAILED',
    "PATIENT_ID",'SAMPLE_ID',"SPECIMEN_PRESERVATION_TYPE", 
             'SAMPLE_TYPE', "MATCHED_STATUS")
  clin.cn<-c("proj_id", "proj_name",
             conVar1, conVar2, "specimen_processing", 'specimen_type',
             'first_therapy_type')
  clin.prm<-distinct(msk0$samp[msk0$samp$SAMPLE_TYPE=='Primary',
                      c('PRIMARY_SITE',clin.c2)])
  clin.met<-distinct(msk0$samp[msk0$samp$SAMPLE_TYPE=='Metastasis',
                      c('METASTATIC_SITE',clin.c2)])
  colnames(clin.prm)[1]<-colnames(clin.met)[1]<-'tumour_histological_type'
  clin<-rbind_list(list(clin.prm, clin.met))
  colnames(clin)[-1]<-clin.cn
  clin$proj_id<-as.character(clin$proj_id)
  unique(clin$proj_id); # table(clin$proj_id)
  clin$proj_id<-gsub(' |Tumor|Tumour|Cancer|Carcinoma|Adenoma', "", clin$proj_id)
  clin$proj_id<-gsub('Skin,Non-|Sarcoma|Stromal|Neuroendocrine','',clin$proj_id)
  tb<-table(clin$proj_id)
  tb.ot<-names(tb)[tb<10]
  clin$proj_id[clin$proj_id%in%tb.ot]<-'Other'
  clin$proj_name[clin$proj_id%in%tb.ot]<-'Other'
  
  #add survial time
  pts<-distinct(msk0$pts[,c("PATIENT_ID",'SEX','SMOKING_HISTORY', 
                            'OS_MONTHS','VITAL_STATUS')])
  colnames(pts)<-c(conVar1, "sex", 'smoking_history','survival_time',
                   'vital_status')
  clin1<-data.frame(full_join(clin,pts, by=conVar1))###check survival_time
  for(i in 1:ncol(clin1)){clin1[,i]<-as.character(clin1[,i])}
  clin1[,"survival_time"]<-as.numeric(clin1[,"survival_time"])
  clin1<-clin1[!is.na(clin1$proj_id),]
  clin1$mut<-0
  clin1$mut[clin1[,conVar2]%in%mut.subj]<-1
  clin1$cna<-0
  clin1$cna[clin1[,conVar2]%in%cna.subj]<-1
  clin1$fus<-0
  clin1$fus[clin1[,conVar2]%in%fus.subj]<-1
  
  #split data by proj_id and proj_name
  kv1<-c(conVar1, conVar2)
  prjid<-1
  uni.prj<-sort(unique(clin1$proj_id))
  prj.sum2<-NULL
  icgc.stid<-nrow(prj.sum)
  for(i in 1:length(uni.prj)){
    pnm1<-uni.prj[i]
    sel1<-clin1$proj_id==uni.prj[i]
    pnm1.nm<-paste(unique(clin1$proj_name[sel1]), collapse=', ')
    clin1s<-data.frame(clin1[sel1,])
    datasum<-getDSum(clin1s, kvar=kv1, ds=NULL)
    samp1<-unique(clin1s[,conVar2])
    mut1s<-mut1[mut1[,conVar2]%in%samp1, ]
    #save data
    ff1<-file.path(topath,db2,i)#output folder name
    if(!dir.exists(ff1)){dir.create(ff1)}
    of.clin<-file.path(ff1, paste0(pnm1,com.clin))
    of.mut1<-file.path(ff1, paste0(pnm1,com.mut1))#including all the mutations
    of.sum<-file.path(ff1, paste0(pnm1,com.sum))#available for heatmap
    write_feather(clin1s, of.clin)
    write_feather(mut1s, of.mut1)
    write_feather(datasum, of.sum)
    
    #---update the summary table---#
    os.rg<-ifunc1(clin1s[,grepl('survival_time', colnames(clin1s))])
    prj.sum2<-rbind(prj.sum2, prj.sum0)
    prj.sum2$id[prjid] <- prjid + icgc.stid
    prj.sum2$db[prjid] <- db2
    prj.sum2$db_id[prjid] <- prjid
    prj.sum2$proj_id[prjid]<-pnm1
    prj.sum2$proj_name[prjid]<-pnm1.nm
    prj.sum2$country[prjid]<-'US'
    prj.sum2$age_range[prjid]<-'NA'
    prj.sum2$OS_range[prjid]<-paste(os.rg,collapse='~')
    prj.sum2$N[prjid]<-length(unique(clin1s[,conVar1]))
    prj.sum2$N_mut[prjid]<-sum(datasum$mut==1)
    prj.sum2$N_mut_exp_array[prjid]<-sum(datasum$mut==1&datasum$gexpa==1)
    prj.sum2$N_mut_exp_seq[prjid]<-sum(datasum$mut==1&datasum$gexps==1)
    prj.sum2$m_mutgene[prjid]<-length(unique(mut$Hugo_Symbol))
    print(prj.sum2[prjid,])
    #---End: update the summary table---#
    prjid<-prjid+1
    print(paste(i, pnm1,'is done.'))
  }
  write.csv(prj.sum2, 
            file=file.path(topath, 
                           gsub('.csv',paste0('_',db2,'.csv'), projfile)), 
            row.names=F) 
}
#---End: MSK-IMPACT---#

#---GENIE---#
if(runDB3){
  #load the GENIE data
  path3file <-dir(path3f); path3file;
  gne<-list()
  for(i in 1:length(path3file)){
    gne[[i]]<-get(load(file.path(path3f, path3file[i])))
  }
  names(gne)<-gsub('genie_|.Rdata','', path3file)
  lapply(gne, dim)
  #pull out mutation data
  mut<-distinct(gne$mut1[,c('Tumor_Sample_Barcode','Hugo_Symbol',
                            'Variant_Classification','HGVSp_Short', 'HGVSc','Chromosome')])
  colnames(mut)[1]<-conVar2
  colnames(mut)[3]<-'mutation_type'
  colnames(mut)[4]<-'aa_mutation'
  colnames(mut)[5]<-'cds_mutation'
  #pull out the samples with mutation test
  mut.subj<-as.character(levels(mut[,conVar2]))
  
  #pull out samples with cna test
  cna<-gne$cna1
  colnames(cna)[-1]<-gsub('.','-', colnames(cna)[-1], fixed=T)
  #pull out samples had cna test
  cna.subj<-colnames(cna)[-1]
  #setup cna data
  cna<-melt(cna, id='Hugo_Symbol');  #proc.time()-t1
  cna<-cna[!is.na(cna$value) & cna$value!=0,]
  table(cna$value)#check different levels of copy number alteration
  colnames(cna)[2]<-conVar2
  cna$mutation_type<-'CNV:amp'
  cna$mutation_type[cna$value<0]<-'CNV:del'
  cna1<-cna[,c(conVar2, 'Hugo_Symbol','mutation_type')]
  cna1[,conVar2]<-gsub('.','-',cna1[,conVar2], fixed=T)
  
  #pull out samples with fus test
  fus.subj<-as.character(levels(gne$fus1$Tumor_Sample_Barcode))
  #setup fusion data
  fus<-distinct(gne$fus1[,c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Fusion')])
  fus$mutation_type<-paste0('FUS:',gsub(' fusion', '', fus$Fusion))
  colnames(fus)[1]<-conVar2
  fus1<-fus[,c(conVar2, 'Hugo_Symbol', 'mutation_type')]

  #add fus1 and cna1 into mut
  mut1<-data.frame(rbind_list(list(mut, cna1, fus1))) #a data.table object

  #set up clinical data
  #get mapping between sample id and gene name
  clin.c1<-c('PATIENT_ID','SAMPLE_ID', 'CANCER_TYPE', 'SEX','PRIMARY_RACE',
             'ETHNICITY','AGE_AT_SEQ_REPORT','SAMPLE_TYPE')
  clin<-distinct(gne$clin[,clin.c1])
  clin$CANCER_TYPE<-as.character(clin$CANCER_TYPE)
  clin$ct<-gsub(' |Cancer|, NOS| Tumor|of|Sarcoma|Carcinoma|Neoplasm', "", 
                clin$CANCER_TYPE)
  clin.c2<-c("ct", "CANCER_TYPE", 
             "PATIENT_ID",'SAMPLE_ID',"SEX", 'AGE_AT_SEQ_REPORT', 
             'SAMPLE_TYPE', "PRIMARY_RACE", 'ETHNICITY')
  clin.cn<-c("proj_id", "proj_name",
             conVar1, conVar2, "sex", 'age_at_enrollment', 'specimen_type',
             'race','ethnicity')
  clin<-clin[,clin.c2]
  colnames(clin)<-clin.cn
  clin$proj_id<-as.character(clin$proj_id)
  tb<-table(clin$proj_id)
  tb.ot<-names(tb)[tb<10]
  clin$proj_id[clin$proj_id%in%tb.ot]<-'Other'
  clin$proj_name[clin$proj_id%in%tb.ot]<-'Other'
  
  #create the folder for data depository
  if(!dir.exists(path3)){dir.create(path3)}
  
  #add survial time
  clin1<-clin
  clin1$survival_time<-clin1$vital_status<-NA
  for(i in 1:ncol(clin1)){clin1[,i]<-as.character(clin1[,i])}
  clin1[,"survival_time"]<-as.numeric(clin1[,"survival_time"])
  clin1<-clin1[!is.na(clin1$proj_id),]
  clin1$mut<-0
  clin1$mut[clin1[,conVar2]%in%mut.subj]<-1
  clin1$cna<-0
  clin1$cna[clin1[,conVar2]%in%cna.subj]<-1
  clin1$fus<-0
  clin1$fus[clin1[,conVar2]%in%fus.subj]<-1
  
  #split data by proj_id and proj_name
  kv1<-c(conVar1, conVar2)
  prjid<-1
  uni.prj<-sort(unique(clin1$proj_id))
  prj.sum3<-NULL
  
  for(i in 1:length(uni.prj)){
    pnm1<-uni.prj[i]
    sel1<-clin1$proj_id==uni.prj[i]
    pnm1.nm<-paste(unique(clin1$proj_name[sel1]), collapse=', ')
    clin1s<-data.frame(clin1[sel1,])
    datasum<-getDSum(clin1s, kvar=kv1, ds=NULL)
    samp1<-unique(clin1s[,conVar2])
    mut1s<-mut1[mut1[,conVar2]%in%samp1, ]
    #save data
    ff1<-file.path(topath,db3,i)#output folder name
    if(!dir.exists(ff1)){dir.create(ff1)}
    of.clin<-file.path(ff1, paste0(pnm1,com.clin))
    of.mut1<-file.path(ff1, paste0(pnm1,com.mut1))#including all the mutations
    of.sum<-file.path(ff1, paste0(pnm1,com.sum))#available for heatmap
    write_feather(clin1s, of.clin)
    write_feather(mut1s, of.mut1)
    write_feather(datasum, of.sum)
    
    #---update the summary table---#
    os.rg<-ifunc1(clin1s[,grepl('survival_time', colnames(clin1s))])
    prj.sum3<-rbind(prj.sum3, prj.sum0)
    prj.sum3$id[prjid] <- prjid
    prj.sum3$db[prjid] <- db3
    prj.sum3$db_id[prjid] <- prjid
    prj.sum3$proj_id[prjid]<-pnm1
    prj.sum3$proj_name[prjid]<-pnm1.nm
    prj.sum3$country[prjid]<-'US'
    prj.sum3$age_range[prjid]<-'NA'
    prj.sum3$OS_range[prjid]<-paste(os.rg,collapse='~')
    prj.sum3$N[prjid]<-length(unique(clin1s[,conVar1]))
    prj.sum3$N_mut[prjid]<-sum(datasum$mut==1)
    prj.sum3$N_mut_exp_array[prjid]<-sum(datasum$mut==1&datasum$gexpa==1)
    prj.sum3$N_mut_exp_seq[prjid]<-sum(datasum$mut==1&datasum$gexps==1)
    prj.sum3$m_mutgene[prjid]<-length(unique(mut$Hugo_Symbol))
    print(prj.sum3[prjid,])
    #---End: update the summary table---#
    prjid<-prjid+1
    print(paste(i, pnm1,'is done.'))
  }
  write.csv(prj.sum3, 
            file=file.path(topath, gsub('.csv',paste0('_',db3,'.csv'),projfile)), 
            row.names=F)
}
#---End: GENIE---#


#---integrate projects_summary file---#
if(TRUE){
  candF<-paste0('projects_summary_',c(db1,db2,db3), '.csv')
  psf<-dir(gsub('/output','',topath))
  psf<-candF[candF%in%psf]
  prj.sum<-NULL
  for(j in 1:length(psf)){
    prj.sum<-rbind(prj.sum,
                   read.csv(file.path(topath,psf[j]), h=T))
  }
  prj.sum$id<-1:nrow(prj.sum)
  write.csv(prj.sum, file=file.path(topath, projfile), row.names=F)
}

#----save platform information---#
if(FALSE){
  write.csv(cnv.plat, file=file.path(topath, 'CNV_plotforms.csv'),row.names=F)
}

Sys.time()-beginTime


#process template
if(FALSE){
  #---setup original mutation files---#
  #---End original mutation data setup---#
  
  #---setup cnv files---#
  #---End CNV data setup---#
  
  #---get sample id for mutation data---#
  #---End sample id setup---#
  
  #---set gene expression array data---#
  #---End gene expression micro array data setup---#
  
  #---set gene expression rnaseq data---#
  #---End gene expression rnaseq data setup---#
  
  #---setup protein data---#
  #---End protein data setup---#
  
  #---import and set up clin files---#
  #setup clin_pts
  kv1<-c('submitted_id', 'submitted_sample_id','submitted_specimen_id')
  datasum<-getDSum(clin, kvar=kv1, ds=NULL)
  #---End clin setup---#
  
  #---save the data summary file---#
  write_feather(datasum, of.sum)
  #---End: save the data summary file---#
  
  #---get ranges for age, os---#
  #---End get ranges for age, os---#
  
  #---update the summary table---#
  prj.sum<-rbind(prj.sum, prj.sum0)
  prj.sum$id[prjid] <- prjid + icgc.stid
  prj.sum$db[prjid] <- db1
  prj.sum$db_id[prjid] <- prjid
  prj.sum$proj_id[prjid]<-pnm1
  prj.sum$proj_name[prjid]<-prj_code$Project_Name[prj_code[,1]==pnm1]
  prj.sum$country[prjid]<-prj_code$Country[prj_code[,1]==pnm1]
  prj.sum$age_range[prjid]<-paste(age.rg,collapse='~')
  prj.sum$OS_range[prjid]<-paste(os.rg,collapse='~')
  prj.sum$N[prjid]<-length(unique(clin$submitted_id))
  prj.sum$N_mut[prjid]<-sum(datasum$mut==1)
  prj.sum$N_mut_exp_array[prjid]<-sum(datasum$mut==1&datasum$gexpa==1)
  prj.sum$N_mut_exp_seq[prjid]<-sum(datasum$mut==1&datasum$gexps==1)
  prj.sum$m_mutgene[prjid]<-length(unique(mut$Hugo_Symbol))
  print(prj.sum[prjid,])
  #---End: update the summary table---#
}
