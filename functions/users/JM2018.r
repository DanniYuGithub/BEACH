print('ppp1')


#all functions
if(TRUE){

library(lavaan) #a package for SEM
library(semPlot)#for SEM plots

library(survminer)#for surv plot with ggsurvplot
library(survival) #for survfit

print2<-function(x, sp='\n'){ #print splitted by sp
  print(strsplit(x, split=sp)[[1]])
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#func: simulate a set of markers in which only 1 baseline marker and 
#multiple chgFromBsl markers are correlated with survival and tumor chg
#
sim.bmk1 <- function(
  outL=NULL,
  n=30, #the number of patients
  bmkType=1, #if 1 cateogrial, then binary variable; else contineous var
  bmkP.norm=data.frame(mu=c(1,3), sd=c(0.35, 0.35)), 
  #parameters for random norm if bmkType!=1
  bmkIns=0.5,#biomarker insidence for group1
  plotbmk=TRUE,#whether plot a histogram of bmk.t1
  lambda=1,  #baseline  hazards
  beta=c(-1,3), #coefficients of the true bmk in TC model and Survival model
  #good candidate is between 1 and 2 (absolute value)
  m.bmk=c(1, #Numb of true chgFromBsl markers
          2, #Numb of false categorical bmks
          2, #Numb of false numeric bmks
          2) #Numb of false chgFromBsl markers
){
  if(TRUE){
    #--------------------------------------------------------------------------#
    #input: parameters for data simulation to create the following objects     #
    # -> bmk: a data.frame with baseline biomarkers in columns                #
    # -> pmk: a list of data.frames for longitudinal                          #
    #          pharmacodynamic (PD) markers, a data.frame is for a PD marker   #
    #rm->pmk.bst: a data.frame of best changes of PD marker for each patient  #
    #                                                                           
    # -> dtte: a data.frame of time-to-event and censored data in two columns  #
    #rm->tchg: a list of data.frames for tumor changes from baseline over time #
    # -> tchg.bst: a data.frame of best tumor change                           #
    #--------------------------------------------------------------------------#
    #output is a list of the ojects described as above.
  }
  
  
  #for robust output in BEACH
  if(is.null(n)|is.null(bmkType)|is.null(bmkP.norm)|
     is.null(bmkIns)|is.null(lambda)|is.null(beta)|
     is.null(m.bmk)){return(NULL)}
  if(length(beta)!=2){
    if(length(beta)==1){
      beta<-rep(beta,2)
    }else{
      print('length(beta) must be 2! Please change')
      return(NULL)
    }
  }
  if(length(m.bmk)!=4){
    print('check m.bmk, 4 integers are required.')
    return(NULL)
  }
  
  
  #adjust samples sizes
  nump <- sum(c(m.bmk,4))
  if(nump>=n){
    print(paste0('n<p! This program increases n from ', n,
                 ' to n+p=', n <- n + nump))
  }
  
  if(plotbmk) par(mfrow=c(2,2))
  
  #~~~1. simulate baseline biomarker~~~#
  bmk.prob<-c(bmkIns, 1-bmkIns)
  bmk.t1 <- sample(x=1:0, size=n, replace=TRUE, prob=bmk.prob)
  bmkIns.est<-mean(bmk.t1)
  if(bmkType!=1){
    if(any(dim(bmkP.norm)!=c(2,2))){return(NULL)}
    bmk.t1<-bmk.t1*rnorm(n, bmkP.norm$mu[1], bmkP.norm$sd[1])+
      (1-bmk.t1)*rnorm(n, bmkP.norm$mu[2], bmkP.norm$sd[2])
    tt1 <- paste0("bmk incidence ",bmkIns,'\n',
                  '(mu=',bmkP.norm$mu[1],', sd=',bmkP.norm$sd[1],
                  '),  (',bmkP.norm$mu[2],', ', bmkP.norm$sd[2], ')')
    if(plotbmk) hist(bmk.t1, main=tt1)
  }else{
    if(plotbmk) hist(bmk.t1, main=paste0("bmk incidence ",bmkIns))
  }
  
  #~~~2. simulate tumor size change (or tumor burden) data~~~#
  #assume the coefficient of bmk.t1 is beta[1]
  #TC ranges from -100 to infinity
  rt1 <- exp(-beta[1]*bmk.t1)
  rt1 <- rt1/max(rt1)/10
  sp1 <-mean(rt1)*50
  tc<-rgamma(n, shape=sp1, rate=rt1)-100 
  if(plotbmk) {
    range(tc); 
    tt2 <- paste0('simmulate tumor change\nfrom rgama with bmk.t1,',
                  ' beta=', beta[1])
    plot(tc~bmk.t1, main=tt2)
  }
  #derive bor from tc
  bor<-rep('PD', n)
  bor[tc<= -30] <- 'PR|CR'
  bor[tc<=20 & tc>-30] <- 'SD'
  bor <- factor(bor, level=c('PD','SD','PR|CR'))
  print(table(bor))
  
  
  #~~~3. simulate survival data~~~#
  #assume weibull distribution
  tc.norm <- (tc-mean(tc))/sd(tc)
  range(sc2 <- exp(beta[2]*bmk.t1 - tc.norm))
  lambda=c(lambda, max(sc2))
  range(tte1 <- rweibull(n, shape=1, scale=lambda[1]* sc2))
  range(tte2 <- rweibull(n, shape=1, scale=lambda[2]))   #censoring time
  tte = pmin(tte1,tte2)  #observed time is max of censored and true
  cnsr = 1*(tte==tte2)   #set to 1 if tte is from censored data
  # head(cbind(tte,tc, cnsr))
  if(plotbmk){
    tt3<-paste0('survial time vs tumor change\npoint size matches bmk.t1')
    plot(tte~tc, pch=cnsr*15+1, cex=bmk.t1-min(bmk.t1)+0.5, main=tt3,col='red')
    legend('topright', legend=c('event','censor'), pch=c(1, 16), col='red')
    if(bmkType==1){
      bmk.bin<-bmk.t1
      leg3<-c('bmk.t1=0','bmk.t1=1')
    }else{
      cut<-round(quantile(bmk.t1, prob=bmkIns.est),3)
      bmk.bin<-1*(bmk.t1>cut)
      leg3<-paste0(c('bmk.t1<=','bmk.t1>'), cut)
    }
    #survival plot
    require(survival)
    sf1<-survfit(Surv(tte, 1-cnsr)~bmk.bin)
    plot(sf1, lty=1:2)
    legend('topright', legend=leg3, lty=1:2)
    cf1 <- survdiff(Surv(tte, 1-cnsr)~bmk.t1)
    pv1 <- round(1 - pchisq(cf1$chisq, 1), 3)
    legend('bottomleft', legend=paste0('p=',pv1), bty='n')
  }
  
  
  #~~~4. simulate the true chgFromBsl markers~~~#
  #for one marker, one change value per patient
  if(m.bmk[1]>=1){
    pmk.t<-data.frame(matrix(NA,nrow=n, ncol=m.bmk[1]))
    colnames(pmk.t)<-paste0('pmk.t',1:m.bmk[1])
    rslp<-0
    while(any(rslp==0)){
      rslp <- round(rnorm(m.bmk[1])) #random slope
    }
    for(j in 1:m.bmk[1]){
      pmk.t[,j] <- rnorm(n, mean=tc*rslp[j], sd=abs(rslp[j])*0.1)
      # boxplot(pmk.t[,j]~bor, ylab=colnames(pmk.t)[j])
    }
    #head(pmk.t)
  }else{pmk.t<-NULL}
  
  #~~~5. simulate the fasle markers~~~#
  if(m.bmk[2]>=1){ #false categorical
    bmk.fc<-data.frame(matrix(NA, nrow=n, ncol=m.bmk[2]))
    colnames(bmk.fc)<-paste0('bmk.fc',1:m.bmk[2])
    r.bi<-round(runif(m.bmk[2], 0.2,0.8),3) #random bmk incidence
    for(j in 1:m.bmk[2]){
      pb<-c(r.bi[j], 1-r.bi[j])
      bmk.fc[,j]<-factor(sample(c(1,0), n, replace=T, prob=pb),levels=c(0,1))
      # boxplot(tc~bmk.fc[,j], xlab=colnames(bmk.fc)[j])
    }
  }else{bmk.fc<-NULL}
  #
  if(m.bmk[3]>=1){ #false numeric
    bmk.fn<-data.frame(matrix(NA, nrow=n, ncol=m.bmk[3]))
    colnames(bmk.fn)<-paste0('bmk.fn',1:m.bmk[3])
    r.mu<-round(rnorm(m.bmk[3])) #random mean and sd
    r.sd<- rchisq(m.bmk[3], 1)
    for(j in 1:m.bmk[3]){
      bmk.fn[,j]<-rnorm(n, mean=r.mu[j], sd=r.sd[j])
      # boxplot(bmk.fn[,j]~bor, ylab=colnames(bmk.fn)[j])
    }
  }else(bmk.fn<-NULL)
  #
  if(m.bmk[4]>=1){ #false chgFromBsl
    pmk.f<-data.frame(matrix(NA, nrow=n, ncol=m.bmk[4]))
    colnames(pmk.f)<-paste0('pmk.f',1:m.bmk[4])
    r.mu<-round(rnorm(m.bmk[4])) #random mean and sd
    r.sd<- rchisq(m.bmk[4], 1)
    for(j in 1:m.bmk[4]){
      pmk.f[,j]<-rnorm(n, mean=r.mu[j], sd=r.sd[j])
      # boxplot(pmk.f[,j]~bor, ylab=colnames(pmk.f)[j])
    }
  }else{pmk.f<-NULL}
  
  #convert response status variable a sa factor
  if(bmkType==1){bmk.t1 <- factor(bmk.t1)}
  dtte<-data.frame(tte,cnsr)
  bmk <- data.frame(bmk.t1, bmk.fc, bmk.fn)
  pmk <- data.frame(pmk.t, pmk.f)
  outL <- list(dtte=dtte, 
               tchg.bst=tc, bor=bor, 
               bmk = bmk, pmk=pmk, 
               slopes=beta)
  outL$all <- data.frame(tc, dtte, bor, bmk, pmk)
  return(outL)
}
#if see error, try to sample again
#sim2 <- sim.bmk1(10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
#plot biomarker in the struction of output from bmkJM.sim1
bmkJM.plot <- function(
  allData=sim1$all, #a data.frame including all variables
  tte.nm='tte',     #variable name for TTE
  cnsr.nm='cnsr',   #variable name for censor or not
  tc.nm='tc',       #variable name for tumor change
  bor.nm = 'bor',   #variable name of best objective response
  bmk.nm=colnames(sim1$bmk),  #variable names for baseline markers
  pmk.nm=colnames(sim1$pmk),  #variable names for longitudinal markers
  #note: all the logitudinal variables must have '.t1', '.t2', ....
  getTotNum=TRUE,   #get the total number of plots and not plot output
  plotAll=TRUE,     #whether plot all graphs or one-by-one
  plotID=NULL       #the index number for a graph
){
  #obj: is a listing object use the same structure as bmkJM.sim1 output
  #outL <- list(dtte=dtte, tchg.bst=tc, bor=bor, 
  #            bmk=bmk, pmk=pmk, slopes=slopes)
  #outL$all <- data.frame(tc, dtte, bor, bmk, pmk)
  
  ask1 <- function(id=NULL){
    if(is.null(id)){
      ot<-readline(prompt='Next plot? y or n: ')
    }else{
      ot<-readline(prompt=paste0('Drow plot ',id,'? y or n: '))
    }
    while(!ot%in%c('y','n')){
      ot<-readline(prompt=paste0(ot, ' is not right. Please input y or n: '))
    }
    return(ot)
  }
  #if plotAll is false, then users will see plots one by one
  
  #pre-define the interactive variable
  nxt <- 'y';
  #pre-define the total number of plots
  totnum <- 0
  #when plotAll==FALSE and getTotNum==FALSE and plotID is specified
  #then only one plot matched to the graph ID will be shown. 
  
  bmk.nm.v <- NULL
  for(i in bmk.nm){
    if(is.numeric(allData[,i])){
      bmk.nm.v <- c(bmk.nm.v, i)
    }
  }
  bmk.nm.c <- bmk.nm[!bmk.nm%in%bmk.nm.v]
  bmk.len.c <- length(bmk.nm.c)
  bmk.len.v <- length(bmk.nm.v)
  
  #1. survival plots grouped by different L or H of a biomarker
  if(all(c(tte.nm, cnsr.nm, bmk.nm)%in%colnames(allData))){
    if (getTotNum){
      totnum <- totnum + bmk.len.c
    }else{
      tte <- allData[,tte.nm]
      cnsr<- allData[,cnsr.nm]
      i <- 1
      while(nxt=="y" & i<=bmk.len.c){
        totnum <- totnum+1
        if(!is.null(plotID) && plotID!=totnum){i <- i+1;next}
        testData <- data.frame(tte=tte, cnsr=cnsr, lev=allData[,bmk.nm.c[i]])
        fit <- survfit(Surv(tte, 1-cnsr) ~ lev, data=testData)
        leg <- paste0(bmk.nm.c[i], ': ')
        # Visualize with survminer
        #ggsurvplot(fit, data=testData, risk.table = TRUE, legend.title=leg)
        print( ggsurvplot(
          fit,                     # survfit object with calculated statistics.
          data = testData,  # data used to fit survival curves. 
          legend.title=leg, 
          risk.table = TRUE,       # show risk table.
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = TRUE,         # show confidence intervals for 
          # point estimaes of survival curves.
          #xlim = c(0,2000),        # present narrower X axis, but not affect
          # survival estimates.
          #break.time.by = 500,     # break X axis in time intervals by 500.
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
          # in legend of risk table
        ))
        if(!plotAll && is.null(plotID)) nxt <- ask1(totnum+1)
        i <- i+1
      }
    }
  }
  
  #2. boxplots of tc grpd by different L or H of a biomarker
  if(all(c(tc.nm, tte.nm, bmk.nm)%in%colnames(allData)) ){
    if (getTotNum){
      totnum <- totnum + bmk.len.c
    }else{
      i <-1; 
      par(mfrow=c(2,2))
      while(nxt=='y' & i<=bmk.len.c){
        totnum <- totnum+1
        if(!is.null(plotID) && plotID!=totnum){i <- i+1;next}
        boxplot(allData[,tc.nm]~allData[,bmk.nm.c[i]],
                ylab=tc.nm, xlab=bmk.nm.c[i])
        boxplot(allData[,tte.nm]~allData[,bmk.nm.c[i]],
                ylab=tte.nm, xlab=bmk.nm.c[i])
        if(!plotAll && is.null(plotID)) nxt <- ask1(totnum+1)
        i <- i+1
      }
      par(mfrow=c(1,1))
    }
  }
  
  #3. scatter plots of tte and tc grpd by different L or H of a biomarker
  if(all(c(tte.nm, cnsr.nm, tc.nm, bmk.nm)%in%colnames(allData))){
    if (getTotNum){
      totnum <- totnum + bmk.len.c
    }else{
      i <-1; 
      par(mfrow=c(1,2))
      while(nxt=='y' & i<=bmk.len.c){
        totnum <- totnum+1
        if(!is.null(plotID) && plotID!=totnum){i <- i+1;next}
        lev <- as.factor(unique(allData[,bmk.nm.c[i]]))
        col1<- as.numeric(lev); names(col1)<-lev;
        pch <- unique(allData[,cnsr.nm])*16
        plot(allData[,tte.nm]~allData[,tc.nm],
             col=allData[,bmk.nm.c[i]], pch=allData[,cnsr.nm]*16,
             ylab=tte.nm, xlab=tc.nm, main=bmk.nm.c[i])
        legend('topright', 
               legend=paste0(#bmk.nm.c[i],'=',
                 rep(lev, each=2), 
                 ', cnsr=', rep(pch/16, 2)),
               pch=rep(pch, 2),
               col=as.numeric(rep(lev, each=2))
        )
        form1<-as.formula(paste0(tte.nm, '~', tc.nm))
        for(k in lev){
          f1<-lm(form1, data=allData[allData[,bmk.nm.c[i]]==k,])
          abline(f1, col=col1[k])
        }
        if(!plotAll && is.null(plotID)) nxt <- ask1(totnum+1)
        i <- i+1
      }
      par(mfrow=c(1,1))
    }
  }
  
  #4. scatter plots of tc vs numeric baseline biomarker
  if(all(c(tc.nm,tte.nm, bmk.nm.v)%in%colnames(allData))){
    if (getTotNum){
      totnum <- totnum + bmk.len.v
    }else{
      i <-1; 
      par(mfrow=c(2,2))
      while(nxt=='y' & i<=bmk.len.v){
        totnum <- totnum+1
        if(!is.null(plotID) && plotID!=totnum){i <- i+1;next}
        plot(allData[,tc.nm]~allData[,bmk.nm.v[i]], pch=16, 
             ylab=tc.nm, xlab=bmk.nm.v[i])
        plot(allData[,tte.nm]~allData[,bmk.nm.v[i]], pch=16, 
             ylab=tte.nm, xlab=bmk.nm.v[i])
        if(!plotAll && is.null(plotID)) nxt <- ask1(totnum+1)
        i <- i+1
      }
    }
  }
  
  #5. scatter plots of tte vs numeric baseline biomarker pch by cnsr
  if(all(c(tte.nm, cnsr.nm, bmk.nm.v)%in%colnames(allData))){
    if (getTotNum){
      totnum <- totnum + bmk.len.v
    }else{
      i<-1;
      pchs <- (allData[,cnsr.nm])*16
      pch1 <- unique(allData[,cnsr.nm])
      par(mfrow=c(2,2))
      while(nxt=='y' & i<=bmk.len.v){
        totnum <- totnum+1
        if(!is.null(plotID) && plotID!=totnum){i <- i+1;next}
        plot(allData[,tte.nm]~allData[,bmk.nm.v[i]], pch=pchs,
             ylab=tte.nm, xlab=bmk.nm.v[i])
        legend('topleft', legend=paste0(cnsr.nm,'=',unique(pchs)/16), xpd=T,
               pch=unique(pchs), col='black', inset=c(0,0))
        form1<-as.formula(paste0(tte.nm, '~', bmk.nm.v[i]))
        for(k in pch1){
          print(k)
          f1<-lm(form1, data=allData[allData[,cnsr.nm]==k,])
          print(f1)
          if(all(!is.na(f1$coefficients)))
          abline(f1, lty=length(pch1)-as.numeric(k))
        }
        i <- i+1
        if(!plotAll && is.null(plotID)) nxt <- ask1(totnum+1)
      }
      par(mfrow=c(1,1))
    }
  }
  
  #6. line plots of longitudinal marker colored by bor level
  if(all(c(bor.nm, pmk.nm)%in%colnames(allData))){
    cols.nm<- as.factor(allData[,bor.nm])
    cols <- as.numeric(cols.nm)
    names(cols) <- cols.nm
    col1.nm <- unique(cols.nm)
    col1 <- as.numeric(col1.nm)
    #note: all the logitudinal variables must have '.t1', '.t2', ....
    pmk.nm.L <- 
      sapply(strsplit(pmk.nm, split='.t'), function(x){x[1]})
    pmk.nm.L1 <- unique(pmk.nm.L)
    if (getTotNum){
      totnum <- totnum + length(pmk.nm.L1)
    }else{
      par(mfrow=c(2,2))
      i<-1; # nxt='y'
      while(nxt=='y' & i<=length(pmk.nm.L1)){
        totnum <- totnum+1
        if(!is.null(plotID) && plotID!=totnum){i <- i+1;next}
        pmk.clnm <- pmk.nm[pmk.nm.L==pmk.nm.L1[i]]
        pmk.1 <- allData[,pmk.clnm]
        if(is.null(ncol(pmk.1))){
          bp<-boxplot(pmk.1~allData[,bor.nm], ylab=pmk.clnm)
        }else{
          xx <- 1:ncol(pmk.1)
          ylim<-range(pmk.1, na.rm=T)
          yy <- as.vector(as.numeric(pmk.1[1,]))
          plot(yy~xx, type='o', col=cols[1], ylim=ylim, pch=16,
               ylab=pmk.nm.L1[i], xlab='time points')
          for(m in 2:nrow(pmk.1)){
            yy <- as.vector(as.numeric(pmk.1[m,]))
            lines(x=xx, y=yy, type='o', pch=16, col=cols[m])
          }
          legend("topleft", legend=col1.nm, text.col=col1,
                 inset=c(0, 0), xpd=T)
        }
        i <- i+1
        if(!plotAll && is.null(plotID)) nxt <- ask1(totnum+1)
      }
      par(mfrow=c(1,1))
    }
  }
  
  if(getTotNum){
    return(totnum)
  }
  #the end of plots
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Set up SEM structures
sem.struct <- function(
  allData=sim1$all, 
  subset=NULL, #T/F or index values to take subset of allData
  model=NULL,
  dtte.nm='tte',
  tc.nm = 'tc',
  bor.nm='bor',
  bmk.nm=colnames(sim1$bmk),
  pmk.nm=colnames(sim1$pmk),
  pmk.cov=NULL,
  cov.nm = NULL,
  getPlot=TRUE,
  est = 'ML',
  rmNode=TRUE, 
  rmTh=c(0.2,0.2), #Threshold of est and std for removed nodes
  fit=NULL,
  shortNM=T
){
  
   if(is.null(fit)){
     if(shortNM){
       colnames(allData)<-gsub('mk.','.',colnames(allData))
       bmk.nm <- gsub('mk.','.', bmk.nm)
       pmk.nm <- gsub('mk.','.', pmk.nm)
     }
    
    #take the subsets
    if(!is.null(subset) & length(subset)==nrow(allData)){
      #and remove all categorical variables
      tmp.bmk.nm<-NULL
      for(i in bmk.nm){
        if(is.numeric(allData[,i])){tmp.bmk.nm<-c(tmp.bmk.nm, i)}
      }
      bmk.nm <- tmp.bmk.nm
      print(bmk.nm)
      nrow1 <- nrow(allData)
      allData <- allData[subset,]
      print(paste('data reduced from', nrow1, 'to', nrow(allData)))
    }
    
    #convert factor into numeric
    all.nm <- c(dtte.nm, tc.nm, bor.nm, bmk.nm, pmk.nm, cov.nm)
    allData2 <- allData[,all.nm]
    for(i  in all.nm){
      if(!is.numeric(allData2[,i])){
        allData2[,i]<-as.numeric(as.factor(allData2[,i]))
      }
    }
    #apply(allData2[,all.nm], 2, is.numeric)
    
    #pre-define covariance.
    lower <- ''
    for(k in 1:length(all.nm)){
      lower <- paste(lower, '\n', 
                     paste(rep('0', k), collapse=' '))
    }
    # print2(lower)
    #crea.cov = getCov(lower, names = all.nm) # Assign variable names
    
    #construct models
    y.nm <- c(dtte.nm, tc.nm)
    mod.y <- paste(y.nm, collapse='+')
    mod.pmk <- paste(pmk.nm, collapse='+'); mod.pmk.l <-length(pmk.nm);
    mod.bmk <- paste(bmk.nm, collapse='+'); mod.bmk.l <-length(bmk.nm);
    cov.bmk <- ''
    for(kb in 1:length(bmk.nm)){
      cov.bmk<-paste(cov.bmk, ifelse(kb==1, '', '\n'),
                     paste(bmk.nm[kb],'~~',
                           paste0(bmk.nm[-kb], collapse='+')))
    }
    
    if(is.null(model)){
      model <- paste0(
        '# latent variables    \n', 
        paste0('eff =~ ', mod.y), '    \n',
        #assume at least one baseline markers are specified, 
        ifelse(mod.bmk.l>1, paste0('bsl =~ ', mod.bmk, '\n'), ''),
        #logitudinal markers are optional
        ifelse(mod.pmk.l>0, 
               ifelse(mod.pmk.l>1, paste0('log =~ ', mod.pmk, '\n'), ''),
               ''),
        '# regressions \n', 
        'eff ~ ', ifelse(mod.bmk.l>1, 'bsl', mod.bmk),
           ifelse(mod.pmk.l>0, 
                  paste0(' + ',  ifelse(mod.pmk.l>1, 'log', mod.pmk)),
                  ''), 
           '\n',
        '# residual covariances    \n',
        ifelse(is.null(dtte.nm)|is.null(tc.nm), '', 
               paste(dtte.nm, '~~', tc.nm)),  '\n',
        ifelse(is.null(pmk.cov), '', 
               paste(pmk.cov, collapse='     \n'))
      )
    }
    print2(model)
    
    fit <- lavaan::sem(model,  data=allData2, 
                       sample.nobs=round(n*.8),
                       #test = "bootstrap",
                       estimator=est)
    #varTable(fit)
    #summary(fit)
  }

  
  #--convert the fit object to semPlotModel structure, not lavaan
  if(class(fit)[1] !=  "semPlotModel"){
    fit.ch <- do.call(semPlotModel, 
                    c(list(fit), modelOpts = list(mplusStd = "std")))
  }else{fit.ch<-fit}
  if(rmNode){
    #print(fit.ch@Vars)
    if(is.null(rmTh)){
      rmTh<-c(quantile(abs(fit.ch@Pars$est), probs=0.25),
              quantile(abs(fit.ch@Pars$std), probs=0.25))
    }
    if(length(rmTh)==1) rmTh <- rep(rmTh,2)
    fit.ch@Pars<-fit.ch@Pars[!is.na(fit.ch@Pars$est)&!is.na(fit.ch@Pars$std),]
    fit.ch@Pars<-fit.ch@Pars[fit.ch@Pars$lhs!=fit.ch@Pars$rhs,]
    fit.ch@Pars<-fit.ch@Pars[abs(fit.ch@Pars$est) > rmTh[1], ]
    fit.ch@Pars<-fit.ch@Pars[abs(fit.ch@Pars$std) > rmTh[2], ]
    
  }
  fit <- fit.ch
  print(fit@Pars)
  
  
  if(getPlot){
    
    if(nrow(fit@Pars)==0){
      x=0; y=0;
      plot(y~x, ylab='', xlab='', axes=F, col='white')
      text(x=x,y=y, labels='0 nodes 0 edges', col='red')
    }else{
    
      # Visualize model
      if(length(unique(allData$bmk.t1))==2){
        bif <- 'bmk.t1'
      }else{bif<-NULL
      }
      bif<-c(bif,colnames(allData)[grepl('.fc',colnames(allData),fixed=T)])
      semPaths(fit,
               style = "lisrel",
               what='std',
               whatLabels='std',
               layout="tree",
               rotation = 3,
               allVars=F,
               intercepts=T,
               residuals=F, 
               thresholds=T,
               nCharNodes=0,
               exoCov=T,
      #         bifactor=bif,
               structural=F,
               optimizeLatRes=T)
      #std: standardized estimate
      #the solid edges represent free parameters estimated from the observed data 
      #and are believed by the investigator to be non-zero.
      #dashed edges represent the fixed parameters  not estimated from the data 
      #and are typically fixed at zero (indicating no relationship between variables)
      #http://userwww.sfsu.edu/efc/classes/biol710/path/SEMwebpage.htm
    }
  }
  
  return(fit)
}



###get selected data names###
get.cnm <- function(dat1, key='bmk.'){
  cnms <- colnames(dat1)
  cnms <- cnms[grepl(key, cnms, fixed=T)]
  return(cnms)
}

###check data names###
check.nm <- function(dat1,   
                     req.vnm=c('tc','tte','cnsr','bmk.', 'pmk.')
                     ){
  cnms <- substring(colnames(dat1), 1, 4)
  mis.vnm <- req.vnm[!req.vnm%in%cnms]
  oo <- 'Pass data names check'
  if(length(mis.vnm)>0){
    oo<-paste0('Fail data names check\nMissed var names: ',
               paste(mis.vnm, collapse=', '))
  }
  return(oo)
}

####plot text###
plottext<-function(tt){
  y=0; x=0; 
  plot(y~x, ylab='', xlab='', axes=F, col='white')
  text(y=y, x=x, labels=tt)
  
}

#convert all binary baseline marker into factor
bin2fac<-function(dat1, varnm){
  for(i in varnam){
    if(length(unique(dat1[,i]))==2){
      dat1[,i] <- factor(dat1[,i])
    }
  }
  return(dat1)
}

} #~~~end of all functions





if(FALSE){
  #BEACH code
  
  
  #^^^^^^^ I. get the data ^^^^^^^^# Figure
  input <- NULL
  #1. Get the data
  input$radio <- c("simulate a new set", "save as default", "use the default", "use uploaded data") [2]
  #2. for data simulation: sample size
  input$text <- ifelse(input$radio=='simulate a new set','30', NULL)
  #3. for data simulation: the true tailoring biomarker type, 1 for binary
  input$slide <- {if(input$radio=='simulate a new set'){1:2}else{NULL}}[2]
  #4. for data simulation: two means for a Normally distributed baseline marker that correlates tc or tte
  input$text2 <- ifelse(input$radio=='simulate a new set','1,3', NULL)
  #5. for data simulation: two std for a Normally distributed baseline marker that correlates tc or tte
  input$text3 <- ifelse(input$radio=='simulate a new set','0.35,0.35', NULL)
  #6. for data simulation: group 1 incidence for baseline marker that correlates tc or tte
  input$text4 <- ifelse(input$radio=='simulate a new set','0.5', NULL)
  #7. for data simulation: coefficients of bmk.t1 from TC model and Survival model
  input$text5 <- ifelse(input$radio=='simulate a new set','-1, 3', NULL)
  #8. for data simulation: Number of changeFromBaseline markers (pmk.t1) correlated with efficacy
  input$text6 <- ifelse(input$radio=='simulate a new set','1', NULL)
  #9. for data simulation: Number of bsl binary markers not correlated with efficacy
  input$text7 <- ifelse(input$radio=='simulate a new set','2', NULL)
  #10. for data simulation: Number of bsl continueous markers not correlated with efficacy
  input$text8 <- ifelse(input$radio=='simulate a new set','2', NULL)
  #11. for data simulation: Number of changeFromBaseline markers not correlated with efficacy
  input$text9 <- ifelse(input$radio=='simulate a new set','2', NULL)
  
  
  if(input$radio=='use the default'){
    load(file='data/sem_sim1_jsm.Rdata')
    allD <<- allD.d
    plottext(check.nm(allD))
  }else if(input$radio=='use uploaded data'){
    allD <- indataset[[1]]
    bmk.nms <- get.cnm(allD, 'bmk.')
    allD <<- bin2fac(allD, bmk.nms)
    plottext(check.nm(allD))
  }else if(input$radio=='save as default'){
    allD.d <- allD
    save(allD.d, file='data/sem_sim1_jsm.Rdata')
    plottext(check.nm(allD.d))
  }else{ #create the 
    sim1 <-   sim.bmk1(
      n=round(as.numeric(input$text)), #the number of patients
      bmkType=input$slide, #if 1 cateogrial, then binary variable; else contineous var
      bmkP.norm=data.frame(mu=as.numeric(strsplit(input$text2, split=',')[[1]]), 
                           sd=as.numeric(strsplit(input$text3, split=',')[[1]])), 
      #parameters for random norm if bmkType!=1
      bmkIns=as.numeric(input$text4),#biomarker insidence for group1
      plotbmk=TRUE,#whether plot a histogram of bmk.t1
      lambda=1,  #baseline  hazards
      #coefficients of the true bmk in TC model and Survival model
      beta=as.numeric(strsplit(input$text5, split=',')[[1]]), 
      #good candidate is between 1 and 2 (absolute value)
      m.bmk=c(as.numeric(input$text6), #Numb of true chgFromBsl markers
              as.numeric(input$text7), #Numb of false categorical bmks
              as.numeric(input$text8), #Numb of false numeric bmks
              as.numeric(input$text9)) #Numb of false chgFromBsl markers
    )
    allD <<- sim1$all
  }
  bmk.nms <<- get.cnm(allD, 'bmk.')
  pmk.nms <<- get.cnm(allD, 'pmk.')
  totPlot<<-bmkJM.plot(
    allData=allD,
    tte.nm='tte', 
    cnsr.nm='cnsr',
    tc.nm='tc',  
    bor.nm = 'bor',
    bmk.nm=bmk.nms,
    pmk.nm=pmk.nms,
    getTotNum = TRUE,
    plotAll=F
  )
  
  
  
  #^^^^^^^ II. see the data ^^^^^^^^# Table
  input <- NULL
  #1. choose variable names
  input$dropdown <- colnames(allD)[1:min(20, ncol(allD))]
  #2. number of rows
  input$slide <- min(50, nrow(allD))
  
  data.frame(allD[1:input$slide,input$dropdown])

  
  #^^^^^^^ III. Visual the data ^^^^^^^^# Figure
  input <- NULL
  #1. plot ID
  input$slide <- seq(1,   totPlot, by=1)[2]
  bmkJM.plot(
    allData=allD,
    tte.nm='tte', 
    cnsr.nm='cnsr',
    tc.nm='tc',  
    bor.nm = 'bor',
    bmk.nm=bmk.nms,
    pmk.nm=pmk.nms,
    getTotNum = F,
    plotAll=F,
    plotID = input$slide
  )
  
  
  #^^^^^^^ IV. analyse the data with strual equation model ^^^^^^^^# Figure
  input <- NULL
  #1. only keep the nodes above threshold
  input$radio <- c('TRUE', 'FALSE')[1]
  #2. thresholds for estimates and standardized estimates
  input$text <- '0.2, 0.2'
  #3. index of baseline markers
  input$slide <- c(0, 1:length(bmk.nms)) [2]
  #4. select baseline markers
  input$dropdown <- if(input$slide>0){bmk.nms[input$slide]}else{bmk.nms}
  #5. select changFromBaseline markers
  input$dropdown2 <- pmk.nms
  #6. shrink nodes names
  input$radio2 <- c('TRUE', 'FALSE')[1]
  
  if(input$slide==0){
    bmk.nm.sel<-input$dropdown
  }else{bmk.nm.sel<-bmk.nms[input$slide]}
  fit0 <- sem.struct(allData=allD, 
                     rmNode=as.logical(input$radio), 
                     rmTh=as.numeric(strsplit(input$text, split=',')[[1]]),
                     bmk.nm=bmk.nm.sel,
                     pmk.nm=input$dropdown2,
                     shortNM=as.logical(input$radio2) )
  title(paste(bmk.nm.sel, collapse=', '))    
}

if(length(indataset)==0){
  load('data/sem_sim1_jsm.Rdata')
  #pre-defined sim1 set is used
  allD <<- sim1$all
  indataset[[1]] <- allD
  print(head(allD))
}
print('ppp2')
