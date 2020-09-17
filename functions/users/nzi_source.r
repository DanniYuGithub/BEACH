if(TRUE){#header
  #*****************************************************************************
  #Eli Lilly and Company - GLOBAL STATISTICAL SCIENCES - PROGRAM
  #CODE NAME (required)                : nzi_source.r
  #PROJECT NAME (required)             : Near-Zero-Informed Bayesian (NZIB)
  #STUDY NUMBER                        : 
  #DESCRIPTION (required)              : 
  #SPECIFICATIONS(optional)            : 
  #VALIDATION TYPE (required)          : 
  #INDEPENDENT REPLICATION (optional)  : 
  #ORIGINAL CODE (required)            : 
  #EXTERNAL CODE FILES THAT ARE NO 
  #RE-USEABLE CODE MODULES             : 
  #SOFTWARE/VERSION# (required)        : R Version 3.2.0
  #INFRASTRUCTURE                      : R studio server
  #DATA INPUT location                 : 
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
  #1.0  Danni Yu                Initial Code created                 01/21/2020
  #2.0  Danni Yu                Last Revision                        06/12/2020
  #3.0  Michelle Carlsen        Peer Review Validation               06/12/2020
  #*****************************************************************************
}
# rm(list=ls())

#program name, track_info
if(TRUE){
  progfile<-'nzi_source.r'
}

#soure functions for NZI Bayesians
if(TRUE){
  library(ggplot2)
  library(gridExtra)
  library(ggrepel)
  library(extraDistr)
  library(cowplot)
  
  #----------------------------------------------------------------------------#
  #### Recommend the maximum P_a is at least ... #####
  recPa <- function(n, dp=6){
    if(is.null(n)){return("0.2")}
    if(is.character(n)){n<-round(as.numeric(n))}
    if(is.na(n)){return("0.2")}
    # p0<-1/n
    # rp<-p0+sqrt(p0*(1-p0)/n)
    rp<-1/n
    for(i in 2:dp){
      rPa<-round(rp, i)
      if(rPa>0) break
    }
    if(i==dp & rPa==0){
      rPa<-as.numeric(paste0('0.',paste(rep(0,dp-1), collapse=''),'1'))
    }
    rPa<-min(rPa,1)
    return(rPa)
  }
  
  
  #----------------------------------------------------------------------------#
  #### Re-set the structure of output #####
  mytab<-function(ttb, ttb.th=0.99){
    colnames(ttb) <- NULL
    rownames(ttb) <- c("Number of Events", #x
                       "Event Rate",        #orr
                       "Superiority", #csr
                       "FPR (type I)", 
                       "FNR (type II)")
    if(all(!is.na(ttb[3,]))){
      return( ttb[,(ttb[3,]<=ttb.th)] )
    }else{return( ttb ) }
  }
  
  
  #----------------------------------------------------------------------------#
  # simulation to estimate the chance of >=1 occurance in N subjects
  # Bayesian approach: 
  nzi_Bayes_sim<-function(
    B=5000, #simulation size 
    N=45,    #sample size
    #p0=NULL,  #a referral rate
    p=0.02,  #a referral rate
    bet=1, a1=NULL #the parameters in Beta distribution
    #plot.p=T, xL=NULL, yL=NULL, bkn=10
  ){
    set.seed(1234567)
    
    if((is.null(p) & is.null(a1))|(is.null(bet))){return(NULL)}
    if(is.null(a1)&!is.null(p)&!is.na(bet)){a1<-p*bet/(1-p)}
    if(!is.null(a1)&is.null(p)&!is.na(bet)){p<-a1/(a1+bet)}

    #using excact probablities for risk estimates
    #upper tail probability in the beta-binomial distribution
    tot.m<-sapply(0:N, s.pbb, n=N, a=a1, b=bet)
    names(tot.m)<-paste0('P[Y>=', 0:N, ']')
    
    #simulate number of events
    if(!is.null(p)){
      tab <- sapply(1:N, function(x, n=N){
        orr <- round(x/n,3) #occurence rate
        csf <- round(1-pbeta(p, a1+x, bet+n-x),3)#under Beta distb only
        pp<-rbeta(B, a1, bet)#randomly generate true rates from Beta distb
        pp1<-pp[pp<=p]; if(length(pp1)==0) pp1<-NA #the true rates < or = bkgrd 
        pp2<-pp[pp>p];  if(length(pp2)==0) pp2<-NA #the true rates > bkgrd
        #approx. Pr( X>= Obs_Incidence x |n, pp1)
        prob.fp <- round(mean(1-pbinom(x-1, n, pp1)), 3)
        #approx. Pr( X < Obs_Incidence x |n, pp2)
        prob.fn <- round(mean(pbinom(x-1, n, pp2)), 3) 
        return(c(x,orr,csf,prob.fp,prob.fn))
      })
      mtt<-mytab(tab, 1.1) 
    }else{mtt<-NULL}
    
    return(list(risk=tot.m[-1], mytab=mtt))
  }

  
  #----------------------------------------------------------------------------#
  # simulation to estimate the chance of >=1 occurance in N
  # Frequentist approach: Binomial distribution
  binom_sim<-function( 
    N=45,     #sample size
    p0=NULL,  #a referral rate
    p=0.02,   #the higher alternative rate or maximumly acceptable rate
    getPlot=F, getPnt=F, getPnt.col='black', getPnt.lty=1,
    out='pp'  #or 'mytab'
  ){
    set.seed(1234567)
    
    if(is.null(p))return(NULL)
    
    o<-1:N #number of events
    
    #simulate number of events
    if(!is.null(p0)){
      tab <- sapply(o, function(x, n=N){
        orr <- round(x/n,3)
        csf <- NA
        prob.fp <- round(1-pbinom(x-1, n, p0), 3)
        prob.fn <- round(pbinom(x-1, n, p), 3) 
        return(c(x,orr,csf,prob.fp,prob.fn))
      })
      mtt<-mytab(tab) 
    }else{mtt<-NULL}
    
    #probability of observing at least #events: P(X>=x)
    if(!is.null(p0)&is.numeric(p0)){
      prob0<-pbinom(o-1, N, p0, lower.tail=F) #P(Y>=y)
      proba<-pbinom(o-1, N, p, lower.tail=F)
      pp<-rbind(prob0, proba)
      colnames(pp)<-o
      rownames(pp)<-paste0('p=', c(p0, p))
      if(getPlot){
        plot(prob0~o, type='l', xlim=c(0,9), ylim=c(0,0.6),
             ylab='Chance observing >= #events', 
             xlab='#events')
        points(proba~o, type='l', col=getPnt.col, lty=getPnt.lty)
        legend('topright', bty='n', legend=paste0(c('p0=','p='), c(p0, p)),
               lty=c(1,2), col=c('black', 'blue'), text.col=c('black', 'blue'))
      }
    }else{
      proba<-pbinom(o-1, N, p, lower.tail=F) #P(Y>=y)
      pp<-proba
      names(pp)<-o
      if(getPlot){
        plot(proba~o, type='l',  col='blue', lty=2)
      }
      if(getPnt){
        points(proba~o, type='l',  col='blue', lty=2)
      }
    }
    return(list(postP=pp, mytab=mtt))
  }
  
  
  #----------------------------------------------------------------------------#
  #check and plot NZI Beta distribution
  fbet<-function(bet, B=5000, p0=0.01, alp=NULL, xlim1=c(0,0.1),
                 ylim1=NULL, is.plot=F, ... ){
    if(is.null(alp)) alp<-p0*bet/(1-p0)
    ss<-rbeta(B, shape1=alp, shape2=bet) 
    if(is.plot){
      hist(ss, freq=F, ylim=ylim1, 
           xlim=xlim1, 
           main=paste0('p=', p0, ', beta=', bet), ...)
    }
    print(c(mean(ss), alp/(alp+bet)))
    oo<-density(ss, from=0, to=1)
    oo$y[oo$y<=0]<-min(oo$y[oo$y>0])
    return(oo)
  }
  #function to generate density values from the prior Beta distribution
  s.dbeta<-function(p,a,b){
    d<-dbeta(p,a,b)
    return(d)
  }
  #function to generate density values from the BetaBinomial distribution
  s.dbb<-function(x, n, a,b){
    d<- extraDistr::dbbinom(x, n, a, b)
    return(d)
  }
  #function to generate upper tail probability in the beta-binomial distribution
  s.pbb<-function(x, n, a,b){
    #P[X>=x]
    pp<- extraDistr::pbbinom(x, n, a, b, FALSE) +
      extraDistr::dbbinom(x, n, a, b)
    return(pp)
  }
  
  
  #----------------------------------------------------------------------------#
  #visualize the prior Beta distribution and the Beta-Binomial distribution
  bbin<-function(
    prior.a="0.5,1,", #parameters in p~Beta(a=, b=) prior distribution
    prior.b="0.5,1",  #parameters in p~Beta(a=, b=) prior distribution
    N11=30,         #planned sample size
    pa=0.2,         #expected average for event rate p
    p.xL='',        #xlim for beta distribution
    x.xL='',        #xlim for beta-binomial distribution
    showBinom=T,    #whether show Biomial distribution in the plot
    drawPlot=TRUE   #if F, only return the ggplot object for predictive distb
  ){
    #clean up paramters
    c2v<-function(ll){as.numeric(unlist(strsplit(ll, split=',', fixed=T)))}
    if(is.character(prior.a)){a<-c2v(prior.a)}else{a<-prior.a}
    if(is.character(prior.b)){b<-c2v(prior.b)}else{b<-prior.b}
    if(is.character(p.xL)){xL.p<-c2v(p.xL)}
    if(is.character(x.xL)){xL.x<-c2v(x.xL)}
    if(is.character(pa)){pa<-as.numeric(pa)}
    if(is.character(N11)){N11<-round(as.numeric(N11))}
    if(length(a)!=length(b)){
      num<-min(length(a), length(b))
      a<-a[1:num]; b<-b[1:num];
    }else{num<-length(a)}
    col2<-rep(c('red', 'purple','magenta','pink', 'brown','valet'),
              each=4)
    col3<-c('black','blue')
    lin2<-rep(3:6,len=length(col2))
    lin3<-1:2;
    
    #function to set x-axis limit
    xyL0<-function(pp, xl){
      if(is.numeric(xl) & length(xl)==2) pp<-pp+ xlim(xl)
      return(pp)
    }
    
    #generate prior distribution
    x1<-0:N11  #use exact probability instead of simulations
    u1<-x1/N11
    ppL<-list() 
    
    legLab1<-legLab2<-NULL #used for math expression in legend
    if(num>0){
      for(i in 1:num){
        #get density values
        ppL[[i]]<- data.frame(p=u1, x=x1,
                              d=s.dbeta(u1, a[i],  b[i]),
                              dx=s.dbb(x1, N11,a[i], b[i]),
                        model=paste0('Beta(a=',a[i],', b=',b[i],')'))
        legLab1<-c(legLab1, #to show math experriosn in legend
                   bquote(paste("Beta(", alpha, "=", .(a[i]),
                                ", ", beta, "=", .(b[i]),')')))
        legLab2<-c(legLab2, 
                   bquote(paste("BetaB(n=", .(N11), ", ", alpha, "=", .(a[i]),
                                ", ", beta, "=", .(b[i]),')')))
      }
      pp <- do.call(rbind, ppL)
    }else{pp<-NULL}
    col1<-rep(col2, len=num); lin1<-rep(lin2, len=num)
    
    #set NZI lable
    model.nzi1<-paste0("NZI Beta(mu=", pa,")")
    
    #combine 3 model data 
    if(is.numeric(pa) && !is.na(pa)){
      tm0<-data.frame(p=u1, x=x1,
                      d=s.dbeta(u1, pa/(1-pa), 1),
                      dx=s.dbb(x1, N11, pa/(1-pa), 1),
                      model=model.nzi1 )
      pp <- rbind(tm0, pp)
      legLab1<-c(bquote(paste("NZI Beta(", mu, "=", .(pa), ")")), legLab1)
      legLab2<-c(bquote(paste("NZIB(n=", .(N11),", ", mu,"=", .(pa), ")")),legLab2)
      col1<-c(col3[1], col1)
      lin1<-c(lin3[1], lin1)
      if(showBinom){
        tm<-data.frame(p=rep(pa,length(x1)), x=x1, d=rep(1, length(x1)),
                       dx=dbinom(x1, N11, pa),
                       model=paste('Binomial(n=',N11,', p=',pa, ')'))
        pp2<-rbind(pp, tm)
        legLab2<-c(legLab2, 
                   bquote(paste('Binomial(n=',.(N11),', p=', .(pa), ')')))
        col1<-c(col1, col3[2])
        lin1<-c(lin1, lin3[2])
      }else{pp2<-pp}
    }else{pp2<-pp}
    
    #generate beta-binomial distribution
    pp2$model<-gsub('NZI Beta(', paste0('NZIB(n=',N11, ', '), pp2$model, fixed=TRUE)
    pp2$model<-gsub('Beta(', paste0('BetaB(n=',N11, ', '), pp2$model, fixed=TRUE)
    pp2$model<-factor(pp2$model, levels=unique(pp2$model))
    pp$model<-factor(pp$model, levels=unique(pp$model))
    
    #clean data for predictive distribution
    pp2<-unique(pp2[,c('x','dx','model')])
  
    #generate distribution plot for prior Beta
    p1<-ggplot(pp, aes(x=p, y=d, group=model,color=model, linetype=model)) + 
      geom_line() + ylab('Density') + 
      theme_bw() + ggtitle('Prior distribution') +
      xlab('p') + 
      scale_color_manual(values=col1, labels=legLab1)+
      scale_linetype_manual(values=lin1, labels=legLab1)
    p1<-xyL0(p1, xL.p)
  
    if(drawPlot){
      #generate distribution plot for Beta-Binomial
      p2<-ggplot(pp2, aes(x=x, y=dx, group=model,color=model, linetype=model)) + 
        geom_line() + ylab('Probability') + 
        theme_bw() + ggtitle('Predictive distribution of events') +
        xlab('n: number of events') + 
        scale_color_manual(values=col1, labels=legLab2)+
        scale_linetype_manual(values=lin1, labels=legLab2)
      p2<-xyL0(p2, xL.x)
      grid.arrange(grobs=list(p1, p2), nrow=1, widths=c(1, 1.1))
    }else{
      #generate distribution plot for Beta-Binomial
      p2<-ggplot(pp2, aes(x=x, y=dx, group=model,color=model, linetype=model)) + 
        geom_line() + 
        theme_bw() + ggtitle('Predictive distribution') +
        scale_color_manual(values=col1, labels=legLab2)+
        scale_linetype_manual(values=lin1, labels=legLab2)
      p2<-xyL0(p2, xL.x) 
      
      return(list(data=pp2, plot=p2))
    }
  }
  
  #----------------------------------------------------------------------------#
  #All-in-one function for NIZB analysis
  nzib <- function( stopAt=NULL, #number of ADR for stopping rule
    N1=45,         #max number of samples
    maxIDR=0.0051, #max expected acceptable occurrence rate
    normIDR=NULL,  #normal occurrence rate, default maxIDR/2 
    plot.type=1,   #1. NZIB; 2. NZIB, Binomial; 3.NZIB, Bayesian; 4.all 3 lines
    bay.prior="c(alpha=0.05, beta=0.8)", #hyper parameter in beta-binomial
    col1=c("black",'red', 'blue'),     #line colors in graph
    rsk.xlim=NULL,    #x-axis range for risk plot
    rsk.ylim=NULL,    #y-axis range for risk plot
    errI.xlim=c(0,5), #x-axis range for type I error plot
    errI.ylim=NULL,   #y-axis range for type I error plot
    errII.xlim=NULL,  #x-axis range for type II error plot
    errII.ylim=NULL,  #y-axis range for type II error plot
    sup.xlim=NULL,    #x-axis range for superiority plot, default=errI.xlim
    sup.ylim=NULL,    #y-axis range for superiority plot
    B1=5000,           #number of simulations
    plot4=FALSE,       #show 4 plots
    goTh="0.05, 0.2", #go/no-go thresholds for type I and type II errors
    noplot=F           #don't create any plot if T
  ){
    #set prior according to plot type
    if(plot.type%in%c(1,2)){bay.prior<-NULL}
    if(plot.type%in%c(2,4)){showbm<-T}else{showbm<-F}
    
    #clean up parameters
    if(is.null(goTh)){
      goTh<-c(0.05, 0.2)
    }else{
      goTh<-as.numeric(strsplit(as.character(goTh), split=',', fixed=T)[[1]])
    }
    if(length(goTh[!is.na(goTh)])<2){goTh<-c(0.05, 0.2)}
    
    if(is.character(B1)){B1<-max(1000, round(as.numeric(B1)))}
    if(is.null(B1)|is.na(B1)){B1<-5000}
    if(is.character(N1)){N1<-max(3, round(as.numeric(N1)))}
    if(!is.null(stopAt) & !is.numeric(stopAt)){
      stopAt<-gsub(' ', '', stopAt, fixed=T)
      if(stopAt!=''){
        stopAt<-unlist(strsplit(as.character(stopAt), split=',', fixed=T))
        stopAt<-as.numeric(stopAt)
        stopAt<-stopAt[!is.na(stopAt)]
        if(length(stopAt)>0){
          stopAt<-pmin(pmax(0, round(stopAt)),N1) #lowest bound
        }
      }
    }
    if(is.character(maxIDR)){maxIDR<-max(min(as.numeric(maxIDR),1),0)}
    if(is.character(normIDR)){normIDR<-max(min(as.numeric(normIDR),maxIDR),0)}
    c2v<-function(xx){
      if(is.character(xx)){xx<-eval(parse(text=paste0('c(',xx,')')))}
      return(xx)
    }
    bay.prior<-c2v(bay.prior)
    rsk.xlim<-c2v(rsk.xlim)
    rsk.ylim<-c2v(rsk.ylim)
    errI.xlim<-c2v(errI.xlim)
    errI.ylim<-c2v(errI.ylim)
    errII.xlim<-c2v(errII.xlim)
    errII.ylim<-c2v(errII.ylim)
    sup.xlim<-c2v(sup.xlim)
    sup.ylim<-c2v(sup.ylim)

    #get stop reference relying number of ADR events
    if(!is.null(stopAt)&is.numeric(stopAt)&length(stopAt)>0){
      sa <- round(stopAt) #require an integer
      saLab<-paste(sa, collapse=',')
    }else{sa<-NULL; saLab<-NULL}
    
    #set default for normIDR
    if(is.null(normIDR)){normIDR <- maxIDR/2}
    
    #get normal expected rate 
    if(plot.type%in%c(1,3)){
      tP1<-NULL
    }else{tP1<-maxIDR}
    
    #order Bayesian hyper parameters
    bay.prior.nm<-c('alpha','beta')
    if(all(bay.prior.nm%in%names(bay.prior))){
      bay.prior<-bay.prior[bay.prior.nm]
    }
    
    #label the methods
    mNm<-c('NZI Bayesian', 'Bayesian', 'Binomial')
    #colors
    if(plot.type==2 & length(col1)>2){
      lab2<-col1[c(1,3)]; lin2<-c(1,3)
    }else{ lab2<-col1; lin2<-1:3 }
    
    #get type I/II errors for 
    b.eg<-eb.eg<-NULL
    #: classic binomial
    if(plot.type%in%c(2,4)){
      b.eg<-binom_sim(N=N1, p0=normIDR, p=maxIDR)
    }
    #: Bayesian
    if(plot.type%in%c(3,4)){
      eb.eg<-nzi_Bayes_sim(B=B1, N=N1,a1=bay.prior[1], bet=bay.prior[2])
    }
    #: NZI Bayesian
    nzi.eg<-nzi_Bayes_sim(B=B1, N=N1, p=maxIDR)
    
    #adjust by plot type
    if(plot.type==4){#show all 3 lines
      eg<-list(nzi.eg, eb.eg, b.eg)
      names(eg)<-mNM0<-mNm
    }else if(plot.type==3){#show lines for Bayesian and NZIB
      eg<-list(nzi.eg, eb.eg)
      names(eg)<-mNM0<-mNm[1:2]
    }else if(plot.type==2){#show lines for binomial and NZIB
      eg<-list(nzi.eg, b.eg)
      names(eg)<-mNM0<-mNm[c(1,3)]
    }else{
      eg<-list(nzi.eg)
      names(eg)<-mNM0<-mNm[1]
    }
    
    #compare type I / II errors
    set.errs<-function(egL=eg){
      numEvts<-fpr<-fnr<-sup<-grp<-NULL
      for(i in 1:length(egL)){
        numEvts<-c(numEvts, egL[[i]]$mytab[1,]) 
        fpr<-c(fpr, egL[[i]]$mytab[4,]) 
        fnr<-c(fnr, egL[[i]]$mytab[5,])
        sup<-c(sup, egL[[i]]$mytab[3,])
        grp<-c(grp, rep(names(egL)[i], ncol(egL[[i]]$mytab)))
      }
      ors<- data.frame( numEvts, fpr, fnr, sup, grp )
      ors$sen<-1-ors$fpr;  
      ors$pwr<-1-ors$fnr
      ors$grp <- factor(ors$grp, levels=names(egL))
      return(ors)
    }
    errs<-set.errs(eg)
    errs.2<-errs[order(errs$grp, errs$numEvts),]
    errs.3<-errs[!is.na(errs$sup),]

    #internal functions for adding stopping thresholds into plot
    xyL<-function(pp, xl, yl, ssa, myLabel, leg='yes'){
      if(is.numeric(xl) & length(xl)==2) pp<-pp+ xlim(xl)
      if(is.numeric(yl) & length(yl)==2) pp<-pp+ ylim(yl)
      if(!is.null(ssa)){
        pp<-pp+ geom_vline(xintercept=ssa,colour='darkgreen', linetype=2)
        pp<-pp+geom_text_repel( data=myLabel[myLabel$numEvts%in%ssa,],
                       aes(label = mylab),
                       box.padding = unit(0.35, "lines"),
                       point.padding = unit(0.3, "lines"),
                       show.legend = FALSE)
      }
      if(leg=='no'){pp<-pp+theme(legend.position='none')}
      
      return(pp)
    }
    
    #1. get risk estimation
    #biomial
    b.rsk <- b.eg$postP
    if(is.matrix(b.rsk)) b.rsk<-b.rsk[paste0('p=',maxIDR),]
    #Bayesian
    eb.rsk <- eb.eg$risk
    #NZIB
    nzi.rsk <- nzi.eg$risk
    #integrate results
    rsko<-cbind(b.rsk, eb.rsk, nzi.rsk)
    rsko<-data.frame(at_least_number_ADR=rownames(rsko), rsko)
    colnames(rsko)[-1]<- paste('risk under', mNM0)
    getid<-function(n1){if(n1<=0){return(NULL)}else{1:n1}}
    rsk<-data.frame(
      numEvts=c(getid(length(b.rsk)), getid(length(eb.rsk)), 
                getid(length(nzi.rsk))),
      risk=c(b.rsk, eb.rsk, nzi.rsk),
      grp=c(rep(mNm[3], length(b.rsk)),
            rep(mNm[2], length(eb.rsk)),
            rep(mNm[1], length(nzi.rsk)) ) )
    rsk$grp <- factor(rsk$grp, levels=mNM0)
    rsk$mylab<-as.character(round(rsk$risk, 3))
    
    if(is.na(normIDR)|plot.type%in%c(1,3)){
      p0lab<-bquote(paste("P(wrong stop | p"<="", .(maxIDR), ")"))
    }else{
      p0lab<-bquote(paste("P(wrong stop | p"<="", .(maxIDR), ", ",
                           p[0], '=', .(normIDR), ")" ))
    }
    
    #1. risk plot
    p1<-ggplot(data=rsk, aes(x=numEvts, y=risk, color=grp, linetype=grp)) + 
      geom_line() + theme_bw() + theme(plot.title = element_text(size = 11)) +
      ggtitle(bquote(paste("Prior study risk profile, ",
                           p[a], '=', .(maxIDR)))) +
      xlab("#IADR") + 
      ylab(expression(paste("P(Y">="#IADR)"))) +
      scale_linetype('') + scale_colour_manual('', values=lab2) 
    p1<-xyL(pp=p1, xl=rsk.xlim, yl=rsk.ylim, ssa=sa, myLabel=rsk, leg='no')
    
    #2. type I
    errs.2$mylab<-as.character(round(errs.2$fpr, 3))
    p2<-ggplot(data=errs.2, aes(x=numEvts, y=fpr, color=grp, linetype=grp)) + 
      geom_hline(yintercept=goTh[1], color='tan3') +
      geom_line() + theme_bw() + theme(plot.title = element_text(size = 11))+
      ggtitle(p0lab) + 
      scale_colour_manual( values=lab2) +
       ylab('Type I error')  
    if(is.null(saLab)){ p2<- p2 + xlab("#IADR")
    }else{ p2<- p2+xlab(bquote(paste("if stopping at #IADR=", .(saLab)))) }  
    
    #graph
    if(plot4){
      p2<-xyL(pp=p2, xl=errI.xlim, yl=errI.ylim, ssa=sa, myLabel=errs.2,
              leg='no') # + labs(colour='For A-E', linetype='For A-E')
      
      #3. power
      errs.2$mylab<-as.character(round(errs.2$fnr, 3))
      p3<-ggplot(data=errs.2, aes(x=numEvts, y=fnr, color=grp, linetype=grp)) + 
        geom_hline(yintercept=goTh[2], color='tan3') +
        geom_line() + theme_bw() + theme(plot.title = element_text(size = 11))+
        ggtitle(paste0("P(wrong go | p>", maxIDR,")")) + 
        ylab('type II error') + 
        xlab(bquote(paste("if going with #IADR=", .(saLab))) ) +
        scale_linetype('For A-E') + 
        scale_colour_manual('For A-E', values=lab2) #+
        # theme(legend.position='none')
      p3<-xyL(pp=p3, xl=errII.xlim, yl=errII.ylim, ssa=sa, myLabel=errs.2, 
              leg='yes') 
      
      #add type I vs type II
      #add bg color
      y.max<-max(errs.2$fpr, na.rm=T);
      y.max2<-max(y.max, goTh[1])
      x.max<-max(1, max(errs.2$fnr, na.rm=T));
      if(y.max<y.max2){
        plg<-data.frame(gp=c(rep('go',4),         rep('stop', 4)),
                        x=c(0,goTh[2],goTh[2],0,  goTh[2],goTh[2], x.max, x.max), 
                        y=c(0,0, y.max2,y.max2,   0,goTh[1],goTh[1], 0) )
        plg$gp<-factor(plg$gp, levels=c('go','stop'))
        plg.col<-c('green','magenta')
      }else{
        plg<-data.frame(gp=c(rep('go',4),   rep('stop', 4),    rep('n/a',4)), 
                      x=c(0,goTh[2],goTh[2],0,  goTh[2],goTh[2], x.max, x.max,
                          goTh[2],goTh[2],x.max,x.max), 
                      y=c(0,0,y.max2,y.max2,    0,goTh[1],goTh[1], 0,
                          goTh[1],y.max2,y.max2, goTh[1]) )
        plg$gp<-factor(plg$gp, levels=c('go', 'stop', 'n/a'))
        plg.col<-c('green', 'magenta', 'gray80')
      }
      #
      p6<-ggplot() + 
        geom_polygon(data=plg, aes(x=x, y=y, group=gp, fill=gp), 
                     colour='white', alpha=0.3, size=1)+
        scale_fill_manual("   ", values = plg.col) +
        geom_line(data=errs.2, aes(x=fnr, y=fpr, color=grp, linetype=grp)) +
        geom_point(data=errs.2, aes(x=fnr, y=fpr, color=grp)) + 
        geom_point(data=errs.2[errs.2$numEvts%in%sa,], 
                   aes(x=fnr, y=fpr, color=grp), size=5, shape=0) + 
        theme_bw() + theme(plot.title = element_text(size = 11))+
        ggtitle("type I vs type II error") + 
        ylab(p0lab) + 
        xlab(paste(paste0("P(wrong go | p>", maxIDR,")"))) +
        scale_linetype('')+
        scale_colour_manual('', values=lab2) +
        theme(legend.position=c(0.95, 0.95), 
              legend.background = element_rect(fill = "transparent"),
              legend.justification=c(1,1), 
              legend.title=element_blank())+
        guides(linetype=F, col=F,
               shape = guide_legend(override.aes = list(size = 3)))
      
      
      #4. superiority prob
      errs.3$mylab<-as.character(round(errs.3$sup, 3))
      p4<-ggplot(data=errs.3, 
             aes(x=numEvts, y=sup, color=grp, linetype=grp)) + 
        geom_line() + theme_bw() + theme(plot.title = element_text(size = 11))+
        ggtitle('Posterior Probability') + 
        ylab(bquote(paste("P(p\'>", .(maxIDR), " | #IADR)"))) + xlab('#IADR') +
        scale_linetype('') + scale_colour_manual('', values=lab2)
      p4<-xyL(pp=p4, xl=sup.xlim, yl=sup.ylim, ssa=sa, myLabel=errs.3, leg='no')
      
      #5. show predictive distribution
      b.out<-bbin(prior.a=bay.prior[1], prior.b=bay.prior[2], 
                     N11=N1, pa=maxIDR, showBinom=showbm, drawPlot=F)
      errs.5<-b.out$data[,c('x','dx', 'model')] 
      #add variables to show text
      errs.5$numEvts<-errs.5$x
      errs.5$mylab<-as.character(round(errs.5$dx,3))
      #revise plot
      p5<- b.out$plot+ ylab('P(Y=#IADR)') +
        labs(colour='For E', linetype='For F')+
        theme(plot.title = element_text(size = 11))+
        xlab('#IADR')
      p5<-xyL(pp=p5, xl=errII.xlim, yl=errII.ylim, ssa=sa, myLabel=errs.5)
        
      #layout the 5 graphs
      if(!noplot){
        x2loc<-0.435; x3loc<-0.28; x4loc <- x3loc*2
        show(ggdraw() + 
              draw_plot(p1, x=0,     y=.5, width=x3loc, height=.5) + 
              draw_plot(p2, x=x3loc, y=.5, width=x4loc-x3loc, height=.5) +
              draw_plot(p3, x=x4loc, y=.5, width=1-x4loc, height=.5) + 
              draw_plot(p6, x=0,     y=0,  width=x3loc, height=.5) + 
              draw_plot(p4, x=x3loc, y=0,  width=x4loc-x3loc, height=.5) + 
              draw_plot(p5, x=x4loc, y=0,  width=1-x4loc, height=.5) + 
              draw_plot_label(c('A','B', 'C', 'D', 'E', 'F'), 
                              x=c(0, x3loc+0.01, x4loc+0.01, 0, x3loc, x4loc), 
                              y=c(1,1,1, 0.5,0.5,0.5), 
                              size=15))
      }
    }else{
      p2<-xyL(pp=p2, xl=errI.xlim, yl=errI.ylim, ssa=sa, myLabel=errs.2)+
         labs(colour='', linetype='')
      #layout the 2 graphs
      if(!noplot){
        x2loc<-0.435; 
        show( ggdraw() + 
                draw_plot(p1, x=0, y=.5, width=x2loc, height=.5) + 
                draw_plot(p2, x=x2loc+0.01, y=.5, width=1-x2loc, height=.5) +
                draw_plot_label(c('A','B'), 
                                x=c(0,x2loc+0.01), y=c(1,1), size=15))
      }
    }
    
    oth<-errs.2[,!colnames(errs.2)%in%'mylab']
    sel<-oth[oth$numEvts%in%sa, c('grp','numEvts', 'fpr', 'fnr')]
    if(nrow(sel)>0){
      sel$note<-"n/a"
      sel$note[sel$fnr<=goTh[2]]<-'go'
      sel$note[sel$fn>goTh[2]&sel$fpr<goTh[1]]<-'stop'
    }
    return(list(risk=rsko, oth=oth, sel=sel))
  }#End nzib function
  
  

} #end for all source functions


if(F){
  # main analysis
  rsk.tb1<- 
    nzib(stopAt="1, 4",   #1. stopping threshold
         N1=32,           #2. max number of samples
       maxIDR=0.01,       #3. max expected acceptable occurrence rate
       plot.type=1, 
       normIDR='0.01',    #4. normal occurrence rate, default maxIDR/2 
       bay.prior=c(alpha=0.2, beta=0.8), #5. hyper parameter in beta-binomial
       rsk.xlim="",       #6. x-axis range for risk plot
       rsk.ylim=NULL,     #7. y-axis range for risk plot
       errI.xlim=NULL,    #8. x-axis range for type I error plot
       plot4=T,           #9. whether show all 4 available plots
       B1=5000,           #10. number of simulations
       goTh="0.05, 0.2",  #11. threshold for type I and II errors
       noplot=F           #don't create any plot
  )
  #input$text:  1. stopping threshold for number of Adverse Drug Reaction (#IADR)
  #input$text2: 2. max number of samples
  #input$text3: 3. max expected acceptable occurence rate of ADR
  #input$text4: 4. expected occurence rate of ADR considered as normal, only useful for frequentist Binomial
  #input$text5: 5. Hyper parameters for traditional Bayesian
  #input$text6: 6. X-axis range in risk plot
  #input$text7: 7. Y-axis range in risk plot
  #input$text8: 8. X-axis range in type I error plot and Superiority plot
  #input$radio: 9. whether show all 4 available plots
  #input$text9: 10. number of simulations
  
  
  
  bbin(
    prior.a="1,0.5, 0.0051", #1. parameters in p~Beta(a=, b=) prior distribution
    prior.b="1,0.5,1", #2. parameters in p~Beta(a=, b=) prior distribution
    N1=45,            #3. planned sample size
    pa='0.01',           #4. expected average for event rate p
    p.xL='',        #5. xlim for beta distribution
    x.xL=''         #6. xlim for beta-binomial distribution
  )
  #input$text:  1. parameters in p~Beta(a=, b=) prior distribution
  #input$text2: 2. parameters in p~Beta(a=, b=) prior distribution
  #input$text3: 3. planned sample size
 
  
}






