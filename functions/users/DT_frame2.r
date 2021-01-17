
#define global variables
if(TRUE){
  fL<<-'https://www.fda.gov/drugs/informationondrugs/approveddrugs/ucm279174.htm'
  FDA.Link <<- fL
}



#Begin 0. ---------------------------------------------------------------------#
if(TRUE){
  #function: get the estimate of p=Pr(CDA+|CDx+), the TRUE positive rate, or sensitivity
  #output: p_hat and Var_p_hat
  TPR <<- function(
    m1,  #number of enrolled patients who are CTA+ and randomly assigned to treatment or control
    m,   #the total number of screened patients with valid CTA results (eligible for bridging study)
    n11, #number of CTA+ and CDx+ out of n1
    n1,  #number of CTA+ patients enrolled in bridging (concordance) study
    n01, #number of CTA- and CDx+ out of n0
    n0   #number of CTA- patients enrolled in bridging (concordance) study
  ){
    if(FALSE){
      phi_hat    =0.09  #Pr(CTA+)
      psi_hat_11 =0.95  #Pr(CDx+|CTA+), sensitivity
      psi_hat_10 =0.95  #Pr(CDx+|CTA-), 1-specificity
    }
    if(m1>m)   {print('Error: m1 > m. So m1 is set to be m.'); m1 <- m}
    if(n1>m1)   {print('Error: n1 > m1. So n1 is set to be m1'); n1 <- m1}
    if(n11>n1) {print('Error: n11 > n1. So n11 is set to be n1.'); n11 <- n1}
    if(n0>(m-n1)) {print('Error: n0 > m-n1. So n0 is set to be m-n1.'); n0 <- m-n1}
    if(n01>n0) {print('Error: n01 > n0. So n01 is set to be n0.'); n01 <- n0}
    
    phi_hat     <- m1/m
    var_phi_hat <- phi_hat * (1-phi_hat) / m
    psi_hat_11     <- n11/n1
    var_psi_hat_11 <- psi_hat_11 * (1-psi_hat_11) / n1
    psi_hat_10     <- n01/n0
    var_psi_hat_10 <- psi_hat_10 * (1-psi_hat_10) / n0
    
    p_hat <- (phi_hat*psi_hat_11)/(phi_hat*psi_hat_11+(1-phi_hat)*psi_hat_10)
    var_p_hat <- (p_hat*(1-p_hat))^2 *
      (  var_phi_hat/((phi_hat*(1-phi_hat))^2) +
           var_psi_hat_11/(psi_hat_11^2) + 
           var_psi_hat_10/(psi_hat_10^2)   )
    se_p_hat <- sqrt(var_p_hat)
    
    return(data.frame(phat=p_hat, var_phat=var_p_hat, 
                      phat_95CI_low=max(0, p_hat-1.96*se_p_hat),
                      phat_95CI_high=min(1, p_hat+1.96*se_p_hat) ))
  } #end of TPR
  
  #function: get estimate of efficacy of briding study
  EFF <<- function(
    delta1_hat=0.5,  #efficacy in subgroup CDA+ and CDx+
    var_delta1_hat,  
    c = 1, #proportion of delta2/delta1, where delta2 is efficacy in CDA- and CDx+
    p_hat,    #p=Pr(CDA+|CDx+)   
    var_p_hat # obtain from the function of PPV   
  ){
    
    #a vector of bridging efficacy given different c
    delta_CDx_pos_hat <- ((1-c)*p_hat+c)*delta1_hat
    var_delta_CDx_pos_hat <- (2*p_hat^2-2*p_hat+1)*var_delta1_hat +
      ( ((1-c)^2)*(delta1_hat^2) + 2*var_delta1_hat ) * var_p_hat
    se_delta_CDx_pos_hat <- sqrt(var_delta_CDx_pos_hat)
    
    return(data.frame(dhat=delta_CDx_pos_hat,
                      var_dhat=var_delta_CDx_pos_hat,
                      dhat_95CI_low=delta_CDx_pos_hat-1.96*se_delta_CDx_pos_hat,
                      dhat_95CI_high=delta_CDx_pos_hat+1.96*se_delta_CDx_pos_hat))
  }#end of EFF
  
  #calculate the variance
  get_var_delta1_hat <<- function(
    orr,  #single-arm objective response rate in subgroup CDA+ and CDx+
    n11,  #number of CTA+ and CDx+ out of n1
    hr,   #two-arm hazard ratio (treatmentHazardRatio/controlHazardRatio)
    nExpectEvt=round(n11*c(0.4, 0.4)), #expected number of events (i.e. death) in 
    #treatment and control
    method='ORR' #or 'logHR'
  ){
    if(method=='ORR'){
      delta1_hat <- orr
      var1 <- delta1_hat*(1-delta1_hat)/n11
      return(data.frame(est=delta1_hat, estVar=var1))
    }else if (method=='logHR'){
      #assume proportional hazards: the hr is consistent over time and
      #any differences are due to random sampling.
      #assume the distribution of the time-to-event measure has an
      #exponential distribution
      #logrank approach: as part of the KM calculation, compute the number
      #of observed events in each group(Oa and Ob), and the number of 
      #expected events assuming a null hypothesis of no difference in 
      #survival (Ea and Eb)
      #http://aac.asm.org/content/48/8/2787.full
      if(sum(nExpectEvt)>n11){
        r1 <- nExpectEvt[1] /sum(Control)
        nr1 <- round(r1*n11)
        nr2 <- n11-nr1
        nExpectEvt <- c(nr1, nr2)
      }
      var_loghr <-sum(1/nExpectEvt)
      return(data.frame(est=log(hr), estVar=var_loghr))
    }else{
      print('Currently only two methods (ORR or logHR) are avaible.')
      var1 <- NA
      return(data.frame(est=NA, estVar=NA))
    }
    
  } #end of ge_var_delta1_hat
  
  #estimate dhat and phat if the efficacy is ORR
  est_ORR <<- function(
    m1=60, m=round(60/0.10), n11=29, n1=30, n01=1, n0=10,
    orr=0.4, c=0.5
  ){    
    tpr1 <- TPR(m1=m1, m=m, n11=n11, n1=n1, n01=n01, n0=n0)
    var_d1 <- get_var_delta1_hat(orr=orr, n11=n11, method='ORR')
    eff1 <- EFF(delta1_hat=var_d1$est[1], var_delta1_hat=var_d1$estVar, c=c, 
                p_hat=tpr1$phat, var_p_hat=tpr1$var_phat)
    eff1$dhat_95CI_low <- max(0, eff1$dhat_95CI_low)
    eff1$dhat_95CI_high <- min(1, eff1$dhat_95CI_high)
    
    return(list(tpr1=tpr1, eff1=eff1))
  }#end of estimate dhat and phat
  
  
  #plot
  plot_dp <<- function(
    tpr0,     #the output for TPR
    eff0,     #the otuput EFF
    xlim0=c(0,1), #range of X-axis
    ylim0=NULL,   #range of Y-axis
    xlab0= 'Concordance Pr(CTA+|CDx+)', 
    ylab0= 'ORR in CDx+ patients',
    add=FALSE #add the est bars to current plot if TRUE
  ){
    if(is.null(ylim0)){
      ylim0 <- range(as.vector(eff0[, c('dhat_95CI_low', 'dhat_95CI_high')]))
    }
    if(is.character(ylim0) & grepl(',', ylim0)){
      y12<-strsplit(ylim0, split=',', fixed=TRUE)[[1]]
      ylim0<-c(as.numeric(y12))
    }
    if(!add){
      plot(eff0$dhat~tpr0$phat, type='b', col='black', pch=16, 
           xlab=xlab0, ylab=ylab0, xlim=xlim0, ylim=ylim0)
    }else{
      points(eff0$dhat~tpr0$phat, type='b', col='black', pch=16,
             xlab=xlab0, ylab=ylab0, xlim=xlim0, ylim=ylim0)
    }
    #add 95%CI for phat
    arrows(x0=tpr0$phat_95CI_low,   y0=eff0$dhat, 
           x1=tpr0$phat_95CI_high, y1=eff0$dhat, 
           code=3, angle=90, length=0.15)
    #add 95%CI for dhat
    arrows(y0=eff0$dhat_95CI_low,   x0=tpr0$phat, 
           y1=eff0$dhat_95CI_high, x1=tpr0$phat, 
           code=3, angle=90, length=0.15)
    
  }#end of plot_dp
  
}#End 0. ---------------------------------------------------------------------#



#Begin 1. ---------------------------------------------------------------------#
if(TRUE){
  #Optimal basket design phase 1a/b study
  #Clinical Benifit and Biological Response are {0, 1} for {no or yes}. 
  #Output@: ranked table {dose, Disease, baselineBiomarker, expectedUtility}
  #Output@: tree plot
  DT_frame <- function(
    blBMK=c('B0','B1'),          #Baseline Biomarkers, required input
    tumorType=c('T0','T1','T2'), #Types of tumors, required input
    dose=c(20, 100),             #numeric values of dose levels, required input
    #parameters
    prior_ti=c(0.1, 0.2,  0.2, 0.3, 0.05, 0.15),#PrioInfo Tumor Incidence: 
    #length is length(blBMK)*length(tumorType)
    #values obtained from prior knowledge
    prior_prop=c(0.1, 0.9),      #proportion of patients in the subgroup of dose 
    #and bmk length is either length(dose) or
    #length(dose)*length(blBMK)*length(tumorType)
    #values obtained from prior knowledge
    prob_stop0 = c(0.75, 0.05),  #Prob(not stop|dose) matching to the levels in 
    #dose length is length(dose), the 
    #proportions obtained from early phase trials
    prob_BR1   = c(0.1, 0.75),   #Prob(BioResp=1|stop=0, dose, tumor, bmk), 
    #length is either length(dose) or 
    #length(dose)*length(blBMK)*length(tumorType)
    #values obtained from early phase trials  
    prob_CB1_BR = c(0.8, 0.1),   #Prob(ClinBenefit=1|BioResp=1) and 
    #Prob(CB=1|BR=0). The lengh is either 2, or
    #2*length(dose)*length(blBMK)*length(tumorType)
    showTree   = TRUE,           #if FASLE only return the tree table
    showProb   = TRUE,           #if TRUE show the probablities on the tree and 
    #return the probability table
    showBar    = TRUE,           #show the barplot of expected U(dose|{T,B})
    #other args for plotting
    th.arrow        = 0.8,       #horizontal space between an arrow and target
    th.utDB    = 1,              #vertical space between dose sign and utility
    topRatio   = 0.2,            #the top ratio of joint p-values (or utilities) 
    #that need to be colored
    topCol     = 'red',          #the color for the top joint p-values  
    payoff     = c(100, -100)    #payoff value for CB=1 and CB=0
  ){
    #an internal function trunc values to 0, 1
    trunc01 <- function(val){
      val[val>1]<-1
      val[val<-0]<-0
      return(val)
    }
    
    #cleanup the payoff values
    if(is.character(payoff)){
      payoff <- as.numeric(strsplit(payoff, split=',', fixed=TRUE)[[1]])
      if(length(payoff)==1){
        payoff<-c(max(0, payoff), min(0, payoff))
      }else{payoff <- payoff[1:2]}
    }
    
    #cleanup the top ratio
    if(is.null(topRatio)) topRatio <- 0.2
    topRatio <- trunc01( as.numeric(topRatio) )
    
    #cleanup the input
    if(length(blBMK)==1){
      blBMK<-strsplit(as.character(blBMK), split=',', fixed=T)[[1]]
    }
    if(length(tumorType)==1){
      tumorType<-strsplit(as.character(tumorType), split=',', fixed=T)[[1]]
    }
    if(length(dose)==1){
      dose<-strsplit(as.character(dose), split=',', fixed=T)[[1]]
    }
    if(length(prior_ti)==1){
      prior_ti<-as.numeric(strsplit(as.character(prior_ti), split=',', fixed=T)[[1]])
    }
    if(length(prior_prop)==1){
      prior_prop<-as.numeric(strsplit(as.character(prior_prop), split=',', fixed=T)[[1]])
    }
    if(length(prob_stop0)==1){
      prob_stop0<-as.numeric(strsplit(as.character(prob_stop0), split=',', fixed=T)[[1]])
    }
    if(length(prob_BR1)==1){
      prob_BR1<-as.numeric(strsplit(as.character(prob_BR1), split=',', fixed=T)[[1]])
    }
    if(length(prob_CB1_BR)==1){
      prob_CB1_BR<-as.numeric(strsplit(as.character(prob_CB1_BR), split=',', fixed=T)[[1]])
    }
    
    #check wether the probablities are matching to the actions
    if(showProb){
      con1 <- length(prior_ti) == length(blBMK)*length(tumorType)
      con2 <- length(prior_prop) == length(dose) | 
        length(prior_prop)==length(dose)*length(blBMK)*length(tumorType)
      con3 <- length(prob_stop0)==length(dose) | 
        length(prob_stop0)==length(dose)*length(blBMK)*length(tumorType)
      con4 <- length(prob_BR1) == length(dose) | 
        length(prob_BR1)==length(dose)*length(blBMK)*length(tumorType)
      con5 <- (length(prob_CB1_BR)==2)| (
        length(prob_CB1_BR)==2*length(dose)*length(blBMK)*length(tumorType)
      )
      if(!all(c(con1, con2, con3, con4, con5))){showProb<-FALSE}
    }
    
    #construct the output matrix
    if(TRUE){
      numL.ClinBenif<- 2
      numL.BioResp  <- 2
      numL.stop      <- 2
      numL.dose     <- length(dose)
      numL.tumorType<- length(tumorType)
      numL.blBMK    <- length(blBMK)
      
      ClinBenif <- rep(c('yes', 'no'), 
                       times=numL.BioResp*numL.stop*numL.dose*numL.tumorType*numL.blBMK)  
      BioResp   <- rep(rep(c('yes', 'no'), each=numL.ClinBenif),
                       times=numL.stop*numL.dose*numL.tumorType*numL.blBMK)
      stop0       <- rep(rep(c('yes', 'no'), each=numL.ClinBenif*numL.BioResp),
                         times=numL.dose*numL.tumorType*numL.blBMK)
      dose.in   <- rep(rep(dose, each=numL.ClinBenif*numL.BioResp*numL.stop),
                       times=numL.tumorType*numL.blBMK)
      BMK.in    <- rep(rep(blBMK, each=numL.ClinBenif*numL.BioResp*numL.stop*numL.dose), 
                       times=numL.tumorType)  
      tumor.in  <- rep(rep(tumorType),
                       each=numL.ClinBenif*numL.BioResp*numL.stop*numL.dose*numL.blBMK)
      numRow <- length(BMK.in)
      if(any(c(length(ClinBenif), length(BioResp), length(stop0), length(dose.in),
               length(tumor.in))!=numRow)){
        stop('Error in level definition: check input value for blBMK, tumorType, dose')
      }
      mat <- cbind(tumorType=tumor.in, blBMK=BMK.in, 
                   dose=dose.in, stop0=stop0, 
                   BioResp=BioResp, ClinBenif=ClinBenif)
      fun1<-function(x, lastCol=NULL){
        if(length(x)<=1) return(x)
        sel <- c(FALSE, x[2:length(x)]==x[1:(length(x)-1)])
        if(!is.null(lastCol)) sel <- sel & (lastCol=="")
        x[sel] <- ""
        return(x)
      }
      
      mat.tab <- mat[mat[,"stop0"]=='no', ]
      m11 <- matrix(fun1(mat.tab[,1]), ncol=1)
      for(i in 2:ncol(mat.tab)){
        m11 <- cbind(m11, fun1(mat.tab[,i], m11[,ncol(m11)]))
      }
      dimnames(m11) <- dimnames(mat.tab)
      mat.tab <- m11
      
      if(is.matrix(mat.tab)) {
        mat.tab[which(mat.tab[,'dose']!=''), 'stop0'] <- 'no'
      }
    }
    
    
    #build arrows' coordinates
    if(showTree){
      
      tot.col <- ncol(mat.tab)
      tot.row <- nrow(mat.tab)
      tot.col2<- tot.col+(tot.col-3)+2 #to enable the fitting lines longer
      text.size1 <- 3/log(tot.row)
      
      th <- th.arrow/log(tot.row)
      
      y.tt <- which(mat.tab[,1]!='')
      pnt <- data.frame(x=c(0, rep(1, length(y.tt))), y=c(1, y.tt), 
                        pch=rep(22,length(y.tt)+1), cex=rep(3,length(y.tt)+1))
      arr <- data.frame(x0=rep(0, length(y.tt))+th*0.8, 
                        y0=rep(1, length(y.tt)), 
                        x1=rep(1, length(y.tt))-th*0.8, 
                        y1=y.tt-th, 
                        lty=rep(1, length(y.tt)),
                        lwd=rep(1, length(y.tt)), 
                        col=rep('gray80', length(y.tt)))  
      y.tt0 <-y.tt;
      dum1 <- 1
      for(i in 2:tot.col){
        y.tt <- which(mat.tab[,i]!='')
        lty1 <- 1;  lwd1 <- 1; col1='gray80';
        if (i>=4) {col1 <- 'black'}
        if (i==5) {lty1<-4; }
        
        if(i <= 3){
          pnt <- rbind(pnt, 
                       data.frame(x=rep(i, length(y.tt)), 
                                  y=y.tt, pch=22, cex=3) )
          
          arr <- rbind(arr, 
                       data.frame(x0=rep(i-1, length(y.tt))+th*0.8, 
                                  y0=rep(which(mat.tab[,i-1]!=''), each=length(y.tt)/length(y.tt0)),
                                  x1=rep(i, length(y.tt))-th*0.8,
                                  y1=y.tt, 
                                  lty=rep(lty1, length(y.tt)),
                                  lwd=rep(lwd1, length(y.tt)),
                                  col=rep(col1, length(y.tt))))
        }else{
          arr <- rbind(arr, 
                       data.frame(x0=rep(i+dum1-2, length(y.tt))+th , 
                                  y0=rep(which(mat.tab[,i-1]!=''), each=length(y.tt)/length(y.tt0)),
                                  x1=rep(i+dum1, length(y.tt))-th,
                                  y1=y.tt, 
                                  lty=rep(lty1, length(y.tt)),
                                  lwd=rep(lwd1, length(y.tt)),
                                  col=rep(col1, length(y.tt))))
          dum1 <- dum1+1
        }
        y.tt0 <-y.tt
      }
      
      
      par(mar=c(0.2, 0.2, 0.2, 0.2), mfrow=c(1,1))
      plot(0~0, col='white', ylim=c(0, tot.row+3), xlim=c(0, tot.col2+ 0.5),
           axes=F, ylab='', xlab='')
      dum2 <- 0;   pos.col<-NULL; 
      th.p <- 0.5; #threshold for probability X-axis position
      #joint probability for Pr(br=1,...), Pr(br=0, ...), Pr(br=1,...), Pr(br=0, ...), etc.
      num.mat.L <- nrow(mat.tab)/2
      #print(num.mat.L)
      j.prob <- rep(1, num.mat.L) 
      for(i in 1:tot.col){
        if(i>3) { dum2 <- dum2+1; th.p <- 1;}
        text(x=i+dum2, y=1:tot.row, labels=mat.tab[,i], cex=text.size1)
        pos.col <- c(pos.col, i+dum2)
        
        if(showProb){
          # add probabilities to the tree and the the expected Utility for each action
          if(i==2){
            text(x=i+dum2-th.p, y=which(mat.tab[,i]!=''), labels=prior_ti, cex=text.size1, col='gray')
            #print(prior_ti)
            j.prob <- j.prob * rep(prior_ti, each=num.mat.L/length(prior_ti))
            #print(j.prob)
          }else if(i==3){
            if(length(prior_prop)==length(dose)) 
              prior_prop <- rep(prior_prop, numL.blBMK*numL.tumorType)
            text(x=i+dum2-th.p, y=which(mat.tab[,i]!=''), labels=prior_prop, cex=text.size1, col='gray')
            j.prob <- j.prob * rep(prior_prop, each=num.mat.L/length(prior_prop))
            #print(j.prob)
          }else if (i==4){
            if(length(prob_stop0)==length(dose))
              prob_stop0 <- rep(prob_stop0, numL.blBMK*numL.tumorType) 
            text(x=i+dum2-th.p, y=which(mat.tab[,i]!=''), labels=prob_stop0, cex=text.size1, col='gray')
            j.prob <- j.prob * rep(prob_stop0, each=num.mat.L/length(prob_stop0))
            #print(j.prob)
          }else if (i==5){
            if(length(prob_BR1)==numL.dose && numL.dose<numL.blBMK*numL.tumorType)
              prob_BR1 <- rep(prob_BR1, numL.blBMK*numL.tumorType)
            text(x=i+dum2-th.p, y=which(mat.tab[,i]=='yes'), labels=prob_BR1, cex=text.size1, col='gray') 
            #the prob are Pr(br=1|...), Pr(br=0|...), Pr(br=1|...), Pr(br=0|...), etc.
            prob_BR <- as.vector(rbind(prob_BR1, 1-prob_BR1))
            j.prob <- j.prob * prob_BR
            #print(j.prob)
          }else if (i==6){
            if(length(prob_CB1_BR)==2)
              prob_CB1_BR <- rep(prob_CB1_BR, numL.blBMK*numL.tumorType)
            #Pr(CB=1|BR=1,...), Pr(CB=1|BR=0,...), Pr(CB=1|BR=1,...), Pr(CB=1|BR=0,...), etc.
            text(x=i+dum2-th.p, y=which(mat.tab[,i]=='yes'), labels=prob_CB1_BR, cex=text.size1, col='gray')          
            #for U(CB=0|....)
            #Pr(CB=0,BR=1|...), Pr(CB=0,BR=0|...),Pr(CB=0,BR=1|...), Pr(CB=0,BR=0|...), etc.
            j.prob0 <- round(payoff[2]*j.prob * trunc01(1-prob_CB1_BR), 3)          
            #for U(CB=1|....)
            #Pr(CB=1,BR=1|...), Pr(CB=1,BR=0|...),Pr(CB=1,BR=1|...), Pr(CB=1,BR=0|...), etc.
            j.prob <- round(payoff[1]*j.prob * prob_CB1_BR, 3)
            
            
            
            #color the top 10%
            num.col <- round(length(j.prob)*topRatio)
            sub.u <- j.prob+j.prob0
            top.p   <- sort(sub.u, decreasing=TRUE)[1:num.col]
            col.p   <- rep('gray', length(j.prob))
            col.p[sub.u%in%top.p] <- topCol
            #print(j.prob)
            text(x=i+dum2+1, y=which(mat.tab[,i]=='yes'), labels=j.prob, 
                 cex=text.size1, col=col.p) 
            text(x=i+dum2+2.5, y=which(mat.tab[,i]=='yes'), labels=j.prob0, 
                 cex=text.size1, col='gray')           
            
            #add the utility into the treat leave
            mat.tab <- cbind(mat.tab, 
                             U=rep("", nrow(mat.tab)),
                             U0=rep("", nrow(mat.tab)),
                             topColor=rep("", nrow(mat.tab)))
            
            #print(prob_CB1_BR)
            #print(j.prob)
            #print(mat.tab)
            
            mat.tab[mat.tab[, i]=='yes', 'U']<-j.prob
            mat.tab[mat.tab[, i]=='yes', 'U0']<-j.prob0
            mat.tab[mat.tab[, i]=='yes', 'topColor']<-col.p
            #get the index for top p values in output table
            wh.top.p<- which(mat.tab[,'topColor']==topCol)
            
          }
          
        }
      }
      points(x=pnt$x, y=pnt$y, pch=pnt$pch, cex=pnt$cex, col='gray80')
      
      arr$col <- as.character(arr$col)
      arrows(x0=arr$x0, y0=arr$y0, x1=arr$x1, y1=arr$y1, length=0.1, lty=arr$lty, lwd=arr$lwd, 
             col=arr$col)
      
      nms <- c('Tumor\nType (T)', 'Baseline\nBiomarker (B)', '\nDose', 
               'Stop due\nto toxicity', 'Biological\nResponse', 
               'Clinical\nBenefit')
      mtext(text=nms, side=3, at=pos.col, padj=1.1)
      
      #It is alway true: P(BR|stop=1)=0
      note1 <- c('P(stop=0|dose)', 'P(BR=1|stop=0, dose, TI)', 'P(CB=1|BR)') 
      text(x=pos.col[-(1:3)]-1, y=tot.row, labels=note1, col='darkgreen', cex=text.size1)
      
      #Assume {dose, stop} is independent from the tumor incident of a biomarker in the tumor type.
      #TI is independent from dose, so P({T,B}|dose)=P({T,B})
      note2 <- c('PriorInfo\nTumorIncidence(TI)', 'Proportion\nDoseLevel')
      text(x=c(pos.col[1]-0.5, pos.col[2]+0.5), y=tot.row, labels=note2, col='blue', cex=text.size1)
      
      if(showProb){
        note3 <- c(paste0('U=', payoff[1], '*Prob(CB=1, BR,\nstop=0, dose, {T,B})'),
                   paste0('U=', payoff[2], '*Prob(CB=0, BR,\nstop=0, dose, {T,B})'))
        pos.note3 <- tot.col2 + c(-1, 0.5)
        text(labels=note3, x=pos.note3, y=rep(tot.row+1, 2), 
             cex=text.size1*0.7, col=c('red', 'magenta'))
        
        #add the expected utility
        mat.tab <- cbind(mat.tab, 
                         U_dTB=rep("", nrow(mat.tab)),
                         U_dTB_color=rep("", nrow(mat.tab)),
                         U_dTB_topColor=rep("", nrow(mat.tab)))
        dum.u <- 0; dum.col <- ''
        for(r in nrow(mat.tab):1){
          if(mat.tab[r,'dose']=="" & mat.tab[r,'U']!=''){
            payoff.pos <-as.numeric(mat.tab[r, 'U'])
            payoff.neg <-as.numeric(mat.tab[r, 'U0'])
            dum.u <- dum.u+ payoff.pos + payoff.neg
            if(mat.tab[r,'topColor']!='') 
              dum.col<-mat.tab[r,'topColor']
          }else if (mat.tab[r,'U']!=''){
            payoff.pos <-as.numeric(mat.tab[r, 'U'])
            payoff.neg <-as.numeric(mat.tab[r, 'U0'])
            dum.u <- dum.u+ payoff.pos + payoff.neg
            mat.tab[r, 'U_dTB'] <- dum.u
            dum.u <- 0
            if(mat.tab[r,'topColor']!='') {
              dum.col<-mat.tab[r,'topColor']
              mat.tab[r, 'U_dTB_color']<-dum.col
            }
          }
        }      
        top.U_dTB <- sort(as.numeric(mat.tab[,"U_dTB"]), decreasing=T)
        top.U_dTB <- top.U_dTB[1:ceiling(length(top.U_dTB)*topRatio)]
        
        u.tDB <- mat.tab[,'U_dTB']; 
        wh.utDB <- which(u.tDB!='')
        mat.tab[wh.utDB, 'U_dTB_topColor']<-'black'
        mat.tab[wh.utDB&u.tDB%in%as.character(top.U_dTB), 'U_dTB_topColor']<-'red'
        u.tDB <- u.tDB[wh.utDB];
        #u.tDB.col<-mat.tab[wh.utDB, 'U_dTB_color']
        u.tDB.col<-mat.tab[wh.utDB, 'U_dTB_topColor']
        note4 <- paste0('U(d,T,B)= ', u.tDB)
        
        if(showBar){
          for(k in 1:length(u.tDB)){
            lines(x=0+c(0, as.numeric(u.tDB[k])), y=rep(wh.utDB[k]-th.utDB, 2), lwd=8, 
                  col=rgb(0, 0, 255, alpha=80, maxColorValue=255) )
          }
          lines(x=c(0,0), y=c(0,nrow(mat.tab)), 
                col=rgb(0, 0, 255, alpha=80, maxColorValue=255))
        }
        
        #u.tDB.col[u.tDB.col=='gray'] <- 'black'
        text(labels=note4, x=3, y=wh.utDB-th.utDB, 
             cex=text.size1, col=u.tDB.col)
        
        
      }
      
    }
    
    
    return(mat.tab)
    
  }
  
  #Obtain combination annotation
  anno_tt_bmk <- function(bmk, tt, other.note=''){
    bmk <- strsplit(as.character(bmk), split=',', fixed=TRUE)[[1]]
    tt <- strsplit(as.character(tt), split=',', fixed=TRUE)[[1]]
    ot <- paste(rep(tt, each=length(bmk)), rep(bmk, length(tt)), sep='_')
    ot <- paste(ot, collapse=',')
    ot <- paste(other.note, ot)
    return(ot)
  }
}
#End 1. -----------------------------------------------------------------------#



#Begin 2. ---------------------------------------------------------------------#
#1. Extension to have the CSF analysis
#2. Generalized to user-defined variables
#3. For discrete variables
if(TRUE){
  #expected loss function for discrete X variables
  #L(theta,a)=0 if x in the range else abs(lev-a)*abs(theta-th)
  #f(theta|data) is a binomial distriubtion
  #Users can define their own expected loss function however the input must
  #be [th, n, p_pos] and the output must be a vector of expected loss under
  #each Bayes decision levels
  my.eLoss <<- function(
    th,             #the vector of thresholds (delta) of decision rule
    n,              #the number of patients in a cohort
    p_pos,          #the posterior probability of responding to drug
    #is the value update by data and affect decsion loss
    d.fun=pbinom    #probability function (lower.tail=T)
  ){
    #defin the distribution as binomial
    #d.fun=pbinom   #the density function of p_pos
    #user can re-define the loss function from here to the end....
    
    th <- sort(as.numeric(th))
    len_a <- length(th)
    
    #~~~data construction~~~#
    #the number of decision levels is the number thresholds plus 1
    xs <- 0:n              #get all numbers in the binomial distribution
    xs_lev <- rep(1, n)    #get decision levels for each number
    for(i in 1:len_a) xs_lev[xs>th[i]]<-i+1
    
    #~~~get the expected loss for each action~~~#
    #a is the action level from 0 to len_a according to theta
    e_loss <- rep(0, len_a+1)
    
    #construct loss elements based on the higher bound for the lowest level
    i<-1
    #els.h<-abs(xs_lev - i)*abs(xs-th[i])*(1-d.fun(xs, size=n, prob=p_pos))
    els.h<-(xs<=th[i])*(xs-th[i])*(d.fun(xs, size=n, prob=p_pos))
    e_loss[i]<-sum(els.h)
    tot.dist <- sum(sapply(n, function(x){
        sum((xs<=x)*abs(xs-x)*(d.fun(xs, size=n, prob=p_pos)))
       }))
    
    #expected loss from level 2 to lev_a
    if(len_a>1){
      for(i in 2:len_a){#if taking the action as level i
        els.h<-abs(xs_lev - i)*abs(xs-th[i])*(1-d.fun(xs, size=n, prob=p_pos))
        els.l<-abs(xs_lev - i)*abs(xs-th[i-1])*(1-d.fun(xs, size=n, prob=p_pos))
        e_loss[i] <- sum(els.h[xs_lev>i]) + sum(els.l[xs_lev<i])
      }
    } 
    
    #for highest level
    #construct loss elements based on the lower bound
    i <- len_a
    els.l<-(xs>th[i])*abs(xs-th[i])*(1-d.fun(xs, size=n, prob=p_pos))
    e_loss[i+1]<-sum(els.l)
    tot.dist2 <- sum(sapply(0, function(x){
        sum((xs>x)*abs(xs-x)*(1-d.fun(xs, size=n, prob=p_pos)))
       }))
    
    #return(e_loss/c(tot.dist, tot.dist2))
    return(e_loss/n)
  }
  
  
  #improved function for Baysian decision theory with Critical Success Factor
  #available to add or remove variables
  #available to specify the Bayse loss function
  #available to select utility bar and the loss bar
  #Clinical Benifit and Biological Response are {0, 1} for {no or yes}. 
  #Output@: ranked table {dose, Disease, baselineBiomarker, expectedUtility}
  #Output@: tree plot
  #function name: Bayesian Decision Theory Utility and Loss
  #fixing the expected utility and make loss function flexible
  BDT_UaL <- function(
    levVars="B1::T1::coh1::BR1,B1::T1::coh2::BR1,B1::T2::coh3::BR2",       
    #variables and levels separated by "::"
    dr_lev="nogo::go,nogo::go,nogo::go",  #order does matter.
    #decision rule labels
    incidence="0.3,0.3,0.1",      #Biomarker incidence in the tumor type 
    #values obtained from prior knowledge
    pBprior=NULL,                 #the hyper parameters "alph, beta"
    #if NULL, then alpha=1+incidence
    #beta=2-incidence
    #if not NULL, the value should be like
    #"1.3 1.7,1.3 1.7, 1.1 1.9," 
    n_ij="10, 10, 10",            #sample sizes for each cohort
    #values obtained from decision makers
    dr_th="0.1,0.5,0.6",          #decision rule threshold (delta)
    #if 3 levels of decision rule such as 
    #dr_lev="go::moreData::nogo," then
    #dr_th="0.9::0.3,"
    drFunc=my.eLoss,              #user-defined Bayes decision loss function
    #input: [th, n, p_pos]
    #output: a vector of Bayes decision loss
    showBar=TRUE,                 #show the barplot of utility & loss
    th.arrow= 0.8,                #horizontal space between an arrow and target
    payoff= c(10, -1)            #payoff value for utility, only two values
    #gain vs lost
  ){
    #an internal function truncates values to 0, 1
    trunc01 <- function(val){
      val[val>1]<-1
      val[val<-0]<-0
      return(val)
    }
    
    #an internal function gets the hierarchical variables
    my.split1 <- function(mylab="", s1=",", s2="::"){
      if(length(mylab)==1 & is.character(mylab)){
        L1 <- strsplit(mylab, split=s1)[[1]]
      }else{L1 <- mylab}
      if(length(L1)==0) return('L1 in my.split1 is missing.')
      if(is.null(s2)) return(L1)
      if(all(is.character(L1))){
        L2 <- strsplit(L1, split=s2)
      }else{
        L2 <- mylab
      }
      return(L2)
      #mat1<-t(matrix(unlist(L2), ncol=length(L2)))
    }
    
    
    #cleanup the input parameters
    if(TRUE){
      #sample size proportions
      if(is.null(n_ij) || all(n_ij=='')|all(incidence=='')|
         all(dr_th=='')) return(NULL)
      n_1 <- as.numeric(unlist(my.split1(n_ij)))
      n_rt<- n_1/sum(n_1)
      #incidences
      if(is.null(incidence)) return(NULL)
      incd<- as.numeric(unlist(my.split1(incidence)))
      
      p.0<-list()
      for(o in 1:length(incd)){
        #p.0[[o]]<-c(1+incd[o], 2-incd[o])
        p.0[[o]]<-c(1, 2)
      }
      
      if(!is.null(pBprior) && pBprior!="~" && 
         gsub(" ", "", pBprior)!=""){
        p.0p<-my.split1(pBprior, s2=" ")
        p.0p<-lapply(p.0p, function(x){x[x!=""]})
        for(o in 1:length(p.0)){
          p.0[[o]] <- p.0[[o]]+as.numeric(p.0p[[o]])*incd[o]
        }
        #note the length p.0 == the leve of plans
      }
      num.p0 <- length(p.0)
      
      #if using default threshold
      if(dr_th=="~"){
        th1 <- list()
        for(o in 1:length(p.0))
          th1[[o]] <- p.0[[o]][1]/sum(p.0[[o]])
      }else{
        th1 <- lapply(my.split1(dr_th),as.numeric)
      }
      #print(th1)
      
      #cleanup the payoff values
      if(is.character(payoff)){
        payoff <- as.numeric(strsplit(payoff, split=',', fixed=TRUE)[[1]])
        if(length(payoff)==1){
          payoff<-c(payoff, 0)
        }else{
          payoff <- payoff[1:2]
        }
      }
    }
    
    #construct the decision tree with user-defined variables
    if(TRUE){
      
      if(is.null(levVars)) return(NULL)
      
      LV0 <- LV <- my.split1(levVars)
      num.var <- length(LV[[1]])
      drLV  <- my.split1(dr_lev)
      num.dr  <- length(drLV[[1]])
      
      if(length(LV)!=length(drLV) | length(LV)!=num.p0){
        #print('Error: lengths of decision rule and layers do not match!')
        return(NULL)
      }
      
      iLV <- E.L <- U <- p.1 <- list()
      for(o in 1:length(LV)){
        th2<-round(th1[[o]]*n_1[o])
        p.1[[o]] <- p.0[[o]][1]/sum(p.0[[o]])
        eL <- drFunc(th=th2,
                     n=n_1[o],
                     p_pos=p.1[[o]] )
        E.L[[o]]<-eL  #expected loss
        
        U[[o]] <- incd[o]*n_rt[o]*p.1[[o]]*payoff[1]+
          incd[o]*n_rt[o]*(1-p.1[[o]])*payoff[2]
        #expected utility
        
        if(o>1){
          o.wh <- which(LV0[[o]]!=LV0[[o-1]])[1]
          if(length(o.wh)==0) o.wh<-1
          no.wh <- which(LV0[[o]]==LV0[[o-1]])
          LV[[o]][ no.wh[no.wh<o.wh] ]<-''
        }
        
        iLV[[o]] <- c(LV[[o]][1:(num.var-1)],
                      paste0(LV[[o]][num.var], ", n=", n_1[o],
                             "\nI=", incd[o],
                             ", U=", round(U[[o]],3),
                             ", p=", round(p.1[[o]],3)  ), 
                      paste0(drLV[[o]][1],": go if r>", th2[1],
                             ", E(L)=", round(eL[1],3)))
        if(num.dr==1) next
        for(h in 2:num.dr){
          if(h==num.dr){
            iLV[[o]] <- c(iLV[[o]], rep('', num.var), 
                          paste0(drLV[[o]][h], ": stop if r<=",th2[h-1],
                                 ", E(L)=", round(eL[h],3))  )
          }else{
            iLV[[o]] <- c(iLV[[o]], rep('', num.var), 
                          paste0(drLV[[o]][h], ": ",th2[h-1],
                                 "<= r <",th2[h],
                                 ", E(L)=", round(eL[h],3))  )
          }
        }
      }
      varMat <- t(matrix(unlist(iLV), nrow=num.var+1))
      #E.L[[o]]: expected Bayes decision loss
      #U[[o]]: utility of plan
    }
    
    #build arrows and coordinates to show the decision tree
    if(TRUE){
      tot.col   <- ncol(varMat)
      max.nchar <- apply(varMat, 2, function(x){max(nchar(x))})
      cex.1char<- 0.1
      tot.row <- nrow(varMat)
      tot.col2<- tot.col+2 #to enable the fitting lines longer
      text.size1 <- 3/log(tot.row)
      th.a <- th.arrow/log(tot.row) #about arrow locaiton
      
      #for Layer 1
      lty1 <- 1;  lwd1 <- 1; col1='gray80';
      y.tt <- which(varMat[,1]!='')
      pnt <- data.frame(x=c(0, rep(1, length(y.tt))), y=c(1, y.tt), 
                        lab=c('', varMat[y.tt,1]),
                        pch=rep(22,length(y.tt)+1), cex=rep(3,length(y.tt)+1))
      arr <- data.frame(x0=rep(0, length(y.tt)), 
                        y0=rep(1, length(y.tt)), 
                        x1=rep(1, length(y.tt)), 
                        y1=y.tt, 
                        lty=rep(1, length(y.tt)),
                        lwd=rep(1, length(y.tt)), 
                        col=rep(col1, length(y.tt)))  
      shf <- sum(max.nchar[1])*cex.1char
      y.tt0 <- y.tt
      for(i in 2:tot.col){
        y.tt <- which(varMat[,i]!='')
        pnt <- rbind(pnt, data.frame(x=rep(i+shf, length(y.tt)), 
                                     y=y.tt, lab=varMat[y.tt,i],
                                     pch=22, cex=3) )
        wh.a1<-which(!y.tt%in%y.tt0)
        y.tt0a <- y.tt
        for(a in wh.a1){
          y.tt0a[a] <- y.tt0a[a-1]
        }
        arr <- rbind(arr,
                     data.frame(x0=rep(i-1+shf, length(y.tt)), 
                                y0=y.tt0a,
                                x1=rep(i+shf, length(y.tt)),
                                y1=y.tt, 
                                lty=rep(lty1, length(y.tt)),
                                lwd=rep(lwd1, length(y.tt)),
                                col=rep(col1, length(y.tt))))
        shf <- sum(max.nchar[1:i])*cex.1char
        y.tt0 <- y.tt
      }
      
      
      par(mar=c(0.2, 0.2, 0.2, 0.2), mfrow=c(1,1))
      plot(y~x, data=pnt, col='gray80', pch=pnt$pch,
           ylim=c(0, tot.row), xlim=c(0, tot.col2+shf+0.5),
           axes=F, ylab='', xlab='')
      text(x=pnt$x, y=pnt$y, labels=pnt$lab, adj=-0.07)
      arrows(x0=arr$x0, y0=arr$y0, x1=arr$x1, y1=arr$y1, 
             length=0.1, lty=arr$lty, lwd=arr$lwd, 
             col=arr$col)
    }
    
    #add barplot of utility and loss
    if(showBar){
      u.x0<-rep(0, length(LV))
      u.y0<-which(varMat[,ncol(varMat)-1]!='')
      u.x1<-unlist(U)
      col.bar1 <- rgb(0, 0, 255, alpha=80, maxColorValue=255)
      abline(v=0, col=col.bar1)
      for(i in 1:length(u.y0)){
        lines(x=c(u.x0[i], u.x1[i]), y=c(u.y0[i], u.y0[i]),
              lwd=8, 
              col=col.bar1)
      }
      
      l.x0<-rep(tot.col2+shf, nrow(varMat))
      l.y0<-1:nrow(varMat)
      l.x1<-l.x0-unlist(E.L)
      col.bar2 <- rgb(255, 0, 0, alpha=80, maxColorValue=255)
      abline(v=tot.col2+shf, col=col.bar2)
      for(i in 1:length(l.y0)){
        lines(x=c(l.x0[i], l.x1[i]), y=c(l.y0[i], l.y0[i]),
              lwd=8, 
              col=col.bar2)
      }
    }
    
    return(list(dat=varMat, BayesLoss=E.L, U=U, p=p.1))
    
  }
  
}
#define global variables for FDA_log analysis
if(TRUE){
  #get a subset
  key.words <<- c('all', 'NSCLC|lung', 
                  'urothelial', 
                  'gastrointestinal', 
                  'msi', 'breast', 'head', 'hcc',
                  'other')
  mySelect<-function(vec1, sp=NULL,
                     ctype='all',
                     fdaLink=FALSE, 
                     allkeys=key.words  ){
    
    if(!is.null(fdaLink)&&as.logical(fdaLink)){
      vec1 <- readLines(fL)
      #vec1 <- vec1[grepl('<li>', vec1)]
      vec1 <- vec1[grepl('approv', vec1)]
      if(length(vec1)==1){
        vec1<-vec1[grepl('href=\"#updates\"', vec1, fixed=T)]
        vec1<-strsplit(split='<li>', vec1, fixed=TRUE)[[1]]
        vec1<-gsub('\t|  ','', vec1)
        vec1 <- vec1[grepl('approv', vec1)]
      }
      if(length(vec1)==0){
        return(data.frame(FDA_log='Fail to reach the link.'))
      }
    }
    if(is.null(sp) || is.null(ctype)){
      return(data.frame(FDA_log='no data'))
    }
    if(!ctype%in%c('all')){
      if(ctype=='other'){
        ctype <- paste(allkeys[!allkeys%in%c('all','other')], collapse='|')
        vec1<-vec1[!grepl(ctype, vec1, ignore.case=T)]
      }else{
        vec1<-vec1[grepl(ctype, vec1, ignore.case=T)]
      }
      if(length(vec1)==0){
        return(data.frame(FDA_log='The disease is not found.'))
      }
    }
    if(sp==''|sp=='~'){
      if(is.data.frame(vec1)){
        dat1 <- data.frame(FDA_log=vec1[,1])
      }else if (is.vector(vec1)){
        dat1<- data.frame(FDA_log=vec1)
      }else{
        dat1<-data.frame(FDA_log='no data')
      }
      return(dat1)
    }
    vec2<- vec1[grepl(
      paste(paste0("(?=.*",strsplit(as.character(sp), split='&')[[1]], ")"), 
            collapse=""), vec1, perl=T, ignore.case=T)]
    dat1<-data.frame(FDA_log=vec2)
    return(dat1)
  }
}
#End 2. -----------------------------------------------------------------------#



#Begin 3. ---------------------------------------------------------------------#
#1.add contineus variables
#2.add benchmark reference
if(TRUE){
  
  #A function for P(y.trt-y.ref<x) given the equal length samples of y.trt, y.ref
  prob.diff <- function(x, y.diff=NULL, y.trt=NULL, y.ref=NULL){
    if(!is.null(y.diff)){
      dif1 <- y.diff[is.finite(y.diff)]
      return( mean(dif1<x, na.rm=T) )
    }else if(!is.null(y.trt) & !is.null(y.ref) & length(y.trt)==length(y.ref)){
      dif1 <- (y.trt-y.ref)
      dif1 <- dif1[is.finite(dif1)]
      return( mean(dif1<x, na.rm=T) )
    }else{return(0)}
  }
  #functions getting the random samples from given distribution
  #rbinom with n and pi
  #rlnorm with meanlog and sdlog
  
  
  #expected loss function for either discrete or continues X variables
  #L(theta,a)=0 if x in the range else abs(lev-a)*abs(theta-th)
  #f(theta|data) is a binomial distriubtion
  #Users can define their own expected loss function however the input must
  #be [th, n, p_pos] and the output must be a vector of expected loss under
  #each Bayes decision levels
  eLoss.diff <<- function(
    th=NULL,              #the vector of thresholds (delta) of decision rule
    sample1,             #a vector of samples under treatment assumption
    sample2,             #a vector of samples under control assumption
    sample1.prob,        #a vector of samples prob under treatment assumption
    sample2.prob,        #a vector of samples prob under control assumption
    len_cutoffs=NULL #the number of cutoffs or thresholds of a decision rule
  ){
    prob.1g2<-NULL
    #create the matrix of difference and the probablity
    if(!is.vector(sample1) | !is.vector(sample1.prob) |
       length(sample1)!=length(sample1.prob)){
      stop('Sample1 input is wrong for eLoss.diff')
    }else if(!is.vector(sample2) | !is.vector(sample2.prob) |
             length(sample2)!=length(sample2.prob)){
      stop('Sample2 input is wrong for eLoss.diff')
    }else{
      n.r <- length(sample1)
      n.c <- length(sample2)
      sam1.mat<-matrix(sample1, nrow=n.r, ncol=n.c)
      sam1prob.mat<-matrix(sample1.prob, nrow=n.r, ncol=n.c)
      sam2.mat<-t(matrix(sample2, nrow=n.c, ncol=n.r))
      sam2prob.mat<-t(matrix(sample2.prob, nrow=n.c, ncol=n.r))
      prob.mat <- sam1prob.mat*sam2prob.mat
      diff.mat <- sam1.mat - sam2.mat
      y2sampleD <- as.vector(diff.mat)
      ord <- order(y2sampleD)
      y2sampleD <- y2sampleD[ord]
      y2sampleD.prob <- as.vector(prob.mat)[ord]
      prob.1g2 <- sum(y2sampleD.prob[y2sampleD>0])
    }
    
    
    #Recomend threshold maximizing the loss difference
    if(is.null(th)||th%in%c("", "~", "-")){
      if(is.null(len_cutoffs)){len_cutoffs <- 1}
      th.prob <- seq(0,1, by=1/(len_cutoffs+1))
      th.prob <- th.prob[c(-1, -length(th.prob))]
      th <- quantile(y2sampleD, th.prob)
      recom.th <- TRUE
    }else{recom.th<-FALSE}
    
    
    #order the thereshold for a decision rule
    th <- sort(as.numeric(th))
    len_a<- length(th)
    n  <- length(y2sampleD)
    
    #~~~data construction~~~#
    #the number of decision levels is the number thresholds plus 1
    xs <- y2sampleD  #get the distribution consits of samples
    xs_lev <- rep(1, n)    #get decision levels for each sample
    for(i in 1:len_a) {
      xs_lev[xs>th[i]]<-i+1
    }
    
    #~~~get the expected loss for each action~~~#
    #a is the action level from 0 to len_a according to theta
    e_loss <- rep(0, len_a+1)
    
    #construct loss elements based on the higher bound for the lowest level
    #i<-1
    els.h<-abs(xs-th[1])*y2sampleD.prob
    e_loss[1]<-sum(els.h[xs_lev==1])
    
    #expected loss from level 2 to lev_a
    i<-2
    while(i>=2 & i<=len_a){#if taking the action as level i
      els.h<-abs(xs-th[i])*y2sampleD.prob
      els.l<-abs(xs-th[i-1])*y2sampleD.prob
      els <- els.h+els.h
      e_loss[i] <- sum(els[xs_lev>=i-1 & xs_lev<i])
      i<-i+1
    }
    
    #for highest level
    #construct loss elements based on the lower bound
    els.l<-abs(xs-th[len_a])*y2sampleD.prob
    e_loss[len_a+1]<-sum(els.l[xs_lev>len_a])
    
    if(recom.th){
      return(list(e_loss=e_loss, th=th, p.1g2=prob.1g2))
    }else{
      return(list(e_loss=e_loss, p.1g2=prob.1g2))  
    }
    
  }
  
  #A revised function of generating random variable from a lognormal distb
  rlnorm2 <- function(B=1000000, m, s){
    location <- log(m^2 / sqrt(s^2 + m^2))
    shape <- sqrt(log(1 + (s^2 / m^2)))
    smp <- rlnorm(n=B, location, shape)
    return(smp)
  }
  
  #improved function for Baysian decision theory with Critical Success Factor
  #take the difference between treatment variable and the benchmark control
  #available to add or remove variables
  #available to specify the Bayse loss function
  #available to select utility bar and the loss bar
  #Clinical Benifit and Biological Response are {0, 1} for {no or yes}. 
  #Output@: ranked table {dose, Disease, baselineBiomarker, expectedUtility}
  #Output@: tree plot
  #function name: Bayesian Decision Theory Utility and Loss
  #fixing the expected utility and make loss function flexible
  BDT_UaL.diff <- function(
    levVars="B1::T1::coh1::BR1,B1::T1::coh2::BR1,B1::T2::coh3::BR2",       
    #variables and levels separated by "::"
    dr_lev="No Go::Go,No Go::Go,No Go::Go",  #order does matter.
    #decision rule labels
    incidence="0.3,0.3,0.1",      #Biomarker incidence in the tumor type 
    #values obtained from prior knowledge
    numRsp=NULL,                 #the hyper parameters "alph, 2-alpha"
    #if NULL, then alpha=1+incidence and beta=2-incidence, 
    #and only trt ORR is onsidered without the benchmark control ORR
    #if not NULL, the value should be like
    #"1.3vs1.3,1.3vs1.2, 1.1vs1" matches "TRT_ORRvsCTR_ORR"
    muTTE=NULL,
    sdTTE=NULL,           #the mean and sd of TTE variables
    #if not NULL, e.g. "5vs5, 9vs2, 10vs8" and "1vs1, 1vs1, 1vs1"
    #matches "TRTvsCTR"
    n_ij="10, 10, 10",            #sample sizes for each cohort
    #or "10vs20,10vs80,10vs10" for "TRTvsCTR"
    #values obtained from decision makers
    dr_th="0.1,0.5,0.6",          #decision rule threshold for ORR
    #if 3 levels of decision rule such as 
    #dr_lev="go::moreData::nogo," then
    #dr_th="0.9::0.3,"
    #input: [th, n, p_pos]
    #output: a vector of Bayes decision loss
    showPlot=TRUE,  showBar=TRUE,      #show the barplot of utility & loss
    onebar=TRUE,                 #only show integrated bar at left
    th.arrow= 0.8,               #horizontal space between an arrow and target
    payoff= c(10, -1),            #payoff value for utility, only two values
    #gain vs lost
    bar.wid=8,          #the width of the horizontal bar
    Bsample=10000        #number of bootstrapping samples
  ){
    set.seed(1234567)
    #a list of difference in ORR between trt and ctr
    RT.diff <- list()
    #a list of difference in tte between trt and ctr
    tte.diff<- list()
    #by default, the difference variable between trt and ctr is not used.
    #it will be TRUE if 'vs' shows in the string of numRsp
    useRspDiff<-FALSE
    #it will be TRUE if 'vs's shows in the string muTTE
    useTteDiff<-FALSE
    
    #an internal function gets the hierarchical variables 
    #seperated by "::" for different levels of a variable
    #seperated by "," for different subgroups
    #seperated by "vs" for treatment vs control
    my.split1 <- function(mylab="", s1=",", s2="::"){
      if(length(mylab)==1 & is.character(mylab)){
        #the input is a string
        L1 <- strsplit(mylab, split=s1)[[1]]
      }else{ #when the input is already a vector
        L1 <- mylab
      }
      if(length(L1)==0) return('L1 in my.split1 is missing.')
      if(is.null(s2)) return(L1)
      if(all(is.character(L1))){
        L2 <- strsplit(L1, split=s2)
      }else{
        L2 <- mylab
      }
      return(L2)
    }
    
    
    #cleanup the input parameters
    if(TRUE){
      #Cleanup sample size parameters
      if(is.null(n_ij) || all(n_ij=='')|all(incidence=='')) return(NULL)
      if(grepl('vs',n_ij)){
        nL <- my.split1(n_ij, s2='vs')
        n_1.trt<- as.numeric(as.vector(sapply(nL,function(x){x[1]})))
        n_1.ctr<- as.numeric(as.vector(sapply(nL,function(x){x[2]})))
        n_1 <- n_1.trt #n_1 alwasy has values for treatment sizes
        if(length(n_1.trt)!=length(n_1.ctr)) return(NULL)
      } else {
        n_1 <- as.numeric(unlist(my.split1(n_ij)))
        n_1.trt <- n_1.ctr <- n_1
      }
      n_rt<- n_1/sum(n_1)
      
      #incidences
      if(is.null(incidence)) return(NULL)
      incd<- as.numeric(unlist(my.split1(incidence)))
      
      #for treatment prior parameters of Beta Variable pi based on incidences
      p.0<-list()
      for(o in 1:length(incd)){
        #p.0[[o]]<-c(1+incd[o], 2-incd[o])
        p.0[[o]]<-c(1, 2)
      }
      
      #posterior parameters of Beta Variable pi updated by trt ORR
      #and get RT.diff
      if(!is.null(numRsp) && numRsp!="~" && gsub(" ", "", numRsp)!=""){
        if(!(grepl("vs", numRsp)) && !is.null(n_1)){
          #when no sample sizes for control or benchmark
          #only the Beta posterrior parameters are updated
          p.0p<-as.numeric(my.split1(numRsp, s2=NULL))
          if(length(incd)!=length(p.0p)) return(NULL)
          if(length(p.0p)!=length(n_1)){
            #use the probability only without the sample sizes
            #if(any(p.0p>1|p.0p<0)) return(NULL)
            for(o in 1:length(p.0)){
              if(p.0p[o]>1){
                p.0[[o]]<-p.0[[o]]+c(p.0p[o]*incd[o], n_1[o]-p.0p[o]*incd[o])
              }else{
                p.0[[o]]<-p.0[[o]]+c(p.0p[o]*incd[o], 1-p.0p[o]*incd[o])
              }
            }
          } else {
            for(o in 1:length(p.0)){
              if(p.0p[o]>1){
                p.0[[o]]<-p.0[[o]]+c(p.0p[o]*incd[o], n_1[o]-p.0p[o]*incd[o])
              }else{
                p.0[[o]]<-p.0[[o]]+c(p.0p[o]*incd[o], 1-p.0p[o]*incd[o])
              }
            }
          }
          #note: p.0 is a listing object. Each element is a vector of
          #the two Beta parameters. The total number of elements matches
          #the number of incidences or cohorts/arms.
        }else{
          #when number of responders in control are inserted in numRsp
          #then the difference between trt and ref is given
          p.0p <- my.split1(numRsp, s2='vs')
          p.0p.trt<-as.numeric(sapply(p.0p, function(x){x[1]}))
          p.0p.ctr<-as.numeric(sapply(p.0p, function(x){x[2]}))
          if(is.null(n_1.ctr)) n_1.ctr <- n_1.trt
          if(!is.null(n_1.ctr) & length(p.0p)==length(n_1.ctr)){
            for(o in 1:length(p.0p)){
              if(p.0p.trt[o]>1){#if the new alpha>1
                p.0[[o]] <- p.0[[o]]+c(p.0p.trt[o]*incd[o], n_1[o]-p.0p.trt[o]*incd[o])
              }else{
                p.0[[o]] <- p.0[[o]]+n_1[o]*c(p.0p.trt[o]*incd[o], 1-p.0p.trt[o]*incd[o])
                #p.0p.trt[o] <- round(p.0p.trt[o], n_1[o])
                p.0p.trt[o] <- p.0p.trt[o]
              }
              if(p.0p.ctr[o]>1){
                prob.ctr<-p.0p.ctr[o]/n_1.ctr[o]
              }else{
                prob.ctr<-p.0p.ctr[o];
                p.0p.ctr[o]<-p.0p.ctr[o]*n_1.ctr[o]
              }
              sam.trt <- sam.ctr <- 0:n_1.trt[o]
              use1 <- 1
              if(use1==1){ #######2020-10-29   use emprical binomial
                sam.trt.prob<-dbinom(sam.trt, size=n_1.trt[o], 
                                   prob=p.0[[o]][1]/sum(p.0[[o]]) )
                sam.ctr.prob<-dbinom(sam.ctr, size=n_1.ctr[o], 
                                   prob=prob.ctr )
              }else if (use1==2){ #use preditive distribution for both trt and cntr
                sam.trt.prob<-dbbinom(sam.trt, size=n_1.trt[o], 
                                  alpha=1+p.0p.trt[o], beta=2-p.0p.trt[o] )
                sam.ctr.prob<-dbbinom(sam.ctr, size=n_1.ctr[o], 
                                  alpha=1+prob.ctr, beta=2-prob.ctr)
              }else{
                sam.trt.prob<-dbeta(sam.trt/n_1[o],  
                                   1+p.0p.trt[o], 2-p.0p.trt[o] )
                sam.ctr.prob<-dbeta(sam.ctr/n_1[o],  
                                   1+prob.ctr, 2-prob.ctr)
              }
              RT.diff[[o]] <- list()
              RT.diff[[o]]$sample1 <- sam.trt
              RT.diff[[o]]$sample2 <- sam.ctr
              RT.diff[[o]]$sample1.prob <- sam.trt.prob
              RT.diff[[o]]$sample2.prob <- sam.ctr.prob
              #print(paste('sum(sam.trt.prob)', sum(sam.trt.prob))); 
              #print(paste(p.0[[o]][1]/sum(p.0[[o]]), ',', prob.ctr))
              #print(quantile(sam.trt.prob))
              #print(quantile(sam.trt))
              #print(quantile(sam.ctr))
              #print(paste('sum(sam.ctr.prob)', sum(sam.ctr.prob)))
            }
            useRspDiff<-TRUE
          }else{return(NULL)}
        }
        #note: RT.diff[[o]] is a listing object for the o_th cohorts/arms. 
        #Four objects are saved in the listting RT.diff[[o]]: 
        #sample1, sample1.prob for treatment and  
        #sample2, sample2.prob for control
      }
      
      #Note:the length p.0 == the level of plans or #cohorts
      num.p0 <- length(p.0)
      
      #When using the default threshold
      if(!is.null(dr_th) && dr_th%in%c("~", "")){
        th1 <- list()
        #use prevalence for the thresholds
        for(o in 1:length(p.0))
          th1[[o]] <- p.0[[o]][1]/sum(p.0[[o]])
      }else{
        #when the thresholds are defined in a string such as "3::5,4::6"
        #then 3 and 4 are the lower cutoff for cohort 1 and 2, 
        #5 and 6 are the higher cutoff for cohort 1 and 2. 
        th1 <- lapply(my.split1(dr_th),as.numeric)
      }
      #note: th1 is a listing object. Each element is a vector of ORR 
      #cutoffs. The number of elements matches the number of incidences.
      
      #cleanup the two parameters for lognormal variables, such as TTE
      #Note: The benchmark control reference must be provided for the analysis
      tte.mu <- tte.sd <- list()
      if(!is.null(muTTE) & !is.null(sdTTE) & 
         grepl('vs', muTTE) & grepl('vs', sdTTE)){
        #a listing object of lognormal means. each element is the two means
        #of trt vs cntr
        tte.mu<-lapply(my.split1(muTTE, s2='vs'), as.numeric)
        tte.sd<-lapply(my.split1(sdTTE, s2='vs'), as.numeric)
        if(length(tte.mu)==length(tte.sd) & 
           all(sapply(tte.mu, length)==2) & 
           all(sapply(tte.sd, length)==2) &
           length(tte.mu)==length(incd)
        ){
          for(o in 1:length(incd)){
            cy.trt <- rlnorm2(B=Bsample, m=tte.mu[[o]][1], s=tte.sd[[o]][1])
            cy.ctr <- rlnorm2(B=Bsample, m=tte.mu[[o]][2], s=tte.sd[[o]][2])
            tte.diff[[o]] <- cy.trt-cy.ctr
            tte.diff[[o]] <- tte.diff[[o]][is.finite(tte.diff[[o]]) & 
                                             !is.na(tte.diff[[o]])]
          }
          useTteDiff<-TRUE
        }
      }
      
      #cleanup the payoff values
      if(is.character(payoff)){
        payoff <- as.numeric(strsplit(payoff, split=',', fixed=TRUE)[[1]])
        if(length(payoff)==1){
          payoff<-c(payoff, 0)
        }else{
          payoff <- payoff[1:2]
        }
      }
    }
    
    #construct the decision tree with user-defined variables
    #calculate Utility and Loss
    if(TRUE){
      
      #A listing object of plans. 
      #Each element is a vector of critical variables.
      if(is.null(levVars)) return(NULL)
      LV0 <- LV <- my.split1(levVars)
      num.var <- length(LV[[1]]) #number of variables
      
      #A listing object of decision rules. 
      #Each element is a vector of decision thresholds
      drLV  <- my.split1(dr_lev) 
      num.dr  <- length(drLV[[1]])
      
      if(length(LV)!=length(drLV) | length(LV)!=num.p0){
        #print('Error: lengths of decision rule and layers do not match!')
        return(NULL)
      }
      
      #lost, utility, response probability pi. 
      #p.2 is the posterior prob for benefit.
      iLV <- E.L <- U <- Utte <- p.1 <-p.2 <- list()
      for(o in 1:length(LV)){
        #convert the ORR_trt decision thresholds into numRsp_trt
        th2.lab<-th2<-round(th1[[o]]*n_1[o])
        
        #eLoss.diff(th, y2sampleD)
        if(useRspDiff){
          #if benchmark ref is available then the threshold th2 change to be
          #the difference. 
          th2 <- round((th1[[o]]-p.0p.ctr[o]/n_1.ctr[o])*n_1[o])
          #print(th2); print(paste("length(RT.diff[[o]])", length(RT.diff[[o]])));
          
          eL <- eLoss.diff(th=th2, 
                           sample1=RT.diff[[o]]$sample1,
                           sample2=RT.diff[[o]]$sample2,
                           sample1.prob=RT.diff[[o]]$sample1.prob,
                           sample2.prob=RT.diff[[o]]$sample2.prob)
          #normalize the loss
#print('sum(eL$e_loss)'); print(sum(eL$e_loss))
          eL$e_loss <- eL$e_loss/sum(eL$e_loss)
          #probability of the difference pi_trt-pi_ctr>0
          #print(paste0('p.1[[',o,']]=', p.1[[o]] <- eL$p.1g2))
          p.1[[o]] <- eL$p.1g2
          eL <- eL$e_loss
        }else{
          #expected response probability pi_trt
          p.1[[o]] <- p.0[[o]][1]/sum(p.0[[o]])
          eL <- my.eLoss(th=th2, n=n_1[o], p_pos=p.1[[o]] )
          #normalized the loss
          #eL <- eL/sum(eL)
          #convert the expected distance to success in the scale of response rate
          #eL <- eL/n_1[o]
        }
        E.L[[o]]<-eL  #expected loss
        
        #expected utility based on ORR
        # U[[o]] <- incd[o]*n_rt[o]*p.1[[o]]*payoff[1]+
        #   incd[o]*n_rt[o]*(1-p.1[[o]])*payoff[2]
        U[[o]] <- incd[o]*p.1[[o]]*payoff[1] + incd[o]*(1-p.1[[o]])*payoff[2]
        
        #expected utility based on TTE
        if(useTteDiff){
          #P(mu_trt-mu_ctr>0|...)
          p.2[[o]] <- mean(tte.diff[[o]]>0)
          # Utte[[o]] <- incd[o]*n_rt[o]*p.2[[o]]*payoff[1]+
          #   incd[o]*n_rt[o]*(1-p.2[[o]])*payoff[2]
          Utte[[o]] <- incd[o]*p.2[[o]]*payoff[1]+
            incd[o]*(1-p.2[[o]])*payoff[2]
        }
        
        if(o>1){
          o.wh <- which(LV0[[o]]!=LV0[[o-1]])[1]
          if(length(o.wh)==0) o.wh<-1
          no.wh <- which(LV0[[o]]==LV0[[o-1]])
          LV[[o]][ no.wh[no.wh<o.wh] ]<-''
        }
        
        iLV[[o]] <- c(LV[[o]][1:(num.var-1)],
                      paste0(LV[[o]][num.var], ", n=", n_1[o],
                             "\nI=", incd[o],
                             ", U_orr=", round(U[[o]],3),
                             ifelse(useRspDiff,
                                    ", p(trt>ref)=", 
                                    ", p="), round(p.1[[o]],3),
                             ifelse(useTteDiff, 
                                    paste0("\nU_tte=",round(Utte[[o]],3),
                                           ", m_trt=", tte.mu[[o]][1],
                                           ", m_ref=", tte.mu[[o]][2],
                                           ", s_trt=", tte.sd[[o]][1],
                                           ", s_ref=", tte.sd[[o]][2]), 
                                    "")), 
                      paste0(drLV[[o]][1],
                             ifelse(useRspDiff, " if d_r<=", " if r<="),
                             th2.lab[1],
                             ", E(L)=", round(eL[1],3)))
        if(num.dr==1) next
        for(h in 2:num.dr){
          if(h==num.dr){
            iLV[[o]] <- c(iLV[[o]], rep('', num.var), 
                          paste0(drLV[[o]][h], 
                                 ifelse(useRspDiff," if d_r>", " if r>"),
                                 th2.lab[h-1],
                                 ", E(L)=", round(eL[h],3))  )
          }else{
            iLV[[o]] <- c(iLV[[o]], rep('', num.var), 
                          paste0(drLV[[o]][h], ": ",th2.lab[h-1],
                                 ifelse(useRspDiff,"<= d_r <" , "<= r <"),
                                 th2.lab[h],
                                 ", E(L)=", round(eL[h],3))  )
          }
        }
      }
      varMat <- t(matrix(unlist(iLV), nrow=num.var+1))
      #E.L[[o]]: expected Bayes decision loss
      #U[[o]]: utility of plan based on ORR
      #Utte[[o]]: utility based on tte
    }
    
    #build arrows and coordinates to show the decision tree
    if(TRUE){
      tot.col   <- ncol(varMat)
      max.nchar <- apply(varMat, 2, function(x){max(nchar(x))})
      cex.1char<- 0.05 #affect distance between text and arrow
      tot.row <- nrow(varMat)
      tot.col2<- tot.col+2 #to enable the fitting lines longer
      text.size1 <- 3/log(tot.row)
      th.a <- th.arrow/log(tot.row) #about arrow locaiton
      
      #for Layer 1
      lty1 <- 1;  lwd1 <- 1; col1='gray80';
      y.tt <- which(varMat[,1]!='')
      pnt <- data.frame(x=c(0, rep(1, length(y.tt))), y=c(1, y.tt), 
                        lab=c('', varMat[y.tt,1]),
                        pch=rep(22,length(y.tt)+1), cex=rep(3,length(y.tt)+1))
      arr <- data.frame(x0=rep(0, length(y.tt)), 
                        y0=rep(1, length(y.tt)), 
                        x1=rep(1, length(y.tt)), 
                        y1=y.tt, 
                        lty=rep(1, length(y.tt)),
                        lwd=rep(1, length(y.tt)), 
                        col=rep(col1, length(y.tt)))  
      shf <- sum(max.nchar[1])*cex.1char
      y.tt0 <- y.tt
      for(i in 2:tot.col){
        y.tt <- which(varMat[,i]!='')
        pnt <- rbind(pnt, data.frame(x=rep(i+shf, length(y.tt)), 
                                     y=y.tt, lab=varMat[y.tt,i],
                                     pch=22, cex=3) )
        wh.a1<-which(!y.tt%in%y.tt0)
        y.tt0a <- y.tt
        for(a in wh.a1){
          y.tt0a[a] <- y.tt0a[a-1]
        }
        arr <- rbind(arr,
                     data.frame(x0=rep(i-0.5+shf, length(y.tt)), 
                                y0=y.tt0a,
                                x1=rep(i+shf, length(y.tt)),
                                y1=y.tt, 
                                lty=rep(lty1, length(y.tt)),
                                lwd=rep(lwd1, length(y.tt)),
                                col=rep(col1, length(y.tt))))
        shf <- sum(max.nchar[1:i])*cex.1char
        y.tt0 <- y.tt
      }
      
      if(showPlot){
        par(mar=c(0.2, 0.2, 0.2, 0.2), mfrow=c(1,1))
        plot(y~x, data=pnt, col='gray80', pch=pnt$pch,
             ylim=c(0, tot.row+0.5), xlim=c(0, tot.col2+shf+0.5),
             axes=F, ylab='', xlab='')
        text(x=pnt$x, y=pnt$y, labels=pnt$lab, adj=-0.07)
        arrows(x0=arr$x0, y0=arr$y0, x1=arr$x1, y1=arr$y1, 
               length=0.1, lty=arr$lty, lwd=arr$lwd, 
               col=arr$col)
        text(x=0, y=tot.row+0.5, labels='Utility', col='blue')
        text(x=tot.col2+shf, y=tot.row+0.5, labels='Risk', col='red')
      }
    }
    
    #add barplot of utility and loss
    if(showBar&showPlot){
      #      bar.wid <- 8
      u.x0<-rep(0, length(LV))
      u.y0<-which(varMat[,ncol(varMat)-1]!='')
      u.x1<-unlist(U)
      if(useTteDiff){utte.x1<-unlist(Utte)}else{utte.x1<-NULL}
      col.bar1 <- rgb(0, 0, 255, alpha=80, maxColorValue=255) #blue
      col.bar1tte <- rgb(0, 100, 0, alpha=80, maxColorValue=255)  #darkgreen
      abline(v=0, col=col.bar1)
      for(i in 1:length(u.y0)){
        if(onebar & useTteDiff){
          uInt.x1<-u.x1[i]*utte.x1[i]
          lines(x=c(u.x0[i], uInt.x1), 
                y=c(u.y0[i], u.y0[i]),
                lwd=bar.wid, col=col.bar1)
          text(x=u.x0[i], y=u.y0[i], labels=round(uInt.x1, 3))
        }else{
          lines(x=c(u.x0[i], u.x1[i]), y=c(u.y0[i], u.y0[i]),
                lwd=bar.wid, 
                col=col.bar1)
          if(useTteDiff){
            lines(x=c(u.x0[i], utte.x1[i]), y=c(u.y0[i], u.y0[i])+0.1,
                  lwd=bar.wid, 
                  col=col.bar1tte)
          }
        }
      }
      
      l.x0<-rep(tot.col2+shf, nrow(varMat))
      l.y0<-1:nrow(varMat)
      l.x1<-l.x0-unlist(E.L)
      col.bar2 <- rgb(255, 0, 0, alpha=80, maxColorValue=255)
      abline(v=tot.col2+shf, col=col.bar2)
      for(i in 1:length(l.y0)){
        lines(x=c(l.x0[i], l.x1[i]), y=c(l.y0[i], l.y0[i]),
              lwd=bar.wid, 
              col=col.bar2)
      }
    }
    
    return(list(dat=varMat, BayesLoss=E.L, U=U, Utte=Utte, p=p.1))
    
  }
  
  
}
#End 3. -----------------------------------------------------------------------#


#Begin 3.1 ---------------------------------------------------------------------#
#1.dynamic CSF with iBDT
if(TRUE){
  iBDT_CSF<-function(
    levInt=paste(c("Cervical::2L_PD-L1_CPS>1%",  "NSCLC::1L_IIIB/IV_PD-L1>5%",
                   "NSCLC::2L_IIIB/IV_PD-L1>5%", "NSCLC::3L_IIIB/IV_PD-L1>5%",
                   "NSCLC::1Lplatinum_IIIB/IV"),collapse=','), 
    #set cohort lables to compare, '::' defines levels
    n1c='40',  #fixed a samples for other parameter change
    orr0="0.146, 0.22, 0.19, 0.18, 0.137", #ORRs of controls match to lev1
    delta.r='0.10',  #fixed an improved ORR effect size
    delta.t='4', #fixed TTE effect size month
    mt0='2.1, 5.4, 2.8, 2.8, 2.7', #tte values of control or benchmarks
    cv='0.8', #coefficient variation for tte
    cutoff='0.4', #decision threshold over ORR or delta ORR
    nsel="20, 40, 60, 80, 100, 120, 140, 160",#simulate different samples
    dorr="0.05, 0.1, 0.15, 0.2", #simulate different ORR effect size
    dtte="1, 2, 3, 4, 5, 6, 7", #simulate different TTE effect size
    dr1r='0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70', 
    #simulate decision threshold
    yellowC="0.25, 0.4", #yellow range for ORR decision thresholds
    showAll=TRUE #if F, only show the iBDT risk vs threshold
  ){
    #set up basic parameters
    if(TRUE){
      #levInt<-paste(lev1, collapse = ',')
      lev1<- unlist(strsplit(levInt, split=',', fixed=T)[[1]])
      ngrp <- length(lev1)
      c2n<-function(x, s=NULL){
        if(is.numeric(x)) return(x)
        if(is.null(s))
          x<-as.numeric(unlist(strsplit(as.character(x), split=',')[[1]]))
        else
          x<-as.numeric(unlist(strsplit(as.character(x), split=',')[[1]]))[s]
        return(x)
      }#covert character to vector
      n1c<-c2n(n1c,1); orr0<-c2n(orr0); delta.r<-c2n(delta.r,1);
      delta.t<-c2n(delta.t,1); mt0<-c2n(mt0); cv<-c2n(cv,1); 
      cutoff<-c2n(cutoff,1); nsel<-c2n(nsel); dorr<-c2n(dorr);
      dtte<-c2n(dtte); dr1r<-c2n(dr1r); yellowC<-c2n(yellowC,1:2)
      
      dr_1Int<-paste0(rep('stop::go', len=ngrp), collapse=',')
      #bmk prevalence internal setup
      ic_1Int<-paste0(rep(1, len=ngrp), collapse=',') 
      n1Int  <- paste(rep(n1c, len=ngrp),collapse=',')
      
      #set up ORR trt vs ctrl labels
      orrV<-paste(orr0+delta.r, orr0, sep='vs')
      orr1Int<-paste(orrV,collapse=', ')
      #orr1Int <- paste(orr0+delta.r, collapse=', ')
      
      #set up TTE trt vs ctrl lables      
      mt1V<-paste(mt0+delta.t, mt0, sep='vs')
      mt1Int <- paste(mt1V,collapse=', ')
      sd0 <- mt0*cv #standard deviation of TTE
      sd1Int <- paste(paste(sd0, sd0, sep='vs'),collapse=',')
      
      #decision rule threshold for ORR
      dr1Int <- paste(rep(cutoff, len=ngrp), collapse=',')
      poInt  <-c(1,0) #convert utility to success rate
    }
    
    p1<-function(nsel.p1=nsel, ngrp.p=ngrp, dlt.r=delta.r, dlt.t=delta.t,
                 lev2=levInt, dr_1=dr_1Int, ic_1=ic_1Int, orr1=orr1Int,
                 mt1=mt1Int, sd1=sd1Int, po=poInt, dr1=dr1Int,
                 lev1p=lev1, orr0p=orr0, mt0p=mt0 ){ #change n
      u1<-u2<-r1<-matrix(NA, nrow=ngrp.p, ncol=length(nsel.p1))
      for(i in 1:length(nsel.p1)){#start looping
        n1  <- paste(rep(nsel.p1[i], len=ngrp.p),collapse=',')
        ibdt1<-BDT_UaL.diff( levVars=lev2, dr_lev=dr_1, incidence=ic_1, 
                             numRsp=orr1, muTTE=mt1,  sdTTE=sd1,  
                             n_ij=n1, dr_th=dr1, showPlot=F, payoff=po)
        u1[,i]<-unlist(ibdt1$U)
        u2[,i]<-unlist(ibdt1$Utte)
        r1[,i]<-sapply(ibdt1$BayesLoss, function(x){x[2]})
      }
      u <- sqrt(u1*u2)
      colnames(u)<-colnames(r1)<-nsel.p1 #start plotting
      plot(u[1,]~nsel.p1, ylim=range(u),  type='o', 
           ylab='iBDT Utility',  xlab='n',
           main=paste0('delta_ORR=', 100*dlt.r, '%, delta_TTE=', dlt.t,'mon'))
      j<-2; while(j>=2 & j<=nrow(u)){
        points(u[j,]~nsel.p1, type='o', col=j)
        j<-j+1
      }
      col1<-1:nrow(u)
      legend('bottomright', 
             legend=paste0(lev1p,'orr_ctrl',orr0p, ' mPFS', mt0p), 
             text.col=col1, bty='n', lty=1, col=col1)
    }
    
    
    p2<-function(dorr.p2=dorr,n.1=n1c, ngrp.p=ngrp, dlt.r=delta.r, dlt.t=delta.t,
                 lev2=levInt, dr_1=dr_1Int, ic_1=ic_1Int,  n1=n1Int,
                 mt1=mt1Int, sd1=sd1Int, po=poInt, dr1=dr1Int, #orr1=orr1Int,
                 lev1p=lev1, orr0p=orr0, mt0p=mt0){ #change delata
      u1<-u2<-r1<-matrix(NA, nrow=ngrp.p, ncol=length(dorr.p2))
      for(i in 1:length(dorr.p2)){#start looping
        orr1<-paste(paste(orr0p+dorr.p2[i], orr0p, sep='vs'), collapse=', ')
        ibdt1<-BDT_UaL.diff( levVars=lev2, dr_lev=dr_1, incidence=ic_1, 
                             numRsp=orr1, muTTE=mt1,  sdTTE=sd1,  n_ij=n1, 
                             dr_th=dr1, showPlot=F, payoff=po)
        u1[,i]<-unlist(ibdt1$U)
        u2[,i]<-unlist(ibdt1$Utte)
        r1[,i]<-sapply(ibdt1$BayesLoss, function(x){x[2]})
      }
      u <- sqrt(u1*u2)
      colnames(u)<-colnames(r1)<-dorr.p2
      plot(u[1,]~dorr.p2, ylim=range(u),  type='o', ylab='iBDT Utility', 
           xlab='ORR effect size', 
           main=paste0('n=',n.1,', delta_TTE=', dlt.t, 'mon'))
      j<-2; while(j>=2 & j<=nrow(u)){
        points(u[j,]~dorr.p2, type='o', col=j)
        j<-j+1
      }
      col1<-1:nrow(u)
      legend('bottomright', 
             legend=paste0(lev1p,'orr_ctrl',orr0p, ' mPFS', mt0p), 
             text.col=col1, bty='n', lty=1, col=col1)
    }
    
    p3<-function(dtte.p3=dtte,n.1=n1c, ngrp.p=ngrp, dlt.r=delta.r, dlt.t=delta.t,
                 lev2=levInt, dr_1=dr_1Int, ic_1=ic_1Int,  n1=n1Int,
                 sd1=sd1Int, po=poInt, dr1=dr1Int, orr1=orr1Int, #mt1=mt1Int, 
                 lev1p=lev1, orr0p=orr0, mt0p=mt0 ){ #change delata TTE
      u1<-u2<-r1<-matrix(NA, nrow=ngrp.p, ncol=length(dtte.p3))
      for(i in 1:length(dtte.p3)){#start looping
        mt1 <- paste(paste(mt0p+dtte.p3[i], mt0p, sep='vs'),collapse=', ')
        ibdt1<-BDT_UaL.diff( levVars=lev2, dr_lev=dr_1, incidence=ic_1, 
                             numRsp=orr1, muTTE=mt1,  sdTTE=sd1,  n_ij=n1, 
                             dr_th=dr1, showPlot=F, payoff=po)
        u1[,i]<-unlist(ibdt1$U)
        u2[,i]<-unlist(ibdt1$Utte)
        r1[,i]<-sapply(ibdt1$BayesLoss, function(x){x[2]})
      }
      u <- sqrt(u1*u2)
      colnames(u)<-colnames(r1)<-dtte.p3
      plot(u[1,]~dtte.p3, ylim=range(u, na.rm=T),  type='o', 
           ylab='iBDT Utility', xlab='TTE effect size', 
           main=paste0('n=',n.1, ', delta_ORR=,', dlt.r*100, '%'))
      j<-2; while(j>=2 & j<=nrow(u)){
        points(u[j,]~dtte.p3, type='o', col=j)
        j<-j+1
      }
    }
    
    p4<-function(dr1r.p4=dr1r, n.1=n1c, ngrp.p=ngrp, dlt.r=delta.r, dlt.t=delta.t,
                 lev2=levInt, dr_1=dr_1Int, ic_1=ic_1Int,  n1=n1Int,
                 sd1=sd1Int, po=poInt, dr1=dr1Int, orr1=orr1Int, mt1=mt1Int, 
                 lev1p=lev1, orr0p=orr0, mt0p=mt0,
                 dc_rg=yellowC){#change CSF with different rate
      u1<-u2<-r1<-matrix(NA, nrow=ngrp.p, ncol=length(dr1r.p4))
      for(i in 1:length(dr1r.p4)){#start looping
        dr1 <- paste(rep(dr1r.p4[i], ngrp.p), collapse=',')
#print('dr1');print(dr1); print(dr_1); print(ic_1); print(orr1);
        ibdt1<-BDT_UaL.diff( levVars=lev2, dr_lev=dr_1, incidence=ic_1, 
                             numRsp=orr1, muTTE=mt1,  sdTTE=sd1,  n_ij=n1, 
                             dr_th=dr1, showPlot=F, payoff=po)
        u1[,i]<-unlist(ibdt1$U)
        u2[,i]<-unlist(ibdt1$Utte)
        r1[,i]<-sapply(ibdt1$BayesLoss, function(x){x[2]})
#print('BayL'); print(ibdt1$BayesLoss)
      }
      u <- sqrt(u1*u2)
      colnames(u)<-colnames(r1)<-dr1r.p4
      plot(0~min(dr1r.p4), ylim=range(r1),  type='o', 
           ylab='iBDT risk = expected distance to success', 
           col='white', xlim=range(dr1r.p4), xlab='ORR decision threshold', 
           main=paste0('n=',n.1,', delta_ORR=', dlt.r*100, '%'))
      polygon(x=c(min(dr1r.p4), dc_rg[1], dc_rg[1], min(dr1r.p4)), 
              y=rep(range(r1), each=2), col='red', border='red')
      polygon(x=c(dc_rg[1], dc_rg[2], dc_rg[2], dc_rg[1]), 
              y=rep(range(r1), each=2), col='yellow', border='yellow')
      polygon(x=c(dc_rg[2], max(dr1r.p4), max(dr1r.p4), dc_rg[2]),
              y=rep(range(r1), each=2), col='green', border='green')
      j<-1; while(j>=1 & j<=nrow(u)){
        points(r1[j,]~dr1r.p4, type='o', col=j)
        j<-j+1
      }
#print('r1:'); print(r1)
      col1<-1:nrow(u)
      # legend('topright', legend=paste0(lev1,', orr_ctrl',orr0, ', mPFS', mt0), 
      # 			 text.col=col1, bty='n', lty=1, col=col1)
    }
    
    if(showAll){
      par(mfrow=c(1,2))
      #p1();  p4(); 
      p2(); p3();
    }else{
      par(mfrow=c(1,2))
      p1(); p4();

    }
  }
  if(F){#BEACH code
    input<-NULL
    
    input$text<-"Cervical::2L_PD-L1_CPS>1%,NSCLC::1L_IIIB/IV_PD-L1>5%,NSCLC::2L_IIIB/IV_PD-L1>5%,NSCLC::3L_IIIB/IV_PD-L1>5%,NSCLC::1Lplatinum_IIIB/IV"
    
    input$text2<-'40'
    input$text3<-"0.146, 0.22, 0.19, 0.18, 0.137"
    input$text4<-'0.10'
    input$text5<-'4'
    input$text6<-'2.1, 5.4, 2.8, 2.8, 2.7'
    input$text7<-'0.8'
    input$text8<-'0.4'
    input$text9<-"20, 40, 60, 80, 100, 120, 140, 160"
    input$text10<-"0.05, 0.1, 0.15, 0.2"
    input$text11<-"1, 2, 3, 4, 5, 6, 7"
    input$text12<-'0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70'
    input$text13<-"0.25, 0.4"
    
    iBDT_CSF(
      levInt=input$text, #set cohort lables to compare, '::' defines levels
      n1c=input$text2,  #fixed a samples for other parameter change
      orr0=input$text3, #ORRs of controls match to lev1
      delta.r=input$text4,  #fixed an improved ORR effect size
      delta.t=input$text5, #fixed TTE effect size month
      mt0=input$text6, #tte values of control or benchmarks
      cv=input$text7, #coefficient variation for tte
      cutoff=input$text8, #decision threshold over ORR or delta ORR
      nsel=input$text9,#simulate different samples
      dorr=input$text10, #simulate different ORR effect size
      dtte=input$text11, #simulate different TTE effect size
      dr1r=input$text12, #simulate decision threshold
      yellowC=input$text13, #yellow range for ORR decision thresholds\
      showAll=F
    )
  }
}
#End 3.1 -----------------------------------------------------------------------#


#Begin 3.2 ---------------------------------------------------------------------#
if(T){ #######2020-10-29
  library(extraDistr)

  #convert a text to num
  c2n<<-function(x, rm.c=c('c','(',')'), sp=',', sp2=';'){
    x<-as.character(x)
    for(jj in rm.c)
      x <- gsub(jj, '', x, fixed=T)
    if(grepl(sp2, x)){
      x3 <- unlist(strsplit(x, split=sp2, fixed=T))
      x2 <- lapply(x3, function(m){
              as.numeric(unlist(strsplit(m, split=sp, fixed=T)))})
    }else{
      x2 <- as.numeric(unlist(strsplit(x, split=sp, fixed=T)))
    }
    return(x2)
  }
  
  #generate a futility table for no go decision
  futTb <<- function(
       n=c(6, 8, 9), #sample sizes
       nr=c(0, 1, 2, 4), #number of observed responders
       trt.orr=c(0.1, 0.2), #target ORRs
       dlt.orr=c(0.05, 0.05), #ORR effect size
       show.type=1, #1=both posterior prob and iBDT risk; 2=iBDT only; 3=post prob
       beta.parm=c(0.2, 1.8), #hyper prior only used in Bayesian Posterior only
       num.dig=3, #number of digits to keep
       direc='1. P(X>r|data)' #for success probability, or'2. P(X<=r|data)' for futility
      ){
     if(is.character(show.type)){show.type<-as.numeric(substring(show.type,1,1))}
     #setup number of responders, sample sizes
     trt.orr <- sort(unique(trt.orr[trt.orr>=0]))
     n <- sort(unique(n[n>0]))
     trt.orrs <- rep(trt.orr, each=length(n))
     trt.id <- rep(1:length(trt.orr), each=length(n))
     dlt.orrs <- rep(dlt.orr, each=length(n))
     ns <- rep(n, times=length(trt.orr))
     if(is.list(beta.parm)){beta.L <- beta.parm}else{beta.L<-NULL}

     #set up rsp labels and fut table
     nr <- sort(unique(nr))
     lab.rsp<-paste0('#rsp=', nr)
     fut<-matrix(NA, nrow=length(n)*length(trt.orr), ncol=length(lab.rsp)) 
     colnames(fut)<-lab.rsp

     show.type<-as.character(show.type)
     for(i in 1:nrow(fut)){
       if(!is.null(beta.L)){beta.parm<-beta.L[[trt.id[i]]]}
       irsk<-NULL #get iBDT risk
       for(j in 1:length(nr)){
         #expected similarity between observed data and target data
         rsk<-BDT_UaL.diff(levVars='my.tst', dr_lev="stop::go", incidence="1", 
              numRsp=paste0(trt.orrs[i],'vs', trt.orrs[i]-dlt.orrs[i]), 
              muTTE=3, sdTTE=0.8, # not useful for rsk result
              n_ij=ns[i], dr_th=nr[j]/ns[i], showPlot=F, payoff=c(1,0))
         irsk <- c(irsk, 1-rsk$BayesLoss[[1]][2])
       }#Similarity to target success?wwith number of responders (m) improved by treatment
       #posterior parameter expectation
       p2 <- 1- pbeta(trt.orrs[i], beta.parm[1]+nr, beta.parm[2]+ns[i]-nr)#Posterior Prob to success
       p4 <- 1- pbbinom(nr, ns[i], 1+trt.orrs[i], 2-trt.orrs[i]) #Predictive Prob to success
       p5 <- pbinom(nr, ns[i], trt.orrs[i])        #false negative P(Y<=nr|targetORR)
       p6 <- 1-pbinom(nr, ns[i], trt.orrs[i]-dlt.orrs[i]) #false positive P(Y>nr|controlORR)
       if(substring(direc, 1,1)!='1'){
          irsk <- 1-irsk
          p2 <- 1-p2 #Posterior Prob to failure
          p4 <- 1-p4 #Predictive Prob to failure
          p5 <- 1-p5 #true negative
          p6 <- 1-p6 #true positive
       }
       p2 <- round(p2, num.dig)
       p4 <- round(p4, num.dig)
       p5 <- round(p5, num.dig)
       p6 <- round(p6, num.dig)
       irsk<-round(irsk, num.dig)
       if(show.type==4){ #FN
         fut[i,] <- p5
       }else if (show.type==5){ #FP
         fut[i,] <- p6
       }else if(show.type==100) {#iBDT posterior for number of responders
         fut[i,] <- paste(p2, 'vs', p4)
         lst.row<-'Posterior Prob Beta vs BetaBinomial'
       }else if(show.type==101) {#iBDT posterior for number of responders
         fut[i,] <- p4
         lst.row<-'Posterior Prob BetaBinomial'
       }else if(show.type==2){ #both table
         fut[i,] <- paste(p2,'vs', irsk)
         lst.row<-'P(ORR>target_orr)  vs  iBDT risk for the expected distance to success. Effect size=0.05'
       }else if(show.type==1){ #iBDT risk
         fut[i,] <- irsk
         lst.row<-'iBDT risk for the expected distance to success. Effect size=0.05'
       }else { #posterior for response rate  #show.type==3
         fut[i,] <- p2
         lst.row<-'P(ORR>target_orr)'
       }
     }
     ot <- data.frame(targetORR=trt.orrs, sampleSize=ns, fut)
     ot[ot==0]<-paste0('<0.', paste(rep('0', num.dig-1), collapse=''), '1')
     colnames(ot)[-(1:2)] <- lab.rsp
     if(show.type%in%c(1,5)){ot <- cbind(controlORR=trt.orrs-dlt.orrs, ot)}
     #ot <- rbind(ot, c('', '', lst.row, rep('', ncol(ot)-3)))
     return(ot)
  }

  #test result
  if(F){

    source('functions/users/DT_frame2.r')

    input<-NULL
    input$text <- '6, 9, 12'
    input$text2 <- '0, 1, 2, 4'
    input$text3 <- '0.1, 0.15'
    input$text4 <- '0.05, 0.05'
    input$text5 <- '0.2, 1.8; 0.3, 1.7;0.4, 16'
#    input$radio <- '0. FN'
    input$radio <- '4. FP'
    input$dropdown<-'3'
    input$radio2<- '1.'

    futTb(n=c2n(input$text), #sample sizes
       nr=c2n(input$text2), #number of observed responders
       trt.orr=c2n(input$text3), #target ORRs
       dlt.orr=c2n(input$text4), #ORR effect size
       show.type=input$radio, #1=both posterior prob and iBDT risk; 2=iBDT only; 3=post prob
       beta.parm=c2n(input$text5), # Hyper parameter in Beta distribution
       num.dig=c2n(input$dropdown),
       direc=input$radio2 )

  }
}
#End 3.2 -----------------------------------------------------------------------#




#Begin 4. ---------------------------------------------------------------------#
#power analysis
if(TRUE){
  #X ~ w*N(muA, sgmS) + (1-w)*N(muB, sgmS)
  #mu = w*muA + (1-w)*muB where muA>=muB
  #sgmS = sqrt(sgm^2 - ee)
  #ee = w*muA^2 + (1-w)*muB^2 - mu^2
  #theta = muA - muB
  BEACH_PowerPlot<<-function(
    n.rg='10~200',
    n.vl=20, #a fixed vertical line in power vs n plot
    theta.rg='0~200',
    theta.vl=50, #a fixed vertical line in power vs theat
    sgm.rg='10~200', 
    sgm.vl=100, #a fixed vertical line in power vs observed standardization
    muA=15, #a default assumption of the lower group mean
    prev='0.15~0.65', #prevalence rate for group A
    alpha=0.05,  #significance level
    simSize=100, #simulation size
    m1=1,           #number of group pairs
    r1=ifelse(as.numeric(m1)<2,  NA,
              paste(rep(1,round(as.numeric(m1))-1),collapse=',') )
    #proportions in the groups with higher means
  ){
    require(pwr)
    xx2<-function(rg1, ss=simSize){
      return(seq(rg1[1], rg1[2], len=ss))
    }
    add1<-function(xL='sample size', LL, v1=n.vl, v1c='n', 
                   LG1=NULL, LG1.pos='topright', 
                   LG2=NULL, LG2.pos='bottomright'){
      plot(y~x, data=LL[[1]], ylab='power', xlab=xL, type='l', lty=2,
           ylim=c(0,1), xlim=range(LL[[1]]$x))
      if(length(LL)>1){
        for(i in 2:length(LL)) lines(y~x, data=LL[[i]], lty=i+1)
      }
      abline(v=v1, col='red')
      mtext(text=paste0(v1c,'=', v1), at=v1, col='red')
      if(!is.null(LG1)){
        legend(LG1.pos,legend=LG1, text.col='red', bty='n')
      }
      if(!is.null(LG2)){
        legend(LG2.pos,legend=LG2, bty='n')
      }
    }
    #convert number of group paris into integer
    m1<-round(as.numeric(m1))
    #get the sample size ratios of other group to group a
    r1<-as.numeric(strsplit(r1, split=',', fixed=T)[[1]])
    if(!is.na(r1)){
      l_r1 <- m1 - length(r1)
      if(l_r1>0){
        r1<-c(r1, rep(1, l_r1)) #set default ratio as 1 if missing
      }
      r1<-r1[1:m1]
      r1<-r1/sum(r1)
      alpha<-round(1-(1-alpha)^(1/m1),3)#set type I error for each test
    }
    
    n.rg<-spF(n.rg); #set total sample size range
    theta.rg<-spF(theta.rg); #set effect size range
    sgm.rg<-spF(sgm.rg); #set standard deviation range
    prev<-spF(prev, s='~|,', noClean=TRUE); 
    prev<-prev[prev>0&prev<1];#set the prevalences of the lowest group mean
    if(alpha>1|alpha<0){alpha<-0.05}; #set default significance level
    
    cv.vl<-sgm.vl/theta.vl; #set coefficient variation selected value
    cv.rg<-c(min(cv.vl,0.01), max(cv.vl+1, 2)); #set range for cv
    d1t<-abs(theta.vl)/sgm.vl; #set the selected effect size
    
    #get the point value
    y1<-nA1<-nB1<-NULL; #note B1 could have >1 values if m>1
    for(i in 1:length(prev)){
      nA1t<-round(n.vl*prev[i]); nB1t<-n.vl-nA1t
      if(is.na(r1)){
        y1t<-pwr.t2n.test(n1=nA1t, n2=nB1t, d=d1t, sig.level=alpha)$power
      }else{
        nB1t1<- round(nB1t*r1)
        nB1t1<-nB1t1[nB1t1>1]
        y1t <- pwr.t2n.test(n1=nA1t, n2=nB1t1, d=d1t, sig.level=alpha)$power
        y1t <- 1-prod(1-y1t)
      }
      y1<-c(y1, y1t)
      nA1<-c(nA1, nA1t)
      nB1<-c(nB1, nB1t) 
    }
    leg1<-paste0(paste0('rA=',prev), paste0(', pwr=',round(y1,2)))
    leg1<-leg1[order(y1, decreasing=T)]
    
    par(mfrow=c(2,2))
    #power vs n
    x<-round(xx2(n.rg, ss=simSize)); pnL<-list(); 
    for(i in 1:length(prev)){
      nAs<-pmax(2, round(x*prev[i])); nBs<-pmax(2, x-nAs); 
      if(is.na(r1)){
        y<-pwr.t2n.test(n1=nAs, n2=nBs, d=d1t, sig.level=alpha)$power
      }else{
        yt<-rep(1, length(x))
        for(r_m in 1:length(r1)){
          nBs1<-round(r1[r_m]*nBs)
          s0<-nBs1>1
          if(all(!s0)){next}
          yt[s0]<-yt[s0]*(1-pwr.t2n.test(n1=nAs[s0], n2=nBs1[s0], d=d1t, 
                                         sig.level=alpha)$power)
        }
        y<-1-yt
      }
      pnL[[i]]<-data.frame(x=x, y=y)
    }
    leg2<-paste0('meanA=',muA, ', diff=',theta.vl,
                 ', sd=', sgm.vl, ', sig.level=', alpha)
    add1(xL='sample size', LL=pnL, v1=n.vl, v1c='n', LG1=leg1,LG2=leg2)
    
    ##power vs cv
    x<-xx2(cv.rg, ss=simSize); pcL<-list(); 
    for(i in 1:length(prev)){
      if(is.na(r1)){
        y<-pwr.t2n.test(n1=nA1[i], n2=nB1[i], d=1/x, sig.level=alpha)$power
      }else{
        yt<-rep(1, length(x))
        for(r_m in 1:length(r1)){
          nBs1<-round(r1[r_m]*nB1[i])
          if(nBs1<2) next
          yt<-yt*(1-pwr.t2n.test(n1=nA1[i], n2=nBs1, d=1/x, sig.level=alpha)$power)
        }
        y<-1-yt
      }
      pcL[[i]]<-data.frame(x=x, y=y)
    }
    leg2<-paste0('meanA=',muA, ', n=',n.vl, ', sig.level=', alpha)
    add1(xL='coefficient of variation (sd/mean)', LL=pcL, v1=round(cv.vl,2), 
         v1c='cv', LG2=leg2, LG2.pos='topright')
    
    #power vs sd
    x<-xx2(sgm.rg, ss=simSize); psL<-list(); 
    muB<-muA+theta.vl
    for(i in 1:length(prev)){
      mmt<-prev[i]*muA^2+(1-prev[i])*muB^2-(prev[i]*muA+(1-prev[i])*muB)^2
      xS<-sqrt(x^2-mmt)
      if(is.na(r1)){
        y<-pwr.t2n.test(n1=nA1[i], n2=nB1[i], 
                        d=abs(theta.vl)/xS, sig.level=alpha)$power
      }else{
        yt<-rep(1, length(x))
        for(r_m in 1:length(r1)){
          nBs1<-round(r1[r_m]*nB1[i])
          if(nBs1<2) next
          yt<-yt*(1-pwr.t2n.test(n1=nA1[i], n2=nBs1, 
                                 d=abs(theta.vl)/xS, sig.level=alpha)$power)
        }
        y<-1-yt
      }
      psL[[i]]<-data.frame(x=x, y=y)
    }
    leg2<-paste0('meanA=',muA, ', diff=',theta.vl,
                 ', n=', n.vl, ', sig.level=', alpha)
    add1(xL='total standard deviation', LL=psL, v1=sgm.vl, v1c='sd', LG2=leg2, 
         LG2.pos='topright')
    
    #power vs diff
    x<-xx2(theta.rg, ss=simSize); pdL<-list(); 
    for(i in 1:length(prev)){
      if(is.na(r1)){
        y<-pwr.t2n.test(n1=nA1[i], n2=nB1[i], 
                        d=abs(x)/sgm.vl, sig.level=alpha)$power
      }else{
        yt<-rep(1, length(x))
        for(r_m in 1:length(r1)){
          nBs1<-round(r1[r_m]*nB1[i])
          if(nBs1<2) next
          yt<-yt*(1-pwr.t2n.test(n1=nA1[i], n2=nBs1, 
                                 d=abs(x)/sgm.vl, sig.level=alpha)$power)
        }
        y<-1-yt
      }
      pdL[[i]]<-data.frame(x=x, y=y)
    }
    leg2<-paste0('meanA=',muA, ', sd=',sgm.vl,
                 ', n=', n.vl, ', sig.level=', alpha)
    add1(xL='absolute difference', LL=pdL, v1=theta.vl, 
         v1c='meanB - meanA', LG2=leg2)
    
    #density plot
    if(F){
      dsL<-list(); mu.rg<-c(muA-3*sgm.vl, muA+theta.vl+3*sgm.vl);
      x<-xx2(mu.rg); yrg<-NULL;
      muB<-muA+theta.vl; 
      sgmS<-NULL
      for(i in 1:length(prev)){
        mmt<-prev[i]*muA^2+(1-prev[i])*muB^2-(prev[i]*muA+(1-prev[i])*muB)^2
        sgmSt<-sqrt(sgm.vl^2-mmt)
        ya<-prev[i]*dnorm(x, mean=muA, sd=sgmSt)
        yb<-(1-prev[i])*dnorm(x, mean=muB, sd=sgmSt)
        sgmS<-c(sgmS, sgmSt)
        dsL[[i]]<-data.frame(x=x, ya=ya, yb=yb)
        yrg<-range(c(yrg, range(ya), range(yb)))
      }
      yrg<-c(0,max(yrg))
      plot(0~0, ylim=yrg, xlim=range(x), col='white', ylab='density', 
           xlab='means')
      for(i in 1:length(dsL)){
        lines(ya~x, data=dsL[[i]], lty=i+1, col='blue')
        lines(yb~x, data=dsL[[i]], lty=i+1, col='magenta')
      }
      abline(v=c(muA,muB), col=c('blue', 'magenta'))
      mtext(at=(muA+muB)/2, text=paste0('muB-muA=',theta.vl))
      #legend("topleft",legend=leg1, lty=(1:length(prev))+1, bty='n')
    }
  }
  #setup parameters
  spF<<-function(x, s="~", ss=simSize, noClean=FALSE){
    if(is.null(x)){
      rg1<-NULL
    }else if(all(grepl(s, x))){
      x<-as.numeric(strsplit(x, split=s)[[1]])
      rg1<-x[!is.na(x)]
    }else{rg1<-as.numeric(x)}
    if(noClean){return(rg1)}
    if(is.null(rg1)||length(rg1)==0){rg1<-c(0, ss)}#defult setup
    if(length(rg1)==1){rg1<-c(rg1, rg1+ss)}else{rg1<-rg1[1:2]}
    return(rg1)
  }
  
}
#End 4. -----------------------------------------------------------------------#

