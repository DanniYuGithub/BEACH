##############################################################################################
#Source Code received from Profession Jason Hsu
#BEACH CD created by Danni Yu
#Initial Date: 2020-11-19
#Last Modify Date: 2021-1-4
##############################################################################################


if(F){
  getwd(); dir();
  source("functions/users/JasonHsu.R")

  lapply(indataset, head)
  varlist <<- colnames(indataset[[1]])

  #"Simultaneous Confidence Intervals for Efficacy"
  input<-NULL
  input$dropdown<-c("Ratio of Median", "Difference of Median") #"Efficacy Measurement"
  input$slide<-5 # "Alpha (in percents)",  min = 1,  max = 50, value = 5)
  input$slide4<-50 #"Proportion of marker negative subjects in population (in percents)", min = 0,  max = 100, value = 50)
  input$dropdown2<-'Subj' #"Patient ID Variable", varlist
  input$dropdown3<-'Trt'  #"Treatment Variable", varlist
  input$dropdown4<-'0' #"Control Value", trt_levels= sort(unique(indataset[[1]][,input$dropdown3]))  
  input$dropdown5<-'M' #"Biomarker Variable", varlist
  input$slide2<-0.95 #"Data filter cut-point c0 (subjects w/ value < C0 are excluded)", 
                    #min=min(indataset[[1]][, input$dropdown5],na.rm=TRUE)
                    #max=as.numeric(quantile(indataset[[1]][, input$dropdown5], 0.50,na.rm=TRUE))
                    #value=min(indatast[[1]][, input$dropdown5],na.rm=TRUE)
  input$slide3<-1.66 #"Marker +/- cut-point c1 (subjects w/ value < C1 are marker -)"
                    #min=min.inSliderC1 
                    #max=max.inSliderC1
                    #value=median(c(CMarker, RxMarker)
                    #get3(mrk=input$dropdown5, slide0=input$slide2,trtVar=input$dropdown3, conVal=input$dropdown4)
  input$dropdown6<-'Y' #"Time to Event Outcome Variable", varlist
  input$dropdown7<-'Delta' #"Censor Indicator Variable", varlist
  input$dropdown8<-'1' #"Censor Value", censor_levels = c("", sort(unique(indataset[[1]][, input$dropdown7])))

  intiPlot(meth=input$dropdown)

  #grepl('Relative Response',input$dropdown, fixed=T)

}

#-----------------------------------------#
#Author: Danni Yu
#Additional functions created to adapt into BEACH. 
if(TRUE){
  library(eha)
  library(survival)
  library(rootSolve)
  library(mvtnorm)
  library(plotrix)

  varlist <<- colnames(indataset[[1]])

  fill.sp<<-function(xx, nch, sep=' '){
    yy<-NULL
    for(i in 1:length(xx)){
      yy<-c(yy, paste0(xx[i], paste(rep(' ', max(0, nch-nchar(xx[i]))), collapse='') ))
    }
    return(yy)
  }

  get3<<-function(mrk=input$dropdown5, slide0=input$slide2, 
                     trtVar=input$dropdown3, conVal=input$dropdown4, 
                     dd=indataset[[1]]){
    trt<-dd[dd[, mrk]>=slide0, trtVar]
    Marker<-dd[dd[, mrk]>=slide0, mrk]
    CMarker<-Marker[trt==conVal]
    RxMarker<-Marker[trt!=conVal]
    ind.C.max<-order(-CMarker,na.last=NA)[2] # na.last=NA for missing values
    ind.Rx.max<-order(-RxMarker,na.last=NA)[2]
    ind.C.min<-order(CMarker,na.last=NA)[2]
    ind.Rx.min<-order(RxMarker,na.last=NA)[2] 

    if(CMarker[order(CMarker,na.last=NA)[2]]==CMarker[order(CMarker,na.last=NA)[3]]){
       min.inSliderC1<-max(CMarker[ind.C.min], RxMarker[ind.Rx.min],na.rm=TRUE)+0.01
       max.inSliderC1<-min(CMarker[ind.C.max], RxMarker[ind.Rx.max],na.rm=TRUE)
     }else{
       min.inSliderC1<-max(CMarker[ind.C.min], RxMarker[ind.Rx.min],na.rm=TRUE)
       max.inSliderC1<-min(CMarker[ind.C.max], RxMarker[ind.Rx.max],na.rm=TRUE)
     }  
     return(list(min=min.inSliderC1, max=max.inSliderC1, 
                 med=median(c(CMarker, RxMarker)), 
                 by=(max.inSliderC1-min.inSliderC1)/20))
  }

  ## ---- model fit reactive object ----##
  ## ---- analysis data reactive object ----
  analysisData <<- function(dd2=input$dropdown2, 
       dd3=input$dropdown3, dd4=input$dropdown4, dd5=input$dropdown5,
       dd6=input$dropdown6, dd7=input$dropdown7, dd8=input$dropdown8, 
       ss2=input$slide2, ss3=input$slide3
    ){
    if (is.null(dd2)|  
          is.null(dd3)|  
          is.null(dd4)|
          is.null(dd5)|
          is.null(dd6)|
          is.null(dd7)|
          is.null(dd8)) { return() }
    tmp <- indataset[[1]]
    # assemble all variable names selected from various inputs
    cols_select <- c(dd2, dd3, dd5, dd6, dd7)
    
    analysis_set <- tmp[, cols_select]
    if(class(analysis_set[, dd3]) != "character") {
      analysis_set[, dd3] <- as.character(analysis_set[, dd3])
    }
    
    ## --------------------------------------------##
    ## -- remove data with Marker value below C0 --##
    ## --------------------------------------------##
    
    if(sum(indataset[[1]][, dd5]<ss2,na.rm=TRUE)>0){
      analysis_set<- analysis_set[indataset[[1]][, dd5]>=ss2,]
    }
    
    analysis_set[, dd5]<-factor(ifelse(analysis_set[, dd5] < ss3, "Neg", "Pos"))
    analysis_set
  }# end of analysisData 

  WeibullFit <<- function( d3=input$dropdown3, d2=input$dropdown2, d4=input$dropdown4,
       d5=input$dropdown5, d6=input$dropdown6, d7=input$dropdown7, d8=input$dropdown8, 
       s2=input$slide2, s3=input$slide3, s4=input$slide4, s0=input$slide, rtLab=F, metd=input$dropdown
   ){
        anaTemp <<- analysisData(dd2=d2, dd3=d3, dd4=d4, dd5=d5,
                       dd6=d6, dd7=d7, dd8=d8, ss2=s2, ss3=s3)
        if (length(unique(indataset[[1]][, d3]))>10 |
              is.null(d2)|  
              is.null(d3)|  
              is.null(d4)|
              is.null(d5)|
              length(levels(factor(anaTemp[, d5])))<2|
              is.null(d6)|
              is.null(d7)|
              is.null(d8)) { return() }
                
        outcomeEventInd <- anaTemp[, d7]
        outcomeEventInd[outcomeEventInd == d8] <- 0
        outcomeEventInd[outcomeEventInd != d8] <- 1    
        
        # define survival object
        response_vec <- Surv(time = as.numeric(anaTemp[, d6]), 
                             event = as.numeric(outcomeEventInd))
        Trt <- anaTemp[, d3]
        M <- anaTemp[, d5]  
        d4 <- as.character(d4)      
    
        # fit Weibull model that reacts to event in user interface 
        fit <- weib(formula = response_vec ~ 
                    relevel(as.factor(Trt), d4)+
                    relevel(as.factor(M),"Neg")+
                    relevel(as.factor(Trt),d4)*relevel(as.factor(M),"Neg"), 
                  data=anaTemp) 
        # log.param are the log scaled parameter estimates from weibull fitting
        log.param.est <- fit$log.coef
        log.param.var <- diag(fit$log.var)
        # param are the original scale (exponential of log.param) parameter estimates
        param.est <- fit$coef
        param.var <- fit$var
        param.var.diag <- diag(param.var)
        
        # get prevalence of marker negative group from slider 
        p.neg = as.double(s4)/100
        
        #-------------------------------------------------------------------------------#
        # compute median survival times and their associated variance covariance matrix #
        #-------------------------------------------------------------------------------#
        comp.median <- median.surv(lambda=param.est[1],k=param.est[2],theta=param.est[3:5],p.neg=p.neg)
        (median.surv.est <- unlist(comp.median))
    
        median.surv.C <- median.surv.est[c(1,2,5)]
        median.surv.Rx <- median.surv.est[c(3,4,6)]
    
        median.surv.estiMates <- matrix(ncol=3, rbind(median.surv.Rx, median.surv.C))
        Negative <- median.surv.estiMates[,1]
        Positive <- median.surv.estiMates[,2]
        Mixture <- median.surv.estiMates[,3]
        Estimates <- data.frame(Negative, Positive, Mixture)
        row.names(Estimates) <- c("Rx", "C")
    
        median.var.est <-var.median.surv.dt.dim6(med.C.neg=median.surv.est[1], med.C.pos=median.surv.est[2], 
                                             med.Rx.neg=median.surv.est[3], med.Rx.pos=median.surv.est[4],
                                             med.C=median.surv.est[5],med.Rx=median.surv.est[6],p.neg=p.neg,
                                             lambda=param.est[1],k=param.est[2],theta=param.est[3:5],
                                             param.var=fit$var)$median.var        
        #----------------------------------------#
        # Calculate ratio of median and their CI #
        #----------------------------------------#  
        #----------------------------------------------------------#
        # begin log scale delta method
        comp.ratio <- log.ratio.median.surv.dt.dim3(median.surv=median.surv.est,median.var=median.var.est)
        log.ratio <- c(comp.ratio$log.ratio.neg,comp.ratio$log.ratio.pos,comp.ratio$log.ratio.mix)
        var.log.ratio <- comp.ratio$var.log.ratio
        
        original.ratio <- c(comp.ratio$ratio.neg,comp.ratio$ratio.pos,comp.ratio$ratio.mix)
        ratio.negative <- original.ratio[1]
        ratio.positive <- original.ratio[2]
        ratio.mixture  <- original.ratio[3]
        Ratios <- data.frame(ratio.negative, ratio.positive, ratio.mixture)
        row.names(Ratios) <- c("Rx/C")
        colnames(Ratios) <- c("Negative", "Positive", "Mixture")
        
        comp.ratio.CI <- log.ratio.simu.CI(seed=12345, alpha=as.double(s0)/100,
                                       log.ratio=log.ratio, var.log.ratio=var.log.ratio)
        ratio.CI.lbd <- comp.ratio.CI$ratio.CI.lbd
        ratio.CI.ubd <- comp.ratio.CI$ratio.CI.ubd
        
        # end of log delta method
        #----------------------------------------------------------#
        # CIs from lg scale delta method
        ratio.CI <- matrix(ncol=3, rbind(ratio.CI.ubd, ratio.CI.lbd))
        
        ratio.CI.negative <- ratio.CI[,1]
        ratio.CI.positive <- ratio.CI[,2]
        ratio.CI.mixture <- ratio.CI[,3]
        Ratio.CI <- data.frame(ratio.CI.negative, ratio.CI.positive, ratio.CI.mixture)
        row.names(Ratio.CI) <- c("Upper","Lower")
        colnames(Ratio.CI) <- c("Negative", "Positive", "Mixture")
        
    
    
        #---------------------------------------------#
        # Calculate difference of median and their CI #
        #---------------------------------------------#
        #----------------------------------------------------------#
        # begin of delta method
        comp.dif <- dif.median.surv.dt.dim3(median.surv=median.surv.est,median.var=median.var.est)
        dif.median <- c(comp.dif$dif.neg, comp.dif$dif.pos,comp.dif$dif.mix)
        var.dif.median <- comp.dif$var.dif
        
        dif.median.negative <- dif.median[1]
        dif.median.positive <- dif.median[2]
        dif.median.mixture <- dif.median[3]
        Diff <- data.frame(dif.median.negative, dif.median.positive, dif.median.mixture)
        row.names(Diff) <- c("Rx-C")
        colnames(Diff) <- c("Negative", "Positive", "Mixture")
        
        comp.dif.CI <- dif.median.simu.CI(seed=12345,alpha=as.double(s0)/100, 
                                     dif.median=dif.median, var.dif.median=var.dif.median)
        median.dif.CI.lbd <- comp.dif.CI$dif.median.CI.lbd
        median.dif.CI.ubd <- comp.dif.CI$dif.median.CI.ubd
    
        # end of delta method
        #----------------------------------------------------------#
        # CIs from delta method
        median.dif.CI <- matrix(ncol=3, rbind(median.dif.CI.ubd, median.dif.CI.lbd))
        
        median.dif.CI.negative <- median.dif.CI[,1]
        median.dif.CI.positive <- median.dif.CI[,2]
        median.dif.CI.mixture <- median.dif.CI[,3]
        Diff.CI <- data.frame(median.dif.CI.negative, median.dif.CI.positive, median.dif.CI.mixture)
        row.names(Diff.CI) <- c("Upper","Lower")
        colnames(Diff.CI) <- c("Negative", "Positive", "Mixture")
        
        #---------------#
        # return retuls #
        #---------------#    
        WeibullResult <<- list(Estimates = Estimates,
                              Ratios = Ratios,
                              Ratio.CI = Ratio.CI,
                              Diff = Diff,
                              Diff.CI = Diff.CI,
                              median.surv.estiMates = median.surv.estiMates)
        if(!rtLab){
          return(WeibullResult)
        }else{
          printDF<-function(x=WeibullResult$Estimates){
            xx<-cbind(` `=rownames(x), round(x, 3))
            xx<-rbind(c(paste(rep(' ', max(nchar(xx[,1]))), collapse=''), colnames(x)), xx)
            p1<-paste(apply(xx, 1, paste, collapse='\t'), collapse='\n')
          }
          tm1<-paste0("Estimated Median Time-to-Event:\n", printDF(WeibullResult$Estimates))

          tm2<-paste0("Estimated Ratio of Median Time-to-Event:\n", printDF(WeibullResult$Ratios))
          tm3<-paste0("Confidence Interval for Ratio of Median Time-to-Event:\n", printDF(WeibullResult$Ratio.CI))

          tm4<-paste0("Estimated Difference of Median Time-to-Event:\n", printDF(WeibullResult$Diff))
          tm5<-paste0("Confidence Interval for Ratio of Difference of Median Time-to-Event:\n", printDF(WeibullResult$Diff.CI))

          if(metd=="Ratio of Median"){
            rt.label<-paste(tm1, tm2, tm3, sep='\n')
          }else{
            rt.label<-paste(tm1, tm4, tm5, sep='\n')
          }
          return(rt.label)
        }
  } #end of WeibullFit


  #--------------------For plot output in BEACH----------------------#
  intiPlot<<-function( meth=input$dropdown, 
      drd3=input$dropdown3, drd2=input$dropdown2, drd4=input$dropdown4,
      drd5=input$dropdown5, drd6=input$dropdown6, drd7=input$dropdown7, drd8=input$dropdown8, 
      sid2=input$slide2, sid3=input$slide3, sid4=input$slide4, sid0=input$slide,
      rt.tb=F
  ){

    if(is.null(meth) ){
      return(NULL)
    }else if(meth=="Ratio of Median"){
      weibTemp <<- WeibullFit(d3=drd3, d2=drd2, d4=drd4,
         d5=drd5, d6=drd6, d7=drd7, d8=drd8, 
         s2=sid2, s3=sid3, s4=sid4, s0=sid0 )
      plotLimit <- max(weibTemp$median.surv.estiMates)*1.3
      if(is.null(weibTemp)){return(NULL)}
      RatPlot(xlim=c(0,plotLimit), ylim=c(0,plotLimit), # range of the data
              xlab=expression(paste(nu[C]," = median time for C")),
              ylab=expression(paste(nu[Rx]," = median time for Rx")), # label of two axes
              med.surv.Rx=weibTemp$median.surv.estiMates[1,], # median survival time for the Rx (sub)groups
              med.surv.C=weibTemp$median.surv.estiMates[2,], # median survival time for the C (sub)groups
              angle.1=as.numeric(weibTemp$Ratio.CI[2,]), # angle.1 (a vector) = the upper bounds of simultaneous CIs for the ratios
              angle.2=as.numeric(weibTemp$Ratio.CI[1,]), # angle.2 (a vector) = the lower bounds of simultaneous CIs for the ratios
              color=c(2,3,4), #the color of the lines/arcs
              pch=c(21,22,23), #the shape of points (to represent the median survival times)
              legend=c(paste(":       ", "g-"), 
                       paste(":       ", "g+"),
                       paste(": {", "g-", "," , "g+", "}")), # the legend of the plot
              axis.1=c(expression(paste(nu[C],"-")),expression(paste(nu[C],"+")), expression(paste(nu[C],"*"))),
              axis.2=c(expression(paste(nu[Rx],"-")),expression(paste(nu[Rx],"+")), expression(paste(nu[Rx],"*")))# labels of the two axes
        )
    }else if(grepl('Relative Response',meth, fixed=T)){
    ####convert colnames to function default #cn0 <- c(input$dropdown6, input$dropdown3, input$dropdown5)
      cn0<-c(drd6, drd3, drd5); cn2 <- c("Y", "Trt", "M")
      dd <- indataset[[1]][,cn0]; colnames(dd)<-cn2
      if(grepl('Logistic', meth, fixed=T)){
        mod <- logistic(dd, 0)
        fit <- logistic.log.RR.dt.dim3(theta=mod$theta, p.neg=sid4/100, param.var=mod$fit.theta.Var)
        log.rr <- c(fit$log.RR.neg, fit$log.RR.pos, fit$log.RR.mix)
        sim <- log.RR.simu.CI(seed=123456, alpha=sid0/100, log.RR=log.rr, var.log.RR=fit$var.log.RR)
      }else{
        mod <- loglinear(dd, 0)
        fit <- loglinear.log.RR.dt.dim3(theta=mod$theta, p.neg=sid4/100, param.var=mod$fit.theta.Var)
        log.rr <- c(fit$log.RR.neg, fit$log.RR.pos, fit$log.RR.mix)
        sim <- log.RR.simu.CI(seed=123456, alpha=sid0/100, log.RR=log.rr, var.log.RR=fit$var.log.RR)
      }
      if(rt.tb){
        ot <- do.call('rbind', sim)
        ot <- rbind(round(log.rr, 3), round(ot, 3))
        nch0 <- max(nchar(ot), na.rm=T)
        gp <- fill.sp(c('g-', 'g+', 'g.mix'), nch0)
        ot <- rbind(gp, ot)
        rownames(ot)[1:2] <- c("group", "RR")
        ot <- cbind(fill.sp(rownames(ot), max(nchar(rownames(ot)), na.rm=T)), ot)
        ot <- paste(apply(ot, 1, paste, collapse= '\t\t'), collapse='\t\t\n')
        return(ot)
      }else{
        lineP(rr=exp(log.rr), s=sim, ylab=meth)#line plot
      }
    }else if(meth=="Difference of Median"){
      weibTemp <<- WeibullFit(d3=drd3, d2=drd2, d4=drd4,
         d5=drd5, d6=drd6, d7=drd7, d8=drd8, 
         s2=sid2, s3=sid3, s4=sid4, s0=sid0 )
      plotLimit <- max(weibTemp$median.surv.estiMates)*1.3
      if(is.null(weibTemp)){return(NULL)}
      DifPlot(xlim=c(0,plotLimit),ylim=c(0,plotLimit), # range of the data
                xlab=expression(paste(nu[C]," = median time for C")),
                ylab=expression(paste(nu[Rx]," = median time for Rx")), # label of two axes
                med.surv.Rx=weibTemp$median.surv.estiMates[1,], # median survival time for the Rx (sub)groups
                med.surv.C=weibTemp$median.surv.estiMates[2,], # median survival time for the C (sub)groups
                UCI=weibTemp$Diff.CI[2,], # UCI (a vector) = the upper bounds of simultaneous CIs for the differences
                LCI=weibTemp$Diff.CI[1,], # LCI (a vector) = the lower bounds of simultaneous CIs for the differences
                color=c(2,3,4), # the color of the lines/arcs
                pch=c(21,22,23), # the shape of points (to represent the median survival times)
                legend=c(paste(":    ", "g-"), 
                         paste(":    ", "g+"),
                         paste(": {", "g-", "," , "g+", "}")), # the legend of the plot
                axis.1=c(expression(paste(nu[C],"-")),expression(paste(nu[C],"+")), expression(paste(nu[C],"*"))),
                axis.2=c(expression(paste(nu[Rx],"-")),expression(paste(nu[Rx],"+")), expression(paste(nu[Rx],"*")))# labels of the two axes
        )
    }else{return(NULL)}
  } 

}

#^^^^^^^^Original name: Weibull_median_RatDiff_2grp.R^^^^^^^^^^^#
#-----------------------------------------#
# Weibull_Median_RatDiff_2grp.R includes 4 functions: 
#  	1. weib: weibull fitting for 2 group case
#   2. median.surv: calculate median survival for each subgroup and combined group
#  	3. var.median.surv.dt.dim6
#  	4. log.ratio.median.surv.dt.dim3
#   5. log.ratio.simu.CI: CI computation is based on log-scale ratio 
#	  6. dif.median.surv.dt.dim3
#   7. dif.median.simu.CI
#
#-----------------------------------------#

weib <- function(formula, data){
  # data contain the following columns: Y, Delta, Trt, M
  # Trt=1: treatment; Trt=0: control
  fit.weib <- weibreg(formula=formula, data=data)
  #fit.weib <- weibreg(Surv(Y,Delta) ~ relevel(as.factor(Trt),"1")+as.factor(M)+relevel(as.factor(Trt),"1")*as.factor(M), data=data)
  fit.LogCoeff <- fit.weib$coef
  fit.Coeff <- exp(fit.LogCoeff) # the coefficient in the "original" scale (exponentiated from the log scale)
  # change the order of the parameters (so that they correspond to "lambda","k","theta1","theta2","theta3")
  fit.Coeff <- fit.Coeff[c(4,5,1:3)]
  names(fit.Coeff)[1:2]<-c("Scale","Shape")
  # change the order of the var parameters accordingly 
  fit.LogVar <- fit.weib$var
  tmp1 <- fit.LogVar[1:3,1:3]
  tmp2 <- fit.LogVar[4:5,4:5]
  tmp3 <- fit.LogVar[1:3,4:5]
  
  fit.LogVar.2 <- rbind(cbind(tmp2,t(tmp3)),cbind(tmp3,tmp1))
  # delta method to compute the Variance matrix for the parameter estimates in the "original" scale
  fit.Var <- diag(fit.Coeff) %*% fit.LogVar.2 %*% diag(fit.Coeff) 
  
  return(list(log.coef=fit.LogCoeff,log.var=fit.LogVar,coef=fit.Coeff,var=fit.Var))
}



# calculate the median survival time for each subgroup and the combined groups
median.surv <- function(lambda,k,theta,p.neg){
  
  # lambda: Weibull distribution's scale parameter 
  # k: Weibull distribution's shape parameter 
  # theta=(theta1,theta2,theta3): =exp(beta) where beta is the slope parameter in the Cox PH model
  # p.neg: Prob(g-) (assume Prob(g-|C)=Prob(g-|Rx))
  
  median.C.neg <- as.numeric(lambda*(log(2))^(1/k))
  median.C.pos <- as.numeric(lambda*(log(2)/theta[2])^(1/k))
  median.Rx.neg <- as.numeric(lambda*(log(2)/theta[1])^(1/k))  
  median.Rx.pos <- as.numeric(lambda*(log(2)/(theta[1]*theta[2]*theta[3]))^(1/k))
  
  f.C <- function(t){
    f=p.neg*exp(-(t/lambda)^k)+(1-p.neg)*exp(-theta[2]*(t/lambda)^k)-0.5
    return(f)
  }
  f.Rx <- function(t){
    f=p.neg*exp(-theta[1]*(t/lambda)^k)+(1-p.neg)*exp(-theta[1]*theta[2]*theta[3]*(t/lambda)^k)-0.5
    return(f)
  }
  median.C <- uniroot(f.C,c(max(median.C.neg,median.C.pos),min(median.C.neg,median.C.pos)))$root
  median.Rx <- uniroot(f.Rx,c(max(median.Rx.neg,median.Rx.pos),min(median.Rx.neg,median.Rx.pos)))$root
  
  return(list(median.C.neg=median.C.neg,median.C.pos=median.C.pos,
              median.Rx.neg=median.Rx.neg,median.Rx.pos=median.Rx.pos,
              median.C=median.C,median.Rx=median.Rx))
  
}

# use the delta method for implicitly defined RVs
# reference paper: Jacques Benichou and Mitchell Gail (The American Statistician, 1989, Vol.43, No.1) 

var.median.surv.dt.dim6 <- function(med.C.neg, med.C.pos, med.Rx.neg, med.Rx.pos, med.C, med.Rx,
                                    p.neg,lambda,k,theta,param.var){
  

  # calculate the variance estimate for the median survival time estimate of all group
  # using delta method for implicitly defined random variables
  
  J.11 <- (-k*exp(-(med.C.neg/lambda)^k)*(med.C.neg/lambda)^(k-1))/lambda 
  J.22 <- (-k*theta[2]*exp(-theta[2]*(med.C.pos/lambda)^k)*(med.C.pos/lambda)^(k-1))/lambda 
  J.33 <- (-k*theta[1]*exp(-theta[1]*(med.Rx.neg/lambda)^k)*(med.Rx.neg/lambda)^(k-1))/lambda 
  J.44 <- (-k*theta[1]*theta[2]*theta[3]*exp(-theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k)*(med.Rx.pos/lambda)^(k-1))/lambda 
  J.55 <- (-p.neg*k*exp(-(med.C/lambda)^k)*(med.C/lambda)^(k-1)
           -(1-p.neg)*k*theta[2]*exp(-theta[2]*(med.C/lambda)^k)*(med.C/lambda)^(k-1))/lambda 
  J.66 <- (-p.neg*k*theta[1]*exp(-theta[1]*(med.Rx/lambda)^k)*(med.Rx/lambda)^(k-1)
           -(1-p.neg)*k*theta[1]*theta[2]*theta[3]*exp(-theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k)*(med.Rx/lambda)^(k-1))/lambda 
  
  J.mat <- diag(c(J.11, J.22, J.33, J.44, J.55, J.66))  
  J.mat.inv <- solve(J.mat)
  
  H.mat <- mat.or.vec(6,5)
  H.mat[1,] <- c(exp(-(med.C.neg/lambda)^k)*k*(med.C.neg/lambda)^(k-1)*(med.C.neg/lambda^2), 
                 -exp(-(med.C.neg/lambda)^k)*(med.C.neg/lambda)^k*log(med.C.neg/lambda),
                 0,
                 0,
                 0)
  H.mat[2,] <- c(exp(-theta[2]*(med.C.pos/lambda)^k)*theta[2]*k*(med.C.pos/lambda)^(k-1)*(med.C.pos/lambda^2), 
                 -exp(-theta[2]*(med.C.pos/lambda)^k)*theta[2]*(med.C.pos/lambda)^k*log(med.C.pos/lambda),
                 0,
                 -exp(-theta[2]*(med.C.pos/lambda)^k)*(med.C.pos/lambda)^k,
                 0)
  H.mat[3,] <- c(exp(-theta[1]*(med.Rx.neg/lambda)^k)*theta[1]*k*(med.Rx.neg/lambda)^(k-1)*(med.Rx.neg/lambda^2),
                 -exp(-theta[1]*(med.Rx.neg/lambda)^k)*theta[1]*(med.Rx.neg/lambda)^k*log(med.Rx.neg/lambda),   
                 -exp(-theta[1]*(med.Rx.neg/lambda)^k)*(med.Rx.neg/lambda)^k,
                 0,
                 0)
  H.mat[4,] <- c(exp(-theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k)*theta[1]*theta[2]*theta[3]*k*(med.Rx.pos/lambda)^(k-1)*(med.Rx.pos/lambda^2),
                 -exp(-theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k)*theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k*log(med.Rx.pos/lambda),   
                 -exp(-theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k)*theta[2]*theta[3]*(med.Rx.pos/lambda)^k,
                 -exp(-theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k)*theta[1]*theta[3]*(med.Rx.pos/lambda)^k,
                 -exp(-theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k)*theta[1]*theta[2]*(med.Rx.pos/lambda)^k)
  
  H.mat[5,] <- c(p.neg*exp(-(med.C/lambda)^k)*k*(med.C/lambda)^(k-1)*(med.C/lambda^2)
                 +(1-p.neg)*exp(-theta[2]*(med.C/lambda)^k)*theta[2]*k*(med.C/lambda)^(k-1)*(med.C/lambda^2), 
                 -p.neg*exp(-(med.C/lambda)^k)*(med.C/lambda)^k*log(med.C/lambda)
                 -(1-p.neg)*exp(-theta[2]*(med.C/lambda)^k)*theta[2]*(med.C/lambda)^k*log(med.C/lambda),
                 0,
                 -(1-p.neg)*exp(-theta[2]*(med.C/lambda)^k)*(med.C/lambda)^k,
                 0)
  
  H.mat[6,] <- c(p.neg*exp(-theta[1]*(med.Rx/lambda)^k)*theta[1]*k*(med.Rx/lambda)^(k-1)*(med.Rx/lambda^2)
                 +(1-p.neg)*exp(-theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k)*theta[1]*theta[2]*theta[3]*k*(med.Rx/lambda)^(k-1)*(med.Rx/lambda^2),
                 -p.neg*exp(-theta[1]*(med.Rx/lambda)^k)*theta[1]*(med.Rx/lambda)^k*log(med.Rx/lambda) 
                 -(1-p.neg)*exp(-theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k)*theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k*log(med.Rx/lambda),   
                 -p.neg*exp(-theta[1]*(med.Rx/lambda)^k)*(med.Rx/lambda)^k-(1-p.neg)*exp(-theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k)*theta[2]*theta[3]*(med.Rx/lambda)^k,
                 -(1-p.neg)*exp(-theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k)*theta[1]*theta[3]*(med.Rx/lambda)^k,
                 -(1-p.neg)*exp(-theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k)*theta[1]*theta[2]*(med.Rx/lambda)^k)
  
  median.var <- J.mat.inv %*% H.mat %*% param.var %*% t(H.mat) %*% t(J.mat.inv)
  
  return(list(median.var = median.var))
}

#-------#
# ratio #
#-------#
# Use ratio (med.Rx.neg/med.C.neg, med.Rx.pos/med.C.pos, med.Rx/med.C) as the efficacy measure

log.ratio.median.surv.dt.dim3 <- function(median.surv,median.var){
  
  median.C.neg <- median.surv["median.C.neg"]
  median.C.pos <- median.surv["median.C.pos"]
  median.Rx.neg <- median.surv["median.Rx.neg"]
  median.Rx.pos <- median.surv["median.Rx.pos"]
  median.C <- median.surv["median.C"]
  median.Rx <- median.surv["median.Rx"]
  
  # ratio for M-, M+, combined
  gamma.neg <- median.Rx.neg/median.C.neg
  gamma.pos <- median.Rx.pos/median.C.pos
  gamma.mix <- median.Rx/median.C
  
  # log(ratio) for M-, M+, combined
  log.gamma.neg <- log(median.Rx.neg/median.C.neg)
  log.gamma.pos <- log(median.Rx.pos/median.C.pos)
  log.gamma.mix <- log(median.Rx/median.C)
  
  
  D.log.gamma <- t(matrix(c(-1/median.C.neg, 0, 1/median.Rx.neg, 0, 0, 0,
                            0, -1/median.C.pos, 0, 1/median.Rx.pos, 0, 0, 
                            0, 0, 0, 0,  -1/median.C, 1/median.Rx),3, 6, byrow =TRUE))
  # var.log.gamma is the variance/covariance matrix for log(ratio)
  var.log.gamma <- t(D.log.gamma) %*% median.var %*% (D.log.gamma)  
  
  return(list(log.ratio.neg=log.gamma.neg, log.ratio.pos=log.gamma.pos, log.ratio.mix=log.gamma.mix, 
              ratio.neg=gamma.neg, ratio.pos=gamma.pos, ratio.mix=gamma.mix,
              var.log.ratio=var.log.gamma))
}

# CI computation is based on log-scale ratio
log.ratio.simu.CI <- function(seed,alpha,log.ratio,var.log.ratio){
  # qmvnorm is simulation based, provide a seed
  set.seed(seed)
  log.CIs <- qmvnorm(1-alpha, corr=cov2cor(var.log.ratio), tail = "both")$quantile*sqrt(diag(var.log.ratio))
  
  ratio.CI.lbd <- exp(log.ratio - log.CIs)
  ratio.CI.ubd <- exp(log.ratio + log.CIs)
  
  return(list(ratio.CI.lbd=ratio.CI.lbd,ratio.CI.ubd=ratio.CI.ubd))
}



#------------#
# difference #
#------------#
# Use difference (med.Rx.neg-med.C.neg, med.Rx.pos-med.C.pos, med.Rx-med.C) as the efficacy measure

dif.median.surv.dt.dim3 <- function(median.surv, median.var){
	
	median.C.neg <- median.surv["median.C.neg"]
	median.C.pos <- median.surv["median.C.pos"]
	median.Rx.neg <- median.surv["median.Rx.neg"]
	median.Rx.pos <- median.surv["median.Rx.pos"]
	median.C <- median.surv["median.C"]
	median.Rx <- median.surv["median.Rx"]
	
	delta.neg <- median.Rx.neg-median.C.neg
	delta.pos <- median.Rx.pos-median.C.pos
	delta.mix <- median.Rx-median.C
	
	D.delta <- t(matrix(c(-1, 0, 1, 0, 0, 0,
							0, -1, 0, 1, 0, 0, 
							0, 0, 0, 0,	-1, 1),3, 6, byrow =TRUE))
	var.delta <- t(D.delta) %*% median.var %*% (D.delta)  
	
	return(list(dif.neg=delta.neg, dif.pos=delta.pos, dif.mix=delta.mix, var.dif=var.delta))
}


# CI computation is based on difference of median

dif.median.simu.CI <- function(seed, alpha, dif.median, var.dif.median){
  # qmvnorm is simulation based, provide a seed
  set.seed(seed)
  CIs <- qmvnorm(1-alpha, corr=cov2cor(var.dif.median), tail = "both")$quantile*sqrt(diag(var.dif.median))
  
  dif.median.CI.lbd <- dif.median - CIs
  dif.median.CI.ubd <- dif.median + CIs
  
  return(list(dif.median.CI.lbd=dif.median.CI.lbd, dif.median.CI.ubd=dif.median.CI.ubd))
}



#^^^^^^^^Original name: App_RateDiffPlot.R^^^^^^^^^^^#
#----------------------#
#-- M&M PLOT for App --#
#----------------------#

# the following is the function for generating the ratio M&M plot
RatPlot <- function(xlim=c(0,20),ylim=c(0,20), # range of the data
                    xlab=expression(paste(m[C]," = median OS for C")),
                    ylab=expression(paste(m[Rx]," = median OS for Rx")), # label of two axes
                    med.surv.Rx=c(10,18,14), # median survival time for the Rx (sub)groups
                    med.surv.C=c(8,10,9), # median survival time for the C (sub)groups
                    angle.1=c(1.65,2.5,2.16), 
                    angle.2=c(0.85,1.1,0.96), 
                    color=c(2,3,4), # the color of the lines/arcs
                    pch=c(21,22,23), # the shape of points (to represent the median survival times)
                    legend=c(": g+ group",": g- group", ": {g+,g-} combined group"), # the legend of the plot
                    axis.1=c(expression(paste(m[C],"+")),expression(paste(m[C],"-")), expression(paste(m[C],"*"))),
                    axis.2=c(expression(paste(m[Rx],"+")),expression(paste(m[Rx],"-")), expression(paste(m[Rx],"*")))# labels of the two axes
){
  
  # calculate the median survival ratios for the (sub)groups  			
  ratio = med.surv.Rx/med.surv.C
  
  par(mar=c(5, 5, 2, 2), cex.lab=1.2, cex.axis=1.2)
  # generate a blank plot with axes labels
  plot(xlim,ylim,type="n",xaxt='n',yaxt='n',xaxs="i",yaxs="i",
       ylab=ylab,
       xlab=xlab,asp=1)
  
  # add the 45 degree diagonal line
  abline(a=0,b=1,lwd=1,lty=2,col="grey") 
  # dots for the observed median survival
  lines(med.surv.C,med.surv.Rx,type="p",lwd=1.5,pch=pch,col="black",bg=color)
  for(i in 1:length(ratio)){
    # add lines and arcs
    abline(a=0,b=ratio[i],lwd=1.5,lty=1,col=color[i])
    draw.arc(0,0,radius=sqrt(med.surv.C[i]^2+med.surv.Rx[i]^2),
             angle1=atan(angle.1[i]),angle2=atan(angle.2[i]),n=0.05,col=color[i],
             lwd=1.5,xaxs="i",yaxs="i",lend=1)
    
    # guidelines (vertical/horizontal) for the observed median survival
    lines(c(med.surv.C[i],med.surv.C[i]),c(0,med.surv.Rx[i]),type="l",lwd=1,lty=3,col=color[i])
    lines(c(0,med.surv.C[i]),c(med.surv.Rx[i],med.surv.Rx[i]),type="l",lwd=1,lty=3,col=color[i])
  }
  
  #axis labels
  axis(1, at=med.surv.C, lab=axis.1)
  axis(2, at=med.surv.Rx, lab=axis.2)
  
  #legend
  legend("topright", legend, cex=1.2,col="black",pt.bg=color, pch=pch)
  
}

lineP <- function(rr=log.rr, s, ylab=meth){ #fit and sim result from logistic or loglinear for RR
    #setup 
    leg <- c(paste(":       ", "g-"),  
             paste(":       ", "g+"), paste(": {", "g-", "," , "g+", "}"))
    rr.x <- 1:length(rr)
    rr.lo<- s$RR.CI.lbd
    rr.hi<- s$RR.CI.ubd
    ylim = range(c(rr.lo, rr.hi), rm.na=T) + c(-0.5, 0.5)
    xlim = c(0,5)
    xlab = ''
    color=c(2,3,4)
    pch=c(21,22,23)
  
    par(mar=c(5, 5, 2, 2), cex.lab=1.2, cex.axis=1.2)
    plot(rr~rr.x, xlim=xlim, ylim=ylim, axes=F, ylab=ylab, xlab=xlab, 
         col='black', bg=color, pch=pch, cex=3)
    arrows(rr.x, rr.lo, rr.x, rr.hi, code=3, angle=90, col=color, lwd=4)
      
    #legend
    legend("topright", leg, cex=1.2, pch=pch, col="black", pt.bg=color)
    box()
    axis(2)
}


DifPlot <- function(xlim=c(0,20),ylim=c(0,20), # range of the data
                    xlab=expression(paste(m[C]," = median OS for C")),
                    ylab=expression(paste(m[Rx]," = median OS for Rx")), # label of two axes
                    med.surv.Rx=c(10,18,14), # median survival time for the Rx (sub)groups
                    med.surv.C=c(8,10,9), # median survival time for the C (sub)groups
                    # radius=c(sqrt(20),sqrt(25),sqrt(30)), # the length of the arcs
                    UCI=c(3,10, 7), # UCI (a vector) = the upper bounds of simultaneous CIs for the differences
                    LCI=c(1, 6, 3), # LCI (a vector) = the lower bounds of simultaneous CIs for the differences
                    color=c(2,3,4), # the color of the lines/CI
                    pch=c(21,22,23), # the shape of points (to represent the median survival times)
                    legend=c(": g+ group",": g- group", ": {g+,g-} combined group"), # the legend of the plot
                    axis.1=c(expression(paste(m[C],"+")),expression(paste(m[C],"-")), expression(paste(m[C],"*"))),
                    axis.2=c(expression(paste(m[Rx],"+")),expression(paste(m[Rx],"-")), expression(paste(m[Rx],"*"))) # labels of the two axes
){
  # calculate the differences of median survival for the (sub)groups				
  dif  = med.surv.Rx-med.surv.C
  
  par(mar=c(5, 5, 2, 1), cex.lab=1.2, cex.axis=1.2)
  # generate a blank plot with axes labels
  plot(xlim,ylim,type="n",xaxt='n',yaxt='n',xaxs="i",yaxs="i",
       ylab=ylab,
       xlab=xlab,asp=1)
  
  # add the 45 degree diagonal line
  abline(a=0,b=1,lwd=1,lty=2,col="grey") 
  # dots for the observed median survival
  lines(med.surv.C, med.surv.Rx,type="p",lwd=1.5,pch=pch,col="black",bg=color)
  
  for(i in 1:length(dif)){
    # add lines for prognostic and predictive
    abline(a=dif[i],b=1,lwd=1.5,lty=1,col=color[i])
    
    # CI for the difference
    lines(c((med.surv.C[i]-(UCI[i]-dif[i])/2), (med.surv.C[i]+(dif[i]-LCI[i])/2)),
          c((med.surv.Rx[i]+(UCI[i]-dif[i])/2), (med.surv.Rx[i]-(dif[i]-LCI[i])/2)), type='l', lwd=1.5, col=color[i], lend=1)
    
    # guidelines for the observed median survival
    lines(c(med.surv.C[i],med.surv.C[i]),c(0,med.surv.Rx[i]),type="l",lwd=1,lty=3,col=color[i])
    lines(c(0,med.surv.C[i]),c(med.surv.Rx[i],med.surv.Rx[i]),type="l",lwd=1,lty=3,col=color[i])
  }
  
  #axis labels
  axis(1, at=med.surv.C, lab=axis.1)
  axis(2, at=med.surv.Rx, lab=axis.2)
    
  #legend
  legend("topright", legend, cex=1.2,col="black",pt.bg=color, pch=pch)
  
}


#-----------------------------------------#
# Binary outcome_logRR.R includes 5 functions: 
#	1. logistic: run logistic model
#     2. logistic.log.RR.dt.dim3: delta method for log RR (logistic model)
#     3. loglinear: run loglinear model
#	4. loglinear.log.RR.dt.dim3:: delta method for log RR (loglinear model)
#	5. log.RR.simu.CI: CI computation
#
# jch 01/01/2018
#-----------------------------------------#

#install required libraries
if(TRUE){
  library(msm) # for deltamethod

  list.of.packages <- c("eha", "survival", "rootSolve", "mvtnorm", "plotrix")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)

  library(eha)
  library(survival)
  library(rootSolve)
  library(mvtnorm)
  library(plotrix)
}


logistic <- function(data, m){
  # m: the total number of marker groups
  # data are generated from "dataLogistic" function
  
  fit.logit <- glm(Y ~ Trt+as.factor(M)+Trt*as.factor(M), family=binomial(logit), data=data)
  beta <- fit.logit$coef	# beta1, ... , beta4
  theta <- exp(beta)	# theta1, ... ,theta4
  
  fit.beta.Var <- vcov(fit.logit)
  #fit.theta.Var <- diag(theta) %*% fit.beta.Var %*% diag(theta) 
  
  #-- applying delta method : using the function "deltamethod" 
  fit.theta.Var <- deltamethod ( list( ~exp(x1), ~exp(x2), ~exp(x3),~exp(x4)), beta, fit.beta.Var, ses=FALSE)		  
  
  return(list(beta=beta, fit.beta.Var=fit.beta.Var, theta=theta, fit.theta.Var=fit.theta.Var))
}

# Use log(response.Rx.neg/response.C.neg, response.Rx.pos/response.C.pos, response.Rx/response.C) as the efficacy measure

logistic.log.RR.dt.dim3 <- function(theta, p.neg, param.var){
	
	# response rate 
	response.C.neg = theta[1] / (1 + theta[1])
	response.C.pos = theta[1] * theta[3] / (1 + theta[1] * theta[3])
	response.Rx.neg = theta[1] * theta[2] / (1 + theta[1] * theta[2])
	response.Rx.pos = theta[1] * theta[2] * theta[3] * theta[4] / (1 + theta[1] * theta[2] * theta[3] * theta[4])
	
	response.C = p.neg * (theta[1] / (1 + theta[1])) + (1-p.neg) * (theta[1] * theta[3] / (1 + theta[1] * theta[3]))
	response.Rx = p.neg * (theta[1] * theta[2] / (1 + theta[1] * theta[2])) + (1-p.neg) * (theta[1] * theta[2] * theta[3] * theta[4] / (1 + theta[1] * theta[2] * theta[3] * theta[4]))
	
	RR.neg <- response.Rx.neg/response.C.neg
	RR.pos <- response.Rx.pos/response.C.pos
	RR.mix <- response.Rx/response.C
	
	# log relative response
	log.RR.neg <- log(response.Rx.neg/response.C.neg)
	log.RR.pos <- log(response.Rx.pos/response.C.pos)
	log.RR.mix <- log(response.Rx/response.C)
	
	# compute the variance estimates for the relative response
	
	#D.RR <- matrix(c(0,	1/theta[2],	0,	0,
	#					0,	1/theta[2],	0,	1/theta[4], 
	#					0,	1/theta[2],	(p.neg*(1-p.neg)*(theta[4]-1))/((p.neg+(1-p.neg)*theta[3])*(p.neg+((1-p.neg)*theta[3]*theta[4]))),	
	#					((1-p.neg)*theta[3])/(p.neg+(1-p.neg)*theta[3]*theta[4])),3, 4, byrow =TRUE)
	#var.log.RR <- D.RR %*% param.var %*% t(D.RR)  	
	
	## extra variable p.neg build up the formula as a string, and convert to a formula.
	form <- sprintf("~ log( ( %f*x2 / (1+x1*x2) + (1- %f)*x2*x3*x4 / (1+x1*x2*x3*x4) ) / ( %f / (1+x1) + (1- %f)*x3/ (1+x1*x3) ) )", p.neg, p.neg, p.neg, p.neg)
	var.log.RR <- deltamethod ( list( ~log( (x2*(1+x1)) / (1+x1*x2) ), ~log((x2*x4*(1+x1*x3)) / (1+x1*x2*x3*x4)), as.formula(form)), theta, param.var, ses=FALSE)
	
	return(list(response.C.neg=response.C.neg, 		response.C.pos=response.C.pos,
					response.Rx.neg=response.Rx.neg,	response.Rx.pos=response.Rx.pos,
					response.C=response.C, response.Rx=response.Rx,
					RR.neg=RR.neg, RR.pos=RR.pos, RR.mix=RR.mix,
					log.RR.neg=log.RR.neg, log.RR.pos=log.RR.pos, log.RR.mix=log.RR.mix, var.log.RR=var.log.RR))
}



loglinear <- function(data, m){
  # m: the total number of marker groups
  # data are generated from "dataLoglinear" function
  
  fit.loglinear <- glm(Y ~ Trt+as.factor(M)+Trt*as.factor(M), family=poisson, data=data)
  beta <- fit.loglinear$coef	# beta1, ... , beta4
  theta <- exp(beta)	# theta1, ... ,theta4
  
  fit.beta.Var <- vcov(fit.loglinear)
  #fit.theta.Var <- diag(theta) %*% fit.beta.Var %*% diag(theta) 
  
  #-- applying delta method : using the function "deltamethod" 
  fit.theta.Var <- deltamethod ( list( ~exp(x1), ~exp(x2), ~exp(x3),~exp(x4)), beta, fit.beta.Var, ses=FALSE)
  
  return(list(beta=beta, fit.beta.Var=fit.beta.Var, theta=theta, fit.theta.Var=fit.theta.Var))
}

# Use log(response.Rx.neg/response.C.neg, response.Rx.pos/response.C.pos, response.Rx/response.C) as the efficacy measure

loglinear.log.RR.dt.dim3 <- function(theta, p.neg, param.var){
	
	# response rate 
	response.C.neg = theta[1]
	response.C.pos = theta[1]*theta[3]
	response.Rx.neg = theta[1]*theta[2]
	response.Rx.pos = theta[1]*theta[2]*theta[3]*theta[4]
	
	response.C = p.neg*theta[1] + (1-p.neg)*theta[1]*theta[3]
	response.Rx = p.neg*theta[1]*theta[2] + (1-p.neg)*theta[1]*theta[2]*theta[3]*theta[4]
	
	RR.neg <- response.Rx.neg/response.C.neg
	RR.pos <- response.Rx.pos/response.C.pos
	RR.mix <- response.Rx/response.C
	
	# log relative response
	log.RR.neg <- log(response.Rx.neg/response.C.neg)
	log.RR.pos <- log(response.Rx.pos/response.C.pos)
	log.RR.mix <- log(response.Rx/response.C)
	
	# compute the variance estimates for the relative response

	#D.RR <- matrix(c(0,	1/theta[2],	0,	0,
	#					0,	1/theta[2],	0,	1/theta[4], 
	#					0,	1/theta[2],	(p.neg*(1-p.neg)*(theta[4]-1))/((p.neg+(1-p.neg)*theta[3])*(p.neg+((1-p.neg)*theta[3]*theta[4]))),	
	#					((1-p.neg)*theta[3])/(p.neg+(1-p.neg)*theta[3]*theta[4])),3, 4, byrow =TRUE)
	#var.log.RR <- D.RR %*% param.var %*% t(D.RR)  	
	
	## extra variable p.neg build up the formula as a string, and convert to a formula.
	form <- sprintf("~ log(( %f * x2 + (1- %f) * x2 * x3 * x4)/( %f + (1- %f) * x3))", p.neg, p.neg, p.neg, p.neg)
	var.log.RR <- deltamethod ( list( ~log(x2), ~log(x2 * x4), as.formula(form)), theta, param.var, ses=FALSE)
	
	return(list(response.C.neg=response.C.neg, 		response.C.pos=response.C.pos,
				response.Rx.neg=response.Rx.neg,	response.Rx.pos=response.Rx.pos,
				response.C=response.C, response.Rx=response.Rx,
				RR.neg=RR.neg, RR.pos=RR.pos, RR.mix=RR.mix,
				log.RR.neg=log.RR.neg, log.RR.pos=log.RR.pos, log.RR.mix=log.RR.mix, var.log.RR=var.log.RR))
}


# CI computation
log.RR.simu.CI <- function(seed, alpha, log.RR, var.log.RR){
  # qmvnorm is simulation based, provide a seed
  set.seed(seed)
  log.CIs <- qmvnorm(1-alpha, corr=cov2cor(var.log.RR), tail = "both")$quantile*sqrt(diag(var.log.RR))

  RR.CI.lbd <- exp(log.RR - log.CIs)
  RR.CI.ubd <- exp(log.RR + log.CIs)
  
  return(list(RR.CI.lbd=RR.CI.lbd, RR.CI.ubd=RR.CI.ubd))
}


#######try functions##############
if(F){

  #load data
  dd <- indataset[[1]]
  head(dd)

  #run model
  mod1 <- logistic(dd, 1)
  fit1 <- logistic.log.RR.dt.dim3(theta=mod1$theta, p.neg=0.2, param.var=mod1$fit.theta.Var)
  sim1 <- log.RR.simu.CI(seed=123456, alpha=0.05, log.RR=fit1$log.RR.neg, var.log.RR=fit1$var.log.RR)

  mod2 <- loglinear(dd, 1)
  fit2 <- loglinear.log.RR.dt.dim3(theta=mod2$theta, p.neg=0.2, param.var=mod2$fit.theta.Var)
  sim2 <- log.RR.simu.CI(seed=123456, alpha=0.05, log.RR=fit2$log.RR.neg, var.log.RR=fit2$var.log.RR)

  #show result
  mod1
  sim1  
  mod2
  sim2



}






