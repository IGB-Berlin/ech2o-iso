
library(ggplot2)
library(grid)
library(hydroGOF)
library(plyr)
library(RColorBrewer)
library(scales)
library(Hmisc)
library(reshape)
library(ggthemes)
library(plotrix)

setwd('C:/Users/s01ks6/OneDrive - University of Aberdeen')

ext<-'LHS10'
d.sim<-paste0('RUNS/ECH2O/BB/Calibration/Results_',ext)

byjob = 1
jobs <- 1:50#; jobs<-jobs[-3]

metric1<-'KGE2012'
metric2<-'MAE'
metric3<-'RMSE'

obs.ref<- c('Streamflow','SWC_Peat','SWC_Gley.L1','SWC_Gley.L2',
            'SWC_Podzol.L1','SWC_Podzol.L2','SWC_ForestB.L1','SWC_ForestB.L2','SWC_HeatherB.L1','SWC_HeatherB.L2',
            'T_ForestA','T_ForestB','NetRad_WS1','NetRad_WS2','NetRad_WS3')

if(byjob==0){
  # Right order (not useful now though...)
  oname<-names(read.csv(file.path(d.sim,paste0(ext,'__',metric1,'.txt')),header=T))
  ind.c<-match(c('Iteration',obs.ref),oname)
  oname<-onae[ind.c]
  
  fits1<-read.csv(file.path(d.sim,paste0(ext,'_',metric1,'.txt')))[,ind.c] # remove saturation area
  fits2<-read.csv(file.path(d.sim,paste0(ext,'_',metric2,'.txt')))[,ind.c]#[,-5] # remove saturation area
  fits3<-read.csv(file.path(d.sim,paste0(ext,'_',metric3,'.txt')))[,ind.c]#[,-5] # remove saturation area
  params<-read.csv(file.path(d.sim,paste0(ext,'_parameters.txt')))
}
if(byjob==1){
  # Right order (not useful now though...)
  oname<-names(read.csv(file.path(d.sim,paste0(ext,'_job',jobs[1],'_',metric1,'.txt')),header=T))
  ind.c<-match(c('Iteration',obs.ref),oname)
  oname<-oname[ind.c]
  
  pname<-t(read.csv(file.path(d.sim,paste0(ext,'_job',jobs[1],'_parameters.txt'))))[1,]
  for(j in jobs){
    # !! REMOVE THE SKIP=1 AFTER CORRECTING OF_MIX-SCRATCH.PY !!
    tmp1<-read.csv(file.path(d.sim,paste0(ext,'_job',j,'_',metric1,'.txt')))[,ind.c]
    tmp2<-read.csv(file.path(d.sim,paste0(ext,'_job',j,'_',metric2,'.txt')))[,ind.c]
    tmp3<-read.csv(file.path(d.sim,paste0(ext,'_job',j,'_',metric3,'.txt')))[,ind.c]
    tmpp<-t(read.csv(file.path(d.sim,paste0(ext,'_job',j,'_parameters.txt')))[,-1])
    print(paste0(j,' ',nrow(tmp1),' ',nrow(tmp2),' ',nrow(tmp3),' ',nrow(tmpp)))
    if(j==jobs[1]){
      fits1<-tmp1
      fits2<-tmp2
      fits3<-tmp3
      params<-tmpp
    } else {
      fits1<-rbind(fits1,tmp1)
      fits2<-rbind(fits2,tmp2)
      fits3<-rbind(fits3,tmp3)
      params<-rbind(params,tmpp)
    }
  }
  names(fits1)<-names(fits2)<-names(fits3)<-oname
  params<-as.data.frame(params)
  names(params)<-pname
  params$Iteration<-fits1$Iteration
  fits1$Sample<-fits2$Sample<-fits3$Sample<-params$Sample<-1:nrow(fits1)
  fits1<-fits1[,c(1,ncol(fits1),2:(ncol(fits1)-1))]
  fits2<-fits2[,c(1,ncol(fits2),2:(ncol(fits2)-1))]
  fits3<-fits3[,c(1,ncol(fits3),2:(ncol(fits3)-1))]
  params<-params[,c(ncol(params)-1:0,1:(ncol(params)-2))]
}

obsnames<-oname[-1]

Nsim<-nrow(fits1)
Nout<-ifelse(byjob==0,ncol(fits1)-3,ncol(fits1)-2)
maximet<-c('KGE2012','NSE','KGE','corr')
minimet<-c('RMSE','MAE','MAEn','RMSEn')

# what to use
#####################
constr = c('Streamflow',
           'SWC_Peat','SWC_Gley.L1','SWC_Gley.L2',
           'SWC_Podzol.L1','SWC_Podzol.L2','SWC_ForestB.L1','SWC_ForestB.L2',
           'SWC_HeatherB.L1','SWC_HeatherB.L2',
           'T_ForestA','T_ForestB',
           'NetRad_WS1','NetRad_WS2','NetRad_WS3',
           #'GWD_DW1','GWD_DW2','GWD_DW3','GWD_DW4',
           #'T_SHeather','T_NHeather','ET_SHeather','ET_NHeather',
           #'dD_stream','d18O_stream',
           #'dD_GW1','d18O_GW1',
           #'dD_GW2','d18O_GW2',
           #'dD_GW3','d18O_GW3',
           #'dD_GW4','d18O_GW4',
           'SWC','Tpine','NR',#'GWD','Heather',#'Iso','Iso2',
           'Hydro')#,'Hydro2','Hydro3')
#if(ext%in%c('lhs2','MC18'))
#  constr<-c(constr,'Hydro+Iso','Hydro+Iso2','Iso','Iso2')
nc<-length(constr)
nbest<-30

# which metrics for which ? MAE for discharge, RMSE for the rest
#wmet <- c(1,2,2,2,2,2,2,2,2,2)
wmet<-list()
wmet[['Streamflow']]<-2
wmet[['SWC_Peat']]<-3
wmet[['SWC_Gley.L1']]<-wmet[['SWC_Gley.L2']]<-3
wmet[['SWC_Podzol.L1']]<-wmet[['SWC_Podzol.L2']]<-3
wmet[['SWC_ForestB.L1']]<-wmet[['SWC_ForestB.L2']]<-3
wmet[['SWC_HeatherB.L1']]<-wmet[['SWC_HeatherB.L2']]<-3
wmet[['T_ForestA']]<-wmet[['T_ForestB']]<-3
wmet[['NetRad_WS1']]<-wmet[['NetRad_WS2']]<-wmet[['NetRad_WS3']]<-3
          #1,1,1,1,#
          # 3,3,3,3,3,3,3, #swc
          # #1,1,1,1, #swc
          # 3,3, #Tpine
          # 3,3,3, #Netrad
          # 1,1,1,1, # GWD
          # #3,3,3,3, # GWD
          # 3,3,3,3, # T and ET heather
          # 3,3, #stream iso
          # 3,3, #iso DW1
          # 3,3, #iso DW2
          # 3,3, #iso DW3
          # 3,3) #iso DW4
          

rlt<-c('>','<','<')
crit2<-c('maxi','mini','mini')
im<-which(crit2=='mini' & rlt=='<')
im2<-ifelse(length(im)>1,im[1],im)
iM<-which(crit2=='maxi' & rlt=='>')
iM2<-ifelse(length(iM)>1,iM[1],iM)

if(length(which(c(metric1,metric2,metric3)%in%maximet))==0){
  rlt<-c('<','<')
  im<-c(1,2)
  im2<-1
  iM<-iM2<-integer(0)
}
if(length(which(c(metric1,metric2,metric3)%in%minimet))==0){
  rlt<-c('>','>')
  im<-im2<-integer(0)
  iM<-c(1,2)
  iM2<-1
}

fits<-evals<-fits1
# Fits array with metrics used for calibration,
# while evals uses KGE2012 (fits1) as a qualitative indicator of range
for(i in 1:Nout){
  if(byjob==0)
    fits[,i+3]<-cbind(fits1[,i+3],fits2[,i+3],fits3[,i+3])[,wmet[[obsnames[i]]]]
  if(byjob==1)
    fits[,i+2]<-cbind(fits1[,i+2],fits2[,i+2],fits3[,i+2])[,wmet[[obsnames[i]]]]
}

#####################
# (Multi-objective) criterion to follow
crit<-list()

## Single-obs
for(i in obsnames)
crit[[i]] <- 
  paste0('which(fits$',i,' ',rlt[wmet[[i]]],' quantile(fits$',i,',exc.q[',wmet[[i]],'],na.rm=T))')

## Multi-criteria

# Function
MCrit = function(char) {
  if(class(char)!='character' | length(char)==1)
    stop('error, incorrect format/class')
  else{
    lch <- length(char)
    conds<-list()
    for(i in char)
      conds[[i]]<-paste0('fits$',i,' ',rlt[wmet[[i]]],' quantile(fits$',i,',exc.q[',wmet[[i]],'],na.rm=T)')
    cmd<-paste0('which(',paste(unlist(conds),collapse=' & '),')')
  }
}

crit[['SWC']]<- MCrit(c('SWC_Peat','SWC_Gley.L1','SWC_Gley.L2','SWC_Podzol.L1','SWC_Podzol.L2',
                        'SWC_ForestB.L1','SWC_ForestB.L2','SWC_HeatherB.L1','SWC_HeatherB.L2'))
  
crit[['Tpine']] <- MCrit(c('T_ForestA','T_ForestB'))

crit[['NR']] <- MCrit(c('NetRad_WS1','NetRad_WS2','NetRad_WS3'))

crit[['Hydro']] <- MCrit(c('Streamflow',
                           'SWC_Peat','SWC_Gley.L1','SWC_Gley.L2','SWC_Podzol.L1','SWC_Podzol.L2',
                           'SWC_ForestB.L1','SWC_ForestB.L2','SWC_HeatherB.L1','SWC_HeatherB.L2',
                           'T_ForestA','T_ForestB',
                           'NetRad_WS1','NetRad_WS2','NetRad_WS3'))


#####################

good.i<- array(dim=c(nbest,nc))
quant.i<-c() 

#ind<-
for(ic in 1:nc){
  print(constr[ic])
  
  #exc.q <- c(0.5,0.5,0.5)           # the quantile to start the iteration
  exc.q <- c(0.3,0.7,0.7)           # the quantile to start the iteration
  Ngood <- Nsim                     # initial number of good simulations equals total runs
  stp<-0.001
  while (Ngood>nbest){
    
    if(length(im)>0){
      # Take the minimize criterion as reference
      exc.q[im] <- exc.q[im] - stp
      #print(exc.q)
      # In case the step is too big to isolate the nbest runs
      while(exc.q[im2]<=0){
        exc.q[im] <- exc.q[im] + stp
        stp<-stp/10
        exc.q[im] <- exc.q[im] - stp}
      # Translate for the maximize criterion
      if(length(iM)>0)
        exc.q[iM] <- 1 - exc.q[im2]
      
    } else {
      # Take the maximize criterion as reference
      exc.q[iM] <- exc.q[iM] + stp
      # In case the step is too big to isolate the nbest runs
      while(exc.q[iM2]>=1){
        exc.q[iM] <- exc.q[iM] - stp
        stp<-stp/10
        exc.q[iM] <- exc.q[iM] + stp}
    }
    
    # Run selection
    eval(parse(text=paste0('goodSims <- ',crit[[constr[ic]]])))
    Ngood <- length(goodSims)
    #cat(Ngood,'...')
    
    # In case the step reduces too much
    if(Ngood<nbest){
      if(length(im)>0){
        exc.q[im] <- exc.q[im] + stp
        if(length(iM)>0)
          exc.q[iM] <- 1 - exc.q[im2]
      } else {
        exc.q[iM] <- exc.q[iM] - stp
      }
      
      stp<-stp/10
      eval(parse(text=paste0('goodSims <- ',crit[[constr[ic]]])))
      Ngood <- length(goodSims)
    }
    #else
    #  cat('[',exc.q,',',Ngood,']...')
  }
  cat(exc.q,'\n')
  if(length(im)>0)
    quant.i<-c(quant.i,exc.q[im2])
  if(length(im)==0)
    quant.i<-c(quant.i,1-exc.q[iM2])
  good.i[,ic]<-goodSims
  #cat('\n')
  
  # Print performances in term of KGE (just an indicator)
  if(ic<=length(obs.ref))
    cat('KGEs:',signif(quantile(evals[goodSims,ic+2],probs=c(0.05,0.5,0.95),na.rm=T),digits=2),'\n')
  
  if(constr[ic]=='SWC')
    for(ic2 in 2:10)
      cat('KGEs:',obs.ref[ic2],signif(quantile(evals[goodSims,ic2+2],probs=c(0.05,0.5,0.95),na.rm=T),digits=2),'\n')
  
  if(constr[ic]=='Tpine')
    for(ic2 in 11:12)
      cat('KGEs:',obs.ref[ic2],signif(quantile(evals[goodSims,ic2+2],probs=c(0.05,0.5,0.95),na.rm=T),digits=2),'\n')
  
  if(constr[ic]=='NR')
    for(ic2 in 13:15)
      cat('KGEs:',obs.ref[ic2],signif(quantile(evals[goodSims,ic2+2],probs=c(0.05,0.5,0.95),na.rm=T),digits=2),'\n')
  
  if(constr[ic]=='Hydro')
    for(ic2 in 1:15)
      cat('KGEs:',obs.ref[ic2],signif(quantile(evals[goodSims,ic2+2],probs=c(0.05,0.5,0.95),na.rm=T),digits=2),'\n')
  
}

# Store job, iterations and samples number
colnames(good.i)<-constr
good.job<-good.iter<-good.sample<-good.i
for(ic in 1:nc){
  #good.job[,ic]<-fits$Job[good.i[,ic]]
  good.iter[,ic]<-fits$Iteration[good.i[,ic]]
  good.sample[,ic]<-fits$Sample[good.i[,ic]]
}
write.table(cbind(quant.i,1-quant.i),
            file=file.path(d.sim,paste0(ext,'_',metric1,'+',metric2,'+',metric3,'.',nbest,'best.quant.txt')),
            sep=',',col.names = c(metric1,metric2), row.names=constr,quote=F)


write.table(good.job,file=file.path(d.sim,paste0(ext,'_',metric1,'+',metric2,'+',metric3,'.',nbest,'best.jobs.txt')),
            sep=',',row.names = F, col.names=constr,quote=F)
write.table(good.iter,file=file.path(d.sim,paste0(ext,'_',metric1,'+',metric2,'+',metric3,'.',nbest,'best.iter.txt')),
            sep=',',row.names = F, col.names=constr,quote=F)
write.table(good.i,file=file.path(d.sim,paste0(ext,'_',metric1,'+',metric2,'+',metric3,'.',nbest,'best.sample.txt')),
            sep=',',row.names = F, col.names=constr,quote=F)

# Output parameters values of best runs
ind.c<-c(1:nc)
#ind.c<-27
d.out<-'RUNS/ECH2O/BB/Iso_test/Input_ensemble'
for(i in ind.c){
  print(constr[i])
  # Parameters sets
  tmp<-subset(params,Sample%in%good.sample[,i])[,-(1:2)]
  write.table(t(tmp),file=file.path(d.out,paste0(ext,'_',metric1,'+',metric2,'+',metric3,'-',constr[i],
                                              '.',nbest,'bestParams.txt')),
              sep=',',col.names = F, row.names=names(tmp),quote=F)
}

