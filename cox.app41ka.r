#    PREDICTION USING COX MODEL WITH TIME-VARYING COVARIATES 
#    JAMES F. TROENDLE   SEP, 2019
#    APPLICATION PROGRAM

#Function to perform numerical integration
#Continues on each interval until further subdivision is not helpful

numint=function(f,a,b,eps=1.0e-06,lim=10) {
  fun1=function(f,a,b,fa,fb,a0,eps,lim,fun) {
    d=(a+b)/2
    h=(b-a)/4
#    cat('Evaluate hazard at',d,'\n')
    fd=f(d)
#    cat('fd=',fd,'\n')
    a1=h*(fa+fd)
    a2=h*(fd+fb)
#    cat('fa=',fa,'\n')
#    cat('fb=',fb,'\n')
#    cat('a0=',a0,'\n')
#    cat('eps=',eps,'\n')
#    cat('lim=',lim,'\n')
    if (abs(a0-a1-a2)<eps || lim==0)
      return(a1+a2)
    else {
      return(fun(f,a,d,fa,fd,a1,eps,lim-1,fun)+
             fun(f,d,b,fd,fb,a2,eps,lim-1,fun))
    }
  }
#  cat('Evaluate hazard at',a,'\n')
  fa=f(a)
#  cat('Evaluate hazard at',b,'\n')
  fb=f(b)
  a0=((fa+fb)*(b-a))/2
#  cat('a=',a,'\n')
#  cat('b=',b,'\n')
#  cat('call fun1','\n')
  fun1(f,a,b,fa,fb,a0,eps,lim,fun1)
}





#---------------------------------------------------------------------------------------------------------------
dyn.load("addapp.dll")
# CALL ADD_APP.F TO ADD OBSERVATIONS
add_app=function(n,maxobs,num_imp_bp,imp_times_bp,
    d1,d2,d3,d4,d5,base_cov,ncov) {
x=.Fortran("addapp",n=as.integer(n),
                    maxobs=as.integer(maxobs),
                    numimpbp=as.integer(num_imp_bp),
                    imptimesbp=as.double(imp_times_bp),
                    d1=as.integer(d1),
                    d2=as.double(d2),
                    d3=as.double(d3),
                    d4=as.integer(d4),
                    d5=as.integer(d5),
                    basecov=as.double(base_cov),
                    ncov=as.integer(ncov)
)
  return(x)
}
#---------------------------------------------------------------------------------------







#---------------------------------------------------------------------------------------
dyn.load("condprob.dll")
# CALL CONDPROB.F TO CALCULATE CONDITIONAL SURVIVAL PROBS
condprob=function(nperson_future,n_future,npredtimes,index_s,future_surv_short,idvec,cond_surv_pred_st) {
x=.Fortran("condprob",nperson=as.integer(nperson_future),
                      n=as.integer(n_future),
                      npredtimes=as.integer(npredtimes),
                      index=as.integer(index_s),
                      future=as.double(future_surv_short),
                      idvec=as.integer(idvec),
                      cs=as.double(cond_surv_pred_st)
)
  return(x)
}
#---------------------------------------------------------------------------------------





#------------------------------------------------------------
#SET PARAMETERS FOR PROGRAM
#------------------------------------------------------------
#Initialize seed for random number generation
seed=99803332 
seedstart=seed
set.seed(seed)
#
#Initialize seed for random number generation in FORTRAN
iseed=765110902
iseedstart=iseed
iseed=-1*iseed
#
# FOR CALLING RAN2 IN FORTRAN
iseed2=123456789
iv=rep(0,times=32)
iy=0
#

library(survival)
library(stats)
library(gtools)
#library(plyr)
#library(pryr)
library(mvtnorm)
#library(pec)
#library(nlme)
library(lme4)
#library(merTools)
#library(JM)
#library(JMbayes)
#library(doBy)
#
#----------------------------------------------------------------------------------------
# DIM_BP IS THE NUMBER OF DIFFERENT BIOMARKERS (<=5)
#----------------------------------------------------------------------------------------
dim_bp=3
#------------------------------------------------------------------------------------
#  DATA LOCATION
#------------------------------------------------------------------------------------
infile_surv="/home/troendlj/accord/dataset_survival.Rda"
infile_bio1="/home/troendlj/accord/dataset_biomarker_HBA1C.Rda"
infile_bio2="/home/troendlj/accord/dataset_biomarker_LDL.Rda"
infile_bio3="/home/troendlj/accord/dataset_biomarker_SBP.Rda"
name_bio=list()
name_bio[[1]]="HBA1C"
name_bio[[2]]="LDL"
name_bio[[3]]="SBP"
var_name_bio=list()
var_name_bio[[1]]="hba1c"
var_name_bio[[2]]="ldl"
var_name_bio[[3]]="sbp"
#------------------------------------------------------------------------------------
#  BASELINE COV NAMES
#------------------------------------------------------------------------------------
ncov=6
base_cov_name=list()
base_cov_name[[1]]="network"
base_cov_name[[2]]="prev_cvd"
base_cov_name[[3]]="glycemia"
base_cov_name[[4]]="bp_lipid"
base_cov_name[[5]]="bp"
base_cov_name[[6]]="lipid"
base_cov_factor=rep(0,times=ncov)
base_cov_factor[1]=7
#------------------------------------------------------------
#------------------------------------------------------------
# PARAMETERS RELATED TO PREDICTION
#------------------------------------------------------------
# S = conditioning time, given survival to time S
# S MUST BE IN SET OF TIMES_BP
# TAU = prediction for survival to time S+TAU
# BP_S = value of Biomarker Process #1 at time S
# BP_TOL = tolerance on value of Biomarker Process #1 at time S
# START_EVAL_ID_NUM
# END_EVAL_ID_NUM
#------------------------------------------------------------
s=1.0
tau=2.0
bp_s=0
bp_tol=.05

start_eval_id_num=1
end_eval_id_num=1

#------------------------------------------------------------
# SPLINE PARAMETERS (MAX 8)
#------------------------------------------------------------
nspline_max=c(4,2,4)
nspline=rep(0,times=dim_bp)
spline_width=c(1,3,1)
#knot1=spline_width
#knot2=2*spline_width
#knot3=3*spline_width
#knot4=4*spline_width
#knot5=5*spline_width
#knot6=6*spline_width
#knot7=7*spline_width

#------------------------------------------------------------
#------------------------------------------------------------
# DESIGN TIMES
#------------------------------------------------------------
design_times_bp=c(0,.33333,.66666,1)
num_design_times_bp=length(design_times_bp)
#------------------------------------------------------------
# IMPUTATIUON TIMES
# BIOMARKER PROCESS IMPUTED UP TO NUM_IMP_BP TIMES
# S_INDEX = INDEX OF S IN IMP_TIMES_BP
# IMP_INC IS INCREMENT OF TIME BETWEEN IMPUTATION TIMES
# CENS_TIME IS TIME AFTER END OF IMPUTATION TIMES WHERE SUTRVIVAL RIGHT CENSORED
# LEFT_CENS_TIME IS TIME BEFORE IMPUTATION TIMES WHERE SUTRVIVAL LEFT CENSORED
#
# CONTINUE TO S+TAU
#------------------------------------------------------------
imp_inc=0.25
digits_imp_inc=2
left_cens_time=s
window_width=0.00001
imp_times_bp=seq(from=left_cens_time+window_width,to=s+tau,by=imp_inc)
imp_times_bp=sort(unique(c(imp_times_bp,s,s+tau)))
num_imp_bp=length(imp_times_bp)
s_index=which(round(imp_times_bp,digits=5)==s)
cens_time=imp_times_bp[num_imp_bp]+imp_inc
num_imp_bp_save=num_imp_bp
imp_times_bp_save=imp_times_bp
#------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# BOOTSTRAP USED FOR ESTIMATING SD OF PREDICTION
#------------------------------------------------------------------------------------
boot_num=0
#
#------------------------------------------------------------------------------------
# num_similar= NUMBER OF MOST SIMILAR PAST SUBJECTS USED TO ESTIMATE BLUP FOR FUTURE SUBJECT 
#------------------------------------------------------------------------------------
num_similar=1
#
#------------------------------------------------------------
#-----------------------------------------------------------
#Indicators of which methods will be run
# ICOX1: Use EBLUP from future subject
# ICOX2: Use EBLUP from most "similar" past subject
#-----------------------------------------------------------
icox1_ind=1
icox2_ind=0
lvcf_ind=1
jm_ind=0
#------------------------------------
#options(warn=1)
options(warn=-1)
z975=qnorm(.975,mean=0,sd=1)

#------------------------------------
#------------------------------------


#---------------------------------------------------------------------------------------------------------------
#
#  READ IN DATA
#
#---------------------------------------------------------------------------------------------------------------

dataset_bio=list()

load(infile_surv, verbose=TRUE)
dataset=dataset_primary[dataset_primary$start < dataset_primary$stop,]
nperson=length(dataset$id)
rm(dataset_primary)
print(summary(dataset))
#------------------------------------
#------------------------------------
# Select only Lipid trial participants
#  Add interaction between glycemic and lipid treatments
#------------------------------------
dataset=dataset[dataset$bp_lipid==0,]
nperson=length(dataset$id)
dataset$glylip=dataset$glycemia*dataset$lipid
cat('------------------------------------------','/n')
cat('After restricting to lipid trial','/n')
print(summary(dataset))

ncov=5
base_cov_name=list()
base_cov_name[[1]]="network"
base_cov_name[[2]]="prev_cvd"
base_cov_name[[3]]="glycemia"
base_cov_name[[4]]="lipid"
base_cov_name[[5]]="glylip"
base_cov_factor=rep(0,times=ncov)
base_cov_factor[1]=7
#------------------------------------
#------------------------------------
#------------------------------------------------------------

j=1
while (j<= dim_bp) {
  if (j==1) {
    load(infile_bio1, verbose=TRUE)
    dataset_bio[[j]]=dataset_HBA1C
    rm(dataset_HBA1C)
    dataset_bio[[j]]=dataset_bio[[j]][!is.na(dataset_bio[[j]]$obstime) & !is.na(dataset_bio[[j]]$hba1c),]
  }
  if (j==2) {
    load(infile_bio2, verbose=TRUE)
    dataset_bio[[j]]=dataset_LDL
    rm(dataset_LDL)
    dataset_bio[[j]]=dataset_bio[[j]][!is.na(dataset_bio[[j]]$obstime) & !is.na(dataset_bio[[j]]$ldl),]
  }
  if (j==3) {
    load(infile_bio3, verbose=TRUE)
    dataset_bio[[j]]=dataset_SBP
    rm(dataset_SBP)
    dataset_bio[[j]]=dataset_bio[[j]][!is.na(dataset_bio[[j]]$obstime) & !is.na(dataset_bio[[j]]$sbp),]
  }
  dataset_bio[[j]]$glylip=dataset_bio[[j]]$glycemia*dataset_bio[[j]]$lipid
  print(summary(dataset_bio[[j]]))
  print(dataset_bio[[j]][1:40,])
  j=j+1
}
#---------------------------------------------------------------------------------------------------------------
#
#  END: READ IN DATA
#
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
n=length(dataset$id)
n_bio=list()
j=1
while (j<= dim_bp) {
  n_bio[[j]]=length(dataset_bio[[j]]$id)
  j=j+1
}

#---------------------------------------------------------------------------------------------------------------
#
#  SAVE DATASETS
#
#---------------------------------------------------------------------------------------------------------------
dataset_save=dataset
dataset_bio_save=dataset_bio
#dataset_lvcf=dataset
#dataset_bio_lvcf=dataset_bio

cat('After reading in data','\n')
cat('# of subjects=',nperson,'\n')
cat('# of observations=',n,'\n')
cat('# of HBA1C observations=',n_bio[[1]],'\n')
cat('# of LDL observations=',n_bio[[2]],'\n')
cat('# of SBP observations=',n_bio[[3]],'\n')
cat('-------------------------------------------------------','\n')



#-------------------------------------------------
#-------------------------------------------------
if (icox1_ind==1 | icox2_ind==1) {
#-------------------------------------------------
# EB COX MODEL METHODS
#-------------------------------------------------

#cat('Dataset before admin censoring=','\n')
#  print(dataset[1:20,])
#cat('# of obs in dataset=',length(dataset$id),'\n')
#cat('Dataset$stop before admin censoring=','\n')
#print(summary(dataset$stop))
#cat('-------------------------------------------------------','\n')
#-----------------------------------------------------------------------------
#  Apply administrative right_censoring to Dataset for integrity of bp_curr term
#-----------------------------------------------------------------------------
#$status=dataset$status*(dataset$stop <= cens_time)+0*(dataset$stop > cens_time)
#dataset$stop=dataset$stop*(dataset$stop <= cens_time)+cens_time*(dataset$stop > cens_time)

#cat('Right censoring time=',cens_time,'\n')
#cat('Dataset after admin censoring=','\n')
#print(dataset[1:10,])
#-----------------------------------------------------------------------------
# Determine imputation times
#-----------------------------------------------------------------------------
#event_times=sort(unique(dataset$stop[dataset$status==1]))

#cat('Ordered Unique Event times','\n')
#print(event_times)

#num_imp_bp=num_imp_bp_save
#imp_times_bp=imp_times_bp_save
#num_imp_bp_new=0
#j=1
#while (j <= num_imp_bp) {
#  if (any(event_times > imp_times_bp[j] & event_times <= imp_times_bp[j]+imp_inc)) {
#    num_imp_bp_new=num_imp_bp_new+1
#    imp_times_bp[num_imp_bp_new]=imp_times_bp[j]
#  }
#  j=j+1
#}
#num_imp_bp=num_imp_bp_new
#imp_times_bp=imp_times_bp[1:num_imp_bp]
#imp_times_bp=sort(unique(c(imp_times_bp,s,s+tau)))
#num_imp_bp=length(imp_times_bp)
#s_index=which(round(imp_times_bp,digits=5)==s)

#cat('------------------------------------------------','\n')
#cat('Imputation times','\n')
#print(imp_times_bp)
#cat('------------------------------------------------','\n')

#-----------------------------------------------------------------------------
# Determine prediction times
#-----------------------------------------------------------------------------
#predtimes=sort(unique(c(imp_times_bp[2:num_imp_bp])))
#predtimes=imp_times_bp[imp_times_bp > s]
#index_s=which(predtimes==s)
#if (length(index_s)==0) {index_s=0}
#cat('s_index=',s_index,'\n')
#cat('--------------------------------------------','\n')
#cat('predtimes=',predtimes,'\n')
#cat('index of predtimes=s is',index_s,'\n')
#cat('--------------------------------------------','\n')
#cat('------------------------------------------------','\n')

#-----------------------------------------------------------------------------
# Create Multiple Observations per Subject dataset
#-----------------------------------------------------------------------------
#maxobs=n*num_imp_bp
#-------------------------------------------------
#d1=c(dataset$id,rep(0,times=maxobs-n))
#d2=c(dataset$start,rep(0,times=maxobs-n))
#d3=c(dataset$stop,rep(0,times=maxobs-n))
#d4=c(dataset$status,rep(0,times=maxobs-n))
#d5=c(rep(1,times=n),rep(0,times=maxobs-n))

#base_cov=rep(0,times=maxobs*ncov)
#dim(base_cov)=c(maxobs,ncov)
#j=1
#while (j<= ncov) {
#  temp=eval(parse(text=paste('dataset$',base_cov_name[[j]],sep='')))
#
#  cat('temp=','\n')
#  print(temp)
#  base_cov[,j]=c(temp,rep(0,times=maxobs-n))
#  j=j+1
#}
#base_cov_id=base_cov[1:nperson,]

#z=add_app(n,maxobs,num_imp_bp,imp_times_bp,
#  d1,d2,d3,d4,d5,base_cov,ncov)

# base_cov=z[[10]]
# dim(base_cov)=c(maxobs,ncov)
# dataset=data.frame(id=z[[5]],start=z[[6]],stop=z[[7]],status=z[[8]],lastobs=z[[9]],
#                    bp_curr1=rep(NA,times=maxobs),bp_curr2=rep(NA,times=maxobs),bp_curr3=rep(NA,times=maxobs),
#                    bp_curr4=rep(NA,times=maxobs),bp_curr5=rep(NA,times=maxobs),
#                    base_cov=base_cov)
# n=z[[1]]
# dataset=dataset[1:n,]
# n=length(dataset$id)
# 
# cat('Obs with start >= stop for last obs for that subject','\n')
# print(dataset[dataset$lastobs==1 & dataset$start >= dataset$stop,1:4])
# cat('# of events=',sum(dataset$status==1),'\n')
# cat('# of subjects according to those read in=',nperson,'\n')
# cat('According to dataset$lastobs # of subjects=',sum(dataset$lastobs==1),'\n')
# cat('According to dataset # of subjects=',length(unique(dataset$id)),'\n')
# cat('# of all observations=',n,'\n')
# cat('-------------------------------------------------------','\n')

#--------------------------------------------------------------------------
#
#  Create multivariate covariates for longitudinal datasets
#
#--------------------------------------------------------------------------
j=1
while (j<= dim_bp) {
  base_cov_bio=rep(0,times=n_bio[[j]]*ncov)
  dim(base_cov_bio)=c(n_bio[[j]],ncov)
  if (j==1) {bptype=rep('TYPE1',times=n_bio[[j]])}
  if (j==2) {bptype=rep('TYPE2',times=n_bio[[j]])}
  if (j==3) {bptype=rep('TYPE3',times=n_bio[[j]])}
  if (j==4) {bptype=rep('TYPE4',times=n_bio[[j]])}
  if (j==5) {bptype=rep('TYPE5',times=n_bio[[j]])}
  temp_bio=eval(parse(text=paste('dataset_bio[[j]]$',var_name_bio[[j]],sep='')))
  jj=1
  while (jj<= ncov) {
    temp=eval(parse(text=paste('dataset_bio[[j]]$',base_cov_name[[jj]],sep='')))
#  cat('temp=','\n')
#  print(temp)
    base_cov_bio[,jj]=temp
    jj=jj+1
  }
  dataset_bio[[j]]=data.frame(id=dataset_bio[[j]]$id,bp_curr=temp_bio,
                      obstime=dataset_bio[[j]]$obstime,bptype=bptype,base_cov=base_cov_bio)
  j=j+1
}
#--------------------------------------------------------------------------
#
#  Add lastobs to longitudinal dataset
#
#--------------------------------------------------------------------------

j=1
while (j <= dim_bp) {
  dataset_bio[[j]]$lastobs=rep(0,times=n_bio[[j]])
  i=1
  while (i <= n_bio[[j]]) {
    id=dataset_bio[[j]]$id[i]
    if (i==n_bio[[j]]) {
      dataset_bio[[j]]$lastobs[i]=1
    } else {
      if (id != dataset_bio[[j]]$id[i+1]) {dataset_bio[[j]]$lastobs[i]=1}
    }  
    i=i+1
  }
#  cat('j=',j,'\n')
#  cat('Longitudinal Dataset with lastobs','\n')
#  print(dataset_bio[[j]][1:50,])
#  cat('-------------------------------------------------------','\n')
  if (j==1) {
    dataset_nonmiss=dataset_bio[[j]]
  } else {
    dataset_nonmiss=rbind(dataset_nonmiss,dataset_bio[[j]])
  }
  j=j+1
}
cat('According to dataset_nonmiss # of subjects=',length(unique(dataset_nonmiss$id)),'\n')
#--------------------------------------------------------------------------
#
#  Restrict longitudinal dataset to those with survival data
#
#--------------------------------------------------------------------------
dataset_nonmiss=dataset_nonmiss[dataset_nonmiss$id %in% unique(dataset$id),]
n_nonmiss=length(dataset_nonmiss$id)
cat('After restriction to those with survival info: dataset_nonmiss # of subjects=',length(unique(dataset_nonmiss$id)),'\n')

cat('Survival Dataset','\n')
print(summary(dataset))
cat('Number of subjects with an event=',sum(dataset$status),'\n')
#cat('SUBJECT id=1819 =','\n')
#print(dataset[dataset$id==1819,])
#cat('---------------------------------------------------------','\n')
cat('Combined Longitudinal Dataset with BPTYPE','\n')
#print(dataset_nonmiss[1:50,])
print(summary(dataset_nonmiss))
dataset_nonmiss_all=unique(dataset_nonmiss$id[dataset_nonmiss$bptype=='TYPE1'])
cat('Number of subjects with at least one BPTYPE=TYPE1',length(dataset_nonmiss_all),'\n')
if (dim_bp >= 2) {dataset_nonmiss_2=unique(dataset_nonmiss$id[dataset_nonmiss$bptype=='TYPE2'])}
if (dim_bp >= 2) {dataset_nonmiss_all=dataset_nonmiss_all[dataset_nonmiss_all %in% dataset_nonmiss_2]}
if (dim_bp >= 3) {dataset_nonmiss_3=unique(dataset_nonmiss$id[dataset_nonmiss$bptype=='TYPE3'])}
if (dim_bp >= 3) {dataset_nonmiss_all=dataset_nonmiss_all[dataset_nonmiss_all %in% dataset_nonmiss_3]}
if (dim_bp >= 4) {dataset_nonmiss_4=unique(dataset_nonmiss$id[dataset_nonmiss$bptype=='TYPE4'])}
if (dim_bp >= 4) {dataset_nonmiss_all=dataset_nonmiss_all[dataset_nonmiss_all %in% dataset_nonmiss_4]}
if (dim_bp >= 5) {dataset_nonmiss_5=unique(dataset_nonmiss$id[dataset_nonmiss$bptype=='TYPE5'])}
if (dim_bp >= 5) {dataset_nonmiss_all=dataset_nonmiss_all[dataset_nonmiss_all %in% dataset_nonmiss_5]}
cat('Number of subjects with at least one of each BPTYPE=',length(dataset_nonmiss_all),'\n')
cat('---------------------------------------------------------','\n')
#cat('Biomarker data for ID=631','\n')
#print(dataset_nonmiss[dataset_nonmiss$id==631,])
#cat('---------------------------------------------------------','\n')

#---------------------------------------------------------------------------------------------------------------
#
#  CREATE SEQUENTIAL ID FROM 1 TO nperson
#
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------
nid=length(unique(dataset$id))
idvec=rep(0,times=nid)
i=1
idcount=0
idprev=0
while (i <=n) {
  idcurr=dataset$id[i]
  if (idcurr != idprev) {
    idcount=idcount+1
    idvec[idcount]=as.numeric(dataset$id[i])
  }
  dataset$id[i]=idcount
  idprev=idcurr
  i=i+1
}
n=length(dataset$id)
#---------------------------
#************************************
# FOR TESTING ONLY
#************************************
#nperson=500
#dataset=dataset[dataset$id<=nperson,]
#n=length(dataset$id)
#idvec=idvec[1:nperson]
#************************************
#---------------------------
#-----------------------------------------------------------------------------
# Define time known to have event
# Define time known to survive to
#-----------------------------------------------------------------------------
dataset_1=dataset[dataset$lastobs==1,]
nonsurvive_time=dataset_1$stop*(dataset_1$status==1)+cens_time*(dataset_1$status!=1)
nonsurvive_time_save=nonsurvive_time
survive_time=dataset_1$stop
#-----------------------------------------------------------------------------
cat('---------------------------------------------','\n')
cat('After sequential IDs: dataset=','\n')
print(summary(dataset))
cat('After sequential IDs: # of IDs in dataset =',length(unique(dataset$id)),'\n')
cat('After sequential IDs: # of IDs in dataset with lastobs=1 =',length(unique(dataset$id[dataset$lastobs==1])),'\n')
cat('---------------------------------------------','\n')

cat('Create sequential IDs in longitudinal dataset_nonmiss','\n')
cat('n_nonmiss=',n_nonmiss,'\n')
i=1
while (i <= n_nonmiss) {
#  cat('i=',i,'\n')
  temp=which(idvec==as.numeric(dataset_nonmiss$id[i]))
  if (length(temp) > 0) {
    dataset_nonmiss$id[i]=temp[1]
  } else {
    dataset_nonmiss$id[i]=NA
  }
  i=i+1
}
dataset_nonmiss=dataset_nonmiss[!is.na(dataset_nonmiss$id),]
n_nonmiss=length(dataset_nonmiss$id)
nonmiss_id=unique(dataset_nonmiss$id)

cat('---------------------------------------------','\n')
cat('After sequential IDs: longitudinal dataset=','\n')
print(summary(dataset_nonmiss))
cat('---------------------------------------------','\n')

#-----------------------------------------------------------------------------
#  Left Censor Dataset for Imputation Models to have sufficient sample sizes
#-----------------------------------------------------------------------------
# used_left_cens_time=imp_times_bp[1]
# dataset=dataset[dataset$stop >= used_left_cens_time-.000001,]
# dataset$start=dataset$start*(dataset$start >= used_left_cens_time)+used_left_cens_time*(dataset$start < used_left_cens_time)
# dataset_id=unique(dataset$id[dataset$lastobs==1])
# n=length(dataset$id)

#cat('Dataset after left censoring=','\n')
#print(dataset[1:50,])
# cat('After left-censoring: # of IDs in dataset=',length(unique(dataset$id)),'\n')
# cat('After left-censoring: # of IDs in dataset with lastobs=1=',length(dataset_id),'\n')
# cat('After left-censoring: # of events=',sum(dataset$status),'\n')
# cat('After left-censoring: survive_times of those in dataset=','\n')
print(summary(survive_time[dataset_id]))
cat('','\n')

#--------------------------------------------------------------------------
#  CALCULATE NUMBER OF OBS PRIOR TO S FOR EACH BIOMARKER FOR EACH SUBJECT  
#--------------------------------------------------------------------------
n_biomarker=rep(0,times=dim_bp*nperson)
dim(n_biomarker)=c(nperson,dim_bp)
i=1
while (i <=n_nonmiss) {
  if (dataset_nonmiss$obstime[i] <= s) {
    if (dataset_nonmiss$bptype[i]=='TYPE1') {
      n_biomarker[dataset_nonmiss$id[i],1]=n_biomarker[dataset_nonmiss$id[i],1]+1
    } else {
      if (dataset_nonmiss$bptype[i]=='TYPE2') {
        n_biomarker[dataset_nonmiss$id[i],2]=n_biomarker[dataset_nonmiss$id[i],2]+1
      } else {
        if (dataset_nonmiss$bptype[i]=='TYPE3') {
          n_biomarker[dataset_nonmiss$id[i],3]=n_biomarker[dataset_nonmiss$id[i],3]+1
        } else {
          if (dataset_nonmiss$bptype[i]=='TYPE4') {
            n_biomarker[dataset_nonmiss$id[i],4]=n_biomarker[dataset_nonmiss$id[i],4]+1
          } else {
            if (dataset_nonmiss$bptype[i]=='TYPE5') {
              n_biomarker[dataset_nonmiss$id[i],5]=n_biomarker[dataset_nonmiss$id[i],5]+1
            }
          }
        }
      }
    }
  }
  i=i+1
}






if (icox1_ind==1 | icox2_ind==1) {
#--------------------------------------------------------------------------
#
# IDENTIFY EVALUABLE SUBJECTS:
#    1) Survive to s with >=3 prior obs for each biomarker
#    2) Event before s+tau or survive to s+tau
#
#--------------------------------------------------------------------------
select_id=unique(dataset_nonmiss$id[survive_time[dataset_nonmiss$id]>s & dataset_nonmiss$bptype=='TYPE1' & dataset_nonmiss$obstime<= s & n_biomarker[dataset_nonmiss$id,1] >=4 ])
cat('Biomarker type =TYPE1','\n')
cat('# of Subjects known to survive to s and have >=4 prior obs biomarker =',length(select_id),'\n')
cat('--------------------------------------------','\n')
if (dim_bp >= 2) {
  select_id_2=unique(dataset_nonmiss$id[survive_time[dataset_nonmiss$id]>s & dataset_nonmiss$bptype=='TYPE2' & dataset_nonmiss$obstime<= s & n_biomarker[dataset_nonmiss$id,2] >=4 ])
  cat('Biomarker type =TYPE2','\n')
  cat('# of Subjects known to survive to s and have >=4 prior obs biomarker =',length(select_id_2),'\n')
  cat('--------------------------------------------','\n')
  select_id=select_id[select_id %in% select_id_2]
}
if (dim_bp >= 3) {
  select_id_3=unique(dataset_nonmiss$id[survive_time[dataset_nonmiss$id]>s & dataset_nonmiss$bptype=='TYPE3' & dataset_nonmiss$obstime<= s & n_biomarker[dataset_nonmiss$id,3] >=4 ])
  cat('Biomarker type =TYPE3','\n')
  cat('# of Subjects known to survive to s and have >=4 prior obs biomarker =',length(select_id_3),'\n')
  cat('--------------------------------------------','\n')
  select_id=select_id[select_id %in% select_id_3]
}
if (dim_bp >= 4) {
  select_id_4=unique(dataset_nonmiss$id[survive_time[dataset_nonmiss$id]>s & dataset_nonmiss$bptype=='TYPE4' & dataset_nonmiss$obstime<= s & n_biomarker[dataset_nonmiss$id,4] >=4 ])
  cat('Biomarker type =TYPE4','\n')
  cat('# of Subjects known to survive to s and have >=4 prior obs biomarker =',length(select_id_4),'\n')
  cat('--------------------------------------------','\n')
  select_id=select_id[select_id %in% select_id_4]
}
if (dim_bp >= 5) {
  select_id_5=unique(dataset_nonmiss$id[survive_time[dataset_nonmiss$id]>s & dataset_nonmiss$bptype=='TYPE5' & dataset_nonmiss$obstime<= s & n_biomarker[dataset_nonmiss$id,5] >=4 ])
  cat('Biomarker type =TYPE5','\n')
  cat('# of Subjects known to survive to s and have >=4 prior obs biomarker =',length(select_id_5),'\n')
  cat('--------------------------------------------','\n')
  select_id=select_id[select_id %in% select_id_5]
}
cat('--------------------------------------------','\n')
cat('# of Subjects known to survive to s and have >=4 prior obs for each biomarker =',length(select_id),'\n')
cat('--------------------------------------------','\n')
cat('--------------------------------------------','\n')
#cat('# of IDs surviving to s with valid EB estimates (EVAL POP)=',length(select_id),'\n')
dataset_eval=dataset[dataset$id %in% select_id,]
dataset_eval_1=dataset[dataset$id %in% select_id & dataset$lastobs==1,]
survive_time_1=survive_time[select_id]
#cat('# of IDs in select_id[[1]]=',length(select_id[[1]]),'\n')
#cat('# of IDs in dataset with ID in select_id[[1]]=',length(unique(dataset$id[dataset$id %in% select_id[[1]]])),'\n')
#cat('# of events in dataset_eval=',sum(dataset_eval$status),'\n')
cat('# of obs in dataset_eval_1=',length(dataset_eval_1$id),'\n')
#cat('Dataset_eval_1 =','\n')
#print(dataset_eval_1[1:50,])
cat('--------------------------------------------','\n')
#cat('--------------------------------------------','\n')
#cat('summary of survive_time for those in EVAL POP=','\n')
#print(summary(survive_time_1))
#cat('dataset_eval=','\n')
#print(dataset_eval[1:50,])
#cat('survive_time_1=','\n')
#print(survive_time_1[1:50])

#eval_id=select_id[[1]][survive_time_1 > s+tau | (survive_time_1 < s+tau & dataset_model_lastobs$status==1)]
#n_eval=length(eval_id)
#dataset_eval=dataset[dataset$id %in% eval_id,]
#dataset_eval_lastobs=dataset[dataset$id %in% eval_id & dataset$lastobs==1,]

eval_id=select_id
n_eval=length(eval_id)

survive_id=select_id[survive_time_1 > s+tau]
event_id=select_id[survive_time_1 < s+tau & dataset_eval_1$status==1]

cat('--------------------------------------------','\n')
cat('--------------------------------------------','\n')
#cat('# of events in dataset_eval=',sum(dataset_eval$status),'\n')
#cat('--------------------------------------------','\n')
#cat('summary of survive_time for those in eval_id=','\n')
#print(summary(survive_time[eval_id]))
cat('********************************************','\n')
cat('--------------------------------------------','\n')
cat('# of IDs surviving to s with valid EB estimates (EVAL POP)=',length(select_id),'\n')
cat('  # of Evaluable IDs surviving to s+tau =',length(survive_id),'\n')
cat('  # of Evaluable IDs with event prior to s+tau =',length(event_id),'\n')
cat('nonsurvive_time[EVAL POP] =','\n')
print(nonsurvive_time[eval_id][1:40])
cat('--------------------------------------------','\n')
cat('********************************************','\n')
cat('--------------------------------------------','\n')

dataset_impute=dataset_eval
n=length(dataset_impute$id)





#---------------------------------------------------------------------
#---------------------------------------------------------------------
KATIE: 
  1) RESTRICT EVAL_ID TO THOSE WITH OBSERVED BIOMARKER VALUES FOR ALL THREE BIOMARKERS AT TIME S+TAU
  2) SKIP TO LINE 1585     LOOP OVER EVALUABLE IDs
#---------------------------------------------------------------------
#---------------------------------------------------------------------





#cat('dataset_nonmiss=','\n')
#print(dataset_nonmiss[1:20,])
#cat('dataset_impute=','\n')
#print(dataset_impute[1:20,])

#cat('dataset before NAs added to Dataset_nonmiss=','\n')
#print(dataset[1:80,])
#cat('dataset_impute before NAs added to Dataset_nonmiss=','\n')
#print(dataset_impute[1:80,])
#--------------------------------------------------------------------------
#
# For Each Obs in Dataset_impute Add One Obs With BP_CURR=NA to Dataset_nonmiss
#
#--------------------------------------------------------------------------
#i=1
#idprev=0
#i_nonmiss=1
#temp_start=0
#while (i <= n) {
#  cat('i=',i,'\n')
#  id=dataset_impute$id[i]
#  cat('id=',id,'\n')
#  if (id != idprev) {
#    while (dataset_nonmiss$id[i_nonmiss] < id) {
#      i_nonmiss=i_nonmiss+1
#    }
#  }
#  cat('i_nonmiss=',i_nonmiss,'\n')
#  cat('temp_start=',temp_start,'\n')
#  if (temp_start==0) {
#    temp=dataset_nonmiss[i_nonmiss,]
#    temp$bp_curr[1]=NA
#    temp$obstime[1]=dataset_impute$start[i]
#    temp_start=1
#  } else { 
#    temp1=dataset_nonmiss[i_nonmiss,]
#    temp1$bp_curr[1]=NA
#    temp1$obstime[1]=dataset_impute$start[i]
#    temp=rbind(temp,temp1)
#  }
#  idprev=id
#  i=i+1
#}
#temp$bptype[is.na(temp$bp_curr)]='TYPE1'
#dataset_nonmiss_imp=rbind(dataset_nonmiss,temp)
#dataset_nonmiss_imp=dataset_nonmiss_imp[order(dataset_nonmiss_imp$id,dataset_nonmiss_imp$obstime),]
#n_nonmiss_imp=length(dataset_nonmiss_imp$id)
#--------------------------------------------------------------------------
#
# For Each Obs in Dataset_impute Add One Obs With BP_CURR=NA to Dataset_nonmiss
#
#--------------------------------------------------------------------------
id=dataset_impute$id
start=dataset_impute$start
lastobs=dataset_impute$lastobs
temp=data.frame(id=id,bp_curr=rep(NA,times=n),obstime=start,bptype=rep('TYPE1',times=n),dataset_impute[,11:(10+ncov)],lastobs=lastobs)

#cat('temp=','\n')
#print(temp[1:20,])

dataset_nonmiss_imp=rbind(dataset_nonmiss,temp)
dataset_nonmiss_imp=dataset_nonmiss_imp[order(dataset_nonmiss_imp$id,dataset_nonmiss_imp$obstime),]
n_nonmiss_imp=length(dataset_nonmiss_imp$id)

dataset_nonmiss_imp_save=dataset_nonmiss_imp

#dataset_nonmiss_save=list()

#cat('dataset_nonmiss_imp=','\n')
#print(dataset_nonmiss_imp[1:40,])

dataset_impute$int1_s=rep(NA,times=n)
dataset_impute$int2_s=rep(NA,times=n)
dataset_impute$int3_s=rep(NA,times=n)
dataset_impute$int4_s=rep(NA,times=n)
dataset_impute$int5_s=rep(NA,times=n)
dataset_impute$slope1_s=rep(NA,times=n)
dataset_impute$slope2_s=rep(NA,times=n)
dataset_impute$slope3_s=rep(NA,times=n)
dataset_impute$slope4_s=rep(NA,times=n)
dataset_impute$slope5_s=rep(NA,times=n)
#--------------------------------------------------------------------------
#
# Fit Longitudinal Models for Biomarker Process 
#
#--------------------------------------------------------------------------
fit1=list()
fit2=list()
fit3=list()
fit4=list()
fit5=list()
id_model1=list()
id_model2=list()
id_model3=list()
id_model4=list()
id_model5=list()

j=1
while (j <= num_imp_bp) {
  dataset_nonmiss_imp=dataset_nonmiss_imp_save
  nspline[1]=max(1,min(floor(imp_times_bp[j]/spline_width[1]),nspline_max[1]))
  if (dim_bp >=2) {nspline[2]=max(1,min(floor(imp_times_bp[j]/spline_width[2]),nspline_max[2]))}
  if (dim_bp >=3) {nspline[3]=max(1,min(floor(imp_times_bp[j]/spline_width[3]),nspline_max[3]))}
  if (dim_bp >=4) {nspline[4]=max(1,min(floor(imp_times_bp[j]/spline_width[4]),nspline_max[4]))}
  if (dim_bp >=5) {nspline[5]=max(1,min(floor(imp_times_bp[j]/spline_width[5]),nspline_max[5]))}
  cat('--------------------------------------------','\n')
  cat('--------------------------------------------','\n')
  cat('Imp model for imp time =',imp_times_bp[j],'\n')
  cat('--------------------------------------------','\n')
  cat('nspline[1] for Biomarker #1 =',nspline[1],'\n')
  cat('nspline[2] for Biomarker #2 =',nspline[2],'\n')
  cat('nspline[3] for Biomarker #3 =',nspline[3],'\n')

#--------------------------------------------------------------------------
# Create Ordered Formulas for Longitudinal Models of Biomarker Process 
#--------------------------------------------------------------------------
  jj=1
  base_names=""
  while (jj <= ncov) {
    if (base_cov_factor[jj]==0) {
      base_names=paste(base_names,"+base_cov.",jj,sep="")
    } else {
      base_names=paste(base_names,"+as.factor(base_cov.",jj,")",sep="")
    }
    jj=jj+1
  }
  names=base_names
  if (nspline[1]>=1) {
    names=paste(names,"+spline111+spline112+spline113",sep="")
  }
  if (nspline[1]>=2) {
    names=paste(names,"+spline121+spline122+spline123",sep="")
  }
  if (nspline[1]>=3) {
    names=paste(names,"+spline131+spline132+spline133",sep="")
  }
  if (nspline[1]>=4) {
    names=paste(names,"+spline141+spline142+spline143",sep="")
  }
  if (nspline[1]>=5) {
    names=paste(names,"+spline151+spline152+spline153",sep="")
  }
  if (nspline[1]>=6) {
    names=paste(names,"+spline161+spline162+spline163",sep="")
  }
  if (nspline[1]>=7) {
    names=paste(names,"+spline171+spline172+spline173",sep="")
  }
  if (nspline[1]>=8) {
    names=paste(names,"+spline181+spline182+spline183",sep="")
  }

  names=paste(names,"+nonsurvive+nonsurvive1+ (obstime | id)",sep="")
#  cat('names =',names,'\n')
  (formula1=as.formula(paste("bp_curr~1+",names,sep="")))

  if (dim_bp >=2) {
    names=base_names
    if (nspline[2]>=1) {
      names=paste(names,"+spline211+spline212+spline213",sep="")
    }
    if (nspline[2]>=2) {
      names=paste(names,"+spline221+spline222+spline223",sep="")
    }
    if (nspline[2]>=3) {
      names=paste(names,"+spline231+spline232+spline233",sep="")
    }
    if (nspline[2]>=4) {
      names=paste(names,"+spline241+spline242+spline243",sep="")
    }
    if (nspline[2]>=5) {
      names=paste(names,"+spline251+spline252+spline253",sep="")
    }
    if (nspline[2]>=6) {
      names=paste(names,"+spline261+spline262+spline263",sep="")
    }
    if (nspline[2]>=7) {
      names=paste(names,"+spline271+spline272+spline273",sep="")
    }
    if (nspline[2]>=8) {
      names=paste(names,"+spline281+spline282+spline283",sep="")
    }

    names=paste(names,"+nonsurvive+nonsurvive1+ (obstime | id)",sep="")
#    cat('names =',names,'\n')
    (formula2=as.formula(paste("bp_curr~1+",names,sep="")))
  }

  if (dim_bp >=3) {
    names=base_names
    if (nspline[3]>=1) {
      names=paste(names,"+spline311+spline312+spline313",sep="")
    }
    if (nspline[3]>=2) {
      names=paste(names,"+spline321+spline322+spline323",sep="")
    }
    if (nspline[3]>=3) {
      names=paste(names,"+spline331+spline332+spline333",sep="")
    }
    if (nspline[3]>=4) {
      names=paste(names,"+spline341+spline342+spline343",sep="")
    }
    if (nspline[3]>=5) {
      names=paste(names,"+spline351+spline352+spline353",sep="")
    }
    if (nspline[3]>=6) {
      names=paste(names,"+spline361+spline362+spline363",sep="")
    }
    if (nspline[3]>=7) {
      names=paste(names,"+spline371+spline372+spline373",sep="")
    }
    if (nspline[3]>=8) {
      names=paste(names,"+spline381+spline382+spline383",sep="")
    }

    names=paste(names,"+nonsurvive+nonsurvive1+ (obstime | id)",sep="")
#    cat('names =',names,'\n')
    (formula3=as.formula(paste("bp_curr~1+",names,sep="")))
  }

  if (dim_bp >=4) {
    names=base_names
    if (nspline[4]>=1) {
      names=paste(names,"+spline411+spline412+spline413",sep="")
    }
    if (nspline[4]>=2) {
      names=paste(names,"+spline421+spline422+spline423",sep="")
    }
    if (nspline[4]>=3) {
      names=paste(names,"+spline431+spline432+spline433",sep="")
    } 
    if (nspline[4]>=4) {
      names=paste(names,"+spline441+spline442+spline443",sep="")
    }
    if (nspline[4]>=5) {
      names=paste(names,"+spline451+spline452+spline453",sep="")
    }
    if (nspline[4]>=6) {
      names=paste(names,"+spline461+spline462+spline463",sep="")
    }
    if (nspline[4]>=7) {
      names=paste(names,"+spline471+spline472+spline473",sep="")
    }
    if (nspline[4]>=8) {
      names=paste(names,"+spline481+spline482+spline483",sep="")
    }

    names=paste(names,"+nonsurvive+nonsurvive1+ (obstime | id)",sep="")
#    cat('names =',names,'\n')
    (formula4=as.formula(paste("bp_curr~1+",names,sep="")))
  }

  if (dim_bp >=5) {
    names=base_names
    if (nspline[5]>=1) {
      names=paste(names,"+spline511+spline512+spline513",sep="")
    }
    if (nspline[5]>=2) {
      names=paste(names,"+spline521+spline522+spline523",sep="")
    }
    if (nspline[5]>=3) {
      names=paste(names,"+spline531+spline532+spline533",sep="")
    } 
    if (nspline[5]>=4) {
      names=paste(names,"+spline541+spline542+spline543",sep="")
    }
    if (nspline[5]>=5) {
      names=paste(names,"+spline551+spline552+spline553",sep="")
    }
    if (nspline[5]>=6) {
      names=paste(names,"+spline561+spline562+spline563",sep="")
    }
    if (nspline[5]>=7) {
      names=paste(names,"+spline571+spline572+spline573",sep="")
    }
    if (nspline[5]>=8) {
      names=paste(names,"+spline581+spline582+spline583",sep="")
    }

    names=paste(names,"+nonsurvive+nonsurvive1+ (obstime | id)",sep="")
#    cat('names =',names,'\n')
    (formula5=as.formula(paste("bp_curr~1+",names,sep="")))
  }

#--------------------------------------------------------------------------
# Spline Creation 
#--------------------------------------------------------------------------
  knot1=spline_width[1]
  knot2=2*spline_width[1]
  knot3=3*spline_width[1]
  knot4=4*spline_width[1]
  knot5=5*spline_width[1]
  knot6=6*spline_width[1]
  knot7=7*spline_width[1]
  dataset_nonmiss_imp$spline111=dataset_nonmiss_imp$obstime
  dataset_nonmiss_imp$spline112=dataset_nonmiss_imp$obstime**2
  dataset_nonmiss_imp$spline113=dataset_nonmiss_imp$obstime**3
  if (nspline[1] >= 2) {
    dataset_nonmiss_imp$spline111=dataset_nonmiss_imp$obstime**1*(dataset_nonmiss_imp$obstime < knot1) +knot1**1*(dataset_nonmiss_imp$obstime >= knot1)
    dataset_nonmiss_imp$spline112=dataset_nonmiss_imp$obstime**2*(dataset_nonmiss_imp$obstime < knot1) +knot1**2*(dataset_nonmiss_imp$obstime >= knot1)
    dataset_nonmiss_imp$spline113=dataset_nonmiss_imp$obstime**3*(dataset_nonmiss_imp$obstime < knot1) +knot1**3*(dataset_nonmiss_imp$obstime >= knot1)
    dataset_nonmiss_imp$spline121=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime >= knot1) 
    dataset_nonmiss_imp$spline122=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime >= knot1) 
    dataset_nonmiss_imp$spline123=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime >= knot1) 
  }
  if (nspline[1] >= 3) {
    dataset_nonmiss_imp$spline121=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**1*(dataset_nonmiss_imp$obstime >= knot2)
    dataset_nonmiss_imp$spline122=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**2*(dataset_nonmiss_imp$obstime >= knot2)
    dataset_nonmiss_imp$spline123=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**3*(dataset_nonmiss_imp$obstime >= knot2)
    dataset_nonmiss_imp$spline131=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime >= knot2) 
    dataset_nonmiss_imp$spline132=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime >= knot2) 
    dataset_nonmiss_imp$spline133=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime >= knot2) 
  }
  if (nspline[1] >= 4) {
    dataset_nonmiss_imp$spline131=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**1*(dataset_nonmiss_imp$obstime >= knot3)
    dataset_nonmiss_imp$spline132=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**2*(dataset_nonmiss_imp$obstime >= knot3)
    dataset_nonmiss_imp$spline133=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**3*(dataset_nonmiss_imp$obstime >= knot3)
    dataset_nonmiss_imp$spline141=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime >= knot3) 
    dataset_nonmiss_imp$spline142=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime >= knot3) 
    dataset_nonmiss_imp$spline143=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime >= knot3) 
  }
  if (nspline[1] >= 5) {
    dataset_nonmiss_imp$spline141=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**1*(dataset_nonmiss_imp$obstime >= knot4)
    dataset_nonmiss_imp$spline142=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**2*(dataset_nonmiss_imp$obstime >= knot4)
    dataset_nonmiss_imp$spline143=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**3*(dataset_nonmiss_imp$obstime >= knot4)
    dataset_nonmiss_imp$spline151=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime >= knot4) 
    dataset_nonmiss_imp$spline152=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime >= knot4) 
    dataset_nonmiss_imp$spline153=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime >= knot4) 
  }
  if (nspline[1] >= 6) {
    dataset_nonmiss_imp$spline151=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**1*(dataset_nonmiss_imp$obstime >= knot5)
    dataset_nonmiss_imp$spline152=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**2*(dataset_nonmiss_imp$obstime >= knot5)
    dataset_nonmiss_imp$spline153=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**3*(dataset_nonmiss_imp$obstime >= knot5)
    dataset_nonmiss_imp$spline161=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime >= knot5) 
    dataset_nonmiss_imp$spline162=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime >= knot5) 
    dataset_nonmiss_imp$spline163=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime >= knot5) 
  }
  if (nspline[1] >= 7) {
    dataset_nonmiss_imp$spline161=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**1*(dataset_nonmiss_imp$obstime >= knot6)
    dataset_nonmiss_imp$spline162=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**2*(dataset_nonmiss_imp$obstime >= knot6)
    dataset_nonmiss_imp$spline163=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**3*(dataset_nonmiss_imp$obstime >= knot6)
    dataset_nonmiss_imp$spline171=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime >= knot6) 
    dataset_nonmiss_imp$spline172=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime >= knot6) 
    dataset_nonmiss_imp$spline173=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime >= knot6) 
  }
  if (nspline[1] >= 8) {
    dataset_nonmiss_imp$spline171=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**1*(dataset_nonmiss_imp$obstime >= knot7)
    dataset_nonmiss_imp$spline172=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**2*(dataset_nonmiss_imp$obstime >= knot7)
    dataset_nonmiss_imp$spline173=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**3*(dataset_nonmiss_imp$obstime >= knot7)
    dataset_nonmiss_imp$spline181=(dataset_nonmiss_imp$obstime-knot7)**1*(dataset_nonmiss_imp$obstime >= knot7) 
    dataset_nonmiss_imp$spline182=(dataset_nonmiss_imp$obstime-knot7)**2*(dataset_nonmiss_imp$obstime >= knot7) 
    dataset_nonmiss_imp$spline183=(dataset_nonmiss_imp$obstime-knot7)**3*(dataset_nonmiss_imp$obstime >= knot7) 
  }

#--------------------------------------------------------------------------
# END: Spline Creation 
#--------------------------------------------------------------------------
  dataset_nonmiss_imp$nonsurvive=(nonsurvive_time[dataset_nonmiss_imp$id] <= imp_times_bp[j])
  dataset_nonmiss_imp$nonsurvive1=dataset_nonmiss_imp$nonsurvive*dataset_nonmiss_imp$obstime
#  dataset_nonmiss_save[[j]]=dataset_nonmiss_imp

  select=which(dataset_nonmiss_imp$obstime <= imp_times_bp[j])
  dataset_nonmiss_select=dataset_nonmiss_imp[select,]

#  cat('dataset_nonmiss_select =','\n')
#  print(dataset_nonmiss_select[1:80,])
  cat('# unique ID in dataset_nonmiss_select =',length(unique(dataset_nonmiss_select$id)),'\n')
#  cat('dataset_nonmiss_imp =','\n')
#  print(dataset_nonmiss_imp[1:80,])
  cat('# unique ID in dataset_nonmiss_imp =',length(unique(dataset_nonmiss_imp$id)),'\n')

#-------------------------------------------------
# FIT LINEAR FIXED EFFECTS
  fit1[[j]]=try(lmer(formula1,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=='TYPE1'))
  id_model1[[j]]=unique(dataset_nonmiss_select$id[dataset_nonmiss_select$bptype=='TYPE1' & !is.na(dataset_nonmiss_select$bp_curr)])
  cat('--------------------------------------------','\n')
  cat('BPTYPE=1','\n')
  print(summary(fit1[[j]]))
  cat('--------------------------------------------','\n')
  dataset_impute$bp_curr1[dataset_impute$start==imp_times_bp[j]]=predict(fit1[[j]],newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr) & dataset_nonmiss_imp$bptype=='TYPE1',])[dataset_impute$start==imp_times_bp[j]]
  if (imp_times_bp[j]==s) {
    dataset_impute$int1_s[dataset_impute$start==imp_times_bp[j]]=random.effects(fit1[[j]])[[1]][which(id_model1[[j]] %in% dataset_impute$id[dataset_impute$start==imp_times_bp[j]]),1]
    dataset_impute$slope1_s[dataset_impute$start==imp_times_bp[j]]=random.effects(fit1[[j]])[[1]][which(id_model1[[j]] %in% dataset_impute$id[dataset_impute$start==imp_times_bp[j]]),2]
  }

#  cat('Biomarker #1: dataset_nonmiss_select for ID=5','\n')
#  print(dataset_nonmiss_select[dataset_nonmiss_select$id==5,])
#  cat('----------------------------------------','\n')

  if (dim_bp>=2) {

#  SPLINE CREATION
    knot1=spline_width[2]
    knot2=2*spline_width[2]
    knot3=3*spline_width[2]
    knot4=4*spline_width[2]
    knot5=5*spline_width[2]
    knot6=6*spline_width[2]
    knot7=7*spline_width[2]
    dataset_nonmiss_imp$spline211=dataset_nonmiss_imp$obstime
    dataset_nonmiss_imp$spline212=dataset_nonmiss_imp$obstime**2
    dataset_nonmiss_imp$spline213=dataset_nonmiss_imp$obstime**3
    if (nspline[2] >= 2) {
      dataset_nonmiss_imp$spline211=dataset_nonmiss_imp$obstime**1*(dataset_nonmiss_imp$obstime < knot1) +knot1**1*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline212=dataset_nonmiss_imp$obstime**2*(dataset_nonmiss_imp$obstime < knot1) +knot1**2*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline213=dataset_nonmiss_imp$obstime**3*(dataset_nonmiss_imp$obstime < knot1) +knot1**3*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline221=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline222=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline223=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime >= knot1) 
    }
    if (nspline[2] >= 3) {
      dataset_nonmiss_imp$spline221=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**1*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline222=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**2*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline223=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**3*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline231=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline232=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline233=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime >= knot2) 
    }
    if (nspline[2] >= 4) {
      dataset_nonmiss_imp$spline231=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**1*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline232=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**2*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline233=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**3*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline241=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline242=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline243=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime >= knot3) 
    }
    if (nspline[2] >= 5) {
      dataset_nonmiss_imp$spline241=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**1*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline242=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**2*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline243=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**3*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline251=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline252=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline253=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime >= knot4) 
    }
    if (nspline[2] >= 6) {
      dataset_nonmiss_imp$spline251=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**1*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline252=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**2*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline253=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**3*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline261=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline262=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline263=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime >= knot5) 
    }
    if (nspline[2] >= 7) {
      dataset_nonmiss_imp$spline261=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**1*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline262=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**2*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline263=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**3*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline271=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline272=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline273=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime >= knot6) 
    }
    if (nspline[2] >= 8) {
      dataset_nonmiss_imp$spline271=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**1*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline272=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**2*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline273=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**3*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline281=(dataset_nonmiss_imp$obstime-knot7)**1*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline282=(dataset_nonmiss_imp$obstime-knot7)**2*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline283=(dataset_nonmiss_imp$obstime-knot7)**3*(dataset_nonmiss_imp$obstime >= knot7) 
    }
#    dataset_nonmiss_save[[j]]=dataset_nonmiss_imp
    select=which(dataset_nonmiss_imp$obstime <= imp_times_bp[j])
    dataset_nonmiss_select=dataset_nonmiss_imp[select,]

    fit2[[j]]=try(lmer(formula2,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=='TYPE2'))
    id_model2[[j]]=unique(dataset_nonmiss_select$id[dataset_nonmiss_select$bptype=='TYPE2' & !is.na(dataset_nonmiss_select$bp_curr)])
    cat('--------------------------------------------','\n')
    cat('BPTYPE=2','\n')
    print(summary(fit2[[j]]))
    cat('--------------------------------------------','\n')
    dataset_impute$bp_curr2[dataset_impute$start==imp_times_bp[j]]=predict(fit2[[j]],newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr) & dataset_nonmiss_imp$bptype=='TYPE1',])[dataset_impute$start==imp_times_bp[j]]
    if (imp_times_bp[j]==s) {
      dataset_impute$int2_s[dataset_impute$start==imp_times_bp[j]]=random.effects(fit2[[j]])[[1]][which(id_model2[[j]] %in% dataset_impute$id[dataset_impute$start==imp_times_bp[j]]),1]
      dataset_impute$slope2_s[dataset_impute$start==imp_times_bp[j]]=random.effects(fit2[[j]])[[1]][which(id_model2[[j]] %in% dataset_impute$id[dataset_impute$start==imp_times_bp[j]]),2]
#      cat('dataset_impute for ID=10229 after model at s=','\n')
#      print(dataset_impute[dataset_impute$id==10229,])
#      cat('dataset_nonmiss_imp for ID=10229=','\n')
#      print(dataset_nonmiss_imp[dataset_nonmiss_imp$id==10229,])
#      cat('dataset_nonmiss_select for ID=10229=','\n')
#      print(dataset_nonmiss_select[dataset_nonmiss_select$id==10229,])
    }
  }

  if (dim_bp>=3) {

#  SPLINE CREATION
    knot1=spline_width[3]
    knot2=2*spline_width[3]
    knot3=3*spline_width[3]
    knot4=4*spline_width[3]
    knot5=5*spline_width[3]
    knot6=6*spline_width[3]
    knot7=7*spline_width[3]
    dataset_nonmiss_imp$spline311=dataset_nonmiss_imp$obstime
    dataset_nonmiss_imp$spline312=dataset_nonmiss_imp$obstime**2
    dataset_nonmiss_imp$spline313=dataset_nonmiss_imp$obstime**3
    if (nspline[3] >= 2) {
      dataset_nonmiss_imp$spline311=dataset_nonmiss_imp$obstime**1*(dataset_nonmiss_imp$obstime < knot1) +knot1**1*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline312=dataset_nonmiss_imp$obstime**2*(dataset_nonmiss_imp$obstime < knot1) +knot1**2*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline313=dataset_nonmiss_imp$obstime**3*(dataset_nonmiss_imp$obstime < knot1) +knot1**3*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline321=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline322=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline323=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime >= knot1) 
    }
    if (nspline[3] >= 3) {
      dataset_nonmiss_imp$spline321=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**1*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline322=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**2*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline323=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**3*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline331=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline332=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline333=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime >= knot2) 
    }
    if (nspline[3] >= 4) {
      dataset_nonmiss_imp$spline331=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**1*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline332=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**2*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline333=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**3*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline341=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline342=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline343=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime >= knot3) 
    }
    if (nspline[3] >= 5) {
      dataset_nonmiss_imp$spline341=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**1*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline342=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**2*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline343=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**3*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline351=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline352=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline353=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime >= knot4) 
    }
    if (nspline[3] >= 6) {
      dataset_nonmiss_imp$spline351=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**1*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline352=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**2*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline353=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**3*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline361=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline362=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline363=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime >= knot5) 
    }
    if (nspline[3] >= 7) {
      dataset_nonmiss_imp$spline361=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**1*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline362=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**2*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline363=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**3*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline371=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline372=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline373=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime >= knot6) 
    }
    if (nspline[3] >= 8) {
      dataset_nonmiss_imp$spline371=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**1*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline372=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**2*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline373=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**3*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline381=(dataset_nonmiss_imp$obstime-knot7)**1*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline382=(dataset_nonmiss_imp$obstime-knot7)**2*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline383=(dataset_nonmiss_imp$obstime-knot7)**3*(dataset_nonmiss_imp$obstime >= knot7) 
    }
#    dataset_nonmiss_save[[j]]=dataset_nonmiss_imp
    select=which(dataset_nonmiss_imp$obstime <= imp_times_bp[j])
    dataset_nonmiss_select=dataset_nonmiss_imp[select,]

    fit3[[j]]=try(lmer(formula3,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=='TYPE3'))
    id_model3[[j]]=unique(dataset_nonmiss_select$id[dataset_nonmiss_select$bptype=='TYPE3' & !is.na(dataset_nonmiss_select$bp_curr)])
    cat('--------------------------------------------','\n')
    cat('BPTYPE=3','\n')
    print(summary(fit3[[j]]))
    cat('--------------------------------------------','\n')
    dataset_impute$bp_curr3[dataset_impute$start==imp_times_bp[j]]=predict(fit3[[j]],newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr) & dataset_nonmiss_imp$bptype=='TYPE1',])[dataset_impute$start==imp_times_bp[j]]
    if (imp_times_bp[j]==s) {
      dataset_impute$int3_s[dataset_impute$start==imp_times_bp[j]]=random.effects(fit3[[j]])[[1]][which(id_model3[[j]] %in% dataset_impute$id[dataset_impute$start==imp_times_bp[j]]),1]
      dataset_impute$slope3_s[dataset_impute$start==imp_times_bp[j]]=random.effects(fit3[[j]])[[1]][which(id_model3[[j]] %in% dataset_impute$id[dataset_impute$start==imp_times_bp[j]]),2]
    }
  }

  if (dim_bp>=4) {

#  SPLINE CREATION
    knot1=spline_width[4]
    knot2=2*spline_width[4]
    knot3=3*spline_width[4]
    knot4=4*spline_width[4]
    knot5=5*spline_width[4]
    knot6=6*spline_width[4]
    knot7=7*spline_width[4]
    dataset_nonmiss_imp$spline411=dataset_nonmiss_imp$obstime
    dataset_nonmiss_imp$spline412=dataset_nonmiss_imp$obstime**2
    dataset_nonmiss_imp$spline413=dataset_nonmiss_imp$obstime**3
    if (nspline[4] >= 2) {
      dataset_nonmiss_imp$spline411=dataset_nonmiss_imp$obstime**1*(dataset_nonmiss_imp$obstime < knot1) +knot1**1*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline412=dataset_nonmiss_imp$obstime**2*(dataset_nonmiss_imp$obstime < knot1) +knot1**2*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline413=dataset_nonmiss_imp$obstime**3*(dataset_nonmiss_imp$obstime < knot1) +knot1**3*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline421=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline422=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline423=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime >= knot1) 
    }
    if (nspline[4] >= 3) {
      dataset_nonmiss_imp$spline421=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**1*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline422=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**2*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline423=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**3*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline431=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline432=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline433=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime >= knot2) 
    }
    if (nspline[4] >= 4) {
      dataset_nonmiss_imp$spline431=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**1*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline432=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**2*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline433=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**3*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline441=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline442=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline443=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime >= knot3) 
    }
    if (nspline[4] >= 5) {
      dataset_nonmiss_imp$spline441=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**1*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline442=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**2*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline443=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**3*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline451=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline452=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline453=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime >= knot4) 
    }
    if (nspline[4] >= 6) {
      dataset_nonmiss_imp$spline451=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**1*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline452=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**2*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline453=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**3*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline461=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline462=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline463=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime >= knot5) 
    }
    if (nspline[4] >= 7) {
      dataset_nonmiss_imp$spline461=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**1*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline462=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**2*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline463=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**3*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline471=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline472=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline473=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime >= knot6) 
    }
    if (nspline[4] >= 8) {
      dataset_nonmiss_imp$spline471=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**1*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline472=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**2*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline473=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**3*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline481=(dataset_nonmiss_imp$obstime-knot7)**1*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline482=(dataset_nonmiss_imp$obstime-knot7)**2*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline483=(dataset_nonmiss_imp$obstime-knot7)**3*(dataset_nonmiss_imp$obstime >= knot7) 
    }
#    dataset_nonmiss_save[[j]]=dataset_nonmiss_imp
    select=which(dataset_nonmiss_imp$obstime <= imp_times_bp[j])
    dataset_nonmiss_select=dataset_nonmiss_imp[select,]

    fit4[[j]]=try(lmer(formula4,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=='TYPE4'))
    id_model4[[j]]=unique(dataset_nonmiss_select$id[dataset_nonmiss_select$bptype=='TYPE4' & !is.na(dataset_nonmiss_select$bp_curr)])
    cat('--------------------------------------------','\n')
    cat('BPTYPE=4','\n')
    print(summary(fit4[[j]]))
    cat('--------------------------------------------','\n')
    dataset_impute$bp_curr4[dataset_impute$start==imp_times_bp[j]]=predict(fit4[[j]],newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr) & dataset_nonmiss_imp$bptype=='TYPE1',])[dataset_impute$start==imp_times_bp[j]]
    if (imp_times_bp[j]==s) {
      dataset_impute$int4_s[dataset_impute$start==imp_times_bp[j]]=random.effects(fit4[[j]])[[1]][which(id_model4[[j]] %in% dataset_impute$id[dataset_impute$start==imp_times_bp[j]]),1]
      dataset_impute$slope4_s[dataset_impute$start==imp_times_bp[j]]=random.effects(fit4[[j]])[[1]][which(id_model4[[j]] %in% dataset_impute$id[dataset_impute$start==imp_times_bp[j]]),2]
    }
  }

  if (dim_bp>=5) {

#  SPLINE CREATION
    knot1=spline_width[5]
    knot2=2*spline_width[5]
    knot3=3*spline_width[5]
    knot4=4*spline_width[5]
    knot5=5*spline_width[5]
    knot6=6*spline_width[5]
    knot7=7*spline_width[5]
    dataset_nonmiss_imp$spline511=dataset_nonmiss_imp$obstime
    dataset_nonmiss_imp$spline512=dataset_nonmiss_imp$obstime**2
    dataset_nonmiss_imp$spline513=dataset_nonmiss_imp$obstime**3
    if (nspline[5] >= 2) {
      dataset_nonmiss_imp$spline511=dataset_nonmiss_imp$obstime**1*(dataset_nonmiss_imp$obstime < knot1) +knot1**1*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline512=dataset_nonmiss_imp$obstime**2*(dataset_nonmiss_imp$obstime < knot1) +knot1**2*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline513=dataset_nonmiss_imp$obstime**3*(dataset_nonmiss_imp$obstime < knot1) +knot1**3*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline521=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline522=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline523=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime >= knot1) 
    }
    if (nspline[5] >= 3) {
      dataset_nonmiss_imp$spline521=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**1*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline522=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**2*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline523=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**3*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline531=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline532=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline533=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime >= knot2) 
    }
    if (nspline[5] >= 4) {
      dataset_nonmiss_imp$spline531=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**1*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline532=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**2*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline533=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**3*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline541=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline542=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline543=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime >= knot3) 
    }
    if (nspline[5] >= 5) {
      dataset_nonmiss_imp$spline541=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**1*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline542=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**2*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline543=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**3*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline551=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline552=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline553=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime >= knot4) 
    }
    if (nspline[5] >= 6) {
      dataset_nonmiss_imp$spline551=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**1*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline552=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**2*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline553=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**3*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline561=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline562=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline563=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime >= knot5) 
    }
    if (nspline[5] >= 7) {
      dataset_nonmiss_imp$spline561=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**1*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline562=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**2*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline563=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**3*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline571=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline572=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline573=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime >= knot6) 
    }
    if (nspline[5] >= 8) {
      dataset_nonmiss_imp$spline571=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**1*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline572=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**2*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline573=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**3*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline581=(dataset_nonmiss_imp$obstime-knot7)**1*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline582=(dataset_nonmiss_imp$obstime-knot7)**2*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline583=(dataset_nonmiss_imp$obstime-knot7)**3*(dataset_nonmiss_imp$obstime >= knot7) 
    }
#    dataset_nonmiss_save[[j]]=dataset_nonmiss_imp
    select=which(dataset_nonmiss_imp$obstime <= imp_times_bp[j])
    dataset_nonmiss_select=dataset_nonmiss_imp[select,]

    fit5[[j]]=try(lmer(formula5,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=='TYPE5'))
    id_model5[[j]]=unique(dataset_nonmiss_select$id[dataset_nonmiss_select$bptype=='TYPE5' & !is.na(dataset_nonmiss_select$bp_curr)])
    cat('--------------------------------------------','\n')
    cat('BPTYPE=5','\n')
    print(summary(fit5[[j]]))
    cat('--------------------------------------------','\n')
    dataset_impute$bp_curr5[dataset_impute$start==imp_times_bp[j]]=predict(fit5[[j]],newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr) & dataset_nonmiss_imp$bptype=='TYPE1',])[dataset_impute$start==imp_times_bp[j]]
    if (imp_times_bp[j]==s) {
      dataset_impute$int5_s[dataset_impute$start==imp_times_bp[j]]=random.effects(fit5[[j]])[[1]][which(id_model5[[j]] %in% dataset_impute$id[dataset_impute$start==imp_times_bp[j]]),1]
      dataset_impute$slope5_s[dataset_impute$start==imp_times_bp[j]]=random.effects(fit5[[j]])[[1]][which(id_model5[[j]] %in% dataset_impute$id[dataset_impute$start==imp_times_bp[j]]),2]
    }
  }

  cat('--------------------------------------------','\n')

  j=j+1
}
cat('Finished fitting longitudinal models','\n')
#cat('nonsurvive_time for first two subjects in dataset_impute','\n')
#print(nonsurvive_time[eval_id[1:2]])
cat('--------------------------------------------','\n')
cat('Before EB Cox model','\n')
cat('dataset_impute after imputation from longitudinal models=','\n')
print(dataset_impute[1:200,])
#--------------------------------------------------------------------------
#
# END: Fit Longitudinal Models for Biomarker Process 
#
#--------------------------------------------------------------------------
#currentseed=.Random.seed[2]
#currentseed=.Random.seed[currentseed+2]
#cat('After Longitudinal Models fit current R seed=',currentseed,'\n')
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#
# Fit EMPIRICAL-BAYES COX MODEL
#
#--------------------------------------------------------------------------
dataset_cox=dataset_impute
dataset_cox=dataset_cox[dataset_cox$start < dataset_cox$stop-.000001,]

cat('EB COX1 Dataset_cox=','\n')
print(summary(dataset_cox))
cat('#subjects in  Dataset_cox=',length(unique(dataset_cox$id)),'\n')
cat('#events in  Dataset_cox=',sum(dataset_cox$status),'\n')
#print(dataset_cox[is.na(dataset_cox$bp_curr1),])

cat('------------------------------------------------------------','\n')
cat('EB Cox model uses apriori FIXED TERMS','\n')
cat('------------------------------------------------------------','\n')
if (dim_bp==1) {
  (ecformula=as.formula(paste("Surv(start,stop,status)~bp_curr1",base_names,sep="")))
}
if (dim_bp==2) {
  (ecformula=as.formula(paste("Surv(start,stop,status)~bp_curr1+bp_curr2",base_names,sep="")))
}
if (dim_bp==3) {
  (ecformula=as.formula(paste("Surv(start,stop,status)~bp_curr1+bp_curr2+bp_curr3",base_names,sep="")))
}
if (dim_bp==4) {
  (ecformula=as.formula(paste("Surv(start,stop,status)~bp_curr1+bp_curr2+bp_curr3+bp_curr4",base_names,sep="")))
}
if (dim_bp==5) {
  (ecformula=as.formula(paste("Surv(start,stop,status)~bp_curr1+bp_curr2+bp_curr3+bp_curr4+bp_curr5",base_names,sep="")))
}
efit=try(coxph(ecformula,data=dataset_cox,x=TRUE))

cat('------------------------------------------------------------','\n')
cat('EB Cox model with ',dim_bp,' biomarker terms is','\n')
print(efit)

coxcoefvalues=abs(coef(summary(efit))[1:dim_bp,1])
coxcoefvalues[is.na(coxcoefvalues)]=0.0
#--------------------------------------------------------------------------
#
#  Loop over Evaluable IDs
#
#--------------------------------------------------------------------------
# dataset_eval=transform(dataset_eval,int1_s=NA,slope1_s=NA,
#                                     int2_s=NA,slope2_s=NA,
#                                     int3_s=NA,slope3_s=NA,
#                                     int4_s=NA,slope4_s=NA,
#                                     int5_s=NA,slope5_s=NA)
# dataset_eval_save=dataset_eval
# 
# icox1_cond_surv_pred=rep(NA,times=n_eval)
# icox2_cond_surv_pred=rep(NA,times=n_eval)
# lvcf_cond_surv_pred=rep(NA,times=n_eval)
# icox1_boot_sd=rep(0,times=n_eval)
# icox2_boot_sd=rep(0,times=n_eval)
eval_id_num=start_eval_id_num
end_eval_id_num=min(end_eval_id_num,length(eval_id))
while (eval_id_num <= end_eval_id_num) {

  cat('','\n')
  cat('','\n')
  cat('----------------------------------------------------------','\n')
  cat('Evaluaton# =',eval_id_num,'\n')
  cat('Evaluate for ID=',eval_id[eval_id_num],'\n')
  cat('----------------------------------------------------------','\n')


  future_dataset_eval=dataset_eval_save[dataset_eval_save$id==eval_id[eval_id_num] & dataset_eval_save$start < s+tau,]
  jj=length(future_dataset_eval$id)
  future_dataset_eval$status=rep(0,times=jj)
  future_dataset_eval$stop[jj]=predtimes[jj]
  while (jj < length(predtimes)) {
    temp=future_dataset_eval[1,]
    temp$start[1]=predtimes[jj]
    temp$stop[1]=predtimes[jj+1]
    future_dataset_eval=rbind(future_dataset_eval,temp)
    jj=jj+1
  }
#  cat('Future dataset for ID=',eval_id[eval_id_num],'\n')
#  print(future_dataset_eval)

#-------------------------------------------------------------------------------
# Re-calculate nonsurvive_time for ID=eval_id[eval_id_num]
#-------------------------------------------------------------------------------
  nonsurvive_time=nonsurvive_time_save
#  nonsurvive_time[eval_id[eval_id_num]]=imp_times_bp[length(imp_times_bp)]
  nonsurvive_time[eval_id[eval_id_num]]=cens_time

#  cat('After entering loop over evaluable: future_dataset_eval for ID==eval_id[eval_id_num]','\n')
#  print(future_dataset_eval)
#  cat('After entering loop over evaluable: nonsurvive_time for ID==eval_id[eval_id_num]','\n')
#  print(nonsurvive_time[eval_id[eval_id_num]])


  if (icox1_ind==1) {
#-----------------------------------------------------
# START ICOX1: METHOD THAT USES BLUP OF FUTURE SUBJECT
#-----------------------------------------------------
    cat('---------------------------------------------','\n')
    cat('START ICOX1','\n')
    cat('---------------------------------------------','\n')
    cat('---------------------------------------------','\n')

#--------------------------------------------------------------------------
#
# Future Patient
# Replace ALL values of Biomarker Process with fitted value from model
#
#--------------------------------------------------------------------------

    future_dataset_impute=future_dataset_eval
    n_future=length(future_dataset_impute$id)

    j=1
    while (j <= num_imp_bp-1) {
      cat('Imp time=',imp_times_bp[j],'\n')
      cat('----------------------------------------','\n')
      cat('----------------------------------------','\n')
#      dataset_nonmiss_imp=dataset_nonmiss_save[[j]]
      dataset_nonmiss_imp=dataset_nonmiss_imp_save
      dataset_nonmiss_save=list()
#-------------------------------------------------------------------------------
#  If not already there add obs for ID=eval_id[eval_id_num] to dataset_nonmiss_imp with bp_curr=NA and obstime=imp_time_bp[j]
#-------------------------------------------------------------------------------
      temp=dataset_nonmiss_imp[dataset_nonmiss_imp$id==eval_id[eval_id_num],][1,]
      temp$bp_curr[1]=NA
      temp$obstime[1]=imp_times_bp[j]
      temp1=dataset_nonmiss_imp[dataset_nonmiss_imp$id==eval_id[eval_id_num] & dataset_nonmiss_imp$obstime==imp_times_bp[j] & is.na(dataset_nonmiss_imp$bp_curr),]
      if (length(temp1$id)==0) {
        dataset_nonmiss_imp=rbind(dataset_nonmiss_imp,temp) 
      }
#-------------------------------------------------------------------------------
# Restrict to data before s for i=eval_id[eval_id_num]
#-------------------------------------------------------------------------------
      dataset_nonmiss_imp=dataset_nonmiss_imp[dataset_nonmiss_imp$id != eval_id[eval_id_num] | dataset_nonmiss_imp$obstime <= s| is.na(dataset_nonmiss_imp$bp_curr),]

  nspline[1]=max(1,min(floor(imp_times_bp[j]/spline_width[1]),nspline_max[1]))
  if (dim_bp >=2) {nspline[2]=max(1,min(floor(imp_times_bp[j]/spline_width[2]),nspline_max[2]))}
  if (dim_bp >=3) {nspline[3]=max(1,min(floor(imp_times_bp[j]/spline_width[3]),nspline_max[3]))}
  if (dim_bp >=4) {nspline[4]=max(1,min(floor(imp_times_bp[j]/spline_width[4]),nspline_max[4]))}
  if (dim_bp >=5) {nspline[5]=max(1,min(floor(imp_times_bp[j]/spline_width[5]),nspline_max[5]))}
  cat('--------------------------------------------','\n')
  cat('--------------------------------------------','\n')
  cat('Imp model for imp time =',imp_times_bp[j],'\n')
  cat('--------------------------------------------','\n')
  cat('nspline[1] for Biomarker #1 =',nspline[1],'\n')
  cat('nspline[2] for Biomarker #2 =',nspline[2],'\n')
  cat('nspline[3] for Biomarker #3 =',nspline[3],'\n')

#--------------------------------------------------------------------------
# Create Ordered Formulas for Longitudinal Models of Biomarker Process 
#--------------------------------------------------------------------------
  jj=1
  base_names=""
  while (jj <= ncov) {
    if (base_cov_factor[jj]==0) {
      base_names=paste(base_names,"+base_cov.",jj,sep="")
    } else {
      base_names=paste(base_names,"+as.factor(base_cov.",jj,")",sep="")
    }
    jj=jj+1
  }
  names=base_names
  if (nspline[1]>=1) {
    names=paste(names,"+spline111+spline112+spline113",sep="")
  }
  if (nspline[1]>=2) {
    names=paste(names,"+spline121+spline122+spline123",sep="")
  }
  if (nspline[1]>=3) {
    names=paste(names,"+spline131+spline132+spline133",sep="")
  }
  if (nspline[1]>=4) {
    names=paste(names,"+spline141+spline142+spline143",sep="")
  }
  if (nspline[1]>=5) {
    names=paste(names,"+spline151+spline152+spline153",sep="")
  }
  if (nspline[1]>=6) {
    names=paste(names,"+spline161+spline162+spline163",sep="")
  }
  if (nspline[1]>=7) {
    names=paste(names,"+spline171+spline172+spline173",sep="")
  }
  if (nspline[1]>=8) {
    names=paste(names,"+spline181+spline182+spline183",sep="")
  }

  names=paste(names,"+nonsurvive+nonsurvive1+ (obstime | id)",sep="")
#  cat('names =',names,'\n')
  (formula1=as.formula(paste("bp_curr~1+",names,sep="")))

  if (dim_bp >=2) {
    names=base_names
    if (nspline[2]>=1) {
      names=paste(names,"+spline211+spline212+spline213",sep="")
    }
    if (nspline[2]>=2) {
      names=paste(names,"+spline221+spline222+spline223",sep="")
    }
    if (nspline[2]>=3) {
      names=paste(names,"+spline231+spline232+spline233",sep="")
    }
    if (nspline[2]>=4) {
      names=paste(names,"+spline241+spline242+spline243",sep="")
    }
    if (nspline[2]>=5) {
      names=paste(names,"+spline251+spline252+spline253",sep="")
    }
    if (nspline[2]>=6) {
      names=paste(names,"+spline261+spline262+spline263",sep="")
    }
    if (nspline[2]>=7) {
      names=paste(names,"+spline271+spline272+spline273",sep="")
    }
    if (nspline[2]>=8) {
      names=paste(names,"+spline281+spline282+spline283",sep="")
    }

    names=paste(names,"+nonsurvive+nonsurvive1+ (obstime | id)",sep="")
#    cat('names =',names,'\n')
    (formula2=as.formula(paste("bp_curr~1+",names,sep="")))
  }

  if (dim_bp >=3) {
    names=base_names
    if (nspline[3]>=1) {
      names=paste(names,"+spline311+spline312+spline313",sep="")
    }
    if (nspline[3]>=2) {
      names=paste(names,"+spline321+spline322+spline323",sep="")
    }
    if (nspline[3]>=3) {
      names=paste(names,"+spline331+spline332+spline333",sep="")
    }
    if (nspline[3]>=4) {
      names=paste(names,"+spline341+spline342+spline343",sep="")
    }
    if (nspline[3]>=5) {
      names=paste(names,"+spline351+spline352+spline353",sep="")
    }
    if (nspline[3]>=6) {
      names=paste(names,"+spline361+spline362+spline363",sep="")
    }
    if (nspline[3]>=7) {
      names=paste(names,"+spline371+spline372+spline373",sep="")
    }
    if (nspline[3]>=8) {
      names=paste(names,"+spline381+spline382+spline383",sep="")
    }

    names=paste(names,"+nonsurvive+nonsurvive1+ (obstime | id)",sep="")
#    cat('names =',names,'\n')
    (formula3=as.formula(paste("bp_curr~1+",names,sep="")))
  }

  if (dim_bp >=4) {
    names=base_names
    if (nspline[4]>=1) {
      names=paste(names,"+spline411+spline412+spline413",sep="")
    }
    if (nspline[4]>=2) {
      names=paste(names,"+spline421+spline422+spline423",sep="")
    }
    if (nspline[4]>=3) {
      names=paste(names,"+spline431+spline432+spline433",sep="")
    } 
    if (nspline[4]>=4) {
      names=paste(names,"+spline441+spline442+spline443",sep="")
    }
    if (nspline[4]>=5) {
      names=paste(names,"+spline451+spline452+spline453",sep="")
    }
    if (nspline[4]>=6) {
      names=paste(names,"+spline461+spline462+spline463",sep="")
    }
    if (nspline[4]>=7) {
      names=paste(names,"+spline471+spline472+spline473",sep="")
    }
    if (nspline[4]>=8) {
      names=paste(names,"+spline481+spline482+spline483",sep="")
    }

    names=paste(names,"+nonsurvive+nonsurvive1+ (obstime | id)",sep="")
#    cat('names =',names,'\n')
    (formula4=as.formula(paste("bp_curr~1+",names,sep="")))
  }

  if (dim_bp >=5) {
    names=base_names
    if (nspline[5]>=1) {
      names=paste(names,"+spline511+spline512+spline513",sep="")
    }
    if (nspline[5]>=2) {
      names=paste(names,"+spline521+spline522+spline523",sep="")
    }
    if (nspline[5]>=3) {
      names=paste(names,"+spline531+spline532+spline533",sep="")
    } 
    if (nspline[5]>=4) {
      names=paste(names,"+spline541+spline542+spline543",sep="")
    }
    if (nspline[5]>=5) {
      names=paste(names,"+spline551+spline552+spline553",sep="")
    }
    if (nspline[5]>=6) {
      names=paste(names,"+spline561+spline562+spline563",sep="")
    }
    if (nspline[5]>=7) {
      names=paste(names,"+spline571+spline572+spline573",sep="")
    }
    if (nspline[5]>=8) {
      names=paste(names,"+spline581+spline582+spline583",sep="")
    }

    names=paste(names,"+nonsurvive+nonsurvive1+ (obstime | id)",sep="")
#    cat('names =',names,'\n')
    (formula5=as.formula(paste("bp_curr~1+",names,sep="")))
  }

#--------------------------------------------------------------------------
# Spline Creation 
#--------------------------------------------------------------------------
  knot1=spline_width[1]
  knot2=2*spline_width[1]
  knot3=3*spline_width[1]
  knot4=4*spline_width[1]
  knot5=5*spline_width[1]
  knot6=6*spline_width[1]
  knot7=7*spline_width[1]
  dataset_nonmiss_imp$spline111=dataset_nonmiss_imp$obstime
  dataset_nonmiss_imp$spline112=dataset_nonmiss_imp$obstime**2
  dataset_nonmiss_imp$spline113=dataset_nonmiss_imp$obstime**3
  if (nspline[1] >= 2) {
    dataset_nonmiss_imp$spline111=dataset_nonmiss_imp$obstime**1*(dataset_nonmiss_imp$obstime < knot1) +knot1**1*(dataset_nonmiss_imp$obstime >= knot1)
    dataset_nonmiss_imp$spline112=dataset_nonmiss_imp$obstime**2*(dataset_nonmiss_imp$obstime < knot1) +knot1**2*(dataset_nonmiss_imp$obstime >= knot1)
    dataset_nonmiss_imp$spline113=dataset_nonmiss_imp$obstime**3*(dataset_nonmiss_imp$obstime < knot1) +knot1**3*(dataset_nonmiss_imp$obstime >= knot1)
    dataset_nonmiss_imp$spline121=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime >= knot1) 
    dataset_nonmiss_imp$spline122=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime >= knot1) 
    dataset_nonmiss_imp$spline123=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime >= knot1) 
  }
  if (nspline[1] >= 3) {
    dataset_nonmiss_imp$spline121=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**1*(dataset_nonmiss_imp$obstime >= knot2)
    dataset_nonmiss_imp$spline122=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**2*(dataset_nonmiss_imp$obstime >= knot2)
    dataset_nonmiss_imp$spline123=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**3*(dataset_nonmiss_imp$obstime >= knot2)
    dataset_nonmiss_imp$spline131=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime >= knot2) 
    dataset_nonmiss_imp$spline132=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime >= knot2) 
    dataset_nonmiss_imp$spline133=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime >= knot2) 
  }
  if (nspline[1] >= 4) {
    dataset_nonmiss_imp$spline131=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**1*(dataset_nonmiss_imp$obstime >= knot3)
    dataset_nonmiss_imp$spline132=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**2*(dataset_nonmiss_imp$obstime >= knot3)
    dataset_nonmiss_imp$spline133=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**3*(dataset_nonmiss_imp$obstime >= knot3)
    dataset_nonmiss_imp$spline141=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime >= knot3) 
    dataset_nonmiss_imp$spline142=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime >= knot3) 
    dataset_nonmiss_imp$spline143=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime >= knot3) 
  }
  if (nspline[1] >= 5) {
    dataset_nonmiss_imp$spline141=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**1*(dataset_nonmiss_imp$obstime >= knot4)
    dataset_nonmiss_imp$spline142=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**2*(dataset_nonmiss_imp$obstime >= knot4)
    dataset_nonmiss_imp$spline143=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**3*(dataset_nonmiss_imp$obstime >= knot4)
    dataset_nonmiss_imp$spline151=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime >= knot4) 
    dataset_nonmiss_imp$spline152=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime >= knot4) 
    dataset_nonmiss_imp$spline153=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime >= knot4) 
  }
  if (nspline[1] >= 6) {
    dataset_nonmiss_imp$spline151=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**1*(dataset_nonmiss_imp$obstime >= knot5)
    dataset_nonmiss_imp$spline152=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**2*(dataset_nonmiss_imp$obstime >= knot5)
    dataset_nonmiss_imp$spline153=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**3*(dataset_nonmiss_imp$obstime >= knot5)
    dataset_nonmiss_imp$spline161=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime >= knot5) 
    dataset_nonmiss_imp$spline162=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime >= knot5) 
    dataset_nonmiss_imp$spline163=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime >= knot5) 
  }
  if (nspline[1] >= 7) {
    dataset_nonmiss_imp$spline161=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**1*(dataset_nonmiss_imp$obstime >= knot6)
    dataset_nonmiss_imp$spline162=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**2*(dataset_nonmiss_imp$obstime >= knot6)
    dataset_nonmiss_imp$spline163=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**3*(dataset_nonmiss_imp$obstime >= knot6)
    dataset_nonmiss_imp$spline171=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime >= knot6) 
    dataset_nonmiss_imp$spline172=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime >= knot6) 
    dataset_nonmiss_imp$spline173=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime >= knot6) 
  }
  if (nspline[1] >= 8) {
    dataset_nonmiss_imp$spline171=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**1*(dataset_nonmiss_imp$obstime >= knot7)
    dataset_nonmiss_imp$spline172=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**2*(dataset_nonmiss_imp$obstime >= knot7)
    dataset_nonmiss_imp$spline173=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**3*(dataset_nonmiss_imp$obstime >= knot7)
    dataset_nonmiss_imp$spline181=(dataset_nonmiss_imp$obstime-knot7)**1*(dataset_nonmiss_imp$obstime >= knot7) 
    dataset_nonmiss_imp$spline182=(dataset_nonmiss_imp$obstime-knot7)**2*(dataset_nonmiss_imp$obstime >= knot7) 
    dataset_nonmiss_imp$spline183=(dataset_nonmiss_imp$obstime-knot7)**3*(dataset_nonmiss_imp$obstime >= knot7) 
  }

#--------------------------------------------------------------------------
# END: Spline Creation 
#--------------------------------------------------------------------------
  dataset_nonmiss_imp$nonsurvive=(nonsurvive_time[dataset_nonmiss_imp$id] <= imp_times_bp[j])
  dataset_nonmiss_imp$nonsurvive1=dataset_nonmiss_imp$nonsurvive*dataset_nonmiss_imp$obstime
  dataset_nonmiss_save[[1]]=dataset_nonmiss_imp

  select=which(dataset_nonmiss_imp$obstime <= imp_times_bp[j])
  dataset_nonmiss_select=dataset_nonmiss_imp[select,]

#  cat('dataset_nonmiss_select =','\n')
#  print(dataset_nonmiss_select[1:80,])
  cat('# unique ID in dataset_nonmiss_select =',length(unique(dataset_nonmiss_select$id)),'\n')
#  cat('dataset_nonmiss_imp =','\n')
#  print(dataset_nonmiss_imp[1:80,])
  cat('# unique ID in dataset_nonmiss_imp =',length(unique(dataset_nonmiss_imp$id)),'\n')

#-------------------------------------------------
# FIT LINEAR FIXED EFFECTS
  fit1[[j]]=try(lmer(formula1,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=='TYPE1'))
  cat('--------------------------------------------','\n')
  cat('BPTYPE=1','\n')
#  print(summary(fit1[[j]]))
  cat('--------------------------------------------','\n')

#  cat('Biomarker #1: dataset_nonmiss_select for ID=',eval_id[eval_id_num],'\n')
#  print(dataset_nonmiss_select[dataset_nonmiss_select$id==eval_id[eval_id_num],])
#  cat('----------------------------------------','\n')

  if (dim_bp>=2) {

#  SPLINE CREATION
    knot1=spline_width[2]
    knot2=2*spline_width[2]
    knot3=3*spline_width[2]
    knot4=4*spline_width[2]
    knot5=5*spline_width[2]
    knot6=6*spline_width[2]
    knot7=7*spline_width[2]
    dataset_nonmiss_imp$spline211=dataset_nonmiss_imp$obstime
    dataset_nonmiss_imp$spline212=dataset_nonmiss_imp$obstime**2
    dataset_nonmiss_imp$spline213=dataset_nonmiss_imp$obstime**3
    if (nspline[2] >= 2) {
      dataset_nonmiss_imp$spline211=dataset_nonmiss_imp$obstime**1*(dataset_nonmiss_imp$obstime < knot1) +knot1**1*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline212=dataset_nonmiss_imp$obstime**2*(dataset_nonmiss_imp$obstime < knot1) +knot1**2*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline213=dataset_nonmiss_imp$obstime**3*(dataset_nonmiss_imp$obstime < knot1) +knot1**3*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline221=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline222=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline223=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime >= knot1) 
    }
    if (nspline[2] >= 3) {
      dataset_nonmiss_imp$spline221=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**1*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline222=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**2*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline223=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**3*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline231=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline232=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline233=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime >= knot2) 
    }
    if (nspline[2] >= 4) {
      dataset_nonmiss_imp$spline231=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**1*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline232=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**2*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline233=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**3*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline241=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline242=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline243=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime >= knot3) 
    }
    if (nspline[2] >= 5) {
      dataset_nonmiss_imp$spline241=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**1*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline242=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**2*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline243=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**3*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline251=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline252=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline253=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime >= knot4) 
    }
    if (nspline[2] >= 6) {
      dataset_nonmiss_imp$spline251=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**1*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline252=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**2*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline253=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**3*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline261=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline262=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline263=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime >= knot5) 
    }
    if (nspline[2] >= 7) {
      dataset_nonmiss_imp$spline261=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**1*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline262=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**2*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline263=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**3*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline271=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline272=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline273=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime >= knot6) 
    }
    if (nspline[2] >= 8) {
      dataset_nonmiss_imp$spline271=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**1*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline272=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**2*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline273=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**3*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline281=(dataset_nonmiss_imp$obstime-knot7)**1*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline282=(dataset_nonmiss_imp$obstime-knot7)**2*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline283=(dataset_nonmiss_imp$obstime-knot7)**3*(dataset_nonmiss_imp$obstime >= knot7) 
    }
    select=which(dataset_nonmiss_imp$obstime <= imp_times_bp[j])
    dataset_nonmiss_select=dataset_nonmiss_imp[select,]
    dataset_nonmiss_save[[2]]=dataset_nonmiss_imp

    fit2[[j]]=try(lmer(formula2,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=='TYPE2'))
    cat('--------------------------------------------','\n')
    cat('BPTYPE=2','\n')
#    print(summary(fit2[[j]]))
    cat('--------------------------------------------','\n')
  }

  if (dim_bp>=3) {

#  SPLINE CREATION
    knot1=spline_width[3]
    knot2=2*spline_width[3]
    knot3=3*spline_width[3]
    knot4=4*spline_width[3]
    knot5=5*spline_width[3]
    knot6=6*spline_width[3]
    knot7=7*spline_width[3]
    dataset_nonmiss_imp$spline311=dataset_nonmiss_imp$obstime
    dataset_nonmiss_imp$spline312=dataset_nonmiss_imp$obstime**2
    dataset_nonmiss_imp$spline313=dataset_nonmiss_imp$obstime**3
    if (nspline[3] >= 2) {
      dataset_nonmiss_imp$spline311=dataset_nonmiss_imp$obstime**1*(dataset_nonmiss_imp$obstime < knot1) +knot1**1*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline312=dataset_nonmiss_imp$obstime**2*(dataset_nonmiss_imp$obstime < knot1) +knot1**2*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline313=dataset_nonmiss_imp$obstime**3*(dataset_nonmiss_imp$obstime < knot1) +knot1**3*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline321=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline322=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline323=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime >= knot1) 
    }
    if (nspline[3] >= 3) {
      dataset_nonmiss_imp$spline321=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**1*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline322=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**2*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline323=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**3*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline331=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline332=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline333=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime >= knot2) 
    }
    if (nspline[3] >= 4) {
      dataset_nonmiss_imp$spline331=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**1*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline332=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**2*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline333=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**3*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline341=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline342=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline343=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime >= knot3) 
    }
    if (nspline[3] >= 5) {
      dataset_nonmiss_imp$spline341=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**1*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline342=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**2*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline343=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**3*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline351=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline352=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline353=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime >= knot4) 
    }
    if (nspline[3] >= 6) {
      dataset_nonmiss_imp$spline351=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**1*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline352=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**2*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline353=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**3*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline361=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline362=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline363=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime >= knot5) 
    }
    if (nspline[3] >= 7) {
      dataset_nonmiss_imp$spline361=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**1*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline362=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**2*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline363=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**3*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline371=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline372=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline373=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime >= knot6) 
    }
    if (nspline[3] >= 8) {
      dataset_nonmiss_imp$spline371=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**1*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline372=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**2*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline373=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**3*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline381=(dataset_nonmiss_imp$obstime-knot7)**1*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline382=(dataset_nonmiss_imp$obstime-knot7)**2*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline383=(dataset_nonmiss_imp$obstime-knot7)**3*(dataset_nonmiss_imp$obstime >= knot7) 
    }
    select=which(dataset_nonmiss_imp$obstime <= imp_times_bp[j])
    dataset_nonmiss_select=dataset_nonmiss_imp[select,]
    dataset_nonmiss_save[[3]]=dataset_nonmiss_imp

    fit3[[j]]=try(lmer(formula3,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=='TYPE3'))
    cat('--------------------------------------------','\n')
    cat('BPTYPE=3','\n')
#    print(summary(fit3[[j]]))
    cat('--------------------------------------------','\n')
  }

  if (dim_bp>=4) {

#  SPLINE CREATION
    knot1=spline_width[4]
    knot2=2*spline_width[4]
    knot3=3*spline_width[4]
    knot4=4*spline_width[4]
    knot5=5*spline_width[4]
    knot6=6*spline_width[4]
    knot7=7*spline_width[4]
    dataset_nonmiss_imp$spline411=dataset_nonmiss_imp$obstime
    dataset_nonmiss_imp$spline412=dataset_nonmiss_imp$obstime**2
    dataset_nonmiss_imp$spline413=dataset_nonmiss_imp$obstime**3
    if (nspline[4] >= 2) {
      dataset_nonmiss_imp$spline411=dataset_nonmiss_imp$obstime**1*(dataset_nonmiss_imp$obstime < knot1) +knot1**1*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline412=dataset_nonmiss_imp$obstime**2*(dataset_nonmiss_imp$obstime < knot1) +knot1**2*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline413=dataset_nonmiss_imp$obstime**3*(dataset_nonmiss_imp$obstime < knot1) +knot1**3*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline421=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline422=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline423=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime >= knot1) 
    }
    if (nspline[4] >= 3) {
      dataset_nonmiss_imp$spline421=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**1*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline422=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**2*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline423=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**3*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline431=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline432=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline433=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime >= knot2) 
    }
    if (nspline[4] >= 4) {
      dataset_nonmiss_imp$spline431=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**1*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline432=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**2*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline433=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**3*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline441=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline442=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline443=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime >= knot3) 
    }
    if (nspline[4] >= 5) {
      dataset_nonmiss_imp$spline441=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**1*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline442=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**2*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline443=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**3*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline451=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline452=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline453=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime >= knot4) 
    }
    if (nspline[4] >= 6) {
      dataset_nonmiss_imp$spline451=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**1*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline452=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**2*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline453=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**3*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline461=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline462=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline463=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime >= knot5) 
    }
    if (nspline[4] >= 7) {
      dataset_nonmiss_imp$spline461=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**1*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline462=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**2*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline463=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**3*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline471=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline472=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline473=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime >= knot6) 
    }
    if (nspline[4] >= 8) {
      dataset_nonmiss_imp$spline471=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**1*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline472=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**2*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline473=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**3*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline481=(dataset_nonmiss_imp$obstime-knot7)**1*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline482=(dataset_nonmiss_imp$obstime-knot7)**2*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline483=(dataset_nonmiss_imp$obstime-knot7)**3*(dataset_nonmiss_imp$obstime >= knot7) 
    }
    select=which(dataset_nonmiss_imp$obstime <= imp_times_bp[j])
    dataset_nonmiss_select=dataset_nonmiss_imp[select,]
    dataset_nonmiss_save[[4]]=dataset_nonmiss_imp

    fit4[[j]]=try(lmer(formula4,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=='TYPE4'))
    cat('--------------------------------------------','\n')
    cat('BPTYPE=4','\n')
#    print(summary(fit4[[j]]))
    cat('--------------------------------------------','\n')
  }

  if (dim_bp>=5) {

#  SPLINE CREATION
    knot1=spline_width[5]
    knot2=2*spline_width[5]
    knot3=3*spline_width[5]
    knot4=4*spline_width[5]
    knot5=5*spline_width[5]
    knot6=6*spline_width[5]
    knot7=7*spline_width[5]
    dataset_nonmiss_imp$spline511=dataset_nonmiss_imp$obstime
    dataset_nonmiss_imp$spline512=dataset_nonmiss_imp$obstime**2
    dataset_nonmiss_imp$spline513=dataset_nonmiss_imp$obstime**3
    if (nspline[5] >= 2) {
      dataset_nonmiss_imp$spline511=dataset_nonmiss_imp$obstime**1*(dataset_nonmiss_imp$obstime < knot1) +knot1**1*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline512=dataset_nonmiss_imp$obstime**2*(dataset_nonmiss_imp$obstime < knot1) +knot1**2*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline513=dataset_nonmiss_imp$obstime**3*(dataset_nonmiss_imp$obstime < knot1) +knot1**3*(dataset_nonmiss_imp$obstime >= knot1)
      dataset_nonmiss_imp$spline521=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline522=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime >= knot1) 
      dataset_nonmiss_imp$spline523=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime >= knot1) 
    }
    if (nspline[5] >= 3) {
      dataset_nonmiss_imp$spline521=(dataset_nonmiss_imp$obstime-knot1)**1*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**1*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline522=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**2*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline523=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**3*(dataset_nonmiss_imp$obstime >= knot2)
      dataset_nonmiss_imp$spline531=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline532=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime >= knot2) 
      dataset_nonmiss_imp$spline533=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime >= knot2) 
    }
    if (nspline[5] >= 4) {
      dataset_nonmiss_imp$spline531=(dataset_nonmiss_imp$obstime-knot2)**1*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**1*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline532=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**2*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline533=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**3*(dataset_nonmiss_imp$obstime >= knot3)
      dataset_nonmiss_imp$spline541=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline542=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime >= knot3) 
      dataset_nonmiss_imp$spline543=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime >= knot3) 
    }
    if (nspline[5] >= 5) {
      dataset_nonmiss_imp$spline541=(dataset_nonmiss_imp$obstime-knot3)**1*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**1*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline542=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**2*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline543=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**3*(dataset_nonmiss_imp$obstime >= knot4)
      dataset_nonmiss_imp$spline551=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline552=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime >= knot4) 
      dataset_nonmiss_imp$spline553=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime >= knot4) 
    }
    if (nspline[5] >= 6) {
      dataset_nonmiss_imp$spline551=(dataset_nonmiss_imp$obstime-knot4)**1*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**1*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline552=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**2*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline553=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**3*(dataset_nonmiss_imp$obstime >= knot5)
      dataset_nonmiss_imp$spline561=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline562=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime >= knot5) 
      dataset_nonmiss_imp$spline563=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime >= knot5) 
    }
    if (nspline[5] >= 7) {
      dataset_nonmiss_imp$spline561=(dataset_nonmiss_imp$obstime-knot5)**1*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**1*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline562=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**2*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline563=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**3*(dataset_nonmiss_imp$obstime >= knot6)
      dataset_nonmiss_imp$spline571=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline572=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime >= knot6) 
      dataset_nonmiss_imp$spline573=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime >= knot6) 
    }
    if (nspline[5] >= 8) {
      dataset_nonmiss_imp$spline571=(dataset_nonmiss_imp$obstime-knot6)**1*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**1*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline572=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**2*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline573=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**3*(dataset_nonmiss_imp$obstime >= knot7)
      dataset_nonmiss_imp$spline581=(dataset_nonmiss_imp$obstime-knot7)**1*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline582=(dataset_nonmiss_imp$obstime-knot7)**2*(dataset_nonmiss_imp$obstime >= knot7) 
      dataset_nonmiss_imp$spline583=(dataset_nonmiss_imp$obstime-knot7)**3*(dataset_nonmiss_imp$obstime >= knot7) 
    }
    select=which(dataset_nonmiss_imp$obstime <= imp_times_bp[j])
    dataset_nonmiss_select=dataset_nonmiss_imp[select,]
    dataset_nonmiss_save[[5]]=dataset_nonmiss_imp

    fit5[[j]]=try(lmer(formula5,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=='TYPE5'))
    cat('--------------------------------------------','\n')
    cat('BPTYPE=5','\n')
#    print(summary(fit5[[j]]))
    cat('--------------------------------------------','\n')
  }


      future_dataset_nonmiss=dataset_nonmiss_save[[1]]
#      future_dataset_nonmiss=future_dataset_nonmiss[future_dataset_nonmiss$id==eval_id[eval_id_num] & future_dataset_nonmiss$obstime <= s,]
      future_dataset_nonmiss=future_dataset_nonmiss[future_dataset_nonmiss$id==eval_id[eval_id_num] & is.na(future_dataset_nonmiss$bp_curr),]
      future_dataset_nonmiss=future_dataset_nonmiss[future_dataset_nonmiss$obstime==imp_times_bp[j] & is.na(future_dataset_nonmiss$bp_curr),]

      future_dataset_nonmiss$obstime[1]=future_dataset_eval$start[1]
      future_dataset_nonmiss$bp_curr[1]=NA
      future_dataset_nonmiss$nonsurvive[1]=FALSE
      future_dataset_nonmiss$nonsurvive1[1]=0
      jj=2
      while (jj <= length(future_dataset_eval$start)) {
        temp=future_dataset_nonmiss[1,]
        temp$obstime[1]=future_dataset_eval$start[jj]
        temp$nonsurvive[1]=FALSE
        temp$nonsurvive1[1]=0
        future_dataset_nonmiss=rbind(future_dataset_nonmiss,temp)
        jj=jj+1
      }
#      cat('Biomarker #1: Future_dataset_nonmiss used for newdata for ID=',eval_id[eval_id_num],'\n')
#      print(future_dataset_nonmiss)
#      cat('----------------------------------------','\n')

      future_dataset_impute$bp_curr1[future_dataset_impute$start==imp_times_bp[j]]=predict(fit1[[j]],newdata=future_dataset_nonmiss)[future_dataset_impute$start==imp_times_bp[j]]
#      cat('Biomarker #1: Future_dataset_impute for ID=',eval_id[eval_id_num],'\n')
#      print(future_dataset_impute)
#      cat('----------------------------------------','\n')
      if (dim_bp >=2) {
        future_dataset_nonmiss=dataset_nonmiss_save[[2]]
        future_dataset_nonmiss=future_dataset_nonmiss[future_dataset_nonmiss$id==eval_id[eval_id_num] & is.na(future_dataset_nonmiss$bp_curr),]
        future_dataset_nonmiss=future_dataset_nonmiss[future_dataset_nonmiss$obstime==imp_times_bp[j] & is.na(future_dataset_nonmiss$bp_curr),]

        future_dataset_nonmiss$obstime[1]=future_dataset_eval$start[1]
        future_dataset_nonmiss$bp_curr[1]=NA
        future_dataset_nonmiss$nonsurvive[1]=FALSE
        future_dataset_nonmiss$nonsurvive1[1]=0
        jj=2
        while (jj <= length(future_dataset_eval$start)) {
          temp=future_dataset_nonmiss[1,]
          temp$obstime[1]=future_dataset_eval$start[jj]
          temp$nonsurvive[1]=FALSE
          temp$nonsurvive1[1]=0
          future_dataset_nonmiss=rbind(future_dataset_nonmiss,temp)
          jj=jj+1
        }
        future_dataset_impute$bp_curr2[future_dataset_impute$start==imp_times_bp[j]]=predict(fit2[[j]],newdata=future_dataset_nonmiss)[future_dataset_impute$start==imp_times_bp[j]]
      }
      if (dim_bp >=3) {
        future_dataset_nonmiss=dataset_nonmiss_save[[3]]
        future_dataset_nonmiss=future_dataset_nonmiss[future_dataset_nonmiss$id==eval_id[eval_id_num] & is.na(future_dataset_nonmiss$bp_curr),]
        future_dataset_nonmiss=future_dataset_nonmiss[future_dataset_nonmiss$obstime==imp_times_bp[j] & is.na(future_dataset_nonmiss$bp_curr),]

        future_dataset_nonmiss$obstime[1]=future_dataset_eval$start[1]
        future_dataset_nonmiss$bp_curr[1]=NA
        future_dataset_nonmiss$nonsurvive[1]=FALSE
        future_dataset_nonmiss$nonsurvive1[1]=0
        jj=2
        while (jj <= length(future_dataset_eval$start)) {
          temp=future_dataset_nonmiss[1,]
          temp$obstime[1]=future_dataset_eval$start[jj]
          temp$nonsurvive[1]=FALSE
          temp$nonsurvive1[1]=0
          future_dataset_nonmiss=rbind(future_dataset_nonmiss,temp)
          jj=jj+1
        }
        future_dataset_impute$bp_curr3[future_dataset_impute$start==imp_times_bp[j]]=predict(fit3[[j]],newdata=future_dataset_nonmiss)[future_dataset_impute$start==imp_times_bp[j]]
      }
      if (dim_bp >=4) {
        future_dataset_nonmiss=dataset_nonmiss_save[[4]]
        future_dataset_nonmiss=future_dataset_nonmiss[future_dataset_nonmiss$id==eval_id[eval_id_num] & is.na(future_dataset_nonmiss$bp_curr),]
        future_dataset_nonmiss=future_dataset_nonmiss[future_dataset_nonmiss$obstime==imp_times_bp[j] & is.na(future_dataset_nonmiss$bp_curr),]

        future_dataset_nonmiss$obstime[1]=future_dataset_eval$start[1]
        future_dataset_nonmiss$bp_curr[1]=NA
        future_dataset_nonmiss$nonsurvive[1]=FALSE
        future_dataset_nonmiss$nonsurvive1[1]=0
        jj=2
        while (jj <= length(future_dataset_eval$start)) {
          temp=future_dataset_nonmiss[1,]
          temp$obstime[1]=future_dataset_eval$start[jj]
          temp$nonsurvive[1]=FALSE
          temp$nonsurvive1[1]=0
          future_dataset_nonmiss=rbind(future_dataset_nonmiss,temp)
          jj=jj+1
        }
        future_dataset_impute$bp_curr4[future_dataset_impute$start==imp_times_bp[j]]=predict(fit4[[j]],newdata=future_dataset_nonmiss)[future_dataset_impute$start==imp_times_bp[j]]
      }
      if (dim_bp >=5) {
        future_dataset_nonmiss=dataset_nonmiss_save[[5]]
        future_dataset_nonmiss=future_dataset_nonmiss[future_dataset_nonmiss$id==eval_id[eval_id_num] & is.na(future_dataset_nonmiss$bp_curr),]
        future_dataset_nonmiss=future_dataset_nonmiss[future_dataset_nonmiss$obstime==imp_times_bp[j] & is.na(future_dataset_nonmiss$bp_curr),]

        future_dataset_nonmiss$obstime[1]=future_dataset_eval$start[1]
        future_dataset_nonmiss$bp_curr[1]=NA
        future_dataset_nonmiss$nonsurvive[1]=FALSE
        future_dataset_nonmiss$nonsurvive1[1]=0
        jj=2
        while (jj <= length(future_dataset_eval$start)) {
          temp=future_dataset_nonmiss[1,]
          temp$obstime[1]=future_dataset_eval$start[jj]
          temp$nonsurvive[1]=FALSE
          temp$nonsurvive1[1]=0
          future_dataset_nonmiss=rbind(future_dataset_nonmiss,temp)
          jj=jj+1
        }
        future_dataset_impute$bp_curr5[future_dataset_impute$start==imp_times_bp[j]]=predict(fit5[[j]],newdata=future_dataset_nonmiss)[future_dataset_impute$start==imp_times_bp[j]]
      }
      j=j+1
    }

    cat('FINAL future_dataset_nonmiss for ID=',eval_id[eval_id_num],'\n')
    print(future_dataset_nonmiss)
    cat('----------------------------------------','\n')
    cat('----------------------------------------','\n')
    cat('nonsurvive_time for ID=',eval_id[eval_id_num],'\n')
    print(nonsurvive_time[eval_id[eval_id_num]])
    cat('--------------------------------------------------','\n')
    cat('After EB Cox','\n')
    cat('After imputation: future_dataset_impute=','\n')
    print(future_dataset_impute)

#--------------------------------------------------------------------------
#
# FORTRAN VERSION: CALCULATE CONDITIONAL SURVIVAL PROBS
#
#--------------------------------------------------------------------------
    cond_surv_pred_st=0
    cat('Call predictSurvProb','\n')
    future_surv=predictSurvProb(efit,newdata=future_dataset_impute,times=predtimes)
    cat('Return from predictSurvProb','\n')
#  if (!all(is.finite(future_surv))) {
#    cat('survival to predtimes for subjects in future_dataset according to EMPIRICAL COX MODEL=','\n')
#    print(future_surv[1:20,])
#  }
#  future_surv[future_surv > 1]=NA
#  future_surv[future_surv < 0]=NA
    npredtimes=length(predtimes)
    future_surv_short=future_surv
#  future_surv_short=future_surv[1:n_future,]
#  if (all(is.finite(future_surv_short))) {
#    cat('All future_surv_short are finite','\n')
#  } else {
#    cat('survival to predtimes for subjects in future_dataset according to EMPIRICAL COX MODEL=','\n')
#    print(future_surv_short)
#  }


#test=future_dataset_impute_long
#test$bp.1[1]=999
#test_surv=predictSurvProb(efit[[k]],newdata=test,times=predtimes)
#test_surv_short=test_surv[1:n_future,]

#cat('survival to predtimes for subjects in test dataset according to EMPIRICAL COX MODEL=','\n')
#print(test_surv_short)


#  idvec=future_dataset_impute_long$id
    idvec=future_dataset_impute$id

#  cat('survival to predtimes for subjects in future_dataset according to EMPIRICAL COX MODEL=','\n')
#  print(future_surv_short)
#  cat('n_future=',n_future,'\n')
#  cat('npredtimes=',npredtimes,'\n')
#  cat('index_s=',index_s,'\n')
#  cat('cond_surv_pred_st=',cond_surv_pred_st,'\n')

    if (all(is.finite(future_surv_short))) {
      cat('Call condprob','\n')
      z=condprob(1,n_future,npredtimes,index_s,future_surv_short,idvec,cond_surv_pred_st)
      cat('return from condprob','\n')

      cond_surv_pred_st=z[[7]]

      cat('FORTRAN: Cond. survival to S+TAU given survival to S for subjects in future_dataset according to EB1 COX MODEL=','\n')
      print(cond_surv_pred_st)

      icox1_cond_surv_pred_imp=mean(cond_surv_pred_st)

      cat('FORTRAN: Avg. Cond. survival to S+TAU given survival to S for subjects in future_dataset according to EB1 COX MODEL=','\n')
      print(icox1_cond_surv_pred_imp)
    } else {
      cat('EB1 COX MODEL DOES NOT GIVE PROB ESTIMATE','\n')
      icox1_cond_surv_pred_imp=NA
    }

    cat('---------------------------------------------------','\n')
    cat('---------------------------------------------------','\n')
    if (is.finite(icox1_cond_surv_pred_imp)) {
      icox1_cond_surv_pred[eval_id_num]=mean(icox1_cond_surv_pred_imp,na.rm=TRUE)
      cat('FINAL Estimate Cond. survival to S+TAU given survival to S for future subj according to EB1 COX MODEL=','\n')
      print(icox1_cond_surv_pred[eval_id_num])
    } else {
      cat('FINAL EB1 COX MODEL DOES NOT GIVE PROB ESTIMATE','\n')
      icox1_cond_surv_pred[eval_id_num]=NA
    }

#eval_id_num=eval_id_num+1
#}


  }
#-----------------------------------------------------
# END ICOX1: METHOD THAT USES EBLUP OF FUTURE SUBJECT
#-----------------------------------------------------










  if (icox2_ind==1) {
#-------------------------------------------------------
# START ICOX2: METHOD THAT USES PERFECT MATCH PREDICTION
#-------------------------------------------------------
    cat('---------------------------------------------','\n')
    cat('START ICOX2','\n')
    cat('---------------------------------------------','\n')
    cat('---------------------------------------------','\n')

#    cat('------------------------------------','\n')
#    cat('dataset_nonmiss[ID==1,]=','\n')
#    print(dataset_nonmiss[dataset_nonmiss$id==1,])
#    cat('------------------------------------','\n')
#    cat('------------------------------------','\n')
#    cat('dataset_nonmiss[ID==10229,]=','\n')
#    print(dataset_nonmiss[dataset_nonmiss$id==10229,])
#    cat('------------------------------------','\n')
    
    temp_obs=list()
    temp_obs[[1]]=dataset_nonmiss[abs(dataset_nonmiss$obstime-design_times_bp[1]) <= window_width,]
    temp=temp_obs[[1]]
    jj=2
    while (jj <= num_design_times_bp) {
      temp_obs[[jj]]=dataset_nonmiss[abs(dataset_nonmiss$obstime-design_times_bp[jj]) <= window_width,]
      temp=rbind(temp,temp_obs[[jj]])
      jj=jj+1
    }
#    cat('------------------------------------','\n')
#    cat('temp[ID==1,]=','\n')
#    print(temp[temp$id==1,])
#    cat('------------------------------------','\n')
    
    current_last_obs_time1=temp$obstime[temp$id==eval_id[eval_id_num] & temp$obstime <= s & temp$bptype=="TYPE1"]
    current_last_obs_time1=current_last_obs_time1[length(current_last_obs_time1)]
    design_last_obs_time1=design_times_bp[design_times_bp <= current_last_obs_time1+window_width]
    design_last_obs_time1=design_last_obs_time1[length(design_last_obs_time1)]
    if (dim_bp >=2) {
      current_last_obs_time2=temp$obstime[temp$id==eval_id[eval_id_num] & temp$obstime <= s & temp$bptype=="TYPE2"]
      current_last_obs_time2=current_last_obs_time2[length(current_last_obs_time2)]
      design_last_obs_time2=design_times_bp[design_times_bp <= current_last_obs_time2+window_width]
      design_last_obs_time2=design_last_obs_time2[length(design_last_obs_time2)]
    }
    if (dim_bp >=3) {
      current_last_obs_time3=temp$obstime[temp$id==eval_id[eval_id_num] & temp$obstime <= s & temp$bptype=="TYPE3"]
      current_last_obs_time3=current_last_obs_time3[length(current_last_obs_time3)]
      design_last_obs_time3=design_times_bp[design_times_bp <= current_last_obs_time3+window_width]
      design_last_obs_time3=design_last_obs_time3[length(design_last_obs_time3)]
    }
    if (dim_bp >=4) {
      current_last_obs_time4=temp$obstime[temp$id==eval_id[eval_id_num] & temp$obstime <= s & temp$bptype=="TYPE4"]
      current_last_obs_time4=current_last_obs_time4[length(current_last_obs_time4)]
      design_last_obs_time4=design_times_bp[design_times_bp <= current_last_obs_time4+window_width]
      design_last_obs_time4=design_last_obs_time4[length(design_last_obs_time4)]
    }
    if (dim_bp >=5) {
      current_last_obs_time5=temp$obstime[temp$id==eval_id[eval_id_num] & temp$obstime <= s & temp$bptype=="TYPE5"]
      current_last_obs_time5=current_last_obs_time5[length(current_last_obs_time5)]
      design_last_obs_time5=design_times_bp[design_times_bp <= current_last_obs_time5+window_width]
      design_last_obs_time5=design_last_obs_time5[length(design_last_obs_time5)]
    }
#    cat('------------------------------------','\n')
#    cat('Current ID=',eval_id[eval_id_num],'\n')
#    cat('Dataset_nonmiss=','\n')
#    print(dataset_nonmiss[dataset_nonmiss$id==eval_id[eval_id_num],])
#    cat('------------------------------------','\n')
#    cat('current_last_obs_time1=',current_last_obs_time1,'\n')
#    cat('design_last_obs_time1=',design_last_obs_time1,'\n')
#    cat('current_last_obs_time2=',current_last_obs_time2,'\n')
#    cat('design_last_obs_time2=',design_last_obs_time2,'\n')
#    cat('current_last_obs_time3=',current_last_obs_time3,'\n')
#    cat('design_last_obs_time3=',design_last_obs_time3,'\n')

#---------------------------------------------------------------------------------------------------
#  CALCULATE NUMBER OF OBSERVED BP <= S, TIME OF LAST MATCHING OBS, AND VALUE OF LAST MATCHING OBS
#---------------------------------------------------------------------------------------------------
    num_bp1=rep(0,times=n_eval)
    num_bp2=rep(0,times=n_eval)
    num_bp3=rep(0,times=n_eval)
    num_bp4=rep(0,times=n_eval)
    num_bp5=rep(0,times=n_eval)
    match_time1=rep(NA,times=n_eval)
    match_time2=rep(NA,times=n_eval)
    match_time3=rep(NA,times=n_eval)
    match_time4=rep(NA,times=n_eval)
    match_time5=rep(NA,times=n_eval)
    last_match1=rep(NA,times=n_eval)
    last_match2=rep(NA,times=n_eval)
    last_match3=rep(NA,times=n_eval)
    last_match4=rep(NA,times=n_eval)
    last_match5=rep(NA,times=n_eval)
    i=1
    while (i <= n_nonmiss) {
      id=dataset_nonmiss$id[i]
      if (id %in% select_id) {
        id=which(select_id %in% id)
        if (dataset_nonmiss$bptype[i]=="TYPE1") {
          if (dataset_nonmiss$obstime[i] <= s) {
            num_bp1[id]=num_bp1[id]+1
          }
          if (abs(dataset_nonmiss$obstime[i]-design_last_obs_time1) <= window_width) {
            last_match1[id]=dataset_nonmiss$bp_curr[i]
            match_time1[id]=dataset_nonmiss$obstime[i]
          }
        }
        if (dataset_nonmiss$bptype[i]=="TYPE2") {
          if (dataset_nonmiss$obstime[i] <= s) {
            num_bp2[id]=num_bp2[id]+1
          }
          if (abs(dataset_nonmiss$obstime[i]-design_last_obs_time2) <= window_width) {
            last_match2[id]=dataset_nonmiss$bp_curr[i]
            match_time2[id]=dataset_nonmiss$obstime[i]
          }
        }
        if (dataset_nonmiss$bptype[i]=="TYPE3") {
          if (dataset_nonmiss$obstime[i] <= s) {
            num_bp3[id]=num_bp3[id]+1
          }
          if (abs(dataset_nonmiss$obstime[i]-design_last_obs_time3) <= window_width) {
            last_match3[id]=dataset_nonmiss$bp_curr[i]
            match_time3[id]=dataset_nonmiss$obstime[i]
          }
        }
        if (dataset_nonmiss$bptype[i]=="TYPE4") {
          if (dataset_nonmiss$obstime[i] <= s) {
            num_bp4[id]=num_bp4[id]+1
          }
          if (abs(dataset_nonmiss$obstime[i]-design_last_obs_time4) <= window_width) {
            last_match4[id]=dataset_nonmiss$bp_curr[i]
            match_time4[id]=dataset_nonmiss$obstime[i]
          }
        }
        if (dataset_nonmiss$bptype[i]=="TYPE5") {
          if (dataset_nonmiss$obstime[i] <= s) {
            num_bp5[id]=num_bp5[id]+1
          }
          if (abs(dataset_nonmiss$obstime[i]-design_last_obs_time5) <= window_width) {
            last_match5[id]=dataset_nonmiss$bp_curr[i]
            match_time5[id]=dataset_nonmiss$obstime[i]
          }
        }
      }
      i=i+1
    }
#    cat('match_time1=','\n')
#    print(match_time1[1:30])
#    cat('last_match1=','\n')
#    print(last_match1[1:30])
#    cat('match_time2=','\n')
#    print(match_time2[1:30])
#    cat('last_match2=','\n')
#    print(last_match2[1:30])
#    cat('match_time3=','\n')
#    print(match_time3[1:30])
#    cat('last_match3=','\n')
#    print(last_match3[1:30])
#    cat('match_time1[eval_id_num]=',match_time1[eval_id_num],'\n')
#    cat('last_match1=','\n')
#    print(last_match1[1:30])
#    cat('last_match1[eval_id_num]=',last_match1[eval_id_num],'\n')
#    cat('last_match2=','\n')
#    print(last_match2[1:30])
#    cat('last_match2[eval_id_num]=',last_match2[eval_id_num],'\n')
#    cat('last_match3=','\n')
#    print(last_match3[1:30])
#    cat('last_match3[eval_id_num]=',last_match3[eval_id_num],'\n')
#    cat('num_bp1=','\n')
#    print(num_bp1[1:10])
#    cat('--------------------------------------','\n')
#    cat('length(match_time1)=',length(match_time1),'\n')
   
#---------------------------------------------------------------------------------------------------
#  ADVANCE INT & SLOPE VALUES IN DATASET_IMPUTE2
#---------------------------------------------------------------------------------------------------
#    cat('dataset_impute before advanced int1_s for ID=10229=','\n')
#    print(dataset_impute[dataset_impute$id==10229,])

    future_dataset_impute$int1_s[1]=dataset_impute$int1_s[dataset_impute$id==eval_id[eval_id_num]]
    future_dataset_impute$int2_s[1]=dataset_impute$int2_s[dataset_impute$id==eval_id[eval_id_num]]
    future_dataset_impute$int3_s[1]=dataset_impute$int3_s[dataset_impute$id==eval_id[eval_id_num]]
    future_dataset_impute$int4_s[1]=dataset_impute$int4_s[dataset_impute$id==eval_id[eval_id_num]]
    future_dataset_impute$int5_s[1]=dataset_impute$int5_s[dataset_impute$id==eval_id[eval_id_num]]
    future_dataset_impute$slope1_s[1]=dataset_impute$slope1_s[dataset_impute$id==eval_id[eval_id_num]]
    future_dataset_impute$slope2_s[1]=dataset_impute$slope2_s[dataset_impute$id==eval_id[eval_id_num]]
    future_dataset_impute$slope3_s[1]=dataset_impute$slope3_s[dataset_impute$id==eval_id[eval_id_num]]
    future_dataset_impute$slope4_s[1]=dataset_impute$slope4_s[dataset_impute$id==eval_id[eval_id_num]]
    future_dataset_impute$slope5_s[1]=dataset_impute$slope5_s[dataset_impute$id==eval_id[eval_id_num]]

    dataset_impute2=dataset_impute[dataset_impute$id != eval_id[eval_id_num] ,]
    dataset_impute2=rbind(dataset_impute2,future_dataset_impute)
    dataset_impute2=dataset_impute2[order(dataset_impute2$id,dataset_impute2$start),]
    n=length(dataset_impute2$id)

    i=1
    idprev=0
    while (i <= n) {
      id=dataset_impute2$id[i]
      if (id!=idprev) {
        int1_s=dataset_impute2$int1_s[i]
        slope1_s=dataset_impute2$slope1_s[i]
        int2_s=dataset_impute2$int2_s[i]
        slope2_s=dataset_impute2$slope2_s[i]
        int3_s=dataset_impute2$int3_s[i]
        slope3_s=dataset_impute2$slope3_s[i]
        int4_s=dataset_impute2$int4_s[i]
        slope4_s=dataset_impute2$slope4_s[i]
        int5_s=dataset_impute2$int5_s[i]
        slope5_s=dataset_impute2$slope5_s[i]
      }
      if (id==idprev) {
        if (is.na(dataset_impute2$int1_s[i])) {
          dataset_impute2$int1_s[i]=int1_s
          dataset_impute2$slope1_s[i]=slope1_s
        } else {
          int1_s=dataset_impute2$int1_s[i]
          slope1_s=dataset_impute2$slope1_s[i]
        }
        if (is.na(dataset_impute2$int2_s[i])) {
          dataset_impute2$int2_s[i]=int2_s
          dataset_impute2$slope2_s[i]=slope2_s
        } else {
          int2_s=dataset_impute2$int2_s[i]
          slope2_s=dataset_impute2$slope2_s[i]
        }
        if (is.na(dataset_impute2$int3_s[i])) {
          dataset_impute2$int3_s[i]=int3_s
          dataset_impute2$slope3_s[i]=slope3_s
        } else {
          int3_s=dataset_impute2$int3_s[i]
          slope3_s=dataset_impute2$slope3_s[i]
        }
        if (is.na(dataset_impute2$int4_s[i])) {
          dataset_impute2$int4_s[i]=int4_s
          dataset_impute2$slope4_s[i]=slope4_s
        } else {
          int4_s=dataset_impute2$int4_s[i]
          slope4_s=dataset_impute2$slope4_s[i]
        }
        if (is.na(dataset_impute2$int5_s[i])) {
          dataset_impute2$int5_s[i]=int5_s
          dataset_impute2$slope5_s[i]=slope5_s
        } else {
          int5_s=dataset_impute2$int5_s[i]
          slope5_s=dataset_impute2$slope5_s[i]
        }
      }
     idprev=id
      i=i+1
    }
#    cat('dataset_impute2 with advanced int1_s for ID=10229=','\n')
#    print(dataset_impute2[dataset_impute2$id==10229,])
#    cat('# of unique ID in dataset_impute2=',length(unique(dataset_impute2$id)),'\n')

#    cat('dataset_impute2=','\n')
#    print(dataset_impute2[1:20,])
 
#---------------------------------------------------------------------------------------------------
#  FIT PRED MODELS
#---------------------------------------------------------------------------------------------------
    j=s_index+1
    while (j <= num_imp_bp) {
      cat('-----------------------','\n')
      cat('j=',j,'\n')
      cat('Time=',imp_times_bp[j],'\n')

      x=(select_id %in% dataset_impute2$id[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]])
#      cat('length(x)=',length(x),'\n')
     
      pred1=dataset_impute2$bp_curr1[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
#      cat('length(pred1)=',length(pred1),'\n')
#      cat('pred1=','\n')
#      print(pred1[1:30])

      match_time1_s=match_time1[x]
#      cat('length(match_time1_s)=',length(match_time1_s),'\n')
#      cat('match_time1_s=','\n')
#      print(match_time1_s[1:30])

      last_match1_s=last_match1[x]
#      cat('length(last_match1_s)=',length(last_match1_s),'\n')
#      cat('last_match1_s=','\n')
#      print(last_match1_s[1:30])

      num_bp1_s=num_bp1[x]
#      cat('length(num_bp1_s)=',length(num_bp1_s),'\n')
#      cat('num_bp1_s=','\n')
#      print(num_bp1_s[1:10])

      int1_s=dataset_impute2$int1_s[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
#      cat('length(int1_s)=',length(int1_s),'\n')
#      cat('int1_s=','\n')
#      print(int1_s[1:10])

      slope1_s=dataset_impute2$slope1_s[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
#      cat('length(slope1_s)=',length(slope1_s),'\n')
#      cat('slope1_s=','\n')
#      print(slope1_s[1:10])

      base_cov=dataset_impute2[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num],(10+1):(10+ncov)]
#      cat('length(base_cov[,1])=',length(base_cov[,1]),'\n')
#      cat('base_cov[,1]=','\n')
#      print(base_cov[1:10,1])
#      cat('length(base_cov[,2])=',length(base_cov[,2]),'\n')
#      cat('base_cov[,2]=','\n')
#      print(base_cov[1:10,2])
#      cat('length(base_cov[,3])=',length(base_cov[,3]),'\n')
#      cat('base_cov[,3]=','\n')
#      print(base_cov[1:10,3])
#      cat('length(base_cov[,4])=',length(base_cov[,4]),'\n')
#      cat('base_cov[,4]=','\n')
#      print(base_cov[1:10,4])
#      cat('length(base_cov[,5])=',length(base_cov[,5]),'\n')
#      cat('base_cov[,5]=','\n')
#      print(base_cov[1:10,5])
#      cat('length(base_cov[,6])=',length(base_cov[,6]),'\n')
#      cat('base_cov[,6]=','\n')
#      print(base_cov[1:10,6])

      num1_int1_s=num_bp1_s*int1_s
      num1_slope1_s=num_bp1_s*slope1_s
      num1_last1_s=num_bp1_s*last_match1_s
      last1_int1_s=last_match1_s*int1_s
      last1_slope1_s=last_match1_s*slope1_s

      if (dim_bp >=2) {
        pred2=dataset_impute2$bp_curr2[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]

#        cat('pred2=','\n')
#        print(pred2[1:30])

        match_time2_s=match_time2[x]

#        cat('match_time2_s=','\n')
#        print(match_time2_s[1:30])

        last_match2_s=last_match2[x]

#        cat('last_match2_s=','\n')
#        print(last_match2_s[1:30])

        num_bp2_s=num_bp2[x]
        int2_s=dataset_impute2$int2_s[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
        slope2_s=dataset_impute2$slope2_s[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
        num2_int2_s=num_bp2_s*int2_s
        num2_slope2_s=num_bp2_s*slope2_s
        num2_last2_s=num_bp2_s*last_match2_s
        last2_int2_s=last_match2_s*int2_s
        last2_slope2_s=last_match2_s*slope2_s
      }

      if (dim_bp >=3) {
        pred3=dataset_impute2$bp_curr3[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
        match_time3_s=match_time3[x]
        last_match3_s=last_match3[x]
        num_bp3_s=num_bp3[x]
        int3_s=dataset_impute2$int3_s[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
        slope3_s=dataset_impute2$slope3_s[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
        num3_int3_s=num_bp3_s*int3_s
        num3_slope3_s=num_bp3_s*slope3_s
        num3_last3_s=num_bp3_s*last_match3_s
        last3_int3_s=last_match3_s*int3_s
        last3_slope3_s=last_match3_s*slope3_s
      }

      if (dim_bp >=4) {
        pred4=dataset_impute2$bp_curr4[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
        match_time4_s=match_time4[x]
        last_match4_s=last_match4[x]
        num_bp4_s=num_bp4[x]
        int4_s=dataset_impute2$int4_s[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
        slope4_s=dataset_impute2$slope4_s[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
        num4_int4_s=num_bp4_s*int4_s
        num4_slope4_s=num_bp4_s*slope4_s
        num4_last4_s=num_bp4_s*last_match4_s
        last4_int4_s=last_match4_s*int4_s
        last4_slope4_s=last_match4_s*slope4_s
      }

      if (dim_bp >=5) {
        pred5=dataset_impute2$bp_curr5[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
        match_time5_s=match_time5[x]
        last_match5_s=last_match5[x]
        num_bp5_s=num_bp5[x]
        int5_s=dataset_impute2$int5_s[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
        slope5_s=dataset_impute2$slope5_s[dataset_impute2$start==imp_times_bp[j] & dataset_impute2$id != eval_id[eval_id_num]]
        num5_int5_s=num_bp5_s*int5_s
        num5_slope5_s=num_bp5_s*slope5_s
        num5_last5_s=num_bp5_s*last_match5_s
        last5_int5_s=last_match5_s*int5_s
        last5_slope5_s=last_match5_s*slope5_s
      }

      if (dim_bp==1) {

        pfit1=lm(pred1~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                 match_time1_s+num_bp1_s+int1_s+slope1_s+last_match1_s+num1_int1_s+num1_slope1_s+
                 num1_last1_s+last1_int1_s+last1_slope1_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=1','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit1))
        z=coef(pfit1)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time1[eval_id_num]+
          z[1+ncov+2]*num_bp1[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match1[eval_id_num]+
          z[1+ncov+6]*num_bp1[eval_id_num]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp1[eval_id_num]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp1[eval_id_num]*last_match1[eval_id_num]+
          z[1+ncov+9]*last_match1[eval_id_num]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match1[eval_id_num]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

        future_dataset_impute$bp_curr1[future_dataset_impute$start==imp_times_bp[j]]=fitted 

      }

      if (dim_bp==2) {

        pfit1=lm(pred1~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time1_s+num_bp1_s+int1_s+slope1_s+last_match1_s+num1_int1_s+num1_slope1_s+
                   num1_last1_s+last1_int1_s+last1_slope1_s+
                   last_match2_s+int2_s+slope2_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=1','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit1))
        z=coef(pfit1)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time1[eval_id_num]+
          z[1+ncov+2]*num_bp1[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match1[eval_id_num]+
          z[1+ncov+6]*num_bp1[eval_id_num]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp1[eval_id_num]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp1[eval_id_num]*last_match1[eval_id_num]+
          z[1+ncov+9]*last_match1[eval_id_num]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match1[eval_id_num]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match2[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

        future_dataset_impute$bp_curr1[future_dataset_impute$start==imp_times_bp[j]]=fitted 

        pfit2=lm(pred2~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time2_s+num_bp2_s+int2_s+slope2_s+last_match2_s+num2_int2_s+num2_slope2_s+
                   num2_last2_s+last2_int2_s+last2_slope2_s+
                   last_match1_s+int1_s+slope1_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=2','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit2))
        z=coef(pfit2)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time2[eval_id_num]+
          z[1+ncov+2]*num_bp2[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match2[eval_id_num]+
          z[1+ncov+6]*num_bp2[eval_id_num]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp2[eval_id_num]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp2[eval_id_num]*last_match2[eval_id_num]+
          z[1+ncov+9]*last_match2[eval_id_num]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match2[eval_id_num]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match1[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

        future_dataset_impute$bp_curr2[future_dataset_impute$start==imp_times_bp[j]]=fitted 

      }

      if (dim_bp==3) {

        pfit1=lm(pred1~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time1_s+num_bp1_s+int1_s+slope1_s+last_match1_s+num1_int1_s+num1_slope1_s+
                   num1_last1_s+last1_int1_s+last1_slope1_s+
                   last_match2_s+int2_s+slope2_s+
                   last_match3_s+int3_s+slope3_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=1','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit1))
        z=coef(pfit1)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
#        cat('after baseline covs fitted=',fitted,'\n')
        fitted=fitted+z[1+ncov+1]*match_time1[eval_id_num]+
          z[1+ncov+2]*num_bp1[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match1[eval_id_num]+
          z[1+ncov+6]*num_bp1[eval_id_num]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp1[eval_id_num]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp1[eval_id_num]*last_match1[eval_id_num]+
          z[1+ncov+9]*last_match1[eval_id_num]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match1[eval_id_num]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match2[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+14]*last_match3[eval_id_num]+
          z[1+ncov+15]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+16]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

#        cat('fitted=',fitted,'\n')
#        cat('dataset_impute2$int1_s[ID=eval_id & start=imp_times_bp]=',dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]],'\n')
#        cat('dataset_impute2$slope1_s[ID=eval_id & start=imp_times_bp]=',dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]],'\n')
#        cat('num_bp1[eval_id_num]=',num_bp1[eval_id_num],'\n')
#        cat('match_time1[eval_id_num]=',match_time1[eval_id_num],'\n')
#        cat('last_match1[eval_id_num]=',last_match1[eval_id_num],'\n')
#        cat('dataset_impute2$int2_s[ID=eval_id & start=imp_times_bp]=',dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]],'\n')
#        cat('dataset_impute2$slope2_s[ID=eval_id & start=imp_times_bp]=',dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]],'\n')
#        cat('last_match2[eval_id_num]=',last_match2[eval_id_num],'\n')
#        cat('dataset_impute2$int3_s[ID=eval_id & start=imp_times_bp]=',dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]],'\n')
#        cat('dataset_impute2$slope3_s[ID=eval_id & start=imp_times_bp]=',dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]],'\n')
#        cat('last_match3[eval_id_num]=',last_match3[eval_id_num],'\n')
#        cat('future_dataset_impute=','\n')
#        print(future_dataset_impute)

        future_dataset_impute$bp_curr1[future_dataset_impute$start==imp_times_bp[j]]=fitted 

        pfit2=lm(pred2~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time2_s+num_bp2_s+int2_s+slope2_s+last_match2_s+num2_int2_s+num2_slope2_s+
                   num2_last2_s+last2_int2_s+last2_slope2_s+
                   last_match1_s+int1_s+slope1_s+
                   last_match3_s+int3_s+slope3_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=2','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit2))
        z=coef(pfit2)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time2[eval_id_num]+
          z[1+ncov+2]*num_bp2[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match2[eval_id_num]+
          z[1+ncov+6]*num_bp2[eval_id_num]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp2[eval_id_num]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp2[eval_id_num]*last_match2[eval_id_num]+
          z[1+ncov+9]*last_match2[eval_id_num]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match2[eval_id_num]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match1[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+14]*last_match3[eval_id_num]+
          z[1+ncov+15]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+16]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

#        cat('fitted=',fitted,'\n')
        future_dataset_impute$bp_curr2[future_dataset_impute$start==imp_times_bp[j]]=fitted 

        pfit3=lm(pred3~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time3_s+num_bp3_s+int3_s+slope3_s+last_match3_s+num3_int3_s+num3_slope3_s+
                   num3_last3_s+last3_int3_s+last3_slope3_s+
                   last_match1_s+int1_s+slope1_s+
                   last_match2_s+int2_s+slope2_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=3','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit3))
        z=coef(pfit3)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time3[eval_id_num]+
          z[1+ncov+2]*num_bp3[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match3[eval_id_num]+
          z[1+ncov+6]*num_bp3[eval_id_num]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp3[eval_id_num]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp3[eval_id_num]*last_match3[eval_id_num]+
          z[1+ncov+9]*last_match3[eval_id_num]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match3[eval_id_num]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match1[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+14]*last_match2[eval_id_num]+
          z[1+ncov+15]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+16]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

#        cat('fitted=',fitted,'\n')
        future_dataset_impute$bp_curr3[future_dataset_impute$start==imp_times_bp[j]]=fitted 

#        cat('future_dataset_impute=','\n')
#        print(future_dataset_impute)

      }

      if (dim_bp==4) {

        pfit1=lm(pred1~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time1_s+num_bp1_s+int1_s+slope1_s+last_match1_s+num1_int1_s+num1_slope1_s+
                   num1_last1_s+last1_int1_s+last1_slope1_s+
                   last_match2_s+int2_s+slope2_s+
                   last_match3_s+int3_s+slope3_s+
                   last_match4_s+int4_s+slope4_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=1','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit1))
        z=coef(pfit1)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time1[eval_id_num]+
          z[1+ncov+2]*num_bp1[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match1[eval_id_num]+
          z[1+ncov+6]*num_bp1[eval_id_num]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp1[eval_id_num]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp1[eval_id_num]*last_match1[eval_id_num]+
          z[1+ncov+9]*last_match1[eval_id_num]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match1[eval_id_num]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match2[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+14]*last_match3[eval_id_num]+
          z[1+ncov+15]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+16]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+17]*last_match4[eval_id_num]+
          z[1+ncov+18]*dataset_impute2$int4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+19]*dataset_impute2$slope4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

        future_dataset_impute$bp_curr1[future_dataset_impute$start==imp_times_bp[j]]=fitted 

        pfit2=lm(pred2~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time2_s+num_bp2_s+int2_s+slope2_s+last_match2_s+num2_int2_s+num2_slope2_s+
                   num2_last2_s+last2_int2_s+last2_slope2_s+
                   last_match1_s+int1_s+slope1_s+
                   last_match3_s+int3_s+slope3_s+
                   last_match4_s+int4_s+slope4_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=2','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit2))
        z=coef(pfit2)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time2[eval_id_num]+
          z[1+ncov+2]*num_bp2[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match2[eval_id_num]+
          z[1+ncov+6]*num_bp2[eval_id_num]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp2[eval_id_num]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp2[eval_id_num]*last_match2[eval_id_num]+
          z[1+ncov+9]*last_match2[eval_id_num]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match2[eval_id_num]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match1[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+14]*last_match3[eval_id_num]+
          z[1+ncov+15]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+16]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+17]*last_match4[eval_id_num]+
          z[1+ncov+18]*dataset_impute2$int4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+19]*dataset_impute2$slope4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

        future_dataset_impute$bp_curr2[future_dataset_impute$start==imp_times_bp[j]]=fitted 

        pfit3=lm(pred3~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time3_s+num_bp3_s+int3_s+slope3_s+last_match3_s+num3_int3_s+num3_slope3_s+
                   num3_last3_s+last3_int3_s+last3_slope3_s+
                   last_match1_s+int1_s+slope1_s+
                   last_match2_s+int2_s+slope2_s+
                   last_match4_s+int4_s+slope4_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=3','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit3))
        z=coef(pfit3)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time3[eval_id_num]+
          z[1+ncov+2]*num_bp3[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match3[eval_id_num]+
          z[1+ncov+6]*num_bp3[eval_id_num]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp3[eval_id_num]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp3[eval_id_num]*last_match3[eval_id_num]+
          z[1+ncov+9]*last_match3[eval_id_num]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match3[eval_id_num]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match1[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+14]*last_match2[eval_id_num]+
          z[1+ncov+15]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+16]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+17]*last_match4[eval_id_num]+
          z[1+ncov+18]*dataset_impute2$int4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+19]*dataset_impute2$slope4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

        future_dataset_impute$bp_curr3[future_dataset_impute$start==imp_times_bp[j]]=fitted 

        pfit4=lm(pred4~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time4_s+num_bp4_s+int4_s+slope4_s+last_match4_s+num4_int4_s+num4_slope4_s+
                   num4_last4_s+last4_int4_s+last4_slope4_s+
                   last_match1_s+int1_s+slope1_s+
                   last_match2_s+int2_s+slope2_s+
                   last_match3_s+int3_s+slope3_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=4','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit4))
        z=coef(pfit4)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time4[eval_id_num]+
          z[1+ncov+2]*num_bp4[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match4[eval_id_num]+
          z[1+ncov+6]*num_bp4[eval_id_num]*dataset_impute2$int4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp4[eval_id_num]*dataset_impute2$slope4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp4[eval_id_num]*last_match4[eval_id_num]+
          z[1+ncov+9]*last_match4[eval_id_num]*dataset_impute2$int4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match4[eval_id_num]*dataset_impute2$slope4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match1[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+14]*last_match2[eval_id_num]+
          z[1+ncov+15]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+16]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+17]*last_match3[eval_id_num]+
          z[1+ncov+18]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+19]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

        future_dataset_impute$bp_curr4[future_dataset_impute$start==imp_times_bp[j]]=fitted 

      }

      if (dim_bp==5) {

        pfit1=lm(pred1~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time1_s+num_bp1_s+int1_s+slope1_s+last_match1_s+num1_int1_s+num1_slope1_s+
                   num1_last1_s+last1_int1_s+last1_slope1_s+
                   last_match2_s+int2_s+slope2_s+
                   last_match3_s+int3_s+slope3_s+
                   last_match4_s+int4_s+slope4_s+
                   last_match5_s+int5_s+slope5_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=1','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit1))
        z=coef(pfit1)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time1[eval_id_num]+
          z[1+ncov+2]*num_bp1[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match1[eval_id_num]+
          z[1+ncov+6]*num_bp1[eval_id_num]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp1[eval_id_num]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp1[eval_id_num]*last_match1[eval_id_num]+
          z[1+ncov+9]*last_match1[eval_id_num]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match1[eval_id_num]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match2[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+14]*last_match3[eval_id_num]+
          z[1+ncov+15]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+16]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+17]*last_match4[eval_id_num]+
          z[1+ncov+18]*dataset_impute2$int4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+19]*dataset_impute2$slope4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+20]*last_match5[eval_id_num]+
          z[1+ncov+21]*dataset_impute2$int5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+22]*dataset_impute2$slope5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

        future_dataset_impute$bp_curr1[future_dataset_impute$start==imp_times_bp[j]]=fitted 

        pfit2=lm(pred2~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time2_s+num_bp2_s+int2_s+slope2_s+last_match2_s+num2_int2_s+num2_slope2_s+
                   num2_last2_s+last2_int2_s+last2_slope2_s+
                   last_match1_s+int1_s+slope1_s+
                   last_match3_s+int3_s+slope3_s+
                   last_match4_s+int4_s+slope4_s+
                   last_match5_s+int5_s+slope5_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=2','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit2))
        z=coef(pfit2)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time2[eval_id_num]+
          z[1+ncov+2]*num_bp2[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match2[eval_id_num]+
          z[1+ncov+6]*num_bp2[eval_id_num]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp2[eval_id_num]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp2[eval_id_num]*last_match2[eval_id_num]+
          z[1+ncov+9]*last_match2[eval_id_num]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match2[eval_id_num]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match1[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+14]*last_match3[eval_id_num]+
          z[1+ncov+15]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+16]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+17]*last_match4[eval_id_num]+
          z[1+ncov+18]*dataset_impute2$int4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+19]*dataset_impute2$slope4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+20]*last_match5[eval_id_num]+
          z[1+ncov+21]*dataset_impute2$int5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+22]*dataset_impute2$slope5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

        future_dataset_impute$bp_curr2[future_dataset_impute$start==imp_times_bp[j]]=fitted 

        pfit3=lm(pred3~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time3_s+num_bp3_s+int3_s+slope3_s+last_match3_s+num3_int3_s+num3_slope3_s+
                   num3_last3_s+last3_int3_s+last3_slope3_s+
                   last_match1_s+int1_s+slope1_s+
                   last_match2_s+int2_s+slope2_s+
                   last_match4_s+int4_s+slope4_s+
                   last_match5_s+int5_s+slope5_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=3','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit3))
        z=coef(pfit3)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time3[eval_id_num]+
          z[1+ncov+2]*num_bp3[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match3[eval_id_num]+
          z[1+ncov+6]*num_bp3[eval_id_num]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp3[eval_id_num]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp3[eval_id_num]*last_match3[eval_id_num]+
          z[1+ncov+9]*last_match3[eval_id_num]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match3[eval_id_num]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match1[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+14]*last_match2[eval_id_num]+
          z[1+ncov+15]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+16]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+17]*last_match4[eval_id_num]+
          z[1+ncov+18]*dataset_impute2$int4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+19]*dataset_impute2$slope4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+20]*last_match5[eval_id_num]+
          z[1+ncov+21]*dataset_impute2$int5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+22]*dataset_impute2$slope5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

        future_dataset_impute$bp_curr3[future_dataset_impute$start==imp_times_bp[j]]=fitted 

        pfit4=lm(pred4~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time4_s+num_bp4_s+int4_s+slope4_s+last_match4_s+num4_int4_s+num4_slope4_s+
                   num4_last4_s+last4_int4_s+last4_slope4_s+
                   last_match1_s+int1_s+slope1_s+
                   last_match2_s+int2_s+slope2_s+
                   last_match3_s+int3_s+slope3_s+
                   last_match5_s+int5_s+slope5_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=4','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit4))
        z=coef(pfit4)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time4[eval_id_num]+
          z[1+ncov+2]*num_bp4[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match4[eval_id_num]+
          z[1+ncov+6]*num_bp4[eval_id_num]*dataset_impute2$int4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp4[eval_id_num]*dataset_impute2$slope4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp4[eval_id_num]*last_match4[eval_id_num]+
          z[1+ncov+9]*last_match4[eval_id_num]*dataset_impute2$int4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match4[eval_id_num]*dataset_impute2$slope4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match1[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+14]*last_match2[eval_id_num]+
          z[1+ncov+15]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+16]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+17]*last_match3[eval_id_num]+
          z[1+ncov+18]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+19]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+20]*last_match5[eval_id_num]+
          z[1+ncov+21]*dataset_impute2$int5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+22]*dataset_impute2$slope5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

        future_dataset_impute$bp_curr4[future_dataset_impute$start==imp_times_bp[j]]=fitted 

        pfit5=lm(pred5~base_cov[,1]+base_cov[,2]+base_cov[,3]+base_cov[,4]+
                   match_time5_s+num_bp5_s+int5_s+slope5_s+last_match5_s+num5_int5_s+num5_slope5_s+
                   num5_last5_s+last5_int5_s+last5_slope5_s+
                   last_match1_s+int1_s+slope1_s+
                   last_match2_s+int2_s+slope2_s+
                   last_match3_s+int3_s+slope3_s+
                   last_match4_s+int4_s+slope4_s)
       
        cat('-------------------------------------------------------','\n')
        cat('-------------------------------------------------------','\n')
        cat('PRED Model for Time t=',imp_times_bp[j],'\n')
        cat('BP DIM=5','\n')
        cat('-------------------------------------------------------','\n')
#        print(summary(pfit5))
        z=coef(pfit5)
        z[is.na(z)]=0
        fitted=z[1]
        jj=1
        while (jj <= ncov) {
          fitted=fitted+z[1+jj]*dataset_impute2[1,10+jj]
          jj=jj+1
        }
        fitted=fitted+z[1+ncov+1]*match_time5[eval_id_num]+
          z[1+ncov+2]*num_bp5[eval_id_num]+
          z[1+ncov+3]*dataset_impute2$int5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+4]*dataset_impute2$slope5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+5]*last_match5[eval_id_num]+
          z[1+ncov+6]*num_bp5[eval_id_num]*dataset_impute2$int5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+7]*num_bp5[eval_id_num]*dataset_impute2$slope5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+8]*num_bp5[eval_id_num]*last_match5[eval_id_num]+
          z[1+ncov+9]*last_match5[eval_id_num]*dataset_impute2$int5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+10]*last_match5[eval_id_num]*dataset_impute2$slope5_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+11]*last_match1[eval_id_num]+
          z[1+ncov+12]*dataset_impute2$int1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+13]*dataset_impute2$slope1_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+14]*last_match2[eval_id_num]+
          z[1+ncov+15]*dataset_impute2$int2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+16]*dataset_impute2$slope2_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+17]*last_match3[eval_id_num]+
          z[1+ncov+18]*dataset_impute2$int3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+19]*dataset_impute2$slope3_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+20]*last_match4[eval_id_num]+
          z[1+ncov+21]*dataset_impute2$int4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]+
          z[1+ncov+22]*dataset_impute2$slope4_s[dataset_impute2$id==eval_id[eval_id_num] & dataset_impute2$start==imp_times_bp[j]]

        future_dataset_impute$bp_curr5[future_dataset_impute$start==imp_times_bp[j]]=fitted 

      }
      cat('-------------------------------------------------------','\n')
      cat('END: PRED Models at Time=',imp_times_bp[j],'\n')
      cat('-------------------------------------------------------','\n')
 
      j=j+1
    }

    cat('Future Dataset with Perfect Match imputed future BP values=','\n')
    print(future_dataset_impute)
#--------------------------------------------------------------------------
#
# FORTRAN VERSION: CALCULATE CONDITIONAL SURVIVAL PROBS
#
#--------------------------------------------------------------------------
    cond_surv_pred_st=rep(0,times=1)
    cat('Call predictSurvProb','\n')
    future_surv=predictSurvProb(efit,newdata=future_dataset_impute,times=predtimes)
    cat('Return from predictSurvProb','\n')
#  if (!all(is.finite(future_surv))) {
#    cat('survival to predtimes for subjects in future_dataset according to EMPIRICAL COX MODEL=','\n')
#    print(future_surv[1:20,])
#  }
#  future_surv[future_surv > 1]=NA
#  future_surv[future_surv < 0]=NA
    npredtimes=length(predtimes)
    future_surv_short=future_surv
#  future_surv_short=future_surv[1:n_future,]
#  if (all(is.finite(future_surv_short))) {
#    cat('All future_surv_short are finite','\n')
#  } else {
#    cat('survival to predtimes for subjects in future_dataset according to EMPIRICAL COX MODEL=','\n')
#    print(future_surv_short)
#  }


#test=future_dataset_impute_long
#test$bp.1[1]=999
#test_surv=predictSurvProb(efit[[k]],newdata=test,times=predtimes)
#test_surv_short=test_surv[1:n_future,]

#cat('survival to predtimes for subjects in test dataset according to EMPIRICAL COX MODEL=','\n')
#print(test_surv_short)


#  idvec=future_dataset_impute_long$id
    idvec=future_dataset_impute$id

#  cat('survival to predtimes for subjects in future_dataset according to EMPIRICAL COX MODEL=','\n')
#  print(future_surv_short)

    if (all(is.finite(future_surv_short))) {
      cat('Call condprob','\n')
      z=condprob(1,n_future,npredtimes,index_s,future_surv_short,idvec,cond_surv_pred_st)
      cat('return from condprob','\n')

      cond_surv_pred_st=z[[7]]

      cat('FORTRAN: Cond. survival to S+TAU given survival to S for subjects in future_dataset according to EB2 COX MODEL=','\n')
      print(cond_surv_pred_st)

      icox2_cond_surv_pred_imp=mean(cond_surv_pred_st)

      cat('FORTRAN: Avg. Cond. survival to S+TAU given survival to S for subjects in future_dataset according to EB2 COX MODEL=','\n')
      print(icox2_cond_surv_pred_imp)
      } else {
      cat('EB2 COX MODEL DOES NOT GIVE PROB ESTIMATE','\n')
      icox2_cond_surv_pred_imp=NA
    }

    cat('---------------------------------------------------','\n')
    cat('---------------------------------------------------','\n')
    if (is.finite(icox2_cond_surv_pred_imp)) {
     icox2_cond_surv_pred[eval_id_num]=mean(icox2_cond_surv_pred_imp,na.rm=TRUE)
      cat('FINAL Estimate Cond. survival to S+TAU given survival to S for future subj according to EB2 COX MODEL=','\n')
      print(icox2_cond_surv_pred[eval_id_num])
    } else {
      cat('FINAL EB2 COX MODEL DOES NOT GIVE PROB ESTIMATE','\n')
      icox2_cond_surv_pred[eval_id_num]=NA
    }


#  eval_id_num=eval_id_num+1
#}


  }
#-----------------------------------------------------
# END ICOX2
#-----------------------------------------------------










#-------------------------------------------------
#-------------------------------------------------
if (lvcf_ind==1) {
#-------------------------------------------------
# LAST VALUE CARRIED FORWARD METHOD
#-------------------------------------------------
#-------------------------------------------------
#--------------------------------------------------------------------------
cat('-------------------------------------------','\n')
cat('START LVCF','\n')
cat('-------------------------------------------','\n')

dataset_lvcf=dataset_eval_save[dataset_eval_save$lastobs==1,]
dataset_lvcf$start=rep(0,times=length(dataset_lvcf$id))
dataset_bio_lvcf=dataset_nonmiss[order(dataset_nonmiss$id,dataset_nonmiss$obstime),]
n_bio=length(dataset_bio_lvcf$id)

#cat('LVCF Dataset=','\n')
#print (dataset_lvcf[1:20,])

#cat('Biomarker Dataset for LVCF=','\n')
#print (dataset_bio_lvcf[1:20,])

#--------------------------------------------------------------------------
#
# Put bp_curr1, bp_curr2, bp_curr3, bp_curr4, bp_curr5 into dataset_lvcf
#
#--------------------------------------------------------------------------
idvec=rep(0,times=n_bio)
startvec=rep(0,times=n_bio)
stopvec=rep(0,times=n_bio)
statusvec=rep(0,times=n_bio)
basevec=rep(0,times=n_bio*ncov)
dim(basevec)=c(n_bio,ncov)
bp_curr1vec=rep(0,times=n_bio)
bp_curr2vec=rep(0,times=n_bio)
bp_curr3vec=rep(0,times=n_bio)
bp_curr4vec=rep(0,times=n_bio)
bp_curr5vec=rep(0,times=n_bio)
n=length(dataset_lvcf$id)
newn=0
inonmiss=1
breakout=0
i=1
while (i <= n) {
#  cat('i=',i,'\n')
#  cat('inonmiss=',inonmiss,'\n')
  newn=newn+1
#  cat('newn=',newn,'\n')
  idvec[newn]=dataset_lvcf$id[i]
  startvec[newn]=dataset_lvcf$start[i]
  stopvec[newn]=dataset_lvcf$stop[i]
  statusvec[newn]=dataset_lvcf$status[i]
  j=1
  while (j <= ncov) {
    temp=eval(parse(text=paste('dataset_lvcf$base_cov.',j,'[i]',sep='')))
    basevec[newn,j]=temp
    j=j+1
  }
  if (newn > 1) {
    bp_curr1vec[newn]=bp_curr1vec[newn-1]
    bp_curr2vec[newn]=bp_curr2vec[newn-1]
    bp_curr3vec[newn]=bp_curr3vec[newn-1]
    bp_curr4vec[newn]=bp_curr4vec[newn-1]
    bp_curr5vec[newn]=bp_curr5vec[newn-1]
  }
  if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE1') {bp_curr1vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
  if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE2') {bp_curr2vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
  if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE3') {bp_curr3vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
  if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE4') {bp_curr4vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
  if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE5') {bp_curr5vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
  if (inonmiss < n_bio) {
    inonmiss=inonmiss+1
  } else {
    breakout=1
  }
#  cat('inonmiss=',inonmiss,'\n')
  while (dataset_bio_lvcf$id[inonmiss]==dataset_lvcf$id[i]) {
    if (breakout==1) {break}
    if (dataset_bio_lvcf$obstime[inonmiss] < dataset_lvcf$stop[i]) {
      stopvec[newn]=dataset_bio_lvcf$obstime[inonmiss]
      statusvec[newn]=0
      newn=newn+1
      idvec[newn]=dataset_bio_lvcf$id[inonmiss]
      startvec[newn]=dataset_bio_lvcf$obstime[inonmiss]
      stopvec[newn]=dataset_lvcf$stop[i]
      statusvec[newn]=dataset_lvcf$status[i]
      j=1
      while (j <= ncov) {
        temp=eval(parse(text=paste('dataset_lvcf$base_cov.',j,'[i]',sep='')))
        basevec[newn,j]=temp
        j=j+1
      }
      if (newn > 1) {
        bp_curr1vec[newn]=bp_curr1vec[newn-1]
        bp_curr2vec[newn]=bp_curr2vec[newn-1]
        bp_curr3vec[newn]=bp_curr3vec[newn-1]
        bp_curr4vec[newn]=bp_curr4vec[newn-1]
        bp_curr5vec[newn]=bp_curr5vec[newn-1]
      }
      if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE1') {bp_curr1vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
      if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE2') {bp_curr2vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
      if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE3') {bp_curr3vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
      if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE4') {bp_curr4vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
      if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE5') {bp_curr5vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
      continue=0
      if (inonmiss < n_bio){continue=(dataset_bio_lvcf$obstime[inonmiss+1] < dataset_bio_lvcf$obstime[inonmiss]+.000001)}
      while (continue==1) {
        inonmiss=inonmiss+1
        if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE1') {bp_curr1vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
        if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE2') {bp_curr2vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
        if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE3') {bp_curr3vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
        if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE4') {bp_curr4vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
        if (dataset_bio_lvcf$bptype[inonmiss]=='TYPE5') {bp_curr5vec[newn]=dataset_bio_lvcf$bp_curr[inonmiss]}
        continue=0
        if (inonmiss < n_bio){continue=(dataset_bio_lvcf$obstime[inonmiss+1] < dataset_bio_lvcf$obstime[inonmiss]+.000001)}
      }
    }
    inonmiss=inonmiss+1
    if (inonmiss > n_bio) {break}
  }
  i=i+1
}
dataset_lvcf=data.frame(id=idvec,start=startvec,stop=stopvec,status=statusvec,
     bp_curr1=bp_curr1vec,bp_curr2=bp_curr2vec,
     bp_curr3=bp_curr3vec,bp_curr4=bp_curr4vec,bp_curr5=bp_curr5vec,
     base_cov=basevec)
#cat('After adding various bp_curr variables LVCF Dataset=','\n')
#print (dataset_lvcf[1:50,])

dataset_lvcf=dataset_lvcf[1:newn,]
dataset_cox=dataset_lvcf[dataset_lvcf$start < dataset_lvcf$stop,]

#
# Fit LVCF Cox Model
#
#--------------------------------------------------------------------------
cat('------------------------------------------------------------','\n')
cat('LVCF Cox model uses apriori FIXED TERMS','\n')
cat('------------------------------------------------------------','\n')
if (dim_bp==1) {
  (ecformula=as.formula(paste("Surv(start,stop,status)~bp_curr1",base_names,sep="")))
}
if (dim_bp==2) {
  (ecformula=as.formula(paste("Surv(start,stop,status)~bp_curr1+bp_curr2",base_names,sep="")))
}
if (dim_bp==3) {
  (ecformula=as.formula(paste("Surv(start,stop,status)~bp_curr1+bp_curr2+bp_curr3",base_names,sep="")))
}
if (dim_bp==4) {
  (ecformula=as.formula(paste("Surv(start,stop,status)~bp_curr1+bp_curr2+bp_curr3+bp_curr4",base_names,sep="")))
}
if (dim_bp==5) {
  (ecformula=as.formula(paste("Surv(start,stop,status)~bp_curr1+bp_curr2+bp_curr3+bp_curr4+bp_curr5",base_names,sep="")))
}
lfit=try(coxph(ecformula,data=dataset_cox,x=TRUE))
cat('------------------------------------------------------------','\n')
cat('LVCF Cox model with ',dim_bp,' biomarker terms is','\n')
print(lfit)
cat('-----------------------------------------------------------------','\n')

#--------------------------------------------------------------------------
#
# Use COX MODEL to estimate survival for subjects in FUTURE_DATASET
#
#--------------------------------------------------------------------------

#cat('eval_id[eval_id_num]=',eval_id[eval_id_num],'\n')
future_dataset_lvcf=dataset_lvcf[dataset_lvcf$id==eval_id[eval_id_num],]

#cat('Before restriction, Future LVCF Dataset=','\n')
#print(future_dataset_lvcf)

future_dataset_lvcf=future_dataset_lvcf[future_dataset_lvcf$start <= s,]
n_future_lvcf=length(future_dataset_lvcf$id)
if (!any(future_dataset_lvcf$stop==s)) {
  temp=future_dataset_lvcf[1,]
  temp$stop[1]=s
  future_dataset_lvcf=rbind(future_dataset_lvcf,temp)
  n_future_lvcf=length(future_dataset_lvcf$id)
#  cat('n_future_lvcf=',n_future_lvcf,'\n')
  future_dataset_lvcf=future_dataset_lvcf[order(future_dataset_lvcf$stop),]
  temp_num=which(future_dataset_lvcf$stop==s)
#  cat('temp_num=',temp_num,'\n')
  if (temp_num>1) {future_dataset_lvcf$start[temp_num]=future_dataset_lvcf$stop[temp_num-1]}
  if (temp_num>1) {future_dataset_lvcf$bp_curr1[temp_num]=future_dataset_lvcf$bp_curr1[temp_num-1]}
  if (temp_num>1) {future_dataset_lvcf$bp_curr2[temp_num]=future_dataset_lvcf$bp_curr2[temp_num-1]}
  if (temp_num>1) {future_dataset_lvcf$bp_curr3[temp_num]=future_dataset_lvcf$bp_curr3[temp_num-1]}
  if (temp_num>1) {future_dataset_lvcf$bp_curr4[temp_num]=future_dataset_lvcf$bp_curr4[temp_num-1]}
  if (temp_num>1) {future_dataset_lvcf$bp_curr5[temp_num]=future_dataset_lvcf$bp_curr5[temp_num-1]}
  if (temp_num < n_future_lvcf) {future_dataset_lvcf$start[n_future_lvcf]=s}
}
future_dataset_lvcf$stop[n_future_lvcf]=s+tau
#cat('After restriction of info, Future LVCF Dataset=','\n')
#print(future_dataset_lvcf)

future_dataset_lvcf=future_dataset_lvcf[future_dataset_lvcf$start < future_dataset_lvcf$stop,]
n_future_lvcf=length(future_dataset_lvcf$id)

#cat('After restriction  for start < stop, Future LVCF Dataset=','\n')
#print(future_dataset_lvcf)

predtimes_lvcf=future_dataset_lvcf$stop
npredtimes=length(predtimes_lvcf)
index_s_lvcf=which(predtimes_lvcf==s)

cat('Predtimes=','\n')
print(predtimes_lvcf)

#--------------------------------------------------------------------------
#
# R VERSION OF LOOP TO CALCULATE CONDITIONAL SURVIVAL PROBS
#
#--------------------------------------------------------------------------
#surv_pred_s=rep(0,times=nperson_future)
#surv_pred_st=rep(0,times=nperson_future)
#cond_surv_pred_st=rep(0,times=nperson_future)
#future_surv=predictSurvProb(lfit,newdata=future_dataset_lvcf_long,times=predtimes)
#id=1
#while (id <= nperson_future) {
#  future_surv_id=future_surv[future_dataset_lvcf$id==id,]
#  cat('ID=',id,'\n')
#  cat('Survival for subject in future_dataset according to EMPIRICAL COX MODEL=','\n')
#  print(future_surv_id)
#  surv_pred_st[id]=future_surv_id[1,1]
#  surv_pred_s[id]=future_surv_id[1,1]
#  i=2
#  while (i <= length(future_surv_id[,1])) {
#    if (future_surv_id[i,i-1] != 0) {
#      surv_pred_st[id]=surv_pred_st[id]*future_surv_id[i,i]/future_surv_id[i,i-1]
#    }
#    if (i <= index_s) {
#      if (future_surv_id[i,i-1] != 0) {
#        surv_pred_s[id]=surv_pred_s[id]*future_surv_id[i,i]/future_surv_id[i,i-1]
#      }
#    }
#    i=i+1
#  }
#  if (surv_pred_s[id] != 0) {
#    cond_surv_pred_st[id]=surv_pred_st[id]/surv_pred_s[id]
#  } else {
#    cond_surv_pred_st[id]=1
#  }
#  id=id+1
#}
#lvcf_cond_surv_pred[irep]=mean(cond_surv_pred_st)
#cat('R: Avg. Cond. survival to S+TAU given survival to S for subjects in future_dataset according to COX MODEL=','\n')
#print(lvcf_cond_surv_pred[irep])
#cat('Survival to S for subjects in future_dataset according to COX MODEL=','\n')
#print(surv_pred_s)
#cat('Survival to S+TAU for subjects in future_dataset according to COX MODEL=','\n')
#print(surv_pred_st)
#cat('Cond. survival to S+TAU given survival to S for subjects in future_dataset according to COX MODEL=','\n')
#print(cond_surv_pred_st)
#--------------------------------------------------------------------------
#
# FORTRAN VERSION OF LOOP TO CALCULATE CONDITIONAL SURVIVAL PROBS
#
#--------------------------------------------------------------------------
cond_surv_pred_st=rep(0,times=1)
cat('Call predictSurvProb','\n')
#future_surv=predictSurvProb(lfit,newdata=future_dataset_lvcf_long,times=predtimes)
future_surv=predictSurvProb(lfit,newdata=future_dataset_lvcf,times=predtimes_lvcf)
cat('Return from predictSurvProb','\n')
#print(future_surv[1:20,1:10])
#if (!all(is.finite(future_surv))) {
#  cat('survival to predtimes for subjects in future_dataset according to LVCF COX MODEL=','\n')
#  print(future_surv[1:20,1:10])
#}
#future_surv[future_surv > 1]=NA
#future_surv[future_surv < 0]=NA
#cat('predtimes=','\n')
#print(predtimes)

npredtimes=length(predtimes_lvcf)
#cat('npredtimes=',npredtimes,'\n')
#cat('n_future=',n_future,'\n')
#cat('dim(future_surv)=','\n')
#print(dim(future_surv))

#future_surv_short=future_surv[1:n_future,]
future_surv_short=future_surv

#cat('survival to predtimes for subjects in future_dataset according to LVCF COX MODEL=','\n')
#print (future_surv_short)

#test=future_dataset_lvcf_long
#test$bp.1[1]=999
#test_surv=predictSurvProb(lfit,newdata=test,times=predtimes)
#test_surv_short=test_surv[1:n_future,]

#cat('survival to predtimes for subjects in test dataset according to COX MODEL=','\n')
#print(test_surv_short)


#idvec=future_dataset_lvcf_long$id
idvec=future_dataset_lvcf$id

if (all(is.finite(future_surv_short))) {

  cat('Call condprob','\n')
  z=condprob(1,n_future_lvcf,npredtimes,index_s_lvcf,future_surv_short,idvec,cond_surv_pred_st)
  cat('Return from condprob','\n')

  cond_surv_pred_st=z[[7]]

#  cat('FORTRAN: Cond. survival to S+TAU given survival to S for subjects in future_dataset according to LVCF COX MODEL=','\n')
#  print(cond_surv_pred_st)

  lvcf_cond_surv_pred[eval_id_num]=mean(cond_surv_pred_st)

  cat('FORTRAN: Avg. Cond. survival to S+TAU given survival to S for subjects in future_dataset according to LVCF COX MODEL=','\n')
  print(lvcf_cond_surv_pred[eval_id_num])
} else {
  cat('LVCF COX MODEL DOES NOT GIVE PROB ESTIMATE','\n')
  lvcf_cond_surv_pred[eval_id_num]=NA
}


#-------------------------------------------------
#-------------------------------------------------
}
#-------------------------------------------------
# END: LAST VALUE CARRIED FORWARD
# END: lvcf_ind=1
#-------------------------------------------------



  eval_id_num=eval_id_num+1
}


}
#-----------------------------------------------------
# END ICOX1 OR ICOX2
#-----------------------------------------------------














#-------------------------------------------------
# START: BOOTSTRAP
#-------------------------------------------------
icox1_cond_surv_pred_boot=rep(0,times=boot_num)
icox2_cond_surv_pred_boot=rep(0,times=boot_num)

#cat('dataset =','\n')
#print(dataset_save[1:30,])
#cat('--------------------------------------','\n')
#cat('dataset nonmiss =','\n')
#print(dataset_nonmiss_save[1:30,])
#cat('--------------------------------------','\n')

iboot=1
while(iboot <= boot_num) {
 cat('---------------------------------','\n')
 cat('Boot#=',iboot,'\n')
#-------------------------------------------------
# GET BOOTSTRAP SAMPLE OF DATASET
#-------------------------------------------------
  idset=c(1:nperson) 
  boot=sample(idset,nperson,replace=TRUE)

#  cat('bootstrap sample=','\n')
#  print(boot)

  dataset_boot=dataset_save[dataset_save$id==boot[1],]
  dataset_boot$oldid=dataset_boot$id
  dataset_boot$id=rep(1,times=length(dataset_boot$id))
  i=2
  while (i <=nperson) {
    sel=dataset_save[dataset_save$id==boot[i],]
    sel$oldid=sel$id
    sel$id=rep(i,times=length(sel$id))
    dataset_boot=rbind(dataset_boot,sel)
    i=i+1
  }
  n_boot=length(dataset_boot$id)
  dataset_boot=subset(dataset_boot,select=-c(oldid))

#  cat('bootstrap dataset =','\n')
#  print(dataset_boot[1:30,])
#  cat('--------------------------------------','\n')

  dataset_nonmiss_boot=dataset_nonmiss_save[dataset_nonmiss_save$id==boot[1],]
  dataset_nonmiss_boot$oldid=dataset_nonmiss_boot$id
  dataset_nonmiss_boot$id=rep(1,times=length(dataset_nonmiss_boot$id))
  i=2
  while (i <=nperson) {
    sel=dataset_nonmiss_save[dataset_nonmiss_save$id==boot[i],]
    sel$oldid=sel$id
    sel$id=rep(i,times=length(sel$id))
    dataset_nonmiss_boot=rbind(dataset_nonmiss_boot,sel)
    i=i+1
  }
  n_boot_nonmiss=length(dataset_nonmiss_boot$id)

  dataset_nonmiss_boot=subset(dataset_nonmiss_boot,select=-c(oldid))

#  cat('bootstrap dataset_nonmiss =','\n')
#  print(dataset_nonmiss_boot[1:30,])
#  cat('--------------------------------------','\n')

#  cat('future dataset =','\n')
#  print(future_dataset)
#  cat('--------------------------------------','\n')
#  cat('future nonmiss dataset =','\n')
#  print(future_dataset_nonmiss_save)
#  cat('--------------------------------------','\n')
#-------------------------------------------------
# CREATE DATASET NONMISS
#-------------------------------------------------
  dataset_boot$id=dataset_boot$id+1
  dataset_boot=smartbind(future_dataset,dataset_boot)
  dataset_nonmiss_boot$id=dataset_nonmiss_boot$id+1
  dataset_nonmiss_boot=smartbind(future_dataset_nonmiss_save,dataset_nonmiss_boot)

#  cat('bootstrap dataset with future subject=','\n')
#  print(dataset_boot[1:30,])
#  cat('--------------------------------------','\n')
#  cat('bootstrap dataset nonmiss with future subject=','\n')
#  print(dataset_nonmiss_boot[1:30,])
#  cat('--------------------------------------','\n')
#-----------------------------------------------------------------------------
#  Apply administrative censoring to Dataset_boot for integrity of bp_curr term
#-----------------------------------------------------------------------------
  dataset_boot$status=dataset_boot$status*(dataset_boot$stop <= cens_time)+0*(dataset_boot$stop > cens_time)
  dataset_boot$stop=dataset_boot$stop*(dataset_boot$stop <= cens_time)+cens_time*(dataset_boot$stop > cens_time)

#  cat('Right censoring time=',cens_time,'\n')
#  cat('Bootstrap Dataset after admin censoring=','\n')
#  print(dataset_boot[1:10,])
#-----------------------------------------------------------------------------
# Determine imputation intervals
#-----------------------------------------------------------------------------
  event_times=sort(unique(dataset_boot$stop[dataset_boot$status==1]))

#  cat('Ordered Unique Event times','\n')
#  print(event_times)


  num_imp_bp=num_imp_bp_save
  imp_times_bp=imp_times_bp_save
  num_imp_bp_new=0
  j=1
  while (j <= num_imp_bp) {
    if (any(event_times > imp_times_bp[j] & event_times <= imp_times_bp[j]+imp_inc)) {
      num_imp_bp_new=num_imp_bp_new+1
      imp_times_bp[num_imp_bp_new]=imp_times_bp[j]
    }
    j=j+1
  }
  num_imp_bp=num_imp_bp_new
  imp_times_bp=imp_times_bp[1:num_imp_bp]
  imp_times_bp=sort(unique(c(imp_times_bp,times_bp,s,s+tau)))
  num_imp_bp=length(imp_times_bp)

#  cat('Imputation times','\n')
#  print(imp_times_bp)
#  cat('------------------------------------------------','\n')
#  cat('------------------------------------------------','\n')

#-----------------------------------------------------------------------------
# Determine prediction times
#-----------------------------------------------------------------------------
  predtimes=imp_times_bp[2:num_imp_bp]
  index_s=which(predtimes==s)

#-----------------------------------------------------------------------------
# Create Multiple Observations per Subject dataset
#-----------------------------------------------------------------------------
  n=nperson+1
  maxobs=(nperson+1)*num_imp_bp
  d1=c(dataset_boot$id,rep(0,times=maxobs-n))
  d2=c(dataset_boot$start,rep(0,times=maxobs-n))
  d3=c(dataset_boot$stop,rep(0,times=maxobs-n))
  d4=c(dataset_boot$status,rep(0,times=maxobs-n))
  d5=c(dataset_boot$baseage,rep(0,times=maxobs-n))
  d6=c(dataset_boot$lastobs,rep(0,times=maxobs-n))

  z=add(n,maxobs,num_imp_bp,imp_times_bp,
   d1,d2,d3,d4,d5,d6)

  dataset_boot=data.frame(id=z[[5]],start=z[[6]],stop=z[[7]],status=z[[8]],baseage=z[[9]],
    bp_curr1=rep(0,times=maxobs),bp_curr2=rep(0,times=maxobs),bp_curr3=rep(0,times=maxobs),
    bp_curr4=rep(0,times=maxobs),bp_curr5=rep(0,times=maxobs),lastobs=z[[10]])
  n=z[[1]]
  dataset_boot=dataset_boot[1:n,]
  dataset_boot=dataset_boot[dataset_boot$start < dataset_boot$stop,]
  n=length(dataset_boot$id)

#  cat('Bootstrap Dataset after adding observations=','\n')
#  print(dataset_boot[1:30,])
#  cat('------------------------------------------------','\n')
#--------------------------------------------------------------------------
#
#  Create survival cohorts from dataset_nonmiss_boot for each imputation time
#
#--------------------------------------------------------------------------
  n_nonmiss=length(dataset_nonmiss_boot$id)
  dataset_nonmiss_boot$lastobs=rep(0,times=n_nonmiss)
  i=1
  while ( i <= n_nonmiss) {
    id=dataset_nonmiss_boot$id[i]
    if (i==n_nonmiss) {
      dataset_nonmiss_boot$lastobs[i]=1
    } else {
      if (id != dataset_nonmiss_boot$id[i+1]) {dataset_nonmiss_boot$lastobs[i]=1}
    }  
    i=i+1
  }

#  cat('Bootstrap Longitudinal Dataset with lastobs','\n')
#  print(dataset_nonmiss_boot[1:20,])
#  cat('------------------------------------------------','\n')

#  cat('Dataset after admin censoring for ID=19=','\n')
#  print(dataset[dataset$id==19,])

  survive_time=dataset_boot$stop[dataset_boot$lastobs==1]
  survive_time[1]=s+tau

#  cat('Last Time Known to be Event Free','\n')
#  print(survive_time)

  future_dataset_nonmiss_save=dataset_nonmiss_boot[dataset_nonmiss_boot$id==1,]
  dataset_nonmiss_save=dataset_nonmiss_boot[dataset_nonmiss_boot$id> 1,]
  dataset_nonmiss_save$id=dataset_nonmiss_save$id-1
  n_nonmiss=length(dataset_nonmiss_boot$id)

  survive_ind=rep(0,times=n_nonmiss*num_imp_bp)
  dim(survive_ind)=c(n_nonmiss,num_imp_bp)
  survive_ind_past=rep(0,times=n_nonmiss*num_imp_bp)
  dim(survive_ind_past)=c(n_nonmiss,num_imp_bp)
  j=1
  while (j <=num_imp_bp) {
    survive_ind[,j]=(round(survive_time[dataset_nonmiss_boot$id],digits=5) >= round(imp_times_bp[j],digits=5))
    survive_ind_past[,j]=(round(survive_time[dataset_nonmiss_boot$id],digits=5) >= round(imp_times_bp[j],digits=5) & dataset_nonmiss_boot$id > nperson_future)
    j=j+1
  }
#  cat('Survive_ind=','\n')
#  print(survive_ind[1:20,])
#-----------------------------------------------------------------------------
#  Left Censor Dataset_boot for Imputation Models to have sufficient sample sizes
#-----------------------------------------------------------------------------
  dataset_boot=dataset_boot[dataset_boot$stop >= left_cens_time-.000001,]
  dataset_boot$start=dataset_boot$start*(dataset_boot$start >= left_cens_time)+left_cens_time*(dataset_boot$start < left_cens_time)

#  cat('Dataset_boot after left censoring=','\n')
#  print(dataset_boot[1:30,])
#  cat('------------------------------------------------','\n')
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
# Fit Imputation Models for Biomarker Process to Boostrap Dataset
# Imputation Models based on Survival Cohorts and Times prior to Current Time
#
#------------------------------------------------------------------------------
  dataset_nonmiss_boot$obstime1=dataset_nonmiss_boot$obstime
  dataset_nonmiss_boot$obstime2=dataset_nonmiss_boot$obstime**2
  dataset_nonmiss_boot$obstime3=dataset_nonmiss_boot$obstime**3
  dataset_nonmiss_boot$id=as.factor(dataset_nonmiss_boot$id)

#  cat('Longitudinal Dataset with nonmissing bp_curr & powers of time =','\n')
#  print(dataset_nonmiss_boot[1:20,])

  if (icox1_ind==1) {
#---------------------------------------------------------------------
# START ICOX1: METHOD THAT USES BLUP OF FUTURE SUBJECT
#---------------------------------------------------------------------
  cat('START ICOX1','\n')
  cat('-----------------------------------------------','\n')
  cat('-----------------------------------------------','\n')

  coef_impute1=list()
  coef_impute2=list()
  coef_impute3=list()
  coef_impute4=list()
  coef_impute5=list()
  random_coef_impute_id1=list()
  random_coef_impute_id2=list()
  random_coef_impute_id3=list()
  random_coef_impute_id4=list()
  random_coef_impute_id5=list()
  select_id=list()
  resid_error_sd=rep(0,times=num_imp_bp)

  j=1
  while (j <= num_imp_bp) {
#    cat('--------------------------------------------','\n')
#    cat('--------------------------------------------','\n')
#    cat('Imp model for imp time =',imp_times_bp[j],'\n')
#    cat('--------------------------------------------','\n')

    select=which(survive_ind[,j]==1 & dataset_nonmiss_boot$obstime <= imp_times_bp[j]+window_width)
    select_id[[j]]=unique(dataset_nonmiss_boot$id[which(survive_ind[,j]==1)])
    n_select=length(select)
    n_select_id=length(select_id[[j]])
    dataset_nonmiss_select=dataset_nonmiss_boot[select,]
    dataset_nonmiss_select$numgroup=as.numeric(dataset_nonmiss_select$group)
#    cat('Subjects known to survive to imp_time =',select_id[[j]],'\n')
#    cat('# of Subjects known to survive to imp_time =',n_select_id,'\n')
    jj=1
    names=""
    while (jj <= 3) {
      names=paste(names,"+obstime",jj,sep="")
      names=paste(names,sep=" ")
      jj=jj+1
    }
    names=paste(names,"+ (obstime | id)",sep="")
    (formula=as.formula(paste("bp_curr~baseage",paste(names,collapse="+"))))
    fit1=lmer(formula,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=="TYPE1")
    coef_impute1[[j]]=fixef(fit1)
    coef_impute1[[j]][is.na(coef_impute1[[j]])]=0
    random_coef_impute_id1[[j]]=random.effects(fit1)[[1]]

    if (dim_bp>=2) {
      fit2=lmer(formula,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=="TYPE2")
      coef_impute2[[j]]=fixef(fit2)
      coef_impute2[[j]][is.na(coef_impute2[[j]])]=0
      random_coef_impute_id2[[j]]=random.effects(fit2)[[1]]
    }

    if (dim_bp>=3) {
      fit3=lmer(formula,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=="TYPE3")
      coef_impute3[[j]]=fixef(fit3)
      coef_impute3[[j]][is.na(coef_impute3[[j]])]=0
      random_coef_impute_id3[[j]]=random.effects(fit3)[[1]]
    }

    if (dim_bp>=4) {
      fit4=lmer(formula,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=="TYPE4")
      coef_impute4[[j]]=fixef(fit4)
      coef_impute4[[j]][is.na(coef_impute4[[j]])]=0
      random_coef_impute_id4[[j]]=random.effects(fit4)[[1]]
    }

    if (dim_bp>=5) {
      fit5=lmer(formula,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=="TYPE5")
      coef_impute5[[j]]=fixef(fit5)
      coef_impute5[[j]][is.na(coef_impute5[[j]])]=0
      random_coef_impute_id5[[j]]=random.effects(fit5)[[1]]
    }

#    resid_error_sd[j]=sigma(fit1)
#    cat('Imp model residual error sd=',resid_error_sd[j],'\n')
#    cat('--------------------------------------------','\n')


    j=j+1
  }
#--------------------------------------------------------------------------
#
# All Patients
# Replace ALL values of Biomarker Process with fitted value from model
#
  dataset_impute=dataset_boot
  n=length(dataset_impute$id)

#    cat('Bootstrap Dataset before imputed BP_curr =','\n')
#    print(dataset_impute[1:50,])

  i=1
  while (i <= n) {
#      cat('-----------------------','\n')
#      cat('i=',i,'\n')
    id=dataset_impute$id[i]
#      cat('id=',id,'\n')
#      cat('start[i]=',dataset_impute[[k]]$start[i],'\n')
#      cat('round(start[i])=',round(dataset_impute[[k]]$start[i],digits=digits_imp_inc),'\n')
#      cat('imp_times_bp=',imp_times_bp[40:45],'\n')
    j=which(round(imp_times_bp,digits=digits_imp_inc)==round(dataset_impute$start[i],digits=digits_imp_inc))
#      cat('j=',j,'\n')
      fit1=coef_impute1[[j]][1]
      fit1=fit1+coef_impute1[[j]][2]*dataset_impute$baseage[i]
      fit1=fit1+coef_impute1[[j]][3]*dataset_impute$start[i]
      fit1=fit1+coef_impute1[[j]][4]*dataset_impute$start[i]**2
      fit1=fit1+coef_impute1[[j]][5]*dataset_impute$start[i]**3

      if (dim_bp>=2) {
        fit2=coef_impute2[[j]][1]
        fit2=fit2+coef_impute2[[j]][2]*dataset_impute$baseage[i]
        fit2=fit2+coef_impute2[[j]][3]*dataset_impute$start[i]
        fit2=fit2+coef_impute2[[j]][4]*dataset_impute$start[i]**2
        fit2=fit2+coef_impute2[[j]][5]*dataset_impute$start[i]**3
      }

      if (dim_bp>=3) {
        fit3=coef_impute3[[j]][1]
        fit3=fit3+coef_impute3[[j]][2]*dataset_impute$baseage[i]
        fit3=fit3+coef_impute3[[j]][3]*dataset_impute$start[i]
        fit3=fit3+coef_impute3[[j]][4]*dataset_impute$start[i]**2
        fit3=fit3+coef_impute3[[j]][5]*dataset_impute$start[i]**3
      }

      if (dim_bp>=4) {
        fit4=coef_impute4[[j]][1]
        fit4=fit4+coef_impute4[[j]][2]*dataset_impute$baseage[i]
        fit4=fit4+coef_impute4[[j]][3]*dataset_impute$start[i]
        fit4=fit4+coef_impute4[[j]][4]*dataset_impute$start[i]**2
        fit4=fit4+coef_impute4[[j]][5]*dataset_impute$start[i]**3
      }

      if (dim_bp>=5) {
        fit5=coef_impute5[[j]][1]
        fit5=fit5+coef_impute5[[j]][2]*dataset_impute$baseage[i]
        fit5=fit5+coef_impute5[[j]][3]*dataset_impute$start[i]
        fit5=fit5+coef_impute5[[j]][4]*dataset_impute$start[i]**2
        fit5=fit5+coef_impute5[[j]][5]*dataset_impute$start[i]**3
      }

      id_num=which(select_id[[j]]==dataset_impute$id[i])

      fit1=fit1+random_coef_impute_id1[[j]][id_num,1]
      fit1=fit1+random_coef_impute_id1[[j]][id_num,2]*dataset_impute$start[i]
      dataset_impute$bp_curr1[i]=fit1

      if (dim_bp>=2) {
        fit2=fit2+random_coef_impute_id2[[j]][id_num,1]
        fit2=fit2+random_coef_impute_id2[[j]][id_num,2]*dataset_impute$start[i]
        dataset_impute$bp_curr2[i]=fit2
      }

      if (dim_bp>=3) {
        fit3=fit3+random_coef_impute_id3[[j]][id_num,1]
        fit3=fit3+random_coef_impute_id3[[j]][id_num,2]*dataset_impute$start[i]
        dataset_impute$bp_curr3[i]=fit3
      }

      if (dim_bp>=4) {
        fit4=fit4+random_coef_impute_id4[[j]][id_num,1]
        fit4=fit4+random_coef_impute_id4[[j]][id_num,2]*dataset_impute$start[i]
        dataset_impute$bp_curr4[i]=fit4
      }

      if (dim_bp>=5) {
        fit5=fit5+random_coef_impute_id5[[j]][id_num,1]
        fit5=fit5+random_coef_impute_id5[[j]][id_num,2]*dataset_impute$start[i]
        dataset_impute$bp_curr5[i]=fit5
      }

    i=i+1
  }

#    cat('Bootstrap Dataset after imputed BP_curr =','\n')
#    print(dataset_impute[[k]][1:50,])
#--------------------------------------------------------------------------
#
# Fit EMPIRICAL BAYES COX MODEL TO BOOTSTRAP DATA
#
#--------------------------------------------------------------------------
  dataset_cox=dataset_impute[[k]][dataset_impute$id > 1,]
  dataset_cox=dataset_cox[dataset_cox$start < dataset_cox$stop-.000001,]

  cat('------------------------------------------------------------','\n')
  cat('EMP Cox model uses apriori FIXED TERMS','\n')
  cat('------------------------------------------------------------','\n')
  if (dim_bp==1) {
    (ecformula=as.formula(paste("Surv(start,stop,status)~baseage+bp_curr1")))
  } 
  if (dim_bp==2) {
    (ecformula=as.formula(paste("Surv(start,stop,status)~baseage+bp_curr1+bp_curr2")))
  } 
  if (dim_bp==3) {
    (ecformula=as.formula(paste("Surv(start,stop,status)~baseage+bp_curr1+bp_curr2+bp_curr3")))
  } 
  if (dim_bp==4) {
    (ecformula=as.formula(paste("Surv(start,stop,status)~baseage+bp_curr1+bp_curr2+bp_curr3+bp_curr4")))
  } 
  if (dim_bp==5) {
    (ecformula=as.formula(paste("Surv(start,stop,status)~baseage+bp_curr1+bp_curr2+bp_curr3+bp_curr4+bp_curr5")))
  } 
  efit=try(coxph(ecformula,data=dataset_cox,x=TRUE))
  cat('------------------------------------------------------------','\n')
  cat('IMP Cox model with ',dim_bp,' biomarker terms is','\n')
  print(efit)
  coxcoefvalues=abs(coef(summary(efit))[2:(dim_bp+1),1])
  coxcoefvalues[is.na(coxcoefvalues)]=0.0
#--------------------------------------------------------------------------
#
# Remove observations past s+tau in FUTURE_DATASET
#
#--------------------------------------------------------------------------
  future_dataset_impute=dataset_impute[dataset_impute$id==1 & round(dataset_impute$start,digits=5) < s+tau,]
  n_future=length(future_dataset_impute$id)
  nperson_future=1

#    cat('Final Future Dataset with imputed BP process =','\n')
#    print(future_dataset_impute[[k]])
#--------------------------------------------------------------------------
#
# FORTRAN VERSION OF LOOP TO CALCULATE CONDITIONAL SURVIVAL PROBS
#
#--------------------------------------------------------------------------
  cond_surv_pred_st=rep(0,times=nperson_future)
#    cat('Call predictSurvProb','\n')
  future_surv=predictSurvProb(efit,newdata=future_dataset_impute,times=predtimes)
#    cat('Return from predictSurvProb','\n')
  npredtimes=length(predtimes)
  future_surv_short=future_surv

  idvec=future_dataset_impute$id

#  cat('survival to predtimes for subjects in Bootstrap future_dataset according to EMPIRICAL COX MODEL=','\n')
#  print(future_surv_short)

  if (all(is.finite(future_surv_short))) {
#    cat('Call condprob','\n')
    z=condprob(nperson_future,n_future,npredtimes,index_s,future_surv_short,idvec,cond_surv_pred_st)
#    cat('return from condprob','\n')

    cond_surv_pred_st=z[[7]]

#      cat('FORTRAN: Cond. survival to S+TAU given survival to S for subjects in Bootstrap future_dataset according to EMPIRICAL COX MODEL=','\n')
#      print(cond_surv_pred_st)

    icox1_cond_surv_pred_imp=mean(cond_surv_pred_st)

#    cat('FORTRAN: Avg. Cond. survival to S+TAU given survival to S for subjects in Bootstrap future_dataset according to EMPIRICAL COX MODEL=','\n')
#      print(icox1_cond_surv_pred_imp[k])
  } else {
#      cat('EMPIRICAL COX MODEL DOES NOT GIVE PROB ESTIMATE FOT BOOTSTRAP DATASET','\n')
    icox1_cond_surv_pred_imp=NA
  }


  cat('---------------------------------------------------','\n')
  cat('---------------------------------------------------','\n')
  if (any(is.finite(icox1_cond_surv_pred_imp))) {
    icox1_cond_surv_pred_boot[iboot]=mean(icox1_cond_surv_pred_imp,na.rm=TRUE)
    cat('FINAL Estimate Cond. survival to S+TAU given survival to S for future subj in Bootstrap according to IMP1 COX MODEL=','\n')
    print(icox1_cond_surv_pred_boot[iboot])
  } else {
    cat('FINAL IMP1 COX MODEL DOES NOT GIVE PROB ESTIMATE FOR BOOTSTRAP','\n')
    icox1_cond_surv_pred_boot[iboot]=NA
  }


  }
#--------------------------------------------------------------------------
# END ICOX1
#--------------------------------------------------------------------------






  if (icox2_ind==1) {
#------------------------------------------------------------------------------------------
# START ICOX2: METHOD THAT USES BLUP FROM NUM_SIMILAR PAST SUBJECTS MOST SIMILAR BY S-OBSERVED VALUES
#------------------------------------------------------------------------------------------
  cat('---------------------------------------------','\n')
  cat('START ICOX2','\n')
  cat('---------------------------------------------','\n')
  cat('---------------------------------------------','\n')

  if (icox1_ind==0) {
    coef_impute1=list()
    coef_impute2=list()
    coef_impute3=list()
    coef_impute4=list()
    coef_impute5=list()
    random_coef_impute_id1=list()
    random_coef_impute_id2=list()
    random_coef_impute_id3=list()
    random_coef_impute_id4=list()
    random_coef_impute_id5=list()
    select_id=list()
  }
  coef_impute1_past=list()
  coef_impute2_past=list()
  coef_impute3_past=list()
  coef_impute4_past=list()
  coef_impute5_past=list()
  random_coef_impute_id1_past=list()
  random_coef_impute_id2_past=list()
  random_coef_impute_id3_past=list()
  random_coef_impute_id4_past=list()
  random_coef_impute_id5_past=list()
  select_id_past=list()
  idmatch2=list()

  j=1
  while (j <= num_imp_bp) {

    if (icox1_ind==0) {
      select=which(survive_ind[,j]==1 & dataset_nonmiss_boot$obstime <= imp_times_bp[j]+window_width)
      select_id[[j]]=unique(dataset_nonmiss_boot$id[which(survive_ind[,j]==1)])
      n_select=length(select)
      n_select_id=length(select_id[[j]])
      dataset_nonmiss_select=dataset_nonmiss_boot[select,]
#      cat('# of Subjects known to survive to imp_time =',n_select_id,'\n')
#    cat('Subjects known to survive to imp_time =',select_id[[j]],'\n')
      jj=1
      names=""
      while (jj <= 3) {
        names=paste(names,"+obstime",jj,sep="")
        names=paste(names,sep=" ")
        jj=jj+1
      }
      names=paste(names,"+ (obstime | id)",sep="")
      (formula=as.formula(paste("bp_curr~baseage",paste(names,collapse="+"))))
#-------------------------------------------------
# TEST LINEAR FIXED EFFECTS
#      cat('--------------------------------------------','\n')
#      cat('--------------------------------------------','\n')
#      cat('Imp model for imp time =',imp_times_bp[j],'\n')
#      cat('--------------------------------------------','\n')
      fit1=lmer(formula,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=="TYPE1")
      coef_impute1[[j]]=fixef(fit1)
      coef_impute1[[j]][is.na(coef_impute1[[j]])]=0
      random_coef_impute_id1[[j]]=random.effects(fit1)[[1]]

      if (dim_bp>=2) {
        fit2=lmer(formula,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=="TYPE2")
        coef_impute2[[j]]=fixef(fit2)
        coef_impute2[[j]][is.na(coef_impute2[[j]])]=0
        random_coef_impute_id2[[j]]=random.effects(fit2)[[1]]
      }

      if (dim_bp>=3) {
        fit3=lmer(formula,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=="TYPE3")
        coef_impute3[[j]]=fixef(fit3)
        coef_impute3[[j]][is.na(coef_impute3[[j]])]=0
        random_coef_impute_id3[[j]]=random.effects(fit3)[[1]]
      }

      if (dim_bp>=4) {
        fit4=lmer(formula,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=="TYPE4")
        coef_impute4[[j]]=fixef(fit4)
        coef_impute4[[j]][is.na(coef_impute4[[j]])]=0
        random_coef_impute_id4[[j]]=random.effects(fit4)[[1]]
      }

      if (dim_bp>=5) {
        fit5=lmer(formula,data=dataset_nonmiss_select, na.action=na.exclude,subset=dataset_nonmiss_select$bptype=="TYPE5")
        coef_impute5[[j]]=fixef(fit5)
        coef_impute5[[j]][is.na(coef_impute5[[j]])]=0
        random_coef_impute_id5[[j]]=random.effects(fit5)[[1]]
      }
    }

    if (j==s_index) {
      select=which(survive_ind_past[,j]==1 & dataset_nonmiss_boot$obstime <= imp_times_bp[j]+window_width)
      select_id_past[[j]]=unique(dataset_nonmiss_boot$id[which(survive_ind_past[,j]==1)])
      n_select=length(select)
      n_select_id_past=length(select_id_past[[j]])
      dataset_nonmiss_select=dataset_nonmiss_boot[select,]
      baseage_select=aggregate(dataset_nonmiss_select$baseage,by=list(dataset_nonmiss_select$id,dataset_nonmiss_select$bptype),FUN=tail,n=1)
      baseage_select=baseage_select$x[baseage_select$Group.2=='TYPE1']

#      cat('baseage_select=','\n')
#      print(baseage_select)
    }

    j=j+1
  }
#  currentseed=.Random.seed[2]
#  currentseed=.Random.seed[currentseed+2]
#  cat('After Longitudinal Models fit current R seed=',currentseed,'\n')

#--------------------------------------------------------------------------
#
# Past Patients
# Replace ALL values of Biomarker Process with fitted value from model
#
  dataset_impute2=dataset_boot
  n=length(dataset_impute2$id)

#  cat('Dataset before imputed BP_curr =','\n')
#  print(dataset_impute[[k]][1:20,])

  i=1
  while (i <= n) {
#    cat('-----------------------','\n')
#    cat('i=',i,'\n')
    id=dataset_impute2$id[i]
#    cat('id=',id,'\n')
#    cat('start[i]=',dataset_impute[[k]]$start[i],'\n')
#    cat('round(start[i])=',round(dataset_impute[[k]]$start[i],digits=digits_imp_inc),'\n')
#    cat('imp_times_bp=',imp_times_bp[40:45],'\n')
    j=which(round(imp_times_bp,digits=digits_imp_inc)==round(dataset_impute2$start[i],digits=digits_imp_inc))
#    cat('j=',j,'\n')
#    cat('time=',imp_times_bp[j],'\n')
    if (id>1 & icox1_ind==0) {
#      if (j > s_index) {j=s_index}
#      cat('coef_impute[[j]]=',coef_impute[[j]],'\n')
#      cat('dataset_impute[[k]]$baseage[i]=',dataset_impute[[k]]$baseage[i],'\n')
#      cat('fit=',fit,'\n')
#-------------------------------------------------
    if (j <= s_index) {
      fit1=coef_impute1[[j]][1]
      fit1=fit1+coef_impute1[[j]][2]*dataset_impute2$baseage[i]
      fit1=fit1+coef_impute1[[j]][3]*dataset_impute2$start[i]
      fit1=fit1+coef_impute1[[j]][4]*dataset_impute2$start[i]**2
      fit1=fit1+coef_impute1[[j]][5]*dataset_impute2$start[i]**3

      if (dim_bp>=2) {
        fit2=coef_impute2[[j]][1]
        fit2=fit2+coef_impute2[[j]][2]*dataset_impute2$baseage[i]
        fit2=fit2+coef_impute2[[j]][3]*dataset_impute2$start[i]
        fit2=fit2+coef_impute2[[j]][4]*dataset_impute2$start[i]**2
        fit2=fit2+coef_impute2[[j]][5]*dataset_impute2$start[i]**3
      }

      if (dim_bp>=3) {
        fit3=coef_impute3[[j]][1]
        fit3=fit3+coef_impute3[[j]][2]*dataset_impute2$baseage[i]
        fit3=fit3+coef_impute3[[j]][3]*dataset_impute2$start[i]
        fit3=fit3+coef_impute3[[j]][4]*dataset_impute2$start[i]**2
        fit3=fit3+coef_impute3[[j]][5]*dataset_impute2$start[i]**3
      }

      if (dim_bp>=4) {
        fit4=coef_impute3[[j]][1]
        fit4=fit4+coef_impute4[[j]][2]*dataset_impute2$baseage[i]
        fit4=fit4+coef_impute4[[j]][3]*dataset_impute2$start[i]
        fit4=fit4+coef_impute4[[j]][4]*dataset_impute2$start[i]**2
        fit4=fit4+coef_impute4[[j]][5]*dataset_impute2$start[i]**3
      }

      if (dim_bp>=5) {
        fit5=coef_impute3[[j]][1]
        fit5=fit5+coef_impute5[[j]][2]*dataset_impute2$baseage[i]
        fit5=fit5+coef_impute5[[j]][3]*dataset_impute2$start[i]
        fit5=fit5+coef_impute5[[j]][4]*dataset_impute2$start[i]**2
        fit5=fit5+coef_impute5[[j]][5]*dataset_impute2$start[i]**3
      }

#      cat('after fixed effects fit1=',fit1,'\n')
#-------------------------------------------------
# TEST AIDS DATA 
#      fit=coef_impute[[j]][1]
#      if (dataset_impute[[k]]$baseage[i]==1) {fit=fit+coef_impute[[j]][3]*dataset_impute[[k]]$start[i]}
#-------------------------------------------------
      id_num=which(select_id[[j]]==dataset_impute2$id[i])
      fit1=fit1+random_coef_impute_id1[[j]][id_num,1]
      fit1=fit1+random_coef_impute_id1[[j]][id_num,2]*dataset_impute2$start[i]
      dataset_impute2$bp_curr1[i]=fit1
      
      if (dim_bp >=2) {
        fit2=fit2+random_coef_impute_id2[[j]][id_num,1]
        fit2=fit2+random_coef_impute_id2[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr2[i]=fit2
      }

      if (dim_bp >=3) {
        fit3=fit3+random_coef_impute_id3[[j]][id_num,1]
        fit3=fit3+random_coef_impute_id3[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr3[i]=fit3
      }

      if (dim_bp >=4) {
        fit4=fit4+random_coef_impute_id4[[j]][id_num,1]
        fit4=fit4+random_coef_impute_id4[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr4[i]=fit4
      }

      if (dim_bp >=5) {
        fit5=fit5+random_coef_impute_id5[[j]][id_num,1]
        fit5=fit5+random_coef_impute_id5[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr5[i]=fit5
      }
    }

    if (j>s_index) {
      id_num=which(select_id[[j]]==dataset_impute2$id[i])
      if (id==1) {
        id_num=which(select_id[[j]] %in% idmatch2[[j]])
      }
      fit1=coef_impute1[[j]][1]
     
      if (id==1) {
        fit1=fit1+coef_impute1[[j]][2]*mean(baseage_id[idmatch2[[j]]])
      } else {
        fit1=fit1+coef_impute1[[j]][2]*dataset_impute2$baseage[i]
      }
      fit1=fit1+coef_impute1[[j]][3]*dataset_impute2$start[i]
      fit1=fit1+coef_impute1[[j]][4]*dataset_impute2$start[i]**2
      fit1=fit1+coef_impute1[[j]][5]*dataset_impute2$start[i]**3

      if (dim_bp>=2) {
        fit2=coef_impute2[[j]][1]
        fit2=fit2+coef_impute2[[j]][2]*dataset_impute2$baseage[i]
        fit2=fit2+coef_impute2[[j]][3]*dataset_impute2$start[i]
        fit2=fit2+coef_impute2[[j]][4]*dataset_impute2$start[i]**2
        fit2=fit2+coef_impute2[[j]][5]*dataset_impute2$start[i]**3
      }

      if (dim_bp>=3) {
        fit3=coef_impute3[[j]][1]
        fit3=fit3+coef_impute3[[j]][2]*dataset_impute2$baseage[i]
        fit3=fit3+coef_impute3[[j]][3]*dataset_impute2$start[i]
        fit3=fit3+coef_impute3[[j]][4]*dataset_impute2$start[i]**2
        fit3=fit3+coef_impute3[[j]][5]*dataset_impute2$start[i]**3
      }

      if (dim_bp>=4) {
        fit4=coef_impute3[[j]][1]
        fit4=fit4+coef_impute4[[j]][2]*dataset_impute2$baseage[i]
        fit4=fit4+coef_impute4[[j]][3]*dataset_impute2$start[i]
        fit4=fit4+coef_impute4[[j]][4]*dataset_impute2$start[i]**2
        fit4=fit4+coef_impute4[[j]][5]*dataset_impute2$start[i]**3
      }

      if (dim_bp>=5) {
        fit5=coef_impute3[[j]][1]
        fit5=fit5+coef_impute5[[j]][2]*dataset_impute2$baseage[i]
        fit5=fit5+coef_impute5[[j]][3]*dataset_impute2$start[i]
        fit5=fit5+coef_impute5[[j]][4]*dataset_impute2$start[i]**2
        fit5=fit5+coef_impute5[[j]][5]*dataset_impute2$start[i]**3
      }

#      cat('after fixed effects fit1=',fit1,'\n')
#-------------------------------------------------
# TEST AIDS DATA 
#      fit=coef_impute[[j]][1]
#      if (dataset_impute[[k]]$baseage[i]==1) {fit=fit+coef_impute[[j]][3]*dataset_impute[[k]]$start[i]}
#-------------------------------------------------
      fit1=fit1+mean(random_coef_impute_id1[[j]][id_num,1])
      fit1=fit1+mean(random_coef_impute_id1[[j]][id_num,2]*dataset_impute2$start[i])
      dataset_impute2$bp_curr1[i]=fit1
      
      if (dim_bp >=2) {
        fit2=fit2+random_coef_impute_id2[[j]][id_num,1]
        fit2=fit2+random_coef_impute_id2[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr2[i]=fit2
      }

      if (dim_bp >=3) {
        fit3=fit3+random_coef_impute_id3[[j]][id_num,1]
        fit3=fit3+random_coef_impute_id3[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr3[i]=fit3
      }

      if (dim_bp >=4) {
        fit4=fit4+random_coef_impute_id4[[j]][id_num,1]
        fit4=fit4+random_coef_impute_id4[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr4[i]=fit4
      }

      if (dim_bp >=5) {
        fit5=fit5+random_coef_impute_id5[[j]][id_num,1]
        fit5=fit5+random_coef_impute_id5[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr5[i]=fit5
      }

    }
    }

    i=i+1
  }
#  cat('Dataset with imputed BP_curr for all patients=','\n')
#  print(dataset_impute[[k]][1:50,])

#--------------------------------------------------------------------------
#
# Fit EMPIRICAL BAYES COX MODEL
#
#--------------------------------------------------------------------------
  if (icox1_ind==0) {
    dataset_cox=dataset_impute2[dataset_impute2$id > nperson_future,]
    dataset_cox=dataset_cox[dataset_cox$start < dataset_cox$stop-.000001,]

#  cat('num_imp_bp=',num_imp_bp,'\n')
#  cat('Dataset_cox=','\n')
#  print(dataset_cox[1:20,])

    cat('------------------------------------------------------------','\n')
    cat('EMP Cox model uses apriori FIXED TERMS','\n')
    cat('------------------------------------------------------------','\n')
    if (dim_bp==1) {
      (ecformula=as.formula(paste("Surv(start,stop,status)~baseage+bp_curr1")))
    }
    if (dim_bp==2) {
      (ecformula=as.formula(paste("Surv(start,stop,status)~baseage+bp_curr1+bp_curr2")))
    }
    if (dim_bp==3) {
      (ecformula=as.formula(paste("Surv(start,stop,status)~baseage+bp_curr1+bp_curr2+bp_curr3")))
    }
    if (dim_bp==4) {
      (ecformula=as.formula(paste("Surv(start,stop,status)~baseage+bp_curr1+bp_curr2+bp_curr3+bp_curr4")))
    }
    if (dim_bp==5) {
      (ecformula=as.formula(paste("Surv(start,stop,status)~baseage+bp_curr1+bp_curr2+bp_curr3+bp_curr4+bp_curr5")))
    }
    efit=try(coxph(ecformula,data=dataset_cox,x=TRUE))
    cat('------------------------------------------------------------','\n')
    cat('EMP Cox model with baseage+bp_curr is','\n')
    print(efit)
    coxcoefvalues=abs(coef(summary(efit))[2:(dim_bp+1),1])
    coxcoefvalues[is.na(coxcoefvalues)]=0.0
  }
#--------------------------------------------------------------------------
# Calculate Matches
#--------------------------------------------------------------------------
  j=s_index+1
  while (j <= num_imp_bp) {

    select=which(survive_ind[,j]==1 & dataset_nonmiss_boot$obstime <= imp_times_bp[j]+window_width)
    select_id[[j]]=unique(dataset_nonmiss_boot$id[which(survive_ind[,j]==1)])
    n_select=length(select)
    n_select_id=length(select_id[[j]])
    dataset_nonmiss_select=dataset_nonmiss_boot[select,]

    z=dataset_nonmiss_select$bp_curr[abs(dataset_nonmiss_select$obstime-s)<=window_width & as.numeric(dataset_nonmiss_select$id)>1 & dataset_nonmiss_select$bptype=='TYPE1']
    zid=dataset_nonmiss_select$id[abs(dataset_nonmiss_select$obstime-s)<=window_width & as.numeric(dataset_nonmiss_select$id)>1 & dataset_nonmiss_select$bptype=='TYPE1']
#      cat('ICOX2 Candidate bp_curr values=','\n')
#      print(z)
    dev=abs(dataset_nonmiss_select$bp_curr[abs(dataset_nonmiss_select$obstime-s)<=window_width & as.numeric(dataset_nonmiss_select$id)>1 & dataset_nonmiss_select$bptype=='TYPE1'] -bp_s)
#      cat('deviations =',dev,'\n')
#      cat('Corresponding IDs=','\n')
#      print(zid)
#      imp_bp_s=coef_impute1[[s_index]][1]
#      imp_bp_s=imp_bp_s+coef_impute1[[s_index]][2]*dataset_nonmiss$baseage[1]
#      imp_bp_s=imp_bp_s+coef_impute1[[s_index]][3]*s
#      imp_bp_s=imp_bp_s+coef_impute1[[s_index]][4]*s**2
#      imp_bp_s=imp_bp_s+coef_impute1[[s_index]][5]*s**3
#      imp_bp_s=imp_bp_s+random_coef_impute_id1[[s_index]][1,1]

#      cat('imp_bp_s0 =',imp_bp_s0,'\n')
#      cat('imp_bp_s1 =',imp_bp_s1,'\n')

#      cat('select_id[[s_index]]=','\n')
#      print(select_id[[s_index]])
#      cat('select_id[[j]]=','\n')
#      print(select_id[[j]])
#      temp=select_id_past[[s_index]] %in% select_id[[j]]
#      cat('temp=','\n')
#      print(temp)
    idtemp=select_id_past[[s_index]][select_id_past[[s_index]] %in% select_id[[j]]]

#      cat('idtemp=','\n')
#      print(idtemp)

    tempwhich=match(select_id[[j]], select_id[[s_index]])
    tempwhich=tempwhich[!is.na(tempwhich)]
    tempwhich=tempwhich[tempwhich != 1]

#      cat('tempwhich= position in select_id[[s_index]] for those matching from select_id[[j]] =','\n')
#      print(tempwhich)
 
    dev_obs=rep(999999,times=length(idtemp))
    id_obs=match(idtemp,zid)
    id_ind=!is.na(id_obs)
    id_obs=id_obs[!is.na(id_obs)]

#      cat('Indecies of IDs with observed value at s and also in jth model=','\n')
#      print(id_obs)
#      id_ind=which(zid %in% tempwhich)
#      cat('Index of IDs with observed value at s and also in jth model=','\n')
#      print(id_ind)
    dev_obs[id_ind]=dev[id_obs]
#      cat('New deviations =',dev_obs,'\n')

#      cat('baseage_select[tempwhich-1]=','\n')
#      print(baseage_select[tempwhich-1])

#      imp_bp_j=coef_impute1[[s_index]][1]
#      imp_bp_j=imp_bp_j+coef_impute1[[s_index]][2]*baseage_select[tempwhich-1]
#      imp_bp_j=imp_bp_j+coef_impute1[[s_index]][3]*s
#      imp_bp_j=imp_bp_j+coef_impute1[[s_index]][4]*s**2
#      imp_bp_j=imp_bp_j+coef_impute1[[s_index]][5]*s**3
#      imp_bp_j=imp_bp_j+random_coef_impute_id1[[s_index]][tempwhich,1]
    imp_bp_s0=random_coef_impute_id1[[s_index]][1,1]
    imp_bp_s1=random_coef_impute_id1[[s_index]][1,2]
    imp_bp_j0=random_coef_impute_id1[[s_index]][tempwhich,1]
    imp_bp_j1=random_coef_impute_id1[[s_index]][tempwhich,2]

#      cat('imp_bp_j0 =','\n')
#      print(imp_bp_j0)

#      dev0=abs(imp_bp_s0-imp_bp_j0)+dev_obs
#      dev1=abs(imp_bp_s1-imp_bp_j1)+dev_obs
#    dev=abs(imp_bp_s0-imp_bp_j0)+abs(imp_bp_s1-imp_bp_j1)+dev_obs
    dev=dev_obs+coxcoefvalues[1]*(abs(imp_bp_s0-imp_bp_j0)+abs(imp_bp_s1-imp_bp_j1))
    if (dim_bp>=2) {
      imp_bp2_s0=random_coef_impute_id2[[s_index]][1,1]
      imp_bp2_s1=random_coef_impute_id2[[s_index]][1,2]
      imp_bp2_j0=random_coef_impute_id2[[s_index]][tempwhich,1]
      imp_bp2_j1=random_coef_impute_id2[[s_index]][tempwhich,2]
      dev=dev+coxcoefvalues[2]*(abs(imp_bp2_s0-imp_bp2_j0)+abs(imp_bp2_s1-imp_bp2_j1))
    }
    if (dim_bp>=3) {
      imp_bp3_s0=random_coef_impute_id3[[s_index]][1,1]
      imp_bp3_s1=random_coef_impute_id3[[s_index]][1,2]
      imp_bp3_j0=random_coef_impute_id3[[s_index]][tempwhich,1]
      imp_bp3_j1=random_coef_impute_id3[[s_index]][tempwhich,2]
      dev=dev+coxcoefvalues[3]*(abs(imp_bp3_s0-imp_bp3_j0)+abs(imp_bp3_s1-imp_bp3_j1))
    }
    if (dim_bp>=4) {
      imp_bp4_s0=random_coef_impute_id4[[s_index]][1,1]
      imp_bp4_s1=random_coef_impute_id4[[s_index]][1,2]
      imp_bp4_j0=random_coef_impute_id4[[s_index]][tempwhich,1]
      imp_bp4_j1=random_coef_impute_id4[[s_index]][tempwhich,2]
      dev=dev+coxcoefvalues[4]*(abs(imp_bp4_s0-imp_bp4_j0)+abs(imp_bp4_s1-imp_bp4_j1))
    }
    if (dim_bp>=5) {
      imp_bp5_s0=random_coef_impute_id5[[s_index]][1,1]
      imp_bp5_s1=random_coef_impute_id5[[s_index]][1,2]
      imp_bp5_j0=random_coef_impute_id5[[s_index]][tempwhich,1]
      imp_bp5_j1=random_coef_impute_id5[[s_index]][tempwhich,2]
      dev=dev+coxcoefvalues[5]*(abs(imp_bp5_s0-imp_bp5_j0)+abs(imp_bp5_s1-imp_bp5_j1))
    }

    orderdev=sort(dev)
    if (length(orderdev) < num_similar) {cutoff=orderdev[length(orderdev)]}
    if (length(orderdev) >= num_similar) {cutoff=orderdev[num_similar]}

    if (dim_bp==1) {
      idmat=idtemp[coxcoefvalues[1]*(abs(imp_bp_s0-imp_bp_j0)+abs(imp_bp_s1-imp_bp_j1))+dev_obs <= cutoff]
    }
    if (dim_bp==2) {
      idmat=idtemp[coxcoefvalues[1]*(abs(imp_bp_s0-imp_bp_j0)+abs(imp_bp_s1-imp_bp_j1))+dev_obs+
                   coxcoefvalues[2]*(abs(imp_bp2_s0-imp_bp2_j0)+abs(imp_bp2_s1-imp_bp2_j1)) <= cutoff]
    }
    if (dim_bp==3) {
      idmat=idtemp[coxcoefvalues[1]*(abs(imp_bp_s0-imp_bp_j0)+abs(imp_bp_s1-imp_bp_j1))+dev_obs+
                   coxcoefvalues[2]*(abs(imp_bp2_s0-imp_bp2_j0)+abs(imp_bp2_s1-imp_bp2_j1))+
                   coxcoefvalues[3]*(abs(imp_bp3_s0-imp_bp3_j0)+abs(imp_bp3_s1-imp_bp3_j1)) <= cutoff]
    }
    if (dim_bp==4) {
      idmat=idtemp[coxcoefvalues[1]*(abs(imp_bp_s0-imp_bp_j0)+abs(imp_bp_s1-imp_bp_j1))+dev_obs+
                   coxcoefvalues[2]*(abs(imp_bp2_s0-imp_bp2_j0)+abs(imp_bp2_s1-imp_bp2_j1))+
                   coxcoefvalues[3]*(abs(imp_bp3_s0-imp_bp3_j0)+abs(imp_bp3_s1-imp_bp3_j1))+
                   coxcoefvalues[4]*(abs(imp_bp4_s0-imp_bp4_j0)+abs(imp_bp4_s1-imp_bp4_j1)) <= cutoff]
    }
    if (dim_bp==5) {
      idmat=idtemp[coxcoefvalues[1]*(abs(imp_bp_s0-imp_bp_j0)+abs(imp_bp_s1-imp_bp_j1))+dev_obs+
                   coxcoefvalues[2]*(abs(imp_bp2_s0-imp_bp2_j0)+abs(imp_bp2_s1-imp_bp2_j1))+
                   coxcoefvalues[3]*(abs(imp_bp3_s0-imp_bp3_j0)+abs(imp_bp3_s1-imp_bp3_j1))+
                   coxcoefvalues[4]*(abs(imp_bp4_s0-imp_bp4_j0)+abs(imp_bp4_s1-imp_bp4_j1))+
                   coxcoefvalues[5]*(abs(imp_bp5_s0-imp_bp5_j0)+abs(imp_bp5_s1-imp_bp5_j1)) <= cutoff]
    }

#      cat('--------------------------------------------','\n')
#      cat('ICOX2 Matched past subjects =',idmat,'\n')
#      cat('--------------------------------------------','\n')
    idmatch2[[j]]=idmat

    j=j+1
  }

#--------------------------------------------------------------------------
#
# Current Patient
# Replace ALL values of Biomarker Process with fitted value from model
#
  i=1
  while (i <= n) {
#    cat('-----------------------','\n')
#    cat('i=',i,'\n')
    id=dataset_impute2$id[i]
#    cat('id=',id,'\n')
#    cat('start[i]=',dataset_impute[[k]]$start[i],'\n')
#    cat('round(start[i])=',round(dataset_impute[[k]]$start[i],digits=digits_imp_inc),'\n')
#    cat('imp_times_bp=',imp_times_bp[40:45],'\n')
    j=which(round(imp_times_bp,digits=digits_imp_inc)==round(dataset_impute2$start[i],digits=digits_imp_inc))
#    cat('j=',j,'\n')
#    cat('time=',imp_times_bp[j],'\n')
    if (id<=nperson_future) {
#      if (j > s_index) {j=s_index}
#      cat('coef_impute[[j]]=',coef_impute[[j]],'\n')
#      cat('dataset_impute[[k]]$baseage[i]=',dataset_impute[[k]]$baseage[i],'\n')
#      cat('fit=',fit,'\n')
#-------------------------------------------------
    if (j <= s_index) {
      fit1=coef_impute1[[j]][1]
      fit1=fit1+coef_impute1[[j]][2]*dataset_impute2$baseage[i]
      fit1=fit1+coef_impute1[[j]][3]*dataset_impute2$start[i]
      fit1=fit1+coef_impute1[[j]][4]*dataset_impute2$start[i]**2
      fit1=fit1+coef_impute1[[j]][5]*dataset_impute2$start[i]**3

      if (dim_bp>=2) {
        fit2=coef_impute2[[j]][1]
        fit2=fit2+coef_impute2[[j]][2]*dataset_impute2$baseage[i]
        fit2=fit2+coef_impute2[[j]][3]*dataset_impute2$start[i]
        fit2=fit2+coef_impute2[[j]][4]*dataset_impute2$start[i]**2
        fit2=fit2+coef_impute2[[j]][5]*dataset_impute2$start[i]**3
      }

      if (dim_bp>=3) {
        fit3=coef_impute3[[j]][1]
        fit3=fit3+coef_impute3[[j]][2]*dataset_impute2$baseage[i]
        fit3=fit3+coef_impute3[[j]][3]*dataset_impute2$start[i]
        fit3=fit3+coef_impute3[[j]][4]*dataset_impute2$start[i]**2
        fit3=fit3+coef_impute3[[j]][5]*dataset_impute2$start[i]**3
      }

      if (dim_bp>=4) {
        fit4=coef_impute3[[j]][1]
        fit4=fit4+coef_impute4[[j]][2]*dataset_impute2$baseage[i]
        fit4=fit4+coef_impute4[[j]][3]*dataset_impute2$start[i]
        fit4=fit4+coef_impute4[[j]][4]*dataset_impute2$start[i]**2
        fit4=fit4+coef_impute4[[j]][5]*dataset_impute2$start[i]**3
      }

      if (dim_bp>=5) {
        fit5=coef_impute3[[j]][1]
        fit5=fit5+coef_impute5[[j]][2]*dataset_impute2$baseage[i]
        fit5=fit5+coef_impute5[[j]][3]*dataset_impute2$start[i]
        fit5=fit5+coef_impute5[[j]][4]*dataset_impute2$start[i]**2
        fit5=fit5+coef_impute5[[j]][5]*dataset_impute2$start[i]**3
      }

#      cat('after fixed effects fit1=',fit1,'\n')
#-------------------------------------------------
# TEST AIDS DATA 
#      fit=coef_impute[[j]][1]
#      if (dataset_impute[[k]]$baseage[i]==1) {fit=fit+coef_impute[[j]][3]*dataset_impute[[k]]$start[i]}
#-------------------------------------------------
      id_num=which(select_id[[j]]==dataset_impute2$id[i])
#      cat('idnum=',id_num,'\n')

      fit1=fit1+random_coef_impute_id1[[j]][id_num,1]
      fit1=fit1+random_coef_impute_id1[[j]][id_num,2]*dataset_impute2$start[i]
      dataset_impute2$bp_curr1[i]=fit1
      
      if (dim_bp >=2) {
        fit2=fit2+random_coef_impute_id2[[j]][id_num,1]
        fit2=fit2+random_coef_impute_id2[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr2[i]=fit2
      }

      if (dim_bp >=3) {
        fit3=fit3+random_coef_impute_id3[[j]][id_num,1]
        fit3=fit3+random_coef_impute_id3[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr3[i]=fit3
      }

      if (dim_bp >=4) {
        fit4=fit4+random_coef_impute_id4[[j]][id_num,1]
        fit4=fit4+random_coef_impute_id4[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr4[i]=fit4
      }

      if (dim_bp >=5) {
        fit5=fit5+random_coef_impute_id5[[j]][id_num,1]
        fit5=fit5+random_coef_impute_id5[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr5[i]=fit5
      }
    }

    if (j>s_index) {
      id_num=which(select_id[[j]]==dataset_impute2$id[i])
      if (id==1) {
        id_num=which(select_id[[j]] %in% idmatch2[[j]])
      }
      fit1=coef_impute1[[j]][1]
     
      if (id==1) {
        fit1=fit1+coef_impute1[[j]][2]*mean(baseage_id[idmatch2[[j]]])
      } else {
        fit1=fit1+coef_impute1[[j]][2]*dataset_impute2$baseage[i]
      }
      fit1=fit1+coef_impute1[[j]][3]*dataset_impute2$start[i]
      fit1=fit1+coef_impute1[[j]][4]*dataset_impute2$start[i]**2
      fit1=fit1+coef_impute1[[j]][5]*dataset_impute2$start[i]**3

      if (dim_bp>=2) {
        fit2=coef_impute2[[j]][1]
        fit2=fit2+coef_impute2[[j]][2]*dataset_impute2$baseage[i]
        fit2=fit2+coef_impute2[[j]][3]*dataset_impute2$start[i]
        fit2=fit2+coef_impute2[[j]][4]*dataset_impute2$start[i]**2
        fit2=fit2+coef_impute2[[j]][5]*dataset_impute2$start[i]**3
      }

      if (dim_bp>=3) {
        fit3=coef_impute3[[j]][1]
        fit3=fit3+coef_impute3[[j]][2]*dataset_impute2$baseage[i]
        fit3=fit3+coef_impute3[[j]][3]*dataset_impute2$start[i]
        fit3=fit3+coef_impute3[[j]][4]*dataset_impute2$start[i]**2
        fit3=fit3+coef_impute3[[j]][5]*dataset_impute2$start[i]**3
      }

      if (dim_bp>=4) {
        fit4=coef_impute3[[j]][1]
        fit4=fit4+coef_impute4[[j]][2]*dataset_impute2$baseage[i]
        fit4=fit4+coef_impute4[[j]][3]*dataset_impute2$start[i]
        fit4=fit4+coef_impute4[[j]][4]*dataset_impute2$start[i]**2
        fit4=fit4+coef_impute4[[j]][5]*dataset_impute2$start[i]**3
      }

      if (dim_bp>=5) {
        fit5=coef_impute3[[j]][1]
        fit5=fit5+coef_impute5[[j]][2]*dataset_impute2$baseage[i]
        fit5=fit5+coef_impute5[[j]][3]*dataset_impute2$start[i]
        fit5=fit5+coef_impute5[[j]][4]*dataset_impute2$start[i]**2
        fit5=fit5+coef_impute5[[j]][5]*dataset_impute2$start[i]**3
      }

#      cat('after fixed effects fit1=',fit1,'\n')
#-------------------------------------------------
# TEST AIDS DATA 
#      fit=coef_impute[[j]][1]
#      if (dataset_impute[[k]]$baseage[i]==1) {fit=fit+coef_impute[[j]][3]*dataset_impute[[k]]$start[i]}
#-------------------------------------------------
      fit1=fit1+mean(random_coef_impute_id1[[j]][id_num,1])
      fit1=fit1+mean(random_coef_impute_id1[[j]][id_num,2]*dataset_impute2$start[i])
      dataset_impute2$bp_curr1[i]=fit1
      
      if (dim_bp >=2) {
        fit2=fit2+random_coef_impute_id2[[j]][id_num,1]
        fit2=fit2+random_coef_impute_id2[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr2[i]=fit2
      }

      if (dim_bp >=3) {
        fit3=fit3+random_coef_impute_id3[[j]][id_num,1]
        fit3=fit3+random_coef_impute_id3[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr3[i]=fit3
      }

      if (dim_bp >=4) {
        fit4=fit4+random_coef_impute_id4[[j]][id_num,1]
        fit4=fit4+random_coef_impute_id4[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr4[i]=fit4
      }

      if (dim_bp >=5) {
        fit5=fit5+random_coef_impute_id5[[j]][id_num,1]
        fit5=fit5+random_coef_impute_id5[[j]][id_num,2]*dataset_impute2$start[i]
        dataset_impute2$bp_curr5[i]=fit5
      }

    }
    }

    i=i+1
  }
#--------------------------------------------------------------------------
#
# Remove observations past s+tau in FUTURE_DATASET
#
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

  future_dataset_impute=dataset_impute2[dataset_impute2$id<=1 & round(dataset_impute2$start,digits=5) < s+tau,]
  n_future=length(future_dataset_impute$id)
  nperson_future=1

#  cat('Final Future Dataset with imputed BP process =','\n')
#  print(future_dataset_impute[[k]][1:20,])

#--------------------------------------------------------------------------
#
# FORTRAN VERSION OF LOOP TO CALCULATE CONDITIONAL SURVIVAL PROBS
#
#--------------------------------------------------------------------------
  cond_surv_pred_st=rep(0,times=nperson_future)
  cat('Call predictSurvProb','\n')
#  future_surv=predictSurvProb(efit[[k]],newdata=future_dataset_impute_long,times=predtimes)
  future_surv=predictSurvProb(efit[[k]],newdata=future_dataset_impute,times=predtimes)
  cat('Return from predictSurvProb','\n')
#  if (!all(is.finite(future_surv))) {
#    cat('survival to predtimes for subjects in future_dataset according to EMPIRICAL COX MODEL=','\n')
#    print(future_surv[1:20,])
#  }
#  future_surv[future_surv > 1]=NA
#  future_surv[future_surv < 0]=NA
  npredtimes=length(predtimes)
  future_surv_short=future_surv
#  future_surv_short=future_surv[1:n_future,]
#  if (all(is.finite(future_surv_short))) {
#    cat('All future_surv_short are finite','\n')
#  } else {
#    cat('survival to predtimes for subjects in future_dataset according to EMPIRICAL COX MODEL=','\n')
#    print(future_surv_short)
#  }


#test=future_dataset_impute_long
#test$bp.1[1]=999
#test_surv=predictSurvProb(efit[[k]],newdata=test,times=predtimes)
#test_surv_short=test_surv[1:n_future,]

#cat('survival to predtimes for subjects in test dataset according to EMPIRICAL COX MODEL=','\n')
#print(test_surv_short)


#  idvec=future_dataset_impute_long$id
  idvec=future_dataset_impute$id

#  cat('survival to predtimes for subjects in future_dataset according to EMPIRICAL COX MODEL=','\n')
#  print(future_surv_short)

  if (all(is.finite(future_surv_short))) {
    cat('Call condprob','\n')
    z=condprob(nperson_future,n_future,npredtimes,index_s,future_surv_short,idvec,cond_surv_pred_st)
    cat('return from condprob','\n')

    cond_surv_pred_st=z[[7]]

    cat('FORTRAN: Cond. survival to S+TAU given survival to S for subjects in future_dataset according to IMP2 COX MODEL=','\n')
    print(cond_surv_pred_st)

    icox2_cond_surv_pred_imp=mean(cond_surv_pred_st)

    cat('FORTRAN: Avg. Cond. survival to S+TAU given survival to S for subjects in future_dataset according to IMP2 COX MODEL=','\n')
    print(icox2_cond_surv_pred_imp)
  } else {
    cat('IMP2 COX MODEL DOES NOT GIVE PROB ESTIMATE','\n')
    icox2_cond_surv_pred_imp=NA
  }

#  k=k+1
#}
#--------------------------------------------------------------------------
# END LOOP OVER IMPUTATIONS
#--------------------------------------------------------------------------

  cat('---------------------------------------------------','\n')
  cat('---------------------------------------------------','\n')
  if (any(is.finite(icox2_cond_surv_pred_imp))) {
   icox2_cond_surv_pred_boot[iboot]=mean(icox2_cond_surv_pred_imp,na.rm=TRUE)
    cat('FINAL Estimate Cond. survival to S+TAU given survival to S for subjects in future_dataset according to IMP2 COX MODEL=','\n')
    print(icox2_cond_surv_pred_boot[iboot])
  } else {
    cat('FINAL IMP2 COX MODEL DOES NOT GIVE PROB ESTIMATE','\n')
    icox2_cond_surv_pred_boot[iboot]=NA
  }
  cat('--------------------------------------------------------','\n')


  }
#---------------------------------------------------------------------------------------------------
# END ICOX2: METHOD THAT USES BLUP FROM NUM-SIMILAR PAST SUBJECTS MOST SIMILAR BY S-OBSERVED VALUES
#---------------------------------------------------------------------------------------------------


  iboot=iboot+1
}
#cat('ICOX1 bootstrap conditional probs=',icox1_cond_surv_pred_boot,'\n')
icox1_boot_sd[eval_id_num]=sd(icox1_cond_surv_pred_boot,na.rm=TRUE)
#cat('ICOX2 bootstrap conditional probs=',icox2_cond_surv_pred_boot,'\n')
icox2_boot_sd[eval_id_num]=sd(icox2_cond_surv_pred_boot,na.rm=TRUE)
cat('--------------------------------------------------------','\n')
cat('--------------------------------------------------------','\n')
cat('SD over bootstrap samples of Cond. survival to S+TAU given survival to S for subjects in bootstrap future_dataset according to IMP1 COX MODEL=','\n')
print(icox1_boot_sd[eval_id_num])
cat('--------------------------------------------------------','\n')
cat('SD over bootstrap samples of Cond. survival to S+TAU given survival to S for subjects in bootstrap future_dataset according to IMP2 COX MODEL=','\n')
print(icox2_boot_sd[eval_id_num])
cat('--------------------------------------------------------','\n')
cat('--------------------------------------------------------','\n')
#-------------------------------------------------
# END: BOOTSTRAP
#-------------------------------------------------

#-------------------------------------------------
#-------------------------------------------------
}
#-------------------------------------------------
# END: IMP COX MODEL
# END: icox1_ind=1 or icox2_ind=1 or icox3_ind=1
#-------------------------------------------------













#-------------------------------------------------
#-------------------------------------------------
if (jm_ind==1) {
#-------------------------------------------------
cat('---------------------------------------------------------------------','\n')
cat('START   JM','\n')
cat('---------------------------------------------------------------------','\n')
#-------------------------------------------------
# JOINT MODEL EXAMPLE Rizopoulos 2010, J of Stat Software
#
#cat('---------------------------------------------------------------------','\n')
#cat('Start JM of Rizoupoulos 2010 for CD4 dataset=','\n')
#cat('---------------------------------------------------------------------','\n')
#cat('aids dataset=','\n')
#print(aids[1:20,])
#cat('---------------------------------------------------------------------','\n')
#cat('Cox model of CD4 dataset with drug and sqrt(cd4) as covariates=','\n')
#td.Cox=coxph(Surv(start,stop,event)~drug+sqrt(CD4),data=aids)
#print(summary(td.Cox))
#cat('---------------------------------------------------------------------','\n')
#fitLME=lme(sqrt(CD4)~obstime+obstime:drug,random=~obstime | patient, data=aids)
#cat('Longitudinal model of CD4 dataset','\n')
#print(summary(fitLME))
#cat('---------------------------------------------------------------------','\n')
#fitSURV=coxph(Surv(Time,death)~drug,data=aids.id,x=TRUE)
#cat('Cox model of CD4 dataset with only drug as covariate=','\n')
#print(summary(fitSURV))
#cat('---------------------------------------------------------------------','\n')
#fit.JM=jointModel(fitLME,fitSURV,timeVar="obstime",method="piecewise-PH-GH")
#print(summary(fit.JM))
#cat('---------------------------------------------------------------------','\n')
#set.seed(123)
#NewData=aids[aids$patient %in% c("7","15","117","303"),]
#NewData$s=rep(s,times=length(NewData$patient))
#cat('---------------------------------------------------------------------','\n')
#cat('NewData dataset=','\n')
#print(NewData)
#cat('---------------------------------------------------------------------','\n')
#predSurv=survfitJM(fit.JM,newdata=NewData,idVar="patient",last.time="Time")
#print(predSurv)
#predSurv=survfitJM(fit.JM,newdata=NewData,idVar="patient",last.time="s",survTimes=s+tau)
#cat('---------------------------------------------------------------------','\n')
#cat('Set s & s+Tau','\n')
#cat('---------------------------------------------------------------------','\n')
#print(predSurv)
#cat('---------------------------------------------------------------------','\n')
#print(predSurv[[1]][[1]])
#print(predSurv[[1]][[1]][2])
#cat('---------------------------------------------------------------------','\n')
#cat('END OF EXAMPLE FROM RIZ (2010)','\n')
#cat('---------------------------------------------------------------------','\n')
#
#-------------------------------------------------

#-------------------------------------------------
# JOINT MODEL
#
#-------------------------------------------------
converge=1
#-------------------------------------------------
#dataset1=dataset[lastobs(~id,data=dataset),]
#dataset1=dataset
#dataset1=subset(dataset1,select=-c(start))
#cat('dataset1=','\n')
#print(dataset1)

#cat('dataset_nonmiss=','\n')
#print(dataset_nonmiss[1:20,])

#data_bp=data_bp[data_bp$time <= dataset1$stop[data_bp$id],]
#cat('data_bp=','\n')
#print(data_bp[1:50,])

#dataset_jm=data.frame(id=dataset$id,start=dataset$start,stop=dataset$stop,status=dataset$status,bp_curr=dataset$bp_curr,baseage=dataset$baseage,
#                         lastobs=dataset$lastobs)
dataset_nonmiss_jm=dataset_nonmiss_save
future_dataset_nonmiss_jm=future_dataset_nonmiss
dataset_nonmiss_jm$id=dataset_nonmiss_jm$id+1
#------------------------------------------------------------------
#  TRY INCLUDING FUTURE PERSON IN LONGITUDINAL DATASET TO FIT JM
#------------------------------------------------------------------
dataset_nonmiss_jm=smartbind(future_dataset_nonmiss_jm,dataset_nonmiss_jm)

dataset_jm1=dataset_save[dataset_save$lastobs==1,c("id","baseage","status","stop","bp_curr1")]
dataset_jm1$id=dataset_jm1$id+1
future_dataset_jm1=future_dataset[future_dataset$lastobs==1,c("id","baseage","status","stop","bp_curr1")]
#------------------------------------------------------------------
#  TRY INCLUDING FUTURE PERSON IN SURVIVAL DATASET TO FIT JM1
#------------------------------------------------------------------
dataset_jm1=smartbind(future_dataset_jm1,dataset_jm1)


#future_dataset_jm=dataset_nonmiss[dataset_nonmiss$id==1 & dataset_nonmiss$obstime <= s+.0001,]
#dataset_jm=dataset_nonmiss[dataset_nonmiss$id != 1,]

dataset_jm=dataset_nonmiss_jm
future_dataset_jm=future_dataset_nonmiss_jm
dataset_jm$obstime2=dataset_jm$obstime**2
dataset_jm$obstime3=dataset_jm$obstime**3
future_dataset_jm$obstime2=future_dataset_jm$obstime**2
future_dataset_jm$obstime3=future_dataset_jm$obstime**3

#cat('dataset_jm for ID=15=','\n')
#print(dataset_jm[dataset_jm$id==15,])
#cat('dataset_jm for ID=16=','\n')
#print(dataset_jm[dataset_jm$id==16,])
#cat('dataset_jm$id','\n')
#print(dataset_jm$id)
#cat('# of unique ids=','\n')
#print(length(unique(dataset_jm$id)))

#cat('future_dataset_jm=','\n')
#print(future_dataset_jm)

#cat('dataset_jm=','\n')
#print(dataset_jm[1:20,])

#cat('dataset_jm1=','\n')
#print(dataset_jm1[1:20,])

fitlme=try(lme(bp_curr1~obstime+baseage,random=~obstime | id,data=dataset_jm, control=lmeControl(opt="optim"),na.action=na.exclude))
#fitlme=try(lme(bp_curr1~obstime+obstime2+obstime3+baseage,random=~obstime | id,data=dataset_jm, control=lmeControl(opt="optim"),na.action=na.exclude))
#fitlme=try(lme(bp_curr~start+baseage,random=~1 | id,data=dataset_jm, control=lmeControl(opt="optim")))
#fitlme=lmer(bp_curr~start+baseage + (start | id),data=dataset_jm)
if (class(fitlme)=="try-error") {
  converge=0
  cat('LME DOES NOT CONVERGE','\n')
} else {
  cat('LME MODEL OF BP=','\n')
  print(fitlme)

#  dataset_jm1=dataset_save[dataset_save$lastobs==1,c("id","baseage","status","stop","bp_curr")]
#  dataset_jm1$id=dataset_jm1$id+1

#  cat('dataset_jm1=','\n')
#  print(dataset_jm1[1:20,])
#  print(summary(dataset_jm1$stop))

  fitcox=coxph(Surv(stop,status)~baseage,data=dataset_jm1,x=TRUE)
  cat('COX MODEL with only baseage=','\n')
  print(fitcox)

  control_jm=vector("list",14)
  control_jm[[2]]=100
  fitjm=try(jointModel(fitlme,fitcox,timeVar="obstime",method="piecewise-PH-aGH"))
#  fitjm=try(jointModelBayes(fitlme,fitcox,timeVar="obstime",param="shared-RE",control=control_jm))
  if (class(fitjm)=="try-error") {
    converge=0
    cat('JOINT MODEL DOES NOT CONVERGE','\n')
  } else {

    cat('JOINT MODEL FINISHED SUCCESSFULLY','\n')
#    print(summary(fitjm))
#    cat('coefficients=','\n')
#    print(fitjm$coefficients)
#    cat('Hessian=','\n')
#    print(fitjm$Hessian)

    n_future_jm=length(future_dataset_jm$id)
#    cat('Future Dataset_jm=','\n')
#    print(future_dataset_jm)

#    cat('Dataset=','\n')
#    print(dataset[1:50,])

#    temp=dataset_jm
#    temp$id=temp$id+nperson_future
#    temp[1:n_future_jm,]=future_dataset_jm
#    jm_future_dataset=temp

#    jm_future_dataset$bp_curr=jm_future_dataset$bp_curr*(jm_future_dataset$bp_curr != -999.)
#    jm_future_dataset$time=jm_future_dataset$start
#    jm_future_dataset=subset(jm_future_dataset,select=-c(start))

#    jm_future_dataset$status=rep(0,times=length(jm_future_dataset$id))
#    jm_future_dataset$svec=rep(s,times=length(jm_future_dataset$id))


#    jm_future_dataset$bp_curr[jm_future_dataset$id > 1]=0

#    jm_future_dataset=jm_future_dataset[jm_future_dataset$id==1,]
#    jm_future_dataset$stop[jm_future_dataset$stop>s]=s

    cat('JM Future Dataset=','\n')
    print(future_dataset_jm[1:n_future_jm,])
    future_dataset_jm$s=rep(s,times=length(future_dataset_jm$id))

#
#  TEST OF JM
#
#    future_dataset_jm=future_dataset_impute[[1]]
#    future_dataset_jm$obstime=future_dataset_jm$start
#    future_dataset_jm$obstime1=future_dataset_jm$start
#    future_dataset_jm$obstime2=future_dataset_jm$start**2
#    future_dataset_jm$obstime3=future_dataset_jm$start**3
#    n_future_jm=length(future_dataset_jm$id)
#    future_dataset_jm$lastobs=rep(0,times=n_future)
#    future_dataset_jm=future_dataset_jm[,c("id","baseage","bp_curr","obstime","lastobs","obstime1","obstime2","obstime3")]
#    cat('Imputed Future Dataset=','\n')
#    print(future_dataset_jm[1:n_future_jm,])

#    cat('s+tau=',s+tau,'\n')

    jmsurv=try(survfitJM(fitjm,newdata=future_dataset_jm,idVar="id",last.time="s",survTimes=c(s+tau)))
#    jmsurv=try(survfitJM(fitjm,newdata=future_dataset_jm,idVar="id",last.time=s,survTimes=c(tau+.01)))
    if (class(jmsurv)=="try-error") {
      converge=0
      cat('JOINT MODEL DOES NOT LEAD TO SUCCESSFUL ESTIMATE OF SURVIVAL','\n')
    } else {
      cat('JM Future Predicted Survival=','\n')
      print(jmsurv[[1]][[1]])
#      print(jmsurv[[1]][[1]][2])
      i=1
      while (i <= nperson_future) {
        cond_surv_pred_st[i]=jmsurv[[1]][[i]][2]
        i=i+1
      }
#      print(cond_surv_pred_st)
      jm_cond_surv_pred[iconv]=mean(cond_surv_pred_st)
      cat('Avg. Cond. survival to S+TAU given survival to S for subjects in future_dataset according to JOINT MODEL=','\n')
      print(jm_cond_surv_pred[iconv])
    }
  }
}
#-------------------------------------------------
#-------------------------------------------------
}
#-------------------------------------------------
# END: JOINT MODEL
#-------------------------------------------------









#-------------------------------------------------
# SAVE RESULTS
#-------------------------------------------------

data_summary=data.frame(id=dataset_eval_1$id, stop=dataset_eval_1$stop,
                        status=dataset_eval_1$status,
                        icox1=icox1_cond_surv_pred,
                        icox2=icox2_cond_surv_pred,
                        lvcf=lvcf_cond_surv_pred)
data_summary1=data_summary[!is.na(data_summary$icox1),]
cat('DATA_SUMMARY of predictions','\n')
print(data_summary1)
print(summary(data_summary1))
save(data_summary1,file="/data/troendlj/accord/accord.pred.ms41.Rda")




warnings()

write_results=function() {
 cat('Predict risk of event based on past history','\n'
 ,'James F. Troendle, January 2021','\n'
 ,'-------------------------------------------------------------------------','\n'
 ,'-------------------------------------------------------------------------','\n'
 ,'Number of Evaluable subjects=     ',n_eval,'\n'
 ,'-------------------------------------------------------------------------','\n'
 ,'Conditioning time:                     S= ',s,'\n'
 ,'Prediction of survival to S+TAU:     TAU= ',tau,'\n'
 ,'-------------------------------------------------------------------------','\n'
 ,'-------------------------------------------------------------------------','\n'
 ,'# of Bootstraps used for estimating SD  =',boot_num,'\n'
 ,'-------------------------------------------------------------------------','\n'
 ,'-------------------------------------------------------------------------','\n'
 ,'Parameters related to ICOX Models','\n'
 ,'---------------------------------','\n'
 ,'Increment for BP updates=                ',imp_inc,'\n'
 ,'Times when BP is updated=                ',imp_times_bp,'\n'
 ,'Survival data left censoring time=       ',left_cens_time,'\n'
 ,'Survival data right censoring time=      ',cens_time,'\n'
 ,'-------------------------------------------------------------------------','\n'
 ,'------------------------------------------------------------------------------------------------------------','\n')
}





write_results()
