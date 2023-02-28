#------------------------------------------------------------
#SET PARAMETERS FOR PROGRAM
#------------------------------------------------------------
#Initialize seed for random number generation
seed=99803332 
seedstart=seed
set.seed(seed)
#


#------------------------------------------------------------
# LOAD R LIBRARIES
#------------------------------------------------------------
library(survival)
library(stats)
library(gtools)
library(mvtnorm)
library(tibble)
library(lme4)
library(dplyr)
library(plyr)
library(pryr)
library(mvtnorm)
library(pec)
library(nlme)
library(merTools)
library(JM)
#library(JMbayes)
library(doBy)
#------------------------------------------------------------------------------------
#  DATA LOCATION
#------------------------------------------------------------------------------------
infile_surv="dataset_survival.Rda"
infile_bio1="dataset_biomarker_HBA1C.Rda"
infile_bio2="dataset_biomarker_LDL.Rda"
infile_bio3="dataset_biomarker_SBP.Rda"

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
# S = conditioning time, given survival to time S
# TAU = prediction for survival to time S+TAU
# START_EVAL_ID_NUM
# END_EVAL_ID_NUM
#------------------------------------------------------------
s=1.0
tau=2.0

#---------------------------------------------------
# PREDICTION FOR SUBSET OF EVALUABLE SUBJECTS
#------------------------------------------------------------
start_eval_id_num=1
end_eval_id_num=2

#------------------------------------------------------------
# BIOMARKER DIMENSION (MAX 5)
#------------------------------------------------------------
dim_bp=3
inc=0.01

#------------------------------------------------------------
# INDICATOR OF METHODS ACTIVE
#------------------------------------------------------------
lmm1_ind=1
lmm2_ind=0

#------------------------------------------------------------
# SPLINE PARAMETERS
#------------------------------------------------------------
nspline_max=c(7,3,7)
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
maxfollow=7
window_width=0.00001
num_design_times_bp=length(design_times_bp)

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# BOOTSTRAP USED FOR ESTIMATING SD OF PREDICTION
#------------------------------------------------------------------------------------
boot_num=0
#
#------------------------------------------------------------------------------------
#------------------------------------
#options(warn=1)
options(warn=-1)
z975=qnorm(.975,mean=0,sd=1)

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

# cat('Before sequential IDs Longitudinal Dataset','\n')
# cat('Dataset_nonmiss =','\n')
# print(dataset_nonmiss[1:50,])
# cat('--------------------------------------------','\n')
# cat('Before sequential IDs Survival Dataset','\n')
# print(summary(dataset))
# print(dataset[1:30,])

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
#  CREATE SEQUENTIAL ID FROM 1 TO NPERSON
#
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------
nid=length(unique(dataset$id))
idvec=rep(0,times=nid)
i=1
idcount=0
idprev=0
while (i <=nid) {
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

#-----------------------------------------------------------------------------
# Define time known to have event
# Define time known to survive to
#-----------------------------------------------------------------------------
dataset_1=dataset
nonsurvive_time=dataset_1$stop*(dataset_1$status==1)+(s+tau+inc)*(dataset_1$status!=1)
nonsurvive_time_save=nonsurvive_time
survive_time=dataset_1$stop
#-----------------------------------------------------------------------------
# cat('---------------------------------------------','\n')
# cat('After sequential IDs: dataset=','\n')
# print(summary(dataset))
# print(dataset[1:30,])
# cat('After sequential IDs: # of IDs in dataset =',length(unique(dataset$id)),'\n')
# cat('---------------------------------------------','\n')

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

# cat('---------------------------------------------','\n')
# cat('After sequential IDs: longitudinal dataset=','\n')
# print(summary(dataset_nonmiss))
# cat('---------------------------------------------','\n')
# cat('Dataset_nonmiss for ID=20 or ID=22 =','\n')
# print(dataset_nonmiss[dataset_nonmiss$id==20 | dataset_nonmiss$id==22,])
# cat('--------------------------------------------','\n')

#--------------------------------------------------------------------------
#  CALCULATE NUMBER OF OBS PRIOR TO S FOR EACH BIOMARKER FOR EACH SUBJECT  
#--------------------------------------------------------------------------
n_biomarker_s=rep(0,times=dim_bp*nperson)
n_biomarker_s_tau=rep(0, times=dim_bp*nperson)
dim(n_biomarker_s)=c(nperson,dim_bp)
dim(n_biomarker_s_tau)=c(nperson,dim_bp)
i=1
while (i <=n_nonmiss) {
  if (dataset_nonmiss$obstime[i] <= s) {
    if (dataset_nonmiss$bptype[i]=='TYPE1') {
      n_biomarker_s[dataset_nonmiss$id[i],1]=n_biomarker_s[dataset_nonmiss$id[i],1]+1
    } else {
      if (dataset_nonmiss$bptype[i]=='TYPE2') {
        n_biomarker_s[dataset_nonmiss$id[i],2]=n_biomarker_s[dataset_nonmiss$id[i],2]+1
      } else {
        if (dataset_nonmiss$bptype[i]=='TYPE3') {
          n_biomarker_s[dataset_nonmiss$id[i],3]=n_biomarker_s[dataset_nonmiss$id[i],3]+1
        } else {
          if (dataset_nonmiss$bptype[i]=='TYPE4') {
            n_biomarker_s[dataset_nonmiss$id[i],4]=n_biomarker_s[dataset_nonmiss$id[i],4]+1
          } else {
            if (dataset_nonmiss$bptype[i]=='TYPE5') {
              n_biomarker_s[dataset_nonmiss$id[i],5]=n_biomarker_s[dataset_nonmiss$id[i],5]+1
            }
          }
        }
      }
    }
  }
  if (abs(dataset_nonmiss$obstime[i]-s-tau) <= window_width) {
    if (dataset_nonmiss$bptype[i]=='TYPE1') {
      n_biomarker_s_tau[dataset_nonmiss$id[i],1]=n_biomarker_s_tau[dataset_nonmiss$id[i],1]+1
    } else {
      if (dataset_nonmiss$bptype[i]=='TYPE2') {
        n_biomarker_s_tau[dataset_nonmiss$id[i],2]=n_biomarker_s_tau[dataset_nonmiss$id[i],2]+1
      } else {
        if (dataset_nonmiss$bptype[i]=='TYPE3') {
          n_biomarker_s_tau[dataset_nonmiss$id[i],3]=n_biomarker_s_tau[dataset_nonmiss$id[i],3]+1
        } else {
          if (dataset_nonmiss$bptype[i]=='TYPE4') {
            n_biomarker_s_tau[dataset_nonmiss$id[i],4]=n_biomarker_s_tau[dataset_nonmiss$id[i],4]+1
          } else {
            if (dataset_nonmiss$bptype[i]=='TYPE5') {
              n_biomarker_s_tau[dataset_nonmiss$id[i],5]=n_biomarker_s_tau[dataset_nonmiss$id[i],5]+1
            }
          }
        }
      }
    }
  }
  i=i+1
}
# cat('survive_time[1:10]=',survive_time[1:10],'\n')
# cat('n_biomarker_s[1:10,1]=',n_biomarker_s[1:10,1],'\n')
# IDENTIFY EVALUABLE SUBJECTS:
#    1) Survive to s with >=4 prior obs for each biomarker
#    2) Observe each biomarker at s+tau
#
#--------------------------------------------------------------------------
select_id=unique(dataset_nonmiss$id[survive_time[dataset_nonmiss$id]>s & dataset_nonmiss$bptype=='TYPE1' & dataset_nonmiss$obstime<= s & n_biomarker_s[dataset_nonmiss$id,1] >=4 & n_biomarker_s_tau[dataset_nonmiss$id, 1] == 1])
cat('Biomarker type =TYPE1','\n')
cat('# of Subjects known to survive to s and have >=4 prior obs biomarker =',length(select_id),'\n')
cat('--------------------------------------------','\n')
if (dim_bp >= 2) {
  select_id_2=unique(dataset_nonmiss$id[survive_time[dataset_nonmiss$id]>s & dataset_nonmiss$bptype=='TYPE2' & dataset_nonmiss$obstime<= s & n_biomarker_s[dataset_nonmiss$id,2] >=4 & n_biomarker_s_tau[dataset_nonmiss$id, 2] == 1])
  cat('Biomarker type =TYPE2','\n')
  cat('# of Subjects known to survive to s and have >=4 prior obs biomarker =',length(select_id_2),'\n')
  cat('--------------------------------------------','\n')
  select_id=select_id[select_id %in% select_id_2]
}
if (dim_bp >= 3) {
  select_id_3=unique(dataset_nonmiss$id[survive_time[dataset_nonmiss$id]>s & dataset_nonmiss$bptype=='TYPE3' & dataset_nonmiss$obstime<= s & n_biomarker_s[dataset_nonmiss$id,3] >=4 & n_biomarker_s_tau[dataset_nonmiss$id, 3] == 1])
  cat('Biomarker type =TYPE3','\n')
  cat('# of Subjects known to survive to s and have >=4 prior obs biomarker =',length(select_id_3),'\n')
  cat('--------------------------------------------','\n')
  select_id=select_id[select_id %in% select_id_3]
}
if (dim_bp >= 4) {
  select_id_4=unique(dataset_nonmiss$id[survive_time[dataset_nonmiss$id]>s & dataset_nonmiss$bptype=='TYPE4' & dataset_nonmiss$obstime<= s & n_biomarker_s[dataset_nonmiss$id,4] >=4 & n_biomarker_s_tau[dataset_nonmiss$id, 4] == 1])
  cat('Biomarker type =TYPE4','\n')
  cat('# of Subjects known to survive to s and have >=4 prior obs biomarker =',length(select_id_4),'\n')
  cat('--------------------------------------------','\n')
  select_id=select_id[select_id %in% select_id_4]
}
if (dim_bp >= 5) {
  select_id_5=unique(dataset_nonmiss$id[survive_time[dataset_nonmiss$id]>s & dataset_nonmiss$bptype=='TYPE5' & dataset_nonmiss$obstime<= s & n_biomarker_s[dataset_nonmiss$id,5] >=4 & n_biomarker_s_tau[dataset_nonmiss$id,5 ] == 1])
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
#dataset_eval_1=dataset[dataset$id %in% select_id & dataset$lastobs==1,]
dataset_eval_1=dataset_eval
survive_time_1=survive_time[select_id]
#cat('# of IDs in select_id[[1]]=',length(select_id[[1]]),'\n')
#cat('# of IDs in dataset with ID in select_id[[1]]=',length(unique(dataset$id[dataset$id %in% select_id[[1]]])),'\n')
#cat('# of events in dataset_eval=',sum(dataset_eval$status),'\n')
#cat('# of obs in dataset_eval_1=',length(dataset_eval_1$id),'\n')
# cat('Dataset_eval_1 =','\n')
# print(dataset_eval_1[1:50,])
# cat('--------------------------------------------','\n')
# cat('Dataset_nonmiss for ID=20 or ID=22 =','\n')
# print(dataset_nonmiss[dataset_nonmiss$id==20 | dataset_nonmiss$id==22,])
# cat('--------------------------------------------','\n')

eval_id=select_id
n_eval=length(eval_id)

survive_id=select_id[survive_time_1 > s+tau-window_width]
event_id=select_id[survive_time_1 < s+tau & dataset_eval_1$status==1]

cat('--------------------------------------------','\n')
cat('--------------------------------------------','\n')
#cat('# of events in dataset_eval=',sum(dataset_eval$status),'\n')
#cat('--------------------------------------------','\n')
#cat('summary of survive_time for those in eval_id=','\n')
#print(summary(survive_time[eval_id]))
cat('********************************************','\n')
cat('--------------------------------------------','\n')
cat('# of Evaluable IDs surviving to s with >=4 obs for each biomarker and obs at s+tau=',length(select_id),'\n')
cat('# of Evaluable IDs surviving to s+tau with each biomarker observed =',length(survive_id),'\n')
cat('# of Evaluable IDs with event prior to s+tau =',length(event_id),'\n')
eval_id=survive_id
dataset_eval=dataset[dataset$id %in% eval_id,]
cat('Final # of IDs surviving to s+tau with  >=4 obs for each biomarker and obs at s+tau (EVAL POP)=',length(eval_id),'\n')
cat('--------------------------------------------','\n')
cat('********************************************','\n')
cat('--------------------------------------------','\n')


cat('nonsurvive_time[EVAL POP] =','\n')
print(nonsurvive_time[eval_id][1:40])
cat('--------------------------------------------','\n')


dataset_impute=dataset_eval
n=length(dataset_impute$id)

#--------------------------------------------------------------
#  SUBSET TO OBS PRIOR TO END OF CLINICAL TRIAL & BEFORE EVENT
#--------------------------------------------------------------
dataset_nonmiss_imp=dataset_nonmiss[dataset_nonmiss$obstime <= survive_time[dataset_nonmiss$id],]
dataset_nonmiss_imp_save=dataset_nonmiss_imp
n_nonmiss=length(dataset_nonmiss_imp$id)

#-----------------------------------------------------------------------
#  FOR EACH BPTYPE
#  FOR EACH SUBJECT IN EVAL_ID ADD AN OBS AT TIME=S+TAU WITH BP_CURR=NA
#-----------------------------------------------------------------------
i=start_eval_id_num
while (i <= end_eval_id_num) {
  cat('i=',i,'\n')
  id=eval_id[i]
  na_row=dataset_nonmiss_imp[dataset_nonmiss_imp$id==id,][1,]
  na_row$bp_curr[1]=NA
  na_row$bptype[1]="TYPE1"
  na_row$obstime[1]=s+tau
  dataset_nonmiss_imp<-rbind(dataset_nonmiss_imp, na_row)
  if (dim_bp >=2) {
    na_row$bptype[1]="TYPE2"
    dataset_nonmiss_imp<-rbind(dataset_nonmiss_imp, na_row)
  }
  if (dim_bp >=3) {
    na_row$bptype[1]="TYPE3"
    dataset_nonmiss_imp<-rbind(dataset_nonmiss_imp, na_row)
  }
  if (dim_bp >=4) {
    na_row$bptype[1]="TYPE4"
    dataset_nonmiss_imp<-rbind(dataset_nonmiss_imp, na_row)
  }
  if (dim_bp >=5) {
    na_row$bptype[1]="TYPE5"
    dataset_nonmiss_imp<-rbind(dataset_nonmiss_imp, na_row)
  }
  i=i+1
}
n_nonmiss_imp=length(dataset_nonmiss_imp$id)


# #--------------------------------------------------------------------------
# ##--------------------------------------------------------------------------
# # Spline Creation 
# #--------------------------------------------------------------------------
nspline[1]=max(1,min(floor(maxfollow/spline_width[1]),nspline_max[1]))
if (dim_bp >=2) {nspline[2]=max(1,min(floor(maxfollow/spline_width[2]),nspline_max[2]))}
if (dim_bp >=3) {nspline[3]=max(1,min(floor(maxfollow/spline_width[3]),nspline_max[3]))}
if (dim_bp >=4) {nspline[4]=max(1,min(floor(maxfollow/spline_width[4]),nspline_max[4]))}
if (dim_bp >=5) {nspline[5]=max(1,min(floor(maxfollow/spline_width[5]),nspline_max[5]))}
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

}
cat('--------------------------------------------','\n')
# 
# #--------------------------------------------------------------------------
# # END: Spline Creation 
# #--------------------------------------------------------------------------
#
# 
#-------------------------------------------------
#-------------------------------------------------
if (lmm1_ind==1 | lmm2_ind==1) {
  #-------------------------------------------------
  # LMM MODEL METHODS
  #-------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Define last time known event free
  #-----------------------------------------------------------------------------
  survive_time=dataset$stop
  nonsurvive_time=dataset$stop*(dataset$status==1)+design_time_max*(dataset$status!=1)
  censor_time=dataset$stop*(dataset$status==0)+design_time_max*(dataset$status!=0)

  #cat('Last Time Known to be Event Free','\n')
  #print(survive_time)
  #cat('Time Known to have Event','\n')
  #print(nonsurvive_time[1:20])

  #cat('------------------------------------------------','\n')
  #cat('------------------------------------------------','\n')

  if (lmm1_ind==1) {
    #----------------------------------------------------------------------
    # START LMM1: METHOD THAT IGNORES INFORMATIVE CENSORING FROM EVENTS
    #----------------------------------------------------------------------
    cat('---------------------------------------------','\n')
    cat('START LMM1','\n')
    cat('---------------------------------------------','\n')
    cat('---------------------------------------------','\n')
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
    
    names=paste(names,"+ (obstime | id)",sep="")
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
      
      names=paste(names,"+ (obstime | id)",sep="")
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
      
      names=paste(names,"+ (obstime | id)",sep="")
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
      
      names=paste(names,"+ (obstime | id)",sep="")
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
      
      names=paste(names,"+ (obstime | id)",sep="")
      #    cat('names =',names,'\n')
      (formula5=as.formula(paste("bp_curr~1+",names,sep="")))
    }
    
#    names=paste(names,"+(obstime | id)",sep="")
#    (formula1=as.formula(paste("bp_curr~baseage",paste(names,collapse="+"))))
#    (formula1=as.formula(paste("bp_curr~1+(obstime | id)")))
    
    
    dataset_nonmiss$id=as.numeric(dataset_nonmiss$id)
    dataset$id=as.numeric(dataset$id)
    n=length(dataset$id)
    
    
    
    dataset_impute=dataset
   
    dataset_impute$bp_curr1=rep(NA,times=n)
    dataset_impute$bp_curr2=rep(NA,times=n)
    dataset_impute$bp_curr3=rep(NA,times=n)
    dataset_impute$bp_curr4=rep(NA,times=n)
    dataset_impute$bp_curr5=rep(NA,times=n)
    # 
    # 
    #-------------------------------------------------
    #  FIT MODEL LMM1
    #-------------------------------------------------
    
    n_nonmiss=length(dataset_nonmiss$id)
    n_nonmiss_imp=length(dataset_nonmiss_imp$id)
    
    cat('--------------------------------------------','\n')
    cat('--------------------------------------------','\n')
    cat('LMM1 model','\n')
    cat('--------------------------------------------','\n')
    
    #select_id=unique(dataset_nonmiss_imp$id)
    n_eval_id=length(eval_id)
    
    cat('#Evaluable Subjects =',n_eval_id,'\n')
#    dataset_nonmiss_imp$id=as.numeric(dataset_nonmiss_imp$id)
    
    
    preds<-data.frame(matrix(NA, nrow=3, ncol=4))
    colnames(preds)=c('id', 'bptype', 'pred', 'real')
    #-------------------------------------------------
    #  LOOP OVER EVAL POP
    #-------------------------------------------------
    i=start_eval_id_num
    while (i<=end_eval_id_num){
        cat('i=',i,'\n')
        id=eval_id[i]
        cat('id=',id,'\n')
        temp_trial<-dataset_nonmiss_imp[!(dataset_nonmiss_imp$id==id & dataset_nonmiss_imp$obstime>s & !is.na(dataset_nonmiss_imp$bp_curr)) & dataset_nonmiss_imp$bptype=='TYPE1',]
        cat('temp_trial[1,]','\n')
        print(temp_trial[1,])
        temp_preds<-preds
        temp_preds$id[1]=id
        temp_preds$id[2]=id
        temp_preds$id[3]=id
        temp_preds$bptype[1]="TYPE1"
        temp_preds$bptype[2]="TYPE2"
        temp_preds$bptype[3]="TYPE3"
        cat('formula1','\n')
        print(formula1)
        cat('length(temp_trial)=',length(temp_trial$id),'\n')
        fit1=try(lmer(formula1,data=temp_trial, na.action=na.exclude))
        cat('TYPE 1 MODEL','\n')
        cat('------------------------------------','\n')
        print(summary(fit1))
        cat('------------------------------------','\n')
#        na_row<- c(id, NA, s+tau, "TYPE1", NA, NA, NA, NA, NA, NA, NA, NA, NA)
#        na_row=temp_trial[temp_trial$id==id,][1,]
#        na_row$id[1]=id
#        na_row$bp_curr[1]=NA
#        na_row$bptype[1]="TYPE1"
#        na_row$obstime[1]=s+tau
#        temp_trial<-rbind(temp_trial, na_row)
#        cat('temp_trial[length(temp_trial$id),]','\n')
#        print(temp_trial[length(temp_trial$id),])
        temp=temp_trial[is.na(temp_trial$bp_curr) & temp_trial$id==id & temp_trial$bptype=="TYPE1",]
        cat('temp=','\n')
        print(temp)
        fitted1_na=predict(fit1,newdata=temp,allow.new.levels=TRUE)
        cat('fitted11_na','\n')
        print(fitted1_na)
        temp_preds$pred[1]<-fitted1_na
        temp_preds$real[1]<-mean(dataset_nonmiss_imp$bp_curr[dataset_nonmiss_imp$id==id & abs(dataset_nonmiss_imp$obstime-s-tau) <= window_width & dataset_nonmiss_imp$bptype=='TYPE1' & !is.na(dataset_nonmiss_imp$bp_curr)])
        temp_preds$diffs[1]<-abs(temp_preds$real[1]-temp_preds$pred[1])
        if (dim_bp >=2) {
          temp_trial<-dataset_nonmiss_imp[!(dataset_nonmiss_imp$id==id & dataset_nonmiss_imp$obstime>s & !is.na(dataset_nonmiss_imp$bp_curr)) & dataset_nonmiss_imp$bptype=='TYPE2',]
          fit2=try(lmer(formula2,data=temp_trial, na.action=na.exclude))
          cat('TYPE 2 MODEL','\n')
          cat('------------------------------------','\n')
          print(summary(fit2))
          cat('------------------------------------','\n')
#          na_row$bptype[1]="TYPE2"
#          temp_trial<-rbind(temp_trial, na_row)
#          temp=temp_trial[is.na(temp_trial$bp_curr),]
          temp=temp_trial[is.na(temp_trial$bp_curr) & temp_trial$id==id & temp_trial$bptype=="TYPE2",]
          cat('temp=','\n')
          print(temp)
          fitted2_na=predict(fit2,newdata=temp,allow.new.levels=TRUE)
          temp_preds$pred[2]<-fitted2_na
          temp_preds$real[2]<-mean(dataset_nonmiss_imp$bp_curr[dataset_nonmiss_imp$id==id & abs(dataset_nonmiss_imp$obstime-s-tau) <= window_width & dataset_nonmiss_imp$bptype=='TYPE2' & !is.na(dataset_nonmiss_imp$bp_curr)])
          temp_preds$diffs[2]<-abs(temp_preds$real[2]-temp_preds$pred[2])
        }
        if (dim_bp >=3) {
          temp_trial<-dataset_nonmiss_imp[!(dataset_nonmiss_imp$id==id & dataset_nonmiss_imp$obstime>s & !is.na(dataset_nonmiss_imp$bp_curr)) & dataset_nonmiss_imp$bptype=='TYPE3',]
          fit3=try(lmer(formula3,data=temp_trial, na.action=na.exclude))
          cat('TYPE 3 MODEL','\n')
          cat('------------------------------------','\n')
          print(summary(fit3))
          cat('------------------------------------','\n')
#          na_row$bptype[1]="TYPE3"
#          temp_trial<-rbind(temp_trial, na_row)
#          temp=temp_trial[is.na(temp_trial$bp_curr),]
          temp=temp_trial[is.na(temp_trial$bp_curr) & temp_trial$id==id & temp_trial$bptype=="TYPE3",]
          fitted3_na=predict(fit3,newdata=temp,allow.new.levels=TRUE)
          temp_preds$pred[3]<-fitted3_na
          temp_preds$real[3]<-mean(dataset_nonmiss_imp$bp_curr[dataset_nonmiss_imp$id==id & abs(dataset_nonmiss_imp$obstime-s-tau) <= window_width & dataset_nonmiss_imp$bptype=='TYPE3' & !is.na(dataset_nonmiss_imp$bp_curr)])
          temp_preds$diffs[3]<-abs(temp_preds$real[3]-temp_preds$pred[3])
        }
        if (dim_bp >=4) {
          temp_trial<-dataset_nonmiss_imp[!(dataset_nonmiss_imp$id==id & dataset_nonmiss_imp$obstime>s & !is.na(dataset_nonmiss_imp$bp_curr)) & dataset_nonmiss_imp$bptype=='TYPE4',]
          fit4=try(lmer(formula4,data=temp_trial, na.action=na.exclude))
#          na_row$bptype[1]="TYPE4"
#          temp_trial<-rbind(temp_trial, na_row)
#          temp=temp_trial[is.na(temp_trial$bp_curr),]
          temp=temp_trial[is.na(temp_trial$bp_curr) & temp_trial$id==id & temp_trial$bptype=="TYPE4",]
          fitted4_na=predict(fit4,newdata=temp,allow.new.levels=TRUE)
          temp_preds$pred[4]<-fitted4_na
          temp_preds$real[4]<-mean(dataset_nonmiss_imp$bp_curr[dataset_nonmiss_imp$id==id & abs(dataset_nonmiss_imp$obstime-s-tau) <= window_width & dataset_nonmiss_imp$bptype=='TYPE4' & !is.na(dataset_nonmiss_imp$bp_curr)])
          temp_preds$diffs[4]<-abs(temp_preds$real[4]-temp_preds$pred[4])
         }
        if (dim_bp >=5) {
          temp_trial<-dataset_nonmiss_imp[!(dataset_nonmiss_imp$id==id & dataset_nonmiss_imp$obstime>s & !is.na(dataset_nonmiss_imp$bp_curr)) & dataset_nonmiss_imp$bptype=='TYPE5',]
          fit5=try(lmer(formula5,data=temp_trial, na.action=na.exclude))
#          na_row$bptype[1]="TYPE5"
#          temp_trial<-rbind(temp_trial, na_row)
#          temp=temp_trial[is.na(temp_trial$bp_curr),]
          temp=temp_trial[is.na(temp_trial$bp_curr) & temp_trial$id==id & temp_trial$bptype=="TYPE5",]
          fitted5_na=predict(fit5,newdata=temp,allow.new.levels=TRUE)
          temp_preds$pred[5]<-fitted5_na
          temp_preds$real[5]<-mean(dataset_nonmiss_imp$bp_curr[dataset_nonmiss_imp$id==id & abs(dataset_nonmiss_imp$obstime-s-tau) <= window_width & dataset_nonmiss_imp$bptype=='TYPE5' & !is.na(dataset_nonmiss_imp$bp_curr)])
          temp_preds$diffs[5]<-abs(temp_preds$real[5]-temp_preds$pred[5])
        }
        
        if (i != start_eval_id_num) {all_preds<-rbind(all_preds, temp_preds)}
        if (i==start_eval_id_num) {all_preds=temp_preds}
      i=i+1
    }
    cat('Predictions=','\n')
    print(all_preds)
    cat('----------------------------------------------', '\n')
    cat('----------------------------------------------', '\n')
    cat('MEASUREMENTS FOR BPTYPE 1', '\n')
    mean_pred_1=mean(all_preds[all_preds$bptype == 'TYPE1', 3])
    sd_pred_1=sd(all_preds[all_preds$bptype == 'TYPE1', 3])
    mean_true_1= mean(all_preds[all_preds$bptype == 'TYPE1', 4])
    sd_true_1= sd(all_preds[all_preds$bptype == 'TYPE1', 4])
    mean_diff_1= mean(all_preds[all_preds$bptype == 'TYPE1', 5])
    sd_diff_1=sd(all_preds[all_preds$bptype == 'TYPE1', 5])
    cat('predicted mean=', mean_pred_1, '\n')
    cat('predicted sd=', sd_pred_1, '\n')
    cat('true mean=', mean_true_1, '\n')
    cat('true sd=', sd_true_1, '\n')
    cat('mean abs difference=', mean_diff_1, '\n')
    cat('sd of differences=', sd_diff_1, '\n')
    if (dim_bp >=2) {
      cat('----------------------------------------------', '\n')
      cat('----------------------------------------------', '\n')
      cat('MEASUREMENTS FOR BPTYPE2', '\n')
      mean_pred_2=mean(all_preds[all_preds$bptype == 'TYPE2', 3])
      sd_pred_2=sd(all_preds[all_preds$bptype == 'TYPE2', 3])
      mean_true_2= mean(all_preds[all_preds$bptype == 'TYPE2', 4])
      sd_true_2= sd(all_preds[all_preds$bptype == 'TYPE2', 4])
      mean_diff_2= mean(all_preds[all_preds$bptype == 'TYPE2', 5])
      sd_diff_2=sd(all_preds[all_preds$bptype == 'TYPE2', 5])
      cat('predicted mean=', mean_pred_2, '\n')
      cat('predicted sd=', sd_pred_2, '\n')
      cat('true mean=', mean_true_2, '\n')
      cat('true sd=', sd_true_2, '\n')
      cat('mean abs difference=', mean_diff_2, '\n')
      cat('sd of differences=', sd_diff_2, '\n')
    }
    if (dim_bp >=3) {
      cat('----------------------------------------------', '\n')
      cat('----------------------------------------------', '\n')
      cat('MEASUREMENTS FOR BPTYPE3', '\n')
      mean_pred_3=mean(all_preds[all_preds$bptype == 'TYPE3', 3])
      sd_pred_3=sd(all_preds[all_preds$bptype == 'TYPE3', 3])
      mean_true_3= mean(all_preds[all_preds$bptype == 'TYPE3', 4])
      sd_true_3= sd(all_preds[all_preds$bptype == 'TYPE3', 4])
      mean_diff_3= mean(all_preds[all_preds$bptype == 'TYPE3', 5])
      sd_diff_3=sd(all_preds[all_preds$bptype == 'TYPE3', 5])
      cat('predicted mean=', mean_pred_3, '\n')
      cat('predicted sd=', sd_pred_3, '\n')
      cat('true mean=', mean_true_3, '\n')
      cat('true sd=', sd_true_3, '\n')
      cat('mean abs difference=', mean_diff_3, '\n')
      cat('sd of differences=', sd_diff_3, '\n')
    }
    
    # fit11=try(lmer(formula1,data=subset(dataset_nonmiss_imp, dataset_nonmiss_imp$bptype=='TYPE1'), na.action=na.exclude))
    # fit12=try(lmer(formula1,data=subset(dataset_nonmiss_imp, dataset_nonmiss_imp$bptype=='TYPE2'), na.action=na.exclude))
    # fit13=try(lmer(formula1,data=subset(dataset_nonmiss_imp, dataset_nonmiss_imp$bptype=='TYPE3'), na.action=na.exclude))
    # cat('Model LMM1 ignoring informative censoring from events=','\n')
    # print(summary(fit11))
    # print(summary(fit12))
    # print(summary(fit13))
    
    # fitted11_na=predict(fit11,newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr),])
    # fitted12_na=predict(fit12,newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr),])
    # fitted13_na=predict(fit13,newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr),])
    
    # All Patients
    # Replace ALL values of Biomarker Process with fitted value from model
    # dataset_impute$bp_curr1[survive_time[dataset_impute$id] >= s+tau | dataset_impute$id==1]=fitted11_na
    # dataset_impute$bp_curr2[survive_time[dataset_impute$id] >= s+tau | dataset_impute$id==1]=fitted12_na
    # dataset_impute$bp_curr3[survive_time[dataset_impute$id] >= s+tau | dataset_impute$id==1]=fitted13_na
    # #cat('Dataset after imputed BP_curr =','\n')
    # #print(dataset_impute[1:40,])
    
    # cat('END LMM1 model ignoring informative censoring from events','\n')
    # cat('--------------------------------------------','\n')
    # 
    # 
    # lmm1_pred=cbind(TYPE1=dataset_impute$bp_curr1[1], TYPE2=dataset_impute$bp_curr2[1], TYPE3=dataset_impute$bp_curr3[1])
    # lmm1_pred=as.data.frame(lmm1_pred)
    # cat('--------------------------------------------','\n')
    # cat('--------------------------------------------','\n')
    # cat('LMM1 predicted value =',lmm1_pred,'\n')
    # cat('--------------------------------------------','\n')
    # cat('True BP value for current subject=',bp_st,'\n')
    # cat('---------------------------------------------------','\n')
    
  }
  #-------------------------------------------------------------------
  # END LMM1: IGNORES INFORMATIVE CENSORING FROM EVENTS
  #-------------------------------------------------------------------
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # if (lmm2_ind==1) {
  #   #--------------------------------------------------------------------------
  #   # START LMM2: METHOD THAT ACCOUNTS FOR INFORMATIVE CENSORING FROM EVENTS
  #   #--------------------------------------------------------------------------
  #   cat('---------------------------------------------','\n')
  #   cat('START LMM2','\n')
  #   cat('---------------------------------------------','\n')
  #   cat('---------------------------------------------','\n')
  #   #-----------------------------------------------------
  #   # CREATE ORDERED FORMULAS FOR LONGITUDINAL MODELS
  #   #-----------------------------------------------------
  #   names=""
  #   if (knot2==1) {
  #     names=paste(names,"+spline11+spline12+spline13",sep="")
  #   }
  #   if (knot2==2) {
  #     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23",sep="")
  #   }
  #   if (knot2==3) {
  #     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33",sep="")
  #   }
  #   if (knot2==4) {
  #     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43",sep="")
  #   }
  #   if (knot2==5) {
  #     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43+spline51+spline52+spline53",sep="")
  #   }
  #   if (knot2==6) {
  #     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43+spline51+spline52+spline53+spline61+spline62+spline63",sep="")
  #   }
  #   if (knot2==7) {
  #     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43+spline51+spline52+spline53+spline61+spline62+spline63+spline71+spline72+spline73",sep="")
  #   }
  #   if (knot2==8) {
  #     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43+spline51+spline52+spline53+spline61+spline62+spline63+spline71+spline72+spline73+spline81+spline82+spline83",sep="")
  #   }
  #   
  #   names=paste(names,"+nonsurvive+nonsurvive1",sep="")
  #   names=paste(names,"+(obstime | id)",sep="")
  #   (formula2=as.formula(paste("bp_curr~baseage",paste(names,collapse="+"))))
  #   
  #   censor_time[1]=s
  #   survive_time[1]=s+tau
  #   dataset_nonmiss$id=as.numeric(dataset_nonmiss$id)
  #   dataset$id=as.numeric(dataset$id)
  #   n=length(dataset$id)
  #   
  #   #--------------------------------------------------------------------------
  #   #  CALCULATE NUMBER OF SUBJECTS AND OBS OF EACH TYPE
  #   #--------------------------------------------------------------------------
  #   num_subj_nonsurv=sum(nonsurvive_time <= s+tau)
  #   num_subj_surv=sum(survive_time >= s+tau)
  #   num_subj_cens=sum(nonsurvive_time > s+tau & survive_time < s+tau)
  #   num_obs_nonsurv=sum(nonsurvive_time[dataset_nonmiss$id] <= s+tau)
  #   num_obs_cens=sum(nonsurvive_time[dataset_nonmiss$id] > s+tau & survive_time[dataset_nonmiss$id] < s+tau)
  #   num_obs_surv=sum(survive_time[dataset_nonmiss$id] >= s+tau & dataset_nonmiss$obstime <= s+tau)
  #   
  #   
  #   cat('number subjects events prior to s+tau=',num_subj_nonsurv,'\n')
  #   cat('number subjects censored prior to s+tau=',num_subj_cens,'\n')
  #   cat('number subjects surviving to s+tau=',num_subj_surv,'\n')
  #   cat('------------------------------------------------------------------','\n')
  #   cat('number observations of nonsurvivors=',num_obs_nonsurv,'\n')
  #   cat('number observations of censored=',num_obs_cens,'\n')
  #   cat('number past observations of survivors=',num_obs_surv,'\n')
  #   
  #   
  #   cat('--------------------------------------------------','\n')
  #   cat('Number of subjects (including current) surviving to s+tau=',length(nonsurvive_time[nonsurvive_time > s+tau]),'\n')
  #   
  #   if (lmm1_ind==0) {
  #     #---------------------------------------------------------------------------------------------------
  #   # FOR EACH ID SURVIVING PAST TIME=S+TAU ADD OBS WITH TIME=S+TAU AND BP_CURR=NA TO DATASET_NONMISS 
  #   #---------------------------------------------------------------------------------------------------
  #   i=1
  #   temp_start=0
  #   idprev=0
  #   while (i <= n_nonmiss) {
  #     #  cat('i=',i,'\n')
  #     id=dataset_nonmiss$id[i]
  #     #  cat('id=',id,'\n')
  #     #  cat('idprev=',idprev,'\n')
  #     if (id != idprev & (survive_time[id] >= s+tau | id==1) ) {
  #       if (temp_start==0) {
  #         temp=dataset_nonmiss[i,]
  #         temp$bp_curr[1]=NA
  #         temp$obstime[1]=s+tau
  #         temp$bptype='ANY'
  #         temp_start=1
  #       } else {
  #         temp1=dataset_nonmiss[i,]
  #         temp1$bp_curr[1]=NA
  #         temp1$obstime[1]=s+tau
  #         temp1$bptype='ANY'
  #         temp=rbind(temp,temp1) 
  #       }
  #     }
  #     idprev=id
  #     i=i+1
  #   }
  #   dataset_nonmiss_imp=rbind(dataset_nonmiss,temp)
  #   dataset_nonmiss_imp=dataset_nonmiss_imp[order(dataset_nonmiss_imp$id,dataset_nonmiss_imp$obstime),]
  #   n_nonmiss_imp=length(dataset_nonmiss_imp$id)
  #   
  # }
  # #cat('Dataset_nonmiss with NAs =','\n')
  # #print(dataset_nonmiss_imp[1:20,])
  # #cat('dataset =','\n')
  # #print(dataset[1:50,1:5])
  # #cat('survive_time=','\n')
  # #print(survive_time[1:50])
  # #cat('censor_time=','\n')
  # #print(censor_time[1:50])
  # 
  # #-------------------------------------------------------------
  # #  SPLINE CREATION
  # #-------------------------------------------------------------
  # dataset_nonmiss_imp$spline11=dataset_nonmiss_imp$obstime
  # dataset_nonmiss_imp$spline12=dataset_nonmiss_imp$obstime**2
  # dataset_nonmiss_imp$spline13=dataset_nonmiss_imp$obstime**3
  # if (knot2 >= 2) {
  #   dataset_nonmiss_imp$spline11=dataset_nonmiss_imp$obstime*(dataset_nonmiss_imp$obstime < knot1) +knot1*(dataset_nonmiss_imp$obstime >= knot1)
  #   dataset_nonmiss_imp$spline12=dataset_nonmiss_imp$obstime**2*(dataset_nonmiss_imp$obstime < knot1) +knot1**2*(dataset_nonmiss_imp$obstime >= knot1)
  #   dataset_nonmiss_imp$spline13=dataset_nonmiss_imp$obstime**3*(dataset_nonmiss_imp$obstime < knot1) +knot1**3*(dataset_nonmiss_imp$obstime >= knot1)
  #   dataset_nonmiss_imp$spline21=(dataset_nonmiss_imp$obstime-knot1)*(dataset_nonmiss_imp$obstime >= knot1) 
  #   dataset_nonmiss_imp$spline22=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime >= knot1) 
  #   dataset_nonmiss_imp$spline23=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime >= knot1) 
  # }
  # if (knot2 >= 3) {
  #   dataset_nonmiss_imp$spline21=(dataset_nonmiss_imp$obstime-knot1)*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)*(dataset_nonmiss_imp$obstime >= knot2) 
  #   dataset_nonmiss_imp$spline22=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**2*(dataset_nonmiss_imp$obstime >= knot2) 
  #   dataset_nonmiss_imp$spline23=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**3*(dataset_nonmiss_imp$obstime >= knot2) 
  #   
  #   dataset_nonmiss_imp$spline31=(dataset_nonmiss_imp$obstime-knot2)*(dataset_nonmiss_imp$obstime >= knot2) 
  #   dataset_nonmiss_imp$spline32=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime >= knot2) 
  #   dataset_nonmiss_imp$spline33=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime >= knot2) 
  # }
  # if (knot2 >= 4) {
  #   dataset_nonmiss_imp$spline31=(dataset_nonmiss_imp$obstime-knot2)*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)*(dataset_nonmiss_imp$obstime >= knot3) 
  #   dataset_nonmiss_imp$spline32=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**2*(dataset_nonmiss_imp$obstime >= knot3) 
  #   dataset_nonmiss_imp$spline33=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**3*(dataset_nonmiss_imp$obstime >= knot3) 
  #   
  #   dataset_nonmiss_imp$spline41=(dataset_nonmiss_imp$obstime-knot3)*(dataset_nonmiss_imp$obstime >= knot3) 
  #   dataset_nonmiss_imp$spline42=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime >= knot3) 
  #   dataset_nonmiss_imp$spline43=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime >= knot3) 
  # }
  # if (knot2 >= 5) {
  #   dataset_nonmiss_imp$spline41=(dataset_nonmiss_imp$obstime-knot3)*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)*(dataset_nonmiss_imp$obstime >= knot4) 
  #   dataset_nonmiss_imp$spline42=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**2*(dataset_nonmiss_imp$obstime >= knot4) 
  #   dataset_nonmiss_imp$spline43=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**3*(dataset_nonmiss_imp$obstime >= knot4) 
  #   
  #   dataset_nonmiss_imp$spline51=(dataset_nonmiss_imp$obstime-knot4)*(dataset_nonmiss_imp$obstime >= knot4) 
  #   dataset_nonmiss_imp$spline52=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime >= knot4) 
  #   dataset_nonmiss_imp$spline53=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime >= knot4) 
  # }
  # if (knot2 >= 6) {
  #   dataset_nonmiss_imp$spline51=(dataset_nonmiss_imp$obstime-knot4)*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)*(dataset_nonmiss_imp$obstime >= knot5) 
  #   dataset_nonmiss_imp$spline52=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**2*(dataset_nonmiss_imp$obstime >= knot5) 
  #   dataset_nonmiss_imp$spline53=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**3*(dataset_nonmiss_imp$obstime >= knot5) 
  #   
  #   dataset_nonmiss_imp$spline61=(dataset_nonmiss_imp$obstime-knot5)*(dataset_nonmiss_imp$obstime >= knot5) 
  #   dataset_nonmiss_imp$spline62=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime >= knot5) 
  #   dataset_nonmiss_imp$spline63=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime >= knot5) 
  # }
  # if (knot2 >= 7) {
  #   dataset_nonmiss_imp$spline61=(dataset_nonmiss_imp$obstime-knot5)*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)*(dataset_nonmiss_imp$obstime >= knot6) 
  #   dataset_nonmiss_imp$spline62=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**2*(dataset_nonmiss_imp$obstime >= knot6) 
  #   dataset_nonmiss_imp$spline63=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**3*(dataset_nonmiss_imp$obstime >= knot6) 
  #   
  #   dataset_nonmiss_imp$spline71=(dataset_nonmiss_imp$obstime-knot6)*(dataset_nonmiss_imp$obstime >= knot6) 
  #   dataset_nonmiss_imp$spline72=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime >= knot6) 
  #   dataset_nonmiss_imp$spline73=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime >= knot6) 
  # }
  # if (knot2 >= 8) {
  #   dataset_nonmiss_imp$spline71=(dataset_nonmiss_imp$obstime-knot6)*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)*(dataset_nonmiss_imp$obstime >= knot7) 
  #   dataset_nonmiss_imp$spline72=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**2*(dataset_nonmiss_imp$obstime >= knot7) 
  #   dataset_nonmiss_imp$spline73=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**3*(dataset_nonmiss_imp$obstime >= knot7) 
  #   
  #   dataset_nonmiss_imp$spline81=(dataset_nonmiss_imp$obstime-knot7)*(dataset_nonmiss_imp$obstime >= knot7) 
  #   dataset_nonmiss_imp$spline82=(dataset_nonmiss_imp$obstime-knot7)**2*(dataset_nonmiss_imp$obstime >= knot7) 
  #   dataset_nonmiss_imp$spline83=(dataset_nonmiss_imp$obstime-knot7)**3*(dataset_nonmiss_imp$obstime >= knot7) 
  # }
  # 
  # #-------------------------------------------------------------
  # #-------------------------------------------------------------
  # 
  # dataset_impute2=dataset
  # n=length(dataset_impute2$id)
  # dataset_impute2$bp_curr1=rep(NA,times=n)
  # dataset_impute2$bp_curr2=rep(NA,times=n)
  # dataset_impute2$bp_curr3=rep(NA,times=n)
  # dataset_impute2$bp_curr4=rep(NA,times=n)
  # dataset_impute2$bp_curr5=rep(NA,times=n)
  # 
  # #-------------------------------------------------
  # #  FIT MODEL LMM2
  # #-------------------------------------------------
  # 
  # n_nonmiss=length(dataset_nonmiss$id)
  # n_nonmiss_imp=length(dataset_nonmiss_imp$id)
  # 
  # cat('--------------------------------------------','\n')
  # cat('--------------------------------------------','\n')
  # cat('LMM2 model','\n')
  # cat('--------------------------------------------','\n')
  
  # 
  # dataset_nonmiss_imp$nonsurvive=(nonsurvive_time[dataset_nonmiss_imp$id] <= s+tau)
  # dataset_nonmiss_imp$nonsurvive1=dataset_nonmiss_imp$nonsurvive*dataset_nonmiss_imp$obstime
  # 
  # select_id=unique(dataset_nonmiss_imp$id)
  # n_select_id=length(select_id)
  # 
  # cat('#Subjects included in model =',n_select_id,'\n')
  # 
  # dataset_nonmiss_imp=dataset_nonmiss_imp[dataset_nonmiss_imp$obstime <= s+tau,]
  # 
  # #-------------------------------------------------
  # # FIT LMM2
  # #-------------------------------------------------
  # fit21=try(lmer(formula2,data=subset(dataset_nonmiss_imp, dataset_nonmiss_imp$bptype=='TYPE1'), na.action=na.exclude))
  # fit22=try(lmer(formula2,data=subset(dataset_nonmiss_imp, dataset_nonmiss_imp$bptype=='TYPE2'), na.action=na.exclude))
  # fit23=try(lmer(formula2,data=subset(dataset_nonmiss_imp, dataset_nonmiss_imp$bptype=='TYPE3'), na.action=na.exclude))
  # cat('Model LMM2=','\n')
  # print(summary(fit21))
  # print(summary(fit22))
  # print(summary(fit23))
  # 
  # fitted21_na=predict(fit21,newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr),])
  # fitted22_na=predict(fit22,newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr),])
  # fitted23_na=predict(fit23,newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr),])
  # 
  # #cat('--------------------------------------------------','\n')
  # #cat('fitted1_na=','\n')
  # #print(fitted1_na[1:20])
  # #cat('Dataset before imputed BP_curr =','\n')
  # #print(dataset_impute[1:40,])
  # 
  # # All Patients
  # # Replace ALL values of Biomarker Process with fitted value from model
  # dataset_impute2$bp_curr1[survive_time[dataset_impute2$id] >= s+tau | dataset_impute2$id==1]=fitted21_na
  # dataset_impute2$bp_curr2[survive_time[dataset_impute2$id] >= s+tau | dataset_impute2$id==1]=fitted22_na
  # dataset_impute2$bp_curr3[survive_time[dataset_impute2$id] >= s+tau | dataset_impute2$id==1]=fitted23_na
  # 
  # #cat('Dataset after imputed BP_curr =','\n')
  # #print(dataset_impute[1:40,])
  # 
  # cat('END LMM2 model','\n')
  # cat('--------------------------------------------','\n')
  # 
  # 
  # lmm2_pred=cbind(TYPE1=dataset_impute2$bp_curr1[1], TYPE2=dataset_impute2$bp_curr2[1], TYPE3=dataset_impute2$bp_curr3[1])
  # }
  
  # cat('True BP value for current subject=',bp_st,'\n')
  # cat('---------------------------------------------------','\n')
  
#}
#-------------------------------------------------------------------
# END LMM1: IGNORES INFORMATIVE CENSORING FROM EVENTS
#-------------------------------------------------------------------
















# if (lmm2_ind==1) {
#   #--------------------------------------------------------------------------
#   # START LMM2: METHOD THAT ACCOUNTS FOR INFORMATIVE CENSORING FROM EVENTS
#   #--------------------------------------------------------------------------
#   cat('---------------------------------------------','\n')
#   cat('START LMM2','\n')
#   cat('---------------------------------------------','\n')
#   cat('---------------------------------------------','\n')
#   #-----------------------------------------------------
#   # CREATE ORDERED FORMULAS FOR LONGITUDINAL MODELS
#   #-----------------------------------------------------
#   names=""
#   if (knot2==1) {
#     names=paste(names,"+spline11+spline12+spline13",sep="")
#   }
#   if (knot2==2) {
#     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23",sep="")
#   }
#   if (knot2==3) {
#     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33",sep="")
#   }
#   if (knot2==4) {
#     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43",sep="")
#   }
#   if (knot2==5) {
#     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43+spline51+spline52+spline53",sep="")
#   }
#   if (knot2==6) {
#     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43+spline51+spline52+spline53+spline61+spline62+spline63",sep="")
#   }
#   if (knot2==7) {
#     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43+spline51+spline52+spline53+spline61+spline62+spline63+spline71+spline72+spline73",sep="")
#   }
#   if (knot2==8) {
#     names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43+spline51+spline52+spline53+spline61+spline62+spline63+spline71+spline72+spline73+spline81+spline82+spline83",sep="")
#   }
#   
#   names=paste(names,"+nonsurvive+nonsurvive1",sep="")
#   names=paste(names,"+(obstime | id)",sep="")
#   (formula2=as.formula(paste("bp_curr~baseage",paste(names,collapse="+"))))
#   
#   censor_time[1]=s
#   survive_time[1]=s+tau
#   dataset_nonmiss$id=as.numeric(dataset_nonmiss$id)
#   dataset$id=as.numeric(dataset$id)
#   n=length(dataset$id)
#   
#   #--------------------------------------------------------------------------
#   #  CALCULATE NUMBER OF SUBJECTS AND OBS OF EACH TYPE
#   #--------------------------------------------------------------------------
#   num_subj_nonsurv=sum(nonsurvive_time <= s+tau)
#   num_subj_surv=sum(survive_time >= s+tau)
#   num_subj_cens=sum(nonsurvive_time > s+tau & survive_time < s+tau)
#   num_obs_nonsurv=sum(nonsurvive_time[dataset_nonmiss$id] <= s+tau)
#   num_obs_cens=sum(nonsurvive_time[dataset_nonmiss$id] > s+tau & survive_time[dataset_nonmiss$id] < s+tau)
#   num_obs_surv=sum(survive_time[dataset_nonmiss$id] >= s+tau & dataset_nonmiss$obstime <= s+tau)
#   
#   
#   cat('number subjects events prior to s+tau=',num_subj_nonsurv,'\n')
#   cat('number subjects censored prior to s+tau=',num_subj_cens,'\n')
#   cat('number subjects surviving to s+tau=',num_subj_surv,'\n')
#   cat('------------------------------------------------------------------','\n')
#   cat('number observations of nonsurvivors=',num_obs_nonsurv,'\n')
#   cat('number observations of censored=',num_obs_cens,'\n')
#   cat('number past observations of survivors=',num_obs_surv,'\n')
#   
#   
#   cat('--------------------------------------------------','\n')
#   cat('Number of subjects (including current) surviving to s+tau=',length(nonsurvive_time[nonsurvive_time > s+tau]),'\n')
#   
#   if (lmm1_ind==0) {
#     #---------------------------------------------------------------------------------------------------
#   # FOR EACH ID SURVIVING PAST TIME=S+TAU ADD OBS WITH TIME=  AND BP_CURR=NA TO DATASET_NONMISS 
#   #---------------------------------------------------------------------------------------------------
#   i=1
#   temp_start=0
#   idprev=0
#   while (i <= n_nonmiss) {
#     #  cat('i=',i,'\n')
#     id=dataset_nonmiss$id[i]
#     #  cat('id=',id,'\n')
#     #  cat('idprev=',idprev,'\n')
#     if (id != idprev & (survive_time[id] >= s+tau | id==1) ) {
#       if (temp_start==0) {
#         temp=dataset_nonmiss[i,]
#         temp$bp_curr[1]=NA
#         temp$obstime[1]=s+tau
#         temp$bptype='ANY'
#         temp_start=1
#       } else {
#         temp1=dataset_nonmiss[i,]
#         temp1$bp_curr[1]=NA
#         temp1$obstime[1]=s+tau
#         temp1$bptype='ANY'
#         temp=rbind(temp,temp1) 
#       }
#     }
#     idprev=id
#     i=i+1
#   }
#   dataset_nonmiss_imp=rbind(dataset_nonmiss,temp)
#   dataset_nonmiss_imp=dataset_nonmiss_imp[order(dataset_nonmiss_imp$id,dataset_nonmiss_imp$obstime),]
#   n_nonmiss_imp=length(dataset_nonmiss_imp$id)
#   
# }
# #cat('Dataset_nonmiss with NAs =','\n')
# #print(dataset_nonmiss_imp[1:20,])
# #cat('dataset =','\n')
# #print(dataset[1:50,1:5])
# #cat('survive_time=','\n')
# #print(survive_time[1:50])
# #cat('censor_time=','\n')
# #print(censor_time[1:50])
# 
# #-------------------------------------------------------------
# #  SPLINE CREATION
# #-------------------------------------------------------------
# dataset_nonmiss_imp$spline11=dataset_nonmiss_imp$obstime
# dataset_nonmiss_imp$spline12=dataset_nonmiss_imp$obstime**2
# dataset_nonmiss_imp$spline13=dataset_nonmiss_imp$obstime**3
# if (knot2 >= 2) {
#   dataset_nonmiss_imp$spline11=dataset_nonmiss_imp$obstime*(dataset_nonmiss_imp$obstime < knot1) +knot1*(dataset_nonmiss_imp$obstime >= knot1)
#   dataset_nonmiss_imp$spline12=dataset_nonmiss_imp$obstime**2*(dataset_nonmiss_imp$obstime < knot1) +knot1**2*(dataset_nonmiss_imp$obstime >= knot1)
#   dataset_nonmiss_imp$spline13=dataset_nonmiss_imp$obstime**3*(dataset_nonmiss_imp$obstime < knot1) +knot1**3*(dataset_nonmiss_imp$obstime >= knot1)
#   dataset_nonmiss_imp$spline21=(dataset_nonmiss_imp$obstime-knot1)*(dataset_nonmiss_imp$obstime >= knot1) 
#   dataset_nonmiss_imp$spline22=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime >= knot1) 
#   dataset_nonmiss_imp$spline23=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime >= knot1) 
# }
# if (knot2 >= 3) {
#   dataset_nonmiss_imp$spline21=(dataset_nonmiss_imp$obstime-knot1)*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)*(dataset_nonmiss_imp$obstime >= knot2) 
#   dataset_nonmiss_imp$spline22=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**2*(dataset_nonmiss_imp$obstime >= knot2) 
#   dataset_nonmiss_imp$spline23=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**3*(dataset_nonmiss_imp$obstime >= knot2) 
#   
#   dataset_nonmiss_imp$spline31=(dataset_nonmiss_imp$obstime-knot2)*(dataset_nonmiss_imp$obstime >= knot2) 
#   dataset_nonmiss_imp$spline32=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime >= knot2) 
#   dataset_nonmiss_imp$spline33=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime >= knot2) 
# }
# if (knot2 >= 4) {
#   dataset_nonmiss_imp$spline31=(dataset_nonmiss_imp$obstime-knot2)*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)*(dataset_nonmiss_imp$obstime >= knot3) 
#   dataset_nonmiss_imp$spline32=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**2*(dataset_nonmiss_imp$obstime >= knot3) 
#   dataset_nonmiss_imp$spline33=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**3*(dataset_nonmiss_imp$obstime >= knot3) 
#   
#   dataset_nonmiss_imp$spline41=(dataset_nonmiss_imp$obstime-knot3)*(dataset_nonmiss_imp$obstime >= knot3) 
#   dataset_nonmiss_imp$spline42=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime >= knot3) 
#   dataset_nonmiss_imp$spline43=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime >= knot3) 
# }
# if (knot2 >= 5) {
#   dataset_nonmiss_imp$spline41=(dataset_nonmiss_imp$obstime-knot3)*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)*(dataset_nonmiss_imp$obstime >= knot4) 
#   dataset_nonmiss_imp$spline42=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**2*(dataset_nonmiss_imp$obstime >= knot4) 
#   dataset_nonmiss_imp$spline43=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**3*(dataset_nonmiss_imp$obstime >= knot4) 
#   
#   dataset_nonmiss_imp$spline51=(dataset_nonmiss_imp$obstime-knot4)*(dataset_nonmiss_imp$obstime >= knot4) 
#   dataset_nonmiss_imp$spline52=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime >= knot4) 
#   dataset_nonmiss_imp$spline53=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime >= knot4) 
# }
# if (knot2 >= 6) {
#   dataset_nonmiss_imp$spline51=(dataset_nonmiss_imp$obstime-knot4)*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)*(dataset_nonmiss_imp$obstime >= knot5) 
#   dataset_nonmiss_imp$spline52=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**2*(dataset_nonmiss_imp$obstime >= knot5) 
#   dataset_nonmiss_imp$spline53=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**3*(dataset_nonmiss_imp$obstime >= knot5) 
#   
#   dataset_nonmiss_imp$spline61=(dataset_nonmiss_imp$obstime-knot5)*(dataset_nonmiss_imp$obstime >= knot5) 
#   dataset_nonmiss_imp$spline62=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime >= knot5) 
#   dataset_nonmiss_imp$spline63=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime >= knot5) 
# }
# if (knot2 >= 7) {
#   dataset_nonmiss_imp$spline61=(dataset_nonmiss_imp$obstime-knot5)*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)*(dataset_nonmiss_imp$obstime >= knot6) 
#   dataset_nonmiss_imp$spline62=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**2*(dataset_nonmiss_imp$obstime >= knot6) 
#   dataset_nonmiss_imp$spline63=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**3*(dataset_nonmiss_imp$obstime >= knot6) 
#   
#   dataset_nonmiss_imp$spline71=(dataset_nonmiss_imp$obstime-knot6)*(dataset_nonmiss_imp$obstime >= knot6) 
#   dataset_nonmiss_imp$spline72=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime >= knot6) 
#   dataset_nonmiss_imp$spline73=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime >= knot6) 
# }
# if (knot2 >= 8) {
#   dataset_nonmiss_imp$spline71=(dataset_nonmiss_imp$obstime-knot6)*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)*(dataset_nonmiss_imp$obstime >= knot7) 
#   dataset_nonmiss_imp$spline72=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**2*(dataset_nonmiss_imp$obstime >= knot7) 
#   dataset_nonmiss_imp$spline73=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**3*(dataset_nonmiss_imp$obstime >= knot7) 
#   
#   dataset_nonmiss_imp$spline81=(dataset_nonmiss_imp$obstime-knot7)*(dataset_nonmiss_imp$obstime >= knot7) 
#   dataset_nonmiss_imp$spline82=(dataset_nonmiss_imp$obstime-knot7)**2*(dataset_nonmiss_imp$obstime >= knot7) 
#   dataset_nonmiss_imp$spline83=(dataset_nonmiss_imp$obstime-knot7)**3*(dataset_nonmiss_imp$obstime >= knot7) 
# }
# 
# #-------------------------------------------------------------
# #-------------------------------------------------------------
# 
# dataset_impute2=dataset
# n=length(dataset_impute2$id)
# dataset_impute2$bp_curr1=rep(NA,times=n)
# dataset_impute2$bp_curr2=rep(NA,times=n)
# dataset_impute2$bp_curr3=rep(NA,times=n)
# dataset_impute2$bp_curr4=rep(NA,times=n)
# dataset_impute2$bp_curr5=rep(NA,times=n)
# 
# #-------------------------------------------------
# #  FIT MODEL LMM2
# #-------------------------------------------------
# 
# n_nonmiss=length(dataset_nonmiss$id)
# n_nonmiss_imp=length(dataset_nonmiss_imp$id)
# 
# cat('--------------------------------------------','\n')
# cat('--------------------------------------------','\n')
# cat('LMM2 model','\n')
# cat('--------------------------------------------','\n')
# 
# dataset_nonmiss_imp$nonsurvive=(nonsurvive_time[dataset_nonmiss_imp$id] <= s+tau)
# dataset_nonmiss_imp$nonsurvive1=dataset_nonmiss_imp$nonsurvive*dataset_nonmiss_imp$obstime
# 
# select_id=unique(dataset_nonmiss_imp$id)
# n_select_id=length(select_id)
# 
# cat('#Subjects included in model =',n_select_id,'\n')
# 
# dataset_nonmiss_imp=dataset_nonmiss_imp[dataset_nonmiss_imp$obstime <= s+tau,]
# 
# #-------------------------------------------------
# # FIT LMM2
# #-------------------------------------------------
# fit21=try(lmer(formula2,data=subset(dataset_nonmiss_imp, dataset_nonmiss_imp$bptype=='TYPE1'), na.action=na.exclude))
# fit22=try(lmer(formula2,data=subset(dataset_nonmiss_imp, dataset_nonmiss_imp$bptype=='TYPE2'), na.action=na.exclude))
# fit23=try(lmer(formula2,data=subset(dataset_nonmiss_imp, dataset_nonmiss_imp$bptype=='TYPE3'), na.action=na.exclude))
# cat('Model LMM2=','\n')
# print(summary(fit21))
# print(summary(fit22))
# print(summary(fit23))
# 
# fitted21_na=predict(fit21,newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr),])
# fitted22_na=predict(fit22,newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr),])
# fitted23_na=predict(fit23,newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr),])
# 
# #cat('--------------------------------------------------','\n')
# #cat('fitted1_na=','\n')
# #print(fitted1_na[1:20])
# #cat('Dataset before imputed BP_curr =','\n')
# #print(dataset_impute[1:40,])
# 
# # All Patients
# # Replace ALL values of Biomarker Process with fitted value from model
# dataset_impute2$bp_curr1[survive_time[dataset_impute2$id] >= s+tau | dataset_impute2$id==1]=fitted21_na
# dataset_impute2$bp_curr2[survive_time[dataset_impute2$id] >= s+tau | dataset_impute2$id==1]=fitted22_na
# dataset_impute2$bp_curr3[survive_time[dataset_impute2$id] >= s+tau | dataset_impute2$id==1]=fitted23_na
# 
# #cat('Dataset after imputed BP_curr =','\n')
# #print(dataset_impute[1:40,])
# 
# cat('END LMM2 model','\n')
# cat('--------------------------------------------','\n')
# 
# 
# lmm2_pred=cbind(TYPE1=dataset_impute2$bp_curr1[1], TYPE2=dataset_impute2$bp_curr2[1], TYPE3=dataset_impute2$bp_curr3[1])
# }

if (lmm2_ind==1) {
  #--------------------------------------------------------------------------
  # START LMM2: METHOD THAT ACCOUNTS FOR INFORMATIVE CENSORING FROM EVENTS
  #--------------------------------------------------------------------------
  cat('---------------------------------------------','\n')
  cat('START LMM2','\n')
  cat('---------------------------------------------','\n')
  cat('---------------------------------------------','\n')
  #-----------------------------------------------------
  # CREATE ORDERED FORMULAS FOR LONGITUDINAL MODELS
  #-----------------------------------------------------
  names=""
  if (knot2==1) {
    names=paste(names,"+spline11+spline12+spline13",sep="")
  }
  if (knot2==2) {
    names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23",sep="")
  }
  if (knot2==3) {
    names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33",sep="")
  }
  if (knot2==4) {
    names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43",sep="")
  }
  if (knot2==5) {
    names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43+spline51+spline52+spline53",sep="")
  }
  if (knot2==6) {
    names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43+spline51+spline52+spline53+spline61+spline62+spline63",sep="")
  }
  if (knot2==7) {
    names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43+spline51+spline52+spline53+spline61+spline62+spline63+spline71+spline72+spline73",sep="")
  }
  if (knot2==8) {
    names=paste(names,"+spline11+spline12+spline13+spline21+spline22+spline23+spline31+spline32+spline33+spline41+spline42+spline43+spline51+spline52+spline53+spline61+spline62+spline63+spline71+spline72+spline73+spline81+spline82+spline83",sep="")
  }
  
  names=paste(names,"+nonsurvive+nonsurvive1",sep="")
  names=paste(names,"+(obstime | id)",sep="")
  (formula2=as.formula(paste("bp_curr~baseage",paste(names,collapse="+"))))
  
  censor_time[1]=s
  survive_time[1]=s+tau
  dataset_nonmiss$id=as.numeric(dataset_nonmiss$id)
  dataset$id=as.numeric(dataset$id)
  n=length(dataset$id)
  
  #--------------------------------------------------------------------------
  #  CALCULATE NUMBER OF SUBJECTS AND OBS OF EACH TYPE
  #--------------------------------------------------------------------------
  num_subj_nonsurv=sum(nonsurvive_time <= s+tau)
  num_subj_surv=sum(survive_time >= s+tau)
  num_subj_cens=sum(nonsurvive_time > s+tau & survive_time < s+tau)
  num_obs_nonsurv=sum(nonsurvive_time[dataset_nonmiss$id] <= s+tau)
  num_obs_cens=sum(nonsurvive_time[dataset_nonmiss$id] > s+tau & survive_time[dataset_nonmiss$id] < s+tau)
  num_obs_surv=sum(survive_time[dataset_nonmiss$id] >= s+tau & dataset_nonmiss$obstime <= s+tau)
  
  
  cat('number subjects events prior to s+tau=',num_subj_nonsurv,'\n')
  cat('number subjects censored prior to s+tau=',num_subj_cens,'\n')
  cat('number subjects surviving to s+tau=',num_subj_surv,'\n')
  cat('------------------------------------------------------------------','\n')
  cat('number observations of nonsurvivors=',num_obs_nonsurv,'\n')
  cat('number observations of censored=',num_obs_cens,'\n')
  cat('number past observations of survivors=',num_obs_surv,'\n')
  
  
  cat('--------------------------------------------------','\n')
  cat('Number of subjects (including current) surviving to s+tau=',length(nonsurvive_time[nonsurvive_time > s+tau]),'\n')
  
  if (lmm1_ind==0) {
    #---------------------------------------------------------------------------------------------------
    # FOR EACH ID SURVIVING PAST TIME=S+TAU ADD OBS WITH TIME=S+TAU AND BP_CURR=NA TO DATASET_NONMISS 
    #---------------------------------------------------------------------------------------------------
    i=1
    temp_start=0
    idprev=0
    while (i <= n_nonmiss) {
      #  cat('i=',i,'\n')
      id=dataset_nonmiss$id[i]
      #  cat('id=',id,'\n')
      #  cat('idprev=',idprev,'\n')
      if (id != idprev & (survive_time[id] >= s+tau | id==1) ) {
        if (temp_start==0) {
          temp=dataset_nonmiss[i,]
          temp$bp_curr[1]=NA
          temp$obstime[1]=s+tau
          temp$bptype='ANY'
          temp_start=1
        } else {
          temp1=dataset_nonmiss[i,]
          temp1$bp_curr[1]=NA
          temp1$obstime[1]=s+tau
          temp1$bptype='ANY'
          temp=rbind(temp,temp1) 
        }
      }
      idprev=id
      i=i+1
    }
    # dataset_nonmiss_imp=rbind(dataset_nonmiss,temp)
    dataset_nonmiss_imp=dataset_nonmiss_imp[order(dataset_nonmiss_imp$id,dataset_nonmiss_imp$obstime),]
    n_nonmiss_imp=length(dataset_nonmiss_imp$id)
    
  }
  #cat('Dataset_nonmiss with NAs =','\n')
  #print(dataset_nonmiss_imp[1:20,])
  #cat('dataset =','\n')
  #print(dataset[1:50,1:5])
  #cat('survive_time=','\n')
  #print(survive_time[1:50])
  #cat('censor_time=','\n')
  #print(censor_time[1:50])
  
  #-------------------------------------------------------------
  #  SPLINE CREATION
  #-------------------------------------------------------------
  dataset_nonmiss_imp$spline11=dataset_nonmiss_imp$obstime
  dataset_nonmiss_imp$spline12=dataset_nonmiss_imp$obstime**2
  dataset_nonmiss_imp$spline13=dataset_nonmiss_imp$obstime**3
  if (knot2 >= 2) {
    dataset_nonmiss_imp$spline11=dataset_nonmiss_imp$obstime*(dataset_nonmiss_imp$obstime < knot1) +knot1*(dataset_nonmiss_imp$obstime >= knot1)
    dataset_nonmiss_imp$spline12=dataset_nonmiss_imp$obstime**2*(dataset_nonmiss_imp$obstime < knot1) +knot1**2*(dataset_nonmiss_imp$obstime >= knot1)
    dataset_nonmiss_imp$spline13=dataset_nonmiss_imp$obstime**3*(dataset_nonmiss_imp$obstime < knot1) +knot1**3*(dataset_nonmiss_imp$obstime >= knot1)
    dataset_nonmiss_imp$spline21=(dataset_nonmiss_imp$obstime-knot1)*(dataset_nonmiss_imp$obstime >= knot1) 
    dataset_nonmiss_imp$spline22=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime >= knot1) 
    dataset_nonmiss_imp$spline23=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime >= knot1) 
  }
  if (knot2 >= 3) {
    dataset_nonmiss_imp$spline21=(dataset_nonmiss_imp$obstime-knot1)*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)*(dataset_nonmiss_imp$obstime >= knot2) 
    dataset_nonmiss_imp$spline22=(dataset_nonmiss_imp$obstime-knot1)**2*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**2*(dataset_nonmiss_imp$obstime >= knot2) 
    dataset_nonmiss_imp$spline23=(dataset_nonmiss_imp$obstime-knot1)**3*(dataset_nonmiss_imp$obstime < knot2 & dataset_nonmiss_imp$obstime >= knot1) +(knot2-knot1)**3*(dataset_nonmiss_imp$obstime >= knot2) 
    
    dataset_nonmiss_imp$spline31=(dataset_nonmiss_imp$obstime-knot2)*(dataset_nonmiss_imp$obstime >= knot2) 
    dataset_nonmiss_imp$spline32=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime >= knot2) 
    dataset_nonmiss_imp$spline33=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime >= knot2) 
  }
  if (knot2 >= 4) {
    dataset_nonmiss_imp$spline31=(dataset_nonmiss_imp$obstime-knot2)*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)*(dataset_nonmiss_imp$obstime >= knot3) 
    dataset_nonmiss_imp$spline32=(dataset_nonmiss_imp$obstime-knot2)**2*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**2*(dataset_nonmiss_imp$obstime >= knot3) 
    dataset_nonmiss_imp$spline33=(dataset_nonmiss_imp$obstime-knot2)**3*(dataset_nonmiss_imp$obstime < knot3 & dataset_nonmiss_imp$obstime >= knot2) +(knot3-knot2)**3*(dataset_nonmiss_imp$obstime >= knot3) 
    
    dataset_nonmiss_imp$spline41=(dataset_nonmiss_imp$obstime-knot3)*(dataset_nonmiss_imp$obstime >= knot3) 
    dataset_nonmiss_imp$spline42=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime >= knot3) 
    dataset_nonmiss_imp$spline43=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime >= knot3) 
  }
  if (knot2 >= 5) {
    dataset_nonmiss_imp$spline41=(dataset_nonmiss_imp$obstime-knot3)*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)*(dataset_nonmiss_imp$obstime >= knot4) 
    dataset_nonmiss_imp$spline42=(dataset_nonmiss_imp$obstime-knot3)**2*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**2*(dataset_nonmiss_imp$obstime >= knot4) 
    dataset_nonmiss_imp$spline43=(dataset_nonmiss_imp$obstime-knot3)**3*(dataset_nonmiss_imp$obstime < knot4 & dataset_nonmiss_imp$obstime >= knot3) +(knot4-knot3)**3*(dataset_nonmiss_imp$obstime >= knot4) 
    
    dataset_nonmiss_imp$spline51=(dataset_nonmiss_imp$obstime-knot4)*(dataset_nonmiss_imp$obstime >= knot4) 
    dataset_nonmiss_imp$spline52=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime >= knot4) 
    dataset_nonmiss_imp$spline53=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime >= knot4) 
  }
  if (knot2 >= 6) {
    dataset_nonmiss_imp$spline51=(dataset_nonmiss_imp$obstime-knot4)*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)*(dataset_nonmiss_imp$obstime >= knot5) 
    dataset_nonmiss_imp$spline52=(dataset_nonmiss_imp$obstime-knot4)**2*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**2*(dataset_nonmiss_imp$obstime >= knot5) 
    dataset_nonmiss_imp$spline53=(dataset_nonmiss_imp$obstime-knot4)**3*(dataset_nonmiss_imp$obstime < knot5 & dataset_nonmiss_imp$obstime >= knot4) +(knot5-knot4)**3*(dataset_nonmiss_imp$obstime >= knot5) 
    
    dataset_nonmiss_imp$spline61=(dataset_nonmiss_imp$obstime-knot5)*(dataset_nonmiss_imp$obstime >= knot5) 
    dataset_nonmiss_imp$spline62=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime >= knot5) 
    dataset_nonmiss_imp$spline63=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime >= knot5) 
  }
  if (knot2 >= 7) {
    dataset_nonmiss_imp$spline61=(dataset_nonmiss_imp$obstime-knot5)*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)*(dataset_nonmiss_imp$obstime >= knot6) 
    dataset_nonmiss_imp$spline62=(dataset_nonmiss_imp$obstime-knot5)**2*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**2*(dataset_nonmiss_imp$obstime >= knot6) 
    dataset_nonmiss_imp$spline63=(dataset_nonmiss_imp$obstime-knot5)**3*(dataset_nonmiss_imp$obstime < knot6 & dataset_nonmiss_imp$obstime >= knot5) +(knot6-knot5)**3*(dataset_nonmiss_imp$obstime >= knot6) 
    
    dataset_nonmiss_imp$spline71=(dataset_nonmiss_imp$obstime-knot6)*(dataset_nonmiss_imp$obstime >= knot6) 
    dataset_nonmiss_imp$spline72=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime >= knot6) 
    dataset_nonmiss_imp$spline73=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime >= knot6) 
  }
  if (knot2 >= 8) {
    dataset_nonmiss_imp$spline71=(dataset_nonmiss_imp$obstime-knot6)*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)*(dataset_nonmiss_imp$obstime >= knot7) 
    dataset_nonmiss_imp$spline72=(dataset_nonmiss_imp$obstime-knot6)**2*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**2*(dataset_nonmiss_imp$obstime >= knot7) 
    dataset_nonmiss_imp$spline73=(dataset_nonmiss_imp$obstime-knot6)**3*(dataset_nonmiss_imp$obstime < knot7 & dataset_nonmiss_imp$obstime >= knot6) +(knot7-knot6)**3*(dataset_nonmiss_imp$obstime >= knot7) 
    
    dataset_nonmiss_imp$spline81=(dataset_nonmiss_imp$obstime-knot7)*(dataset_nonmiss_imp$obstime >= knot7) 
    dataset_nonmiss_imp$spline82=(dataset_nonmiss_imp$obstime-knot7)**2*(dataset_nonmiss_imp$obstime >= knot7) 
    dataset_nonmiss_imp$spline83=(dataset_nonmiss_imp$obstime-knot7)**3*(dataset_nonmiss_imp$obstime >= knot7) 
  }
  
  #-------------------------------------------------------------
  #-------------------------------------------------------------
  
  dataset_impute2=dataset
  n=length(dataset_impute2$id)
  dataset_impute2$bp_curr1=rep(NA,times=n)
  dataset_impute2$bp_curr2=rep(NA,times=n)
  dataset_impute2$bp_curr3=rep(NA,times=n)
  dataset_impute2$bp_curr4=rep(NA,times=n)
  dataset_impute2$bp_curr5=rep(NA,times=n)
  
  #-------------------------------------------------
  #  FIT MODEL LMM2
  #-------------------------------------------------
  
  n_nonmiss=length(dataset_nonmiss$id)
  n_nonmiss_imp=length(dataset_nonmiss_imp$id)
  
  cat('--------------------------------------------','\n')
  cat('--------------------------------------------','\n')
  cat('LMM2 model','\n')
  cat('--------------------------------------------','\n')
  
  dataset_nonmiss_imp$nonsurvive=(nonsurvive_time[dataset_nonmiss_imp$id] <= s+tau)
  dataset_nonmiss_imp$nonsurvive1=dataset_nonmiss_imp$nonsurvive*dataset_nonmiss_imp$obstime
  
  select_id=unique(dataset_nonmiss_imp$id)
  n_select_id=length(select_id)
  
  cat('#Subjects included in model =',n_select_id,'\n')
  
  dataset_nonmiss_imp=dataset_nonmiss_imp[dataset_nonmiss_imp$obstime <= s+tau,]
  
  #-------------------------------------------------
  # FIT LMM2
  #-------------------------------------------------
  fit21=try(lmer(formula2,data=subset(dataset_nonmiss_imp, dataset_nonmiss_imp$bptype=='TYPE1'), na.action=na.exclude))
  fit22=try(lmer(formula2,data=subset(dataset_nonmiss_imp, dataset_nonmiss_imp$bptype=='TYPE2'), na.action=na.exclude))
  fit23=try(lmer(formula2,data=subset(dataset_nonmiss_imp, dataset_nonmiss_imp$bptype=='TYPE3'), na.action=na.exclude))
  cat('Model LMM2=','\n')
  print(summary(fit21))
  print(summary(fit22))
  print(summary(fit23))
  
  fitted21_na=predict(fit21,newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr),])
  fitted22_na=predict(fit22,newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr),])
  fitted23_na=predict(fit23,newdata=dataset_nonmiss_imp[is.na(dataset_nonmiss_imp$bp_curr),])
  
  #cat('--------------------------------------------------','\n')
  #cat('fitted1_na=','\n')
  #print(fitted1_na[1:20])
  #cat('Dataset before imputed BP_curr =','\n')
  #print(dataset_impute[1:40,])
  
  # All Patients
  # Replace ALL values of Biomarker Process with fitted value from model
  dataset_impute2$bp_curr1[survive_time[dataset_impute2$id] >= s+tau | dataset_impute2$id==1]=fitted21_na
  dataset_impute2$bp_curr2[survive_time[dataset_impute2$id] >= s+tau | dataset_impute2$id==1]=fitted22_na
  dataset_impute2$bp_curr3[survive_time[dataset_impute2$id] >= s+tau | dataset_impute2$id==1]=fitted23_na
  
  #cat('Dataset after imputed BP_curr =','\n')
  #print(dataset_impute[1:40,])
  
  cat('END LMM2 model','\n')
  cat('--------------------------------------------','\n')
  
  
  lmm2_pred=cbind(TYPE1=dataset_impute2$bp_curr1[1], TYPE2=dataset_impute2$bp_curr2[1], TYPE3=dataset_impute2$bp_curr3[1])
}

}
# END: LMM1_IND==1 | LMM2_IND==1
