
construct_outcomes = function(summary,
                              obs,
                              vax){
  
  summary=data.frame(summary)
  obs=data.frame(obs)
  vax=data.frame(vax)
  vax[,"vaccine_date"]=as.Date(vax[,"vaccine_date"])
  summary[,"admission_date"]=as.Date(summary[,"admission_date"])
  
  ### Clean vax data
  vax = vax %>%
    group_by(patient_num)%>%
    slice_min(vaccine_date) %>%
    ungroup()
  vax=data.frame(vax)
  
  ### Only select patients who have >=1 Pfizer or Moderna dose
  vax=vax%>%
    filter(vaccine_type %in% c("Pfizer","Moderna"))
  pat.keep=vax[,"patient_num"]
  obs=obs%>%
    filter(patient_num %in% pat.keep)
  summary=summary%>%
    filter(patient_num %in% pat.keep)
  
  obs=left_join(obs,
                dplyr::select(summary,patient_num,admission_date,death_date),
                by="patient_num")
  
  obs = obs %>%
    mutate("obs_date"=admission_date+days_since_admission)
  
  obs=left_join(obs,
                dplyr::select(vax,patient_num,vaccine_date),
                by="patient_num")
  
  obs = obs %>%
    mutate("days_since_vaccine"=as.numeric(obs_date-vaccine_date))
  
  ### Construct pre-vax infection flag
  tmp=obs%>%
    filter(concept_code %in% c("covidpos",'U07.1'))%>%
    mutate("inf_pre"=obs_date<vaccine_date)
  
  pat.inf.pre=unique(tmp$patient_num[tmp$inf_pre==TRUE])
  res=data.frame("patient_num"=pat.keep,
                 "pre_vaccine_infection"=ifelse(pat.keep %in% pat.inf.pre,1,0))
  
  
  ### Construct post-vax infection flag
  tmp=obs%>%
    filter(concept_code %in% c("covidpos",'U07.1'),
           days_since_vaccine>0)%>%
    group_by(patient_num)%>%
    slice_min(days_since_vaccine)%>%
    ungroup()
  tmp=data.frame(tmp)
  
  pat.inf.post90=unique(tmp$patient_num[tmp$days_since_vaccine<=90])
  pat.inf.post180=unique(tmp$patient_num[tmp$days_since_vaccine<=180])
  pat.inf.post270=unique(tmp$patient_num[tmp$days_since_vaccine<=270])
  pat.inf.post999=unique(tmp$patient_num[tmp$days_since_vaccine<=999])
  
  
  res=cbind.data.frame(res,
                       "post_vaccine_infection_3month"=ifelse(pat.keep %in% pat.inf.post90,1,0),
                       "post_vaccine_infection_6month"=ifelse(pat.keep %in% pat.inf.post180,1,0),
                       "post_vaccine_infection_9month"=ifelse(pat.keep %in% pat.inf.post270,1,0),
                       "post_vaccine_infection_999"=ifelse(pat.keep %in% pat.inf.post999,1,0))
  
  
  ### Construct hospitalization flag among those infected after vaccine
  tmp=obs%>%
    filter(concept_code %in% c("covidpos",'U07.1'),
           days_since_vaccine>0,
           grepl("PosAdm|U071Adm",cohort))%>%
    mutate("days_from_vaccine_to_hosp"=as.numeric(admission_date-vaccine_date))%>%
    filter(days_from_vaccine_to_hosp>0)
  
  tmp=data.frame(tmp)
  
  pat.hosp.post90=unique(tmp$patient_num[tmp$days_from_vaccine_to_hosp<=90])
  pat.hosp.post180=unique(tmp$patient_num[tmp$days_from_vaccine_to_hosp<=180])
  pat.hosp.post270=unique(tmp$patient_num[tmp$days_from_vaccine_to_hosp<=270])
  pat.hosp.post999=unique(tmp$patient_num[tmp$days_from_vaccine_to_hosp<=999])
  
  
  res=cbind.data.frame(res,
                       "hosp_3month"=ifelse(pat.keep %in% pat.hosp.post90,1,0),
                       "hosp_6month"=ifelse(pat.keep %in% pat.hosp.post180,1,0),
                       "hosp_9month"=ifelse(pat.keep %in% pat.hosp.post270,1,0),
                       "hosp_999"=ifelse(pat.keep %in% pat.hosp.post999,1,0))
  
  ### Construct death flag among those infected after vaccine
  tmp=obs%>%
    filter(concept_code %in% c("covidpos",'U07.1'),
           days_since_vaccine>0,
           as.numeric(as.Date(death_date)-obs_date)>0)%>%
    mutate("days_from_vaccine_to_death"=as.numeric(as.Date(death_date)-vaccine_date))%>%
    filter(days_from_vaccine_to_death>0)
  
  tmp=data.frame(tmp)
  
  pat.death.post90=unique(tmp$patient_num[tmp$days_from_vaccine_to_death<=90])
  pat.death.post180=unique(tmp$patient_num[tmp$days_from_vaccine_to_death<=180])
  pat.death.post270=unique(tmp$patient_num[tmp$days_from_vaccine_to_death<=270])
  pat.death.post999=unique(tmp$patient_num[tmp$days_from_vaccine_to_death<=999])
  
  
  res=cbind.data.frame(res,
                       "death_3month"=ifelse(pat.keep %in% pat.death.post90,1,0),
                       "death_6month"=ifelse(pat.keep %in% pat.death.post180,1,0),
                       "death_9month"=ifelse(pat.keep %in% pat.death.post270,1,0),
                       "death_999"=ifelse(pat.keep %in% pat.death.post999,1,0))
  
  
  return(res)
  
  
}




create.table = function(input.path){
  #### Read the data
  dat.loc.race = fread(paste0(input.path, '/LocalPatientRace.csv'))
  dat.loc.sum = fread(paste0(input.path, '/LocalPatientSummary.csv'))
  dat.loc.patobs = fread(paste0(input.path, '/LocalPatientObservations.csv'))
  dat.loc.vac = fread(paste0(input.path, '/LocalPatientVaccine.csv'))
  
  # Create outcomes
  res = construct_outcomes(dat.loc.sum, dat.loc.patobs, dat.loc.vac)
  
  #### Select unique patients and vaccine info / Only keep Moderna(1) and Pfizer(0). Save date for now as t0
  dat.imp = dat.loc.vac %>% filter(vaccine_type %in% c('Moderna', 'Pfizer'))
  dat.imp$A[which(dat.imp$vaccine_type == 'Moderna')] = 1
  dat.imp$A[which(dat.imp$vaccine_type == 'Pfizer')] = 0
  dat.imp$vaccine_date = substr(as.character(dat.imp$vaccine_date), 1, 7)
  
  #### Select demographics info / Age, gender, race
  dat.imp = left_join(dat.imp, dat.loc.sum[, c('cohort', 'patient_num', 'age', 'sex')], by = c('cohort', 'patient_num'))
  dat.imp$age_50 = dat.imp$age_70 = dat.imp$age_80 = 0
  dat.imp$age_50[dat.imp$age<50] = 1
  dat.imp$age_70[dat.imp$age<70] = 1
  dat.imp$age_80[dat.imp$age<80] = 1
  
  dat.imp$gender_male[dat.imp$sex == 'male'] = 1
  
  dat.loc.race$race = 0
  dat.loc.race$race[dat.loc.race$race_4ce == 'white'] = 1
  
  dat.imp.save = left_join(dat.imp[, c('cohort', 'patient_num','A', 'gender_male', 'age_50', 'age_70', 'age_80')],
                           dat.loc.race[, c('cohort', 'patient_num', 'race')], by = c('cohort', 'patient_num'))
  
  #### Select covariates info
  icd.asthma = c('J45')# Asthma
  icd.bronch = c('J20', 'J21', 'J40', 'J41', 'J42')# Bronchitis
  icd.copd = c('J41.0', 'J41.1', 'J41.8', 'J42', 'J43.0', 'J43.1','J43.2', 'J43.8', 'J43.9', 'J44.0', 'J44.1', 'J44.9')# COPD
  icd.ami = c('I21, I22')# AMI
  icd.chd7 = c('I25.110','I25.700', 'I25.710', 'I25.720', 'I25.730','I25.750', 'I25.760', 'I25.790',
               'I25.810','I25.811', 'I25.812')
  icd.chd6 = c('I20', 'I21', 'I24', 'I25.10', 'I25.41', 'I25.42', 'I25.82', 'I25.83', 'I25.84', 'I25.89')# CHD
  icd.chd5 = c('I25.2', 'I25.3','I25.5', 'I25.9')
  icd.hf5 = c('I11.0', 'I13.0', 'I13.2', 'I50.9', 'I50.1')# HF
  icd.hf6 = c('I50.20', 'I50.21', 'I50.22', 'I50.23', 'I50.30', 'I50.31', 'I50.32',
              'I50.33', 'I50.40', 'I50.41', 'I50.42', 'I50.43', 'I50.82', 'I50.83', 'I50.84', 'I50.89')
  icd.hf7 = c('I50.814','I50.810', 'I50.811', 'I50.812', 'I50.813')
  icd.hypert = c('I10', 'I11', 'I12', 'I13', 'I15', 'I16')# Hypertension
  icd.hypert5 = c('I11.0','I10.0','I11.9','I12.0','I12.9','I13.0','I13.1','I13.9','I15.0','I15.1','I15.9')
  icd.diabetes = c('E11.9','E11.8','E11.7')# Diabetes T2
  icd.ckd = c('I12.0', 'I13.1', 'N03.2', 'N03.3', 'N03.4', 'N03.5', 'N03.6', 'N03.7', 'N05.2', 'N05.3',
              'N05.4', 'N05.5', 'N05.6', 'N05.7', 'N25.0',
              'N18.1','N18.2', 'N18.3', 'N18.4', 'N18.5','N18.9', 'Z49.0', 'Z49.1', 'Z49.2', 'Z94.0', 'Z99.2')# CKD
  icd.cancer = c(paste0('C0', c(0:9)), paste0('C', c(10:41)), 'C43', 'C45', paste0('C', c(47:86)), 'C88', paste0('C', c(90:96)), 'C7A', 'C7B')# Cancer
  icd.covid = c('U07.1')
  dat.loc.obs = left_join(dat.loc.patobs, 
                          dat.loc.sum[,c('cohort', 'patient_num', 'admission_date')], 
                          by = c('cohort', 'patient_num'))
  # dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 5) %in% icd.covid] = 'infected'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 3) %in% icd.asthma] = 'asthma'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 3) %in% icd.bronch] = 'bronchitis'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 5) %in% icd.copd] = 'copd'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 3) %in% icd.ami] = 'ami'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 7) %in% icd.chd7] = 'chd'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 6) %in% icd.chd6] = 'chd'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 5) %in% icd.chd5] = 'chd'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 7) %in% icd.hf7] = 'hf'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 6) %in% icd.hf6] = 'hf'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 5) %in% icd.hf5] = 'hf'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 3) %in% icd.hypert] = 'hypertension'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 5) %in% icd.hypert5] = 'hypertension'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 5) %in% icd.diabetes] = 'diabetes'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 5) %in% icd.ckd] = 'ckd'
  dat.loc.obs$cov[dat.loc.obs$concept_code %like% paste(icd.cancer, collapse="|") ] = 'cancer'
  
  
  ## delete rows with  othe codes
  dat.loc.obs = na.omit(dat.loc.obs)
  dat.loc.obs = left_join(dat.loc.obs, dat.loc.vac[,c('cohort', 'patient_num', 'vaccine_date', 'vaccine_type')], by = c('cohort', 'patient_num'))
  
  ## delete patients with no vaccine
  dat.loc.obs = dat.loc.obs %>% filter(vaccine_type %in% c('Moderna', 'Pfizer'))
  
  ## create binary code
  dat.loc.obs$code_date = dat.loc.obs$admission_date + dat.loc.obs$days_since_admission
  dat.loc.obs$code_days = as.numeric(difftime(dat.loc.obs$vaccine_date, dat.loc.obs$code_date, units = 'days'))
  dat.loc.obs = dat.loc.obs %>% select(cohort, patient_num, cov, code_days) %>%
    filter(code_days > 0) %>% group_by(cohort, patient_num, cov) %>% arrange(code_days) %>% slice(1)
  
  dat.loc.obs$cov_index = 0
  dat.loc.obs$cov_index[dat.loc.obs$code_days < 365] = 1
  dat.loc.obs$cov_index[dat.loc.obs$code_days >= 365] = 0
  
  dat.cov.save = spread(dat.loc.obs[,c('cohort', 'patient_num', 'cov', 'cov_index')], key = cov, value = cov_index, fill = 0)
  
 
  
  #### combine the datasets
  dat.input = left_join(dat.imp.save, dat.cov.save, by = c('cohort', 'patient_num'))
  dat.input[is.na(dat.input)] = 0 ## 0 when don't have any cov
  dat.input = left_join(dat.input, res, by = c('patient_num'))
  
  # make data split reproducible
  set.seed(2022)
  
  # use 20% of dataset for tree and 80% for FACE
  train = dat.input %>% dplyr::sample_frac(0.20)
  test  = dplyr::anti_join(dat.input, train, by = 'patient_num')
  
  write.csv(train, file = paste0(output.path, '/train_', siteid, '.csv'), row.names = F)
  write.csv(test, file = paste0(output.path, '/test_', siteid, '.csv'), row.names = F)
  
  
  return(list(X.train = train[, c(4:18)],
              A.train = train[, 3],
              Y.train = train[, c(19:30)],
              X.test = test[,c(4:18)],
              A.test = test[, 3],
              Y.test = test[, c(19:30)]))
  
}


