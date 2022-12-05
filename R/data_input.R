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
  #input.path="/Users/harrisonzhang/Dropbox (HMS)/4CE/Vaccine Study/Phase2.2CATE-main/R/"
  dat.loc.race = fread(paste0(input.path, '/LocalPatientRace.csv'))
  dat.loc.sum = fread(paste0(input.path, '/LocalPatientSummary.csv'))
  dat.loc.patobs = fread(paste0(input.path, '/LocalPatientObservations.csv'))
  dat.loc.vac = fread(paste0(input.path, '/LocalPatientVaccine.csv'))
  
  #### Map observatons to PheCodes
  dat.icd=dplyr::filter(dat.loc.patobs,concept_type%in%c("DIAG-ICD10","DIAG-ICD9"))
  dat.phecode=left_join(dat.icd, icd10.phecode.map[,c("concept_code","phecode","description")],
                    by="concept_code")
  dat.phecode=dplyr::filter(dat.phecode,!is.na(dat.phecode$phecode))
  
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
  dat.loc.obs = left_join(dat.phecode, 
                          dat.loc.sum[,c('cohort', 'patient_num', 'admission_date')], 
                          by = c('cohort', 'patient_num'))
  dat.loc.obs=dat.loc.obs[!duplicated(dat.loc.obs),]
  ### Identify relevant covariates of interest
  cov.str=c("asthma","bronchitis","emphysema","hypertension",
            "CKD","heart failure","obesity")
  icd.cancer = c(paste0('C0', c(0:9)), paste0('C', c(10:41)), 'C43', 'C45', paste0('C', c(47:86)), 'C88', paste0('C', c(90:96)), 'C7A', 'C7B')
  phecode.cancer = icd10.phecode.map$phecode[icd10.phecode.map$concept_code %in% icd.cancer];phecode.cancer=unique(phecode.cancer)
  phecode.chd="411"
  phecode.t2d="250.2"
  phecode.copd="496"
  dat.loc.obs = dat.loc.obs %>%
    dplyr::filter(grepl(paste0(cov.str,collapse = "|"),
                        description,
                        ignore.case = T) | phecode %in% c(phecode.cancer,phecode.chd,phecode.t2d,phecode.copd))
  
  ### HARRISON: remember to correct for COPD using Emphysema and chronic bronchitis 496.2 496.21 496.1
  dat.loc.obs$cov = dat.loc.obs$description
  dat.loc.obs$cov = tolower(dat.loc.obs$cov)
  dat.loc.obs$cov[grepl("asthma",dat.loc.obs$cov,ignore.case = T)]="asthma"
  dat.loc.obs$cov[grepl("bronchitis",dat.loc.obs$cov,ignore.case = T)]="bronchitis"
  dat.loc.obs$cov[dat.loc.obs$phecode %in% c("496.2","496.21","496.1","496")]="copd"
  dat.loc.obs$cov[dat.loc.obs$phecode %in% phecode.chd]="chd"
  dat.loc.obs$cov[grepl("heart failure",dat.loc.obs$cov,ignore.case = T)]="hf"
  dat.loc.obs$cov[grepl("hypertension",dat.loc.obs$cov,ignore.case = T)]="hypertension"
  dat.loc.obs$cov[dat.loc.obs$phecode %in% phecode.t2d]="diabetes"
  dat.loc.obs$cov[grepl("CKD",dat.loc.obs$cov,ignore.case = T)]="ckd"
  dat.loc.obs$cov[dat.loc.obs$phecode %in% phecode.cancer]="cancer"
  dat.loc.obs$cov[grepl("obesity",dat.loc.obs$cov,ignore.case = T)]="obesity"
  
  dat.loc.obs = dat.loc.obs %>%
    dplyr::filter(cov %in% c("asthma","bronchitis","copd","chd","hf","hypertension","diabetes","ckd","cancer","obesity"))
  
  ## delete rows with other codes
  #dat.loc.obs = na.omit(dat.loc.obs)
  dat.loc.obs = left_join(dat.loc.obs, dat.loc.vac[,c('cohort', 'patient_num', 'vaccine_date', 'vaccine_type')], by = c('cohort', 'patient_num'))
  
  ## delete patients with no vaccine
  dat.loc.obs = dat.loc.obs %>% filter(vaccine_type %in% c('Moderna', 'Pfizer'))
  
  ## create binary code
  dat.loc.obs$code_date = as.Date(dat.loc.obs$admission_date) + dat.loc.obs$days_since_admission
  dat.loc.obs$code_days = as.numeric(difftime(dat.loc.obs$vaccine_date, dat.loc.obs$code_date, units = 'days'))
  dat.loc.obs = dat.loc.obs %>% dplyr::select(cohort, patient_num, cov, code_days) %>%
    dplyr::filter(code_days > 0) %>% dplyr::group_by(cohort, patient_num, cov) %>% dplyr::arrange(code_days) %>% dplyr::slice(1)
  
  dat.loc.obs$cov_index = 1
  #dat.loc.obs$cov_index[dat.loc.obs$code_days < 365] = 1
  #dat.loc.obs$cov_index[dat.loc.obs$code_days >= 365] = 0
  
  dat.cov.save = tidyr::spread(dat.loc.obs[,c('cohort', 'patient_num', 'cov', 'cov_index')], key = cov, value = cov_index, fill = 0)
  
  #### combine the datasets
  dat.input = left_join(dat.imp.save, dat.cov.save, by = c('cohort', 'patient_num'))
  dat.input[is.na(dat.input)] = 0 ## 0 when don't have any cov
  dat.input = left_join(dat.input, res, by = c('patient_num'))
  
  dat.input = na.omit(dat.input)

  # Summary statistics for the processed data
  colMeans(dat.input[,c((3:31))])
  
  # make data split reproducible
  set.seed(2022)
  
  # use 20% of dataset for tree and 80% for FACE
  train = dat.input %>% dplyr::sample_frac(0.20)
  test  = dplyr::anti_join(dat.input, train, by = 'patient_num')
  
  write.csv(train, file = paste0(output.path, '/train_', siteid, '.csv'), row.names = F)
  write.csv(test, file = paste0(output.path, '/test_', siteid, '.csv'), row.names = F)
  
  
  return(list(X.train = train[, c("gender_male","age_50","age_70","age_80","race",
                                  "asthma","bronchitis","copd","chd","hf","hypertension","diabetes","ckd","cancer","obesity",
                                  "pre_vaccine_infection")],
              A.train = train[, "A"],
              Y.train = train[, c("post_vaccine_infection_3month","post_vaccine_infection_6month",
                                  "post_vaccine_infection_9month","post_vaccine_infection_999","hosp_3month",
                                  "hosp_6month","hosp_9month", "hosp_999","death_3month","death_6month",
                                  "death_9month","death_999")],
              X.test = test[,c("gender_male","age_50","age_70","age_80","race",
                               "asthma","bronchitis","copd","chd","hf","hypertension","diabetes","ckd","cancer","obesity",
                               "pre_vaccine_infection")],
              A.test = test[, "A"],
              Y.test = test[, c("post_vaccine_infection_3month","post_vaccine_infection_6month",
                                "post_vaccine_infection_9month","post_vaccine_infection_999","hosp_3month",
                                "hosp_6month","hosp_9month", "hosp_999","death_3month","death_6month",
                                "death_9month","death_999")]))
  
}


