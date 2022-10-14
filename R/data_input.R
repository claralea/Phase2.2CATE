library(dplyr)
library(data.table)
library(tidyr)



# input.path: character, path to input data
# Y.type : character, outcome type ('severe', 'icu', 'dead')

#' @import dplyr
#' @import data.table
#' @import tidyr

create.table = function(input.path){
  #### Read the data
  dat.loc.race = fread(paste0(input.path, '/LocalPatientRace.csv'))
  dat.loc.sum = fread(paste0(input.path, '/LocalPatientSummary.csv'))
  dat.loc.patobs = fread(paste0(input.path, '/LocalPatientObservations.csv'))
  dat.loc.vac = fread(paste0(input.path, '/LocalPatientVaccine.csv'))

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
  icd.chd = c('I20', 'I21', 'I24', 'I25.10', 'I25.110', 'I25.2', 'I25.3', 'I25.41', 'I25.42', 'I25.5', 'I25.700', 'I25.710', 'I25.720', 'I25.730', 'I25.750', 'I25.760', 'I25.790', 'I25.810', 'I25.811', 'I25.812', 'I25.82', 'I25.83', 'I25.84', 'I25.89', 'I25.9')# CHD
  icd.hf = c('I11.0', 'I13.0', 'I13.2', 'I50.20', 'I50.21', 'I50.22', 'I50.23', 'I50.30', 'I50.31', 'I50.32', 'I50.33', 'I50.40', 'I50.41', 'I50.42', 'I50.43', 'I50.814', 'I50.9', 'I50.1', 'I50.810', 'I50.811', 'I50.812', 'I50.813', 'I50.82', 'I50.83', 'I50.84', 'I50.89')# HF
  icd.hypert = c('I10', 'I11', 'I12', 'I13', 'I15', 'I16')# Hypertension
  icd.diabetes = c('E08', 'E10', 'E11', 'E13')# Diabetes
  icd.ckd = c('I12.0', 'I13.1', 'N03.2', 'N03.3', 'N03.4', 'N03.5', 'N03.6', 'N03.7', 'N18', 'N19', 'N05.2', 'N05.3', 'N05.4', 'N05.5', 'N05.6', 'N05.7', 'N25.0', 'Z49.0', 'Z49.1', 'Z49.2', 'Z94.0', 'Z99.2')# CKD
  icd.cancer = c(paste0('C0', c(0:9)), paste0('C', c(10:41)), 'C43', 'C45', paste0('C', c(47:86)), 'C88', paste0('C', c(90:96)), 'C7A', 'C7B')# Cancer

  dat.loc.obs = left_join(dat.loc.patobs, dat.loc.sum[,c('cohort', 'patient_num', 'admission_date')], by = c('cohort', 'patient_num'))

  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 3) %in% icd.asthma] = 'asthma'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 3) %in% icd.bronch] = 'bronchitis'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 4) %in% icd.bronch] = 'copd'
  dat.loc.obs$cov[substr(dat.loc.obs$concept_code, 1, 3) %in% icd.ami] = 'ami'
  dat.loc.obs$cov[dat.loc.obs$concept_code %like% paste(icd.chd, collapse="|") ] = 'chd'
  dat.loc.obs$cov[dat.loc.obs$concept_code %like% paste(icd.hf, collapse="|") ] = 'hf'
  dat.loc.obs$cov[dat.loc.obs$concept_code %like% paste(icd.hypert, collapse="|") ] = 'hypertension'
  dat.loc.obs$cov[dat.loc.obs$concept_code %like% paste(icd.diabetes, collapse="|") ] = 'diabetes'
  dat.loc.obs$cov[dat.loc.obs$concept_code %like% paste(icd.ckd, collapse="|") ] = 'ckd'
  dat.loc.obs$cov[dat.loc.obs$concept_code %like% paste(icd.cancer, collapse="|") ] = 'cancer'

  dat.loc.obs$cov[dat.loc.obs$concept_code %like% "I21.A" ] = NA

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

  #### Select outcome info
  ## 9 outcomes: svere case, hospitalization, death.
  ## binary indicators for the event within 3 months, 6 months, and 9 months of the first dose.

  dat.out = dat.loc.sum[,c('cohort', 'patient_num', 'severe_date', 'severe', 'icu_date', 'icu', 'death_date', 'dead')]
  dat.out = left_join(dat.out, dat.loc.vac[,c('cohort', 'patient_num', 'vaccine_date', 'vaccine_type')], by = c('cohort', 'patient_num'))

  dat.out = dat.out %>% filter(vaccine_type %in% c('Moderna', 'Pfizer'))

  dat.out$severe_days = as.numeric(difftime(dat.out$severe_date, dat.out$vaccine_date, units = 'days'))
  dat.out$icu_days = as.numeric(difftime(dat.out$icu_date, dat.out$vaccine_date, units = 'days'))
  dat.out$dead_days = as.numeric(difftime(dat.out$death_date, dat.out$vaccine_date, units = 'days'))

  dat.out$severe_3[dat.out$severe_days < 90 & dat.out$severe_days >=0] = 1
  dat.out$severe_6[dat.out$severe_days < 180 & dat.out$severe_days >=0] = 1
  dat.out$severe_9[dat.out$severe_days < 270 & dat.out$severe_days >=0] = 1

  dat.out$icu_3[dat.out$icu_days < 90 & dat.out$icu_days >=0] = 1
  dat.out$icu_6[dat.out$icu_days < 180 & dat.out$icu_days >=0] = 1
  dat.out$icu_9[dat.out$icu_days < 270 & dat.out$icu_days >=0] = 1

  dat.out$dead_3[dat.out$dead_days < 90 & dat.out$dead_days >=0] = 1
  dat.out$dead_6[dat.out$dead_days < 180 & dat.out$dead_days >=0] = 1
  dat.out$dead_9[dat.out$dead_days < 270 & dat.out$dead_days >=0] = 1

  dat.out.save = dat.out %>% select(cohort, patient_num, severe_3, severe_6, severe_9, icu_3, icu_6, icu_9, dead_3, dead_6, dead_9)
  dat.out.save[is.na(dat.out.save)] = 0

  #### combine the datasets
  dat.input = left_join(dat.imp.save, dat.cov.save, by = c('cohort', 'patient_num'))
  dat.input[is.na(dat.input)] = 0 ## 0 when don't have any cov
  dat.input = left_join(dat.input, dat.out.save, by = c('cohort', 'patient_num'))

  return(list(X = dat.input[, c(4:17)],
              A = dat.input[, 3],
              Y.all = dat.input[, c(18:26)]))
}






