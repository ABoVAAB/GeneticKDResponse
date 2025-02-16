cat("\014")  
rm(list = ls()) 
gc()  

invisible(source("Load GeneticKDResponse Data.R"))

packages <- c("dplyr", "tibble", "table1", "stats", "MatrixModels", "binom", "glmglrt")
invisible(lapply(packages, library, character.only = TRUE))

# Table 1
table1(~ age_at_kd + sex + age_at_epilepsy_onset + epilepsy_type + num_seizure_types + num_previous_ASMs + num_concomitant_ASMs + development + motor_function + gastrostomy + kd_ratio_3mo + bohb_3mo + seizure_type_fas_fias + seizure_type_fbtcs + seizure_type_spasms + seizure_type_gtc + seizure_type_tonic + seizure_type_atonic + seizure_type_myoclonic + response_3mo + response_6mo + response_12mo + response_24mo + genetic_group + functional_group + gene_name + variant | response_3mo_50, data)
table1(~ age_at_kd + sex + age_at_epilepsy_onset + epilepsy_type + num_seizure_types + num_previous_ASMs + num_concomitant_ASMs + development + motor_function + gastrostomy + kd_ratio_3mo + bohb_3mo + seizure_type_fas_fias + seizure_type_fbtcs + seizure_type_spasms + seizure_type_gtc + seizure_type_tonic + seizure_type_atonic + seizure_type_myoclonic + response_3mo + response_6mo + response_12mo + response_24mo + genetic_group + functional_group + gene_name + variant | response_3mo_90, data)
table1(~ age_at_kd + sex + age_at_epilepsy_onset + epilepsy_type + num_seizure_types + num_previous_ASMs + num_concomitant_ASMs + development + motor_function + gastrostomy + kd_ratio_3mo + bohb_3mo + seizure_type_fas_fias + seizure_type_fbtcs + seizure_type_spasms + seizure_type_gtc + seizure_type_tonic + seizure_type_atonic + seizure_type_myoclonic + response_3mo + response_6mo + response_12mo + response_24mo + genetic_group + functional_group + gene_name + variant | response_3mo_100, data)
table1(~ age_at_kd + sex + age_at_epilepsy_onset + epilepsy_type + num_seizure_types + num_previous_ASMs + num_concomitant_ASMs + development + motor_function + gastrostomy + kd_ratio_3mo + bohb_3mo + seizure_type_fas_fias + seizure_type_fbtcs + seizure_type_spasms + seizure_type_gtc + seizure_type_tonic + seizure_type_atonic + seizure_type_myoclonic + response_3mo + response_6mo + response_12mo + response_24mo + genetic_group + functional_group + gene_name + variant | genetic_group, data)

#table 2
invisible(source("logbin.R"))

prev <- function(y_var, x_var, data) {
  
  if (!(y_var %in% names(data))) stop(paste("Variable", y_var, "not found in data"))
  if (!(x_var %in% names(data))) stop(paste("Variable", x_var, "not found in data"))
  
  tab <- xtabs(as.formula(paste(" ~ ", y_var, " + ", x_var)), data = data)
  
  binom_results <- binom.profile(t(tab)[, 2], t(tab)[, 1] + t(tab)[, 2])
  
  formatted_results <- paste0(
    sprintf("%3.1f", binom_results$mean * 100), " (", 
    sprintf("%3.1f", binom_results$lower * 100), "-", 
    sprintf("%3.1f", binom_results$upper * 100), ")"
  )
  
  out <- cbind(colnames(tab), binom_results, formatted_results)
  colnames(out)[ncol(out)] <- "Prevalence (%)"
  
  return(out)
}

d0<-cbind(data,
          model.matrix(~ sex - 1, data = data),
          model.matrix(~ epilepsy_type - 1, data = data),
          model.matrix(~ development - 1, data = data),
          model.matrix(~ motor_function - 1, data = data),
          model.matrix(~ addNA(gastrostomy) - 1, data = data),
          model.matrix(~ addNA(seizure_type_fas_fias) - 1, data = data),
          model.matrix(~ addNA(seizure_type_fbtcs) - 1, data = data),
          model.matrix(~ addNA(seizure_type_spasms) - 1, data = data),
          model.matrix(~ addNA(seizure_type_gtc) - 1, data = data),
          model.matrix(~ addNA(seizure_type_tonic) - 1, data = data),
          model.matrix(~ addNA(seizure_type_atonic) - 1, data = data),
          model.matrix(~ addNA(seizure_type_myoclonic) - 1, data = data),
          model.matrix(~ addNA(genetic_group) - 1, data = data),
          model.matrix(~ functional_group - 1, data = data),
          model.matrix(~ addNA(gene_name) - 1, data = data),
          model.matrix(~ addNA(variant) - 1, data = data),
          model.matrix(~ response_3mo_50 - 1, data = data),
          model.matrix(~ response_3mo_90 - 1, data = data),
          model.matrix(~ response_3mo_100 - 1, data = data))

d0$functional_group_1_4 <- d0$"functional_groupIon Channel"+d0$"functional_groupNeurotransmitter receptors"+d0$"functional_groupTransporters"+d0$"functional_groupSynapse-related"           
d0$functional_group_5_6 <- d0$"functional_groupCell growth, division, and proliferation-related"+d0$"functional_groupCell structural integrity and or homeostasis"    
d0$functional_group_7_8 <- d0$"functional_groupMetabolism-related"+d0$"functional_groupEnergy homeostasis and mitochondrial function"

logbin("response_3mo_50",c("age_at_kd"),d0)$table

prev("response_3mo_50","sex",d0)
logbin("response_3mo_50",c("sexFemale"),d0)$table

logbin("response_3mo_50",c("age_at_epilepsy_onset"),d0)$table

prev("response_3mo_50","epilepsy_type",d0)
logbin("response_3mo_50",c("epilepsy_typeFocal", "epilepsy_typeCombined"),d0)$table

logbin("response_3mo_50",c("num_seizure_types"),d0)$table

prev("response_3mo_50","seizure_type_fas_fias",d0)
prev("response_3mo_50","seizure_type_fbtcs",d0)
prev("response_3mo_50","seizure_type_spasms",d0)
prev("response_3mo_50","seizure_type_gtc",d0)
prev("response_3mo_50","seizure_type_tonic",d0)
prev("response_3mo_50","seizure_type_atonic",d0)
prev("response_3mo_50","seizure_type_myoclonic",d0)
logbin("response_3mo_50",c("addNA(seizure_type_fas_fias)Yes", "addNA(seizure_type_fbtcs)Yes", "addNA(seizure_type_spasms)Yes", "addNA(seizure_type_gtc)Yes", "addNA(seizure_type_tonic)Yes", "addNA(seizure_type_atonic)Yes", "addNA(seizure_type_myoclonic)Yes"),d0)$table

logbin("response_3mo_50",c("num_previous_ASMs"),d0)$table

logbin("response_3mo_50",c("num_concomitant_ASMs"),d0)$table

prev("response_3mo_50","gastrostomy",d0)
logbin("response_3mo_50",c("addNA(gastrostomy)Yes", "addNA(gastrostomy)NA"),d0)$table

prev("response_3mo_50","development",d0)
logbin("response_3mo_50",c("developmentID","developmentnot_defined"),d0)$table

prev("response_3mo_50","motor_function",d0)
logbin("response_3mo_50",c("motor_functionmotor_dysfunction", "motor_functionnot_defined"),d0)$table

logbin("response_3mo_50",c("kd_ratio_3mo"),d0)$table

logbin("response_3mo_50",c("bohb_3mo"),d0)$table

logbin("response_3mo_50",c("age_at_epilepsy_onset", "epilepsy_typeFocal", "epilepsy_typeCombined", "addNA(seizure_type_fas_fias)Yes", "addNA(seizure_type_fbtcs)Yes", "addNA(seizure_type_spasms)Yes", "addNA(seizure_type_gtc)Yes", "addNA(seizure_type_tonic)Yes", "addNA(seizure_type_atonic)Yes", "addNA(seizure_type_myoclonic)Yes", "num_previous_ASMs"),d0)$table

#Table 4

d0$fixed<-1

prev("response_3mo_50","fixed",d0)
timepoints <- c("3mo", "6mo", "12mo", "24mo")
thresholds <- c("50", "90", "100")

for (t in timepoints) {
  for (th in thresholds) {
    response_var <- paste0("response_", t, "_", th)
    print(response_var)
    print(prev(response_var,"fixed",d0))
    print(prev(response_var,"genetic_group",d0))
    cat("\n")
  }
}

#table 3
d0<-d0[d0$"addNA(genetic_group)NA"!=1,]
d0<-d0[d0$"addNA(genetic_group)not_done_genetics"!=1,]
d0$functional_group<-factor(d0$functional_group,levels=unique(d0$functional_group))

#outcumel<-"addNA(genetic_group)gen_fynd"
outcumel<-"response_3mo_50"
#outcumes<-"genetic_group"
outcumes<-"response_3mo_50"

prev("genetic_group","seizure_type_fas_fias",d0)
logbin("addNA(genetic_group)gen_fynd",c("addNA(seizure_type_fas_fias)Yes"),d0)$table 

prev("genetic_group","seizure_type_fbtcs",d0)
logbin("addNA(genetic_group)gen_fynd",c("addNA(seizure_type_fbtcs)Yes"),d0)$table                  

prev("genetic_group","seizure_type_spasms",d0)
logbin("addNA(genetic_group)gen_fynd",c("addNA(seizure_type_spasms)Yes"),d0)$table              

prev("genetic_group","seizure_type_gtc",d0)
logbin("addNA(genetic_group)gen_fynd",c("addNA(seizure_type_gtc)Yes"),d0)$table            

prev("genetic_group","seizure_type_tonic",d0)
logbin("addNA(genetic_group)gen_fynd",c("addNA(seizure_type_tonic)Yes"),d0)$table              

prev("genetic_group","seizure_type_atonic",d0)
logbin("addNA(genetic_group)gen_fynd",c("addNA(seizure_type_atonic)Yes"),d0)$table            

prev("genetic_group","seizure_type_myoclonic",d0)
logbin("addNA(genetic_group)gen_fynd",c("addNA(seizure_type_myoclonic)Yes"),d0)$table 

logbin("addNA(genetic_group)gen_fynd",c("addNA(seizure_type_fas_fias)Yes","addNA(seizure_type_fbtcs)Yes","addNA(seizure_type_spasms)Yes","addNA(seizure_type_gtc)Yes","addNA(seizure_type_tonic)Yes","addNA(seizure_type_atonic)Yes","addNA(seizure_type_myoclonic)Yes"),d0)$table 

#table 6
prev("response_3mo_50","functional_group",d0)
logbin("response_3mo_50",c("functional_groupIon Channel", "functional_groupNeurotransmitter receptors", "functional_groupTransporters", "functional_groupSynapse-related", "functional_groupCell growth, division, and proliferation-related", "functional_groupCell structural integrity and or homeostasis", "functional_groupMetabolism-related", "functional_groupEnergy homeostasis and mitochondrial function", "functional_groupstrukturell", "functional_groupnot_done_genetics"),d0)$table
