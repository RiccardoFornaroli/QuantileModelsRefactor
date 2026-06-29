
rm(list=ls())
### reading custom functions
source("Functions.R")
### reading packages (if not installed, automatically install)
package.list <- c("Hmisc","quantreg","fmsb","Amelia","sm")
tmp.install <- which(lapply(package.list, require, character.only = TRUE)==FALSE)
if(length(tmp.install)>0) install.packages(package.list[tmp.install])
lapply(package.list, require, character.only = TRUE)

### loading data
dati<-read.table("METRICS_2026_LB_OK.csv",header=T, sep = ",", na.string=c("NA"))
summary(dati)
names(dati)
dati[,c(6,7)]<-rm_outlier_15iqr(dati[,c(6,7)]) ## calcola i quantili di velocità e depth
dati<-dati[complete.cases(dati),]
dati$SUB<-factor(dati$SUB,levels=levels(factor(dati$SUB))[c(1,3,2)]) ##subtrato diventa categorico
summary(dati)

#######TRASFORMO LE METRICHE CHE NECESSITANO
##########################
names(dati)
##logit +1 delle abbondanze delle famiglie
Metrics_or<-dati[c(10:length(dati))]
fine_tassonomia<-which(names(Metrics_or)=="Viviparidae")
Metrics_or[,1:fine_tassonomia]<-log10(Metrics_or[, 1:fine_tassonomia] + 1)
Metrics<-Metrics_or
summary(Metrics)

##logit delle proporzioni
eps <- 0.001
Metrics$logit_EPT_prop <- log((Metrics$EPT_prop + eps) / ((1 - Metrics$EPT_prop) + eps))
Metrics$logit_OCH_prop <- log((Metrics$OCH_prop + eps) / ((1 - Metrics$OCH_prop) + eps))
Metrics$logit_EPT_EPTOCH <- log((Metrics$EPT_EPTOCH + eps) / ((1 - Metrics$EPT_EPTOCH) + eps))

Metrics$EPT_abu<-log10(Metrics$EPT_abu+1)
Metrics$OCH_abu<-log10(Metrics$OCH_abu+1)
Metrics$ABUNDANCE<-log10(Metrics$ABUNDANCE+1)
summary(Metrics)

#######SISTEMO LE VARIABILI IDRAULICHE
#####variabili idrauliche
Variables<-dati[c(6,7)]
Variables<-abs(Variables)                ##voglio solo valori positivi
Variables$VEL[Variables$VEL==0]<-0.001   #velocità zero diventano 0.001
Variables_ST<-Variables
Variables_ST<-as.data.frame(Variables_ST)
summary(Variables_ST)
GROUP<-factor(dati$GROUP)
SUB<-dati$SUB

# # #CHOICE OF QUANTILES
# taus_s<-seq(0.90,0.99,0.01)
taus_s<-seq(0.02,0.98,0.02)
taus_m<-list((2:10)/100,(45:55)/100,(90:98)/100)
names(taus_m)<-c("Floor","Median","Ceiling")
############################################################
#####################CHOICE OF MODELS#######################
############################################################

Models<-list(
null<-function(fit_data){rq(VAR ~ 1 , tau=taus, method="sfn", data= fit_data)},
SBS_GRP<-function(fit_data){rq(VAR ~ 1 + SUB + GROUP, tau=taus, method="sfn", data= fit_data)},
LIN_SBS_GRP<-function(fit_data){rq(VAR ~ 1 + INVARy + GROUP + SUB, tau=taus, method="sfn", data= fit_data)},
LOG_SBS_GRP<-function(fit_data){rq(VAR ~ 1 + log10(INVARy) + GROUP + SUB, tau=taus, method="sfn", data= fit_data)},
EXP_SBS_GRP<-function(fit_data){rq(VAR ~ 1 + exp(INVARy) + GROUP + SUB, tau=taus, method="sfn", data= fit_data)},
QUA_SBS_GRP<-function(fit_data){rq(VAR ~ 1 + poly(INVARy,2) + GROUP + SUB, tau=taus, method="sfn", data= fit_data)},
LIN_SBS_GRP_INT<-function(fit_data){rq(VAR ~ 1 + INVARy * GROUP + SUB, tau=taus, method="sfn", data= fit_data)},
LOG_SBS_GRP_INT<-function(fit_data){rq(VAR ~ 1 + log10(INVARy) * GROUP + SUB, tau=taus, method="sfn", data= fit_data)},
EXP_SBS_GRP_INT<-function(fit_data){rq(VAR ~ 1 + exp(INVARy) * GROUP + SUB, tau=taus, method="sfn", data= fit_data)},
QUA_SBS_GRP_INT<-function(fit_data){rq(VAR ~ 1 + poly(INVARy,2) * GROUP + SUB, tau=taus, method="sfn", data= fit_data)})

names(Models)<-c("null","SBS_GRP",
"LIN_SBS_GRP","LOG_SBS_GRP","EXP_SBS_GRP","QUA_SBS_GRP",
"LIN_SBS_GRP_INT","LOG_SBS_GRP_INT","EXP_SBS_GRP_INT","QUA_SBS_GRP_INT")
NULL_Models<-c(1:2)
NON_NULL_Models<-(1:length(Models))[(1:length(Models))%nin% NULL_Models]


######## UNIVARIATE MODELS FITTING FOR SHAPE SELECTION ##########
taus<-taus_s
wiSHAPE<-list()
# for (i in 12:13){
for (i in 1:length(Metrics)){
VAR<-Metrics[,i]
metric<-names(Metrics[i])
wiINVAR<-list()
for (y in 1:ncol(Variables_ST)){
print(Sys.time())
print(paste("Metric ",i,"/",length(Metrics)," Variable",y,"/",ncol(Variables_ST),sep=""))
perc<-((((i-1)*ncol(Variables_ST))+y)/(length(Metrics)*ncol(Variables_ST)))*100
print(paste("Progress ",round(perc,1),"%"))
INVARy<-Variables_ST[,y]
fit_data<-data.frame(VAR,INVARy,GROUP,SUB)
fit_data<-fit_data[complete.cases(fit_data),]
fitted_model<-list()
for (m in NON_NULL_Models){
fit<-Models[[m]](fit_data)
fitted_model<-lappend(fitted_model,fit)
}
names(fitted_model)<-names(Models)[NON_NULL_Models]
meanWI<-mean_wi(fitted_model,taus)
wiINVAR<-lappend(wiINVAR,meanWI)
}
names(wiINVAR)<-colnames(Variables_ST)
wiSHAPE<-lappend(wiSHAPE,wiINVAR)
}
names(wiSHAPE)<-names(Metrics)
save.image("wiSHAPE.RData")
load("wiSHAPE.RData" ) 


# Variables selection: the variable will always be included in the subsequent calculation, 
# this is just to select the shape
wi_selection<-list()
wi_selected<-list()
for (v in 1:ncol(Metrics)){
selection<-as.data.frame(matrix(ncol=5,nrow=ncol(Variables_ST),NA))
names(selection)<-c("Index","Variable","Selected","Model","Wi")
selection$Variable<-colnames(Variables_ST)
selection$Index<-1:ncol(Variables_ST)
for (i in 1:ncol(Variables_ST)){
    # null_models<-which(rownames(wiSHAPE[[v]][[i]])%in%names(Models)[NULL_Models])
# if ((null_models[1]>1)&(all(wiSHAPE[[v]][[i]][null_models,1]<0.05))) {
# if ((null_models[1]>1)) {
   selection$Selected[i]<-"YES"
# } else {
   # selection$Selected[i]<-"NO"
# }
selection$Model[i]<-rownames(wiSHAPE[[v]][[i]])[1]
selection$Wi[i]<-wiSHAPE[[v]][[i]][1,1]
index_selected<-which(selection$Selected=="YES")
selected<-selection[selection$Selected=="YES",]
}
wi_selection<-lappend(wi_selection,selection)
wi_selected<-lappend(wi_selected,selected)
}
names(wi_selection)<-names(Metrics)
names(wi_selected)<-names(Metrics)


var_list<-list()
for (v in 1:ncol(Metrics)){
if (nrow(wi_selected[[v]])>1){
Variable<-try(vif_func(in_frame=Variables_ST[,wi_selected[[v]]$Index],thresh=1.5,trace=F))
if (Variable[1]=="Error in if (vif_max < thresh) break : argument is of length zero\n"){
Variable<-wi_selected[[v]]$Variable[which.max(wi_selected[[v]]$Wi)]
}
Model<-wi_selected[[v]]$Model[wi_selected[[v]]$Variable %in% Variable]
Shape<-substr(Model, start = 1, stop = 4)
Shape<-gsub('[M]', '', Shape)
Shape<-gsub('[_]', '', Shape)
sel<-cbind(Variable,Model,Shape)
var_list<-lappend(var_list,sel)
} else {
if (nrow(wi_selected[[v]])>0){
Variable<-colnames(Variables_ST)[wi_selected[[v]]$Index]
Model<-wi_selected[[v]]$Model[wi_selected[[v]]$Variable %in% Variable]
Shape<-substr(Model, start = 1, stop = 4)
Shape<-gsub('[M]', '', Shape)
Shape<-gsub('[_]', '', Shape)
sel<-cbind(Variable,Model,Shape)
var_list<-lappend(var_list,sel)
} else {
var_list<-lappend(var_list,NA)
}}}
names(var_list)<-names(Metrics)
var_list
save.image("wiSHAPE.RData")
load("wiSHAPE.RData" ) 

#### Multivariate models fitting and evaluation #################
var_final<-list()
for (v in 1:length(var_list)){
# for (v in 2:2){
taus<-taus_s
print(Sys.time())
print(paste("Metric ",v,"/",length(Metrics)))
perc<-(v/length(Metrics))*100
print(paste("Progress ",round(perc,1),"%"))
var_met_final<-list()
var_met_final<-lappend(var_met_final,names(var_list)[v])
if (all(!is.na(var_list[[v]]))){
VAR<-Metrics[[v]]
invar<-c()
for (i in 1:nrow(var_list[[v]])){
if (var_list[[v]][i,3]=="LIN"){
inv<-as.character(var_list[[v]][i,1])
invar<-append(invar,inv)
} else {
if (var_list[[v]][i,3]=="LOG"){
inv<-paste("log10(",var_list[[v]][i,1],")",sep="")
invar<-append(invar,inv)
} else {
if (var_list[[v]][i,3]=="EXP"){
inv<-paste("exp(",var_list[[v]][i,1],")",sep="")
invar<-append(invar,inv)
} else {
if (var_list[[v]][i,3]=="QUA"){
inv<-paste("poly(",var_list[[v]][i,1],",2)",sep="")
invar<-append(invar,inv)
}}}}}
nvar<-length(invar)
var_met_final<-lappend(var_met_final,wi_selected[[v]])
var_met_final<-lappend(var_met_final,invar)
fit_data<-as.data.frame(cbind(VAR,Variables_ST))
invar_int<-paste("(",invar,"):GROUP",sep="")
invar_int_iniz<-paste("(",invar,"):GROUP",sep="")

# backward stepwise selection
sel<-list()
full_formula<-as.formula(paste("VAR ~ 1 + ",paste(invar,collapse=" + ")," + ",paste(invar_int,collapse=" + "), " + GROUP + SUB"))
full<-rq(full_formula, tau=taus, method="sfn", data= fit_data)
mod<-list(full)
names(mod)<-c("start")
meanWI<-mean_wi(mod,taus)

# selection interactions
while (rownames(meanWI)[1]!="full" & length(invar_int)>1){
for (b in 1:length(invar_int)){
res_formula<-as.formula(paste("VAR ~ 1 + ",paste(invar,collapse=" + ")," + ",paste(invar_int[-b],collapse=" + "), " + GROUP + SUB"))
res<-rq(res_formula, tau=taus, method="sfn", data= fit_data)
mod<-lappend(mod,res)
}
names(mod)<-c("full",invar_int)
meanWI<-mean_wi(mod,taus)
invar_int<-invar_int[invar_int!=rownames(meanWI)[1]]
full_formula<-as.formula(paste("VAR ~ 1 + ",paste(invar,collapse=" + ")," + ",paste(invar_int,collapse=" + "), " + GROUP + SUB"))
full<-rq(full_formula, tau=taus, method="sfn", data= fit_data)
mod<-list(full)
sel<-lappend(sel,meanWI)
}
if (length(invar_int)==1) {
res_formula<-as.formula(paste("VAR ~ 1 + ",paste(invar,collapse=" + "), " + GROUP + SUB"))
last_res<-rq(res_formula, tau=taus, method="sfn", data= fit_data)
mod<-list(full,last_res)
names(mod)<-c("full","last_res")
meanWI<-mean_wi(mod,taus)
ifelse(rownames(meanWI)[1]!="full",invar_int<-invar_int[-1],invar_int<-invar_int)
sel<-lappend(sel,meanWI)
}
ifelse(length(invar_int)==length(invar_int_iniz),all_int<-T,all_int<-F)
invar_to_mantain<-which(invar_int_iniz %in% invar_int)
print("Interaction selected")

# selection variable
# se viente tenuta l'interazione ovviamnete viente tenuta la variabile
# salta il pezzo seguente quando il ciclo sopra tiente tutte le interazioni
sel_i<-list()
full_formula<-as.formula(paste("VAR ~ 1 + ",paste(invar,collapse=" + ")," + ",paste(invar_int,collapse=" + "), " + GROUP + SUB"))
full<-rq(full_formula, tau=taus, method="sfn", data= fit_data)
mod<-list(full)
names(mod)<-c("start")
meanWI<-mean_wi(mod,taus)
while (rownames(meanWI)[1]!="full" & (length(invar)-length(invar_to_mantain))>1 & all_int==F){
for (b in 1:length(invar)){
res_formula<-as.formula(paste("VAR ~ 1 + ",paste(invar[-b],collapse=" + ")," + ",paste(invar_int,collapse=" + "), " + GROUP + SUB"))
res<-rq(res_formula, tau=taus, method="sfn", data= fit_data)
mod<-lappend(mod,res)
}
names(mod)<-c("full",invar)
meanWI<-mean_wi(mod,taus)
invar<-invar[invar!=rownames(meanWI)[1]]
full_formula<-as.formula(paste("VAR ~ 1 + ",paste(invar,collapse=" + ")," + ",paste(invar_int,collapse=" + "), " + GROUP + SUB"))
full<-rq(full_formula, tau=taus, method="sfn", data= fit_data)
mod<-list(full)
sel_i<-lappend(sel_i,meanWI)
}
if ((length(invar)-length(invar_to_mantain))==1) {
res_formula<-as.formula(paste("VAR ~ 1 + ",paste(invar[invar_to_mantain],collapse=" + ")," + ",paste(invar_int,collapse=" + "), "  + GROUP + SUB"))
last_res<-rq(res_formula, tau=taus, method="sfn", data= fit_data)
mod<-list(full,last_res)
names(mod)<-c("full","last_res")
meanWI<-mean_wi(mod,taus)
ifelse(rownames(meanWI)[1]!="full",invar<-invar[invar_to_mantain],invar<-invar)
sel_i<-lappend(sel_i,meanWI)
}
print("Variable selected")
# evaluation of grouping factor for intercept
full_formula<-as.formula(paste("VAR ~ 1 + ",paste(invar,collapse=" + ")," + ",paste(invar_int,collapse=" + "), " + GROUP + SUB"))
full<-rq(full_formula, tau=taus, method="sfn", data= fit_data)
res_grp_formula<-as.formula(paste("VAR ~ 1 + ",paste(invar,collapse=" + ")," + ",paste(invar_int,collapse=" + "), " + SUB"))
res_group<-rq(res_grp_formula, tau=taus, method="sfn", data= fit_data)
res_sub_formula<-as.formula(paste("VAR ~ 1 + ",paste(invar,collapse=" + ")," + ",paste(invar_int,collapse=" + "), " + GROUP"))
res_sub<-rq(res_sub_formula, tau=taus, method="sfn", data= fit_data)
mod<-list(full,res_group,res_sub)
names(mod)<-c("full","res_group","res_sub")
meanWI<-mean_wi(mod,taus)
if(rownames(meanWI)[1]=="full") {
full<-full
sig<-"full"
} else if (rownames(meanWI)[1]=="res_group"){
full_formula<-as.formula(paste("VAR ~ 1 + ",paste(invar,collapse=" + ")," + ",paste(invar_int,collapse=" + "), " + SUB"))
full<-rq(full_formula, tau=taus, method="sfn", data= fit_data)
res_sub_formula<-as.formula(paste("VAR ~ ",paste(invar,collapse=" + ")," + ",paste(invar_int,collapse=" + "), " + 1"))
res_sub<-rq(res_sub_formula, tau=taus, method="sfn", data= fit_data)
mod<-list(full,res_sub)
names(mod)<-c("full","res_sub")
meanWI<-mean_wi(mod,taus)
if(rownames(meanWI)[1]=="full") {
full<-full
sig<-"no_group"
} else {
full<-res_sub
full_formula<-res_sub_formula
sig<-"no"
}
} else if (rownames(meanWI)[1]=="res_sub"){
full_formula<-as.formula(paste("VAR ~ 1 + ",paste(invar,collapse=" + ")," + ",paste(invar_int,collapse=" + "), " + GROUP"))
full<-rq(full_formula, tau=taus, method="sfn", data= fit_data)
res_grp_formula<-as.formula(paste("VAR ~ ",paste(invar,collapse=" + ")," + ",paste(invar_int,collapse=" + "), " + 1"))
res_group<-rq(res_grp_formula, tau=taus, method="sfn", data= fit_data)
mod<-list(full,res_group)
names(mod)<-c("full","res_group")
meanWI<-mean_wi(mod,taus)
if(rownames(meanWI)[1]=="full"){
full<-full
sig<-"no_sub"
} else {
full<-res_group
full_formula<-res_grp_formula
sig<-"no"
}}
print("Grouping selected")
}
var_met_final<-lappend(var_met_final,invar)
var_met_final<-lappend(var_met_final,invar_int)
var_met_final<-lappend(var_met_final,sig)
var_met_final<-lappend(var_met_final,full_formula)
# probabilities

if (sig=="full"){
res_formula<-as.formula("VAR ~ 1  + GROUP + SUB")
} else if (sig=="no_group"){
res_formula<-as.formula("VAR ~ 1  + SUB")
} else if (sig=="no_sub"){
res_formula<-as.formula("VAR ~ 1  + GROUP")
} else if (sig=="no"){
res_formula<-as.formula("VAR ~ 1 ")
}
res<-rq(res_formula, tau=taus, method="sfn", data= fit_data)
mod<-list(full,res)
names(mod)<-c("full","res")
meanWI<-mean_wi(mod,taus)
var_met_final<-lappend(var_met_final,meanWI)
var_met_final<-lappend(var_met_final,full)
# evaluation of multiple quantiles groups
for (i in 1:length(taus_m)){
taus_i<-taus_m[[i]]
full<-rq(full_formula, tau=taus_i, method="sfn", data= fit_data)
res<-rq(res_formula, tau=taus_i, method="sfn", data= fit_data)
mod<-list(full,res)
names(mod)<-c("full","res")
meanWI<-mean_wi(mod,taus_i)
var_met_final<-lappend(var_met_final,meanWI)
}
names(var_met_final)<-c("Metric","Selected Variables","Invar_Inizial","Invar_Selected","Invar_int_Selected","Intercept_grouping","Formula","Significance","Model","Wi_floor","Wi_median","Wi_ceiling")
var_final<-lappend(var_final,var_met_final)
}
names(var_final)<-names(Metrics)



save.image("var_final_original.RData")
load("var_final_original.RData") 

var_final



################################################################################
#### Pairing variables and models shape for sebsequent analyses#################
################################################################################

# results<-as.data.frame(matrix(ncol=9,nrow=length(Metrics),NA))
# names(results)<-c("Metrics","Variables","Interaction","Intercept_grouping","Wi","Floor","Median","Ceiling","Formula")
# for (i in 1:length(Metrics)){
# results$Metrics[i]<-names(Metrics)[i]
# results$Variables[i]<-paste(var_final[[i]]$Invar_Selected, collapse=" ")
# results$Interaction[i]<-paste(var_final[[i]]$Invar_int_Selected, collapse=" ")
# results$Intercept_grouping[i]<-var_final[[i]]$Intercept_grouping
# results$Wi[i]<-var_final[[i]]$Significance["full",]
# results$Floor[i]<-ifelse(var_final[[i]]$Wi_floor["full",]>0.999,"YES","NO")
# results$Median[i]<-ifelse(var_final[[i]]$Wi_median["full",]>0.999,"YES","NO")
# results$Ceiling[i]<-ifelse(var_final[[i]]$Wi_ceiling["full",]>0.999,"YES","NO")
# results$Formula[i]<-toString(var_final[[i]]$Formula)
# }
# head(results)




############################################################################
# Prepariamo la nuova matrice a 12 colonne
results <- as.data.frame(matrix(ncol=12, nrow=length(Metrics), NA))
names(results) <- c("Metrics", "Variables", "Interaction", "Intercept_grouping", 
                    "Wifull", "Floor_Sig", "Wi_Floor", "Median_Sig", "Wi_Median", "Ceiling_Sig", "Wi_Ceiling", "Formula")

for (i in 1:length(Metrics)){
  results$Metrics[i] <- names(Metrics)[i]
  
  # Variabili e interazioni
  results$Variables[i]   <- paste(var_final[[i]]$Invar_Selected, collapse=" ")
  results$Interaction[i] <- paste(var_final[[i]]$Invar_int_Selected, collapse=" ")
  results$Intercept_grouping[i] <- var_final[[i]]$Intercept_grouping
  
  # 1. Il Wifull generale del modello (quello che mi chiedevi ora)
  results$Wifull[i] <- var_final[[i]]$Significance["full", ][1]
  
  # 2. Quantile basso (Floor - 0.05)
  results$Wi_Floor[i]  <- var_final[[i]]$Wi_floor["full", ][1]
  results$Floor_Sig[i] <- ifelse(results$Wi_Floor[i] > 0.999, "YES", "NO")
  
  # 3. Mediana (Median - 0.50)
  results$Wi_Median[i]  <- var_final[[i]]$Wi_median["full", ][1]
  results$Median_Sig[i] <- ifelse(results$Wi_Median[i] > 0.999, "YES", "NO")
  
  # 4. Quantile alto (Ceiling - 0.95)
  results$Wi_Ceiling[i]  <- var_final[[i]]$Wi_ceiling["full", ][1]
  results$Ceiling_Sig[i] <- ifelse(results$Wi_Ceiling[i] > 0.999, "YES", "NO")
  
  # Formula pulita in testo
  raw_formula <- var_final[[i]]$Formula
  results$Formula[i] <- paste(deparse(raw_formula), collapse = "")
}


head(results)
tail(results,15)
##write.table(results,"results_original_LB.csv", sep = ",", dec = ".", row.names = F, col.names = TRUE)
results


# Crea una lista vuota che conterrà i modelli pronti per i grafici
modelli_pronti <- list()
# Ciclo su tutte le metriche per salvare i modelli reali
for (i in 1:length(Metrics)) {
  m <- names(Metrics)[i]
  # Controlliamo se la metrica è significativa in almeno un quantile
  if (results$Floor_Sig[i] == "YES" || results$Median_Sig[i] == "YES" || results$Ceiling_Sig[i] == "YES") {
      # Prepariamo i dati per questa specifica metrica
    fit_data$VAR <- Metrics[, m]
      # Convertiamo la stringa di testo della formula in una formula di R
    vera_formula <- as.formula(results$Formula[i])
    # Fittiamo il modello rq sui 3 quantili che servono per il grafico
    # Usiamo tryCatch per evitare blocchi in caso di errori matematici
    modelli_pronti[[m]] <- tryCatch({
      rq(vera_formula, tau = c(0.05, 0.50, 0.95), method = "sfn", data = fit_data)
    }, error = function(e) { NULL })
  }
}

# SALVATAGGIO DEFINITIVO DI TUTTO L'AMBIENTE
save(results, modelli_pronti, file = "Tutto_Pronto_Per_Grafici.RData")

