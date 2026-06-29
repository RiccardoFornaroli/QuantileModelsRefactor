# FUNCTION TO COMPUTE MEAN WI
mean_wi<-function(models,taus){
AICC<-sapply(models,AICc<-function(x){aicc(x)})
dimnames(AICC)<-list(taus,names(models))
minAIC<-do.call(pmin, as.data.frame(AICC))
deltai<-apply(AICC,2,"-",minAIC)
denominatorWI<-apply(deltai,1,function(x){sum(exp(x*(-0.5)))})
numeratorWI<-apply(deltai,2,function(x){(exp(x*(-0.5)))})
wi<-apply(numeratorWI,2,"/",denominatorWI)
meanWI<-as.data.frame(rev(sort(apply(wi,2,mean))))
colnames(meanWI)<-"Wi"
return(meanWI)
}

#FUNCTION TO PLOT GLM E GLMM
plot_fit<-function(m,focal_var,inter_var=NULL,RE=NULL,n=20,n_core=4){
  require(arm)  
  dat<-model.frame(m)
  #turn all character variable to factor
  dat<-as.data.frame(lapply(dat,function(x){
    if(is.character(x)){
      as.factor(x)
    }
    else{x}
  }))
  #make a sequence from the focal variable
  x1<-list(seq(min(dat[,focal_var]),max(dat[,focal_var]),length=n))
  #grab the names and unique values of the interacting variables
  isInter<-which(names(dat)%in%inter_var)
  if(length(isInter)==1){
    x2<-list(unique(dat[,isInter]))
    names(x2)<-inter_var
  }
  if(length(isInter)>1){
    x2<-lapply(dat[,isInter],unique)
  }
  if(length(isInter)==0){
    x2<-NULL
  }
  #all_var<-x1
  #add the focal variable to this list
  all_var<-c(x1,x2)
  #expand.grid on it
  names(all_var)[1]<-focal_var
  all_var<-expand.grid(all_var)
  
  #remove varying variables and non-predictors
  dat_red<-dat[,-c(1,which(names(dat)%in%c(focal_var,inter_var,RE,"X.weights."))),drop=FALSE]
  if(dim(dat_red)[2]==0){
    new_dat<-all_var
  }  else{
    fixed<-lapply(dat_red,function(x) if(is.numeric(x)) mean(x) else factor(levels(x)[1],levels = levels(x)))
    #the number of rows in the new_dat frame
    fixed<-lapply(fixed,rep,dim(all_var)[1])
    #create the new_dat frame starting with the varying focal variable and potential interactions
    new_dat<-cbind(all_var,as.data.frame(fixed)) 
    #get the name of the variable to average over, debug for conditions where no variables are to be avergaed over
    name_f<-names(dat_red)[sapply(dat_red,function(x) ifelse(is.factor(x),TRUE,FALSE))]
  }  
    
  
  #get the predicted values
  cl<-class(m)[1]
  if(cl=="lm"){
    pred<-predict(m,newdata = new_dat,se.fit=TRUE)
  }
  
  if(cl=="glm" | cl=="negbin"){
    #predicted values on the link scale
    pred<-predict(m,newdata=new_dat,type="link",se.fit=TRUE)
  }
  # if(cl=="glmerMod" | cl=="lmerMod"){
    # pred<-list(fit=predict(m,newdata=new_dat,type="link",re.form=~0))
    # #for bootstrapped CI
    # new_dat<-cbind(new_dat,rep(0,dim(new_dat)[1]))
    # names(new_dat)[dim(new_dat)[2]]<-as.character(formula(m)[[2]])
    # mm<-model.matrix(formula(m,fixed.only=TRUE),new_dat)
  # }
  # #average over potential categorical variables  
  # if(length(name_f)>0){
    # if(cl=="glmerMod" | cl=="lmerMod"){
      # coef_f<-lapply(name_f,function(x) fixef(m)[grep(paste0("^",x,"\\w+$"),names(fixef(m)))])
    # }
    # else{
      # coef_f<-lapply(name_f,function(x) coef(m)[grep(paste0("^",x,"\\w+$"),names(coef(m)))])
    # }    
    # pred$fit<-pred$fit+sum(unlist(lapply(coef_f,function(x) mean(c(0,x)))))
  # }
  #to get the back-transform values get the inverse link function
  linkinv<-family(m)$linkinv
  
  #get the back transformed prediction together with the 95% CI for LM and GLM
  if(cl=="glm" | cl=="lm"){
    pred$pred<-linkinv(pred$fit)
    pred$LC<-linkinv(pred$fit-1.96*pred$se.fit)
    pred$UC<-linkinv(pred$fit+1.96*pred$se.fit)
  }
  
  # # for GLMM need to use bootstrapped CI, see ?predict.merMod
  # if(cl=="glmerMod" | cl=="lmerMod"){
    # pred$pred<-linkinv(pred$fit)
    # predFun<-function(.) mm%*%fixef(.)
    # bb<-bootMer(m,FUN=predFun,nsim=200,parallel="multicore",ncpus=n_core) #do this 200 times
    # bb$t<-apply(bb$t,1,function(x) linkinv(x))
    # # as we did this 200 times the 95% CI will be bordered by the 5th and 195th value
    # bb_se<-apply(bb$t,1,function(x) x[order(x)][c(5,195)])
    # pred$LC<-bb_se[1,]
    # pred$UC<-bb_se[2,] 
  # }
  
  #the output
  out<-as.data.frame(cbind(new_dat[,1:(length(inter_var)+1)],pred$LC,pred$pred,pred$UC))
  names(out)<-c(names(new_dat)[1:(length(inter_var)+1)],"LC","Pred","UC")
  return(out)
}

#FUNCTION TO ELABORATE MULTIPLE iButton DATA
iButtonMulti <- function(filenames, TYPE = "GIORNALIERO"){
ldf <- lapply(filenames, read.csv2)
res <- lapply(ldf, summary)
OUT<-list()
for(i in 1:length(ldf)){
int<-iButton(ldf[[i]], TYPE = TYPE)
OUT<-lappend(OUT,int)
}
OUT1<-OUT[[1]]
if (length(OUT)>1) {
for(i in 2:length(OUT)){
OUT1<-rbind(OUT1,OUT[[i]])
}
}
DATES<-as.POSIXct(OUT1[,1])
OUT1<-OUT1[order(DATES),]
OUT2<-data.frame(unique(OUT1$DATE))
for(i in 2:ncol(OUT1)){
int<-as.data.frame(tapply(OUT1[,i],OUT1$DATE, mean))
int$DATE<-rownames(int)
int<-int[order(int$DATE),]
int<-int[,1]
OUT2<-cbind(OUT2,int)
}
names(OUT2)<-names(OUT1)
RESULT<-OUT2
}

#FUNCTION TO ELABORATE iButton DATA
iButton <- function(dati, TYPE = "GIORNALIERO"){
#NUMERO DI MUSURE
n<-nrow(dati)
# INTERVALLO DI CALCOLO (SECONDI)
interval = 60
#CONVERSIONE DEL TESTO IN DATE
dati$Time<-strptime(dati$DateTime, format("%d/%m/%y %H.%M.%S"))
dati$day<-format(dati$Time, "%Y-%m-%d")
dati$day_hour<-format(dati$Time, "%Y-%m-%d %H")
# CREAZIONE DEGLI ISTANTI DI CUI SI VUOLE CALCOLARE LA TEMPERATURA
epoch<-strptime("1970-01-01 00:00:00", format("%Y-%m-%d %H:%M:%S"))
startTIME<-as.integer(as.POSIXct(dati$Time[1]))
stopTIME<-as.integer(as.POSIXct(dati$Time[n]))
INTERVAL_SIM<-seq(startTIME,stopTIME,interval)
#CREAZIONE DELLE SPLINE PER CALCOLO TEMPERATURA IN CONTINUO
spline_T<-smooth.spline(dati$Time, dati$Value , df=n-1)
# RESTITUZIONE DEI DATI ISTANTANEI CALCOLATI TRAMITE SPLINE
OUTspline<-predict(spline_T,INTERVAL_SIM)
OUTspline$day<-factor(format(as.POSIXct(OUTspline$x, origin = epoch),format("%Y-%m-%d")))
OUTspline$day_hour<-factor(format(as.POSIXct(OUTspline$x, origin = epoch),format("%Y-%m-%d %H")))
OUTspline$day_minute<-format(as.POSIXct(OUTspline$x, origin = epoch),format("%Y-%m-%d %H:%M"))
if (TYPE == "GIORNALIERO") {
# CREAZIONE DI UNA MATRICE CON LE STATISTICHE GIORNALIERE
DATE<-levels(OUTspline$day)
MEAN<-round(as.vector(tapply(OUTspline$y,OUTspline$day,mean)),2)
MEDIAN<-round(as.vector(tapply(OUTspline$y,OUTspline$day,median)),2)
MIN<-round(as.vector(tapply(OUTspline$y,OUTspline$day,min)),2)
MAX<-round(as.vector(tapply(OUTspline$y,OUTspline$day,max)),2)
PER10<-round(as.vector(tapply(OUTspline$y,OUTspline$day,quantile,probs = 0.1)),2)
PER90<-round(as.vector(tapply(OUTspline$y,OUTspline$day,quantile,probs = 0.9)),2)
STD<-round(as.vector(tapply(OUTspline$y,OUTspline$day,sd)),3)
M_MEAN<-round(as.vector(tapply(dati$Value,dati$day,mean)),2)
M_MIN<-round(as.vector(tapply(dati$Value,dati$day,min)),2)
M_MAX<-round(as.vector(tapply(dati$Value,dati$day,max)),2)
M_STD<-round(as.vector(tapply(dati$Value,dati$day,sd)),3)
RESULT<-data.frame(DATE,MEAN,MEDIAN,MIN,MAX,PER10,PER90,STD,M_MEAN,M_MIN,M_MAX,M_STD)
}
if (TYPE == "ORARIO") {
# CREAZIONE DI UNA MATRICE CON LE STATISTICHE ORARIE
DATE_HOUR<-levels(OUTspline$day_hour)
MEAN<-round(as.vector(tapply(OUTspline$y,OUTspline$day_hour,mean)),2)
MEDIAN<-round(as.vector(tapply(OUTspline$y,OUTspline$day_hour,median)),2)
MIN<-round(as.vector(tapply(OUTspline$y,OUTspline$day_hour,min)),2)
MAX<-round(as.vector(tapply(OUTspline$y,OUTspline$day_hour,max)),2)
PER10<-round(as.vector(tapply(OUTspline$y,OUTspline$day_hour,quantile,probs = 0.1)),2)
PER90<-round(as.vector(tapply(OUTspline$y,OUTspline$day_hour,quantile,probs = 0.9)),2)
STD<-round(as.vector(tapply(OUTspline$y,OUTspline$day_hour,sd)),3)
RESULT<-data.frame(DATE_HOUR,MEAN,MEDIAN,MIN,MAX,PER10,PER90,STD)
}
if (TYPE == "ISTANTANEO") {
# CREAZIONE DI UNA MATRICE CON I VALORI ISTANTANEI
DATE_MINUTE<-OUTspline$day_minute
TEMPERATURE<-round(OUTspline$y,2)
RESULT<-data.frame(DATE_MINUTE,TEMPERATURE)
}
RESULT
}

#FUNCTION TO APPEND ELEMENTS TO LIST
lappend <- function (lst, ...){
lst <- c(lst, list(...))
  return(lst)
}

#FUNCTION TO COMPUTE AICC
aicc<-function(x) AIC(x)+((2*length(coef(x))*(length(coef(x))+1))/(length(resid(x))-length(coef(x))-1))
lappend <- function (lst, ...){
lst <- c(lst, list(...))
  return(lst)
}

#FUNCTION TO COMPUTE VIF AND VARIABLES SELELCTION
vif_func<-function(in_frame,thresh=10,trace=T){
 
if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
#get initial vif value for all comparisons of variables
vif_init<-NULL
for(val in names(in_frame)){
form_in<-formula(paste(val,' ~ .'))
vif_init<-rbind(vif_init,c(val,VIF(lm(form_in,data=in_frame))))
}
vif_max<-max(as.numeric(vif_init[,2]))
 
if(vif_max < thresh){
if(trace==T){ #print output of each iteration
prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
cat('\n')
cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
}
return(names(in_frame))
}
else{
 
in_dat<-in_frame
 
#backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
while(vif_max >= thresh){
vif_vals<-NULL
 
for(val in names(in_dat)){
form_in<-formula(paste(val,' ~ .'))
vif_add<-VIF(lm(form_in,data=in_dat))
vif_vals<-rbind(vif_vals,c(val,vif_add))
}
max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2])))[1]
 
vif_max<-as.numeric(vif_vals[max_row,2])
 
if(vif_max<thresh) break
if(trace==T){ #print output of each iteration
prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
cat('\n')
cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
flush.console()
}
 
in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
 
}
 
return(names(in_dat))
}
} 

#FUNCTION TO REMOVE OUTLIERS
rm_outlier_15iqr<- function(data){
  data_new <- apply(data, 2, function(var) {
	upperq <- quantile(var)[4] 
	lowerq <- quantile(var)[2]
	iqr <- upperq - lowerq
	uppert <- (iqr * 1.5) + upperq
	lowert <- lowerq - (iqr * 1.5)
	sapply(var, function(y) {
if (y > uppert){
y <- NA
} else if (y < lowert){
y <- NA
} else {y <- y}
})
})
}

#FUNCTION TO REMOVE EXTREME OUTLIERS
rm_outlier_3iqr<- function(data){
  data_new <- apply(data, 2, function(var) {
	upperq <- quantile(var)[4] 
	lowerq <- quantile(var)[2]
	iqr <- upperq - lowerq
	uppert <- (iqr * 3) + upperq
	lowert <- lowerq - (iqr * 3)
	sapply(var, function(y) {
if (y > uppert){
y <- NA
} else if (y < lowert){
y <- NA
} else {y <- y}
})
})
}

#FUNCTION TO SUBSTITUTE EXTREME OUTLIERS
sub_outlier_3iqr<- function(data){
  data_new <- apply(data, 2, function(var) {
	upperq <- quantile(var)[4] 
	lowerq <- quantile(var)[2]
	iqr <- upperq - lowerq
	uppert <- (iqr * 3) + upperq
	lowert <- lowerq - (iqr * 3)
	sapply(var, function(y) {
if (y > uppert){
y <- uppert
} else if (y < lowert){
y <- lowert
} else {y <- y}
})
})
}


ASPT_ind <- function(data){
	BMWP_table <- unique(data[data$BMWP_SCORE > 0,c(4,5)])
	ASPT <- NA
	if (length(BMWP_table[,2]) > 0) {
	ASPT <- (round(sum(BMWP_table$BMWP_SCORE)/length(BMWP_table[,2]),3))
	} else {
	ASPT <- 2
	}	
	return(ASPT)
	}

Shannon_ind <- function(data){
	pi <- data[,2]/sum(data[,2])
	lnpi <- log(data[,2]/sum(data[,2]))
	return(round(-sum(pi*lnpi),3))
	}

LIFE_ind <- function(data){
score_1<-ifelse(data$ab_cat[data$LIFE==1]==1,9,ifelse(data$ab_cat[data$LIFE==1]==2,10,ifelse(data$ab_cat[data$LIFE==1]==3,11,12)))
score_2<-ifelse(data$ab_cat[data$LIFE==2]==1,8,ifelse(data$ab_cat[data$LIFE==2]==2,9,ifelse(data$ab_cat[data$LIFE==2]==3,10,11)))
score_3<-ifelse(data$ab_cat[data$LIFE==3]==1,7,ifelse(data$ab_cat[data$LIFE==3]==2,7,ifelse(data$ab_cat[data$LIFE==3]==3,7,7)))
score_4<-ifelse(data$ab_cat[data$LIFE==4]==1,6,ifelse(data$ab_cat[data$LIFE==4]==2,5,ifelse(data$ab_cat[data$LIFE==4]==3,4,3)))
score_5<-ifelse(data$ab_cat[data$LIFE==5]==1,5,ifelse(data$ab_cat[data$LIFE==5]==2,4,ifelse(data$ab_cat[data$LIFE==5]==3,3,2)))
score_6<-ifelse(data$ab_cat[data$LIFE==6]==1,4,ifelse(data$ab_cat[data$LIFE==6]==2,3,ifelse(data$ab_cat[data$LIFE==6]==3,2,1)))
n<-length(data$LIFE[data$LIFE>0])
LIFE<-(sum(score_1)+sum(score_2)+sum(score_3)+sum(score_4)+sum(score_5)+sum(score_6))/n
return(round(LIFE,3))
}

PSI_ind <- function(data){
score_A<-ifelse(data$ab_cat[data$FSSR=="A"]==1,2,ifelse(data$ab_cat[data$FSSR=="A"]==2,3,ifelse(data$ab_cat[data$FSSR=="A"]==3,4,5)))
score_B<-ifelse(data$ab_cat[data$FSSR=="B"]==1,1,ifelse(data$ab_cat[data$FSSR=="B"]==2,2,ifelse(data$ab_cat[data$FSSR=="B"]==3,3,4)))
score_C<-ifelse(data$ab_cat[data$FSSR=="C"]==1,1,ifelse(data$ab_cat[data$FSSR=="C"]==2,2,ifelse(data$ab_cat[data$FSSR=="C"]==3,3,4)))
score_D<-ifelse(data$ab_cat[data$FSSR=="D"]==1,2,ifelse(data$ab_cat[data$FSSR=="D"]==2,3,ifelse(data$ab_cat[data$FSSR=="D"]==3,4,5)))
PSI<-((sum(score_A)+sum(score_B))/(sum(score_A)+sum(score_B)+sum(score_C)+sum(score_D)))*100
return(round(PSI,3))
}

METRICS <- function(data){
	# load reference data
	base <- read.csv("TaxaList.csv", header=T, sep=";", dec=".");
	

	# format data
	data <- data[which(data[,2]>0),]
	data$ab_cat<-ifelse(data[,2]>10000,5,ifelse(data[,2]>1000,4,ifelse(data[,2]>100,3,ifelse(data[,2]>10,2,1))))
	data$GROUP <- base$GROUP[match(data[,1], base$TAXA)]
	data[,5:31] <- base[match(data[,1],base$TAXA),3:29]
	data<-droplevels(data)

	# metrics calculation
	Ind_Tot <- sum(data[,2])
	Ind_Tot_feed<-sum(data[!is.na(data$MIN),2])
	Ind_Tot_drift<-sum(data[!is.na(data$Accidental),2])
	Ind_Tot_FHG<-sum(data[!is.na(data$SWS),2])
	Ind_Tot_volt<-sum(data[!is.na(data$minoreuguale),2])
	N_BMWP_TAXA <- length(levels(data$BMWP_TAXA))
	ASPT <- ASPT_ind(data)
	log_SelEPTD <- round(log10(sum(data[data$selEPTD==1,2])+1),3)
	GOLD <- round(1-(sum(data[data$GOLD==1,2])/Ind_Tot),3)
	N_Fam <- length(data[,1])
	N_EPT_Fam <- length(data[data$EPT==1,1])
	Shannon <- Shannon_ind(data)
	LIFE<-LIFE_ind(data)
	PSI<-PSI_ind(data)
	Ind_EPT <- sum(data[data$EPT==1,2])
	LOG_Ind_Tot<-log10(sum(data[,2])+1)
	LOG_Ind_EPT <-log10(sum(data[data$EPT==1,2])+1)

	GRA<-sum(data[!is.na(data$MIN),2]*data[!is.na(data$MIN),12])/Ind_Tot_feed
	MIN<-sum(data[!is.na(data$MIN),2]*data[!is.na(data$MIN),13])/Ind_Tot_feed
	XYL<-sum(data[!is.na(data$MIN),2]*data[!is.na(data$MIN),14])/Ind_Tot_feed
	SHR<-sum(data[!is.na(data$MIN),2]*data[!is.na(data$MIN),15])/Ind_Tot_feed
	GAT<-sum(data[!is.na(data$MIN),2]*data[!is.na(data$MIN),16])/Ind_Tot_feed
	AFF<-sum(data[!is.na(data$MIN),2]*data[!is.na(data$MIN),17])/Ind_Tot_feed
	PFF<-sum(data[!is.na(data$MIN),2]*data[!is.na(data$MIN),18])/Ind_Tot_feed
	PRE<-sum(data[!is.na(data$MIN),2]*data[!is.na(data$MIN),19])/Ind_Tot_feed
	PAR<-sum(data[!is.na(data$MIN),2]*data[!is.na(data$MIN),20])/Ind_Tot_feed
	OTH<-sum(data[!is.na(data$MIN),2]*data[!is.na(data$MIN),21])/Ind_Tot_feed
	ACCIDENTAL<-sum(data[!is.na(data$Accidental),2]*data[!is.na(data$Accidental),22])/Ind_Tot_drift
	BEHAVIORAL<-sum(data[!is.na(data$Accidental),2]*data[!is.na(data$Accidental),23])/Ind_Tot_drift
	SWS<-sum(data[!is.na(data$SWS),2]*data[!is.na(data$SWS),24])/Ind_Tot_FHG
	SWD<-sum(data[!is.na(data$SWS),2]*data[!is.na(data$SWS),25])/Ind_Tot_FHG
	BUB<-sum(data[!is.na(data$SWS),2]*data[!is.na(data$SWS),26])/Ind_Tot_FHG
	SPW<-sum(data[!is.na(data$SWS),2]*data[!is.na(data$SWS),27])/Ind_Tot_FHG
	SES<-sum(data[!is.na(data$SWS),2]*data[!is.na(data$SWS),28])/Ind_Tot_FHG
	OTH1<-sum(data[!is.na(data$SWS),2]*data[!is.na(data$SWS),29])/Ind_Tot_FHG
	MINOREUGUALE<-sum(data[!is.na(data$minoreuguale),2]*data[!is.na(data$minoreuguale),30])/Ind_Tot_volt
	MAGGIORE<-sum(data[!is.na(data$minoreuguale),2]*data[!is.na(data$minoreuguale),31])/Ind_Tot_volt

	PR<-(MIN+GRA)/(SHR+GAT+PFF+AFF)
	CPOMFPOM<-(SHR+MIN)/(GAT+PFF+AFF)
	SPOMBPOM<-(AFF+PFF)/GAT
    HABITATSTABILITY<-(GRA+AFF+PFF)/(SHR+MIN+GAT)
	TOPDOWN<-(PRE+PAR)/(GRA+AFF+PFF+MIN+XYL+SHR+GAT+OTH)
	LIFECYCLE<-(MAGGIORE)/(MINOREUGUALE)
	BENTHICFOOD<-(SES+SPW)/(BUB+SWD+SWS)
	DRIFTFOOD<-(BEHAVIORAL)/(ACCIDENTAL)
	# produce output: table with all output variables
	out <- cbind(Ind_Tot,log_SelEPTD,GOLD,N_Fam,N_EPT_Fam,Shannon,ASPT,
	LIFE,PSI,Ind_EPT,LOG_Ind_Tot,LOG_Ind_EPT,GRA,MIN,XYL,SHR,GAT,AFF,
	PFF,PRE,PAR,OTH,ACCIDENTAL,BEHAVIORAL,SWS,SWD,BUB,SPW,SES,OTH1,
	MINOREUGUALE,MAGGIORE,PR,CPOMFPOM,SPOMBPOM,HABITATSTABILITY,TOPDOWN,
	LIFECYCLE,BENTHICFOOD,DRIFTFOOD)
	return(out)
}



