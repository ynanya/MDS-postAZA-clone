
library(clinfun)
library(patchwork)
library(MASS)


##  make a comprehensive table for pre- and post- samples 
both_evaluate_sample<-subset(clini_dat, !is.na(before_filename) & !is.na(after_filename))


mut_dat %>% filter(Sample %in% c(both_evaluate_sample$before_filename, both_evaluate_sample$after_filename) &  ID!=0 )-> mut_dat.tmp
mut_dat.tmp %>% filter(Gene.refGene != "DDX41" | (Gene.refGene=="DDX41" & DDX41=="s")) -> mut_dat.tmp
mut_dat.tmp$timing <- ifelse(mut_dat.tmp$Sample %in% both_evaluate_sample$before_filename, "pre", "post")
for(m in 1:nrow(mut_dat.tmp)){
  mut_dat.tmp$patient[m] <- subset(clini_dat, before_filename == mut_dat.tmp$Sample[m] | after_filename == mut_dat.tmp$Sample[m])$patientID[1]
  mut_dat.tmp$aVAF[m] <- ifelse(mut_dat.tmp$aVAF[m]>1,1,mut_dat.tmp$aVAF[m])
}
mut_dat.tmp$variantID <- paste(mut_dat.tmp$patient, mut_dat.tmp$ID, sep="_")
data.frame(variantID=mut_dat.tmp$variantID, aVAF=mut_dat.tmp$aVAF, timing=mut_dat.tmp$timing, gene=mut_dat.tmp$Gene.refGene, patient=mut_dat.tmp$patient)  %>% 
  pivot_wider(names_from = timing, values_from = aVAF) -> mut_dat.wide


## if needed, add clinical information here
data.frame(patient=both_evaluate_sample$patientID, response=both_evaluate_sample$response_best3, 
           TP53=both_evaluate_sample$TP53, blast=both_evaluate_sample$blast) -> df_response


mut_dat.wide %>% left_join(df_response, by="patient") -> mut_dat.wide
mut_dat.wide$pre <- ifelse(is.na(mut_dat.wide$pre),0, mut_dat.wide$pre)
mut_dat.wide$post <- ifelse(is.na(mut_dat.wide$post),0, mut_dat.wide$post)
mut_dat.wide$delta <- with(mut_dat.wide, post-pre)


## if needed, add clinical information here
data.frame(patient=clini_dat$patientID, OS=clini_dat$days_OS, event=clini_dat$days_censor, IPSSR_score=clini_dat$total, 
           days_tx=clini_dat$days_transplant, event_tx=clini_dat$transplant_censor, TP53_multi=clini_dat$TP53_multi, DDX41=(clini_dat$DDX41>0))->tmp.df
mut_dat.wide %>% left_join(tmp.df, by="patient") -> mut_dat.wide


## set flag_pre_max 
mut_dat.wide$flag_pre_max<-0
for(m in 1:nrow(mut_dat.wide)){
  patient_<-mut_dat.wide$patient[m]
  value_<-mut_dat.wide$pre[m]
  
  subset(mut_dat.wide, patient==patient_)$pre ->tmp.values
  if(mut_dat.wide$pre[m] == max(tmp.values) & max(tmp.values)>0){
    if(patient_ == "MDS104P26" & mut_dat.wide$gene[m]=="ASXL1"){mut_dat.wide$flag_pre_max[m]<-0}
    else if(patient_ == "MDS231N49" & mut_dat.wide$gene[m]=="ASXL1"){mut_dat.wide$flag_pre_max[m]<-0}
    else if(patient_ == "MDS274S35" & mut_dat.wide$gene[m]=="TET2"){mut_dat.wide$flag_pre_max[m]<-0}
    else if(patient_ == "MDS486R44" & mut_dat.wide$gene[m]=="ASXL1"){mut_dat.wide$flag_pre_max[m]<-0}
    else if(patient_ == "MDS795B46" & mut_dat.wide$gene[m]=="NRAS"){mut_dat.wide$flag_pre_max[m]<-0}
    else{mut_dat.wide$flag_pre_max[m]<-1}
  }
}

subset(mut_dat.wide,flag_pre_max==1) -> mut_dat.wide.pre_max
mut_dat.wide.pre_max %>% group_by(gene) %>% summarise(n()) -> gene_count 
gene_count %>% arrange(desc(`n()`)) %>% filter(`n()`>5) -> gene_count





## set flag_post_max 
mut_dat.wide$flag_post_max<-0
for(m in 1:nrow(mut_dat.wide)){
  patient_<-mut_dat.wide$patient[m]
  value_<-mut_dat.wide$post[m]
  
  subset(mut_dat.wide, patient==patient_)$post ->tmp.postvalues
  subset(mut_dat.wide, patient==patient_)$pre ->tmp.prevalues
  subset(mut_dat.wide, patient==patient_) -> submut_dat
  submut_dat$post -> tmp.postvalues
  submut_dat$pre -> tmp.prevalues
  
  
  if(mut_dat.wide$post[m] == max(tmp.postvalues)){
    if(patient_ == "MDS231N49" & mut_dat.wide$gene[m]=="ASXL1"){mut_dat.wide$flag_post_max[m]<-0
    }else if(sum(tmp.postvalues==max(tmp.postvalues))>1){
      if(mut_dat.wide$pre[m]==max(tmp.prevalues)){mut_dat.wide$flag_post_max[m]<-1}
    }else if(sum(tmp.postvalues==max(tmp.postvalues))==1){mut_dat.wide$flag_post_max[m]<-1}
  }
}

mut_dat.wide %>% filter(flag_post_max==1) %>% distinct(patient, .keep_all=TRUE) -> mut_dat.wide.post_max

mut_dat.wide.post_max %>% group_by(gene) %>% summarise(n()) -> gene_count 
gene_count %>% arrange(desc(`n()`)) %>% filter(`n()`>5) -> gene_count




## make response2 (mCR/PR combined)
mut_dat.wide$response2 <- ifelse(mut_dat.wide$response %in% c("mCR","PR"), "mCR/PR", as.character(mut_dat.wide$response))
mut_dat.wide$response2 <- factor(mut_dat.wide$response2, levels= c("CR", "mCR/PR", "SD+HI", "SD-HI", "PD"), ordered = T)

mut_dat.wide$response3 <- ifelse(mut_dat.wide$response %in% c("CR","mCR","PR","SD+HI"), "OR+", ifelse(mut_dat.wide$response %in% c("PD","SD-HI"),"OR-","NA"))
mut_dat.wide$response3 <- factor(mut_dat.wide$response3, levels=c("OR+","OR-"))

mut_dat.wide$response4 <- ifelse(mut_dat.wide$response %in% c("CR","mCR","PR"), "OR+", ifelse(mut_dat.wide$response %in% c("PD","SD+HI", "SD-HI"),"OR-","NA"))
mut_dat.wide$response4 <- factor(mut_dat.wide$response4, levels=c("OR+","OR-"))

## IPSS-R stratification
mut_dat.wide$IPSSR <- ifelse(mut_dat.wide$IPSSR_score<=4.5,"Very Low/Low/Int", ifelse(mut_dat.wide$IPSSR_score<=6,"High",ifelse(mut_dat.wide$IPSSR_score>6,"Very High","NA")))
mut_dat.wide$IPSSR <- factor(mut_dat.wide$IPSSR, levels=c("Very Low/Low/Int", "High","Very High"), ordered = T)


mut_dat.wide %>% pivot_longer(cols=c(pre, post), names_to = "timing") -> mut_dat.long
mut_dat.long$timing <-factor(mut_dat.long$timing, levels=c("pre","post"), ordered = T)












##  IPSS-R
res <-coxph(Surv(OS/30, event)~IPSSR_score, data=subset(mut_dat.wide, flag_post_max==1))
summary(res)

res <-survfit(Surv(OS/30, event)~IPSSR, data=subset(mut_dat.wide, flag_post_max==1))
ggsurvplot(res,  surv.median.line = "h",  pval = TRUE,
           xlab="Time from AZA treatment (months)",title="IPSS-R, post-treatment+ ", ylab="Overall survival",
           risk.table=T, risk.table.height=0.2 ,palette ="npg",
           font.title=c(15, "bold"), legend.labs = c("Very Low/Low/Int", "High","Very High"),
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())


##  IPSS-R + genetic
res <-coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41, data=subset(mut_dat.wide, flag_post_max==1))
summary(res)


library(compareC)
with(subset(mut_dat.wide, flag_post_max==1),
  compareC(OS, event, (-1)*IPSSR_score, (-1)*(IPSSR_score*(0.21502)+TP53_multi*(0.51373)+DDX41*(-1.03281)))
)
#Improve by incorporating genetic data


##  IPSS-R + response3 (CR/mCR/PR/SD+HI vs SD-HI,PD)

res <-survfit(Surv(OS/30, event)~ response3, data=subset(mut_dat.wide, flag_post_max==1))
ggsurvplot(res,  surv.median.line = "h",  pval = TRUE,
           xlab="Time from AZA treatment (months)",title="IPSS-R, post-treatment+ ", ylab="Overall survival",
           risk.table=T, risk.table.height=0.2 ,palette ="npg",
           font.title=c(15, "bold"), 
           #legend.labs = c("Very Low/Low/Int", "High","Very High"),
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())


res <-coxph(Surv(OS/30, event)~IPSSR_score + response3, data=subset(mut_dat.wide, flag_post_max==1))
summary(res)

############################################################################################################
##  IPSS-R + genetic + postVAF
res <-coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41 + post + response4, data=subset(mut_dat.wide, flag_post_max==1))
summary(res)
forest_model(res)
############################################################################################################

############################################################################################################ 2021/4/3
mut_dat.wide$score<-with(mut_dat.wide, IPSSR_score*1 + TP53_multi*2 - 3*DDX41 +3*post + (response4=="OR-"))
subset(mut_dat.wide, flag_post_max==1)-> mut_dat.wide_nonredundant
hist(mut_dat.wide_nonredundant$score, breaks=seq(0,16,0.5))

mut_dat.wide_nonredundant$score_stratify <- ifelse(mut_dat.wide_nonredundant$score<7, "Low", ifelse(mut_dat.wide_nonredundant$score < 11, "Int", "High"))

res<-survfit(Surv(OS/30, event)~score_stratify, data=mut_dat.wide_nonredundant)
ggsurvplot(res,  surv.median.line = "none",  pval = TRUE,
           xlab="Time from AZA treatment (months)",title="OS by full model", ylab="Overall survival", palette="npg",
           risk.table=T, risk.table.height=0.2,
           font.title=c(15, "bold"), legend.labs = c("High", "Int", "Low"),
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())
surv_median(res)
############################################################################################################ 



############################################################################################################
##  IPSS-R + genetic + postVAF, DDX41+ cases
res <-coxph(Surv(OS/30, event)~IPSSR_score + post + response4, data=subset(mut_dat.wide, DDX41>0 & flag_post_max==1))
summary(res)
forest_model(res)
############################################################################################################




######################################################################################################################
######################################################################################################################
######################################################################################################################
### transplanted cases
mut_dat.wide %>% filter(flag_post_max==1 & event_tx ==1) -> tmp.df

res.cut <- surv_cutpoint(tmp.df, time = "OS", event = "event", variables = c("post"))
# optimal cutpoint in surv_cutpoint function of maxstat package.  
summary(res.cut)
post_median=res.cut$cutpoint[1,1]
vaf_threshold=post_median

res<-survfit(Surv(OS/30, event)~(post>vaf_threshold),data=tmp.df)
ggsurvplot(res,  surv.median.line = "hv",  pval = TRUE, conf.int=F, palette="jco",
           xlab="Time from AZA treatment (months)",title="OS by transplantation, post AZA\n time from aza ", ylab="Overall survival",
           risk.table=T, risk.table.height=0.2,
           font.title=c(15, "bold"), 
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())
surv_median(res)
res<-coxph(Surv(OS/30, event)~(post>vaf_threshold),data=tmp.df)
summary(res)


res<-survfit(Surv(OS/30, event)~(post>vaf_threshold),data=tmp.df)
ggsurvplot(res,  surv.median.line = "hv",  pval = TRUE, conf.int=F, palette="npg",
           xlab="Time from transplant (months)",title="OS after azacitidine", ylab="Overall survival",
           risk.table=T, risk.table.height=0.2,
           font.title=c(15, "bold"), 
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())
surv_median(res)



tmp.df %>% mutate(OS_tx = OS - days_tx) -> tmp.df
res.cut <- surv_cutpoint(tmp.df, time = "OS_tx"  , event = "event", variables = c("post"))
# optimal cutpoint in surv_cutpoint function of maxstat package.  
summary(res.cut)
post_median=res.cut$cutpoint[1,1]
vaf_threshold=post_median

res<-survfit(Surv(OS_tx/30, event)~(post>vaf_threshold),data=tmp.df)
ggsurvplot(res,  surv.median.line = "hv",  pval = TRUE, conf.int=F, palette="npg",
           xlab="Time from transplant (months)",title="OS after transplantation", ylab="Overall survival",
           risk.table=T, risk.table.height=0.2,
           font.title=c(15, "bold"), 
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())
surv_median(res)
summary(coxph(Surv(OS_tx/30, event)~(post>vaf_threshold),data=tmp.df))


res<-survfit(Surv(OS/30, event)~(post>vaf_threshold),data=tmp.df)
ggsurvplot(res,  surv.median.line = "hv",  pval = TRUE, conf.int=F, palette="npg",
           xlab="Time from transplant (months)",title="OS after azacitidine", ylab="Overall survival",
           risk.table=T, risk.table.height=0.2,
           font.title=c(15, "bold"), 
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())
surv_median(res)



######################################################################################################################
######################################################################################################################
######################################################################################################################







######################################################################################################################################################
##  IG vs I
anova(coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41, data=dataset), coxph(Surv(OS/30, event)~IPSSR_score , data=dataset))

##  IR vs I
anova(coxph(Surv(OS/30, event)~IPSSR_score + response4, data=dataset), coxph(Surv(OS/30, event)~IPSSR_score , data=dataset))

##  IGR vs I
anova(coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41 + response4, data=dataset), coxph(Surv(OS/30, event)~IPSSR_score , data=dataset))

##  IGRP vs I
anova(coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41 + response4 + post, data=dataset), coxph(Surv(OS/30, event)~IPSSR_score , data=dataset))

## IGR vs IG
anova(coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41 + response4, data=dataset), coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41, data=dataset))

## IGRP vs IG
anova(coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41 + response4 + post, data=dataset), coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41, data=dataset))

## IGR vs IR
anova(coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41 + response4, data=dataset), coxph(Surv(OS/30, event)~IPSSR_score +response4, data=dataset))

## IGRP vs IR
anova(coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41 + response4 + post, data=dataset), coxph(Surv(OS/30, event)~IPSSR_score  +response4, data=dataset))

## IGRP vs IGR
anova(coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41 + response4 + post, data=dataset), coxph(Surv(OS/30, event)~IPSSR_score  + TP53_multi + DDX41 +response4, data=dataset))




res_I<-coxph(Surv(OS/30, event)~IPSSR_score, data=dataset)
res_IG<-coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41, data=dataset)
res_IR<-coxph(Surv(OS/30, event)~IPSSR_score + response4, data=dataset)
res_IGR<-coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41+ response4, data=dataset)
res_IGRP<-coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41+ response4 +post, data=dataset)

loglik = c(res_I$loglik[2], res_IG$loglik[2], res_IR$loglik[2], res_IGR$loglik[2], res_IGRP$loglik[2])
loglik.df <- data.frame(model=c("I","IG","IR","IGR","IGRP"), loglik=loglik)
loglik.df$model<- factor(loglik.df$model, levels=c("I","IG","IR","IGR","IGRP"),ordered = T)


ggplot(data=loglik.df, aes(x=model, y=900+as.numeric(loglik), fill=model)) + geom_bar(stat="identity") + themePPP  + scale_fill_npg() + ylim(c(0,60))
barplot(loglik.df$loglik, ylim=c(-880, -840))


mut_dat.wide$IPSSR<-
  with(mut_dat.wide,
       ifelse(IPSSR_score<=4.5, "VeryLow/Low/Intermediate", ifelse(IPSSR_score<=6.0,"High", "VeryHigh"))
  )
mut_dat.wide$IPSSR<-factor(mut_dat.wide$IPSSR, levels=c("VeryLow/Low/Intermediate", "High","VeryHigh"), order=T)

#library(survminer)
#library(prodlim)
#library(Publish)




### Figure 3D
res.cut <- surv_cutpoint(subset(mut_dat.wide, flag_post_max==1), time = "OS", event = "event", variables = c("post"))
# optimal cutpoint in surv_cutpoint function of maxstat package.  
summary(res.cut)
post_median=res.cut$cutpoint[1,1]


res<-survfit(Surv(OS/30, event)~post>post_median, data=subset(mut_dat.wide, flag_post_max==1))
ggsurvplot(res,  surv.median.line = "n",  pval = TRUE, palette="jco",
           xlab="Time from AZA treatment (months)",title="OS by Max(VAFpost)", ylab="Overall survival",
           risk.table=T, risk.table.height=0.2,
           font.title=c(15, "bold"), legend.labs = c("Max(VAFpost)<=0.6888", "Max(VAFpost)>0.6888"),
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())



res<-survfit(Surv(OS/30, event)~IPSSR, data=subset(mut_dat.wide, flag_post_max==1))
ggsurvplot(res,  surv.median.line = "n",  pval = TRUE, palette="jco",
           xlab="Time from AZA treatment (months)",title="OS by IPSS-R, cases with post samples", ylab="Overall survival",
           risk.table=T, risk.table.height=0.2,
           #palette =c("red","darkgreen","blue"),
           font.title=c(15, "bold"), legend.labs = c("VeryLow/Low/Int", "High","VeryHigh"),
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())


res<-survfit(Surv(OS/30, event)~IPSSR+(post>post_median), data=subset(mut_dat.wide, flag_post_max==1))
g<-ggsurvplot(res,  surv.median.line = "h",  pval = FALSE,
              palette = c("#0073C2","#0073C3", "#EFC000", "#EFC001", "#868686", "#868687"),
              #palette="jco",
              linetype=c(1,2,1,2,1,2), censor=F,
              xlab="Time from AZA treatment (months)",title="OS by IPSS-R, cases with post samples", ylab="Overall survival",
              risk.table=T, risk.table.height=0.3,
              font.title=c(15, "bold"), legend.labs = c("VeryLow/Low/int, post<=threshold","VeryLow/Low/int, post>threshold", "High, post=<threshold","High, post>threshold",
                                                        "VeryHigh, post<=threshold", "VeryHigh, post>threshold"),
              font.x =c(13, "bold"), font.y =c(13, "bold"),
              font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())
g

res1<-survdiff(Surv(OS/30, event)~(post>post_median), data= subset(mut_dat.wide, flag_post_max==1 & IPSSR=="VeryLow/Low/Intermediate"))
p_val1 = 1-pchisq(res1$chisq, 1)
res2<-survdiff(Surv(OS/30, event)~(post>post_median), data= subset(mut_dat.wide, flag_post_max==1 & IPSSR=="High"))
p_val2 = 1-pchisq(res2$chisq, 1)
res3<-survdiff(Surv(OS/30, event)~(post>post_median), data= subset(mut_dat.wide, flag_post_max==1 & IPSSR=="VeryHigh"))
p_val3 = 1-pchisq(res3$chisq, 1)

res4<-survdiff(Surv(OS/30, event)~(IPSSR), data= subset(mut_dat.wide, flag_post_max==1 & ((IPSSR=="VeryLow/Low/Intermediate" & post>post_median) | (IPSSR=="High" & post<=post_median)) ) )
p_val4 = 1-pchisq(res4$chisq, 1)
res5<-survdiff(Surv(OS/30, event)~(IPSSR), data= subset(mut_dat.wide, flag_post_max==1 & ((IPSSR=="High" & post>post_median) | (IPSSR=="VeryHigh" & post<=post_median)) ))
p_val5 = 1-pchisq(res5$chisq, 1)

plot(res, xlab="Time from AZA treatment (months)",main="OS by IPSS-R and max(VAFpost)", ylab="Overall survival",
     col=c("#0073C2","#0073C2", "#EFC000", "#EFC000", "#868686", "#868686"), lty=c(1,2,1,2,1,2), lwd=2, mark.t=F, xaxp=c(0,125,5)
)

lines(c(30,40),c(1.00,1.00), col="#0073C2", lty=1, lwd=2);    text(43,1.00,"VeryLow/Low/Intermediate", adj=0)
lines(c(30,40),c(0.95,0.95), col="#EFC000", lty=1, lwd=2); text(43,0.95,"High", adj=0)
lines(c(30,40),c(0.90,0.90), col="#868686", lty=1, lwd=2);   text(43,0.90,"VeryHigh", adj=0)

lines(c(90,100),c(1.00,1.00), col="black", lty=1, lwd=2); text(105,1.00,"Clone size<=median", adj=0)
lines(c(90,100),c(0.95,0.95), col="black", lty=2, lwd=2); text(105,0.95,"Clone size> median", adj=0)

lines(c(60,70),c(0.8,0.8), col="#0073C2", lty=1, lwd=2);      lines(c(80,90),c(0.8,0.8), col="#0073C2", lty=2, lwd=2);      text(75,0.8,"vs");    text(95,0.80, paste0("p=", round(p_val1,5)), adj=0)
lines(c(60,70),c(0.75,0.75), col="#EFC000", lty=1, lwd=2); lines(c(80,90),c(0.75,0.75), col="#EFC000", lty=2, lwd=2); text(75,0.75,"vs");   text(95,0.75, paste0("p=", round(p_val2,4)), adj=0)
lines(c(60,70),c(0.7,0.7), col="#868686", lty=1, lwd=2);     lines(c(80,90),c(0.7,0.7), col="#868686", lty=2, lwd=2);     text(75,0.7,"vs");    text(95,0.70, paste0("p=", round(p_val3,4)), adj=0)
lines(c(60,70),c(0.65,0.65), col="#0073C2", lty=2, lwd=2);    lines(c(80,90),c(0.65,0.65), col="#EFC000", lty=1, lwd=2); text(75,0.65,"vs");   text(95,0.65, paste0("p=", round(p_val4,4)), adj=0)
lines(c(60,70),c(0.6,0.6), col="#EFC000", lty=2, lwd=2);   lines(c(80,90),c(0.6,0.6), col="#868686", lty=1, lwd=2);     text(75,0.6,"vs");    text(95,0.60, paste0("p=", round(p_val5,4)), adj=0)



### Figure 3D_improve

res<-survfit(Surv(OS/30, event)~IPSSR, data=subset(mut_dat.wide, flag_post_max==1 & IPSSR %in% c("VeryLow/Low/Intermediate", "High")))
ggsurvplot(res,  surv.median.line = "n",  pval = TRUE,
           xlab="Time from AZA treatment (months)",title="OS by IPSS-R, cases with post samples", ylab="Overall survival",
           risk.table=T, risk.table.height=0.2, linetype=1,  censor=F,
           palette =c("#0073C2","#EFC000"),
           font.title=c(15, "bold"), legend.labs = c("VeryLow/Low/Int", "High"),
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())

res<-survfit(Surv(OS/30, event)~IPSSR+(post>post_median), data=subset(mut_dat.wide, flag_post_max==1 & IPSSR %in% c("VeryLow/Low/Intermediate")))
ggsurvplot(res,  surv.median.line = "n",  pval = FALSE,
           xlab="Time from AZA treatment (months)",title="OS by IPSS-R, cases with post samples", ylab="Overall survival",
           risk.table=T, risk.table.height=0.2, linetype=c(2,3), censor=F,
           palette =c("#0073C2","#0073C2"),
           font.title=c(15, "bold"), legend.labs = c("VeryLow/Low/Int", "High"),
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())


res<-survfit(Surv(OS/30, event)~IPSSR, data=subset(mut_dat.wide, flag_post_max==1 & IPSSR %in% c("VeryHigh", "High")))
ggsurvplot(res,  surv.median.line = "n",  pval = TRUE,
           xlab="Time from AZA treatment (months)",title="OS by IPSS-R, cases with post samples", ylab="Overall survival",
           risk.table=T, risk.table.height=0.2, linetype=1,  censor=F,
           palette =c("#EFC000","#868686"),
           font.title=c(15, "bold"), legend.labs = c("VeryLow/Low/Int", "High"),
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())

res<-survfit(Surv(OS/30, event)~IPSSR+(post>post_median), data=subset(mut_dat.wide, flag_post_max==1 & IPSSR %in% c("High")))
ggsurvplot(res,  surv.median.line = "n",  pval = FALSE,
           xlab="Time from AZA treatment (months)",title="OS by IPSS-R, cases with post samples", ylab="Overall survival",
           risk.table=T, risk.table.height=0.2, linetype=c(2,3), censor=F,
           palette =c("#EFC000","#EFC000"),
           font.title=c(15, "bold"), legend.labs = c("VeryLow/Low/Int", "High"),
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())

res<-survfit(Surv(OS/30, event)~IPSSR+(post>post_median), data=subset(mut_dat.wide, flag_post_max==1 & IPSSR %in% c("VeryHigh")))
ggsurvplot(res,  surv.median.line = "n",  pval = FALSE,
           xlab="Time from AZA treatment (months)",title="OS by IPSS-R, cases with post samples", ylab="Overall survival",
           risk.table=T, risk.table.height=0.2, linetype=c(2,3), censor=F,
           palette =c("#868686","#868686"),
           font.title=c(15, "bold"), legend.labs = c("VeryLow/Low/Int", "High"),
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())







res<-survdiff(Surv(days_OS, days_censor)~DDX41>0,data= subset(clini_dat, cohort2=="Test" & IPSS_R_calc=="VeryHigh"))
res


library(patchwork)
mut_dat.wide %>% filter(response %in% c("CR","mCR","PR","SD+HI","SD-HI","PD")) ->mut_dat.wide
g<-ggpaired(mut_dat.wide, cond1="pre", cond2="post", fill= "condition", palette = "jco", 
            line.size = 0.1, line.color="lightgray", facet.by = "response2", 
            width=0.6, title="Change of MCF by response", xlab="Response",ylab="MCF",
            panel.labs = list(response=c("CR","mCR","PR","SD+HI","SD-HI","PD")),
            ncol=6)+stat_compare_means(paired = TRUE, label.y = 1.05)
g1<-ggpar(g,
          legend = "right", legend.title = "Timing",
          font.legend = c(13,"plain", "black"), font.x = c(13, "plain", "black"), ylab=c(0,1))


g1


mut_dat.wide %>% pivot_longer(cols=c("pre", "post") , names_to = "timing", values_to = "VAF") -> mut_dat.tmp
mut_dat.tmp$timing<-factor(mut_dat.tmp$timing, levels=c("pre","post"), ordered = T)
mut_dat.tmp$flags <- ifelse(mut_dat.tmp$flag_pre_max==1 & mut_dat.tmp$flag_post_max==1, "pre_post",
                          ifelse(mut_dat.tmp$flag_pre_max==1 & mut_dat.tmp$flag_post_max==0, "pre",
                                 ifelse(mut_dat.tmp$flag_pre_max==0 & mut_dat.tmp$flag_post_max==1, "post", "none") )  )
mut_dat.tmp$flags <- factor(mut_dat.tmp$flags,levels=c("none","pre","post", "pre_post"),ordered = T)
mut_dat.tmp$flag_TP53 <- ifelse(mut_dat.tmp$TP53_multi==1 & mut_dat.tmp$gene=="TP53","TP53 in TP53multi",
                                ifelse(mut_dat.tmp$TP53_multi==1 & mut_dat.tmp$gene!="TP53","Other in TP53_multi","Others"))
ggplot(mut_dat.tmp, aes(x=timing, y=VAF))+ geom_boxplot(aes(fill=timing), alpha=0.8, col="grey")+ 
  geom_point() + geom_line(aes(group=variantID, col=factor(flags), linetype=factor(TP53_multi)))+ facet_grid(.~response2)+
  scale_fill_jco()+scale_color_jco()+themePPP

g1<-ggplot(subset(mut_dat.tmp,TP53_multi==1 ), aes(x=timing, y=VAF))+ geom_boxplot(aes(fill=timing), alpha=0.6, col="grey")+ 
  geom_point() + geom_line(aes(group=variantID, col=factor(flag_TP53) ))+ facet_grid(.~response2)+
  scale_fill_jco()+scale_color_jco()+themePPP
g2<-ggplot(subset(mut_dat.tmp,TP53_multi==0 ), aes(x=timing, y=VAF))+ geom_boxplot(aes(fill=timing), alpha=0.6, col="grey")+ 
  geom_point() + geom_line(aes(group=variantID, col=factor(flags) ))+ facet_grid(.~response2)+
  scale_fill_jco()+scale_color_jco()+themePPP
g1+g2+plot_layout(ncol=1)

ggplot(mut_dat.tmp, aes(x=timing, y=VAF))+ geom_boxplot(aes(fill=timing), alpha=0.8, col="grey")+ 
  geom_point() + geom_line(aes(group=variantID, col=factor(TP53==0), linetype=factor(TP53_multi) ))+ facet_grid(.~response2)+
  scale_fill_jco()+scale_color_jco()+themePPP

ggplot(subset(mut_dat.tmp, flags %in% c("pre_post","pre") ), aes(x=timing, y=VAF))+ geom_boxplot(aes(fill=timing), alpha=0.6, col="grey")+ 
  geom_point() + geom_line(aes(group=variantID, col=factor(flags) ))+ facet_grid(.~response2)+
  scale_fill_jco()+scale_color_jco()+themePPP

ggplot(subset(mut_dat.tmp, flags %in% c("pre_post","post") ), aes(x=timing, y=VAF))+ geom_boxplot(aes(fill=timing), alpha=0.6, col="grey")+ 
  geom_point() + geom_line(aes(group=variantID, col=factor(flags) ))+ facet_grid(.~response2)+
  scale_fill_jco()+scale_color_jco()+themePPP


## Figure 3A 
ggplot(subset(mut_dat.tmp, flags %in% c("pre_post","post") ), aes(x=timing, y=VAF))+ geom_boxplot(aes(fill=timing), alpha=1, col="black")+ 
  geom_point() + geom_line(aes(group=variantID ))+ facet_grid(.~response2)+
  scale_fill_jco()+scale_color_jco()+themePPP + stat_compare_means(method="t.test", paired=TRUE,label = "p.format", label.y = 1.05)

ggplot(subset(mut_dat.tmp, flags %in% c("none") ), aes(x=timing, y=VAF))+ geom_boxplot(aes(fill=timing), alpha=0.6, col="grey")+ 
  geom_point() + geom_line(aes(group=variantID, col=factor(flags) ))+ facet_grid(.~response2)+
  scale_fill_jco()+scale_color_jco()+themePPP

ggplot( subset(mut_dat.tmp, flags %in% c("pre","none")), aes(x=timing, y=VAF))+ geom_boxplot(aes(fill=timing), alpha=0.6, col="grey")+ 
  geom_point() + geom_line(aes(group=variantID, col=factor(flag_post_max) ))+ facet_grid(.~response2)+
  scale_fill_jco()+scale_color_jco()+themePPP

ggplot(mut_dat.tmp, aes(x=timing, y=VAF))+ geom_boxplot(aes(fill=timing), alpha=0.6, col="grey")+ 
  geom_point() + geom_line(aes(group=variantID, col=factor(flag_post_max) ))+ facet_grid(.~response2)+
  scale_fill_jco()+scale_color_jco()+themePPP +stat_compare_means(label = "p.format")


## Figure 3A candidate2
response_col<-c("red2","indianred1","cornflowerblue","forestgreen","wheat4")
ggplot(subset(mut_dat.tmp, flag_pre_max==1 & timing=="pre"), aes(x=response2, y=VAF, fill=response2))+ geom_boxplot(aes(fill=response2), show.legend = FALSE)+ 
  geom_point(show.legend = FALSE) +   scale_fill_jco()+scale_fill_manual(values = response_col)+themePPP 
with(subset(mut_dat.tmp, flag_pre_max==1 & timing=="pre" & !is.na(response2)), jonckheere.test(VAF, response2, alternative = "increasing"))
with(subset(mut_dat.tmp, flag_pre_max==1 & timing=="pre" & !is.na(response2)), anova(aov(VAF ~ response2)))

ggplot(subset(mut_dat.tmp, flag_post_max==1 & timing=="post"), aes(x=response2, y=VAF, fill=response2))+ geom_boxplot(aes(fill=response2), show.legend = FALSE)+ 
  geom_point(show.legend = FALSE) +   scale_fill_jco()+scale_fill_manual(values = response_col)+themePPP 
with(subset(mut_dat.tmp, flag_post_max==1 & timing=="post" & !is.na(response2)),
     jonckheere.test(VAF, response2, alternative = "increasing", nper=10000))


### cross validation to see the c-index
dataset=subset(mut_dat.wide, flag_post_max==1)
set.seed(1234)




response_col<-c("red2","indianred1","cornflowerblue","forestgreen","wheat4")
g<-ggplot(dataset, aes(x=response2, y=post, fill=response2))
g<-g+geom_violin()
#g<-g+geom_jitter(size = 1, width =0.3)
g<-g+scale_fill_manual(values = response_col)
g<-g+themePPP+ theme(legend.position = 'na')
g<-g+xlab("Response")+ylab("MaxVAFpost")+ggtitle("Group2 ")
g



### I: IPSSR,  G: genetic, R:IWG response, P: post clone

c_test_I<-c()
c_valid_I<-c()

c_test_IG<-c()
c_valid_IG<-c()

c_test_IR<-c()
c_valid_IR<-c()

c_test_IGR<-c()
c_valid_IGR<-c()

c_test_IGRP<-c()
c_valid_IGRP<-c()

for(r in 1:1000){
  dataset$random=runif(nrow(dataset))
  testset=subset(dataset, random<=0.75)
  validset=subset(dataset, random>0.75)
  
  
  
  
  
  ## IPSS-R only
  res <-coxph(Surv(OS/30, event)~IPSSR_score, data=testset)
  c_test_I[r] <-res$concordance[6]
  
  df <- data.frame(t=validset$OS/30, e=validset$event, score=predict(res, validset, type="risk"))
  c_valid_I[r] <-estC(df$t, df$e, (-1)*df$score)
  
  
  
  
  
  ## IPSS-R + IWG response + only
  res <-coxph(Surv(OS/30, event)~IPSSR_score + response4 , data=testset)
  c_test_IR[r] <-res$concordance[6]
  
  df <- data.frame(t=validset$OS/30, e=validset$event, score=predict(res, validset, type="risk"))
  c_valid_IR[r] <-estC(df$t, df$e, (-1)*df$score)
  
  
  
  
  ## IPSS-R + Genetic + only
  res <-coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41 , data=testset)
  c_test_IG[r] <-res$concordance[6]
  
  df <- data.frame(t=validset$OS/30, e=validset$event, score=predict(res, validset, type="risk"))
  c_valid_IG[r] <-estC(df$t, df$e, (-1)*df$score)
  
  
  
  
  ## IPSS-R + genetic + IWG response only
  res <-coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41 + response4, data=testset)
  c_test_IGR[r] <-res$concordance[6]
  
  df <- data.frame(t=validset$OS/30, e=validset$event, score=predict(res, validset, type="risk"))
  c_valid_IGR[r] <-estC(df$t, df$e, (-1)*df$score)
  
  
  
  
  
  ## IPSS-R + genetic + IWG response + post
  res <-coxph(Surv(OS/30, event)~IPSSR_score + TP53_multi + DDX41 + post + response4, data=testset)
  c_test_IGRP[r] <-res$concordance[6]
  
  df <- data.frame(t=validset$OS/30, e=validset$event, score=predict(res, validset, type="risk"))
  c_valid_IGRP[r] <-estC(df$t, df$e, (-1)*df$score)
  
}





##ggplot density plot

df_test<-data.frame(I=c_test_I, IG=c_test_IG, IR=c_test_IR, IGR=c_test_IGR, IGRP=c_test_IGRP)
df_test_long <- df_test %>% pivot_longer(cols=c(I, IR, IG, IGR, IGRP),names_to = "model", values_to = "c-index")
df_test_long$model <-factor(df_test_long$model, levels=c("I", "IG", "IR", "IGR", "IGRP"),ordered = T)
g<-ggplot(data=df_test_long, aes(x=`c-index`, fill=model))
#g<-g+geom_histogram(position="identity", alpha=0.5)
g<-g + geom_density(aes(color = model, alpha = 0.2), show.legend = T)
g<-g+themePPP
g<-g+scale_fill_jco()
g<-g+scale_color_jco()
g<-g+xlim(c(0.6,0.8))
g<-g+xlab("C-index") + ylab("Frequency") + ggtitle("Histogram of C-index, test")
g1<-g

df_valid<-data.frame(I=c_valid_I, IG=c_valid_IG, IR=c_valid_IR, IGR=c_valid_IGR,IGRP=c_valid_IGRP)
df_valid_long <- df_valid %>% pivot_longer(cols=c(I, IR, IG, IGR, IGRP),names_to = "model", values_to = "c-index")
df_valid_long$model <-factor(df_valid_long$model, levels=c("I", "IG", "IR", "IGR", "IGRP"),ordered = T)
g<-ggplot(data=df_valid_long, aes(x=`c-index`, fill=model))
#g<-g+geom_histogram(position="identity", alpha=0.5)
g<-g + geom_density(aes(color = model, alpha = 0.2), show.legend = T)
g<-g+themePPP
g<-g+scale_fill_jco()
g<-g+scale_color_jco()
g<-g+xlim(c(0.5,0.85))
g<-g+xlab("C-index") + ylab("Frequency") + ggtitle("Histogram of C-index, validation")
g2<-g

g1+g2+plot_layout(ncol=2)



##ggplot box plot

g<-ggplot(data=df_test_long, aes(x=model, y=as.numeric(`c-index`), fill=model))
g<-g + geom_boxplot(outlier.shape = NA)
g<-g+themePPP + theme(legend.position = 'none')
g<-g+scale_fill_npg()
g<-g+scale_color_npg()
g<-g+ylim(c(0.45,0.9))
g<-g+xlab("C-index") + ylab("Frequency") + ggtitle("Histogram of C-index, test")
g1<-g

g<-ggplot(data=df_valid_long, aes(x=model, y=as.numeric(`c-index`), fill=model))
g<-g + geom_boxplot(outlier.shape = NA)
g<-g+themePPP + theme(legend.position = 'none')
g<-g+scale_fill_npg()
g<-g+scale_color_npg()
g<-g+ylim(c(0.5,0.8))
g<-g+xlab("C-index") + ylab("Frequency") + ggtitle("Histogram of C-index, test")
g2<-g

g1+g2+plot_layout(ncol=2)


df_valid_long %>% group_by(model) %>% summarize(mean(`c-index`))
df_valid_long %>% group_by(model) %>% summarize(median(`c-index`))


### post max(VAF) and response

response_col<-c("red2","indianred1","cornflowerblue","forestgreen","wheat4")
g<-ggplot(data=subset(mut_dat.wide, !is.na(response2) & flag_post_max==1), aes(x=response2, y=post, fill=response2))
g<-g + geom_boxplot() + themePPP + theme(legend.position = "none") + scale_fill_manual(values=response_col)
g<-g+xlab("Response") + ylab("max(Post)") + ggtitle("Post-treatment max(VAF)")

with(subset(mut_dat.wide, !is.na(response2) & flag_post_max==1),
     jonckheere.test(post, response2, alternative = "increasing", nper=10000))
g<-g+annotate("text",x=5,y=1.1, label="p=0.0001", size=6)
g



### post max(VAF) and response in TP53
response_col<-c("red2","indianred1","cornflowerblue","forestgreen","wheat4")
g<-ggplot(data=subset(mut_dat.wide, !is.na(response2) & flag_post_max==1 & gene %in% c("TP53")), aes(x=response2, y=post, fill=response2))
g<-g + geom_boxplot() + themePPP + theme(legend.position = "none") + scale_fill_manual(values=response_col)
g<-g+xlab("Response") + ylab("max(Post)") + ggtitle("Post-treatment max(VAF) \nTP53")

with(subset(mut_dat.wide, !is.na(response2) & flag_post_max==1 & gene %in% c("TP53")),
     jonckheere.test(post, response2, alternative = "increasing", nper=10000))
g<-g+annotate("text",x=5,y=1.1, label="p=0.0001", size=6)
g

### post max(VAF) and response in splicing factors
response_col<-c("red2","indianred1","cornflowerblue","forestgreen","wheat4")
g<-ggplot(data=subset(mut_dat.wide, !is.na(response2) & flag_post_max==1 & gene %in% c("SF3B1","SRSF2","U2AF1")), aes(x=response2, y=post, fill=response2))
g<-g + geom_boxplot() + themePPP + theme(legend.position = "none") + scale_fill_manual(values=response_col)
g<-g+xlab("Response") + ylab("max(Post)") + ggtitle("Post-treatment max(VAF) \nSplicing factors (SF3B1, SRSF2, U2AF1)")

with(subset(mut_dat.wide, !is.na(response2) & flag_post_max==1 & gene %in% c("SF3B1","SRSF2","U2AF1")),
     jonckheere.test(post, response2, alternative = "increasing", nper=10000))
g<-g+annotate("text",x=5,y=1.1, label="p=0.0001", size=6)
g

### post max(VAF) and response in CH-related
response_col<-c("red2","indianred1","cornflowerblue","forestgreen","wheat4")
g<-ggplot(data=subset(mut_dat.wide, !is.na(response2) & flag_post_max==1 & gene %in% c("ASXL1","TET2","DNMT3A")), aes(x=response2, y=post, fill=response2))
g<-g + geom_boxplot() + themePPP + theme(legend.position = "none") + scale_fill_manual(values=response_col)
g<-g+xlab("Response") + ylab("max(Post)") + ggtitle("Post-treatment max(VAF) \nCH-related genes (ASXL1,TET2,DNMT3A)")

with(subset(mut_dat.wide, !is.na(response2) & flag_post_max==1 & gene %in% c("ASXL1","TET2","DNMT3A")),
     jonckheere.test(post, response2, alternative = "increasing", nper=10000))
g<-g+annotate("text",x=5,y=1.1, label="p=0.012", size=6)
g


### post max(VAF) and response in various genes
response_col<-c("red2","indianred1","cornflowerblue","forestgreen","wheat4")
subset(mut_dat.wide, !is.na(response2) & flag_post_max==1 & gene %in% c("TP53","SF3B1","SRSF2","U2AF1","ASXL1", "TET2","DNMT3A")) -> tmp.df
tmp.df$gene_group <-ifelse(tmp.df$gene=="TP53","TP53", ifelse(tmp.df$gene %in% c("SF3B1","SRSF2","U2AF1"), "Splicing F", ifelse(tmp.df$gene %in% c("ASXL1", "TET2","DNMT3A"),"CH-related","NA")))
tmp.df$gene_group <- factor(tmp.df$gene_group, levels=c("TP53","Splicing F","CH-related"), ordered = T) 
g<-ggplot(data=subset(tmp.df, response2 %in% c("CR","mCR/PR")), aes(x=gene_group, y=post, fill=response2))
g<-g + geom_boxplot() + themePPP + scale_fill_manual(values=response_col)
g<-g+xlab("Response") + ylab("max(Post)") + ggtitle("Post-treatment max(VAF)")

#with(subset(mut_dat.wide, !is.na(response2) & flag_post_max==1 & gene %in% c("TET2","DNMT3A")), jonckheere.test(post, response2, alternative = "increasing", nper=10000))
#g<-g+annotate("text",x=5,y=1.1, label="p=0.012", size=6)
g




### post max(VAF) and response by gene
mut_dat.wide %>% filter(flag_post_max==1) %>% distinct(patient, .keep_all=TRUE) -> mut_dat.wide.post_max
mut_dat.wide.post_max %>% group_by(gene) %>% summarise(n()) -> gene_count 
gene_count %>% arrange(desc(`n()`)) %>% filter(`n()`>5) -> gene_count

response_col<-c("red2","indianred1","cornflowerblue","forestgreen","wheat4")
g<-ggplot(data=subset(mut_dat.wide, !is.na(response2) & gene %in% subset(gene_count, `n()`>5)$gene & flag_post_max==1), aes(x=response2, y=post, fill=response2))
g<-g + geom_boxplot() + facet_wrap(~ gene, ncol=4) + themePPP + theme(legend.position = "none") + scale_fill_manual(values=response_col)
g<-g+xlab("Response") + ylab("max(Post)") + ggtitle("Post-treatment max(VAF), genes n>=7")
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g<-g+theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
g










########         Max(VAF), p= prob for H0: delta==0
plot(c(0, 1),c(-0.3,6),col="white", xlab="Median(MAX_VAF)", ylab="-Log10 (p-value)", main="Median(MAX_postVAF)")
lines(c(0, 1),c(0,0),col="gray",lty=1)
lines(c(0.5, 0.5),c(0,5),col="gray",lty=1)

lines(c(0, 1), rep(log10(0.05)*(-1),2), col="lightgray", lty=2)
#lines(c(-0.5, 0.5), rep(log10(0.1)*(-1),2), col="lightgray", lty=2)
pvalues<-c()


genes <-unique(subset(mut_dat.wide, flag_post_max==1)$gene)
for(g in 1:length(genes)){
  pvalues[g]<-"NA"
  mut_dat.wide %>% filter(gene==genes[g] & flag_post_max==1) -> tmp.df
  print(as.character(genes[g]))
  if(nrow(tmp.df)>4){
    #res<-t.test(tmp.df$pre, tmp.df$post, paired = T)
    res<-jonckheere.test(tmp.df$post, tmp.df$response2, alternative="increasing", nperm=10000)
    y=(-log(res$p.value))/log(10); pvalues[g]<-res$p.value
    x=(median(tmp.df$post))
    r=sqrt(nrow(tmp.df)/20000)
    b_col<-rgb(1,0,0,0.8)
    col <- rgb(1,0,0,0.2)
    draw.circle(x,y,radius=r,nv=100, border=b_col, col =col, lty=1,lwd=1)
    text(x+r,y,as.character(genes[g]),adj=0)
  }
}

draw.circle(0.8,5.6,radius=sqrt(50/20000),nv=100, border=b_col, col =col, lty=1,lwd=1); text(0.9,5.6,"n=50")
draw.circle(0.8,5,radius=sqrt(25/20000),nv=100, border=b_col, col =col, lty=1,lwd=1); text(0.9,5,"n=25")
draw.circle(0.8,4.6,radius=sqrt(10/20000),nv=100, border=b_col, col =col, lty=1,lwd=1); text(0.9,4.6,"n=10")
draw.circle(0.8,4.2,radius=sqrt(5/20000),nv=100, border=b_col, col =col, lty=1,lwd=1); text(0.9,4.2,"n=5")

q.values <- p.adjust(pvalues, method = "BH")
data.frame(gene=genes, q.values=q.values, p.values=pvalues)->tmp.df
subset(tmp.df, q.values<0.10)
subset(tmp.df, p.values<0.05)



