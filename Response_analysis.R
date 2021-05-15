setwd(PWD)

library(RColorBrewer)
library(tidyverse)
library(readxl)
library(forestmodel)
library(plotrix)
library(exact2x2)
library(ggsci)
library(glmnet)
library(coefplot)
library(rpubs)


#simple profile of initial response
table(clini_dat$HI_E, clini_dat$cohort2)
table(clini_dat$HI_P, clini_dat$cohort2)
table(clini_dat$HI_N, clini_dat$cohort2)
table(clini_dat$HI_any, clini_dat$cohort2) ### IH_any == 1 when either of HI is (-1,1)
table(clini_dat$HI_any2, clini_dat$cohort2) ### IH_any == 1 when either of HI is 1




## response profile
library(ggalluvial)


data.frame(best=subset(clini_dat, cohort2=="Test")$response_best,
           last=subset(clini_dat, cohort2=="Test")$response_last) -> response.df

response.df %>% group_by(best,last) %>% summarise(count=n()) ->response.df2
response.df2 %>% filter(!is.na(best) ) ->response.df2


themePPP = theme_classic() + theme(text=element_text(size=15,colour="black", family=sans),
                                   axis.text=element_text(size=15,colour="black"),
                                   axis.title=element_text(size=15,colour="black"),
                                   strip.text.x = element_text(size=15),
                                   legend.text=element_text(size=15))
                                   
response_col<-c("red2","indianred1","cornflowerblue","wheat4")

g<-ggplot(response.df2, aes(y=count, axis1=best, axis2=last))+
  geom_alluvium(aes(fill = last),
                width = 1/4,knot.pos = 0, reverse = FALSE)+
  geom_stratum(width = 1/4, reverse = FALSE) +
  #geom_label(stat="stratum", infer.label=TRUE)+
  geom_text(stat = "stratum", infer.label = TRUE, reverse = FALSE) +
  scale_fill_manual(values=response_col)+themePPP+
  ggtitle("Best and last response")
g                                   
                                   






###  HI factor analysis, dichotomous analysis
cohort_for_test="Test"; select="none"; endpoint="HI_E"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1
cohort_for_test="Validation"; select="none"; endpoint="HI_E"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1

cohort_for_test="Test"; select="none"; endpoint="HI_N"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.15
cohort_for_test="Validation"; select="none"; endpoint="HI_N"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1

cohort_for_test="Test"; select="none"; endpoint="HI_P"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1
cohort_for_test="Validation"; select="none"; endpoint="HI_P"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1

cohort_for_test="Test"; select="none"; endpoint="HI_any"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1
cohort_for_test="Validation"; select="none"; endpoint="HI_any"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1

cohort_for_test="Test"; select="none"; endpoint="HI_any2"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1
cohort_for_test="Validation"; select="none"; endpoint="HI_any2"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1




data_both$CR <- ifelse(data_both$response_best %in% c("CR") ,1,ifelse(is.na(data_both$response_best), NA, 0))
cohort_for_test="Test"; select="none"; endpoint="CR"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1
cohort_for_test="Validation"; select="none"; endpoint="CR"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1

cohort_for_test="Test"; select="TP53_multi"; endpoint="CR"; N_cutoff_for_multi=4; p_cutoff_for_multi=0.1
cohort_for_test="Test"; select="noTP53_multi"; endpoint="CR"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1



data_both$Remission <- ifelse(data_both$response_best %in% c("CR", "mCR", "PR") ,1,ifelse(is.na(data_both$response_best), NA, 0))
cohort_for_test="Test"; select="none"; endpoint="Remission"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1
cohort_for_test="Validation"; select="none"; endpoint="Remission"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1


#data_both$benefit <- ifelse(data_both$response_best2 %in% c("CR", "mCR", "PR", "SD+HI"),1,ifelse(is.na(data_both$response_best), NA, 0))
#cohort_for_test="Test"; select="none"; endpoint="benefit"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.2
#cohort_for_test="Validation"; select="none"; endpoint="benefit"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1



data_both$benefit <- ifelse(data_both$response_best3 %in% c("CR", "mCR", "PR", "SD+HI"),1,ifelse(is.na(data_both$response_best), NA, 0))
cohort_for_test="Test"; select="none"; endpoint="benefit"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1
cohort_for_test="Validation"; select="none"; endpoint="benefit"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1

cohort_for_test="Test"; select="TP53_multi"; endpoint="benefit"; N_cutoff_for_multi=4; p_cutoff_for_multi=0.1
cohort_for_test="Test"; select="noTP53_multi"; endpoint="benefit"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1


data_both$PD <- ifelse(data_both$response_best %in% c("PD") ,1,ifelse(is.na(data_both$response_best), NA, 0))
cohort_for_test="Test"; select="none"; endpoint="PD"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1
cohort_for_test="Validation"; select="none"; endpoint="PD"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1


##early PD
data_both$PD_early <- ifelse(data_both$response_best %in% c("PD") & data_both$days_progress_failure<83 ,1,ifelse(is.na(data_both$response_best), NA, 0))
cohort_for_test="Test"; select="none"; endpoint="PD_early"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1
cohort_for_test="Validation"; select="none"; endpoint="PD_early"; N_cutoff_for_multi=0.05; p_cutoff_for_multi=0.1


{
  ## select cohort (Test / Validation) 
  dat_both <-subset(data_both, cohort2==cohort_for_test)
  
  if(select=="TP53_multi"){
    dat_both <-subset(dat_both, TP53_multi==1)
  }else if(select=="noTP53_multi"){
    dat_both <-subset(dat_both, TP53_multi==0)
  }else if(select=="none"){
    dat_both <- dat_both
  }
  
  ## select outcome parameter
  command_txt = paste0("dat_both <- subset(dat_both,", endpoint, " %in% c(0,1))")
  eval(parse(text=command_txt))
  
  command_txt = paste0("dat_both %>% mutate(response = as.numeric(as.character(dat_both$", endpoint,"))) -> dat_both")
  eval(parse(text=command_txt))
  
  
  dat_both$response <-factor(dat_both$response, levels=c(0,1))
  
  
  ## response profile by gene mutation => for graph of response by gene mutation
  response_ta<-c()
  for(g in c(clinic_list_univariate, genetic_list2)){
    
    ta<-c()
    ta$gene <- g
    
    table<-matrix(c(0,0,0,0),nrow=2)
    colnames(table) <- c("TRUE","FALSE")
    rownames(table) <- c(1,0)
    
    paste0('table["0","TRUE"]<-nrow(subset(dat_both, response==0 &', g, ">0))") -> command_txt; eval(parse(text=command_txt))
    paste0('table["1","TRUE"]<-nrow(subset(dat_both, response==1 &', g, ">0))") -> command_txt; eval(parse(text=command_txt))
    paste0('table["0","FALSE"]<-nrow(subset(dat_both, response==0 &', g, "==0))") -> command_txt; eval(parse(text=command_txt))
    paste0('table["1","FALSE"]<-nrow(subset(dat_both, response==1 &', g, "==0))") -> command_txt; eval(parse(text=command_txt))

    ta$par1_res0 <- table["0","TRUE"]
    ta$par1_res1 <- table["1","TRUE"]
    ta$par0_res0 <- table["0","FALSE"]
    ta$par0_res1 <- table["1","FALSE"]
    
    fisher.exact(table)->res
    
    ta$p_value <- res$p.value
    ta$odds <- res$estimate[[1]]
    ta$lower <- res$conf.int[1]
    ta$upper <- res$conf.int[2]
    
    response_ta <-rbind(response_ta, ta)
    
  }
  write.table(response_ta, paste0("out/",endpoint,"_",cohort_for_test,".txt"), sep="\t", col.names=T, row.names = F)
  
  
  
  #outpath=paste0("out/", endpoint,"_",cohort_for_test, ".pdf")
  #pdf(outpath, height=12, width=9)

  ##draw univariate
  as.data.frame(response_ta) %>% filter(as.numeric(par1_res0) + as.numeric(par1_res1) > 
                N_cutoff_for_multi*(as.numeric(par1_res0)+as.numeric(par1_res1)+as.numeric(par0_res0)+as.numeric(par0_res1))) -> tmp.df_draw
  tmp.df_draw$x <- log(as.numeric(tmp.df_draw$odds))/log(2)
  tmp.df_draw$y <- -log(as.numeric(tmp.df_draw$p_value))/log(10)
  
  
  theta <- seq(-pi, pi, length=100)
  
  title=paste0(endpoint," ",cohort_for_test)
  
  odds <- as.numeric(tmp.df_draw$odds)
  odds <- odds[odds>0]
  xlim = max(abs(log(odds) / log(2)), na.rm=T)*1.3
  ylim = (log(min(as.numeric(tmp.df_draw$p_value),na.rm = T)) / log(10)) *(-1.1)

  plot(c(-1*xlim, xlim),c(-0.3,ylim),col="white", xlab="Log2(Odds ratio)", ylab="-Log10(P-value)", main=title)
  #plot(c(-1*xlim, xlim),c(-0.3,2),col="white", xlab="Log2(Odds ratio)", ylab="-Log10(P-value)", main=title)
  lines(c(-xlim, xlim),c(0,0),col="darkgray",lty=1)
  lines(c(0, 0),c(0,40),col="darkgray",lty=1)
  lines(c(-xlim, xlim), rep(log10(0.05)*(-1),2), col="lightgray", lty=2)
  lines(c(-xlim, xlim), rep(log10(p_cutoff_for_multi)*(-1),2), col="lightgray", lty=2)
  
  
  
  for(p in 1:nrow(tmp.df_draw)){
    x<-tmp.df_draw$x[p]
    y<-tmp.df_draw$y[p]
    gene<-tmp.df_draw$gene[p]
    
    if(gene %in% (colnames(clini_dat_clinic_univariate))){
      r<-sqrt(as.numeric(tmp.df_draw$par1_res0[p]) + as.numeric(tmp.df_draw$par1_res1[p]))/40
      b_col<-rgb(1,0,0,0.8)
      col <- rgb(1,0,0,0.2)
    }else if(gene %in% c(genetic_list_SNV, genetic_allele_list)){
      r<-sqrt(as.numeric(tmp.df_draw$par1_res0[p]) + as.numeric(tmp.df_draw$par1_res1[p]))/40
      b_col<-rgb(0,1,0,0.8)
      col <- rgb(0,1,0,0.2)
    }else if(gene %in% genetic_list_CNV){
      r<-sqrt(as.numeric(tmp.df_draw$par1_res0[p]) + as.numeric(tmp.df_draw$par1_res1[p]))/40
      b_col<-rgb(0,0,1,0.8)
      col <- rgb(0,0,1,0.2)
    }
    
    sub("del_upd_", "del/upd ", gene) -> gene
    
    draw.circle(x,y,radius=r,nv=100, border=b_col, col =col, lty=1,lwd=1)
    if( y>(log10(p_cutoff_for_multi)*(-1))){text(x+r+0.1, y, gene, adj=0, cex=0.7)}
  }
  
  
  ### for legend
  
  plot(c(-1*xlim, xlim),c(-0.3,ylim),col="white", xlab="Log2(Odds ratio)", ylab="-Log10(P-value)", main=title)
  #plot(c(-1*xlim, xlim),c(-0.3,2),col="white", xlab="Log2(Odds ratio)", ylab="-Log10(P-value)", main=title)
  
  
  r<-0.1
  b_col<-rgb(1,0,0,0.8)
  col <- rgb(1,0,0,0.2)
  draw.circle(0,0.1,radius=r,nv=100, border=b_col, col =col, lty=1,lwd=1)
  text(0.3, 0.1, "Clinical factors", adj=0)
  
  r<-0.1
  b_col<-rgb(0,1,0,0.8)
  col <- rgb(0,1,0,0.2)
  draw.circle(0,0.2,radius=r,nv=100, border=b_col, col =col, lty=1,lwd=1)
  text(0.3, 0.2, "Mutations", adj=0)
  
  r<-0.1
  b_col<-rgb(0,0,1,0.8)
  col <- rgb(0,0,1,0.2)
  draw.circle(0,0.3,radius=r,nv=100, border=b_col, col =col, lty=1,lwd=1)
  text(0.3, 0.3, "CNAs", adj=0)
  
  r<- sqrt(10)/40
  b_col<-rgb(0,0,0,0.8)
  col <- rgb(0,0,0,0.2)
  draw.circle(0,0.4,radius=r,nv=100, border=b_col, col =col, lty=1,lwd=1)
  text(0.3, 0.4, "n = 10", adj=0)
  
  r<- sqrt(20)/40
  b_col<-rgb(0,0,0,0.8)
  col <- rgb(0,0,0,0.2)
  draw.circle(0,0.45,radius=r,nv=100, border=b_col, col =col, lty=1,lwd=1)
  text(0.3, 0.45, "n = 20", adj=0)  
  
  r<- sqrt(40)/40
  b_col<-rgb(0,0,0,0.8)
  col <- rgb(0,0,0,0.2)
  draw.circle(0,0.5,radius=r,nv=100, border=b_col, col =col, lty=1,lwd=1)
  text(0.3, 0.5, "n = 40", adj=0) 
  
  r<- sqrt(80)/40
  b_col<-rgb(0,0,0,0.8)
  col <- rgb(0,0,0,0.2)
  draw.circle(0,0.6,radius=r,nv=100, border=b_col, col =col, lty=1,lwd=1)
  text(0.3, 0.6, "n = 80", adj=0) 
  
  r<- sqrt(160)/40
  b_col<-rgb(0,0,0,0.8)
  col <- rgb(0,0,0,0.2)
  draw.circle(0,0.7,radius=r,nv=100, border=b_col, col =col, lty=1,lwd=1)
  text(0.3, 0.7, "n = 160", adj=0) 
  #dev.off()
  
  
  # carry the result on to multivariate analysis
 
  
  #dev.off()
  
  
  
  ## elastic net
  
  #explanation parameters
  dat_both[,colnames(dat_both) %in% as.character(tmp.df2$gene)] ->x
  dat_both$response ->y
  mse.df<-NULL
  
  alpha <- seq(0.01, 0.99, 0.01)
  for(i in 1:length(alpha)){
    m<- cv.glmnet(as.matrix(x), y, family="binomial", nfolds = 5, alpha=alpha[i])
    mse.df<-rbind(mse.df, data.frame(alpha = alpha[i], mse = min(m$cvm)))
  }
  
  best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
  m<- cv.glmnet(as.matrix(x), y, family="binomial", nfolds = 5, alpha=best.alpha)
  best.lambda <- m$lambda.min
  
  model<- glmnet(as.matrix(x), y, family="binomial", nfolds = 5, lambda = best.lambda, alpha=best.alpha)
  coefplot(model, sort="magnitude", cex=2)
  

  #dev.off()
  
  
  coef(model, s="lambda.min") ->coefdat
  cbind(name=coefdat@Dimnames[[1]][coefdat@i+1], coef=coefdat@x) ->df.coef
  write.table(df.coef, "out/elasticnet.coeff.txt", sep="\t")
  
}


panels<-list(list(width = 0.03),
             list(width = 0.1, display = ~variable, fontface = "bold", heading = "Variable"),
             list(width = 0.1, display = ~level),
             list(width = 0.05, display = ~n, hjust = 1, heading = "N"),
             list(width = 0.03, item = "vline", hjust = 0.5),
             list(width = 0.30, item = "forest", hjust = 0.5, heading = "Odds ratio", linetype = "dashed",
                  line_x = 0),
             list(width = 0.03, item = "vline", hjust = 0.5),
             list(width = 0.12, display = ~ifelse(reference, "Reference", sprintf("%0.2f (%0.2f, %0.2f)",
                                                                                  trans(estimate), trans(conf.low), trans(conf.high))), display_na = NA),
             list(width = 0.03, item = "vline", hjust = 0.5),
             list(width = 0.05,
                  display = ~ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
                  display_na = NA, hjust = 1, heading = "p"),
             list(width = 0.03)
)





###########################################################################################################################################
res <-glm(CR ~  age_elder + HB_lower + isKTrisk_higher + TP53_multi + ASXL1 + EZH2 + STAG2 + del_7_7q + del_5_5q + del_upd_17p, binomial(), data=dat_both)
(car::vif(res) ->tmp)
tmp %>% as.data.frame()-> tmp2
colnames(tmp2)<-"gvif"
tmp2 %>% mutate(item=row.names(tmp2)) -> tmp2
tmp2 %>% ggplot(aes(x=item,y=gvif))+geom_bar(stat="identity")+themePPP+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ylim(c(0,5))
## remove LOH17p 

res <-glm(CR ~  age_elder + HB_lower + isKTrisk_higher + TP53_multi + ASXL1 + EZH2 + STAG2 + del_7_7q + del_5_5q , binomial(), data=dat_both)
(car::vif(res) ->tmp)
tmp %>% as.data.frame()-> tmp2
colnames(tmp2)<-"gvif"
tmp2 %>% mutate(item=row.names(tmp2)) -> tmp2
tmp2 %>% ggplot(aes(x=item,y=gvif))+geom_bar(stat="identity")+themePPP+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ylim(c(0,5))

summary(step(res))->res2
command_txt=paste0("res2 <- glm(", res2$call[2], ", binomial(), data=dat_both)")
command_txt
eval(parse(text=command_txt))

forest_model(res)
forest_model(res2, panels)
#forest_model(res2, format_options = list(color="black", text_size=4.5), panels )





panels<-list(list(width = 0.03),
             list(width = 0.1, display = ~variable, fontface = "bold", heading = "Variable"),
             list(width = 0.1, display = ~level),
             list(width = 0.03, item = "vline", hjust = 0.5),
             list(width = 0.30, item = "forest", hjust = 0.5, heading = "Odds ratio", linetype = "dashed",
                  line_x = 0),
             list(width = 0.03, item = "vline", hjust = 0.5),
             list(width = 0.14, display = ~ifelse(reference, "Reference", sprintf("%0.2f (%0.2f, %0.2f)",
                                                                                  trans(estimate), trans(conf.low), trans(conf.high))), display_na = NA),
             list(width = 0.03, item = "vline", hjust = 0.5),
             list(width = 0.05,
                  display = ~ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
                  display_na = NA, hjust = 1, heading = "p"),
             list(width = 0.03)
)


forest_model(res)
forest_model(res2, panels)
#forest_model(res2, format_options = list(color="black", text_size=4.5), panels )


data_both$CR <- ifelse(data_both$response_best %in% c("CR") ,1,ifelse(is.na(data_both$response_best), NA, 0))
dat_both <-subset(data_both, cohort2=="Test")
dat_both <-subset(dat_both, CR %in% c(0,1))

res <-glm(CR ~  HB_lower + isKTrisk_higher + TP53 + ASXL1 + del_5_5q + del_upd_17p, binomial(), data=dat_both)
car::vif(res)
res <-glm(CR ~  HB_lower + isKTrisk_higher + TP53 + ASXL1 + del_5_5q  , binomial(), data=dat_both)
car::vif(res)
summary(step(res))->res2
command_txt=paste0("res2 <- glm(", res2$call[2], ", binomial(), data=dat_both)")
command_txt
eval(parse(text=command_txt))

forest_model(res)
forest_model(res2, panels)
#forest_model(res2, format_options = list(color="black", text_size=4.5), panels )



## response duration by PFS
dat_both %>% filter(CR==1) %>% select(TP53_mono, TP53_multi, days_best, days_progress_failure, days_OS, days_censor) -> tmp.df

tmp.df$censor_duration <-ifelse(is.na(tmp.df$days_progress_failure), tmp.df$days_censor,1)
tmp.df$days_duration <-ifelse(is.na(tmp.df$days_progress_failure), tmp.df$days_OS-tmp.df$days_best, tmp.df$days_progress_failure-tmp.df$days_best)
res<-survfit(Surv(days_duration/30, censor_duration)~ (TP53_multi+TP53_mono)>0, data=tmp.df)
ggsurvplot(res,  surv.median.line = "hv",  pval = TRUE,
           xlab="Time from CR (months)",title="PFS in CR cases", ylab="PFS",
           risk.table=T, risk.table.height=0.25, palette="npg",
           font.title=c(15, "bold"), legend.labs = c("TP53 wt", "TP53 mut" ),
           font.x =c(13, "bold"), font.y =c(13, "bold"),
           font.tickslab=c(13), font.legend=c(13), tables.theme = theme_cleantable())
surv_median(res)
res<-survdiff(Surv(days_duration/30, censor_duration)~ (TP53_multi+TP53_mono)>0, data=tmp.df)
res


###########################################################################################################################################









###########################################################################################################################################
# validation cohort for CR
data_both$CR <- ifelse(data_both$response_best %in% c("CR") ,1,ifelse(is.na(data_both$response_best), NA, 0))
dat_both <-subset(data_both, cohort2=="Validation")
dat_both <-subset(dat_both, CR %in% c(0,1))
res <-glm(CR ~  HB_lower + isKTrisk_higher + TP53 + ASXL1 +  del_5_5q + del_upd_17p, binomial(), data=dat_both)
car::vif(res)
forest_model(res)
res$coefficients


library(metafor)
a <- c(9,6,6,7,5,1,4,3,1,0)
n1 <- c(42,27,27,40,29,7,28,22,15,7)
c <- c(4,7,7,6,8,12,9,10,12,13)
n0 <- c(45,60,60,47,58,80,59,65,72,80)
dat <- data.frame(a,n1,c,n0)
dat.escalc <- escalc(measure="OR", ai=a, n1i=n1, ci=c, n2i=n0, data=dat)
res.reml <- rma.uni(yi, vi, method="REML", data=dat.escalc)
forest(res.reml)
###########################################################################################################################################


########  specific calculation  for PD cases
### after univariate analysis with "Test", "PD" 

data_both$response <- ifelse(data_both$response_best %in% c("PD") ,1,ifelse(is.na(data_both$response_best), NA, 0))
dat_both <-subset(data_both, cohort2=="Test")
dat_both <-subset(dat_both, response %in% c(0,1))

#res <-glm(as.numeric(response) ~ ANC_lower + HB_lower + isKTrisk_higher + STAG2 + del_12p + del_7_7q, binomial(), data=dat_both)
res <-glm(as.numeric(response) ~ age_elder + ANC_lower + HB_lower + isKTrisk_higher + STAG2 + del_12p + del_7_7q, binomial(), data=dat_both)
summary(step(res))->res2

command_txt=paste0("res2 <- glm(", res2$call[2], ",binomial(), data=dat_both)")
command_txt
eval(parse(text=command_txt))

forest_model(res)
#forest_model(res2, format_options = list(color="black", text_size=4.5) )
panels<-list(list(width = 0.03),
             list(width = 0.1, display = ~variable, fontface = "bold", heading = "Variable"),
             list(width = 0.1, display = ~level),
             list(width = 0.03, item = "vline", hjust = 0.5),
             list(width = 0.30, item = "forest", hjust = 0.5, heading = "Odds ratio", linetype = "dashed",
                  line_x = 0),
             list(width = 0.03, item = "vline", hjust = 0.5),
             list(width = 0.14, display = ~ifelse(reference, "Reference", sprintf("%0.2f (%0.2f, %0.2f)",
                                                                                  trans(estimate), trans(conf.low), trans(conf.high))), display_na = NA),
             list(width = 0.03, item = "vline", hjust = 0.5),
             list(width = 0.05,
                  display = ~ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
                  display_na = NA, hjust = 1, heading = "p"),
             list(width = 0.03)
)
forest_model(res2,panels)



##########################  TP53 and response, Test 


###  TP53 multi vs others
response_col<-c("red2","indianred1","orange","cornflowerblue","forestgreen","wheat4")
dat_both<-subset(data_both, cohort2=="Test")
ta<-table(dat_both$response_best3, dat_both$TP53_multi)
data.frame(ta) ->ta
colnames(ta) <-c("Response","TP53 multihit","number")
ta$Response <-factor(ta$Response, levels=rev(levels(ta$Response)))

g<-ggplot(ta, aes(x=`TP53 multihit`, y=number, fill=Response))
g<-g+geom_bar(stat="identity", position = "fill") + themePPP + theme(legend.position="bottom")
g<-g+xlab("TP53 multihit")+ylab("Ratio")+ggtitle("Response by TP53 status")
g<-g + scale_fill_manual(values=rev(response_col))
g


###  TP53 mut vs others
dat_both<-subset(data_both, cohort2=="Test")
ta<-table(dat_both$response_best3, dat_both$TP53>0)
data.frame(ta) ->ta
colnames(ta) <-c("Response","TP53","number")
ta$Response <-factor(ta$Response, levels=rev(levels(ta$Response)))

g<-ggplot(ta, aes(x=TP53, y=number, fill=Response))
g<-g+geom_bar(stat="identity", position = "fill") + themePPP + theme(legend.position="bottom")
g<-g+xlab("TP53 mutation")+ylab("Ratio")+ggtitle("Response by TP53 status")
g<-g + scale_fill_manual(values=rev(response_col))
g

############################################################################################
###  TP53 allelic status
dat_both<-subset(clini_dat, cohort2=="Test")
ta<-table(dat_both$response_best3, dat_both$TP53_allelic_status2)
data.frame(ta) ->ta
colnames(ta) <-c("Response","TP53 allelic status","number")
ta$Response <-factor(ta$Response, levels=(levels(ta$Response)))
ta$`TP53 allelic status`<-factor(ta$`TP53 allelic status`, levels=c("multihit","1hit","wt"), ordered = T)

g<-ggplot(ta, aes(x=`TP53 allelic status`, y=number, fill=Response))
g<-g+geom_bar(stat="identity", position = "fill") + themePPP + theme(legend.position="bottom")
g<-g+xlab("TP53 status")+ylab("Ratio")+ggtitle("Response by TP53 status")
g<-g + scale_fill_manual(values=(response_col))
g
############################################################################################










