

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



### Appendix table A2

dim(mut_dat_pre)
mut_dat_pre %>% filter(Sample %in% subset(clini_dat, cohort2=="Test")$before_filename) -> mut_dat_pre_discovery
mut_dat_pre %>% filter(Sample %in% subset(clini_dat, cohort2=="Validation")$before_filename) -> mut_dat_pre_validation
mut_dat_post %>% filter(Sample %in% subset(clini_dat, cohort2=="Test")$after_filename) -> mut_dat_post_discovery
mut_dat_post %>% filter(Sample %in% subset(clini_dat, cohort2=="Validation")$after_filename) -> mut_dat_post_validation

pre_convert = data.frame(patient=clini_dat$patientID, Sample=clini_dat$before_filename)
rbind(
mut_dat_pre_discovery[,c(1,8,96,113)] %>% left_join(pre_convert, by="Sample") %>% mutate(mode="Pre_discovery"),
mut_dat_pre_validation[,c(1,8,96,113)] %>% left_join(pre_convert, by="Sample") %>% mutate(mode="Pre_validation")) -> tmp.df1

post_convert = data.frame(patient=clini_dat$patientID, Sample=clini_dat$after_filename)
rbind(
  mut_dat_post_discovery[,c(1,8,96,113)] %>% left_join(post_convert, by="Sample") %>% mutate(mode="Post_discovery"),
  mut_dat_post_validation[,c(1,8,96,113)] %>% left_join(post_convert, by="Sample") %>% mutate(mode="Post_validation")) -> tmp.df2

rbind(tmp.df1, tmp.df2) -> tmp.df
tmp.df$aVAF <-ifelse(tmp.df$aVAF>1,1,tmp.df$aVAF)

data.frame(patient=clini_dat$patientID, cohort=clini_dat$cohort2, annonymous=paste0("UPN_#",c(1:nrow(clini_dat)))) -> convert
  
  
write.table(tmp.df %>% left_join(convert, by="patient"), "out/mut_dat_simple.txt", sep="\t", row.names = F, col.names = T)





##################################################### Figure 1 Clinical response


## response profile by gene mutation => for graph of response by gene mutation
themePPP = theme_classic() + theme(text=element_text(size=15,colour="black", family=sans),
                                   axis.text=element_text(size=15,colour="black"),
                                   axis.title=element_text(size=15,colour="black"),
                                   strip.text.x = element_text(size=15),
                                   legend.text=element_text(size=15)
)
response_col<-c("red2","indianred1","orange","cornflowerblue","forestgreen","wheat4")

draw<-function(){
  dat_both <-subset(data_both, cohort2==cohort_for_test)
  
  response_ta<-c()
  genes<-c()
  #for(g in c(genetic_list_SNV, genetic_list_CNV)){
  for(g in c(genetic_list2)){
    
    command_txt = paste0("subset(dat_both, ", g, ">0 ) -> sub_clini_dat")
    eval(parse(text=command_txt))
    
    table(sub_clini_dat$response_best3)->ta
    genes <- c(genes,g)
    response_ta <-rbind(response_ta, ta)
    
  }
  total<-apply(response_ta,1,sum)
  cbind(genes, response_ta,total) ->response_ta
  as.data.frame(response_ta)->response_ta
  response_ta %>% arrange(desc(as.numeric(as.character(total)))) -> response_ta
  response_ta$genes <-factor(response_ta$genes, levels=response_ta$genes, order=T)
  response_ta %>% pivot_longer(col = -c(genes,total), names_to="response", values_to="number") ->response_ta2
  response_ta2$response <-factor(response_ta2$response, levels=c("PD","SD-HI","SD+HI","PR","mCR","CR"),order=T)
  response_ta2$number  <-as.numeric(as.character(response_ta2$number))
  response_ta2$total  <-as.numeric(as.character(response_ta2$total))
  
  response_ta2 %>% filter(total> 10) -> response_ta2
  
  
  g<-ggplot(data=response_ta2, aes(x=reorder(genes, (-1)*total), y=number, fill=response))
  g<-g+ geom_bar(stat = "identity")+themePPP
  g<-g+xlab("Genes")+ylab("Number")+ggtitle(paste("Response by genes,", cohort_for_test," cohort",sep=""))
  g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                        vjust = 0.5 , size = rel(1.0), face = "italic"))
  g1<-g + scale_fill_manual(values=rev(response_col))
  g1
  
  g<-ggplot(data=subset(response_ta2, total>5), aes(x=reorder(genes,  (-1)*total), y=number, fill=response))
  g<-g+ geom_bar(stat = "identity", position="fill")+themePPP
  g<-g+xlab("Genes")+ylab("Ratio")+ggtitle(paste("Response by genes,", cohort_for_test," cohort",sep=""))
  g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                        vjust = 0.5 , size = rel(1.0), face = "italic"))
  g<-g + scale_fill_manual(values=rev(response_col))
  g
}

cohort_for_test="Test"; draw()
cohort_for_test="Validation"; draw()





### Figure 1


###  response by gene, TP53-mutated with highlight
cohort_for_test="Validation"
cohort_for_test="Test"
dat_both <-subset(data_both, cohort2==cohort_for_test)

response_ta<-c()
response_ta_withTP53<-c()
genes<-c()


#for(g in genetic_list_SNV){
for(g in genetic_list[1:66]){
  
  command_txt = paste0("subset(dat_both, ", g, ">0 ) -> sub_clini_dat")
  eval(parse(text=command_txt))
  
  table(sub_clini_dat$response_best3)->ta
  genes <- c(genes,g)
  response_ta <-rbind(response_ta, ta)
  
  
  sub_clini_dat %>% filter(TP53>0) -> sub_clini_dat_TP53
  #sub_clini_dat %>% filter(TP53_multi==TRUE) -> sub_clini_dat_TP53
  
  table(sub_clini_dat_TP53$response_best3)->ta
  response_ta_withTP53 <-rbind(response_ta_withTP53, ta)
  
  
}
total<-apply(response_ta,1,sum)
cbind(genes, response_ta,total) ->response_ta
as.data.frame(response_ta)->response_ta
response_ta$CR<-as.numeric(as.character(response_ta$CR));
response_ta$mCR<-as.numeric(as.character(response_ta$mCR));
response_ta$PR<-as.numeric(as.character(response_ta$PR))
response_ta$`SD+HI`<-as.numeric(as.character(response_ta$`SD+HI`))
response_ta$`SD-HI`<-as.numeric(as.character(response_ta$`SD-HI`))
response_ta$PD<-as.numeric(as.character(response_ta$PD))


total<-apply(response_ta_withTP53,1,sum)
cbind(genes, response_ta_withTP53,total) ->response_ta_withTP53
as.data.frame(response_ta_withTP53)->response_ta_withTP53
response_ta_withTP53$CR<-as.numeric(as.character(response_ta_withTP53$CR))
response_ta_withTP53$mCR<-as.numeric(as.character(response_ta_withTP53$mCR))
response_ta_withTP53$PR<-as.numeric(as.character(response_ta_withTP53$PR))
response_ta_withTP53$`SD+HI`<-as.numeric(as.character(response_ta_withTP53$`SD+HI`))
response_ta_withTP53$`SD-HI`<-as.numeric(as.character(response_ta_withTP53$`SD-HI`))
response_ta_withTP53$PD<-as.numeric(as.character(response_ta_withTP53$PD))



response_ta_withTP53$CR_without <- response_ta$CR - response_ta_withTP53$CR
response_ta_withTP53$mCR_without <- response_ta$mCR - response_ta_withTP53$mCR
response_ta_withTP53$PR_without <- response_ta$PR - response_ta_withTP53$PR
response_ta_withTP53$`SD+HI_without` <- response_ta$`SD+HI` - response_ta_withTP53$`SD+HI`
response_ta_withTP53$`SD-HI_without` <- response_ta$`SD-HI` - response_ta_withTP53$`SD-HI`
response_ta_withTP53$PD_without <- response_ta$PD - response_ta_withTP53$PD



response_ta_withTP53$total <-with(response_ta_withTP53,
                                  CR + mCR + PR + `SD+HI` +`SD-HI`+ PD +
                                    CR_without + mCR_without + PR_without +`SD+HI_without` +`SD-HI_without`+PD_without)




response_ta_withTP53 %>% arrange(desc(as.numeric(as.character(total)))) -> response_ta_withTP53
response_ta_withTP53$genes <-factor(response_ta_withTP53$genes, levels=response_ta_withTP53$genes, order=T)
colnames(response_ta_withTP53) <-c("genes","CR with TP53","mCR with TP53","PR with TP53","SD+HI with TP53",
                                  "SD-HI with TP53","PD with TP53", "total", 
                                  "CR without TP53","mCR without TP53","PR without TP53","SD+HI without TP53",
                                  "SD-HI without TP53","PD without TP53")
response_ta_withTP53 %>% pivot_longer(col = -c(genes,total), names_to="response", values_to="number") ->response_ta2_withTP53
response_ta2_withTP53$response <-factor(response_ta2_withTP53$response, 
                        levels=c("PD with TP53","PD without TP53","SD-HI with TP53","SD-HI without TP53", "SD+HI with TP53",
                                 "SD+HI without TP53", "PR with TP53", "PR without TP53",
                        "mCR with TP53","mCR without TP53", "CR with TP53", "CR without TP53"),order=T)
response_ta2_withTP53$number  <-as.numeric(as.character(response_ta2_withTP53$number))
response_ta2_withTP53$total  <-as.numeric(as.character(response_ta2_withTP53$total))

response_ta2_withTP53 %>% filter(total> 10) -> response_ta2_withTP53

response_col<-c("red2","firebrick3","indianred1","coral2","orange","darkorange1",
                "cornflowerblue","dodgerblue4","forestgreen","darkgreen","wheat4","gray36")


g<-ggplot(data=response_ta2_withTP53, aes(x=reorder(genes, (-1)*total), y=number, fill=response))
g<-g+ geom_bar(stat = "identity")+themePPP
g<-g+xlab("Genes")+ylab("Number")+ggtitle(paste("Response by genes,", cohort_for_test," cohort",sep=""))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                      vjust = 0.5 , size = rel(1.0), face = "italic"))
g<-g + scale_fill_manual(values=rev(response_col))
g<-g+theme(element_line(size=1) )
g



g<-ggplot(data=response_ta2_withTP53, aes(x=reorder(genes, (-1)*total), y=number, fill=response))
g<-g+ geom_bar(stat = "identity", position="fill")+themePPP
g<-g+xlab("Genes")+ylab("Number")+ggtitle(paste("Response by genes,", cohort_for_test," cohort",sep=""))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                      vjust = 0.5 , size = rel(1.0), face = "italic"))
g<-g + scale_fill_manual(values=rev(response_col))
g<-g+theme(element_line(size=1) )
g









###  response by gene, TP53-mutated with highlight
cohort_for_test="Validation"
cohort_for_test="Test"
dat_both <-subset(data_both, cohort2==cohort_for_test)

response_ta<-c()
response_ta_withTP53<-c()
genes<-c()
for(g in genetic_list2[c(1:67)]){
  
  command_txt = paste0("subset(dat_both, ", g, ">0 ) -> sub_clini_dat")
  eval(parse(text=command_txt))
  
  table(sub_clini_dat$response_best3)->ta
  genes <- c(genes,g)
  response_ta <-rbind(response_ta, ta)
  
  sub_clini_dat %>% filter(TP53_multi==TRUE) -> sub_clini_dat_TP53
  table(sub_clini_dat_TP53$response_best3)->ta
  response_ta_withTP53 <-rbind(response_ta_withTP53, ta)
  
  
}
total<-apply(response_ta,1,sum)
cbind(genes, response_ta,total) ->response_ta
as.data.frame(response_ta)->response_ta
response_ta$CR<-as.numeric(as.character(response_ta$CR));
response_ta$mCR<-as.numeric(as.character(response_ta$mCR));
response_ta$PR<-as.numeric(as.character(response_ta$PR))
response_ta$`SD+HI`<-as.numeric(as.character(response_ta$`SD+HI`))
response_ta$`SD-HI`<-as.numeric(as.character(response_ta$`SD-HI`))
response_ta$PD<-as.numeric(as.character(response_ta$PD))
response_ta$total<-as.numeric(as.character(response_ta$total))



cbind(genes, response_ta_withTP53,total) ->response_ta_withTP53
as.data.frame(response_ta_withTP53)->response_ta_withTP53
response_ta_withTP53$CR<-as.numeric(as.character(response_ta_withTP53$CR))
response_ta_withTP53$mCR<-as.numeric(as.character(response_ta_withTP53$mCR))
response_ta_withTP53$PR<-as.numeric(as.character(response_ta_withTP53$PR))
response_ta_withTP53$`SD+HI`<-as.numeric(as.character(response_ta_withTP53$`SD+HI`))
response_ta_withTP53$`SD-HI`<-as.numeric(as.character(response_ta_withTP53$`SD-HI`))
response_ta_withTP53$PD<-as.numeric(as.character(response_ta_withTP53$PD))
response_ta_withTP53$total<-as.numeric(as.character(response_ta_withTP53$total))


response_ta_withoutTP53 <- response_ta
response_ta_withoutTP53[,c(2:7)] <- response_ta[,c(2:7)]- response_ta_withTP53[,c(2:7)]

response_ta_withTP53 %>% mutate(isTP53=1)-> response_ta_withTP53
response_ta_withoutTP53 %>% mutate(isTP53=0)-> response_ta_withoutTP53

response_ta_withTP53 %>% pivot_longer(cols=c(CR, mCR, PR, `SD+HI`,`SD-HI`, PD), names_to = "response", values_to = "N") -> response_ta_withTP53.long
response_ta_withoutTP53 %>% pivot_longer(cols=c(CR, mCR, PR, `SD+HI`,`SD-HI`, PD), names_to = "response", values_to = "N") -> response_ta_withoutTP53.long

rbind(response_ta_withTP53.long, response_ta_withoutTP53.long) ->response_ta.long
response_ta.long$response <- factor(response_ta.long$response, levels=c("CR","mCR","PR","SD+HI","SD-HI","PD"), ordered = T)



response_ta.long %>% filter(total>10) -> response_ta.long

response_col<-c("red2","indianred1","orange", "cornflowerblue","forestgreen","wheat4")

response_col<-c("firebrick3","coral2","darkorange1",
                "dodgerblue4","darkgreen","gray36")

g<-ggplot(data=response_ta.long, aes(x=reorder(genes, (-1)*total), y=N, fill=rev(response), alpha=factor(isTP53)))
g<-g+ geom_bar(stat = "identity")+themePPP
g<-g+xlab("Genes")+ylab("Number")+ggtitle(paste("Response by genes,", cohort_for_test," cohort",sep=""))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                      vjust = 0.5 , size = rel(1.0), face = "italic"))
g<-g + scale_fill_manual(values=rev(response_col))
g<-g+theme(element_line(size=1) )
g<-g+scale_alpha_manual(values = c(1, 0.75))
g


response_ta.long$response_TP53<- paste0(response_ta.long$response,"_",response_ta.long$isTP53)
response_ta.long[,c(1,2,5,6)] %>% pivot_wider(names_from = response_TP53, values_from = N) -> tmp.df
tmp.df %>% arrange(desc(total))->tmp.df


response_col<-c("red2","indianred1","orange", "cornflowerblue","forestgreen","wheat4")
barplot(t(tmp.df[,c(9,3,10,4,11,5,12,6,13,7,14,8)]),names.arg = tmp.df$genes, las=2,
        col=rep(response_col, each=2), density=rep(c(100,80),5), angle=rep(c(0,30),5))

legend("topright", legend=c("CR","mCR","PR","SD+HI","SD-HI","PD"), col=response_col, pch=15)
legend(14,70, legend=c("In multihit TP53 case", "Not in multhit TP53 case"), col=c("black","black"), pch=0)

t(tmp.df[,c(9,3,10,4,11,5,12,6,13,7,14,8)])->m
barplot(sweep(m,2,100/colSums(m),"*"),names.arg = tmp.df$genes, las=2, 
        col=rep(response_col, each=2), density=rep(c(100,80),5), angle=rep(c(0,30),5))



# mutated genes with max VAF in pre-sample
data_both$maxgebe_pre <-c()
data_both$maxVAF_pre <-c()
for (p in 1:nrow(data_both)){
  sampleID = rownames(data_both)[p]
  sub_mut_dat <- subset(mut_dat, Sample == sampleID)
  sub_mut_dat$aVAF <- ifelse(sub_mut_dat$aVAF>1,1,sub_mut_dat$aVAF)
  
  #print(sampleID)
  #print(sub_mut_dat$aVAF)
  
  if(nrow(sub_mut_dat)>0 & sum(!is.na(sub_mut_dat$aVAF))>0){
    n=which.max(sub_mut_dat$aVAF)
    data_both$maxgene_pre[p] <- sub_mut_dat$Gene.refGene[n]
    data_both$maxVAF_pre[p] <- sub_mut_dat$aVAF[n]
  }else if(nrow(sub_mut_dat)==1){
    data_both$maxgene_pre[p] <- sub_mut_dat$Gene.refGene[1]
    data_both$maxVAF_pre[p] <- sub_mut_dat$aVAF[1]
  }else{
    data_both$maxgene_pre[p] <- "NA"
    data_both$maxVAF_pre[p] <- "NA"
  }
}

data_both$maxgene_pre
data_both$maxVAF_pre

table(data_both$TP53>0, data_both$maxgene_pre=="TP53")

ta<-with(subset(data_both, cohort2=="Test"),
  table(maxgene_pre, response_best3)
)
total=apply(ta,1,sum)
genes=rownames(ta)
cbind(genes,ta, total) ->ta
as.data.frame(ta) ->ta
ta<-ta[genes!="NA",]
ta %>% arrange(desc(as.numeric(as.character(total)))) ->  ta
ta$genes <-factor(ta$genes, levels=ta$genes, order=T)
ta %>% pivot_longer(col = -c(genes,total), names_to="response", values_to="number") ->ta2
ta2$response <-factor(ta2$response, levels=c("PD","SD-HI","SD+HI","PR","mCR","CR"),order=T)
ta2$number <-as.numeric(as.character(ta2$number))
ta2$total <-as.numeric(as.character(ta2$total))

response_col<-c("red2","indianred1","orange","cornflowerblue","forestgreen","wheat4")
g<-ggplot(data=subset(ta2, total>5), aes(x=genes, y=number, fill=response))
g<-g+ geom_bar(stat = "identity")+themePPP
g<-g+xlab("Genes")+ylab("Number")+ggtitle("Response by maximum clone in pre-treatment samples,\n discovery cohort (n>5)")
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                      vjust = 0.5 , size = rel(1.0), face = "italic"))
g<-g + scale_fill_manual(values=rev(response_col))
g


response_col<-c("red2","indianred1","orange","cornflowerblue","forestgreen","wheat4")
g<-ggplot(data=subset(ta2, total>5), aes(x=genes, y=number, fill=response))
g<-g+ geom_bar(stat = "identity", position="fill")+themePPP
g<-g+xlab("Genes")+ylab("Number")+ggtitle("Response by maximum clone in pre-treatment samples,\n discovery cohort (n>5)")
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                      vjust = 0.5 , size = rel(1.0), face = "italic"))
g<-g + scale_fill_manual(values=rev(response_col))
g



