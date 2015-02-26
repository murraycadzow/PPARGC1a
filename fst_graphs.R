#setwd to Analysis/Dec2014/Fst
library(ggplot2)
library(scales)
get_fst = function(fst, pop1, pop2){
  a= read.table(fst, header=TRUE)
  b=a[a[,3]!= "NaN",]
  b$POP1= as.factor(pop1)
  b$POP2= as.factor(pop2)

  return(b)
 }

bp_stat=function(fst){
  a=data.frame(
    snpwindow=fst[fst["BIN_START"] == 23815662, "MEAN_FST"],
    lower = quantile(fst$MEAN_FST, 0.025),
    mean = mean(fst$MEAN_FST),
    upper= quantile(fst$MEAN_FST, 0.975)
  )
  return(a)
}

window_stat=function(fst){
  a=data.frame(
    snpwindow=fst[fst["BIN_START"] < 23815662 & fst[,"BIN_END"] > 23815662, "MEAN_FST"],
    lower = quantile(fst$MEAN_FST, 0.025),
    mean = mean(fst$MEAN_FST),
    upper= quantile(fst$MEAN_FST, 0.975)
  )
  return(a)
}

Fst = get_fst(fst="axiom/AXIOM_CEU_5000000.windowed.weir.fst", pop1="NZ Maori", pop2="CEU")
Fst = rbind(Fst, get_fst(fst="axiom/AXIOM_CHB_5000000.windowed.weir.fst", pop1="NZ Maori", pop2="CHB") )
Fst = rbind(Fst, get_fst(fst="axiom/AXIOM_CHS_5000000.windowed.weir.fst", pop1="NZ Maori", pop2="CHS") )
Fst = rbind(Fst, get_fst(fst="axiom/AXIOM_GBR_5000000.windowed.weir.fst", pop1="NZ Maori", pop2="GBR") )
Fst = rbind(Fst, get_fst(fst="axiom/AXIOM_YRI_5000000.windowed.weir.fst", pop1="NZ Maori", pop2="YRI"))
Fst = rbind(Fst, get_fst(fst="omni/OMNI_CEU_5000000.windowed.weir.fst", pop1="Samoan", pop2="CEU"))
Fst = rbind(Fst, get_fst(fst="omni/OMNI_CHB_5000000.windowed.weir.fst", pop1="Samoan", pop2="CHB") )
Fst = rbind(Fst, get_fst(fst="omni/OMNI_CHS_5000000.windowed.weir.fst", pop1="Samoan", pop2="CHS") )
Fst = rbind(Fst, get_fst(fst="omni/OMNI_GBR_5000000.windowed.weir.fst", pop1="Samoan", pop2="GBR") )
Fst = rbind(Fst, get_fst(fst="omni/OMNI_YRI_5000000.windowed.weir.fst", pop1="Samoan", pop2="YRI"))

Fst[Fst$MEAN_FST < 0, "MEAN_FST" ] = 0
library(dplyr)
s = ddply(Fst, POP2 ~ POP1 ,summarise,q1=quantile(MEAN_FST,0.025), q2=quantile(MEAN_FST,0.975), m= mean(MEAN_FST))

#ppargc1a_fst = subset(Fst, BIN_START > 23756664 & BIN_START < 23905712)
p = ggplot(data = Fst, aes(x=BIN_START, y=MEAN_FST), ylim=c(min(Fst$MEAN_FST -0.5), max(Fst$MEAN_FST + 0.5)))  + geom_point(shape = 16, colour= alpha("black", 1/5)) + facet_grid(POP1~POP2) + geom_hline(aes(yintercept= q1, colour ="quantile"), data=s,)  + geom_hline(aes(yintercept= q2, colour ="quantile"), data=s)  +geom_hline(aes(yintercept= m, colour="mean"), data=s)+ scale_colour_manual("",breaks=c("Threshold","Mean"),values=c("blue", "purple")) +  xlab("Position (Mbp)") +scale_x_continuous(breaks=c(0e6,50e6,100e6,150e6), labels=c("0","50", "100", "150"))+ ylab("Fst") + ggtitle("Fst (5 Mbp Windows) Chromosome 4") + geom_vline(xintercept = 23815662, color="red", linetype="longdash")
png("fst_5mbp.png")
p
dev.off()



Fst = get_fst(fst="axiom/AXIOM_CEU_1.windowed.weir.fst", pop1="NZ Maori", pop2="CEU")
Fst = rbind(Fst, get_fst(fst="axiom/AXIOM_CHB_1.windowed.weir.fst", pop1="NZ Maori", pop2="CHB") )
Fst = rbind(Fst, get_fst(fst="axiom/AXIOM_CHS_1.windowed.weir.fst", pop1="NZ Maori", pop2="CHS") )
Fst = rbind(Fst, get_fst(fst="axiom/AXIOM_GBR_1.windowed.weir.fst", pop1="NZ Maori", pop2="GBR") )
Fst = rbind(Fst, get_fst(fst="axiom/AXIOM_YRI_1.windowed.weir.fst", pop1="NZ Maori", pop2="YRI"))
Fst = rbind(Fst, get_fst(fst="omni/OMNI_CEU_1.windowed.weir.fst", pop1="Samoan", pop2="CEU"))
Fst = rbind(Fst, get_fst(fst="omni/OMNI_CHB_1.windowed.weir.fst", pop1="Samoan", pop2="CHB") )
Fst = rbind(Fst, get_fst(fst="omni/OMNI_CHS_1.windowed.weir.fst", pop1="Samoan", pop2="CHS") )
Fst = rbind(Fst, get_fst(fst="omni/OMNI_GBR_1.windowed.weir.fst", pop1="Samoan", pop2="GBR") )
Fst = rbind(Fst, get_fst(fst="omni/OMNI_YRI_1.windowed.weir.fst", pop1="Samoan", pop2="YRI"))
Fst[Fst$MEAN_FST < 0, "MEAN_FST" ] = 0

s = ddply(Fst, POP2 ~ POP1 ,summarise,q1=quantile(MEAN_FST,0.025), q2=quantile(MEAN_FST,0.975), m= mean(MEAN_FST))
rs8192678_fst = subset(Fst, BIN_START == 23815662)
p = ggplot(data = Fst, aes(x=BIN_START, y=MEAN_FST), ylim=c(min(Fst$MEAN_FST -0.5), max(Fst$MEAN_FST + 0.5)))  + geom_point(shape = 16, colour= alpha("black", 1/5)) +  geom_point(data=rs8192678_fst, aes(y = MEAN_FST, x=BIN_START,colour="rs8192678")) + facet_grid(POP1~POP2) + geom_hline(aes(yintercept= q1, colour ="quantile"), data=s,)  + geom_hline(aes(yintercept= q2, colour ="quantile"), data=s)  +geom_hline(aes(yintercept= m, colour="mean"), data=s)+ scale_colour_manual("",breaks=c("Threshold", "rs8192678", "Mean"),values=c("blue",  "purple","red")) +  xlab("Position (Mbp)") +scale_x_continuous(breaks=c(0e6,50e6,100e6,150e6), labels=c("0","50", "100", "150"))+ ylab("Fst") + ggtitle("Fst (1 bp Windows) Chromosome 4") + geom_vline(xintercept = 23815662, color="red", linetype="longdash")
png("fst_1bp.png")
p
dev.off()