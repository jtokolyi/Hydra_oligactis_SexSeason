library(ggplot2); library(ggpubr); library(gridExtra); library(gtable); library(grid); library(cowplot);
library(nlme); library(lme4); library(EMAtools)

x <- read.csv("M28seasonal_full_dataset2.csv",sep="\t")

####################################################################################################
########################## IMPORTANT NOTE ##########################################################
## Sex data for two strains are added here to the data table
## These strains contained sex-changed individuals (i.e. both males and females in the same strain),
## therefore, they were not included in the original dataset uploaded to Figshare (https://doi.org/10.6084/m9.figshare.12727673.v1)
## However, for the purpose of this analysis we have no reason to assume that they should bias our results,
## therefore we included them.

x$Sex[x$Strain=="M28/2019autumn/1/5"] <- c("F","F","F","F","F","M","M","M","M")
x$SexStart[x$Strain=="M28/2019autumn/1/5"] <- c(27,27,30,27,27,30,27,51,30)

x$Sex[x$Strain=="M28/2019autumn/9/1"] <- c("F","F","M","F","F","F","F","F","M")
x$SexStart[x$Strain=="M28/2019autumn/9/1"] <- c(27,27,20,27,27,27,27,27,20)
#####################################################################################################

x$Time <- factor(as.character(x$Time), levels=c("2018spring","2018autumn","2019spring","2019autumn"))
x$Season <- factor(as.character(x$Season), levels=c("Spring","Autumn"))
x$Year <- factor(as.character(x$Year), levels=c("2018","2019"))

x$repr.mode <- ifelse(x$Sex=="ASEX","asex","sex")
x$repr.mode.binomial <- ifelse(x$Sex=="ASEX",0,1)
x$wild <- ifelse(grepl("A/1$",x$PolypID), "yes", "no")
x$Strain <- as.character(x$Strain)
x$logSexStart <- log(x$SexStart)

m <- x[x$Sex=="M",]
f <- x[x$Sex=="F",]
mf <- rbind(m,f)

size.mod <- lme(Body_size ~ Season, random= ~1|Strain, data=x[which(!is.na(x$Body_size)),])
summary(size.mod)

##################### Analyze proportion of sexual individuals ###############################
asm <- glmer(repr.mode.binomial~Season+PolypAge+Year+wild+(1|Strain), data=x,family="binomial")
summary(asm)
exp(asm@beta) ## Odds ratio
asm.red1 <- glmer(repr.mode.binomial~Season+PolypAge+Year+(1|Strain), data=x,family="binomial")
summary(asm.red1)
asm.red2 <- glmer(repr.mode.binomial~Season+PolypAge+(1|Strain), data=x,family="binomial")
summary(asm.red2)
exp(asm.red2@beta)

##################### Analyze sex start times ################################################
mf.comb <- lme(logSexStart~Sex*Season+Sex*PolypAge+Sex*Year+Sex*wild,random=~1|Strain, data=mf)
summary(mf.comb)
lme.dscore(mf.comb,type="nlme") # Cohen's d
mf.comb.red1 <- lme(log(SexStart)~Sex*Season+Sex*PolypAge+Year+Sex*wild,random=~1|Strain, data=mf)
summary(mf.comb.red1)
mf.comb.red2 <- lme(log(SexStart)~Sex+Season+Sex*PolypAge+Year+Sex*wild,random=~1|Strain, data=mf)
summary(mf.comb.red2)
mf.comb.red3 <- lme(log(SexStart)~Sex+Season+PolypAge+Year+Sex*wild,random=~1|Strain, data=mf)
summary(mf.comb.red3)
mf.comb.red4 <- lme(log(SexStart)~Sex+Season+PolypAge+Year+wild,random=~1|Strain, data=mf)
summary(mf.comb.red4)
mf.comb.red5 <- lme(log(SexStart)~Sex+Season+PolypAge+Year,random=~1|Strain, data=mf)
summary(mf.comb.red5)
mf.comb.red6 <- lme(log(SexStart)~Sex+Season+PolypAge,random=~1|Strain, data=mf)
summary(mf.comb.red6)
lme.dscore(mf.comb.red6,type="nlme")

##############################################################################################
################### Repeat with wild-collected individuals removed ###########################

asm <- glmer(repr.mode.binomial~Season+PolypAge+Year+(1|Strain), data=x[x$wild=="no",],family="binomial")
summary(asm)
exp(asm@beta)
asm.red1 <- glmer(repr.mode.binomial~Season+PolypAge+(1|Strain), data=x[x$wild=="no",],family="binomial")
summary(asm.red1)
exp(asm.red1@beta)

mf.comb <- lme(log(SexStart)~Sex*Season+Sex*PolypAge+Sex*Year,random=~1|Strain, data=mf[mf$wild=="no",])
summary(mf.comb)
lme.dscore(mf.comb,type="nlme") # Cohen's d
mf.comb.red1 <- lme(log(SexStart)~Sex*Season+Sex*PolypAge+Year,random=~1|Strain, data=mf[mf$wild=="no",])
summary(mf.comb.red1)
mf.comb.red2 <- lme(log(SexStart)~Sex+Season+Sex*PolypAge+Year,random=~1|Strain, data=mf[mf$wild=="no",])
summary(mf.comb.red2)
mf.comb.red3 <- lme(log(SexStart)~Sex+Season+PolypAge+Year,random=~1|Strain, data=mf[mf$wild=="no",])
summary(mf.comb.red3)
mf.comb.red4 <- lme(log(SexStart)~Sex+Season+PolypAge,random=~1|Strain, data=mf[mf$wild=="no",])
summary(mf.comb.red4)
lme.dscore(mf.comb.red4,type="nlme") # Cohen's d

full.output <- rbind(round(summary(asm)[[10]],3), round(summary(mm)[[20]][,(-3)], 3), round(summary(fm)[[20]][,(-3)], 3))
write.table(cbind(row.names(full.output),
                  paste(full.output[,1],"(",full.output[,2],")"),
                  paste(full.output[,3],"(",full.output[,4],")")),
            file="table1full.csv",col.names=F, row.names=F, quote=F, sep="\t")

red.output <- rbind(round(summary(asm.red2)[[10]],3), round(summary(mm.red2)[[20]][,(-3)], 3), round(summary(fm.red2)[[20]][,(-3)], 3))
write.table(cbind(row.names(red.output),
                  paste(red.output[,1],"(",red.output[,2],")"),
                  paste(red.output[,3],"(",red.output[,4],")")),
            file="table1red.csv",col.names=F, row.names=F, quote=F, sep="\t")

#################### Plot sex proportion

ggplot(x, aes(fill=repr.mode, x=Season)) + geom_bar(color="black") + facet_wrap(~Year)+
    theme_classic() +
    scale_fill_manual(labels=c("Sexual","Asexual"),breaks=c("sex","asex"),values=c("white", "grey"))+
    labs(fill="Reproductive mode")+
    theme(strip.background = element_blank(), strip.text.x=element_text(size=20, face="bold", colour="red"),
          text=element_text(size=15),legend.position="top")+
    ylab("No. polyps")
ggsave(filename="seasonal_fig1.tiff",width=6, height=6, units="in",dpi=600)

mp1 <- ggplot(m, aes(x=SexStart,fill=Season)) + geom_histogram(aes(y=stat(count/sum(count))), color="black",alpha=1,binwidth=7) +
    facet_wrap(~Year) + theme_classic() + 
    scale_fill_manual(labels=c("Spring","Autumn"),breaks=c("Spring","Autumn"),values=c("#7FC97F", "#BEAED4"))+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
        theme(strip.background = element_blank(), strip.text.x=element_blank(),plot.margin=unit(c(0,0,0,0), "cm"),
              text=element_text(size=12), legend.justification=c("top"),
              axis.title.y = element_text(margin=margin(t=0,b=0,l=0,r=3)))+
        xlab("Appearance of testes (days after cooling)")+ylab("Prop. of polyps")+
        scale_x_continuous(breaks=seq(0,150,25))+xlim(0,150)+ylim(0,0.25)
mp2 <- ggplot(m, aes(y=SexStart,fill=Season))+geom_boxplot()+coord_flip()+
    facet_wrap(~Year) +
    scale_fill_manual(labels=c("Spring","Autumn"),breaks=c("Spring","Autumn"),values=c("#7FC97F", "#BEAED4"))+
    theme_classic() +
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1),text=element_text(size=12)) +
        theme(strip.background = element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank(),
              strip.text.x=element_text(size=12, face="bold", colour="red"),
              plot.margin = unit(c(0,0,-0.5,0) , units = "cm" ),
              legend.position="none",
              plot.tag=element_text(face="bold"),plot.tag.position=c(-0.05,0.75))+
        ylab("")+labs(tag="a") + scale_x_continuous(breaks=seq(0,150,25)) + ylim(0,150)

fp1 <- ggplot(f, aes(x=SexStart,fill=Season)) + geom_histogram(aes(y=stat(count/sum(count))),color="black",alpha=1,binwidth=7) +
    facet_wrap(~Year) + theme_classic() + 
    scale_fill_manual(labels=c("Spring","Autumn"),breaks=c("Spring","Autumn"),values=c("#7FC97F", "#BEAED4"))+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
        theme(strip.background = element_blank(), strip.text.x=element_blank(),plot.margin=unit(c(0,0,0,0), "cm"),
              text=element_text(size=12), legend.justification=c("top"),
              axis.title.y = element_text(margin=margin(t=0,b=0,l=0,r=10)))+
        xlab("Appearance of eggs (days after cooling)")+ylab("Prop. of polyps")+
        scale_x_continuous(breaks=seq(0,150,25)) + xlim(0,150) + ylim(0,0.25)
fp2 <- ggplot(f, aes(y=SexStart,fill=Season))+geom_boxplot()+coord_flip()+
    facet_wrap(~Year) +
    scale_fill_manual(labels=c("Spring","Autumn"),breaks=c("Spring","Autumn"),values=c("#7FC97F", "#BEAED4"))+
    theme_classic() +
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1),text=element_text(size=12)) +
        theme(strip.background = element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank(),
              plot.margin = unit(c(0,0,-0.5,0) , units = "cm" ),
              strip.text.x=element_text(size=12, face="bold", colour="red"),
              plot.tag=element_text(face="bold"),plot.tag.position=c(-0.05,0.75),
              legend.position="none")+
        ylab("")+labs(tag="b") + scale_x_continuous(breaks=seq(0,150,25)) + ylim(0,150)
tiff(file="seasonal_fig2_mod.tiff",width=6, height=6, units="in", res=600, compression="lzw")
plot_grid(mp2,mp1,fp2,fp1,rel_heights = c(0.9, 2, 0.9, 2), align = "v", ncol = 1, nrow = 4, axis="lr")
dev.off()

