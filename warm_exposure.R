library(readxl); library(ggplot2)

x=as.data.frame(read_excel("Warm_exposure2.xlsx",col_types=c("text","date","text","text","text","text","date","text","text")))
x <- x[-which(x$Final_condition=="dead"),]

x$Sex <- ifelse(grepl("C2/7",x$ID),"Male strain","Female strain")
x$Lineage <- sapply(lapply(strsplit(x$ID,"-"), "[", 1:2),paste,collapse="-")

x$repr.mode.binary <- ifelse(x$Final_condition=="asexual", 0, 1)

x$ExpGroup <- factor(x$ExpGroup, levels=c("control", "warmexposed_1week", "warmexposed_4weeks"))

tiff("warm_exposure.tif",width=6, height=6, units="in",res=600,compression="lzw")
ggplot(x, aes(fill=Final_condition, x=ExpGroup)) + geom_bar(color="black") + facet_wrap(~Sex)+
    theme_classic() +
    scale_fill_manual(labels=c("Sexual","Asexual"),breaks=c("sexual","asexual"),values=c("white", "grey"))+
    labs(fill="Reproductive mode")+
    theme(strip.background = element_blank(), strip.text.x=element_text(size=20, face="bold", colour="red"),
          text=element_text(size=15),legend.position="top")+
    ylab("No. polyps")+xlab("Warm exposure") + scale_x_discrete(labels=c("None", "1 week","4 weeks"))
dev.off()

fisher.test(table(x$Final_condition[x$Sex=="Male strain"], x$ExpGroup[x$Sex=="Male strain"])) # P<0.001
fisher.test(table(x$Final_condition[x$Sex=="Female strain"], x$ExpGroup[x$Sex=="Female strain"])) # P<0.001

length(as.Date(x$Final_date[x$Sex=="Female strain" & x$Final_condition=="sexual" & x$ExpGroup=="warmexposed_1week"]) - as.Date("2019-07-15"))
median(as.Date(x$Final_date[x$Sex=="Female strain" & x$Final_condition=="sexual" & x$ExpGroup=="warmexposed_1week"]) - as.Date("2019-07-15"))

median(as.Date(x$Final_date[x$Sex=="Male strain" & x$Final_condition=="sexual" & x$ExpGroup=="warmexposed_4weeks"]) - as.Date("2019-08-06"))
length(as.Date(x$Final_date[x$Sex=="Male strain" & x$Final_condition=="sexual" & x$ExpGroup=="warmexposed_4weeks"]) - as.Date("2019-08-06"))

median(as.Date(x$Final_date[x$Sex=="Female strain" & x$Final_condition=="sexual" & x$ExpGroup=="warmexposed_4weeks"]) - as.Date("2019-08-06"))
length(as.Date(x$Final_date[x$Sex=="Female strain" & x$Final_condition=="sexual" & x$ExpGroup=="warmexposed_4weeks"]) - as.Date("2019-08-06"))
