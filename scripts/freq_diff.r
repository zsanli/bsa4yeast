library(lattice)
args = commandArgs(trailingOnly=TRUE)
snps<-read.table(args[1],sep=" ",header=F)
colnames(snps)<-c('chr','position','frequency')
snps$chr<-as.roman(unlist(lapply(as.character(snps$chr),as.roman)))
a=as.roman(as.character(snps$chr))
b<-as.numeric(a)
fb = factor(b)
rb = factor(b,labels=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"))
snps$chr<-rb
pdf(args[2],width=25,height=25)
xyplot(frequency ~position|chr,pch='.',data=snps,type='p',
       panel = function(x, y) {
           panel.xyplot(x, y,pch='.',type='p')
			panel.loess(x, y,col='red',lty=1)
           panel.abline(a=0,col='red',lty='longdash')
       })
dev.off()