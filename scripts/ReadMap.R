args<-commandArgs(TRUE)
file<-args[2]
out<-args[3]
rt<-read.table(file)  
pdf(out,width=15,height=10)
par(bg="white");
m=levels(as.factor(rt[,1]));
m=m[grep("[_|M]",m,invert=T)];
m
m=gsub("chr","",m)
m=m[order(as.numeric(m))]
m
rt[,1]=gsub("chr","",rt[,1])
rt<-cbind(rt,apply(rt,1,function(x) as.numeric(x[4])/as.numeric(x[6]) ))
par(oma=c(3,2,1,3),mfrow=c(length(m),1),mar=c(0,5,0,0),mgp=c(3,1,0));
for(i in m){
plot(rt[rt[,1]==i,2],rt[rt[,1]==i,7],type="h",col="royalblue",ylab="",xlim=c(1,max(rt[,2])),ylim=c(0,1),axes=F,yaxt="n",xaxt="n",cex.lab=1.2,xlab="",lwd=1,las=1,xaxs="i");
par(srt=0);
#text(-1/30*(max(rt[,2])-min(rt[,2])),(max(rt[rt[,1]==i,7])-min(rt[rt[,1]==i,7])),paste("chr",i,sep=""),xpd=T,font=1,cex=1.5);
text(-1/30*(max(rt[,2])-min(rt[,2])),0.75,paste("chr",i,sep=""),xpd=T,font=1,cex=1.5);
grid(col="grey");
abline(h=0,col="brown");
par(new=T);
plot(rt[rt[,1]==i,2],rt[rt[,1]==i,7],type="h",col="royalblue2",ylab="",xlim=c(1,max(rt[,2])),ylim=c(0,1),axes=F,yaxt="n",xaxt="n",cex.lab=1.2,xlab="",lwd=1,las=1,xaxs="i");
par(new=T);
plot(rt[rt[,1]==i,2],rt[rt[,1]==i,8],type="l",col="red",ylab="",xlim=c(1,max(rt[,2])),ylim=c(0,1),axes=F,yaxt="n",xaxt="n",cex.lab=1.2,xlab="",lwd=1,las=1,xaxs="i");};
axis(1,xaxp=c(1,max(rt[,2]),20))

