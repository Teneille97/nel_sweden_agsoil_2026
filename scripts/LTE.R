#tnel edits
##LTE

library(here)

Long=read.csv2(here("csv_files", "LTEnitrogen2.csv"))
head(Long)
unique(Long$LTE)

M2=subset(Long,Long$LTE=="M2-1957")
M4=subset(Long,Long$LTE=="M4-1957")
M6=subset(Long,Long$LTE=="M6-1957")
R94=subset(Long,Long$LTE=="R94-1966")

M2_Soil=subset(M2,M2$Temperature...C.=="Soil")
M2_325=subset(M2,M2$Temperature...C.=="325")
M2_400=subset(M2,M2$Temperature...C.=="400")

M4_Soil=subset(M4,M4$Temperature...C.=="Soil")
M4_325=subset(M4,M4$Temperature...C.=="325")
M4_400=subset(M4,M4$Temperature...C.=="400")

M6_Soil=subset(M6,M6$Temperature...C.=="Soil")
M6_325=subset(M6,M6$Temperature...C.=="325")
M6_400=subset(M6,M6$Temperature...C.=="400")

R94_Soil=subset(R94,R94$Temperature...C.=="Soil")
R94_175=subset(R94,R94$Temperature...C.=="175")
R94_250=subset(R94,R94$Temperature...C.=="250")
R94_325=subset(R94,R94$Temperature...C.=="325")
R94_400=subset(R94,R94$Temperature...C.=="400")
R94_445=subset(R94,R94$Temperature...C.=="445")
R94_500=subset(R94,R94$Temperature...C.=="500")


##Plot C and N====
pdf(file="LTE longCN.pdf", encoding='WinAnsi.enc', width=12, height=12.5)
par(mfrow=c(4,3),mar=c(4.14,4.6,0.7,1.2))
plot(M2_Soil$Cgkg~M2_Soil$Year,pch=19,col="brown",cex.lab=1.5,ylim=c(0,30),xlab="",ylab=expression(paste("Total organic carbon (TOC) [g ", kg^-1,"]")))
points(M2_325$Cgkg~M2_325$Year,pch=19,col="blue")
points(M2_400$Cgkg~M2_400$Year,pch=19,col="red")
legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","red"),pch = 19,bty = "n",cex=1.2)
legend("topleft", "(A) Site M2",bty="n", cex=1.5)

plot(M2_Soil$TONgkg~M2_Soil$Year,pch=19,col="brown",cex.lab=1.5,ylim=c(0,2.6),xlab="",ylab=expression(paste("Total organic nitrogen (TON) [g ", kg^-1,"]")))
points(M2_325$TONgkg~M2_325$Year,pch=19,col="blue")
points(M2_400$TONgkg~M2_400$Year,pch=19,col="red")
legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","red"),pch = 19,bty = "n",cex=1.2)
legend("topleft", "(B) Site M2",bty="n", cex=1.5)

plot(M2_Soil$MolarCN~M2_Soil$Year,pch=19,col="brown",cex.lab=1.5,ylim=c(0,20),xlab="",ylab=expression(paste("Molar TOC:TON ratio")))
points(M2_325$MolarCN~M2_325$Year,pch=19,col="blue")
points(M2_400$MolarCN~M2_400$Year,pch=19,col="red")
legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","red"),pch = 19,bty = "n",cex=1.2)
legend("topleft", "(C) Site M2",bty="n", cex=1.5)

plot(M4_Soil$Cgkg~M4_Soil$Year,pch=19,col="brown",cex.lab=1.5,ylim=c(0,30),xlab="",ylab=expression(paste("Total organic carbon (TOC) [g ", kg^-1,"]")))
points(M4_325$Cgkg~M4_325$Year,pch=19,col="blue")
points(M4_400$Cgkg~M4_400$Year,pch=19,col="red")
legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","red"),pch = 19,bty = "n",cex=1.2)
legend("topleft", "(D) Site M4",bty="n", cex=1.5)

plot(M4_Soil$TONgkg~M4_Soil$Year,pch=19,col="brown",cex.lab=1.5,ylim=c(0,2.6),xlab="",ylab=expression(paste("Total organic nitrogen (TON) [g ", kg^-1,"]")))
points(M4_325$TONgkg~M4_325$Year,pch=19,col="blue")
points(M4_400$TONgkg~M4_400$Year,pch=19,col="red")
legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","red"),pch = 19,bty = "n",cex=1.2)
legend("topleft", "(E) Site M4",bty="n", cex=1.5)

plot(M4_Soil$MolarCN~M4_Soil$Year,pch=19,col="brown",cex.lab=1.5,ylim=c(0,20),xlab="",ylab=expression(paste("Molar TOC:TON ratio")))
points(M4_325$MolarCN~M4_325$Year,pch=19,col="blue")
points(M4_400$MolarCN~M4_400$Year,pch=19,col="red")
legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","red"),pch = 19,bty = "n",cex=1.2)
legend("topleft", "(F) Site M4",bty="n", cex=1.5)

plot(M6_Soil$Cgkg~M6_Soil$Year,pch=19,col="brown",cex.lab=1.5,ylim=c(0,30),xlab="",ylab=expression(paste("Total organic carbon (TOC) [g ", kg^-1,"]")))
points(M6_325$Cgkg~M6_325$Year,pch=19,col="blue")
points(M6_400$Cgkg~M6_400$Year,pch=19,col="red")
legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","red"),pch = 19,bty = "n",cex=1.2)
legend("topleft", "(G) Site M6",bty="n", cex=1.5)

plot(M6_Soil$TONgkg~M6_Soil$Year,pch=19,col="brown",cex.lab=1.5,ylim=c(0,2.6),xlab="",ylab=expression(paste("Total organic nitrogen (TON) [g ", kg^-1,"]")))
points(M6_325$TONgkg~M6_325$Year,pch=19,col="blue")
points(M6_400$TONgkg~M6_400$Year,pch=19,col="red")
legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","red"),pch = 19,bty = "n",cex=1.2)
legend("topleft", "(H) Site M6",bty="n", cex=1.5)

plot(M6_Soil$MolarCN~M6_Soil$Year,pch=19,col="brown",cex.lab=1.5,ylim=c(0,20),xlab="",ylab=expression(paste("Molar TOC:TON ratio")))
points(M6_325$MolarCN~M6_325$Year,pch=19,col="blue")
points(M6_400$MolarCN~M6_400$Year,pch=19,col="red")
legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","red"),pch = 19,bty = "n",cex=1.2)
legend("topleft", "(I) Site M6",bty="n", cex=1.5)

plot(R94_Soil$Cgkg~R94_Soil$Year,pch=19,col="brown",cex.lab=1.5,ylim=c(0,30),xlab="Calendar year",ylab=expression(paste("Total organic carbon (TOC) [g ", kg^-1,"]")))
#points(R94_175$Cgkg~R94_175$Year,pch=19,col="blue")
#points(R94_250$Cgkg~R94_250$Year,pch=19,col="green")
points(R94_325$Cgkg~R94_325$Year,pch=19,col="blue")
points(R94_400$Cgkg~R94_400$Year,pch=19,col="red")
#points(R94_445$Cgkg~R94_445$Year,pch=19,col="orange")
#points(R94_500$Cgkg~R94_500$Year,pch=19,col="pink")
legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","red"),pch = 19,bty = "n",cex=1.2)
legend("topleft", "(J) Site R94",bty="n", cex=1.5)

plot(R94_Soil$TONgkg~R94_Soil$Year,pch=19,col="brown",cex.lab=1.5,ylim=c(0,2.6),xlab="Calendar year",ylab=expression(paste("Total organic nitrogen (TON) [g ", kg^-1,"]")))
#points(R94_175$TONgkg~R94_175$Year,pch=19,col="blue")
#points(R94_250$TONgkg~R94_250$Year,pch=19,col="green")
points(R94_325$TONgkg~R94_325$Year,pch=19,col="blue")
points(R94_400$TONgkg~R94_400$Year,pch=19,col="red")
#points(R94_445$TONgkg~R94_445$Year,pch=19,col="orange")
#points(R94_500$TONgkg~R94_500$Year,pch=19,col="pink")
legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","red"),pch = 19,bty = "n",cex=1.2)
legend("topleft", "(K) Site R94",bty="n", cex=1.5)

plot(R94_Soil$MolarCN~R94_Soil$Year,pch=19,col="brown",cex.lab=1.5,ylim=c(0,20),xlab="Calendar year",ylab=expression(paste("Molar TOC:TON ratio")))
#points(R94_175$MolarCN~R94_175$Year,pch=19,col="blue")
#points(R94_250$MolarCN~R94_250$Year,pch=19,col="green")
points(R94_325$MolarCN~R94_325$Year,pch=19,col="blue")
points(R94_400$MolarCN~R94_400$Year,pch=19,col="red")
#points(R94_445$MolarCN~R94_445$Year,pch=19,col="orange")
#points(R94_500$MolarCN~R94_500$Year,pch=19,col="pink")
legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","red"),pch = 19,bty = "n",cex=1.2)
legend("topleft", "(L) Site R94",bty="n", cex=1.5)
dev.off()


##14C====
pdf(file="//storage.slu.se/Home$/mesp0003/Desktop/LTE long delta14C.pdf", encoding='WinAnsi.enc', width=5, height=12.5)
par(mfrow=c(4,1),mar=c(4.14,4.6,0.7,1.2))
plot(M2_Soil$d14C~M2_Soil$Year,pch=19,col="brown",xlab="", cex.lab=1.2, ylim=c(-650,150),ylab=expression(paste(Delta^14, "C (\u2030)")))
points(M2_325$d14C~M2_325$Year,pch=19,col="blue")
points(M2_400$d14C~M2_400$Year,pch=19,col="orange")
legend("bottomright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","orange"),pch = 19,bty = "n",cex=1)
legend("topleft", "(A) Site M2",bty="n", cex=1.3)


plot(M4_Soil$d14C~M4_Soil$Year,pch=19,col="brown",xlab="",cex.lab=1.2, ylim=c(-650,150), ylab=expression(paste(Delta^14, "C (\u2030)")))
points(M4_325$d14C~M4_325$Year,pch=19,col="blue")
points(M4_400$d14C~M4_400$Year,pch=19,col="orange")
#legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","orange"),pch = 19,bty = "n",cex=1)
legend("topleft", "(B) Site M4",bty="n", cex=1.3)


plot(M6_Soil$d14C~M6_Soil$Year,pch=19,col="brown",xlab="",cex.lab=1.2, ylim=c(-650,150),ylab=expression(paste(Delta^14, "C (\u2030)")))
points(M6_325$d14C~M6_325$Year,pch=19,col="blue")
points(M6_400$d14C~M6_400$Year,pch=19,col="orange")
#legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","orange"),pch = 19,bty = "n",cex=1)
legend("topleft", "(C) Site M6",bty="n", cex=1.3)

plot(R94_Soil$d14C~R94_Soil$Year,pch=19,col="brown",xlab="",cex.lab=1.2, ylim=c(-650,150),ylab=expression(paste(Delta^14, "C (\u2030)")))
#points(R94_175$d14C~R94_175$Year,pch=19,col="blue")
#points(R94_250$d14C~R94_250$Year,pch=19,col="green")
points(R94_325$d14C~R94_325$Year,pch=19,col="blue")
points(R94_400$d14C~R94_400$Year,pch=19,col="orange")
#points(R94_445$d14C~R94_445$Year,pch=19,col="orange")
#points(R94_500$d14C~R94_500$Year,pch=19,col="pink")
#legend("topright", legend = c("Soil", "325°C", "400°C") , col = c("brown","blue","orange"),pch = 19,bty = "n",cex=1)
legend("topleft", "(D) Site R96",bty="n", cex=1.3)
dev.off()
