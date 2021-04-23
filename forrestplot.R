# box plot
HRa_quantile
ev_table
HR_a

dev.new(width=6,height=3,pointsize=9)
par(mar=c(3,20,5.5,9)+.1, mgp=c(1.75,0.5,0), tcl=-0.4)
xl <- c(0.3,2)
yl <- c(0.5,3)
#
plot(xl,yl,axes=F,type="n",xaxs="i",xlab="",ylab="",log="x")
segments(1,par('usr')[3],1,11)
segments(HR_a[1],par('usr')[3],HR_a[1],11, col="#374E55FF", lty=3)
segments(HRa_quantile[2,],3:1,HRa_quantile[3,],3:1)
points(HRa_quantile[1,],3:1, pch=15, col="#00A1D5FF", cex=1/(log(HRa_quantile[3,])-log(HRa_quantile[1,])))
axis(1,at=c(0.3, 0.6,1, 1.4, 2), labels=c(0.3, 0.6,1, 1.4, 2), cex.axis=.8, lwd=0.75, lwd.tick=0.75)
axis(1,at=c(seq(0.06,0.09,by=0.01),seq(0.3,0.9,by=0.1),seq(1.1,1.9,by=0.1),3,4), labels=NA, tcl=-0.25, lwd=0.75, lwd.tick=0.75)
mtext("Adjusted Hazard Ratio", side=1, line=1.75, cex=1)
#
mtext("pH at enrollment",side=2, line=19, at=3.25, font=2, las=1, adj=0)
mtext(c("< 7.22","7.22-7.27", "7.27-7.30"),side=2, line=17, at=3:1, las=1, adj=0)
#
mtext("Early", side=2, line=9, at=3.5, font=2, las=1, adj=0.5)
mtext("Delayed", side=2, line=3.5, at=3.5, font=2, las=1, adj=0.5)
#
mtext("Events/N",side=2, line=9, at=3.25, font=2, las=1, adj=0.5)
mtext("Events/N",side=2, line=3.5, at=3.25, font=2, las=1, adj=0.5)
mtext(paste(ev_table[,1], ev_table[,2], sep="/"),side=2, line=9, at=3:1, las=1, adj=0.5)
mtext(paste(ev_table[,3], ev_table[,4], sep="/"),side=2, line=3.25, at=3:1, las=1, adj=0.5)
#
mtext("HR (95% CI)",side=4, line=4, at=3.5, font=2, las=1, adj=0.5)

format.ci <- function(x, dig=2) paste0(format(round(x[1], digits=dig), nsmall = dig), " (", format(round(x[2], digits=dig), nsmall = dig), " to ", format(round(x[3], digits=dig), nsmall = dig), ")")

mtext(apply(HRa_quantile, 2, format.ci), side=4, las=1, at=3:1, line=4, adj=0.5)
mtext(substitute(paste("Interaction ", italic(P), " = 0.46")), side=4, las=1, at=0.5, line=4, adj=0.5)

dev.copy2pdf(file="forrest.pdf")      
