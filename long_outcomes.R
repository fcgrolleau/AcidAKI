akikidbsg <- read.csv("~/Desktop/TEAM METHODS/phd/Trials Data/akiki data/SGakiki.csv", sep=";", comment.char="#")
long_akiki<-data.frame(
        etude="akiki",
        id=substring(akikidbsg$subject_id, 5, nchar(akikidbsg$subject_id)-3),
        bras=akikidbsg$inc_bras,
        phmin_j0=akikidbsg$j0_car_phmin,
        phmin_j1=akikidbsg$j1_car_phmin,
        phmin_j2=akikidbsg$j2_car_phmin,
        phmin_j3=akikidbsg$j3_car_phmin,
        phmin_j4=akikidbsg$j4_car_phmin,
        phmin_j5=akikidbsg$j5_car_phmin,
        phmin_j6=akikidbsg$j6_car_phmin,
        phmin_j7=akikidbsg$j7_car_phmin,
        
        bicarmin_j0=akikidbsg$j0_car_hco3min,
        bicarmin_j1=akikidbsg$j1_car_hco3min,
        bicarmin_j2=akikidbsg$j2_car_hco3min,
        bicarmin_j3=akikidbsg$j3_car_hco3min,
        bicarmin_j4=akikidbsg$j4_car_hco3min,
        bicarmin_j5=akikidbsg$j5_car_hco3min,
        bicarmin_j6=akikidbsg$j6_car_hco3min,
        bicarmin_j7=akikidbsg$j7_car_hco3min,
        
        pco2mon_j0=akikidbsg$j0_car_paco2min,
        pco2mon_j1=akikidbsg$j1_car_paco2min,
        pco2mon_j2=akikidbsg$j2_car_paco2min,
        pco2mon_j3=akikidbsg$j3_car_paco2min,
        pco2mon_j4=akikidbsg$j4_car_paco2min,
        pco2mon_j5=akikidbsg$j5_car_paco2min,
        pco2mon_j6=akikidbsg$j6_car_paco2min,
        pco2mon_j7=akikidbsg$j7_car_paco2min,
        
        bicarven_j0=akikidbsg$j0_ele_q2,
        bicarven_j1=akikidbsg$j1_ele_q2,
        bicarven_j2=akikidbsg$j2_ele_q2,
        bicarven_j3=akikidbsg$j3_ele_q2,
        bicarven_j4=akikidbsg$j4_ele_q2,
        bicarven_j5=akikidbsg$j5_ele_q2,
        bicarven_j6=akikidbsg$j6_ele_q2,
        bicarven_j7=akikidbsg$j7_ele_q2
)

longout <- long_akiki[long_akiki$id %in% imp1$id[imp1$etude=="akiki"],]

int_mask <- longout[,(ncol(longout)-7):ncol(longout)]=="OUI"

int <- list()
for (i in 1:8) {
int[[i]] <- unlist(tapply(int_mask[,i], longout$bras, function(x) return(data.frame(Ev=sum(x, na.rm = T), N=sum(!is.na(x)))) ))
}

int <- t(sapply(int, c))
tab_int <- cbind(int[,3:4], int[,1:2])
rownames(tab_int) <- paste0("j", 0:7)
tab_int
write.csv(tab_int, "interv.csv")

reshape2::melt(longout[,2:12], id.vars="id",  measure.vars="bras")

library(tidyr)
melted <- longout[,2:11] %>% 
        pivot_longer(
                cols = starts_with("phmin_j"), 
                names_to = "jour", 
                values_to = "val",
                values_drop_na = F
        )

### pH box plot
        cx <- 1
        dx <- -3
        eq <- .33
        dev.new(width=6, height=6, pointsize=9)
        par(mfrow=c(1,1),mar=c(5,4,2,1)+.1, mgp=c(3,0.75,0), las=1, oma=c(1,1,0,0))
        boxplot(melted$val[melted$bras=="STRATEGIE PRECOCE"]~ melted$jour[melted$bras=="STRATEGIE PRECOCE"], boxwex=0.3, col="#00A1D5FF", axes=F, ylab="Lowest pH", xlab="")
        boxplot(melted$val[melted$bras=="STRATEGIE D ATTENTE"]~ melted$jour[melted$bras=="STRATEGIE D ATTENTE"], boxwex=0.3, col="#DF8F44FF", add=T, axes= F, at=1:8+eq)
        axis(1,at=1:8+eq/2, labels = c(0:7))
        #mtext(c(10,12,14), at=c(10,12,14), side=1, line=0.75, cex=cx)
        axis(2,at=seq(6.7,7.6,by=.1))
        mtext("Day", side=1, line=2, cex=cx)
        #mtext("A", line=0, side=3, at=dx, adj=0, font=2)
        legend(6.5,6.8,c("Early","Delayed"), fill=c("#00A1D5FF", "#DF8F44FF"), bty="n")
        mtext(table(is.na(melted$val[melted$bras=="STRATEGIE PRECOCE"]), melted$jour[melted$bras=="STRATEGIE PRECOCE"])[1,], side=1, at=1:8+eq/2, line=3.5, cex=cx)
        mtext(table(is.na(melted$val[melted$bras=="STRATEGIE D ATTENTE"]), melted$jour[melted$bras=="STRATEGIE D ATTENTE"])[1,], side=1, at=1:8+eq/2, line=4.5, cex=cx)
        mtext("N evaluated", side=1, line=2.5, at=-1, cex=cx,  adj=0)
        mtext("Early", side=1, line=3.5, at=-1, cex=cx,  adj=0)
        mtext("Delayed", side=1, line=4.5, at=-1, cex=cx,  adj=0)
 
dev.copy2pdf(file="long_ph.pdf")

### HCO3 box plot
        
melted_bic <- longout[,c(rep(T, 3), grepl('^bicarmin_j', names(longout))[4:ncol(longout)])] %>% 
        pivot_longer(
                cols = starts_with("bicarmin_j"), 
                names_to = "jour", 
                values_to = "val",
                values_drop_na = F
                ) 
cx <- 1
dx <- -3
dev.new(width=6, height=6, pointsize=9)
par(mfrow=c(1,1),mar=c(5,4,2,1)+.1, mgp=c(3,0.75,0), las=1, oma=c(1,1,0,0))
boxplot(melted_bic$val[melted_bic$bras=="STRATEGIE PRECOCE"]~ melted_bic$jour[melted_bic$bras=="STRATEGIE PRECOCE"], boxwex=0.3, col="#00A1D5FF", ylim=c(0,45), axes=F, ylab="Lowest serum bicarbonate (mmol/L)", xlab="")
boxplot(melted_bic$val[melted_bic$bras=="STRATEGIE D ATTENTE"]~ melted_bic$jour[melted_bic$bras=="STRATEGIE D ATTENTE"], boxwex=0.3, col="#DF8F44FF", add=T, axes= F, at=1:8+eq)
axis(1,at=1:8+eq/2, labels = c(0:7))
#mtext(c(10,12,14), at=c(10,12,14), side=1, line=0.75, cex=cx)
axis(2,at=seq(0,45,by=5))
mtext("Day", side=1, line=2, cex=cx)
#mtext("A", line=0, side=3, at=dx, adj=0, font=2)
legend(6.5,4,c("Early","Delayed"), fill=c("#00A1D5FF", "#DF8F44FF"), bty="n")
mtext(table(is.na(melted_bic$val[melted_bic$bras=="STRATEGIE PRECOCE"]), melted_bic$jour[melted_bic$bras=="STRATEGIE PRECOCE"])[1,], side=1, at=1:8+eq/2, line=3.5, cex=cx)
mtext(table(is.na(melted_bic$val[melted_bic$bras=="STRATEGIE D ATTENTE"]), melted_bic$jour[melted_bic$bras=="STRATEGIE D ATTENTE"])[1,], side=1, at=1:8+eq/2, line=4.5, cex=cx)
mtext("N evaluated", side=1, line=2.5, at=-1, cex=cx,  adj=0)
mtext("Early", side=1, line=3.5, at=-1, cex=cx,  adj=0)
mtext("Delayed", side=1, line=4.5, at=-1, cex=cx,  adj=0)

dev.copy2pdf(file="long_bicar.pdf")

### PCO2 box plot

melted_co2 <- longout[,c(rep(T, 3), grepl('^pco2mon_j', names(longout))[4:ncol(longout)])] %>% 
        pivot_longer(
                cols = starts_with("pco2mon_j"), 
                names_to = "jour", 
                values_to = "val",
                values_drop_na = F
        ) 
cx <- 1
dx <- -3
dev.new(width=6, height=6, pointsize=9)
par(mfrow=c(1,1),mar=c(5,4,2,1)+.1, mgp=c(3,0.75,0), las=1, oma=c(1,1,0,0))
boxplot(melted_co2$val[melted_co2$bras=="STRATEGIE PRECOCE"]~ melted_co2$jour[melted_co2$bras=="STRATEGIE PRECOCE"], boxwex=0.3, col="#00A1D5FF", ylim=c(10,80), axes=F, ylab="Lowest PaCO2 (mmHg)", xlab="")
boxplot(melted_co2$val[melted_co2$bras=="STRATEGIE D ATTENTE"]~ melted_co2$jour[melted_co2$bras=="STRATEGIE D ATTENTE"], boxwex=0.3, col="#DF8F44FF", add=T, axes= F, at=1:8+eq)
axis(1,at=1:8+eq/2, labels = c(0:7))
#mtext(c(10,12,14), at=c(10,12,14), side=1, line=0.75, cex=cx)
axis(2,at=seq(10,80,by=10))
mtext("Day", side=1, line=2, cex=cx)
#mtext("A", line=0, side=3, at=dx, adj=0, font=2)
legend(6.5,15,c("Early","Delayed"), fill=c("#00A1D5FF", "#DF8F44FF"), bty="n")
mtext(table(is.na(melted_co2$val[melted_co2$bras=="STRATEGIE PRECOCE"]), melted_co2$jour[melted_co2$bras=="STRATEGIE PRECOCE"])[1,], side=1, at=1:8+eq/2, line=3.5, cex=cx)
mtext(table(is.na(melted_co2$val[melted_co2$bras=="STRATEGIE D ATTENTE"]), melted_co2$jour[melted_co2$bras=="STRATEGIE D ATTENTE"])[1,], side=1, at=1:8+eq/2, line=4.5, cex=cx)
mtext("N evaluated", side=1, line=2.5, at=-1, cex=cx,  adj=0)
mtext("Early", side=1, line=3.5, at=-1, cex=cx,  adj=0)
mtext("Delayed", side=1, line=4.5, at=-1, cex=cx,  adj=0)

dev.copy2pdf(file="long_pco2.pdf")
