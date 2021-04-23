## Start from data prep
rm(list=ls()[ls()!="rcts"])

factor_var <- which(sapply(sapply(rcts, class), function(x) x[1])=="character")

for(i in factor_var[-1]){
        rcts[,i] <- as.factor(rcts[,i])
}

### données des CJS gaz du sang et intervention pour prévenir l'acidose dans AKIKI 

###
library(mice)
# Utilisons mice avec l'argument maxit=0 pour produire une matrice de prédiction sans faire aucune imputation
imp <- mice(rcts, maxit=0)

# Stockons cette matrice de prédiction dans l'objet preM
predM <- imp$predictorMatrix

# Définissons les variables que l'on ne souhaite pas utiliser pour les imputations
nonpredictors <- c("id","naiss","datehospi", "datehospirea", "datehospirea",
                   "datesortierea","etatsortierea","derniere.nouvelles.censureJ60",
                   "etat.censureJ60", "date.rando", "date.eer", "eer.initiation",
                   "type.eer", "rrtfree", "ventilfree", "vasofree", "reprisediurese", "lastrrt", "surviesanseer")

# Attribuons la valeur zéro aux colonnes qui leur correspondent dans la matrice de prédiction
predM[,nonpredictors]<-0

# Récupérons également les méthodes d'imputation que mice prévois d'utiliser pour chaque variable
meth <-  imp$method

# A visée illustraive, choissons par exemple de ne pas imputer les ids
meth["hospitalmortality"]<-""

# Définissons une graine pour la reproductibilité
set.seed(123)

# Utilisons mice avec la nouvelle matrice de prediction
# L'argument m=4 signifie que nous souhaitons imputer 4 datasets complets
# Cette commande parfois peut être longue à exécuter
impdata<-mice(rcts, m=4, predictorMatrix = predM, method = meth)

save(impdata, file="imprcts.RData")

# creation d'un dataset imputé  akiki& ideal-icu complet pour lissage sur pH
rcts_imp1 <- complete(impdata, 1)
rcts_imp1$ph_e_p <- apply(sapply(with(impdata, ph_e)$analyses,c),1,mean)

pco2_calc <- with(impdata, bicar_e/(.03*10^(ph_e-6.1)) )
pco2_critera <- apply(sapply(pco2_calc$analyses, c),1, mean)<=40
ph_critera <- rcts$ph_e<=7.30
included <- as.logical(pco2_critera*ph_critera)

imp_included <- filter(impdata, included)
# Stokons le premier dataset imuputé dans l'objet imp1 avec la commande complete
imp1 <- complete(imp_included, 1)

library(survival)
library(survminer)
survfit.obj<-survfit(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~bras, data = imp1)
summary(survfit.obj)

splots<-ggsurvplot(survfit.obj,
                   ggtheme = theme_survminer(),
                   axes.offset=FALSE,
                   title    = "",
                   legend.title = "",
                   legend.labs = c("Delayed \nstrategy", "Early \nstrategy"),
                   #legend = c(0.1, 0.15),
                   palette = c("#0099B4FF", "#AD002AFF"),
                   conf.int = F,
                   pval=F, pval.method = F,
                   risk.table = TRUE, 
                   break.time.by = 10,
                   tables.theme = theme_cleantable(), 
                   tables.y.text = T,
                   tables.y.text.col= FALSE,
                   xlim=c(0,62))

wl <- 4.23
dev.new(width=wl, height=wl*1.5, pointsize=7)
splots + labs(x  = "Days since RRT initiation", y = "Proportion Surviving (%)")

dev.copy2pdf(file="fig1.0.pdf")

library(gtools)
set.seed(984)
imp1$ph_quantile <- quantcut(imp1$ph_e+rnorm(nrow(imp1), mean=0, sd=10^-6), q=3)

##
mod <- coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = imp1)
HR_main <- exp(c(coef(mod), confint(mod)))

HR_quantile <- sapply(levels(imp1$ph_quantile), function(x) {
mod <- coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE"), data = imp1, subset=ph_quantile==x)
exp(c(coef(mod)[1], confint(mod)[1,]))
})

## adjusted HR
mod <- coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE")+age+sofa_e, data = imp1)
HR_a <- exp(c(coef(mod)[1], confint(mod)[1,]))

HRa_quantile <- sapply(levels(imp1$ph_quantile), function(x) {
        mod <- coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE")+age+sofa_e, data = imp1, subset=ph_quantile==x)
        exp(c(coef(mod)[1], confint(mod)[1,]))
})
##
ev_table <- cbind(
matrix(unlist(tapply(imp1$etat.censureJ60[imp1$bras=="STRATEGIE PRECOCE"], imp1$ph_quantile[imp1$bras=="STRATEGIE PRECOCE"], function(x) {
        c(sum(x),length(x)) })
), 3, 2, byrow=T)
,
matrix(unlist(tapply(imp1$etat.censureJ60[imp1$bras=="STRATEGIE D ATTENTE"], imp1$ph_quantile[imp1$bras=="STRATEGIE D ATTENTE"], function(x) {
        c(sum(x),length(x)) })
), 3, 2, byrow=T)
)

rownames(ev_table) <- levels(imp1$ph_quantile)
colnames(ev_table) <- c("precoce_ev","precoce_N", "tardif_ev","tardif_N")

##

d60_km_death <- c(
1-tail(weightedsurvivaltables(survfit.obj$time, survfit.obj$surv, group="A"),1)[2],
1-tail(weightedsurvivaltables(survfit.obj$time, survfit.obj$surv, group="B"),1)[2]
)
names(d60_km_death) <- c("delayed", "early")

library(boot)
res <- boot(imp1, ard_km, R=999, strata=as.numeric(imp1$bras=="STRATEGIE PRECOCE"))
boot.ci(res)
d60_km_death <- c(d60_km_death, d60_km_death[2]-d60_km_death[1],boot.ci(res)$bca[c(4,5)])
names(d60_km_death) <- c("delayed", "early", "ard early-delayed","l_ci", "u_ci")

d60_km_death <- c(d60_km_death, d60_km_death[2]/d60_km_death[1])
res <- boot(imp1, rr_km, R=999, strata=as.numeric(imp1$bras=="STRATEGIE PRECOCE"))
boot.ci(res)
d60_km_death <- c(d60_km_death, boot.ci(res)$bca[c(4,5)])
names(d60_km_death) <- c("delayed", "early", "ard early-delayed","l_ci", "u_ci", "rr early/late", "rr_l_ci", "rr_u_ci")
d60_km_death

# export primary outcome results
write.csv(round(cbind(t(as.data.frame(tapply(imp1$etat.censureJ60, imp1$bras, function(x) table(x)["1"] ))), t(as.data.frame(d60_km_death))),2), "po.csv")
#
ARD_quantile <- sapply(levels(imp1$ph_quantile), function(x) {
        survfit.obj<-survfit(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~bras, data = imp1[imp1$ph_quantile==x,])
        d60_km <- c(
                1-tail(weightedsurvivaltables(survfit.obj$time, survfit.obj$surv, group="A"),1)[2],
                1-tail(weightedsurvivaltables(survfit.obj$time, survfit.obj$surv, group="B"),1)[2]
        )
        d60_km <- c(d60_km,d60_km[2]-d60_km[1])
        res <- boot(imp1[imp1$ph_quantile==x,], ard_km, R=999, strata=as.numeric(imp1[imp1$ph_quantile==x,]$bras=="STRATEGIE PRECOCE"))
        d60_km <- c(d60_km, boot.ci(res)$bca[c(4,5)])
        names(d60_km) <- c("delayed", "early", "ard early-delayed","l_ci", "u_ci")
        d60_km
})

RR_quantile <- sapply(levels(imp1$ph_quantile), function(x) {
        survfit.obj<-survfit(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~bras, data = imp1[imp1$ph_quantile==x,])
        d60_km <- c(
                1-tail(weightedsurvivaltables(survfit.obj$time, survfit.obj$surv, group="A"),1)[2],
                1-tail(weightedsurvivaltables(survfit.obj$time, survfit.obj$surv, group="B"),1)[2]
        )
        d60_km <- c(d60_km,d60_km[2]/d60_km[1])
        res <- boot(imp1[imp1$ph_quantile==x,], rr_km, R=999, strata=as.numeric(imp1[imp1$ph_quantile==x,]$bras=="STRATEGIE PRECOCE"))
        d60_km <- c(d60_km, boot.ci(res)$bca[c(4,5)])
        names(d60_km) <- c("delayed", "early", "rr early/delayed","l_ci", "u_ci")
        d60_km
})

##### Smoothing for event rate plot
# create interaction variables
probs.int.t <- with(rcts_imp1, ifelse(bras=="STRATEGIE PRECOCE", ph_e_p-mean(ph_e_p), 0))
probs.c<-with(rcts_imp1, ph_e_p-mean(ph_e_p))

prec <- 100 # n of mesures along ph_e_p prec as in precision
seq_probs.c <- modelr::seq_range(probs.c, n=prec)

nkn <- 1
kn<-quantile(probs.c, probs = seq(0, 1, length=nkn+2))[-c(1,nkn+2)]
rg <-range(probs.c)
library(splines)

fit.int.er<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE')
                  + ns(probs.c, knots = kn, Boundary.knots=rg)
                  + ns(probs.int.t, knots = kn, Boundary.knots=rg),
                  data = rcts_imp1, x=TRUE)

ypred <- list()
ypred$fit<-as.numeric(1-tail(survfit(fit.int.er, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=seq_probs.c, probs.int.t=seq_probs.c))$surv,1))
ypred$se<-as.numeric(tail(survfit(fit.int.er, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=seq_probs.c, probs.int.t=seq_probs.c))$std.err,1))

ypred2 <- list()
ypred2$fit<-as.numeric(1-tail(survfit(fit.int.er, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=seq_probs.c, probs.int.t=0))$surv,1))
ypred2$se<-as.numeric(tail(survfit(fit.int.er, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=seq_probs.c, probs.int.t=0))$std.err,1))

yy.er <- ypred$fit + outer(ypred$se, c(0, -qnorm(.975), qnorm(.975)), '*')
yy2.er<- ypred2$fit + outer(ypred2$se, c(0, -qnorm(.975), qnorm(.975)), '*')

wl <- 4.23
dev.new(width=wl, height=wl*1.5, pointsize=7)
xl <- seq_probs.c+mean(rcts_imp1$ph_e_p)>7.15 & seq_probs.c+mean(rcts_imp1$ph_e_p)<7.55 

matplot((seq_probs.c+mean(rcts_imp1$ph_e_p))[xl], (as.matrix(cbind(yy.er,yy2.er)))[xl,], type='l',
        bty="n", xaxt="n", yaxt="n",
        lty=c(1,3,3), col=c(rep("#B24745FF",3), rep("#374E55FF",3)), lwd=2,
        xlab="pH at enrollment", ylab="60 Days Death Rate", ylim=c(0,1), xlim=c(7.15, 7.57))

legend("bottomleft", box.lty=0,
       lty=1, lwd=2, col=c("#374E55FF", "#B24745FF"), c("Delayed", "Early"))

axis(1, at=seq(7.15, 7.55, by=.05))
axis(2, at=c(seq(0, 1, by=.2),1),las=1)
#dev.copy2pdf(file="eventrateplot.pdf")

##### HR Smoothing

ypred<-predict(fit.int.er, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=seq_probs.c, probs.int.t=seq_probs.c), se=TRUE)
ypred2<-predict(fit.int.er, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=seq_probs.c, probs.int.t=0), se=TRUE)

yy <- ypred$fit + outer(ypred$se, c(0, -qnorm(.975), qnorm(.975)), '*')
yy2<- ypred2$fit + outer(ypred2$se, c(0, -qnorm(.975), qnorm(.975)), '*')

matplot(seq_probs.c+mean(rcts_imp1$ph_e_p), exp(matrix(cbind(yy,yy2), ncol=6)), type='l', lty=c(1,2,2), col=c(1,1,1,2,2,2),
        lwd=2, log='y', xlim=c(7.15, 7.57),
        xlab="pH at enrollment", ylab="exp(LP)")


y0.precoce_<-data.frame(bras = 1, ns(seq_probs.c, knots = kn, Boundary.knots=rg), ns(seq_probs.c, knots = kn, Boundary.knots=rg)) 
y0.tardif_<-data.frame(bras = 0, ns(seq_probs.c, knots = kn, Boundary.knots=rg), ns(0, knots = kn, Boundary.knots=rg))


l<-sapply(1:prec, function(i){
        matrix(y0.precoce_[i,] - y0.tardif_[i,]) })

l<-matrix(as.numeric(l), nrow(l), prec)

fits<-as.numeric(t(l) %*% coef(fit.int.er))
se_s<-sapply(1:length(seq_probs.c), function(i){
        l = as.numeric(y0.precoce_[i,] - y0.tardif_[i,])
        t(l) %*% vcov(fit.int.er) %*% l
})

yy.hr <- fits + outer(sqrt(se_s), c(0, -qnorm(.975), qnorm(.975)), '*')

wl <- 4.23
dev.new(width=wl, height=wl*1.5, pointsize=7)

plot(1, bty="n", type="n", xaxt="n", yaxt="n", log='y',
     xlim=c(7.15, 7.57), ylim=c(.7, 3.5),
     xlab="pH at enrollment", ylab="Hazard Ratio")

axis(1, at=seq(7.15, 7.55, by=.05))
axis(2, at=c(seq(.7,3.5, by=.7),1),las=1)
lines(c(0,7.57), y=c(1,1), lwd=1)
arrows(7.56, c(1.05,1/1.05), 7.56, c(1.05+.3, 1/(1.05+.3)), length = 0.07, xpd=TRUE)
mtext(text="favors delayed", side=4,line=0, at=1.05, cex=3/4, adj=0)
mtext(text="favors early", side=4,line=-0, at=2-1.05, cex=3/4, adj=1)

xl <- seq_probs.c+mean(rcts_imp1$ph_e_p)>7.15 & seq_probs.c+mean(rcts_imp1$ph_e_p)<7.55 
points((seq_probs.c+mean(rcts_imp1$ph_e_p))[xl], exp(yy.hr[xl,1]), type = 'l', lwd=2, col="#79af97")
polygon(c((seq_probs.c+mean(rcts_imp1$ph_e_p))[xl], rev((seq_probs.c+mean(rcts_imp1$ph_e_p))[xl])), c(exp(yy.hr[xl,2]), rev(exp(yy.hr[xl,3]))),
        col= rgb(55,78,85, alpha = 38, maxColorValue=255), border=NA)
#dev.copy2pdf(file="hrplot.pdf")

##### ARD Smoothing
probs.int.t <- with(rcts_imp1, ifelse(bras=="STRATEGIE PRECOCE", ph_e_p-mean(ph_e_p), 0))
probs.c<-with(rcts_imp1, ph_e_p-mean(ph_e_p))

nkn <- 1
kn<-quantile(probs.c, probs = seq(0, 1, length=nkn+2))[-c(1,nkn+2)]

a1<-ns(probs.c, knots = kn, Boundary.knots=rg)[,1]
a2<-ns(probs.c, knots = kn, Boundary.knots=rg)[,2]

b1<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,1]
b2<-ns(probs.int.t, knots = kn, Boundary.knots=rg)[,2]


fit.int.ard<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + a1+a2 + b1+b2, data = rcts_imp1, x=TRUE)
#fit.int.ard<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE')*ph_e_p, data = rcts_imp1, x=TRUE)

#y0.precoce_<-data.frame(bras = "STRATEGIE PRECOCE", ns(seq_probs.c, knots = kn, Boundary.knots=rg), ns(seq_probs.c, knots = kn, Boundary.knots=rg))

seq_ph <- modelr::seq_range(rcts_imp1$ph_e_p, n=prec)
y0.precoce_<-data.frame(bras = "STRATEGIE PRECOCE", ns(seq_probs.c, knots = kn, Boundary.knots=rg), ns(seq_probs.c, knots = kn, Boundary.knots=rg))
names(y0.precoce_) <- c("bras","a1", "a2", "b1", "b2")

y0.tardif_<-data.frame(bras = "STRATEGIE D ATTENTE", ns(seq_probs.c, knots = kn, Boundary.knots=rg), ns(0, knots = kn, Boundary.knots=rg))
names(y0.tardif_) <- c("bras","a1", "a2", "b1", "b2")

precoced60d<-1-tail(survfit(fit.int.ard, newdata=y0.precoce_)$surv,1)

tardifd60d<-1-tail(survfit(fit.int.ard, newdata=y0.tardif_)$surv,1)

se_s_pooled<-ardbootfunction(dataset=rcts_imp1, nboot = 1000)


wl <- 4.23
dev.new(width=wl, height=wl*1.5, pointsize=7)

xl <- seq_ph>7.15 & seq_ph<7.55 

plot(1, xaxt="n", yaxt="n", main="", xlab="pH at enrollment", ylab="Absolute Risk Reduction",
     xlim=c(7.15, 7.57), ylim=c(-.5,.5), las=1, bty="n")

axis(1, at=seq(7.15, 7.55, by=.05))
axis(2, at=seq(-.5,.5, by=.1),las=1)
lines(c(0,7.57), y=c(0,0), lwd=1)

arrows(7.56, c(0.02,-0.02), 7.56, c(0.02+.1, -0.02-.1), length = 0.07, xpd=TRUE)
mtext(text="favors delayed", side=4,line=0, at=0.02, cex=3/4, adj=0)
mtext(text="favors early", side=4,line=-0, at=-0.02, cex=3/4, adj=1)

points(seq_ph[xl], (precoced60d-tardifd60d)[xl], type="l", lwd=2, col="#79af97")
lines(seq_probs.c, (precoced60d-tardifd60d_)[xl]-qnorm(.975)*se_s_pooled[xl], lty=2)
lines(seq_probs.c, (precoced60d_pooled-tardifd60d_pooled)+qnorm(.975)*se_s_pooled, lty=2)

polygon(c(seq_ph[xl], rev(seq_ph[xl])), c((precoced60d-tardifd60d)[xl]+qnorm(.975)*se_s_pooled[xl],
                                                                          rev((precoced60d-tardifd60d)[xl]-qnorm(.975)*se_s_pooled[xl])),
        col= rgb(55,78,85, alpha = 38, maxColorValue=255), border=NA)
#dev.copy2pdf(file="ardplot.pdf")

#### Hospital mortality (no imputation of outcomes, complete case analysis here)
temp <- imp1
hospitalmortality_val<-c(
                sum(temp[temp$bras=="STRATEGIE PRECOCE","hospitalmortality"]=="Décédé", na.rm = T),
                
                sum(temp[temp$bras=="STRATEGIE PRECOCE","hospitalmortality"]=="Décédé", na.rm = T)/length(temp[temp$bras=="STRATEGIE PRECOCE","hospitalmortality"]),
                
                sum(temp[temp$bras=="STRATEGIE D ATTENTE","hospitalmortality"]=="Décédé", na.rm = T),
                sum(temp[temp$bras=="STRATEGIE D ATTENTE","hospitalmortality"]=="Décédé", na.rm = T)/length(temp[temp$bras=="STRATEGIE D ATTENTE","hospitalmortality"]),
                
                (sum(temp[temp$bras=="STRATEGIE PRECOCE","hospitalmortality"]=="Décédé", na.rm = T)/length(temp[temp$bras=="STRATEGIE PRECOCE","hospitalmortality"]))-
                (sum(temp[temp$bras=="STRATEGIE D ATTENTE","hospitalmortality"]=="Décédé", na.rm = T)/length(temp[temp$bras=="STRATEGIE D ATTENTE","hospitalmortality"])),
                
                (sum(temp[temp$bras=="STRATEGIE PRECOCE","hospitalmortality"]=="Décédé", na.rm = T)/length(temp[temp$bras=="STRATEGIE PRECOCE","hospitalmortality"]))/
                (sum(temp[temp$bras=="STRATEGIE D ATTENTE","hospitalmortality"]=="Décédé", na.rm = T)/length(temp[temp$bras=="STRATEGIE D ATTENTE","hospitalmortality"]))
        )
res_ard <- boot(temp, binout_ard_boot, R=999)
res_rr <- boot(temp, binout_rr_boot, R=999)

hospitalmortality_val <- c(hospitalmortality_val, boot.ci(res_ard)$bca[c(4,5)], boot.ci(res_rr)$bca[c(4,5)])
names(hospitalmortality_val) <- c("precoce_ev", "precoce_ra", "tardif_ev", "tardif_ra", "precoce-tardif", "precoce/tardif", "ard_lci", "ard_uci", "rr_lci", "rr_uci")
hospitalmortality_val

#### ICU mortality (no missing data here : no imputation needed)
temp <- imp1
icumortality_val<-c(
        sum(temp[temp$bras=="STRATEGIE PRECOCE","etatsortierea"]=="Décédé", na.rm = T),
        mean(temp[temp$bras=="STRATEGIE PRECOCE","etatsortierea"]=="Décédé", na.rm = T),
        sum(temp[temp$bras=="STRATEGIE D ATTENTE","etatsortierea"]=="Décédé", na.rm = T),
        mean(temp[temp$bras=="STRATEGIE D ATTENTE","etatsortierea"]=="Décédé", na.rm = T),
        mean(temp[temp$bras=="STRATEGIE PRECOCE","etatsortierea"]=="Décédé", na.rm = T)-mean(temp[temp$bras=="STRATEGIE D ATTENTE","etatsortierea"]=="Décédé", na.rm = T),
        mean(temp[temp$bras=="STRATEGIE PRECOCE","etatsortierea"]=="Décédé", na.rm = T)/mean(temp[temp$bras=="STRATEGIE D ATTENTE","etatsortierea"]=="Décédé", na.rm = T)
)
res_ard <- boot(temp, binout_ard_boot, R=999, cj="etatsortierea", event="Décédé")
res_rr <- boot(temp, binout_rr_boot, R=999, cj="etatsortierea", event="Décédé")

icumortality_val <- c(icumortality_val, boot.ci(res_ard)$bca[c(4,5)], boot.ci(res_rr)$bca[c(4,5)])
names(icumortality_val) <- c("precoce_ev", "precoce_ra", "tardif_ev", "tardif_ra", "precoce-tardif", "precoce/tardif", "ard_lci", "ard_uci", "rr_lci", "rr_uci")
icumortality_val

### export table hospital and ICU mortality
mortl_exp <- rbind(hospitalmortality_val, icumortality_val)
mortl_exp <- mortl_exp[,!(colnames(mortl_exp) %in% c("precoce-tardif", "ard_lci", "ard_uci"))]
mortl_exp <- round(mortl_exp,3)
mortl_exp[,2] <- mortl_exp[,2]*100
mortl_exp[,4] <- mortl_exp[,4]*100

mortl_exp <- cbind(
paste0(mortl_exp[,1], " (", mortl_exp[,2], "%)"),
paste0(mortl_exp[,3], " (", mortl_exp[,4], "%)"),
format(round(mortl_exp[,5],2), nsmall = 2),
paste0("(",format(round(mortl_exp[,6],2), nsmall = 2), " to ", format(round(mortl_exp[,7],2), nsmall = 2), ")")
)

write.csv(mortl_exp, "mortable.csv")

#### Survie sans EER
rr <- data.frame(val=rep(NA, imp_included$m), var=rep(NA, imp_included$m))

for(i in 1:imp_included$m) {
temp <- complete(imp_included, i)
x <- tapply(temp$surviesanseer, temp$bras, function(x) c(sum(x=="OUI"), mean(x=="OUI")))
rr$nprecoce[i] <- x[["STRATEGIE PRECOCE"]][1]
rr$ntardif[i]<- x[["STRATEGIE D ATTENTE"]][1]
rr$val[i] <- x[["STRATEGIE PRECOCE"]][2]/x[["STRATEGIE D ATTENTE"]][2]
res <- boot(temp, binout_rr_boot, R=999, cj="surviesanseer", event="OUI")
rr$var[i] <- var(res$t)
}
rubinvar <- rubinr(thetas = rr$val, vars = rr$var)

surviesanseer_pool <- c(mean(rr$val), mean(rr$val)+c(-1,1)*sqrt(rubinvar)*qnorm(.975))

#### Ventilator free days
ventilfree_val <- list()
for (i in 1:imp_included$m) {
temp <- complete(imp_included, i)
ventilfree_val[[i]]<-c(
        mean(temp[temp$bras=="STRATEGIE PRECOCE",]$ventilfree), sd(temp[temp$bras=="STRATEGIE PRECOCE",]$ventilfree),
        mean(temp[temp$bras=="STRATEGIE D ATTENTE",]$ventilfree), sd(temp[temp$bras=="STRATEGIE D ATTENTE",]$ventilfree),
        mean(temp[temp$bras=="STRATEGIE PRECOCE",]$ventilfree)-mean(temp[temp$bras=="STRATEGIE D ATTENTE",]$ventilfree)
)
res <- boot(temp, mdif, R=999, cj="ventilfree", strata=as.numeric(imp1$bras=="STRATEGIE PRECOCE"))
ventilfree_val[[i]] <- c(ventilfree_val[[i]], var(res$t))
names(ventilfree_val[[i]]) <- c("precoce","sd_precoce", "tardif","sd_tardif", "precoce-tardif", "var")
}
ventilfree_val <- sapply(ventilfree_val, c)
vars <- ventilfree_val["var",]
ventilfree_val_pool <- apply(ventilfree_val, 1, mean)

rubinvar <- rubinr(thetas = ventilfree_val["precoce-tardif",], vars = ventilfree_val["var",])
ventilfree_val_pool <- c(ventilfree_val_pool, rubinvar, ventilfree_val_pool["precoce-tardif"]+c(-1,1)*qnorm(.975)*sqrt(rubinvar))
names(ventilfree_val_pool)[c(7:9)] <- c("rubinvar", "lci", "uci")
ventilfree_val_pool
vfd <- t(round(ventilfree_val_pool[-c(6,7)],2))

#### RRT free days
rrtfree_val <- list()
for (i in 1:imp_included$m) {
        temp <- complete(imp_included, i)
        rrtfree_val[[i]]<-c(
                mean(temp[temp$bras=="STRATEGIE PRECOCE",]$rrtfree), sd(temp[temp$bras=="STRATEGIE PRECOCE",]$rrtfree),
                mean(temp[temp$bras=="STRATEGIE D ATTENTE",]$rrtfree), sd(temp[temp$bras=="STRATEGIE D ATTENTE",]$rrtfree),
                mean(temp[temp$bras=="STRATEGIE PRECOCE",]$rrtfree)-mean(temp[temp$bras=="STRATEGIE D ATTENTE",]$rrtfree)
        )
        res <- boot(temp, mdif, R=999,  cj="rrtfree", strata=as.numeric(imp1$bras=="STRATEGIE PRECOCE"))
        rrtfree_val[[i]] <- c(rrtfree_val[[i]], var(res$t))
        names(rrtfree_val[[i]]) <- c("precoce","sd_precoce", "tardif","sd_tardif", "precoce-tardif", "var")
}
rrtfree_val <- sapply(rrtfree_val, c)
vars <- rrtfree_val["var",]
rrtfree_val_pool <- apply(rrtfree_val, 1, mean)

rubinvar <- rubinr(thetas = rrtfree_val["precoce-tardif",], vars = rrtfree_val["var",])
rrtfree_val_pool <- c(rrtfree_val_pool, rubinvar, rrtfree_val_pool["precoce-tardif"]+c(-1,1)*qnorm(.975)*sqrt(rubinvar))
names(rrtfree_val_pool)[c(7:9)] <- c("rubinvar", "lci", "uci")
rrtfree_val_pool
rfd <- t(round(rrtfree_val_pool[-c(6,7)],1))

#### Vaso free days
vasofree_val <- list()
for (i in 1:imp_included$m) {
        temp <- complete(imp_included, i)
        vasofree_val[[i]]<-c(
                mean(temp[temp$bras=="STRATEGIE PRECOCE",]$vasofree), sd(temp[temp$bras=="STRATEGIE PRECOCE",]$vasofree),
                mean(temp[temp$bras=="STRATEGIE D ATTENTE",]$vasofree), sd(temp[temp$bras=="STRATEGIE D ATTENTE",]$vasofree),
                mean(temp[temp$bras=="STRATEGIE PRECOCE",]$vasofree)-mean(temp[temp$bras=="STRATEGIE D ATTENTE",]$vasofree)
        )
        res <- boot(temp, mdif, R=999,  cj="vasofree", strata=as.numeric(imp1$bras=="STRATEGIE PRECOCE"))
        vasofree_val[[i]] <- c(vasofree_val[[i]], var(res$t))
        names(vasofree_val[[i]]) <- c("precoce","sd_precoce", "tardif","sd_tardif", "precoce-tardif", "var")
}
vasofree_val <- sapply(vasofree_val, c)
vars <- vasofree_val["var",]
vasofree_val_pool <- apply(vasofree_val, 1, mean)

rubinvar <- rubinr(thetas = vasofree_val["precoce-tardif",], vars = vasofree_val["var",])
vasofree_val_pool <- c(vasofree_val_pool, rubinvar, vasofree_val_pool["precoce-tardif"]+c(-1,1)*qnorm(.975)*sqrt(rubinvar))
names(vasofree_val_pool)[c(7:9)] <- c("rubinvar", "lci", "uci")
vasofree_val_pool
vsfd <- t(round(vasofree_val_pool[-c(6,7)],1))

### export free days outcome table
fd_table <- rbind(rfd,vfd,vsfd)
rownames(fd_table) <- c("RRT", "Ventilator", "Vasopressors")
fd_exp <- cbind(
        paste0(format(fd_table[,1],nsmall = 2), " (", format(fd_table[,2],nsmall = 2), ")"),
        paste0(format(fd_table[,3],nsmall = 2), " (", format(fd_table[,4],nsmall = 2), ")"),
        format(fd_table[,5],nsmall = 2),
        paste0("(", format(fd_table[,6],nsmall = 2), " to ", format(fd_table[,7],nsmall = 2), ")")
      )

write.csv(fd_exp, "freedaysoutcomes.csv")

#### Length of ICU stay
imp1$icustay <- as.numeric(round(with(imp1, difftime(datesortierea, date.rando, unit="days"))))

icustay <- tapply(imp1$icustay, imp1$bras, mean)
icustay_sd <- tapply(imp1$icustay, imp1$bras, sd)
icustay_pool <- icustay["STRATEGIE PRECOCE"]-icustay["STRATEGIE D ATTENTE"]
res <- boot(imp1, mdif, cj="icustay", R=999)

icustay_pool <- c(icustay["STRATEGIE PRECOCE"],icustay_sd["STRATEGIE PRECOCE"], icustay["STRATEGIE D ATTENTE"],icustay_sd["STRATEGIE D ATTENTE"], icustay_pool, boot.ci(res)$bca[c(4,5)])
names(icustay_pool) <- c("precoce","precoce_sd", "tardif","tardif_sd", "md", "lci", "uci")

icustay_exp <- round(icustay_pool, 2)

icustay_exp <- c(paste0(format(icustay_exp[1], nsmall = 2), " (", format(icustay_exp[2], nsmall = 2), ")"),
  paste0(format(icustay_exp[3], nsmall = 2), " (", format(icustay_exp[4], nsmall = 2), ")"),
  format(icustay_exp[5], nsmall = 2),
  paste0("(", format(icustay_exp[6], nsmall = 2), " to ", format(icustay_exp[7], nsmall = 2), ")")
  )

#write.csv(t(icustay_exp), "icustay.csv")

### table(s) one
library(tableone)

imp1_table1 <- data.frame(
        arm=imp1$bras,
        study=imp1$etude,
        age=imp1$age,
        female=imp1$sexe=="Féminin",
        baseline_creat=imp1$creat_b,
        hypertension=imp1$hypertension,
        dibatetes=imp1$diabetes!="Aucun diabète OU diabète\ntraité par régime uniquement",
        cirrhosis=imp1$cirrhosis=="OUI" | imp1$cirrhosis=="Modérée (sans hypertension\nportale, en incluant les\nhépatites chroniques)",
        respiratory_disease=imp1$respiratory_disease,
        cancer=imp1$cancer!="NON",
        aides=imp1$aides,
        immunosupressive_drug=imp1$immunosupressive_drug,
        organ_graft=imp1$organ_graft,
        hemopathy=imp1$hemopathy,
        sofa_e=imp1$sofa_e,
        kidney_sofa=imp1$kidney_sofa,
        hemodynamic_sofa=imp1$hemodynamic_sofa,
        bilirubin_sofa=imp1$bilirubin_sofa,
        platelet_sofa=imp1$platelet_sofa,
        gcs_sofa=imp1$gcs_sofa,
        respi_sofa=imp1$respi_sofa, 
        weight_e=imp1$weight_e,
        creat_e=imp1$creat_e,
        urea_e=imp1$urea_e,
        pot_e=imp1$pot_e,
        bicar_e=imp1$bicar_e,
        ph_e=imp1$ph_e
)


table1<- CreateTableOne(vars=colnames(imp1_table1)[-c(1,2)], strata="arm", data=imp1_table1, test=TRUE)
printed <- print(table1,smd=FALSE)
write.csv(printed, file = "table1.csv")

table1<- CreateTableOne(vars=colnames(imp1_table1)[-c(1,2)], strata="arm", data=imp1_table1[imp1_table1$study=="akiki",], test=TRUE)
printed <- print(table1,smd=FALSE)
write.csv(printed, file = "table1_akiki.csv")

table1<- CreateTableOne(vars=colnames(imp1_table1)[-c(1,2)], strata="arm", data=imp1_table1[imp1_table1$study=="idealicu",], test=TRUE)
printed <- print(table1,smd=FALSE)
write.csv(printed, file = "table1_ideal.csv")
