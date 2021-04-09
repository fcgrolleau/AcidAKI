load("explanatorydata.Rdata")
rm(list=ls()[ls()!="rcts"])


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
                   "type.eer", "rrtfree", "ventilfree", "vasofree", "reprisediurese", "lastrrt")

# Attribuons la valeur zéro aux colonnes qui leur correspondent dans la matrice de prédiction
predM[,nonpredictors]<-0

# Récupérons également les méthodes d'imputation que mice prévois d'utiliser pour chaque variable
meth <-  imp$method

# A visée illustraive, choissons par exemple de ne pas imputer les ids
meth["id"]<-""

# Définissons une graine pour la reproductibilité
set.seed(123)

# Utilisons mice avec la nouvelle matrice de prediction
# L'argument m=4 signifie que nous souhaitons imputer 4 datasets complets
# Cette commande parfois peut être longue à exécuter
impdata<-mice(rcts, m=4, predictorMatrix = predM, method = meth)

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
                   title    = "Patients With pH < 7.30 and PCO2 < 40mmHg",
                   legend.title = "",
                   legend.labs = c("Delayed strategy", "Early strategy"),
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
HRa_quantile <- sapply(levels(imp1$ph_quantile), function(x) {
        mod <- coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=="STRATEGIE PRECOCE")+age+sofa_e, data = imp1, subset=ph_quantile==x)
        exp(c(coef(mod)[1], confint(mod)[1,]))
})
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
                mean(temp[temp$bras=="STRATEGIE PRECOCE","hospitalmortality"]=="Décédé", na.rm = T),
                mean(temp[temp$bras=="STRATEGIE D ATTENTE","hospitalmortality"]=="Décédé", na.rm = T),
                mean(temp[temp$bras=="STRATEGIE PRECOCE","hospitalmortality"]=="Décédé", na.rm = T)-mean(temp[temp$bras=="STRATEGIE D ATTENTE","hospitalmortality"]=="Décédé", na.rm = T),
                mean(temp[temp$bras=="STRATEGIE PRECOCE","hospitalmortality"]=="Décédé", na.rm = T)/mean(temp[temp$bras=="STRATEGIE D ATTENTE","hospitalmortality"]=="Décédé", na.rm = T)
        )
res_ard <- boot(temp, binout_ard_boot, R=999)
res_rr <- boot(temp, binout_rr_boot, R=999)

hospitalmortality_val <- c(hospitalmortality_val, boot.ci(res_ard)$bca[c(4,5)], boot.ci(res_rr)$bca[c(4,5)])
names(hospitalmortality_val) <- c("precoce", "tardif", "precoce-tardif", "precoce/tardif", "ard_lci", "ard_uci", "rr_lci", "rr_uci")
hospitalmortality_val

#### ICU mortality (no missing data here : no imputation needed)
temp <- imp1
icumortality_val<-c(
        mean(temp[temp$bras=="STRATEGIE PRECOCE","etatsortierea"]=="Décédé", na.rm = T),
        mean(temp[temp$bras=="STRATEGIE D ATTENTE","etatsortierea"]=="Décédé", na.rm = T),
        mean(temp[temp$bras=="STRATEGIE PRECOCE","etatsortierea"]=="Décédé", na.rm = T)-mean(temp[temp$bras=="STRATEGIE D ATTENTE","etatsortierea"]=="Décédé", na.rm = T),
        mean(temp[temp$bras=="STRATEGIE PRECOCE","etatsortierea"]=="Décédé", na.rm = T)/mean(temp[temp$bras=="STRATEGIE D ATTENTE","etatsortierea"]=="Décédé", na.rm = T)
)
res_ard <- boot(temp, binout_ard_boot, R=999, cj="etatsortierea", event="Décédé")
res_rr <- boot(temp, binout_rr_boot, R=999, cj="etatsortierea", event="Décédé")

icumortality_val <- c(icumortality_val, boot.ci(res_ard)$bca[c(4,5)], boot.ci(res_rr)$bca[c(4,5)])
names(icumortality_val) <- c("precoce", "tardif", "precoce-tardif", "precoce/tardif", "ard_lci", "ard_uci", "rr_lci", "rr_uci")
icumortality_val

#### Ventilator free days
ventilfree_val <- list()
for (i in 1:imp_included$m) {
temp <- complete(imp_included, i)
ventilfree_val[[i]]<-c(
        mean(temp[temp$bras=="STRATEGIE PRECOCE",]$ventilfree),
        mean(temp[temp$bras=="STRATEGIE D ATTENTE",]$ventilfree),
        mean(temp[temp$bras=="STRATEGIE PRECOCE",]$ventilfree)-mean(temp[temp$bras=="STRATEGIE D ATTENTE",]$ventilfree)
)
res <- boot(temp, mdif, R=999, strata=as.numeric(imp1$bras=="STRATEGIE PRECOCE"))
ventilfree_val[[i]] <- c(ventilfree_val[[i]], var(res$t))
names(ventilfree_val[[i]]) <- c("precoce", "tardif", "precoce-tardif", "var")
}
ventilfree_val <- sapply(ventilfree_val, c)
vars <- ventilfree_val["var",]
ventilfree_val_pool <- apply(ventilfree_val, 1, mean)

rubinvar <- rubinr(thetas = ventilfree_val["precoce-tardif",], vars = ventilfree_val["var",])
ventilfree_val_pool <- c(ventilfree_val_pool, rubinvar, ventilfree_val_pool[3]+c(-1,1)*qnorm(.975)*sqrt(rubinvar))
names(ventilfree_val_pool)[c(5:7)] <- c("rubinvar", "lci", "uci")
ventilfree_val_pool

#### RRT free days
rrtfree_val <- list()
for (i in 1:imp_included$m) {
        temp <- complete(imp_included, i)
        rrtfree_val[[i]]<-c(
                mean(temp[temp$bras=="STRATEGIE PRECOCE",]$rrtfree),
                mean(temp[temp$bras=="STRATEGIE D ATTENTE",]$rrtfree),
                mean(temp[temp$bras=="STRATEGIE PRECOCE",]$rrtfree)-mean(temp[temp$bras=="STRATEGIE D ATTENTE",]$rrtfree)
        )
        res <- boot(temp, mdif, R=999, strata=as.numeric(imp1$bras=="STRATEGIE PRECOCE"))
        rrtfree_val[[i]] <- c(rrtfree_val[[i]], var(res$t))
        names(rrtfree_val[[i]]) <- c("precoce", "tardif", "precoce-tardif", "var")
}
rrtfree_val <- sapply(rrtfree_val, c)
vars <- rrtfree_val["var",]
rrtfree_val_pool <- apply(rrtfree_val, 1, mean)

rubinvar <- rubinr(thetas = rrtfree_val["precoce-tardif",], vars = rrtfree_val["var",])
rrtfree_val_pool <- c(rrtfree_val_pool, rubinvar, rrtfree_val_pool[3]+c(-1,1)*qnorm(.975)*sqrt(rubinvar))
names(rrtfree_val_pool)[c(5:7)] <- c("rubinvar", "lci", "uci")
rrtfree_val_pool

#### Vaso free days
vasofree_val <- list()
for (i in 1:imp_included$m) {
        temp <- complete(imp_included, i)
        vasofree_val[[i]]<-c(
                mean(temp[temp$bras=="STRATEGIE PRECOCE",]$rrtfree),
                mean(temp[temp$bras=="STRATEGIE D ATTENTE",]$rrtfree),
                mean(temp[temp$bras=="STRATEGIE PRECOCE",]$rrtfree)-mean(temp[temp$bras=="STRATEGIE D ATTENTE",]$rrtfree)
        )
        res <- boot(temp, mdif, R=999, strata=as.numeric(imp1$bras=="STRATEGIE PRECOCE"))
        vasofree_val[[i]] <- c(vasofree_val[[i]], var(res$t))
        names(vasofree_val[[i]]) <- c("precoce", "tardif", "precoce-tardif", "var")
}
vasofree_val <- sapply(vasofree_val, c)
vars <- vasofree_val["var",]
vasofree_val_pool <- apply(vasofree_val, 1, mean)

rubinvar <- rubinr(thetas = vasofree_val["precoce-tardif",], vars = vasofree_val["var",])
vasofree_val_pool <- c(vasofree_val_pool, rubinvar, vasofree_val_pool[3]+c(-1,1)*qnorm(.975)*sqrt(rubinvar))
names(vasofree_val_pool)[c(5:7)] <- c("rubinvar", "lci", "uci")
vasofree_val_pool
