ard_km <- function(d, i=1:nrow(d)) {
        z<-d[i,]
        survfit.boot<-survfit(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~bras, data = z)
        return(
                1-tail(weightedsurvivaltables(survfit.boot$time, survfit.boot$surv, group="B"),1)[2]-
                (1-tail(weightedsurvivaltables(survfit.boot$time, survfit.boot$surv, group="A"),1)[2])
        )
}

rr_km <- function(d, i=1:nrow(d)) {
        z<-d[i,]
        survfit.boot<-survfit(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~bras, data = z)
        return(
                (1-tail(weightedsurvivaltables(survfit.boot$time, survfit.boot$surv, group="B"),1)[2])/
                        (1-tail(weightedsurvivaltables(survfit.boot$time, survfit.boot$surv, group="A"),1)[2])
        )
}


ardbootfunction<-function(dataset, nboot=1000) {
        
        ardmat<-matrix(999, nrow=nboot, ncol = length(seq_ph))
        
        probs.int.p <- with(dataset, ifelse(bras=="STRATEGIE PRECOCE", seq_ph-7.300908, 0))
        probs.c<-dataset$ph_e_p-7.300908
        rg <-range(probs.c)
        kn<- -0.0009076725 
        
        seq_ph <- modelr::seq_range(rcts_imp1$ph_e_p, n=prec)
        
        set.seed(1000)
        for (i in 1:nboot) {
                bootindex<-c(sample(which(dataset$bras!="STRATEGIE PRECOCE"), replace=TRUE), sample(which(dataset$bras=="STRATEGIE PRECOCE"), replace=TRUE))
                zdata<-dataset[bootindex,]
                zdata$probs.int.p <- with(zdata, ifelse(bras=="STRATEGIE PRECOCE", ph_e_p-7.300908, 0))
                zdata$probs.c<-zdata$ph_e_p-7.300908
                
                zdata$a1<-ns(zdata$probs.c, knots = kn, Boundary.knots=rg)[,1]
                zdata$a2<-ns(zdata$probs.c, knots = kn, Boundary.knots=rg)[,2]

                zdata$b1<-ns(zdata$probs.int.p, knots = kn, Boundary.knots=rg)[,1]
                zdata$b2<-ns(zdata$probs.int.p, knots = kn, Boundary.knots=rg)[,2]

                
                fit.int.ard<-coxph(Surv(derniere.nouvelles.censureJ60, etat.censureJ60)~I(bras=='STRATEGIE PRECOCE') + a1+a2 + b1+b2, data = zdata, x=TRUE)
                
                y0.precoce_<-data.frame(bras = "STRATEGIE PRECOCE", ns(seq_ph-7.300908, knots = kn, Boundary.knots=rg), ns(seq_ph, knots = kn, Boundary.knots=rg))
                names(y0.precoce_) <- c("bras","a1", "a2", "b1", "b2")
                
                y0.tardif_<-data.frame(bras = "STRATEGIE D ATTENTE", ns(seq_ph-7.300908, knots = kn, Boundary.knots=rg), ns(0, knots = kn, Boundary.knots=rg))
                names(y0.tardif_) <- c("bras","a1", "a2", "b1", "b2")
                
                temp_precoced60d<-1-tail(survfit(fit.int.ard, newdata=y0.precoce_)$surv,1)
                
                temp_tardifd60d<-1-tail(survfit(fit.int.ard, newdata=y0.tardif_)$surv,1)
                ####
                
                ardmat[i,]<-temp_precoced60d-temp_tardifd60d
        }
        apply(ardmat, 2, function(x) sd(x))
}

mdif <- function(d, i=1:nrow(d), cj="ventilfree") {
        z <- d[i,]
        return(mean(z[z$bras=="STRATEGIE PRECOCE",cj])-mean(z[z$bras=="STRATEGIE D ATTENTE",cj]))
}

binout_ard_boot <- function(d, i=1:nrow(d), cj="hospitalmortality", event="Décédé") {
        z<-d[i,]
        ARD_val<- sum(z[z$bras=="STRATEGIE PRECOCE",cj]==event, na.rm=T)/length(z[z$bras=="STRATEGIE PRECOCE",cj])
        -
                sum(z[z$bras=="STRATEGIE D ATTENTE",cj]==event, na.rm=T)/length(z[z$bras=="STRATEGIE D ATTENTE",cj])
        return(ARD_val)
}

binout_rr_boot <- function(d, i=1:nrow(d), cj="hospitalmortality", event="Décédé") {
        z<-d[i,]
        RR_val<- (sum(z[z$bras=="STRATEGIE PRECOCE",cj]==event, na.rm=T)/length(z[z$bras=="STRATEGIE PRECOCE",cj]))/(sum(z[z$bras=="STRATEGIE D ATTENTE",cj]==event, na.rm=T)/length(z[z$bras=="STRATEGIE D ATTENTE",cj]))
        return(RR_val)
}

rubinr <-  function(thetas, vars) {
        theta <- mean(thetas)
        w <- mean(vars)
        b <- var(thetas)
        return(w+(1+1/length(vars))*b)
}
