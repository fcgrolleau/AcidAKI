### RMST in weighted populations


weightedsurvivaltables<-function(timeobj, survobj, group="A"){
        for (i in 1:length(timeobj)) {
                if(timeobj[i+1]<timeobj[i]) break
                nbtp<-i+1
        }
        groupA<-cbind(timeobj, survobj)[1:nbtp,]
        groupB<-cbind(timeobj, survobj)[as.numeric(nbtp+1):as.numeric(length(timeobj)),]
        
        if(group=="B") { toprint<-groupB }
        else { toprint<-groupA }
        toprint
}


rmst<-function(timeobj, survobj, group) {
x<-weightedsurvivaltables(timeobj, survobj, group)        
subarea<-vector()
for (i in 1:(nrow(x))-1) {
        subarea[i]<-(x[i+1,1]-x[i,1])*x[i+1,2]
        rmstval<-sum(subarea)
}
rmstval
}

### example
# rmst(survfit.obj$time, survfit.obj$surv, group="A")-rmst(survfit.obj$time, survfit.obj$surv, group="B")


