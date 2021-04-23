load("explanatorydata.Rdata")
rm(list=ls()[ls()!="rcts"])

akikidbsg <- read.csv("~/Desktop/TEAM METHODS/phd/Trials Data/akiki data/SGakiki.csv", sep=";", comment.char="#")
newoutcomes_akiki<-data.frame(
        etude="akiki",
        id=substring(akikidbsg$subject_id, 5, nchar(akikidbsg$subject_id)-3),
        eerdep=akikidbsg$bil_depeerj28,
        saps2=NA,
        saps3=akikidbsg$adm_saps_tot
)

library("readxl")
pre_idealnewoutcomes<- read_excel("~/Desktop/TEAM METHODS/phd/Trials Data/ideal icu data/Ideal_ICU_FGrolleau bis.xlsx")
pre_idealnewoutcomes<-data.frame(pre_idealnewoutcomes)
pre_idealnewoutcomes<-pre_idealnewoutcomes[-c(1:2),]

newoutcomes_ideal<-data.frame(
        etude="idealicu",
        id=pre_idealnewoutcomes$"I_1",
        eerdep=pre_idealnewoutcomes$"Poursuite",
        saps2=as.numeric(pre_idealnewoutcomes$"IGSII..first"),
        saps3=NA
)

newoutcomes_ideal$eerdep[newoutcomes_ideal$eerdep=="Oui"]<-"OUI"
newoutcomes_ideal$eerdep[newoutcomes_ideal$eerdep=="Non"]<-"NON"
newoutcomes_ideal$eerdep[newoutcomes_ideal$eerdep=="."]<-NA

newoutcomes<-rbind(newoutcomes_akiki, newoutcomes_ideal)

###
rcts$d28mortality <- "initialization"
for (i in 1:nrow(rcts)){
if(with(rcts, derniere.nouvelles.censureJ60<=28 & etat.censureJ60==1)[i]==T){ 
        rcts$d28mortality[i] <- "Décédé"
        } else if(with(rcts, derniere.nouvelles.censureJ60<=28 & etat.censureJ60==0)[i]==T){
                rcts$d28mortality[i] <- "Perdu de vue"
                } else if(with(rcts, derniere.nouvelles.censureJ60>28)[i]==T) {
                        rcts$d28mortality[i] <- "Survivant"}
        }
###
newoutcomes<-merge(newoutcomes, rcts[,c("id","d28mortality")], by="id")

newoutcomes$surviesanseer<-"initialization"
for (i in 1:nrow(newoutcomes)) {
        if(newoutcomes$d28mortality[i]=="Décédé" & is.na(newoutcomes$d28mortality[i])==FALSE) { newoutcomes$surviesanseer[i] <- "NON" }
        else { 
                if(newoutcomes$eerdep[i]=="OUI" & is.na(newoutcomes$eerdep[i])==FALSE) { newoutcomes$surviesanseer[i] <- "NON" }
                else {
                        if(newoutcomes$eerdep[i]=="NON" & is.na(newoutcomes$eerdep[i])==FALSE) { newoutcomes$surviesanseer[i] <- "OUI" }                                
                }
        }
}
newoutcomes$surviesanseer[newoutcomes$surviesanseer=="initialization"]<-NA

rcts <- merge(rcts, newoutcomes[,c("id", "surviesanseer")], by="id")



