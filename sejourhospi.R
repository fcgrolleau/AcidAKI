### Recuperer les dates de sorties de l'hospital dans les deux bases

library("readxl")
prep_ideal_hospi_disch<- read_excel("~/Desktop/TEAM METHODS/phd/Trials Data/ideal icu data/Ideal_ICU_FGrolleau.xlsx")
prep_ideal_hospi_disch<-data.frame(prep_ideal_hospi_disch)
ideal_hospi_disch<-data.frame(
        etude="idealicu",
        id=prep_ideal_hospi_disch$"I_1",
        datesortiehosp=prep_ideal_hospi_disch$BSHOP_1
)

akikidbsg <- read.csv("~/Desktop/TEAM METHODS/phd/Trials Data/akiki data/SGakiki.csv", sep=";", comment.char="#")
akiki_hospi_disch<-data.frame(
        etude="akiki",
        id=substring(akikidbsg$subject_id, 5, nchar(akikidbsg$subject_id)-3),
        datesortiehosp=as.POSIXct(akikidbsg$bil_sortiedate, format="%d/%m/%y"))

premerge <- rbind(ideal_hospi_disch, akiki_hospi_disch)

### Selectionner uniquement les patients inclus dans acidAKI

acidsortie <- premerge[premerge$id %in% imp1$id,]

# Merge avec imp1 (impossible d'imputer des dates donc procedure d'imputation non réalisée à nouveau)
acidsortie <- merge(imp1, acidsortie, by="id")

acidsortie$los_h <- as.numeric(with(acidsortie, difftime(datesortiehosp, date.rando, unit="days")))
hospstay <- tapply(acidsortie$los_h, acidsortie$bras, function(x) mean(x, na.rm=T))
hospstay_sd <- tapply(acidsortie$los_h, acidsortie$bras, function(x) sd(x, na.rm=T))
hospstay_pool <- hospstay["STRATEGIE PRECOCE"]-hospstay["STRATEGIE D ATTENTE"]

library(boot)
res <- boot(acidsortie[!is.na(acidsortie$los_h),], mdif, cj="los_h", R=999)

hospstay_pool <- c(hospstay["STRATEGIE PRECOCE"], hospstay_sd["STRATEGIE PRECOCE"], hospstay["STRATEGIE D ATTENTE"], hospstay_sd["STRATEGIE D ATTENTE"], hospstay_pool, boot.ci(res)$bca[c(4,5)])
names(hospstay_pool) <- c("precoce","precoce_sd", "tardif","tardif_sd", "md", "lci", "uci")

hospstay_exp <- round(hospstay_pool, 2)

hospstay_exp <- c(paste0(format(hospstay_exp[1], nsmall = 2), " (", format(hospstay_exp[2], nsmall = 2), ")"),
                 paste0(format(hospstay_exp[3], nsmall = 2), " (", format(hospstay_exp[4], nsmall = 2), ")"),
                 format(hospstay_exp[5], nsmall = 2),
                 paste0("(", format(hospstay_exp[6], nsmall = 2), " to ", format(hospstay_exp[7], nsmall = 2), ")")
)

#write.csv(t(hospstay_exp), "hospstay.csv")


