library(reshape)
library(tidyr)
library(cAIC4)
library(nlme)

###############################
## Reformat raw data for lme
###############################

x=read.csv('result4plot.csv')
xmelt=melt(x, id.vars =c("Species","Sample.ID"))
write.csv(xmelt,'Orobanchaceae_O2K_lme.csv')

#Manual adjustment format

###############################
## lme analysis
###############################

#Read in data

data=read.csv('Orobanchaceae_O2K_lme.csv')
View(data)
data$Life_history=as.factor(data$Life_history)

#Do some test regression to see if CIV needs to be log transformed
test1 <-  lme(value ~ Life_history, random = ~1|Sample_ID, data= subset(data, FCF_type == "CIV"))
plot(test1)
qqnorm(resid(test1))
qqline(resid(test1))

#log transform CIV to better meet normal residuals ####
test2 <-  lme(value ~ Life_history, random = ~1|Sample_ID, data= subset(data, FCF_type == "CIV"))
plot(test2)
qqnorm(resid(test2))
qqline(resid(test2))


data=read.csv('Orobanchaceae_O2K_lme_logCIV.csv')
#All FCF pass the Fligner-Killeen test for homogeneity of group variances 
fligner.test(value~Life_history ,data= subset(data, FCF_type == "CI"))
fligner.test(value~Life_history ,data= subset(data, FCF_type == "CII"))
fligner.test(value~Life_history ,data= subset(data, FCF_type == "CIV"))
fligner.test(value~Life_history ,data= subset(data, FCF_type == "AOX"))
fligner.test(value~Life_history ,data= subset(data, FCF_type == "DHex"))
fligner.test(value~Life_history ,data= subset(data, FCF_type == "OXPHOS_efficiency"))

# functions to fit models
fit_mods = function(FCF){
  m1 <- lme(value ~ Life_history, random = ~1|Species, data = subset(data, FCF_type == FCF))
  m2 <- lme(value ~ Life_history, random = ~1|Field2lab_time, data = subset(data, FCF_type == FCF))
  m3 <- lme(value ~ Life_history, random = ~1|Species/Field2lab_time, data = subset(data, FCF_type == FCF))
  
  p1 = cAIC(m1)$caic
  p2 = cAIC(m2)$caic
  p3 = cAIC(m3)$caic
  
  return(paste("m1 =", round(p1,2), ", m2 =", round(p2,2),", m3 =", round(p3,2)))
}


data %>% 
  group_by(FCF_type) %>% 
  nest(data = c(Species, Sample_ID, Life_history, value)) %>% 
  mutate(fit_mods) %>% 
  select(-data) %>% 
  unnest()
  
#Nest different populations within the species does not improve the model

#### FCF by Life history Models ####

aox.mod <- lme(value ~ Life_history, random = ~1|Species, data= subset(data, FCF_type == "AOX"))
#check model diagnostics
summary(aox.mod)
plot(aox.mod)
qqnorm(resid(aox.mod))
qqline(resid(aox.mod))
#p-value=0.23

dhex.mod <- lme(value ~ Life_history, random = ~1|Species, data= subset(data, FCF_type == "DHex"))
#check model diagnostics
summary(dhex.mod)
plot(dhex.mod)
qqnorm(resid(dhex.mod))
qqline(resid(dhex.mod))
#p-value=0.35

CI.mod <- lme(value ~ Life_history, random = ~1|Species, data= subset(data, FCF_type == "CI"))
#check model diagnostics
summary(CI.mod)
plot(CI.mod)
qqnorm(resid(CI.mod))
qqline(resid(CI.mod))
#p-value=0.89

CII.mod <- lme(value ~ Life_history, random = ~1|Species, data= subset(data, FCF_type == "CII"))
#check model diagnostics
plot(CII.mod)
summary(CII.mod)
qqnorm(resid(CII.mod))
qqline(resid(CII.mod))
#p=0.02

CIV.mod <- lme(value ~ Life_history, random = ~1|Species, data= subset(data, FCF_type == "CIV"))
#check model diagnostics
summary(CIV.mod)
plot(CIV.mod) 
qqnorm(resid(CIV.mod)) 
qqline(resid(CIV.mod))
#p=0.25

OXPHOS.mod <- lme(value ~ Life_history, random = ~1|Species, data= subset(data, FCF_type == "OXPHOS_efficiency"))
#check model diagnostics
plot(OXPHOS.mod)
qqnorm(resid(OXPHOS.mod))
qqline(resid(OXPHOS.mod))
summary(OXPHOS.mod)
#p=0.6124

######################################
#Correlation between FCFs
x$CII.flux.control..1.C.D.=as.numeric(x$CII.flux.control..1.C.D.)
x$Excess.capacity.of.CIV..G.F.1.=as.numeric(x$Excess.capacity.of.CIV..G.F.1.)
x$OXPHOS.coupling.efficiency..1.A.B.=as.numeric(x$OXPHOS.coupling.efficiency..1.A.B.)
x$AOX.flux.control.MAX..1.preASC.prenPG.=as.numeric(x$AOX.flux.control.MAX..1.preASC.prenPG.)
x$CI.flux.control.MAX..1.pre.nPG.D.=as.numeric(x$CI.flux.control.MAX..1.pre.nPG.D.)
x$DHex.flux.control..1.B.C.=as.numeric(x$DHex.flux.control..1.B.C.)

a=x[x$habit=='holo',]
b=x[x$habit=='hemi',]

maxres.mod <- lme(MAX/protein ~ habit, random = ~1|Species, data= x,na.action = na.omit)

Linear mixed-effects model fit by REML
  Data: x 
       AIC     BIC   logLik
  544.1381 550.001 -268.069

Random effects:
 Formula: ~1 | Species
        (Intercept) Residual
StdDev:    813.1369  862.834

Fixed effects:  MAX/protein ~ habit 
               Value Std.Error DF  t-value p-value
(Intercept) 1732.153  447.1897 28 3.873420  0.0006
habitholo   1135.848  770.7791  4 1.473636  0.2146
 Correlation: 
          (Intr)
habitholo -0.58 

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-1.82843657 -0.42434388 -0.05276825  0.39373463  1.59846807 

Number of Observations: 34
Number of Groups: 6 


AOX.CII.global <-lme(AOX.flux.control.MAX..1.preASC.prenPG. ~CII.flux.control..1.C.D., random = ~1|Species, data= x,na.action = na.omit)
summary(AOX.CII.global)
Linear mixed-effects model fit by REML
  Data: x 
        AIC       BIC   logLik
  -20.95396 -13.90915 14.47698

Random effects:
 Formula: ~1 | Species
        (Intercept)  Residual
StdDev:   0.1491936 0.1379706

Fixed effects:  AOX.flux.control.MAX..1.preASC.prenPG. ~ CII.flux.control..1.C.D. 
                             Value  Std.Error DF  t-value p-value
(Intercept)              0.2528502 0.06077173 35 4.160654   2e-04
CII.flux.control..1.C.D. 0.1171381 0.11133618 35 1.052112   3e-01
 Correlation: 
                         (Intr)
CII.flux.control..1.C.D. -0.456

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-1.97089122 -0.55178825 -0.08147687  0.65927514  2.15234012 

Number of Observations: 45
Number of Groups: 9 

AOX.CII.holo <-lme(AOX.flux.control.MAX..1.preASC.prenPG. ~CII.flux.control..1.C.D., random = ~1|Species, data= a,na.action = na.omit)
AOX.CII.hemi <-lme(AOX.flux.control.MAX..1.preASC.prenPG. ~CII.flux.control..1.C.D., random = ~1|Species, data= b,na.action = na.omit)


DEX.CII.global<-lme(DHex.flux.control..1.B.C.~CII.flux.control..1.C.D., random = ~1|Species, data= x,na.action = na.omit)
summary(DEX.CII.global)
Linear mixed-effects model fit by REML
  Data: x 
       AIC       BIC  logLik
  -35.3692 -28.14255 21.6846

Random effects:
 Formula: ~1 | Species
        (Intercept)  Residual
StdDev:   0.1347996 0.1189129

Fixed effects:  DHex.flux.control..1.B.C. ~ CII.flux.control..1.C.D. 
                             Value  Std.Error DF  t-value p-value
(Intercept)              0.4532918 0.05320405 37 8.519874  0.0000
CII.flux.control..1.C.D. 0.2414967 0.08234610 37 2.932704  0.0057
 Correlation: 
                         (Intr)
CII.flux.control..1.C.D. -0.413

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-2.32905508 -0.51546294  0.09331725  0.59541506  1.57344511 

Number of Observations: 47
Number of Groups: 9 

AOX.DEX.global<-lme(AOX.flux.control.MAX..1.preASC.prenPG. ~ DHex.flux.control..1.B.C., random = ~1|Species, data= x,na.action = na.omit)
summary(AOX.DEX.global)
Linear mixed-effects model fit by REML
  Data: x 
        AIC       BIC   logLik
  -21.76913 -14.72433 14.88456

Random effects:
 Formula: ~1 | Species
        (Intercept)  Residual
StdDev:   0.1516199 0.1375128

Fixed effects:  AOX.flux.control.MAX..1.preASC.prenPG. ~ DHex.flux.control..1.B.C. 
                              Value Std.Error DF  t-value p-value
(Intercept)               0.1913710 0.1006342 35 1.901650  0.0655
DHex.flux.control..1.B.C. 0.1771373 0.1649733 35 1.073733  0.2903
 Correlation: 
                          (Intr)
DHex.flux.control..1.B.C. -0.839

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-2.18822165 -0.61657176 -0.02999418  0.46890314  2.44821417 

Number of Observations: 45
Number of Groups: 9 


AOX.DEX.helo<-lme(AOX.flux.control.MAX..1.preASC.prenPG. ~ DHex.flux.control..1.B.C., random = ~1|Species, data= a,na.action = na.omit)




CIV.CII.global <-lme(Excess.capacity.of.CIV..G.F.1.~ CII.flux.control..1.C.D.,random = ~1|Species, data= x,na.action = na.omit)
Linear mixed-effects model fit by REML
  Data: x 
       AIC      BIC    logLik
  173.2469 180.2917 -82.62346

Random effects:
 Formula: ~1 | Species
        (Intercept) Residual
StdDev:    1.389506 1.326159

Fixed effects:  Excess.capacity.of.CIV..G.F.1. ~ CII.flux.control..1.C.D. 
                             Value Std.Error DF   t-value p-value
(Intercept)               2.615088 0.5715167 35  4.575699  0.0001
CII.flux.control..1.C.D. -1.444625 1.0664150 35 -1.354656  0.1842
 Correlation: 
                         (Intr)
CII.flux.control..1.C.D. -0.464

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-1.43383610 -0.43424227 -0.09628033  0.11514910  3.19969987 

Number of Observations: 45
Number of Groups: 9 

OXPHOS.CII.global <-lme(OXPHOS.coupling.efficiency..1.A.B. ~ CII.flux.control..1.C.D.,random = ~1|Species, data= x,na.action = na.omit)
summary(OXPHOS.CII.global)
Linear mixed-effects model fit by REML
  Data: x 
        AIC       BIC   logLik
  -48.97186 -41.74521 28.48593

Random effects:
 Formula: ~1 | Species
        (Intercept)  Residual
StdDev:   0.1273256 0.1005039

Fixed effects:  OXPHOS.coupling.efficiency..1.A.B. ~ CII.flux.control..1.C.D. 
                              Value  Std.Error DF   t-value p-value
(Intercept)               0.3680404 0.04886823 37  7.531281  0.0000
CII.flux.control..1.C.D. -0.0304394 0.07015931 37 -0.433861  0.6669
 Correlation: 
                         (Intr)
CII.flux.control..1.C.D. -0.383

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-2.06272558 -0.49696598  0.05078506  0.64025323  1.76593329 

Number of Observations: 47
Number of Groups: 9 

OXPHOS.AOX.global <-lme(OXPHOS.coupling.efficiency..1.A.B. ~ AOX.flux.control.MAX..1.preASC.prenPG. ,random = ~1|Species, data= x,na.action = na.omit)
Linear mixed-effects model fit by REML
  Data: x 
        AIC       BIC   logLik
  -49.38623 -42.34143 28.69311

Random effects:
 Formula: ~1 | Species
        (Intercept)   Residual
StdDev:   0.1298404 0.09640437

Fixed effects:  OXPHOS.coupling.efficiency..1.A.B. ~ AOX.flux.control.MAX..1.preASC.prenPG. 
                                            Value  Std.Error DF   t-value
(Intercept)                             0.3749747 0.05502334 35  6.814829
AOX.flux.control.MAX..1.preASC.prenPG. -0.0362796 0.10822064 35 -0.335237
                                       p-value
(Intercept)                             0.0000
AOX.flux.control.MAX..1.preASC.prenPG.  0.7394
 Correlation: 
                                       (Intr)
AOX.flux.control.MAX..1.preASC.prenPG. -0.555

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-2.15417819 -0.46181593  0.09213853  0.62793016  1.77896057 

Number of Observations: 45
Number of Groups: 9 

CI.CII.global <-lme(CI.flux.control.MAX..1.pre.nPG.D. ~ CII.flux.control..1.C.D.,random = ~1|Species, data= x,na.action = na.omit)
Linear mixed-effects model fit by REML
  Data: x 
        AIC       BIC   logLik
  -58.97716 -51.93236 33.48858

Random effects:
 Formula: ~1 | Species
        (Intercept)  Residual
StdDev:  0.01553153 0.1040365

Fixed effects:  CI.flux.control.MAX..1.pre.nPG.D. ~ CII.flux.control..1.C.D. 
                               Value  Std.Error DF   t-value p-value
(Intercept)               0.19871975 0.02239485 35  8.873456  0.0000
CII.flux.control..1.C.D. -0.07147155 0.06402001 35 -1.116394  0.2719
 Correlation: 
                         (Intr)
CII.flux.control..1.C.D. -0.681

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-2.20635344 -0.71472630 -0.09179925  0.69632265  2.20841023 

Number of Observations: 45
Number of Groups: 9 

CI.CII.habit <-lme(CI.flux.control.MAX..1.pre.nPG.D. ~ CII.flux.control..1.C.D.*habit,random = ~1|Species, data= x,na.action = na.omit)
summary(CI.CII.habit)
Linear mixed-effects model fit by REML
  Data: x 
        AIC       BIC   logLik
  -49.80431 -39.52288 30.90216

Random effects:
 Formula: ~1 | Species
        (Intercept) Residual
StdDev:  0.03193243 0.102237

Fixed effects:  CI.flux.control.MAX..1.pre.nPG.D. ~ CII.flux.control..1.C.D. *      habit 
                                         Value  Std.Error DF   t-value
(Intercept)                         0.20296483 0.02953592 34  6.871797
CII.flux.control..1.C.D.           -0.18601012 0.10816490 34 -1.719690
habitholo                          -0.00204541 0.05458627  7 -0.037471
CII.flux.control..1.C.D.:habitholo  0.14582257 0.14812247 34  0.984473
                                   p-value
(Intercept)                         0.0000
CII.flux.control..1.C.D.            0.0946
habitholo                           0.9712
CII.flux.control..1.C.D.:habitholo  0.3318
 Correlation: 
                                   (Intr) CII.f...1.C.D. habthl
CII.flux.control..1.C.D.           -0.507                      
habitholo                          -0.541  0.274               
CII.flux.control..1.C.D.:habitholo  0.370 -0.730         -0.659

Standardized Within-Group Residuals:
       Min         Q1        Med         Q3        Max 
-2.4257437 -0.6339922 -0.1169118  0.5165227  2.0646125 

Number of Observations: 45
Number of Groups: 9 


MAX.CII.global <-lme(MAX/protein ~ CII.flux.control..1.C.D.,random = ~1|Species, data= x,na.action = na.omit)
summary(MAX.CII.global)
Linear mixed-effects model fit by REML
  Data: x 
       AIC      BIC    logLik
  528.0222 533.7581 -260.0111

Random effects:
 Formula: ~1 | Species
        (Intercept) Residual
StdDev:    1265.927 812.0569

Fixed effects:  MAX/protein ~ CII.flux.control..1.C.D. 
                             Value Std.Error DF   t-value p-value
(Intercept)               2439.331  567.4127 26  4.299043  0.0002
CII.flux.control..1.C.D. -1374.068  838.1618 26 -1.639382  0.1132
 Correlation: 
                         (Intr)
CII.flux.control..1.C.D. -0.326

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-2.06127741 -0.42140035 -0.05119637  0.35786062  1.86309225 

Number of Observations: 33
Number of Groups: 6 

MAX.AOX.global <-lme(MAX/protein ~ AOX.flux.control.MAX..1.preASC.prenPG. ,random = ~1|Species, data= x,na.action = na.omit)
summary(MAX.AOX.global)
Linear mixed-effects model fit by REML
  Data: x 
       AIC      BIC    logLik
  527.5922 533.3281 -259.7961

Random effects:
 Formula: ~1 | Species
        (Intercept) Residual
StdDev:    916.5727 861.1807

Fixed effects:  MAX/protein ~ AOX.flux.control.MAX..1.preASC.prenPG. 
                                           Value Std.Error DF   t-value
(Intercept)                             2538.860  516.9344 26  4.911379
AOX.flux.control.MAX..1.preASC.prenPG. -1542.016 1229.3527 26 -1.254332
                                       p-value
(Intercept)                             0.0000
AOX.flux.control.MAX..1.preASC.prenPG.  0.2209
 Correlation: 
                                       (Intr)
AOX.flux.control.MAX..1.preASC.prenPG. -0.624

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-1.95739398 -0.45791610 -0.02073471  0.23963445  1.69283101 

Number of Observations: 33
Number of Groups: 6 

MAX.DEX.global <-lme(MAX/protein ~ DHex.flux.control..1.B.C. ,random = ~1|Species, data= x,na.action = na.omit)
> summary(MAX.DEX.global)
Linear mixed-effects model fit by REML
  Data: x 
       AIC      BIC    logLik
  526.6598 532.3957 -259.3299

Random effects:
 Formula: ~1 | Species
        (Intercept) Residual
StdDev:    537.8382 914.3203

Fixed effects:  MAX/protein ~ DHex.flux.control..1.B.C. 
                              Value Std.Error DF   t-value p-value
(Intercept)                359.8225  840.4802 26 0.4281154  0.6721
DHex.flux.control..1.B.C. 3066.3895 1378.0226 26 2.2252098  0.0350
 Correlation: 
                          (Intr)
DHex.flux.control..1.B.C. -0.946

Standardized Within-Group Residuals:
       Min         Q1        Med         Q3        Max 
-1.7444162 -0.4572101 -0.1305728  0.2594659  2.2149007 

Number of Observations: 33
Number of Groups: 6 


max.habit<-lme(MAX/protein ~habit, ,random = ~1|Species, data= x,na.action = na.omit)
summary(max.habit)
Linear mixed-effects model fit by REML
  Data: x 
       AIC     BIC   logLik
  544.1381 550.001 -268.069

Random effects:
 Formula: ~1 | Species
        (Intercept) Residual
StdDev:    813.1369  862.834

Fixed effects:  MAX/protein ~ habit 
               Value Std.Error DF  t-value p-value
(Intercept) 1732.153  447.1897 28 3.873420  0.0006
habitholo   1135.848  770.7791  4 1.473636  0.2146
 Correlation: 
          (Intr)
habitholo -0.58 

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-1.82843657 -0.42434388 -0.05276825  0.39373463  1.59846807 

Number of Observations: 34
Number of Groups: 6 


###################
pdf('DexCII.pdf',width = 4,height = 3)
ggplot(x, aes(x = CII.flux.control..1.C.D., y = DHex.flux.control..1.B.C., color = habit)) +
	geom_point() +                                # Scatter plot
	geom_smooth(method = "lm", se = TRUE) +       
	labs(x= "CII", y = "Dex") +
	theme_minimal()
dev.off()


pdf('AOXCII.pdf',width = 4,height = 3)
ggplot(x, aes(x = CII.flux.control..1.C.D., y = AOX.flux.control.MAX..1.preASC.prenPG., color = habit)) +
	geom_point() +                                # Scatter plot
	geom_smooth(method = "lm", se = TRUE) +       
	labs(x= "CII", y = "AOX") +
	theme_minimal()
dev.off()

pdf('AOXDex.pdf',width = 4,height = 3)
ggplot(x, aes(x =DHex.flux.control..1.B.C., y = AOX.flux.control.MAX..1.preASC.prenPG., color = habit)) +
	geom_point() +                                # Scatter plot
	geom_smooth(method = "lm", se = TRUE) +       
	labs(x= "Dex", y = "AOX") +
	theme_minimal()
dev.off()

pdf('CICII.pdf',width = 4,height = 3)
ggplot(x, aes(x =CII.flux.control..1.C.D., y = CI.flux.control.MAX..1.pre.nPG.D., color = habit)) +
	geom_point() +                                # Scatter plot
	geom_smooth(method = "lm", se = TRUE) +       
	labs(x= "CII", y = "CI") +
	theme_minimal()
dev.off()

pdf('OXPHOSCII.pdf',width = 4,height = 3)
ggplot(x, aes(x =CII.flux.control..1.C.D., y = OXPHOS.coupling.efficiency..1.A.B., color = habit)) +
	geom_point() +                                # Scatter plot
	geom_smooth(method = "lm", se = TRUE) +       
	labs(x= "CII", y = "OXPHOS") +
	theme_minimal()
dev.off()

pdf('CIVCII.pdf',width = 4,height = 3)
ggplot(x, aes(x =CII.flux.control..1.C.D., y = Excess.capacity.of.CIV..G.F.1., color = habit)) +
	geom_point() +                                # Scatter plot
	geom_smooth(method = "lm", se = TRUE) +       
	labs(x= "CII", y = "CIV") +
	theme_minimal()
dev.off()

pdf('CIVCI.pdf',width = 4,height = 3)
ggplot(x, aes(x =log(Excess.capacity.of.CIV..G.F.1.), y = CI.flux.control.MAX..1.pre.nPG.D., color = habit)) +
	geom_point() +                                # Scatter plot
	geom_smooth(method = "lm", se = TRUE) +       
	labs(x= "CI", y = "CIV") +
	theme_minimal()
dev.off()


pdf('MaxCII.pdf',width = 4,height = 3)
ggplot(x, aes(x =CII.flux.control..1.C.D., y = MAX/protein, color = habit)) +
	geom_point() +                                # Scatter plot
	geom_smooth(method = "lm", se = TRUE) +       
	labs(x= "CII", y = "max resp") +
	theme_minimal()
dev.off()



pdf('MaxAOX.pdf',width = 4,height = 3)
ggplot(x, aes(x =AOX.flux.control.MAX..1.preASC.prenPG., y = MAX/protein, color = habit)) +
	geom_point() +                                # Scatter plot
	geom_smooth(method = "lm", se = TRUE) +       
	labs(x= "AOX", y = "max resp") +
	theme_minimal()
dev.off()


pdf('MaxDex.pdf',width = 4,height = 3)
ggplot(x, aes(x =DHex.flux.control..1.B.C., y = MAX/protein, color = habit)) +
	geom_point() +                                # Scatter plot
	geom_smooth(method = "lm", se = TRUE) +       
	labs(x= "Dex", y = "max resp") +
	theme_minimal()
dev.off()