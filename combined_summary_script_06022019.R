## ----------------------------- Microarray, all CYPs:
## 04/29/2016
## ------------------ 
library(MASS)
library(glmnet)
library(multcomp) ## multiT, multinormal
library(sgof)	## Beta-binomial SGoF procedure, for correcting positive correlation

## caculation of the numerator of the sobol main index
t.gs <- function(i, coef, mse, rhom){
  rho <- rhom[i,-i]
  return ( ( coef[i]+ matrix(coef[-i],nr=1) %*% ( matrix(mse[-i],nc=1) * matrix(rho,nc=1) ) / mse[i] )^2 * mse[i]^2  )
}

scan.one.CYP <- function(CYP, input.genes){
  ## fit the sample with GLM via Iteratively Reweighted Least Square (IRLS)
  fit <- glm(CYP~., data=data.frame(input.genes))
  glm.coef <- fit$coef
  glm.p <- summary(fit)$coef[,4]; 
  glm.p.order <- order(glm.p[-1], na.last=TRUE)
  
  ## calculate emprical estimate of sobol ranking using glm.coef
  coef.est <- glm.coef[-1]
  mse.est <- apply(input.genes, 2, function(x){var(x, na.rm=T)^0.5})
  rhom.est <- cor(input.genes, use="pairwise.complete.obs")
  
  est.gs <- apply(matrix(1:length(coef.est), nc=1), 1, function(i){t.gs(i, coef.est, mse.est, rhom.est)})  
  est.gs.order <- order(est.gs, decreasing=TRUE, na.last=TRUE)
  
  my.list <- list(glm.coef.with.intercept=glm.coef, glm.p.with.intercept=glm.p, input.se.est=mse.est, input.rho.est=rhom.est, est.gs=est.gs, est.gs.order=est.gs.order, glm.p.order=glm.p.order)
  
  return(my.list)
}

## ----------------------------- code for scanning for subset effect using polynomial with degree k=3:
SI <- function(y, input){
  fit.1 <- glm(y~., data=data.frame(input));	glm.coef <- fit.1$coef;
  SI <- var( apply(input, 1, function(x){matrix(x, nr=1) %*% matrix(glm.coef[-1], nc=1)}), na.rm=TRUE )
  ##return(SI)
  return(c(SI, fit.1$deviance / fit.1$null.deviance, fit.1$null.deviance))
}

## k >=3
comb.m <- function(input, k=3){
  d <- dim(input)[2]
  if(d>=3){
    combn.index <- apply(matrix(2:(min(d,k)), nr=1),2, function(x){combn(c(1:d), x)} ) 
    
    combn.list <- lapply(combn.index, function(x){cc<-input[,x[1,]]; for(j in 2:dim(x)[1]){cc <- cc*input[,x[j,]]}; return(list(cc))})
    combn.input <- do.call(cbind, lapply(combn.list, function(x){x[[1]]}))
  }
  if(d==2){combn.input <- input[,1]*input[,2]}
  return(combn.input)
}

SI.comb <- function(y, input, k=3){
  main.new.m <- apply(input, 2, function(x){ apply(matrix(1:k),1, function(i){x^i}) })
  main.SI <- apply(main.new.m, 2, function(x){SI(y,matrix(x, nr=dim(input)[1]))})
  
  comb.new.m <- cbind(matrix(main.new.m, nr=dim(input)[1]), comb.m(input, k))
  comb.SI <- SI(y, comb.new.m)
  
  ##return(c(comb.SI, main.SI, comb.SI-sum(main.SI)))
  return(cbind(comb.SI, main.SI))
}

SI.main.total <- function(i, y, input, k=3){
  main.input <- apply(matrix(c(1:k), nc=1), 1, function(j){input[,i]^j})
  main.SI <- SI(y, main.input)
  
  total.input.1 <- apply(input[,-i], 2, function(x){ apply(matrix(1:k),1, function(i){x^i}) })
  total.input <- cbind(matrix(total.input.1, nr=dim(input)[1]), comb.m(input[,-i], k))
  total.SI <- SI(y, input[,-i])
  return(cbind(total.SI, main.SI))
}

test.sig.against.noise <- function(CYP, input.genes, sig.gene.index){
  sig.genes <- input.genes[,sig.gene.index]
  noise.mu.est <- quantile(apply(sig.genes,2,function(x){mean(x,na.rm=T)}),c(1:10)/10)[c(1,9)]
  noise.se.est <- quantile(apply(sig.genes,2,function(x){var(x,na.rm=T)^0.5}),c(1:10)/10)[c(1,9)]
  noise.rho.est <- quantile( cor(sig.genes, use="pairwise.complete.obs"),c(1:10)/10)[c(1,9)]
  
  d <- min(length(sig.gene.index),30)
  noise.rho <- matrix(runif(1, min=0, max=noise.rho.est[2]),nr=d,nc=d)
  noise.rho[lower.tri(noise.rho)] <- t(noise.rho)[lower.tri(noise.rho)]
  noise.rho[row(noise.rho)==col(noise.rho)] <- 1
  noise.se <- matrix(runif(d, min=noise.se.est[1], max=noise.se.est[2]),nc=1)
  noise.cov <- noise.se %*% t(noise.se) * noise.rho
  noise.mu <- matrix(runif(d, min=noise.mu.est[1], max=noise.se.est[2]),nc=1)
  
  noise <- mvrnorm(length(CYP), noise.mu, noise.cov)
  new.input <- cbind(sig.genes, noise)
  return( list( scan.one.CYP(CYP, new.input), length(sig.gene.index) , names(input.genes)[sig.gene.index] ))
}

##----------------------------------------------------------------------------------------------------
##-------------------- code for variable selection via deviance estimate from IRLS -------------------

## Null Deviance = 2(LL(Saturated Model) - LL(Null Model)) on df = df_Sat - df_Null
## Residual Deviance = 2(LL(Saturated Model) - LL(Proposed Model)) df = df_Sat - df_Res

interaction3gene <- function(y,gene){
  xm <- data.frame(gene)
  m1 <- glm(y ~ .,data=xm); ##print(summary(m1));
  m2 <- glm(y ~ .^4, data=xm);##print(summary(m2));
  return(c(m1$deviance, m2$deviance, (m1$deviance-m2$deviance)/m1$deviance))
  
}

##--------------------------------------- combine data files ---------------------------------
data <- read.table("U:/Data/Danxin/CYP3A4_Microarray/CYP3A4_Microarray.csv", sep=",", header=T)
data[1:5,1:10]

data.2 <- read.table("U:/Data/Danxin/CYP3A4_Microarray/moreCYPs.csv", sep=",", header=T)
data.2[1:5,1:10]

dim(data.2) # [1] 427  43
dim(data)  # [1] 427  82

data[1:5,1:10]
#           X     AHR   AHR.1   AHRR    ARNT  ARNT.1   CEBPA  CEBPB  CEBPD   CEBPG
# 1 GSM242213 -0.6493 -0.6612 0.2826 -0.1619 -0.0474 -0.2233 0.4964 0.1738  0.0059
# 2 GSM242214      NA -0.3682 0.0336 -0.0940  0.0750 -0.0048 0.2513 0.1068 -0.2858
# 3 GSM242215 -0.2992 -0.3068 0.3418 -0.0732 -0.0996  0.1988 0.5854 0.2192  0.0251
# 4 GSM242216 -0.3256 -0.3363 0.2179 -0.0157  0.1440  0.0129 0.5688 0.0287  0.0670
# 5 GSM242217 -0.4974 -0.5120 0.1424 -0.0809  0.1111  0.1000 0.2296 0.2814 -0.1105

data.2[1:5,1:10]
#           CYP4F11  CYP1A1  CYP3A7  CYP2U1 CYP51A1 CYP2D6  CYP2C9  CYP2B6 CYP2C19  CYP1B1
# GSM242213  0.1840  0.2021 -0.6217 -0.0769 -0.0534 0.4206  0.1318  0.1976  0.3618  0.4384
# GSM242214  0.0527 -0.2326 -1.4490  0.0717  0.2282 0.3601  0.0240 -0.2499  0.2895  0.2078
# GSM242215  0.3110  0.3020  0.1535 -0.2346  0.3859 0.6635  0.2906 -0.1127  0.4624 -0.3988
# GSM242216  0.0970 -0.0447 -0.9564  0.1200 -0.1218 0.5671 -0.1834 -0.8589  0.0688  0.1541
# GSM242217  0.1087 -0.1591 -0.0728 -0.1518 -0.2899 0.4622  0.3110  0.3884  0.4838 -0.1837

table(data[,1]==rownames(data.2))
# TRUE 
# 427 

data.with.more.CYPs <- cbind(data, data.2)
dim(data.with.more.CYPs)
# [1] 427 125
write.csv(data.with.more.CYPs, "U:/Data/Danxin/CYP3A4_Microarray/CYP3A4_Microarray_with_more_CYPs.csv")

##------------------------------------- load 01 07 2016 post-analysis summary
# save(CYPs.sorted, input.genes, file="U:/Data/Danxin/CYP3A4_Microarray/CYP3A4_Microarray_input_CYPs_separated.Rdata")
#save(CYPs.sorted, input.genes, CYPs.sorted.m, input.genes.m, scan.all.CYPs, scan.all.CYPs.against.noise, sig.est.gs, sig.est.gs.with.noise, file="U:/Data/Danxin/CYP3A4_Microarray/Sobol_selection_result1.Rdata")
#save(CYPs.sorted, input.genes, CYPs.sorted.m, input.genes.m, scan.all.CYPs, scan.all.CYPs.against.noise, sig.est.gs, sig.est.gs.with.noise, scan.all.CYPs.against.noise.insig.genes, insig.est.gs.with.noise, file="U:/Data/Danxin/CYP3A4_Microarray/Sobol_selection_result2.Rdata")
load("U:/Data/Danxin/CYP3A4_Microarray/Sobol_selection_result2.Rdata")
ls()
[1] "CYPs.sorted"           ## response data frame                 
[2] "CYPs.sorted.m"         ## response matrix                 
[3] "input.genes"           ## inputs data frame                 
[4] "input.genes.m"         ## inputs matrix                 
[5] "insig.est.gs.with.noise"               ## above average picks from the pool of irrelevant genes and artificial noise
[6] "scan.all.CYPs"                         ## model fitting without any artificial noise, gene-gene comparison  
[7] "scan.all.CYPs.against.noise"           ## signal genes tested against artificial noise
[8] "scan.all.CYPs.against.noise.insig.genes"     ## irrelevant genes tested against artificial noise
[9] "sig.est.gs"                             ## separation of signal genes and irrelavent genes
[10] "sig.est.gs.with.noise"                  ## above average picks from the pool of signal genes and artificial noise

names(scan.all.CYPs[[1]])
[1] "glm.coef.with.intercept" "glm.p.with.intercept"   
[3] "input.se.est"            "input.rho.est"          
[5] "est.gs"                  "est.gs.order"           
[7] "glm.p.order"      

dim(CYPs.sorted) # [1] 427  46 ## this is "data.2" plus CYP3A4 columns 
dim(input.genes) # [1] 427  78 ## this is "data" exclude CYP3A4 columns

CYPs.sorted.m <- as.matrix(CYPs.sorted)
dim(CYPs.sorted.m)  # [1] 427  46
input.genes.m <- as.matrix(input.genes)
dim(input.genes.m)   # [1] 427  78


##---------------------- analysis code starts from here:
ptm <- proc.time()
scan.all.CYPs <- apply(CYPs.sorted.m, 2, function(y){scan.one.CYP(y,input.genes.m)})
proc.time()-ptm

names(scan.all.CYPs[[1]])

##--------- variable selection via fdr-adjusted p-values ---------------
glm.p.c <- do.call(rbind, lapply(scan.all.CYPs, function(x){x[["glm.p.with.intercept"]]}))
glm.p.c.fdr <- t(apply(glm.p.c, 1, function(x){p.adjust(x,"fdr")}))
## number of variables selected by controlling fdr at 0.05
table(apply(glm.p.c.fdr, 1, function(x){length(which(x<=0.05))}))

##----------- fdr adjusted p-value selection (fdr=0.05) ------------
sig.glm.p.fdr <- apply(glm.p.c.fdr, 1, function(x){which(x<=0.05)})
##sig.glm.p.fdr.v <- unlist(sig.glm.p.fdr)
##length(sig.glm.p.fdr.v)
sig.glm.p.fdr

glm.p.c.fdr.order <- t(apply(glm.p.c.fdr,1,function(x){order(x,decreasing=FALSE, na.last=TRUE)}))
dim(glm.p.c.fdr.order)
input.gene.names <- names(input.genes)
names <- c("intercept",input.gene.names)
apply(glm.p.c.fdr.order[,1:10],1,function(x){names[x]})


##----------
##---------- variable selection via Sobol indices ------------------------
est.gs.c <- do.call(rbind, lapply(scan.all.CYPs, function(x){x[["est.gs"]]}))

##---------- mean cut off of sobol indices  ----------------------
##---- this cufoff criteria works as long as the number of important genes is less than unimportant genes
##----------------------------------------------------------------
sig.est.gs <- apply(est.gs.c, 1, function(x){which( x > mean(x) )})
length(sig.est.gs)
sig.est.gs[1:5]
unlist(lapply(sig.est.gs, function(x){length(x)}))
dim(input.genes)

ptm <- proc.time()
scan.all.CYPs.against.noise <- apply(matrix(1:46,nr=1), 2, function(i){test.sig.against.noise(CYPs.sorted.m[,i],input.genes,sig.est.gs[[i]])})
proc.time()-ptm

##----------
##---------- noise cut off of sobol indices  ----------------------
sig.est.gs.with.noise <- lapply(scan.all.CYPs.against.noise, function(x){gs <- x[[1]][["est.gs"]];return( which(gs[c(1:x[[2]])] > max(gs[-c(1:x[[2]])] ) )); })

sig.est.gs.with.noise[[1]]
sig.est.gs[[1]]

sig.est.gs.with.noise[[2]]
sig.est.gs[[2]]

plot(scan.all.CYPs.against.noise[[1]][[1]][["est.gs"]])
scan.all.CYPs.against.noise[[1]][[3]]

##------------- confirm no signal genes --------------
index <- c(1:78)
ptm <- proc.time()
scan.all.CYPs.against.noise.insig.genes <- apply(matrix(1:46,nr=1), 2, function(i){test.sig.against.noise(CYPs.sorted.m[,i],input.genes,index[-sig.est.gs[[i]] ] )})
proc.time()-ptm

insig.est.gs.with.noise <- lapply(scan.all.CYPs.against.noise.insig.genes, function(x){gs <- x[[1]][["est.gs"]];return( which(gs[c(1:x[[2]])] > max(gs[-c(1:x[[2]])] ) )); })

plot(scan.all.CYPs.against.noise.insig.genes[[1]][[1]][["est.gs"]])
insig.est.gs.with.noise[[1]]
index[-sig.est.gs[[1]] ]
insig.est.gs.with.noise[[2]]
index[-sig.est.gs[[2]] ]

plot(scan.all.CYPs.against.noise[[2]][[1]][["est.gs"]])
dev.new()
plot(scan.all.CYPs.against.noise.insig.genes[[2]][[1]][["est.gs"]])

##-------------- input normality test
input.normality.test <- unlist( apply(input.genes.m, 2, function(x){t <- shapiro.test(x); return(t$p.value)}) )
length(which(input.normality.test<=0.0001))
hist(input.genes[,3],breaks=30)

dim(CYPs.sorted) # [1] 427  46
dim(input.genes) # [1] 427  78

CYPs.sorted.m <- as.matrix(CYPs.sorted)
dim(CYPs.sorted.m)  # [1] 427  46
input.genes.m <- as.matrix(input.genes)
dim(input.genes.m)   # [1] 427  78

CYPs.sorted.m[1:3,1:5]
#            CYP1A1 CYP1A1.1  CYP1B1 CYP1B1.1 CYP1B1.2
# GSM242213  0.2021   0.2235  0.4384  -0.0338   0.4809
# GSM242214 -0.2326  -0.9000  0.2078  -0.1180   0.2592
# GSM242215  0.3020  -0.1735 -0.3988   0.0397  -0.3320

input.genes[1:3,1:5]
#               AHR   AHR.1   AHRR    ARNT  ARNT.1
# GSM242213 -0.6493 -0.6612 0.2826 -0.1619 -0.0474
# GSM242214      NA -0.3682 0.0336 -0.0940  0.0750
# GSM242215 -0.2992 -0.3068 0.3418 -0.0732 -0.0996

ptm <- proc.time()
scan.all.CYPs <- apply(CYPs.sorted.m, 2, function(y){scan.one.CYP(y,input.genes.m)})
proc.time()-ptm
# user  system elapsed 
# 3.96    0.07    4.02 

names(scan.all.CYPs[[1]])
# [1] "glm.coef.with.intercept" "glm.p.with.intercept"    "input.se.est"           
# [4] "input.rho.est"           "est.gs"                  "est.gs.order"           
# [7] "glm.p.order"            

##---------
glm.p.c <- do.call(rbind, lapply(scan.all.CYPs, function(x){x[["glm.p.with.intercept"]]}))
glm.p.c.fdr <- t(apply(glm.p.c, 1, function(x){p.adjust(x,"fdr")}))

## number of variables selected by controlling fdr at 0.05
table(apply(glm.p.c.fdr, 1, function(x){length(which(x<=0.05))}))
#  0  1  2  4 
# 38  6  1  1 

##----------- fdr adjusted p-value selection (fdr=0.05) ------------
sig.glm.p.fdr <- apply(glm.p.c.fdr, 1, function(x){which(x<=0.05)})
##sig.glm.p.fdr.v <- unlist(sig.glm.p.fdr)
##length(sig.glm.p.fdr.v)
sig.glm.p.fdr

glm.p.c.fdr.order <- t(apply(glm.p.c.fdr,1,function(x){order(x,decreasing=FALSE, na.last=TRUE)}))
dim(glm.p.c.fdr.order)   # [1] 46 79

input.gene.names <- names(input.genes)
names <- c("intercept",input.gene.names)
apply(glm.p.c.fdr.order[,1:10],1,function(x){names[x]})
#       CYP1A1    CYP1A1.1  CYP1B1      CYP1B1.1    CYP1B1.2    CYP2A6      CYP2A6.1    CYP2B6     
# [1,] "ONECUT1" "ONECUT1" "intercept" "NR0B1"     "intercept" "intercept" "intercept" "intercept"
# [2,] "NR1I3.1" "RXRA"    "RXRA.1"    "intercept" "RXRA.1"    "AHR.1"     "AHR"       "NFE2L2"   
# [3,] "AHRR"    "ARNT"    "AHR"       "PPARD.1"   "AHR"       "AHRR"      "AHR.1"     "AHR.1"    
# [4,] "HNF4G.1" "CEBPG"   "AHR.1"     "RXRG"      "AHR.1"     "CEBPD"     "AHRR"      "NCOR2.1"  
# [5,] "NCOR2"   "DBP"     "ARNT"      "VDR.1"     "ARNT"      "CEBPG"     "ARNT.1"    "AHRR"     
# [6,] "NCOR2.2" "HNF4G.1" "ARNT.1"    "CEBPD"     "ARNT.1"    "CEBPG.1"   "CEBPB"     "AHR"      
# [7,] "NR0B1"   "NCOA3.1" "CEBPB"     "FOXA1"     "CEBPB"     "ESR1"      "CEBPD"     "HNF4A.3"  
# [8,] "NR1I2"   "NCOR2"   "CEBPG.1"   "FOXA1.1"   "CEBPG.1"   "ESR1.2"    "CEBPG"     "HNF4A.2"  
# [9,] "NR1I2.1" "NR1H2"   "ESR1"      "NCOA1"     "DBP"       "FOXA1"     "CEBPG.1"   "NCOA1"    
# [10,] "NR1I3"   "PPARA"   "ESR1.2"    "NCOR2.1"   "ESR1"      "HNF4A"     "DBP"       "DBP"      
# 
#       CYP2B6.1    CYP2C18   CYP2C18.1 CYP2C19     CYP2C19.1   CYP2C8      CYP2C8.1    CYP2C9     
# [1,] "NFE2L2"    "NR1H3"   "NR1H3"   "intercept" "intercept" "intercept" "intercept" "intercept"
# [2,] "intercept" "CEBPA"   "NR1I3.1" "NR3C1"     "AHR"       "DBP"       "DBP"       "AHRR"     
# [3,] "AHR.1"     "NR1I3.1" "CEBPA"   "CEBPA"     "AHR.1"     "HNF4A.2"   "FOXA1"     "CEBPA"    
# [4,] "NCOR2.1"   "HNF4A.1" "HNF4A.1" "NCOA3.1"   "AHRR"      "HNF4A.3"   "NCOA3.1"   "CEBPD"    
# [5,] "AHRR"      "NCOR2.1" "NR1H4"   "NCOR2"     "ARNT"      "NCOA1"     "NR1H2"     "HNF4A.2"  
# [6,] "DBP"       "NR0B1"   "NCOR2.1" "NCOR2.2"   "ARNT.1"    "AHR.1"     "NR3C1"     "HNF4A.3"  
# [7,] "NR2F1"     "NR1H4"   "NR0B1"   "NR1H2"     "CEBPA"     "NFE2L2"    "PGRMC1"    "NCOA3.1"  
# [8,] "AHR"       "NCOR2"   "NCOR2"   "PGRMC1"    "CEBPB"     "NR1H2"     "PPARG"     "NCOR2"    
# [9,] "HNF4A.2"   "DBP"     "NR2F1.1" "NCOA1"     "CEBPD"     "NR2F1"     "HNF4A.2"   "NCOR2.1"  
# [10,] "HNF4A.3"   "FOXA3"   "FOXA3"   "PPARG.2"   "CEBPG"     "NR2F1.1"   "HNF4A.3"   "NCOR2.2"  
# 
#       CYP2C9.1  CYP2D6    CYP2D6.1    CYP2E1    CYP2E1.1  CYP2F1      CYP2F1.1    CYP2J2    CYP2J2.1 
# [1,] "NR3C1"   "ARNT"    "intercept" "CEBPG"   "CEBPG.1" "HNF4G"     "intercept" "NR1I3.1" "NCOR1"  
# [2,] "AHRR"    "CEBPA"   "AHR"       "CEBPG.1" "FOXA1"   "NR0B1"     "AHR"       "NCOR1"   "NR1I3.1"
# [3,] "CEBPD"   "ESR1"    "AHR.1"     "FOXA1"   "FOXA3"   "USF1"      "AHR.1"     "NR1I3"   "NR1I3"  
# [4,] "NR1H2"   "ESR1.1"  "AHRR"      "FOXA1.1" "HNF4A.1" "intercept" "AHRR"      "PPARG.1" "CEBPD"  
# [5,] "NR1I2"   "ESR1.2"  "ARNT"      "FOXA3"   "NCOR1"   "NCOR2.1"   "ARNT"      "PPARG.2" "PPARA.2"
# [6,] "PGRMC1"  "FOXA1.1" "ARNT.1"    "HNF4A.1" "NR3C1"   "RXRB"      "ARNT.1"    "NCOA3"   "PPARG.1"
# [7,] "PPARG.2" "FOXA2"   "CEBPA"     "HNF4G.1" "PPARA.5" "ONECUT1"   "CEBPA"     "PPARA.2" "PPARG.2"
# [8,] "RXRA"    "HNF4G"   "CEBPB"     "NCOR1"   "PPARD"   "NR2F1"     "CEBPB"     "AHRR"    "ARNT"   
# [9,] "CEBPA"   "HNF4G.1" "CEBPD"     "NCOR2.1" "PPARD.2" "NFE2L2"    "CEBPD"     "ARNT"    "NR1I2.1"
# [10,] "RXRG.1"  "NCOA1"   "CEBPG"     "NR0B1"   "NR1I3.1" "NR1I3.1"   "CEBPG"     "CEBPD"   "AHRR"   
# 
#       CYP2R1    CYP2R1.1  CYP2U1    CYP2U1.1  CYP3A4      CYP3A4.1    CYP3A4.2    CYP3A43    
# [1,] "HNF4G"   "NR0B1"   "HNF4G"   "NCOR2.1" "intercept" "intercept" "intercept" "intercept"
# [2,] "AHRR"    "VDR"     "THRB"    "THRB"    "NFE2L2"    "NFE2L2"    "RXRB.1"    "ESR1"     
# [3,] "ARNT"    "VDR.2"   "FOXA3"   "HNF4G"   "NR0B2"     "RXRB.1"    "NFE2L2"    "ESR1.2"   
# [4,] "CEBPA"   "CEBPA"   "NCOR1"   "YY1"     "NR1D2.1"   "AHRR"      "NR2F1"     "FOXA2"    
# [5,] "CEBPD"   "HNF4A"   "RXRG"    "NR0B1"   "NR2F1"     "ARNT"      "NR2F1.1"   "HNF4G"    
# [6,] "CEBPG.1" "HNF4G"   "ESR1.2"  "PPARG"   "NR2F1.1"   "DBP"       "DBP"       "NCOA1"    
# [7,] "DBP"     "HNF4G.1" "NR1I3"   "AHR"     "PPARD.2"   "ESR1"      "NCOA3.1"   "NR1I3"    
# [8,] "ESR1.1"  "NR1H3"   "CEBPG.1" "AHR.1"   "RXRB.1"    "ESR1.2"    "PPARD.2"   "PPARG"    
# [9,] "FOXA1"   "ONECUT1" "ESR1"    "NR3C1"   "AHRR"      "FOXA2"     "FOXA2"     "VDR.2"    
# [10,] "HNF4A"   "PPARA.3" "NR1H2"   "RXRB"    "ARNT"      "NCOA3.1"   "HNF4A.2"   "RXRA.1"   
# 
#       CYP3A5    CYP3A5.1  CYP3A7      CYP3A7.1    CYP4F11   CYP4F11.1 CYP4F12   CYP4F12.1 CYP51A1    
# [1,] "NR1D2.1" "NR1H2.1" "intercept" "intercept" "NR1H2"   "NR1H2"   "NR3C1"   "CEBPA"   "NR2F1"    
# [2,] "NR1H2.1" "NCOA3"   "RXRB.1"    "AHRR"      "NR1I3.1" "RXRB.1"  "ONECUT1" "NR1I2.1" "NR2F1.1"  
# [3,] "RXRA"    "RXRA"    "AHRR"      "ARNT"      "NR0B1"   "THRA"    "PPARA.4" "NR3C1"   "intercept"
# [4,] "NCOA1"   "PPARD.1" "ARNT"      "ESR1"      "CEBPG.1" "CEBPD"   "THRB"    "NR1I2"   "ESR1"     
# [5,] "PPARD.1" "NCOA1"   "DBP"       "ESR1.2"    "FOXA1"   "HNF4A.2" "NR1I2"   "NR2F2"   "ESR1.2"   
# [6,] "FOXA3"   "NR1D2.1" "ESR1"      "FOXA1"     "FOXA3"   "HNF4A.3" "NR1I2.1" "NR2F2.1" "FOXA2"    
# [7,] "NCOA3"   "USF1"    "ESR1.2"    "FOXA2"     "HNF4G.1" "RXRB"    "USF1"    "ONECUT1" "NCOA1"    
# [8,] "NFE2L2"  "NCOR1"   "FOXA1"     "NFE2L2"    "NCOR1"   "NR1H4"   "NR1D2.1" "CEBPD"   "NCOR2.1"  
# [9,] "ARNT"    "RXRG.1"  "FOXA2"     "NR0B2"     "NCOR2.1" "FOXA3"   "NR2F2.1" "DBP"     "NCOR2.2"  
# [10,] "ESR1.1"  "PPARA.4" "HNF4A.2"   "RXRB.1"    "NR1I3"   "HNF4G"   "VDR.1"   "THRA"    "NR0B1"    
# 
#       CYP51A1.1 CYP51A1.2 CYP7A1      CYP7A1.1   
# [1,] "CEBPB"   "FOXA2"   "intercept" "intercept"
# [2,] "CEBPG.1" "NCOA1"   "NR1D2"     "AHR"      
# [3,] "NCOR1"   "NCOR2.1" "PPARA.1"   "NR1D2"    
# [4,] "NR1D2"   "NR0B1"   "VDR.2"     "NR5A2"    
# [5,] "PPARA.2" "NR1H2.1" "FOXA3"     "PPARA.1"  
# [6,] "PPARA.5" "NR2F1"   "NCOR2.2"   "PPARA.5"  
# [7,] "RXRG.1"  "NR2F1.1" "PPARA.5"   "VDR.2"    
# [8,] "YY1"     "PPARD.1" "NCOR2"     "AHR.1"    
# [9,] "CEBPG"   "PPARG.1" "VDR"       "FOXA3"    
# [10,] "NR5A2"   "PPARG.2" "AHR"       "HNF4G"    

sig.est.gs[1:5]
# $CYP1A1
# [1]  1  2  5  9 10 11 12 14 21 22 28 29 32 33 34 37 40 42 43 44 45 53 57 58 59 60 61 62 66 71 76
# 
# $CYP1A1.1
# [1]  5  9 10 11 12 14 21 22 25 28 42 43 44 45 53 57 58 59 60 61 62 66 72 77
# 
# $CYP1B1
# [1]  4  6 11 15 19 21 22 25 35 41 42 43 44 45 54 56 57 58 59 61 68 73 75 77 78
# 
# $CYP1B1.1
# [1]  1  2  4  8  9 10 12 14 17 18 22 30 31 32 34 38 42 43 44 45 52 53 55 56 60 61 62 69 72 73 74 76
# 
# $CYP1B1.2
# [1]  4  6 11 15 19 20 21 22 25 28 35 41 42 43 44 45 48 49 56 57 58 59 61 68 73 75 77 78

unlist(lapply(sig.est.gs, function(x){length(x)}))
# CYP1A1  CYP1A1.1    CYP1B1  CYP1B1.1  CYP1B1.2    CYP2A6  CYP2A6.1    CYP2B6  CYP2B6.1   CYP2C18 CYP2C18.1 
# 31        24        25        32        28        25        24        28        29        23        23 
# CYP2C19 CYP2C19.1    CYP2C8  CYP2C8.1    CYP2C9  CYP2C9.1    CYP2D6  CYP2D6.1    CYP2E1  CYP2E1.1    CYP2F1 
# 25        29        29        27        25        27        32        29        26        28        10 
# CYP2F1.1    CYP2J2  CYP2J2.1    CYP2R1  CYP2R1.1    CYP2U1  CYP2U1.1    CYP3A4  CYP3A4.1  CYP3A4.2   CYP3A43 
# 27        24        25        24         9        32        29        25        26        26        28 
# CYP3A5  CYP3A5.1    CYP3A7  CYP3A7.1   CYP4F11 CYP4F11.1   CYP4F12 CYP4F12.1   CYP51A1 CYP51A1.1 CYP51A1.2 
# 25        23        27        23        29        29        25        27        33        19        33 
# CYP7A1  CYP7A1.1 
# 29        28 

dim(input.genes) # [1] 427  78

##---------------------------------
ptm <- proc.time()
scan.all.CYPs.against.noise <- apply(matrix(1:46,nr=1), 2, function(i){test.sig.against.noise(CYPs.sorted.m[,i],input.genes.m,sig.est.gs[[i]])})
proc.time()-ptm
# user  system elapsed 
# 3.01    0.15    3.18 

names(scan.all.CYPs.against.noise[[1]][[1]])
# [1] "glm.coef.with.intercept" "glm.p.with.intercept"    "input.se.est"            "input.rho.est"          
# [5] "est.gs"                  "est.gs.order"            "glm.p.order"            

##------------- confirm no signal genes --------------
index <- c(1:78)
ptm <- proc.time()
scan.all.CYPs.against.noise.insig.genes <- apply(matrix(1:46,nr=1), 2, function(i){test.sig.against.noise(CYPs.sorted.m[,i],input.genes,index[-sig.est.gs[[i]] ] )})
proc.time()-ptm
# user  system elapsed 
# 6.02    0.17    6.19 

insig.est.gs.with.noise <- lapply(scan.all.CYPs.against.noise.insig.genes, function(x){gs <- x[[1]][["est.gs"]];return( which(gs[c(1:x[[2]])] > max(gs[-c(1:x[[2]])] ) )); })

plot(scan.all.CYPs.against.noise.insig.genes[[1]][[1]][["est.gs"]])

##----------------------------------- Example of analyzing one CYP at a time: 
## y1 is the response variable=CYP expression of interest
## y1 <- CYPs.sorted.m[,30] ## CYP3A4
index.2m.all <- combn(c(1:78), 2)
dim(index.2m.all)

gene.names.2m.all <- apply(index.2m.all, 2, function(i){gene.names[i]})
all.gene.pairs.3 <- apply(index.2m.all, 2, function(i){SI.comb(y1, input.genes.m[,i], k=3)})
order.2m.3 <- order(-all.gene.pairs.3[2,], all.gene.pairs.3[1,],decreasing=TRUE)

gene.names.2m.all[,order.2m.3][,1:30]
all.gene.pairs.3[,order.2m.3][,1:30]

##------------------------------------------------------- 
index.3m.all <- combn(c(1:78), 3)
dim(index.3m.all)

gene.names.3m.all <- apply(index.3m.all, 2, function(i){gene.names[i]})
all.gene.triplets.3 <- apply(index.3m.all, 2, function(i){SI.comb(y1, input.genes.m[,i], k=3)})
order.3m.3 <- order(-all.gene.triplets.3[2,], all.gene.triplets.3[1,],decreasing=TRUE)

gene.names.3m.all[,order.3m.3][,1:30]
all.gene.triplets.3[,order.3m.3][,1:30]

table(gene.names.2m.all[,order.2m.3][,1:100])
table(gene.names.3m.all[,order.3m.3][,1:100])
table(gene.names.2m.all[,order.2m.3][,1:200])
table(gene.names.3m.all[,order.3m.3][,1:200])
# save.image("U:\\Data\\Danxin\\CYP3A4_Microarray\\gene_triplets_picks_02_23_2016_fullPoly3.Rdata")

dev.new();hist(all.gene.triplets.3[2,], breaks=1200)
dev.new();hist(all.gene.triplets.3[1,order.3m.3][1:200], breaks=50)
##-------------------------------------------------------
## ----------------------------- code for scanning for subset effect using polynomial with degree k=3 and stepwise model selection by AIC:
## need to impute all missing values before using the following code

library(MASS)

SI <- function(y, input){
  colnames(input) <- NULL
  m.full <- glm(y~., data=data.frame(input));
  fit.step <- stepAIC(m.full, direction="both", trace=F)	
  m.best <- formula(fit.step)
  fit.1 <- glm(m.best, data=data.frame(input));	glm.coef <- fit.1$coef;
  v.select <- input[, as.integer(sub("X", "", attr(terms(fit.step),"term.labels") ))]
  SI <- var( matrix(v.select, nr=dim(input)[1]) %*% matrix(glm.coef[-1], nc=1) , na.rm=TRUE )
  return(c(SI, fit.1$deviance / fit.1$null.deviance, fit.1$null.deviance))
}

## k >=3
comb.m <- function(input, k=3){
  d <- dim(input)[2]
  if(d>=3){
    combn.index <- apply(matrix(2:(min(d,k)), nr=1),2, function(x){combn(c(1:d), x)} ) 
    
    combn.list <- lapply(combn.index, function(x){cc<-input[,x[1,]]; for(j in 2:dim(x)[1]){cc <- cc*input[,x[j,]]}; return(list(cc))})
    combn.input <- do.call(cbind, lapply(combn.list, function(x){x[[1]]}))
  }
  if(d==2){combn.input <- input[,1]*input[,2]}
  return(combn.input)
}

## code tested up to slection of 4comb with k=3:
comb.m <- function(input, k=3){
  d <- dim(input)[2]
  combn.index <- apply(matrix(2:(min(d,k)), nr=1),2, function(x){list(combn(c(1:d), x))} ) 
  
  combn.list <- lapply(combn.index, function(x){cc<-input[,x[[1]][1,]]; 
  for(j in 2:dim(x[[1]])[1]){cc <- cc*input[,x[[1]][j,]]}; 
  return(list(cc));
  })
  combn.input <- do.call(cbind, lapply(combn.list, function(x){x[[1]]}))
  
  return(combn.input)
}


SI.comb <- function(y, input, k=3){
  main.new.m <- apply(input, 2, function(x){ apply(matrix(1:k),1, function(i){x^i}) })
  main.SI <- apply(main.new.m, 2, function(x){SI(y,matrix(x, nr=dim(input)[1]))})
  
  comb.new.m <- cbind(matrix(main.new.m, nr=dim(input)[1]), comb.m(input, k))
  comb.SI <- SI(y, comb.new.m)
  
  return(cbind(comb.SI, main.SI))
}


##------------------------- 03 01 2016
SI.main.total <- function(i, y, input, k=3){
  main.input <- apply(matrix(c(1:k), nc=1), 1, function(j){input[,i]^j})
  main.SI <- SI(y, main.input)
  
  total.input.1 <- apply(input[,-i], 2, function(x){ apply(matrix(1:k),1, function(i){x^i}) })
  total.input <- cbind(matrix(total.input.1, nr=dim(input)[1]), comb.m(input[,-i], k))
  total.SI <- SI(y, input[,-i])
  return(cbind(total.SI, main.SI))
}

all.gene.main.total.SI.3.stepwise <- apply(matrix(1:78, nr=1), 2, function(i){SI.main.total(i, y1, input, k=3)})
save.image("U:\\Data\\Danxin\\CYP3A4_Microarray\\All_single_gene_main_total_SI_03_01_2016_Poly3stepwise.Rdata")

summary(all.gene.main.total.SI.3.stepwise[1,])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3477  0.3610  0.3620  0.3608  0.3620  0.3724 
which(all.gene.main.total.SI.3.stepwise[1,]>0.372)  # [1] 59

M.T.SI <- data.frame(gene.names, t(all.gene.main.total.SI.3.stepwise))
dim(M.T.SI) # [1] 78  7
names(M.T.SI) <- c("Gene", "T.SI", "T.ResDev.Perc", "T.NullDev", "M.SI", "M.ResDev.Perc", "M.NullDev")

M.T.SI[1:3,]
#    Gene      T.SI T.ResDev.Perc T.NullDev        M.SI M.ResDev.Perc M.NullDev
# 1   AHR 0.3619527     0.2035096  184.8586 0.068721840     0.8439045  184.8586
# 2 AHR.1 0.3619527     0.2035096  184.8586 0.073308624     0.8338740  184.8586
# 3  AHRR 0.3619527     0.2035096  184.8586 0.008624813     0.9842717  184.8586

order.1m.3.stepwise <- order(-M.T.SI$M.ResDev.Perc, M.T.SI$M.SI,decreasing=TRUE)
M.T.SI.report <- M.T.SI[,-c(2:4)]
M.T.SI.report[order.1m.3.stepwise[1:30],]

M.SI.single <- M.T.SI.report[order.1m.3.stepwise,]
M.T.SI.report[order.1m.3.stepwise[1:30],]
#       Gene       M.SI M.ResDev.Perc M.NullDev
# 12    ESR1 0.22414915     0.4922417  184.8586
# 14  ESR1.2 0.22208708     0.4965189  184.8586
# 53  PGRMC1 0.18834937     0.5704080  184.8586
# 44   NR1I3 0.18398708     0.5865769  184.8586
# 45 NR1I3.1 0.18525503     0.5876062  184.8586
# 42   NR1I2 0.14112393     0.6805029  184.8586
# 33  NFE2L2 0.14040139     0.6813296  184.8586
# 43 NR1I2.1 0.13823927     0.6867121  184.8586
# 32 NCOR2.2 0.13500244     0.6952524  184.8586
# 17   FOXA2 0.12946166     0.7069641  184.8586
# 29   NCOR1 0.12675841     0.7122877  184.8586
# 21 HNF4A.2 0.10978822     0.7519701  184.8586
# 55 PPARA.1 0.11699504     0.7594395  184.8586
# 30   NCOR2 0.10743229     0.7601134  184.8586
# 22 HNF4A.3 0.10260127     0.7650642  184.8586
# 74    USF1 0.10706339     0.7667116  184.8586
# 40   NR1H3 0.08129362     0.8169474  184.8586
# 2    AHR.1 0.07330862     0.8338740  184.8586
# 73    THRB 0.09963322     0.8401774  184.8586
# 67  RXRA.1 0.06920875     0.8416139  184.8586
# 1      AHR 0.06872184     0.8439045  184.8586
# 4     ARNT 0.06796497     0.8464045  184.8586
# 38   NR1H2 0.06957239     0.8473280  184.8586
# 9    CEBPG 0.06718041     0.8495164  184.8586
# 60   PPARD 0.05940257     0.8675627  184.8586
# 51   NR5A2 0.05662914     0.8720955  184.8586
# 8    CEBPD 0.05115976     0.8836267  184.8586
# 15   FOXA1 0.04735714     0.8937608  184.8586
# 62 PPARD.2 0.04597317     0.8972174  184.8586
# 35   NR0B2 0.04275165     0.9035985  184.8586

order.1m.3.stepwise <- order(-M.T.SI$T.ResDev.Perc, M.T.SI$T.SI,decreasing=TRUE)
M.T.SI.report <- M.T.SI[,-c(5:7)]
M.T.SI.report[order.1m.3.stepwise[1:30],]

T.SI.single <- M.T.SI.report[order.1m.3.stepwise,]
M.T.SI.report[order.1m.3.stepwise[1:30],]
# Gene      T.SI T.ResDev.Perc T.NullDev
# 4     ARNT 0.3646555     0.2027955  184.8586
# 57 PPARA.3 0.3649793     0.2030767  184.8586
# 43 NR1I2.1 0.3628726     0.2032585  184.8586
# 73    THRB 0.3586479     0.2034640  184.8586
# 1      AHR 0.3619527     0.2035096  184.8586
# 2    AHR.1 0.3619527     0.2035096  184.8586
# 3     AHRR 0.3619527     0.2035096  184.8586
# 7    CEBPB 0.3619527     0.2035096  184.8586
# 8    CEBPD 0.3619527     0.2035096  184.8586
# 14  ESR1.2 0.3619527     0.2035096  184.8586
# 16 FOXA1.1 0.3619527     0.2035096  184.8586
# 18   FOXA3 0.3619527     0.2035096  184.8586
# 19   HNF4A 0.3619527     0.2035096  184.8586
# 22 HNF4A.3 0.3619527     0.2035096  184.8586
# 25   NCOA1 0.3619527     0.2035096  184.8586
# 26   NCOA2 0.3619527     0.2035096  184.8586
# 27   NCOA3 0.3619527     0.2035096  184.8586
# 29   NCOR1 0.3619527     0.2035096  184.8586
# 31 NCOR2.1 0.3619527     0.2035096  184.8586
# 35   NR0B2 0.3619527     0.2035096  184.8586
# 36   NR1D2 0.3619527     0.2035096  184.8586
# 42   NR1I2 0.3619527     0.2035096  184.8586
# 45 NR1I3.1 0.3619527     0.2035096  184.8586
# 46   NR2F1 0.3619527     0.2035096  184.8586
# 47 NR2F1.1 0.3619527     0.2035096  184.8586
# 48   NR2F2 0.3619527     0.2035096  184.8586
# 50   NR3C1 0.3619527     0.2035096  184.8586
# 51   NR5A2 0.3619527     0.2035096  184.8586
# 53  PGRMC1 0.3619527     0.2035096  184.8586
# 54   PPARA 0.3619527     0.2035096  184.8586

##------------------------------------------------------- 03 21 2016 Generate csv files for Gephi:
dim(all.gene.pairs.3.stepwise)  # [1]    9 3003
all.gene.pairs.3.stepwise[,1:5]
#              [,1]         [,2]         [,3]         [,4]         [,5]
# [1,]   0.07554767 6.872184e-02   0.11255070   0.09972268 7.435926e-02
# [2,]   0.82902043 8.439045e-01   0.74339131   0.77291290 8.299962e-01
# [3,] 184.85859840 1.848586e+02 184.85859840 184.85859840 1.848586e+02
# [4,]   0.06872184 6.872184e-02   0.06872184   0.06872184 6.872184e-02
# [5,]   0.84390450 8.439045e-01   0.84390450   0.84390450 8.439045e-01
# [6,] 184.85859840 1.848586e+02 184.85859840 184.85859840 1.848586e+02
# [7,]   0.07330862 8.624813e-03   0.06796497   0.03873406 7.946078e-03
# [8,]   0.83387399 9.842717e-01   0.84640452   0.91345141 9.819607e-01
# [9,] 184.85859840 1.848586e+02 184.85859840 184.85859840 1.848586e+02

dim(gene.names.2m.all)  ## [1]    2 3003
choose(78,2)  # [1] 3003
gene.names.2m.all[,1:5]
#       [,1]    [,2]   [,3]   [,4]     [,5]   
# [1,] "AHR"   "AHR"  "AHR"  "AHR"    "AHR"  
# [2,] "AHR.1" "AHRR" "ARNT" "ARNT.1" "CEBPA"

Gephi.gene.pairs.edge <- data.frame(t(gene.names.2m.all), t(all.gene.pairs.3.stepwise[1:3,]))
dim(Gephi.gene.pairs.edge)  ## [1] 3003    5
Gephi.gene.pairs.edge[1:5,]
#    X1     X2       X1.1      X2.1       X3
# 1 AHR  AHR.1 0.07554767 0.8290204 184.8586
# 2 AHR   AHRR 0.06872184 0.8439045 184.8586
# 3 AHR   ARNT 0.11255070 0.7433913 184.8586
# 4 AHR ARNT.1 0.09972268 0.7729129 184.8586
# 5 AHR  CEBPA 0.07435926 0.8299962 184.8586
write.csv(Gephi.gene.pairs.edge, file="U:/Data/Gephi.gene.pairs.edge.csv")

Gephi.gene.pairs.edge$Label <- apply(gene.names.2m.all, 2, function(x){paste(x[1],"-",x[2])})
Gephi.gene.pairs.edge[1:5,]
#    X1     X2       X1.1      X2.1       X3        Label
# 1 AHR  AHR.1 0.07554767 0.8290204 184.8586  AHR - AHR.1
# 2 AHR   AHRR 0.06872184 0.8439045 184.8586   AHR - AHRR
# 3 AHR   ARNT 0.11255070 0.7433913 184.8586   AHR - ARNT
# 4 AHR ARNT.1 0.09972268 0.7729129 184.8586 AHR - ARNT.1
# 5 AHR  CEBPA 0.07435926 0.8299962 184.8586  AHR - CEBPA
write.csv(Gephi.gene.pairs.edge, file="U:/Data/Gephi_gene_pairs_edge.csv")

Gephi.gene.pairs.node <- data.frame(c(gene.names.2m.all[1,], gene.names.2m.all[2,]), rbind(t(all.gene.pairs.3.stepwise[4:6,]), t(all.gene.pairs.3.stepwise[7:9,])))
dim(Gephi.gene.pairs.node)  # [1] 6006    4
Gephi.gene.pairs.node.unique <- Gephi.gene.pairs.node[!duplicated(Gephi.gene.pairs.node[,1]),]
dim(Gephi.gene.pairs.node.unique)  # [1] 78  4
write.csv(Gephi.gene.pairs.node.unique, file="U:/Data/Gephi_gene_pairs_node.csv")

##-----------
dim(all.gene.triplets.3.stepwise)
all.gene.triplets.3.stepwise[,1:5]
dim(gene.names.3m.all)
choose(78,3)

gene.names.3m.all[,1:5]
Gephi.gene.triplets.edge <- data.frame(t(gene.names.3m.all), t(all.gene.triplets.3.stepwise[1:3,]))
dim(Gephi.gene.triplets.edge)

Gephi.gene.triplets.edge[1:5,]

Gephi.gene.triplets.edge$Label <- apply(gene.names.3m.all, 2, function(x){paste(x[1],"-",x[2], "-", x[3])})
Gephi.gene.triplets.edge[1:5,]
write.csv(Gephi.gene.triplets.edge, file="U:/Data/Gephi_gene_triplets_edge.csv")

Gephi.gene.triplets.edge.1 <- Gephi.gene.triplets.edge
Gephi.gene.triplets.edge.2 <- Gephi.gene.triplets.edge
Gephi.gene.triplets.edge.2[,1] <- Gephi.gene.triplets.edge.1[,2]
Gephi.gene.triplets.edge.2[,2] <- Gephi.gene.triplets.edge.1[,3]
write.csv(rbind(Gephi.gene.triplets.edge.1, Gephi.gene.triplets.edge.2), file="U:/Data/Gephi_gene_triplets_edge_in_pairs.csv")

##--------------
dim(all.gene.4comb.3.stepwise)
all.gene.4comb.3.stepwise[,1:5]
dim(gene.names.4m.all)

gene.names.4m.all[,1:5]
Gephi.gene.4comb.edge <- data.frame(t(gene.names.4m.all), t(all.gene.4comb.3.stepwise[1:3,]))
dim(Gephi.gene.4comb.edge)

dim(Gephi.gene.4comb.edge) # [1] 194580      7

Gephi.gene.4comb.edge[1:5,]
Gephi.gene.4comb.edge$Label <- apply(gene.names.4m.all, 2, function(x){paste(x[1],"-",x[2], "-", x[3], "-", x[4])})
Gephi.gene.4comb.edge[1:5,]
##write.csv(Gephi.gene.triplets.edge, file="U:/Data/Gephi_gene_triplets_edge.csv")
quantile(all.gene.4comb.3.stepwise[2,], (0:100)/100)

quantile(all.gene.4comb.3.stepwise[2,], (0:100)/100)
#          0%        1%        2%        3%        4%        5%        6%        7%        8%        9%       10%       11% 
#   0.3050708 0.3591035 0.3681105 0.3746085 0.3801594 0.3849702 0.3895254 0.3938849 0.3980541 0.4022310 0.4062780 0.4102412 
#         12%       13%       14%       15%       16%       17%       18%       19%       20%       21%       22%       23% 
#   0.4139544 0.4173212 0.4206390 0.4237780 0.4268039 0.4298029 0.4326973 0.4356101 0.4384056 0.4411423 0.4439477 0.4469685 
#         24%       25%       26%       27%       28%       29%       30%       31%       32%       33%       34%       35% 
#   0.4499103 0.4529346 0.4559858 0.4591132 0.4621087 0.4653899 0.4686134 0.4720120 0.4751040 0.4782437 0.4813679 0.4843980 
#         36%       37%       38%       39%       40%       41%       42%       43%       44%       45%       46%       47% 
#   0.4875293 0.4906413 0.4937513 0.4969185 0.4999888 0.5030339 0.5060292 0.5091500 0.5121302 0.5153868 0.5185882 0.5217753 
#         48%       49%       50%       51%       52%       53%       54%       55%       56%       57%       58%       59% 
#   0.5251140 0.5284943 0.5318367 0.5354842 0.5390414 0.5426846 0.5464264 0.5500144 0.5537387 0.5577006 0.5618032 0.5657288 
#         60%       61%       62%       63%       64%       65%       66%       67%       68%       69%       70%       71% 
#   0.5699786 0.5741031 0.5781108 0.5821493 0.5861424 0.5901255 0.5941462 0.5980342 0.6018650 0.6057665 0.6097853 0.6136810 
#         72%       73%       74%       75%       76%       77%       78%       79%       80%       81%       82%       83% 
#   0.6175909 0.6216001 0.6257048 0.6295785 0.6336936 0.6376393 0.6418813 0.6463276 0.6510995 0.6561015 0.6617206 0.6676337 
#         84%       85%       86%       87%       88%       89%       90%       91%       92%       93%       94%       95% 
#   0.6740480 0.6805090 0.6870766 0.6939313 0.7012260 0.7088078 0.7171415 0.7260902 0.7363856 0.7465898 0.7578890 0.7703898 
#         96%       97%       98%       99%      100% 
#   0.7850554 0.8021205 0.8247854 0.8549862 0.9693799 

length(which(all.gene.4comb.3.stepwise[2,]<0.4))  # [1] 16506

Gephi.gene.4comb.edge.full <- Gephi.gene.4comb.edge
Gephi.gene.4comb.edge <- Gephi.gene.4comb.edge.full[which(Gephi.gene.4comb.edge.full[,6]<0.4),]
dim(Gephi.gene.4comb.edge)

Gephi.gene.4comb.edge[1:5,]

Gephi.gene.4comb.edge.1 <- Gephi.gene.4comb.edge
Gephi.gene.4comb.edge.2 <- Gephi.gene.4comb.edge
Gephi.gene.4comb.edge.3 <- Gephi.gene.4comb.edge

Gephi.gene.4comb.edge.2[,1] <- Gephi.gene.4comb.edge.1[,2]
Gephi.gene.4comb.edge.2[,2] <- Gephi.gene.4comb.edge.1[,3]
Gephi.gene.4comb.edge.3[,1] <- Gephi.gene.4comb.edge.1[,3]
Gephi.gene.4comb.edge.3[,2] <- Gephi.gene.4comb.edge.1[,4]

write.csv(rbind(Gephi.gene.4comb.edge.1, Gephi.gene.4comb.edge.2, Gephi.gene.4comb.edge.3), file="U:/Data/Gephi_gene_4comb_edge_in_pairs.csv")
save.image("U:\\Data\\Danxin\\CYP3A4_Microarray\\gene_comb_picks_03_21_2016_Gephi_format.Rdata")



## ---------------------------- GTEx RNAseq:
## 12/27/2016
SI <- function(y, input){
  colnames(input) <- NULL
  m.full <- glm(y~., data=data.frame(input));
  fit.step <- stepAIC(m.full, direction="both", trace=F)
  m.best <- formula(fit.step)
  fit.1 <- glm(m.best, data=data.frame(input));glm.coef <- fit.1$coef;
  v.select <- input[, as.integer(sub("X", "", attr(terms(fit.step),"term.labels") ))]
  SI <- var( matrix(v.select, nr=dim(input)[1]) %*% matrix(glm.coef[-1], nc=1) , na.rm=TRUE )
  return(c(SI, fit.1$deviance / fit.1$null.deviance, fit.1$null.deviance))
}

## code tested up to slection of 4comb with k=3:
comb.m <- function(input, k=3){
  d <- dim(input)[2]
  combn.index <- apply(matrix(2:(min(d,k)), nr=1),2, function(x){list(combn(c(1:d), x))} ) 
  
  combn.list <- lapply(combn.index, function(x){cc<-input[,x[[1]][1,]]; 
  for(j in 2:dim(x[[1]])[1]){cc <- cc*input[,x[[1]][j,]]}; 
  return(list(cc));
  })
  combn.input <- do.call(cbind, lapply(combn.list, function(x){x[[1]]}))
  
  return(combn.input)
}

####### load data:
data <- read.table("U:/Data/Danxin/consider_gender_effects/CYPs_GTEX.csv", sep=",", header=T)
dim(data) # [1]  22 140
data[,1:4]
#     X               Name Description GTEX.1192X.1026.SM.5H12P
# 1   1  ENSG00000134716.5      CYP2J2               9.84567261
# 2   2  ENSG00000138061.7      CYP1B1              33.74108124
# 3   3 ENSG00000155016.13      CYP2U1               3.64712715
# 4   4 ENSG00000001630.11     CYP51A1               2.56612277
# 5   5  ENSG00000106258.9      CYP3A5              82.95661926
# 6   6  ENSG00000160870.8      CYP3A7               9.41562748
# 7   7 ENSG00000160868.10      CYP3A4              13.41103554
# 8   8 ENSG00000021461.12     CYP3A43               0.09146577
# 9   9  ENSG00000167910.3      CYP7A1               0.05567182
# 10 10  ENSG00000108242.8     CYP2C18              17.10964394
# 11 11  ENSG00000165841.5     CYP2C19               3.71595001
# 12 12  ENSG00000138109.9      CYP2C9              40.32092285
# 13 13  ENSG00000138115.9      CYP2C8              44.09194565
# 14 14  ENSG00000130649.5      CYP2E1              84.64115143
# 15 15  ENSG00000186104.6      CYP2R1               5.11535072
# 16 16  ENSG00000140465.9      CYP1A1               0.22710375
# 17 17  ENSG00000140505.6      CYP1A2               0.06020230
# 18 18 ENSG00000186204.10     CYP4F12               1.92371059
# 19 19 ENSG00000171903.12     CYP4F11               5.77074957
# 20 20  ENSG00000255974.2      CYP2A6              14.17587471
# 21 21  ENSG00000197408.4      CYP2B6               1.21407962
# 22 22 ENSG00000100197.16      CYP2D6              21.47479630

CYP.name <- as.character(data[,3])
CYP.name
# [1] "CYP2J2"  "CYP1B1"  "CYP2U1"  "CYP51A1" "CYP3A5"  "CYP3A7"  "CYP3A4" 
# [8] "CYP3A43" "CYP7A1"  "CYP2C18" "CYP2C19" "CYP2C9"  "CYP2C8"  "CYP2E1" 
# [15] "CYP2R1"  "CYP1A1"  "CYP1A2"  "CYP4F12" "CYP4F11" "CYP2A6"  "CYP2B6" 
# [22] "CYP2D6" 
CYP3A4 <- data[7,]

CYP.ENSG <- as.character(data[,2])

rpkm <- data[,4:140]
rpkm[,1:3]
rpkm <- t(rpkm)
rpkm[1:3,]

rpkm <- t(rpkm)
rpkm[1:3,]

colnames(rpkm) <- CYP.name
colnames(rpkm)

save(CYP.ENSG, rpkm, file="U:/Data/Danxin/consider_gender_effects/CYPs_GTEx.RData")

#------- need to match up two datasets to re-run the same analysis for 22 different CYPs
csv.table <- read.table("U:/Data/Danxin/consider_gender_effects/CYP3A4_GTEx_Liver_subject_info_matched_final.csv", sep=",", header=T)

names(csv.table)
# [1] "SAMPID"   "SMTS"     "SMTSD"    "GENDER"   "AGE"      "RACE"     "BMI"      "SMTORMVE" "AHR"      "AHRR"     "ARNT"    
# [12] "CEBPA"    "CEBPB"    "CEBPD"    "CEBPG"    "DBP"      "ESR1"     "FOXA1"    "FOXA2"    "FOXA3"    "HNF4A"    "HNF4G"   
# [23] "NCOA1"    "NCOA2"    "NCOA3"    "NCOR1"    "NCOR2"    "NFE2L2"   "NR0B1"    "NR0B2"    "NR1D2"    "NR1H2"    "NR1H3"   
# [34] "NR1H4"    "NR1I2"    "NR1I3"    "NR2F1"    "NR2F2"    "NR3C1"    "NR5A2"    "ONECUT1"  "PGRMC1"   "PPARA"    "PPARD"   
# [45] "PPARG"    "RXRA"     "RXRB"     "RXRG"     "THRA"     "THRB"     "USF1"     "VDR"      "YY1"      "y"       

csv.table$SAMPID[1:5]
# [1] GTEX-1192X-1026-SM-5H12P GTEX-11DXY-0526-SM-5EGGQ GTEX-11DXZ-0126-SM-5EGGY GTEX-11EQ9-0526-SM-5A5JZ GTEX-11GSP-0626-SM-5986T
# 137 Levels: GTEX-1192X-1026-SM-5H12P GTEX-11DXY-0526-SM-5EGGQ GTEX-11DXZ-0126-SM-5EGGY ... GTEX-ZZPU-0426-SM-5GZYH

old.sampleID <- as.character(csv.table$SAMPID)
old.sampleID[1:5]
# [1] "GTEX-1192X-1026-SM-5H12P" "GTEX-11DXY-0526-SM-5EGGQ" "GTEX-11DXZ-0126-SM-5EGGY" "GTEX-11EQ9-0526-SM-5A5JZ"
# [5] "GTEX-11GSP-0626-SM-5986T"

new.sampleID <- rownames(rpkm)
new.sampleID[1:5]
# [1] "GTEX.1192X.1026.SM.5H12P" "GTEX.11DXY.0526.SM.5EGGQ" "GTEX.11DXZ.0126.SM.5EGGY" "GTEX.11EQ9.0526.SM.5A5JZ"
# [5] "GTEX.11GSP.0626.SM.5986T"
length(old.sampleID) # [1] 137
length(new.sampleID) # [1] 137

new.reorder <- apply(matrix(old.sampleID, nc=1),1,function(x){which(new.sampleID==x)})

r.old.sampleID <- gsub("-","", old.sampleID)
r.new.sampleID <- gsub("[.]","", new.sampleID)

r.new.sampleID <- gsub("[.]","", new.sampleID)
r.new.sampleID[1:5]
# [1] "GTEX1192X1026SM5H12P" "GTEX11DXY0526SM5EGGQ" "GTEX11DXZ0126SM5EGGY" "GTEX11EQ90526SM5A5JZ" "GTEX11GSP0626SM5986T"
table(r.old.sampleID==r.new.sampleID)
# TRUE 
# 137 
## --------------------- two datasets can be directly combined together ----------------------

dim(rpkm) # [1] 137  22
dim(csv.table[,9:53]) # [1] 137  45
c.rpkm <- cbind(rpkm, csv.table[,9:53])

table(log(csv.table$y) >= log.CYP3A4)
# TRUE 
# 137 

plot(log(csv.table$y), log.CYP3A4)
abline(0,1)
summary(log(csv.table$y)-log.CYP3A4)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2940  0.5143  0.5855  0.5748  0.6426  0.7223

log.gene.tb <- log(gene.tb)
log.gene.tb[which(is.na(log.gene.tb))] <- -745

##--------------------------- pre-processing is done ------------------------------

y.m <- c.rpkm[1:22]
dim(c.rpkm) #[1] 137  67
gene.tb <- c.rpkm[,23:67]
which(names(gene.tb)=="CYP3A4") # integer(0)
which(names(gene.tb)=="y") # integer(0)
dim(gene.tb) # [1] 137  45

names(csv.table)[8] # [1] "SMTORMVE"
which(names(csv.table)=="y") # [1] 54

apply(y.m, 2, summary)
#         CYP2J2  CYP1B1 CYP2U1 CYP51A1  CYP3A5   CYP3A7    CYP3A4 CYP3A43   CYP7A1 CYP2C18  CYP2C19
# Min.     3.402  0.3208 0.2355  0.1134   4.566   0.2066    0.1847  0.0000   0.0000   6.124   0.5102
# 1st Qu. 15.160  1.0530 1.0990  1.0500  49.080   3.8240   43.0400  0.4127   0.1109  27.620   3.7620
# Median  21.310  1.6640 1.3750  2.4550  99.110   9.4160  210.4000  1.5090   0.9797  35.370   6.8940
# Mean    21.190  2.7410 1.4770  3.9200 105.800  23.1300  325.9000  4.1290   9.6560  36.600  18.0600
# 3rd Qu. 26.540  2.5880 1.7510  5.9160 148.100  17.8700  396.8000  5.8720   8.9300  45.360  18.1300
# Max.    53.800 48.5300 3.7460 29.2500 352.400 495.0000 4881.0000 26.1900 169.6000  75.790 184.4000
# 
#           CYP2C9    CYP2C8  CYP2E1 CYP2R1    CYP1A1  CYP1A2 CYP4F12 CYP4F11    CYP2A6    CYP2B6  CYP2D6
# Min.       1.634    0.3413   10.98 0.5232   0.05738   0.000   1.581   5.771    0.1009 1.617e-02   6.289
# 1st Qu.  128.000  162.6000 1096.00 1.7690   1.07400   5.964   7.180  18.010   63.0400 2.301e+01 104.000
# Median   251.200  389.8000 1590.00 2.2030   2.15000  25.580  14.200  21.680  284.5000 5.769e+01 179.000
# Mean     276.300  444.7000 1614.00 2.3460  12.55000  59.100  19.980  22.530  423.5000 1.183e+02 179.100
# 3rd Qu.  398.100  615.2000 2032.00 2.7670  10.62000  77.950  26.620  26.350  639.4000 1.268e+02 219.700
# Max.    1020.000 2054.0000 4045.00 5.1150 336.80000 442.900  91.630  51.350 1880.0000 1.689e+03 722.400

save(y.m, gene.tb, log.gene.tb, csv.table, file="U:/Data/Danxin/consider_gender_effects/21CYPsResults/data.Rdata")

##-------------------------- gene 1 : CYP2J2 -----------------------------------------------

log.CYP2J2 <- log(y.m[,1])
hist(log.CYP2J2, breaks=20) 


index.m <- combn(c(1:45),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})

dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP2J2, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2J2 Residual Deviance")

summary(Residual.Deviance)
summary(Residual.Deviance.Deduction)

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.199  18.100  21.010  20.920  23.900  31.500 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0005906 0.0355200 0.0647400 0.0779500 0.1096000 0.3746000 

large.interaction.index <- which((Residual.Deviance<=10) & (Residual.Deviance.Deduction>=0.30) )

goodfit.index <- which((Residual.Deviance<=9) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]

length(goodfit.index) # [1] 4

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
#      [,1]                [,2]                [,3]                [,4]               
# [1,] "34"                "19"                "27"                "34"               
# [2,] "35"                "34"                "34"                "35"               
# [3,] "38"                "35"                "36"                "36"               
# [4,] "PGRMC1"            "NCOR2"             "NR1I2"             "PGRMC1"           
# [5,] "PPARA"             "PGRMC1"            "PGRMC1"            "PPARA"            
# [6,] "RXRA"              "PPARA"             "PPARD"             "PPARD"            
# [7,] "11.2600128821749"  "11.4689041703849"  "11.6351063469732"  "11.6413525221764" 
# [8,] "8.19931950865824"  "8.88589702693438"  "8.90553187550164"  "8.91705666101246" 
# [9,] "0.271819704430526" "0.225218303778347" "0.234598154075455" "0.234018844114053"


goodfit.index <- which((Residual.Deviance<=9.5) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])

goodfitr[,order(goodfitr[8,])]
#       [,1]                [,2]                [,3]                [,4]                [,5]               
# [1,] "34"                "19"                "27"                "34"                "19"               
# [2,] "35"                "34"                "34"                "35"                "34"               
# [3,] "38"                "35"                "36"                "36"                "39"               
# [4,] "PGRMC1"            "NCOR2"             "NR1I2"             "PGRMC1"            "NCOR2"            
# [5,] "PPARA"             "PGRMC1"            "PGRMC1"            "PPARA"             "PGRMC1"           
# [6,] "RXRA"              "PPARA"             "PPARD"             "PPARD"             "RXRB"             
# [7,] "11.2600128821749"  "11.4689041703849"  "11.6351063469732"  "11.6413525221764"  "10.7641348449992" 
# [8,] "8.19931950865824"  "8.88589702693438"  "8.90553187550164"  "8.91705666101246"  "9.07783691492821" 
# [9,] "0.271819704430526" "0.225218303778347" "0.234598154075455" "0.234018844114053" "0.156658937699432"
#       [,6]               [,7]                [,8]                [,9]               [,10]              
# [1,] "13"               "22"                "24"                "19"               "25"               
# [2,] "34"               "34"                "27"                "27"               "27"               
# [3,] "36"               "36"                "34"                "34"               "34"               
# [4,] "HNF4A"            "NR0B2"             "NR1H2"             "NCOR2"            "NR1H3"            
# [5,] "PGRMC1"           "PGRMC1"            "NR1I2"             "NR1I2"            "NR1I2"            
# [6,] "PPARD"            "PPARD"             "PGRMC1"            "PGRMC1"           "PGRMC1"           
# [7,] "12.3269792577911" "11.9773845489999"  "13.2478526365686"  "12.2183166360337" "12.857310289428"  
# [8,] "9.1786311980431"  "9.27879172330136"  "9.28577332785505"  "9.34071040650534" "9.34689566391791" 
# [9,] "0.25540304675682" "0.225307354427716" "0.299073322855123" "0.23551576827219" "0.273028693131604"
#     [,11]               [,12]               [,13]               [,14]              
# [1,] "17"                "33"                "12"                "10"               
# [2,] "34"                "34"                "19"                "34"               
# [3,] "36"                "35"                "34"                "36"               
# [4,] "NCOA3"             "ONECUT1"           "FOXA3"             "FOXA1"            
# [5,] "PGRMC1"            "PGRMC1"            "NCOR2"             "PGRMC1"           
# [6,] "PPARD"             "PPARA"             "PGRMC1"            "PPARD"            
# [7,] "11.9859985931635"  "12.6220470331962"  "11.8407671148905"  "12.0949507871581" 
# [8,] "9.40687888038402"  "9.4837273070951"   "9.49459539774369"  "9.49480037165693" 
# [9,] "0.215177708618332" "0.248637936290939" "0.198143557286621" "0.214978172400824"

##------------------------------------------- log.gene.tb------------------------
log.gene.tb <- gene.tb
log.gene.tb <- log(log.gene.tb)
log.gene.tb[which(gene.tb==0,arr.ind=T)] <- -750

log.gene.tb <- gene.tb
log.gene.tb[which(gene.tb==0,arr.ind=T)] <- exp(-750)
log.gene.tb <- log(log.gene.tb)
log.gene.tb[1:3, 1:5]
#                               AHR        AHRR     ARNT      CEBPA    CEBPB
# GTEX.1192X.1026.SM.5H12P 2.816572  0.06412156 3.110143 -0.1358014 4.192281
# GTEX.11DXY.0526.SM.5EGGQ 1.647608 -1.70630681 2.775849  3.0242490 4.559045
# GTEX.11DXZ.0126.SM.5EGGY 2.679415 -0.52598833 2.737580  2.8909461 3.993440

exp(-745.0001)>0 # [1] FALSE
exp(-745)>0 # [1] TRUE

hist(gene.tb[,21], breaks=20)
hist(gene.tb[,40], breaks=20)
gene.tb <- gene.tb[,-c(21,40)]
dim(gene.tb) # [1] 137  43

log.gene.tb <- log(gene.tb)

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(log.gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP2J2, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2J2 Residual Deviance")

summary(Residual.Deviance)
summary(Residual.Deviance.Deduction)

large.interaction.index <- which((Residual.Deviance<=10) & (Residual.Deviance.Deduction>=0.30) )

goodfit.index <- which((Residual.Deviance<=8.75) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]               
# [1,] "26"                "19"                "26"                "27"                "33"               
# [2,] "33"                "33"                "31"                "30"                "34"               
# [3,] "34"                "38"                "33"                "33"                "37"               
# [4,] "NR1I2"             "NCOR2"             "NR1I2"             "NR1I3"             "PGRMC1"           
# [5,] "PGRMC1"            "PGRMC1"            "NR5A2"             "NR3C1"             "PPARA"            
# [6,] "PPARA"             "RXRB"              "PGRMC1"            "PGRMC1"            "RXRA"             
# [7,] "10.4529145185069"  "9.8967996423686"   "11.2647525460398"  "10.2503347006783"  "10.1725430126267" 
# [8,] "8.17540328406413"  "8.39387221507176"  "8.41014105874738"  "8.51932344138356"  "8.55161395420898" 
# [9,] "0.217882890978436" "0.151859942770059" "0.253410936070318" "0.168873632895147" "0.159343544323748"
# [,6]                [,7]                [,8]                [,9]                [,10]              
# [1,] "21"                "21"                "19"                "21"                "19"               
# [2,] "27"                "33"                "33"                "24"                "24"               
# [3,] "33"                "35"                "34"                "33"                "33"               
# [4,] "NR0B2"             "NR0B2"             "NCOR2"             "NR0B2"             "NCOR2"            
# [5,] "NR1I3"             "PGRMC1"            "PGRMC1"            "NR1H3"             "NR1H3"            
# [6,] "PGRMC1"            "PPARD"             "PPARA"             "PGRMC1"            "PGRMC1"           
# [7,] "10.5407155527472"  "10.9388030916448"  "10.681125947483"   "10.8642735701907"  "9.87752123892421" 
# [8,] "8.55851099816877"  "8.61828246637009"  "8.67163939858836"  "8.7305857388473"   "8.75310120257948" 
# [9,] "0.188052181529723" "0.212136611824303" "0.188134337032899" "0.196394891711658" "0.113836255994443"
# [,11]                [,12]               [,13]               [,14]               [,15]              
# [1,] "27"                 "26"                "19"                "24"                "27"               
# [2,] "33"                 "33"                "21"                "30"                "33"               
# [3,] "36"                 "38"                "33"                "33"                "34"               
# [4,] "NR1I3"              "NR1I2"             "NCOR2"             "NR1H3"             "NR1I3"            
# [5,] "PGRMC1"             "PGRMC1"            "NR0B2"             "NR3C1"             "PGRMC1"           
# [6,] "PPARG"              "RXRB"              "PGRMC1"            "PGRMC1"            "PPARA"            
# [7,] "9.44886621293677"   "10.4868066247862"  "11.6093087484901"  "10.8177761612624"  "9.99788008313004" 
# [8,] "8.79683000404175"   "8.8347647068918"   "8.85896006499287"  "8.88504196272295"  "8.88913583137568" 
# [9,] "0.0690068198872681" "0.157535270459715" "0.236908910175632" "0.178662801829865" "0.110897934615679"


##-------------------------- gene 2 : CYP1B1 -----------------------------------------------
log.CYP1B1 <- log(y.m[,2])
hist(log.CYP1B1, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP1B1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP1B1 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.96   68.88   72.67   72.69   76.40   88.87 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001338 0.0262800 0.0458000 0.0509000 0.0699800 0.2008000

large.interaction.index <- which((Residual.Deviance<=53) & (Residual.Deviance.Deduction>=0.19) )

goodfit.index <- which((Residual.Deviance<=58) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                 [,2]                 [,3]                [,4]                 [,5]                [,6]                
# [1,] "10"                 "16"                 "16"                "16"                 "15"                "8"                 
# [2,] "16"                 "25"                 "30"                "18"                 "27"                "16"                
# [3,] "42"                 "42"                 "42"                "42"                 "34"                "42"                
# [4,] "FOXA1"              "NCOA2"              "NCOA2"             "NCOA2"              "NCOA1"             "DBP"               
# [5,] "NCOA2"              "NR1H4"              "NR3C1"             "NCOR1"              "NR1I3"             "NCOA2"             
# [6,] "VDR"                "VDR"                "VDR"               "VDR"                "PPARA"             "VDR"               
# [7,] "56.7776838717444"   "56.8874227573163"   "58.0515830060861"  "58.0975907317527"   "65.1494391589049"  "59.7982453603893"  
# [8,] "51.9648945624058"   "53.6057294731687"   "54.1004811200956"  "54.622744244934"    "55.2547010028438"  "55.5358073841277"  
# [9,] "0.0847655096359732" "0.0576875014736987" "0.068061914617838" "0.0598105092320047" "0.151877564623802" "0.0712803185206002"
# [,7]                [,8]               [,9]                 [,10]                [,11]                [,12]              
# [1,] "2"                 "16"               "16"                 "16"                 "16"                 "1"                
# [2,] "16"                "35"               "27"                 "24"                 "21"                 "6"                
# [3,] "42"                "42"               "42"                 "42"                 "42"                 "27"               
# [4,] "AHRR"              "NCOA2"            "NCOA2"              "NCOA2"              "NCOA2"              "AHR"              
# [5,] "NCOA2"             "PPARD"            "NR1I3"              "NR1H3"              "NR0B2"              "CEBPD"            
# [6,] "VDR"               "VDR"              "VDR"                "VDR"                "VDR"                "NR1I3"            
# [7,] "63.4454138387781"  "62.2201710575954" "58.0686275901189"   "58.9861639784969"   "60.9274110639721"   "67.0112370059785" 
# [8,] "55.7243240059275"  "55.8732859146959" "55.881390154159"    "56.0266650405357"   "56.6657505433988"   "56.7209322150653" 
# [9,] "0.121696579243233" "0.10200687389664" "0.0376664220721483" "0.0501727649053434" "0.0699465223641081" "0.153560883975251"
# [,13]               [,14]               [,15]               [,16]                [,17]               [,18]               
# [1,] "6"                 "10"                "15"                "3"                  "15"                "16"                
# [2,] "15"                "13"                "26"                "16"                 "30"                "21"                
# [3,] "27"                "15"                "34"                "42"                 "34"                "27"                
# [4,] "CEBPD"             "FOXA1"             "NCOA1"             "ARNT"               "NCOA1"             "NCOA2"             
# [5,] "NCOA1"             "HNF4A"             "NR1I2"             "NCOA2"              "NR3C1"             "NR0B2"             
# [6,] "NR1I3"             "NCOA1"             "PPARA"             "VDR"                "PPARA"             "NR1I3"             
# [7,] "65.8501697112475"  "64.9573993355014"  "66.7558349224869"  "60.8745594088475"   "65.6065217094524"  "62.692251709453"   
# [8,] "56.8157699653445"  "56.9985005793879"  "57.0302265630037"  "57.0701524039154"   "57.091992405111"   "57.1342463848636"  
# [9,] "0.137196301627145" "0.122524898433914" "0.145689262530775" "0.0624958446003832" "0.129781751607702" "0.0886553788233349"

##------------------------------------------- log.gene.tb------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP1B1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP1B1 Residual Deviance")

summary(Residual.Deviance)
summary(Residual.Deviance.Deduction)

goodfit.index <- which((Residual.Deviance<=58) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                 [,3]                 [,4]                [,5]                [,6]               
# [1,] "10"                "16"                 "10"                 "10"                "10"                "10"               
# [2,] "13"                "25"                 "16"                 "15"                "27"                "13"               
# [3,] "15"                "27"                 "27"                 "27"                "37"                "20"               
# [4,] "FOXA1"             "NCOA2"              "FOXA1"              "FOXA1"             "FOXA1"             "FOXA1"            
# [5,] "HNF4A"             "NR1H4"              "NCOA2"              "NCOA1"             "NR1I3"             "HNF4A"            
# [6,] "NCOA1"             "NR1I3"              "NR1I3"              "NR1I3"             "RXRA"              "NFE2L2"           
# [7,] "61.3364610664631"  "60.5793393792343"   "59.9904538478305"   "63.5734149964616"  "63.6110233251787"  "62.9995585827998" 
# [8,] "53.3449486065618"  "55.4375938741546"   "55.549619125012"    "55.8270840158424"  "55.9769522574456"  "55.9862195443137" 
# [9,] "0.130289754592165" "0.0848762227810327" "0.0740256897219516" "0.121848590028558" "0.120011763192488" "0.111323621883295"
# [,7]                [,8]                 [,9]                 [,10]               [,11]               [,12]              
# [1,] "16"                "2"                  "16"                 "3"                 "4"                 "10"               
# [2,] "21"                "16"                 "21"                 "15"                "10"                "15"               
# [3,] "42"                "27"                 "27"                 "43"                "27"                "26"               
# [4,] "NCOA2"             "AHRR"               "NCOA2"              "ARNT"              "CEBPA"             "FOXA1"            
# [5,] "NR0B2"             "NCOA2"              "NR0B2"              "NCOA1"             "FOXA1"             "NCOA1"            
# [6,] "VDR"               "NR1I3"              "NR1I3"              "YY1"               "NR1I3"             "NR1I2"            
# [7,] "67.1753125488058"  "61.3318826450355"   "58.1608799591775"   "67.0849206175614"  "66.9073625631107"  "65.8964532043302" 
# [8,] "56.0662830329879"  "56.0832289802394"   "56.3150092841274"   "56.3192201391606"  "56.6056018579478"  "56.6143908028943" 
# [9,] "0.165373693017754" "0.0855778991030368" "0.0317373237190649" "0.160478694456152" "0.153970509530184" "0.140858300410407"
# [,13]               [,14]               [,15]               [,16]                [,17]               [,18]               
# [1,] "15"                "8"                 "3"                 "16"                 "10"                "10"                
# [2,] "34"                "15"                "15"                "21"                 "13"                "15"                
# [3,] "43"                "29"                "31"                "35"                 "16"                "34"                
# [4,] "NCOA1"             "DBP"               "ARNT"              "NCOA2"              "FOXA1"             "FOXA1"             
# [5,] "PPARA"             "NCOA1"             "NCOA1"             "NR0B2"              "HNF4A"             "NCOA1"             
# [6,] "YY1"               "NR2F2"             "NR5A2"             "PPARD"              "NCOA2"             "PPARA"             
# [7,] "67.2981772995524"  "66.1362227181462"  "70.8935618092645"  "60.71648278955"     "60.7681063650404"  "61.8097044169949"  
# [8,] "56.7892651048216"  "56.8453574848531"  "56.9269995447638"  "56.965573888807"    "56.9863805759568"  "56.9955543010114"  
# [9,] "0.156154484659433" "0.140480735842568" "0.197007484291409" "0.0617774404644631" "0.062232082177555" "0.0778866386984353"

##-------------------------- gene 3 : CYP2U1 -----------------------------------------------

log.CYP2U1 <- log(y.m[,3])
hist(log.CYP2U1, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP2U1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2U1 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11.22   16.95   18.29   18.23   19.54   23.49 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0008912 0.0352100 0.0617800 0.0713300 0.0969700 0.3910000 

large.interaction.index <- which((Residual.Deviance<=13) & (Residual.Deviance.Deduction>=0.35) )
goodfit.index <- which((Residual.Deviance<=12.85) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]               
# [1,] "16"                "10"                "13"                "1"                 "14"               
# [2,] "25"                "16"                "16"                "16"                "16"               
# [3,] "31"                "31"                "31"                "31"                "31"               
# [4,] "NCOA2"             "FOXA1"             "HNF4A"             "AHR"               "HNF4G"            
# [5,] "NR1H4"             "NCOA2"             "NCOA2"             "NCOA2"             "NCOA2"            
# [6,] "NR5A2"             "NR5A2"             "NR5A2"             "NR5A2"             "NR5A2"            
# [7,] "18.4270472643013"  "18.5298353257378"  "18.2941532620092"  "15.9957995542092"  "17.511177535427"  
# [8,] "11.2222692872628"  "11.6773814352897"  "11.856963036372"   "11.9197461040861"  "12.015458630334"  
# [9,] "0.390989281880026" "0.369806518513961" "0.351871449497753" "0.254820238044964" "0.313840625164959"
# [,6]                [,7]               [,8]                [,9]                [,10]              
# [1,] "16"                "16"               "16"                "1"                 "16"               
# [2,] "31"                "29"               "31"                "5"                 "31"               
# [3,] "37"                "41"               "41"                "16"                "42"               
# [4,] "NCOA2"             "NCOA2"            "NCOA2"             "AHR"               "NCOA2"            
# [5,] "NR5A2"             "NR2F2"            "NR5A2"             "CEBPB"             "NR5A2"            
# [6,] "RXRA"              "USF1"             "USF1"              "NCOA2"             "VDR"              
# [7,] "17.830666495482"   "15.4153326290959" "16.0648927664291"  "16.926174507308"   "18.5270220910963" 
# [8,] "12.2734211974488"  "12.3709095599209" "12.5364008824923"  "12.5847163379536"  "12.667247193796"  
# [9,] "0.311667839193857" "0.19749318048634" "0.219639927588582" "0.256493761628174" "0.316282609719371"
# [,11]               [,12]               [,13]               [,14]               [,15]              
# [1,] "16"                "1"                 "1"                 "16"                "3"                
# [2,] "29"                "6"                 "16"                "31"                "16"               
# [3,] "31"                "16"                "29"                "40"                "31"               
# [4,] "NCOA2"             "AHR"               "AHR"               "NCOA2"             "ARNT"             
# [5,] "NR2F2"             "CEBPD"             "NCOA2"             "NR5A2"             "NCOA2"            
# [6,] "NR5A2"             "NCOA2"             "NR2F2"             "THRB"              "NR5A2"            
# [7,] "17.1259270347792"  "16.8624806371205"  "16.360391252987"   "15.6656047905574"  "18.0094233067053" 
# [8,] "12.6758411783039"  "12.7099783749578"  "12.7465768117747"  "12.7659633667034"  "12.8350565433762" 
# [9,] "0.259844961819477" "0.246256903211594" "0.220888020667138" "0.185096040824533" "0.287314406197703"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP2U1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2U1 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10.92   16.25   17.84   17.81   19.42   23.47 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003415 0.0334000 0.0576900 0.0671100 0.0899300 0.4195000

large.interaction.index <- which((Residual.Deviance<=13) & (Residual.Deviance.Deduction>=0.37) )
goodfit.index <- which((Residual.Deviance<=12.15) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]               [,4]                [,5]               
# [1,] "16"                "10"                "13"               "13"                "16"               
# [2,] "29"                "16"                "16"               "16"                "25"               
# [3,] "32"                "31"                "31"               "29"                "31"               
# [4,] "NCOA2"             "FOXA1"             "HNF4A"            "HNF4A"             "NCOA2"            
# [5,] "NR2F2"             "NCOA2"             "NCOA2"            "NCOA2"             "NR1H4"            
# [6,] "ONECUT1"           "NR5A2"             "NR5A2"            "NR2F2"             "NR5A2"            
# [7,] "14.4606535952933"  "19.2152802438602"  "18.5569289811173" "14.2878828875629"  "19.0229193943342" 
# [8,] "10.9154421747873"  "11.1536654709923"  "11.2388824142503" "11.3441474540021"  "11.4220768761037" 
# [9,] "0.245162599127601" "0.419541878679796" "0.39435655405663" "0.206030204525487" "0.399562357421032"
# [,6]                [,7]                [,8]                [,9]                [,10]              
# [1,] "16"                "14"                "16"                "10"                "24"               
# [2,] "24"                "16"                "29"                "16"                "29"               
# [3,] "29"                "29"                "37"                "29"                "31"               
# [4,] "NCOA2"             "HNF4G"             "NCOA2"             "FOXA1"             "NR1H3"            
# [5,] "NR1H3"             "NCOA2"             "NR2F2"             "NCOA2"             "NR2F2"            
# [6,] "NR2F2"             "NR2F2"             "RXRA"              "NR2F2"             "NR5A2"            
# [7,] "13.2268558729581"  "14.2131987635538"  "14.4361731959835"  "14.5247046259724"  "14.3964940752935" 
# [8,] "11.5075067099727"  "11.7953318876118"  "11.8018807678019"  "11.8070795980663"  "11.8907442941973" 
# [9,] "0.129989256668358" "0.170114195696891" "0.182478582960929" "0.187103634661634" "0.174052777571486"
# [,11]                [,12]               [,13]               [,14]                [,15]              
# [1,] "16"                 "16"                "1"                 "4"                  "13"               
# [2,] "29"                 "27"                "16"                "16"                 "16"               
# [3,] "41"                 "29"                "29"                "29"                 "24"               
# [4,] "NCOA2"              "NCOA2"             "AHR"               "CEBPA"              "HNF4A"            
# [5,] "NR2F2"              "NR1I3"             "NCOA2"             "NCOA2"              "NCOA2"            
# [6,] "USF1"               "NR2F2"             "NR2F2"             "NR2F2"              "NR1H3"            
# [7,] "12.7239558121541"   "13.4923017257435"  "13.7058699051026"  "13.2823861778728"   "16.6705621793128" 
# [8,] "11.9542988553304"   "12.0483792414818"  "12.0799590115298"  "12.0848050543371"   "12.1499255205427" 
# [9,] "0.0604888108844646" "0.107018247413385" "0.118628799545769" "0.0901631007785876" "0.271174817630323"

##-------------------------- gene 4 : CYP51A1 -----------------------------------------------

log.CYP51A1 <- log(y.m[,4])
hist(log.CYP51A1, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP51A1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP51A1 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 74.99  125.90  141.60  138.90  153.50  174.20 
summary(Residual.Deviance.Deduction)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000262 0.022870 0.040420 0.049180 0.067470 0.277100

large.interaction.index <- which((Residual.Deviance<=76) & (Residual.Deviance.Deduction>=0.25) )
goodfit.index <- which((Residual.Deviance<=86.45) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                 [,4]                 [,5]                
# [1,] "25"                "24"                "5"                  "24"                 "24"                
# [2,] "27"                "25"                "24"                 "33"                 "33"                
# [3,] "41"                "41"                "33"                 "41"                 "36"                
# [4,] "NR1H4"             "NR1H3"             "CEBPB"              "NR1H3"              "NR1H3"             
# [5,] "NR1I3"             "NR1H4"             "NR1H3"              "PGRMC1"             "PGRMC1"            
# [6,] "USF1"              "USF1"              "PGRMC1"             "USF1"               "PPARG"             
# [7,] "93.1101882512362"  "91.222640530564"   "89.240271764333"    "92.0273315305412"   "85.3158323672242"  
# [8,] "74.9885292822944"  "75.8618489873061"  "81.8401765289599"   "83.1906341999495"   "83.7281177649804"  
# [9,] "0.194625951351797" "0.168387929289454" "0.0829232709523263" "0.0960225313895914" "0.0186098471783036"
# [,6]                 [,7]                [,8]                [,9]                 [,10]               
# [1,] "6"                  "9"                 "25"                "9"                  "24"                
# [2,] "24"                 "27"                "40"                "24"                 "33"                
# [3,] "33"                 "33"                "41"                "33"                 "43"                
# [4,] "CEBPD"              "ESR1"              "NR1H4"             "ESR1"               "NR1H3"             
# [5,] "NR1H3"              "NR1I3"             "THRB"              "NR1H3"              "PGRMC1"            
# [6,] "PGRMC1"             "PGRMC1"            "USF1"              "PGRMC1"             "YY1"               
# [7,] "92.0139534980167"   "95.1521869340548"  "116.037437273444"  "87.9245162276148"   "88.6423958052026"  
# [8,] "83.8249703112054"   "83.8700314407914"  "83.8890812321849"  "83.907097091145"    "84.5830162276058"  
# [9,] "0.0889971887468977" "0.118569586856501" "0.277051586079938" "0.0456916831486422" "0.0457950119773111"
# [,11]               [,12]               [,13]               [,14]                [,15]               
# [1,] "8"                 "25"                "5"                 "24"                 "24"                
# [2,] "27"                "26"                "24"                "25"                 "33"                
# [3,] "33"                "41"                "25"                "33"                 "37"                
# [4,] "DBP"               "NR1H4"             "CEBPB"             "NR1H3"              "NR1H3"             
# [5,] "NR1I3"             "NR1I2"             "NR1H3"             "NR1H4"              "PGRMC1"            
# [6,] "PGRMC1"            "USF1"              "NR1H4"             "PGRMC1"             "RXRA"              
# [7,] "95.8938310352443"  "109.925590252081"  "96.6680597724405"  "89.787775448199"    "86.5997011721304"  
# [8,] "84.7236142996392"  "85.1663317883616"  "85.5106732465394"  "86.1065892765016"   "86.2466729005929"  
# [9,] "0.116485248477555" "0.225236529610089" "0.115419576560923" "0.0409987456902883" "0.0040765529991348"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP51A1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP51A1 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 70.97  122.30  139.10  136.60  152.30  173.40 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003474 0.0229400 0.0393700 0.0467200 0.0632500 0.2296000 

goodfit.index <- which((Residual.Deviance<=83.45) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]                
# [1,] "25"                "24"                "27"                "25"                "27"                
# [2,] "27"                "25"                "32"                "35"                "33"                
# [3,] "41"                "41"                "41"                "41"                "41"                
# [4,] "NR1H4"             "NR1H3"             "NR1I3"             "NR1H4"             "NR1I3"             
# [5,] "NR1I3"             "NR1H4"             "ONECUT1"           "PPARD"             "PGRMC1"            
# [6,] "USF1"              "USF1"              "USF1"              "USF1"              "USF1"              
# [7,] "87.327926235946"   "89.9213961124129"  "92.4191448288906"  "98.5160032445749"  "87.0844930044211"  
# [8,] "70.967727312081"   "72.4426848026814"  "76.7279156821123"  "79.3293468244404"  "79.695895125139"   
# [9,] "0.187342121003336" "0.194377668334697" "0.169783318984717" "0.194756748022977" "0.0848440132608571"
# [,6]               [,7]                [,8]                [,9]                [,10]              
# [1,] "8"                "24"                "9"                 "8"                 "8"                
# [2,] "31"               "33"                "25"                "40"                "32"               
# [3,] "41"               "41"                "41"                "41"                "41"               
# [4,] "DBP"              "NR1H3"             "ESR1"              "DBP"               "DBP"              
# [5,] "NR5A2"            "PGRMC1"            "NR1H4"             "THRB"              "ONECUT1"          
# [6,] "USF1"             "USF1"              "USF1"              "USF1"              "USF1"             
# [7,] "102.782178623978" "93.7137053930226"  "98.3803150465164"  "103.11685338399"   "98.4597366278861" 
# [8,] "79.9105604208387" "80.3105512445131"  "80.5161069389454"  "81.195830985319"   "81.5371890807384" 
# [9,] "0.22252513528453" "0.143022347609653" "0.181583156133668" "0.212584283551017" "0.171872768775565"
# [,11]               [,12]              [,13]               [,14]               [,15]               
# [1,] "5"                 "8"                "27"                "10"                "8"                 
# [2,] "24"                "25"               "31"                "27"                "36"                
# [3,] "33"                "41"               "41"                "41"                "41"                
# [4,] "CEBPB"             "DBP"              "NR1I3"             "FOXA1"             "DBP"               
# [5,] "NR1H3"             "NR1H4"            "NR5A2"             "NR1I3"             "PPARG"             
# [6,] "PGRMC1"            "USF1"             "USF1"              "USF1"              "USF1"              
# [7,] "92.956490918226"   "97.3283577456726" "96.5995649744729"  "94.1021980711291"  "91.4717355524245"  
# [8,] "81.8690603640992"  "82.1628185145717" "82.3676642541808"  "83.3577896283856"  "83.3735814221267"  
# [9,] "0.119275485171664" "0.15581829984977" "0.147328828282539" "0.114178081521774" "0.0885317642809624"

##-------------------------- gene 5 : CYP3A5 -----------------------------------------------
names(y.m)[5]

log.CYP3A5 <- log(y.m[,5])
hist(log.CYP3A5, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP3A5, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP3A5 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 38.02   60.69   67.34   66.82   73.43   89.42 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0004957 0.0295000 0.0496900 0.0565400 0.0767700 0.2515000 

large.interaction.index <- which((Residual.Deviance<=40) & (Residual.Deviance.Deduction>=0.22) )
goodfit.index <- which((Residual.Deviance<=43) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]               
# [1,] "16"                "16"                "16"                "12"                "16"               
# [2,] "25"                "26"                "26"                "27"                "25"               
# [3,] "26"                "37"                "27"                "32"                "27"               
# [4,] "NCOA2"             "NCOA2"             "NCOA2"             "FOXA3"             "NCOA2"            
# [5,] "NR1H4"             "NR1I2"             "NR1I2"             "NR1I3"             "NR1H4"            
# [6,] "NR1I2"             "RXRA"              "NR1I3"             "ONECUT1"           "NR1I3"            
# [7,] "49.3904011907123"  "47.1940025547328"  "50.0493846920179"  "51.4013807930181"  "46.2792976801333" 
# [8,] "38.0174418821181"  "40.4022509899435"  "40.5795270150006"  "40.6855584108891"  "40.8452363695662" 
# [9,] "0.230266590965309" "0.143911327650429" "0.189210271720436" "0.208473434308685" "0.117418836995441"
# [,6]                [,7]                [,8]                 [,9]                [,10]              
# [1,] "9"                 "10"                "16"                 "20"                "7"                
# [2,] "16"                "16"                "26"                 "26"                "26"               
# [3,] "26"                "26"                "43"                 "32"                "32"               
# [4,] "ESR1"              "FOXA1"             "NCOA2"              "NFE2L2"            "CEBPG"            
# [5,] "NCOA2"             "NCOA2"             "NR1I2"              "NR1I2"             "NR1I2"            
# [6,] "NR1I2"             "NR1I2"             "YY1"                "ONECUT1"           "ONECUT1"          
# [7,] "51.6564111593348"  "49.9164939653989"  "46.2084497652409"   "52.5585717690582"  "55.3221273640936" 
# [8,] "41.5859842848314"  "42.0986918033085"  "42.4922590417214"   "42.5325930879088"  "42.5576781787733" 
# [9,] "0.194950184275115" "0.156617613558948" "0.0804223197791602" "0.190758202585935" "0.230729543376255"
# [,11]               [,12]               [,13]               [,14]               [,15]              
# [1,] "16"                "14"                "19"                "27"                "16"               
# [2,] "23"                "16"                "20"                "37"                "26"               
# [3,] "26"                "26"                "26"                "41"                "32"               
# [4,] "NCOA2"             "HNF4G"             "NCOR2"             "NR1I3"             "NCOA2"            
# [5,] "NR1H2"             "NCOA2"             "NFE2L2"            "RXRA"              "NR1I2"            
# [6,] "NR1I2"             "NR1I2"             "NR1I2"             "USF1"              "ONECUT1"          
# [7,] "52.1353395771465"  "48.3659331247863"  "52.5110610610009"  "51.5428002591519"  "51.561335710515"  
# [8,] "42.6680316125674"  "42.7178646580397"  "42.729360417216"   "42.9198465141711"  "42.931503560025"  
# [9,] "0.181590990705451" "0.116777824841595" "0.186278860989339" "0.167296959063642" "0.167370220952793"

##------------------------------------------- log.gene.tb ------------------------

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP3A5, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP3A5 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 39.13   58.16   66.96   65.96   73.72   89.44 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001632 0.0314400 0.0528700 0.0594800 0.0799500 0.2577000

large.interaction.index <- which((Residual.Deviance<=40) & (Residual.Deviance.Deduction>=0.22) )
goodfit.index <- which((Residual.Deviance<=41.5) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]                [,6]               
# [1,] "20"                "26"                "16"                "9"                 "26"                "11"               
# [2,] "27"                "31"                "25"                "12"                "31"                "26"               
# [3,] "32"                "37"                "27"                "26"                "41"                "31"               
# [4,] "NFE2L2"            "NR1I2"             "NCOA2"             "ESR1"              "NR1I2"             "FOXA2"            
# [5,] "NR1I3"             "NR5A2"             "NR1H4"             "FOXA3"             "NR5A2"             "NR1I2"            
# [6,] "ONECUT1"           "RXRA"              "NR1I3"             "NR1I2"             "USF1"              "NR5A2"            
# [7,] "50.3094986832575"  "44.3191723248665"  "46.3730877907552"  "47.7936120406468"  "50.5748408318114"  "50.2982834903233" 
# [8,] "39.1266804039619"  "39.2973969742901"  "39.5267728671185"  "39.7203532179055"  "40.252026586762"   "40.3458536591671" 
# [9,] "0.222280455420581" "0.113309321612913" "0.147635519862915" "0.168919202337653" "0.204109673412089" "0.197868180393688"
# [,7]                [,8]                 [,9]                [,10]               [,11]               [,12]              
# [1,] "9"                 "9"                  "13"                "3"                 "9"                 "9"                
# [2,] "14"                "26"                 "16"                "20"                "26"                "27"               
# [3,] "26"                "31"                 "26"                "27"                "40"                "40"               
# [4,] "ESR1"              "ESR1"               "HNF4A"             "ARNT"              "ESR1"              "ESR1"             
# [5,] "HNF4G"             "NR1I2"              "NCOA2"             "NFE2L2"            "NR1I2"             "NR1I3"            
# [6,] "NR1I2"             "NR5A2"              "NR1I2"             "NR1I3"             "THRB"              "THRB"             
# [7,] "46.033440390671"   "44.044260835964"    "49.6513651492515"  "51.3825692893036"  "46.9120729021286"  "49.5797535388118" 
# [8,] "40.3643869449979"  "40.3970701922748"   "40.4337104817471"  "40.5397483454995"  "40.6700007629395"  "41.0174704593912" 
# [9,] "0.123150765998841" "0.0828073981596054" "0.185647557520245" "0.211021385146292" "0.133058970815716" "0.172697169071604"
# [,13]               [,14]                [,15]               [,16]               [,17]              
# [1,] "9"                 "26"                 "9"                 "16"                "13"               
# [2,] "26"                "31"                 "16"                "26"                "26"               
# [3,] "34"                "43"                 "26"                "31"                "31"               
# [4,] "ESR1"              "NR1I2"              "ESR1"              "NCOA2"             "HNF4A"            
# [5,] "NR1I2"             "NR5A2"              "NCOA2"             "NR1I2"             "NR1I2"            
# [6,] "PPARA"             "YY1"                "NR1I2"             "NR5A2"             "NR5A2"            
# [7,] "47.5085690177624"  "44.8308467376205"   "46.6564197991455"  "49.9538857315641"  "48.1781381660854" 
# [8,] "41.2386101054092"  "41.2715888909503"   "41.2910772597211"  "41.3779640790682"  "41.4310735721838" 
# [9,] "0.131975326598639" "0.0793930542401154" "0.114996876368184" "0.171676768021211" "0.140044112344947"

##------------------------------- table partition according to the two peaks -----------------------------

index.low <- which(Residual.Deviance<=57)
par(mfrow=c(1,2))
hist(Residual.Deviance[index.low], breaks=52, ylim=c(0,300),xlim=c(39, 72)); 
hist(Residual.Deviance[-index.low], breaks=100, ylim=c(0,300));

length(index.low)/(length(Residual.Deviance[-index.low])) # [1] 0.2987792
pg.tb <- table(gene.names.m[,index.low])
npg.tb <- table(gene.names.m[,-index.low])

dim(pg.tb)
dim(npg.tb)

pg.tb[order(pg.tb)]
npg.tb[order(npg.tb)]

dim(pg.tb) # [1] 43
dim(npg.tb) # [1] 42

pg.tb[order(pg.tb)]
# NR2F1     VDR    AHRR   CEBPG   HNF4G   NR2F2   CEBPA   NCOA1   PPARG   FOXA3   NR3C1     AHR   PPARA   CEBPD    THRA   PPARD   CEBPB 
# 120     120     122     124     124     125     126     127     128     129     130     131     131     132     133     134     135 
# THRB  NFE2L2    ARNT   HNF4A   FOXA1   FOXA2   NR1H4   NR5A2  PGRMC1    RXRB   NCOA3    USF1   NCOR1   NR1D2 ONECUT1   NCOA2   NR1H2 
# 135     137     138     138     139     139     139     139     143     144     146     146     148     149     156     165     172 
# RXRA   NR0B2   NCOR2     DBP     YY1   NR1H3    ESR1   NR1I2   NR1I3 
# 184     186     190     231     232     273     857     859     861 

npg.tb[order(npg.tb)]
# NR1I2    ESR1   NR1H3     YY1     DBP   NCOR2   NR0B2    RXRA   NR1H2   NCOA2 ONECUT1   NR1D2   NCOR1   NCOA3    USF1    RXRB  PGRMC1 
# 2       4     588     629     630     671     675     677     689     696     705     712     713     715     715     717     718 
# FOXA1   FOXA2   NR1H4   NR5A2    ARNT   HNF4A  NFE2L2   CEBPB    THRB   PPARD    THRA   CEBPD     AHR   PPARA   NR3C1   FOXA3   PPARG 
# 722     722     722     722     723     723     724     726     726     727     728     729     730     730     731     732     733 
# NCOA1   CEBPA   NR2F2   CEBPG   HNF4G    AHRR   NR2F1     VDR 
# 734     735     736     737     737     739     741     741 

# There seem to have three clear predictor for CYP3A5:  ESR1 , NR1I2 ,  NR1I3. And a few other potential predictors which are also enriched in triplet models of CYP3A5 with relatively good fit: DBP , YY1 ,  NR1H3.

##-------------------------- gene 6 : CYP3A7 -----------------------------------------------

log.CYP3A7 <- log(y.m[,6])
hist(log.CYP3A7, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP3A7, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP3A7 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 163.9   209.2   218.1   217.9   227.5   257.1 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002147 0.0220900 0.0395900 0.0452200 0.0623600 0.2023000

goodfit.index <- which((Residual.Deviance<=174.05) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                 [,3]                 [,4]                 [,5]                
# [1,] "21"                "25"                 "21"                 "25"                 "25"                
# [2,] "34"                "34"                 "25"                 "33"                 "33"                
# [3,] "37"                "37"                 "33"                 "34"                 "37"                
# [4,] "NR0B2"             "NR1H4"              "NR0B2"              "NR1H4"              "NR1H4"             
# [5,] "PPARA"             "PPARA"              "NR1H4"              "PGRMC1"             "PGRMC1"            
# [6,] "RXRA"              "RXRA"               "PGRMC1"             "PPARA"              "RXRA"              
# [7,] "182.548001352624"  "173.312905103602"   "183.29979621521"    "171.421104010749"   "180.543913785068"  
# [8,] "163.919011525183"  "166.228835168036"   "166.442564932936"   "167.748129048243"   "169.928653664724"  
# [9,] "0.102049815332987" "0.0408744515091425" "0.0919653574654396" "0.0214266206235316" "0.0587960009163255"
# [,6]                [,7]               [,8]                [,9]                [,10]               
# [1,] "25"                "9"                "25"                "24"                "25"                
# [2,] "33"                "21"               "33"                "25"                "30"                
# [3,] "38"                "38"               "35"                "34"                "34"                
# [4,] "NR1H4"             "ESR1"             "NR1H4"             "NR1H3"             "NR1H4"             
# [5,] "PGRMC1"            "NR0B2"            "PGRMC1"            "NR1H4"             "NR3C1"             
# [6,] "RXRB"              "RXRB"             "PPARD"             "PPARA"             "PPARA"             
# [7,] "175.446557938051"  "211.839095220831" "178.480513671484"  "199.221425371552"  "188.529647519615"  
# [8,] "170.456524955709"  "171.368068127436" "171.407604955416"  "172.153202778262"  "172.908631948523"  
# [9,] "0.028441897299028" "0.19104607226162" "0.039628464590183" "0.135870037787388" "0.0828570772640251"
# [,11]                [,12]                [,13]                [,14]                [,15]               
# [1,] "10"                 "9"                  "25"                 "25"                 "19"                
# [2,] "25"                 "25"                 "28"                 "32"                 "25"                
# [3,] "33"                 "34"                 "33"                 "33"                 "33"                
# [4,] "FOXA1"              "ESR1"               "NR1H4"              "NR1H4"              "NCOR2"             
# [5,] "NR1H4"              "NR1H4"              "NR2F1"              "ONECUT1"            "NR1H4"             
# [6,] "PGRMC1"             "PPARA"              "PGRMC1"             "PGRMC1"             "PGRMC1"            
# [7,] "182.89058665686"    "186.471975736371"   "185.522826339865"   "184.902367682222"   "184.081593444647"  
# [8,] "172.964259729468"   "173.054588047568"   "173.886357933015"   "173.961988862503"   "174.018594759824"  
# [9,] "0.0542746737753948" "0.0719539096200293" "0.0627225697043491" "0.0591684084788001" "0.0546659690223099"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP3A7, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP3A7 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 158.6   205.4   215.2   215.4   224.7   256.6 
summary(Residual.Deviance.Deduction)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000117 0.019870 0.034320 0.041000 0.054530 0.212400 

goodfit.index <- which((Residual.Deviance<=170.25) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                 [,2]                 [,3]                [,4]                [,5]                
# [1,] "25"                 "21"                 "9"                 "19"                "13"                
# [2,] "34"                 "34"                 "21"                "21"                "34"                
# [3,] "37"                 "37"                 "38"                "29"                "37"                
# [4,] "NR1H4"              "NR0B2"              "ESR1"              "NCOR2"             "HNF4A"             
# [5,] "PPARA"              "PPARA"              "NR0B2"             "NR0B2"             "PPARA"             
# [6,] "RXRA"               "RXRA"               "RXRB"              "NR2F2"             "RXRA"              
# [7,] "174.791642056638"   "172.958845103128"   "201.656528175621"  "192.66932982476"   "182.881049902627"  
# [8,] "158.635143535741"   "161.058997655226"   "161.358395001652"  "161.921699233176"  "167.018386878062"  
# [9,] "0.0924329008572485" "0.0688016125501244" "0.199835500187095" "0.159587572238666" "0.0867375981984498"
# [,6]                [,7]                [,8]               [,9]                 [,10]             
# [1,] "10"                "13"                "13"               "34"                 "4"               
# [2,] "19"                "31"                "19"               "37"                 "19"              
# [3,] "21"                "37"                "21"               "41"                 "21"              
# [4,] "FOXA1"             "HNF4A"             "HNF4A"            "PPARA"              "CEBPA"           
# [5,] "NCOR2"             "NR5A2"             "NCOR2"            "RXRA"               "NCOR2"           
# [6,] "NR0B2"             "RXRA"              "NR0B2"            "USF1"               "NR0B2"           
# [7,] "204.580638321151"  "187.037756435038"  "206.203488833929" "182.753862492776"   "207.460189860171"
# [8,] "167.159693127259"  "167.472447072031"  "167.74063514647"  "168.798627136352"   "168.988168225806"
# [9,] "0.182915379974271" "0.104606202169682" "0.18652862715837" "0.0763608230549759" "0.18544291153062"
# [,11]               [,12]               [,13]               [,14]               [,15]              
# [1,] "19"                "19"                "19"                "9"                 "21"               
# [2,] "21"                "21"                "21"                "21"                "37"               
# [3,] "39"                "25"                "38"                "29"                "41"               
# [4,] "NCOR2"             "NCOR2"             "NCOR2"             "ESR1"              "NR0B2"            
# [5,] "NR0B2"             "NR0B2"             "NR0B2"             "NR0B2"             "RXRA"             
# [6,] "THRA"              "NR1H4"             "RXRB"              "NR2F2"             "USF1"             
# [7,] "203.126315806509"  "201.699822147564"  "196.063301864292"  "193.431747780103"  "190.034905856981" 
# [8,] "169.103001518783"  "169.348157266475"  "169.449188323142"  "169.511481838592"  "169.656637275891" 
# [9,] "0.167498308393161" "0.160395108615514" "0.135742453014338" "0.123662564268943" "0.107234344602074"

index.low <- which(Residual.Deviance<=190)
par(mfrow=c(1,2))
hist(Residual.Deviance[index.low], breaks=100, ylim=c(0,150)); 
hist(Residual.Deviance[-index.low], breaks=100, ylim=c(0,150));

pg.tb <- table(gene.names.m[,index.low])
npg.tb <- table(gene.names.m[,-index.low])

dim(pg.tb)
dim(npg.tb)

pg.tb[order(pg.tb)]
npg.tb[order(npg.tb)]

dim(pg.tb) # [1] 43
dim(npg.tb) # [1] 43

pg.tb[order(pg.tb)]
# CEBPG    AHRR   NCOA3     VDR   FOXA2   NR2F1     YY1   FOXA3   NCOR1   PPARG    ARNT   CEBPA   CEBPB    THRA   NR1D2    THRB 
# 4       6       6       6       8       8       8       9       9      10      11      11      11      11      12      12 
# AHR   NCOA1 ONECUT1   HNF4G   NR1H2   CEBPD     DBP   NCOA2   NR3C1   FOXA1   PPARD  NFE2L2   NR2F2   NR1H3   NR5A2   HNF4A 
# 13      13      13      14      14      17      18      18      23      30      30      31      33      34      36      41 
# NR1I2  PGRMC1    RXRB   NR1I3   NCOR2    USF1   PPARA   NR1H4    ESR1    RXRA   NR0B2 
# 48      49      60      69      71      73      85     110     112     119     166 

npg.tb[order(npg.tb)]
# NR0B2    RXRA    ESR1   NR1H4   PPARA    USF1   NCOR2   NR1I3    RXRB  PGRMC1   NR1I2   HNF4A   NR5A2   NR1H3   NR2F2  NFE2L2 
# 695     742     749     751     776     788     790     792     801     812     813     820     825     827     828     830 
# FOXA1   PPARD   NR3C1     DBP   NCOA2   CEBPD   HNF4G   NR1H2     AHR   NCOA1 ONECUT1   NR1D2    THRB    ARNT   CEBPA   CEBPB 
# 831     831     838     843     843     844     847     847     848     848     848     849     849     850     850     850 
# THRA   PPARG   FOXA3   NCOR1   FOXA2   NR2F1     YY1    AHRR   NCOA3     VDR   CEBPG 
# 850     851     852     852     853     853     853     855     855     855     857

length(index.low)/(length(Residual.Deviance[-index.low])) # [1] 0.04169832
sum(pg.tb)/sum(npg.tb) # [1] 0.04169832
# For CYP3A7, no clear predictors.

##-------------------------- gene 7 : CYP3A4 -----------------------------------------------
names(y.m)[7]

log.CYP3A4 <- log(y.m[,7])
hist(log.CYP3A4, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP3A4, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP3A4 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 208.6   309.6   344.9   341.4   371.8   441.3 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001631 0.0228200 0.0384400 0.0445200 0.0600300 0.2071000 

large.interaction.index <- which((Residual.Deviance<=210) & (Residual.Deviance.Deduction>=0.18) )
goodfit.index <- which((Residual.Deviance<=235) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]               [,2]                [,3]                [,4]               [,5]                [,6]               
# [1,] "30"               "9"                 "7"                 "26"               "9"                 "25"               
# [2,] "33"               "30"                "9"                 "30"               "10"                "30"               
# [3,] "36"               "36"                "36"                "41"               "36"                "33"               
# [4,] "NR3C1"            "ESR1"              "CEBPG"             "NR1I2"            "ESR1"              "NR1H4"            
# [5,] "PGRMC1"           "NR3C1"             "ESR1"              "NR3C1"            "FOXA1"             "NR3C1"            
# [6,] "PPARG"            "PPARG"             "PPARG"             "USF1"             "PPARG"             "PGRMC1"           
# [7,] "250.183297420385" "256.542074631837"  "272.897063192169"  "274.766468014465" "262.162221112819"  "254.741027191427" 
# [8,] "208.612765660851" "210.073561880869"  "216.376137467962"  "219.311586372445" "224.282753902015"  "227.645506731365" 
# [9,] "0.16616030002068" "0.181134080316667" "0.207114452105354" "0.20182550673942" "0.144488656870597" "0.106364965073728"
# [,7]                 [,8]                [,9]                [,10]              [,11]               [,12]              
# [1,] "25"                 "9"                 "9"                 "9"                "9"                 "24"               
# [2,] "27"                 "27"                "33"                "26"               "36"                "27"               
# [3,] "33"                 "36"                "36"                "36"               "37"                "30"               
# [4,] "NR1H4"              "ESR1"              "ESR1"              "ESR1"             "ESR1"              "NR1H3"            
# [5,] "NR1I3"              "NR1I3"             "PGRMC1"            "NR1I2"            "PPARG"             "NR1I3"            
# [6,] "PGRMC1"             "PPARG"             "PPARG"             "PPARG"            "RXRA"              "NR3C1"            
# [7,] "248.169441644195"   "265.739586080009"  "263.33841454834"   "266.136472314743" "266.034387436546"  "272.934448482713" 
# [8,] "227.918860438221"   "230.120055163562"  "230.626227894521"  "231.256486307265" "231.852047103107"  "232.374408276468" 
# [9,] "0.0815998177366584" "0.134039235335162" "0.124221096682475" "0.13106052584265" "0.128488428367522" "0.148607258745551"
# [,13]               [,14]              
# [1,] "11"                "3"                
# [2,] "24"                "9"                
# [3,] "27"                "36"               
# [4,] "FOXA2"             "ARNT"             
# [5,] "NR1H3"             "ESR1"             
# [6,] "NR1I3"             "PPARG"            
# [7,] "259.265072626372"  "277.285990268979" 
# [8,] "233.501317753629"  "234.452794161627" 
# [9,] "0.099372254857758" "0.154472990380086"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP3A4, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP3A4 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 204.6   297.4   346.5   338.4   373.4   440.5 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003584 0.0240500 0.0427300 0.0492300 0.0673500 0.2293000 

large.interaction.index <- which((Residual.Deviance<=215) & (Residual.Deviance.Deduction>=0.18) )
goodfit.index <- which((Residual.Deviance<=219) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                 [,3]                [,4]                [,5]                [,6]                
# [1,] "26"                "11"                 "30"                "30"                "19"                "11"                
# [2,] "30"                "23"                 "33"                "36"                "30"                "19"                
# [3,] "41"                "30"                 "36"                "41"                "36"                "30"                
# [4,] "NR1I2"             "FOXA2"              "NR3C1"             "NR3C1"             "NCOR2"             "FOXA2"             
# [5,] "NR3C1"             "NR1H2"              "PGRMC1"            "PPARG"             "NR3C1"             "NCOR2"             
# [6,] "USF1"              "NR3C1"              "PPARG"             "USF1"              "PPARG"             "NR3C1"             
# [7,] "248.246746699263"  "230.228799417419"   "243.333564930568"  "246.328992795615"  "244.534910078687"  "232.672657373352"  
# [8,] "204.617789046017"  "212.115444694724"   "212.116718689302"  "213.472740405425"  "213.531082210318"  "213.727262174147"  
# [9,] "0.175748356155098" "0.0786754514141135" "0.128288287109809" "0.133383618457984" "0.126786919129023" "0.0814251034611474"
# [,7]                 [,8]                [,9]                 [,10]                [,11]               [,12]              
# [1,] "9"                  "9"                 "26"                 "24"                 "9"                 "9"                
# [2,] "30"                 "10"                "30"                 "26"                 "26"                "26"               
# [3,] "36"                 "30"                "36"                 "30"                 "36"                "29"               
# [4,] "ESR1"               "ESR1"              "NR1I2"              "NR1H3"              "ESR1"              "ESR1"             
# [5,] "NR3C1"              "FOXA1"             "NR3C1"              "NR1I2"              "NR1I2"             "NR1I2"            
# [6,] "PPARG"              "NR3C1"             "PPARG"              "NR3C1"              "PPARG"             "NR2F2"            
# [7,] "232.731248906286"   "257.336549263378"  "239.279200793499"   "238.155151118352"   "240.218512720581"  "250.880677136908" 
# [8,] "215.283804341399"   "215.428252962785"  "216.478768339889"   "217.027266037899"   "217.144170200444"  "217.856704411095" 
# [9,] "0.0749682075221101" "0.162854038497657" "0.0952879831510627" "0.0887147936176828" "0.096055638088879" "0.131632189065689"
# [,13]               [,14]              
# [1,] "9"                 "7"                
# [2,] "10"                "9"                
# [3,] "26"                "11"               
# [4,] "ESR1"              "CEBPG"            
# [5,] "FOXA1"             "ESR1"             
# [6,] "NR1I2"             "FOXA2"            
# [7,] "254.560839262617"  "247.143662811586" 
# [8,] "218.193999443141"  "218.814678858836" 
# [9,] "0.142861093343421" "0.114625572958134"

##------------------------------- table partition according to the two peaks -----------------------------
index.low <- which(Residual.Deviance<=308)
par(mfrow=c(1,2))
hist(Residual.Deviance[index.low], breaks=100, ylim=c(0,150)); 
hist(Residual.Deviance[-index.low], breaks=100, ylim=c(0,150));

pg.tb <- table(gene.names.m[,index.low])
npg.tb <- table(gene.names.m[,-index.low])

dim(pg.tb)
dim(npg.tb)

pg.tb[order(pg.tb)]
# THRA     VDR   NCOR1 ONECUT1    AHRR   NR1D2    RXRB   FOXA3    ARNT   CEBPA   NR5A2   HNF4G   NR1H3   CEBPB   NR2F2   NCOA1   NR2F1 
# 151     157     159     164     165     166     166     167     168     168     168     170     170     171     171     172     172 
# NFE2L2   CEBPD     AHR   PPARA    USF1    THRB   NCOA2   NCOA3   PPARD     YY1   CEBPG   NCOR2   FOXA1   HNF4A   NR0B2   NR1H4    RXRA 
# 173     174     175     175     184     187     188     188     189     199     201     202     211     211     216     225     237 
# NR1H2   PPARG     DBP  PGRMC1   FOXA2   NR1I3   NR3C1    ESR1   NR1I2 
# 247     247     269     324     365     834     851     861     861 

npg.tb[order(npg.tb)]
# NR3C1   NR1I3   FOXA2  PGRMC1     DBP   NR1H2   PPARG    RXRA   NR1H4   NR0B2   FOXA1   HNF4A   NCOR2   CEBPG     YY1   PPARD   NCOA2 
# 10      27     496     537     592     614     614     624     636     645     650     650     659     660     662     672     673 
# NCOA3    THRB    USF1     AHR   PPARA   CEBPD  NFE2L2   NCOA1   NR2F1   CEBPB   NR2F2   HNF4G   NR1H3    ARNT   CEBPA   NR5A2   FOXA3 
# 673     674     677     686     686     687     688     689     689     690     690     691     691     693     693     693     694 
# NR1D2    RXRB    AHRR ONECUT1   NCOR1     VDR    THRA 
# 695     695     696     697     702     704     710

# For CYP3A4, four clear predictors: NR1I3   NR3C1    ESR1   NR1I2; three potential predictors: DBP  PGRMC1   FOXA2.

##-------------------------- gene 8 : CYP3A43 -----------------------------------------------
which(y.m[,8]==0) # [1] 64
log.CYP3A43 <- y.m[,8]
log.CYP3A43[which(log.CYP3A43==0)] <- exp(-175)
log.CYP3A43 <- log(log.CYP3A43)

log.CYP3A43 <- log(y.m[,8]); 
log.CYP3A43 <- log.CYP3A43[-which(y.m[,8]==0)]

hist(log.CYP3A43, breaks=30) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

gene.tb.tmp <- gene.tb
gene.tb <- gene.tb[-which(y.m[,8]==0), ]

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP3A43, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP3A43 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 115.3   224.0   265.6   262.5   302.6   388.2 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002112 0.0217500 0.0388000 0.0475800 0.0641000 0.2877000 

goodfit.index <- which((Residual.Deviance<=136.85) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]               
# [1,] "9"                 "9"                 "9"                 "9"                 "9"                
# [2,] "14"                "10"                "14"                "25"                "19"               
# [3,] "27"                "27"                "26"                "27"                "27"               
# [4,] "ESR1"              "ESR1"              "ESR1"              "ESR1"              "ESR1"             
# [5,] "HNF4G"             "FOXA1"             "HNF4G"             "NR1H4"             "NCOR2"            
# [6,] "NR1I3"             "NR1I3"             "NR1I2"             "NR1I3"             "NR1I3"            
# [7,] "135.232013619857"  "137.801246037878"  "139.400765833948"  "142.919188896414"  "145.920474668944" 
# [8,] "115.323760799344"  "118.046565710056"  "124.749238927675"  "127.844931546221"  "130.718677393503" 
# [9,] "0.147215531941097" "0.143356325837517" "0.105103632814509" "0.105473991747315" "0.104178644634555"
# [,6]                [,7]                [,8]                 [,9]                 [,10]               
# [1,] "25"                "9"                 "14"                 "19"                 "9"                 
# [2,] "27"                "12"                "27"                 "25"                 "27"                
# [3,] "33"                "27"                "33"                 "27"                 "38"                
# [4,] "NR1H4"             "ESR1"              "HNF4G"              "NCOR2"              "ESR1"              
# [5,] "NR1I3"             "FOXA3"             "NR1I3"              "NR1H4"              "NR1I3"             
# [6,] "PGRMC1"            "NR1I3"             "PGRMC1"             "NR1I3"              "RXRB"              
# [7,] "137.03693981363"   "155.054759925837"  "140.973845904435"   "141.288489484158"   "145.902962370503"  
# [8,] "131.378381540512"  "131.448033563596"  "132.969601878388"   "133.290734239117"   "134.961705042862"  
# [9,] "0.041292211288531" "0.152247672844951" "0.0567782199222492" "0.0566058514337636" "0.0749899601068894"
# [,11]               [,12]                [,13]                [,14]               [,15]               
# [1,] "9"                 "14"                 "9"                  "9"                 "10"                
# [2,] "27"                "27"                 "15"                 "10"                "19"                
# [3,] "36"                "41"                 "27"                 "26"                "27"                
# [4,] "ESR1"              "HNF4G"              "ESR1"               "ESR1"              "FOXA1"             
# [5,] "NR1I3"             "NR1I3"              "NCOA1"              "FOXA1"             "NCOR2"             
# [6,] "PPARG"             "USF1"               "NR1I3"              "NR1I2"             "NR1I3"             
# [7,] "153.552656365158"  "150.467022920427"   "150.608794885044"   "152.836606359709"  "139.720750220491"  
# [8,] "135.390185075639"  "135.556735829294"   "135.876784749737"   "136.049331261965"  "136.371616541215"  
# [9,] "0.118281713383894" "0.0990933880510059" "0.0978164000751135" "0.109838051875046" "0.0239701953646115"

##------------------------------------------- log.gene.tb ------------------------
log.gene.tb.tmp <- log.gene.tb
log.gene.tb <- log.gene.tb[-which(y.m[,8]==0),]

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP3A43, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=150,main="log.CYP3A43 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 107.3   215.9   262.9   258.8   304.3   390.3 
summary(Residual.Deviance.Deduction)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000466 0.025860 0.044120 0.050840 0.069190 0.263800 

goodfit.index <- which((Residual.Deviance<=126.6) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                 [,3]                 [,4]                 [,5]               
# [1,] "9"                 "9"                  "9"                  "9"                  "9"                
# [2,] "14"                "14"                 "10"                 "10"                 "26"               
# [3,] "26"                "27"                 "27"                 "26"                 "40"               
# [4,] "ESR1"              "ESR1"               "ESR1"               "ESR1"               "ESR1"             
# [5,] "HNF4G"             "HNF4G"              "FOXA1"              "FOXA1"              "NR1I2"            
# [6,] "NR1I2"             "NR1I3"              "NR1I3"              "NR1I2"              "THRB"             
# [7,] "122.015714621123"  "121.112146206021"   "124.915755199661"   "130.138861261476"   "138.573392914908" 
# [8,] "107.276406927142"  "112.794077699632"   "114.684954481833"   "117.326472306232"   "121.158753281663" 
# [9,] "0.120798437641814" "0.0686807126036691" "0.0819016040168423" "0.0984516756259392" "0.125670875677692"
# [,6]                 [,7]                [,8]                 [,9]                [,10]              
# [1,] "9"                  "9"                 "9"                  "9"                 "9"                
# [2,] "25"                 "12"                "26"                 "12"                "26"               
# [3,] "27"                 "26"                "37"                 "27"                "43"               
# [4,] "ESR1"               "ESR1"              "ESR1"               "ESR1"              "ESR1"             
# [5,] "NR1H4"              "FOXA3"             "NR1I2"              "FOXA3"             "NR1I2"            
# [6,] "NR1I3"              "NR1I2"             "RXRA"               "NR1I3"             "YY1"              
# [7,] "127.122210751185"   "144.576091457944"  "133.985821689491"   "137.936070941569"  "133.714388262087" 
# [8,] "121.611125149547"   "123.087256682916"  "123.757873041795"   "124.087895760391"  "124.406272069384" 
# [9,] "0.0433526570146318" "0.148633391305081" "0.0763360519697281" "0.100395604185682" "0.069611926687044"
# [,11]                [,12]               [,13]                [,14]                [,15]              
# [1,] "9"                  "9"                 "9"                  "9"                  "9"                
# [2,] "27"                 "26"                "25"                 "21"                 "24"               
# [3,] "40"                 "36"                "26"                 "27"                 "27"               
# [4,] "ESR1"               "ESR1"              "ESR1"               "ESR1"               "ESR1"             
# [5,] "NR1I3"              "NR1I2"             "NR1H4"              "NR0B2"              "NR1H3"            
# [6,] "THRB"               "PPARG"             "NR1I2"              "NR1I3"              "NR1I3"            
# [7,] "137.289126976304"   "140.850465826202"  "134.316172914018"   "134.68118547279"    "141.098637547569" 
# [8,] "124.56551995672"    "125.777785605617"  "125.988827608056"   "126.495960553822"   "126.55474923564"  
# [9,] "0.0926774559632855" "0.107011930220978" "0.0619980835166693" "0.0607748208499503" "0.103076036485647"

gene.tb <- gene.tb.tmp
log.gene.tb <- log.gene.tb.tmp

##------------------------------- table partition according to the two peaks -----------------------------
index.low <- which(Residual.Deviance<=185)

par(mfrow=c(1,2))
hist(Residual.Deviance[index.low], breaks=38, ylim=c(0,200),xlim=c(107, 307)); 
hist(Residual.Deviance[-index.low], breaks=100, ylim=c(0,200));

length(index.low)/(length(Residual.Deviance[-index.low])) # [1] 0.1920216

pg.tb <- table(gene.names.m[,index.low])
npg.tb <- table(gene.names.m[,-index.low])

dim(pg.tb)
dim(npg.tb)

pg.tb[order(pg.tb)]
npg.tb[order(npg.tb)]

dim(pg.tb) # [1] 43
dim(npg.tb) # [1] 41

pg.tb[order(pg.tb)]
# NR2F1    THRA     VDR    AHRR  NFE2L2   NCOA2   CEBPB   FOXA3   NCOR1 ONECUT1   NR5A2    ARNT   CEBPD   FOXA2 
# 84      84      84      85      85      86      87      87      87      87      89      90      90      90 
# NCOA1   NCOA3   NR0B2    RXRB   FOXA1   NR2F2   PPARG   NR1H4   PPARA    THRB   CEBPA   CEBPG    USF1     AHR 
# 90      90      90      90      93      93      93      94      94      94      95      95      96     100 
# HNF4A   NR1D2   NR1H3   HNF4G  PGRMC1   NR1H2     YY1   PPARD    RXRA     DBP   NCOR2   NR3C1   NR1I2    ESR1 
# 100     101     101     102     102     105     111     113     124     133     138     145     345     861 
# NR1I3 
# 861 

npg.tb[order(npg.tb)]
# NR1I2   NR3C1   NCOR2     DBP    RXRA   PPARD     YY1   NR1H2   HNF4G  PGRMC1   NR1D2   NR1H3     AHR   HNF4A 
# 516     716     723     728     737     748     750     756     759     759     760     760     761     761 
# USF1   CEBPA   CEBPG   NR1H4   PPARA    THRB   FOXA1   NR2F2   PPARG    ARNT   CEBPD   FOXA2   NCOA1   NCOA3 
# 765     766     766     767     767     767     768     768     768     771     771     771     771     771 
# NR0B2    RXRB   NR5A2   CEBPB   FOXA3   NCOR1 ONECUT1   NCOA2    AHRR  NFE2L2   NR2F1    THRA     VDR 
# 771     771     772     774     774     774     774     775     776     776     777     777     777  

# For CYP3A43, two clear predictors: ESR1 and NR1I3; one potential predictor: NR1I2

##-------------------------- gene 9 : CYP7A1 -----------------------------------------------
which(y.m[,9]==0)

log.CYP7A1 <- log(y.m[,9])
hist(log.CYP7A1, breaks=20) 

log.CYP7A1 <- log.CYP7A1[-which(y.m[,9]==0)]

gene.tb.tmp <- gene.tb
log.gene.tb.tmp <- log.gene.tb

gene.tb <- gene.tb[-which(y.m[,9]==0),]
log.gene.tb <- log.gene.tb[-which(y.m[,9]==0),]

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP7A1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP7A1 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 271.9   482.4   557.4   551.4   626.0   758.8 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001749 0.0220400 0.0399600 0.0497000 0.0664600 0.2993000

goodfit.index <- which((Residual.Deviance<=312.23) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]               [,3]                [,4]                [,5]               
# [1,] "9"                 "9"                "9"                 "9"                 "8"                
# [2,] "10"                "10"               "10"                "10"                "9"                
# [3,] "18"                "23"               "27"                "26"                "10"               
# [4,] "ESR1"              "ESR1"             "ESR1"              "ESR1"              "DBP"              
# [5,] "FOXA1"             "FOXA1"            "FOXA1"             "FOXA1"             "ESR1"             
# [6,] "NCOR1"             "NR1H2"            "NR1I3"             "NR1I2"             "FOXA1"            
# [7,] "368.959799258165"  "320.719143773953" "335.312944476444"  "346.692422665362"  "364.304841012686" 
# [8,] "271.942961186376"  "281.689106069146" "286.780502336189"  "291.034120411907"  "299.632605843449" 
# [9,] "0.262946907134197" "0.12169537884622" "0.144737753014676" "0.160540867393511" "0.177522305192166"
# [,6]                [,7]                [,8]                 [,9]                 [,10]              
# [1,] "9"                 "9"                 "9"                  "9"                  "9"                
# [2,] "10"                "23"                "23"                 "23"                 "10"               
# [3,] "22"                "26"                "27"                 "30"                 "30"               
# [4,] "ESR1"              "ESR1"              "ESR1"               "ESR1"               "ESR1"             
# [5,] "FOXA1"             "NR1H2"             "NR1H2"              "NR1H2"              "FOXA1"            
# [6,] "NR1D2"             "NR1I2"             "NR1I3"              "NR3C1"              "NR3C1"            
# [7,] "378.11976615253"   "350.521198203846"  "336.954217558459"   "337.019975030638"   "366.668407047475" 
# [8,] "301.314419143379"  "303.650898381685"  "303.892247582394"   "304.49413688521"    "304.596021655027" 
# [9,] "0.203124390429693" "0.133716020778018" "0.0981200657336477" "0.0965101197413346" "0.169287520275538"
# [,11]              [,12]               [,13]               [,14]               [,15]              
# [1,] "8"                "9"                 "9"                 "9"                 "9"                
# [2,] "9"                "10"                "23"                "22"                "23"               
# [3,] "23"               "40"                "31"                "23"                "25"               
# [4,] "DBP"              "ESR1"              "ESR1"              "ESR1"              "ESR1"             
# [5,] "ESR1"             "FOXA1"             "NR1H2"             "NR1D2"             "NR1H2"            
# [6,] "NR1H2"            "THRB"              "NR5A2"             "NR1H2"             "NR1H4"            
# [7,] "347.130038853001" "365.885917642739"  "378.961773094527"  "366.433842510334"  "346.804121521316" 
# [8,] "305.940660759358" "306.47344891417"   "308.895346077191"  "308.922233947626"  "312.165208742013" 
# [9,] "0.11865691090792" "0.162379763373624" "0.184890487621449" "0.156949500539338" "0.099880337717308"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP7A1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP7A1 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 257.4   452.7   549.1   537.8   620.2   762.6 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003026 0.0253300 0.0432800 0.0519500 0.0687600 0.3435000 

goodfit.index <- which((Residual.Deviance<=286.5) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                 [,3]                [,4]               [,5]                
# [1,] "9"                 "8"                  "9"                 "6"                "6"                 
# [2,] "10"                "9"                  "10"                "9"                "8"                 
# [3,] "23"                "23"                 "27"                "37"               "9"                 
# [4,] "ESR1"              "DBP"                "ESR1"              "CEBPD"            "CEBPD"             
# [5,] "FOXA1"             "ESR1"               "FOXA1"             "ESR1"             "DBP"               
# [6,] "NR1H2"             "NR1H2"              "NR1I3"             "RXRA"             "ESR1"              
# [7,] "325.949844800767"  "299.030124247445"   "324.125463424116"  "332.019203730439" "299.647075248856"  
# [8,] "257.437293658076"  "275.876883948405"   "275.965461307993"  "277.731296256676" "278.394608753973"  
# [9,] "0.210193538164033" "0.0774277854356936" "0.148584445070598" "0.16350833585469" "0.0709249922670957"
# [,6]                [,7]               [,8]                [,9]                [,10]              
# [1,] "4"                 "6"                "9"                 "9"                 "9"                
# [2,] "9"                 "9"                "23"                "10"                "10"               
# [3,] "23"                "10"               "25"                "26"                "40"               
# [4,] "CEBPA"             "CEBPD"            "ESR1"              "ESR1"              "ESR1"             
# [5,] "ESR1"              "ESR1"             "NR1H2"             "FOXA1"             "FOXA1"            
# [6,] "NR1H2"             "FOXA1"            "NR1H4"             "NR1I2"             "THRB"             
# [7,] "332.354380903285"  "321.891416699373" "331.444496268238"  "329.18284554174"   "339.985575244409" 
# [8,] "279.672484233297"  "280.12105760974"  "281.154123157784"  "281.948061925791"  "282.357865967015" 
# [9,] "0.158511214826799" "0.12976537093763" "0.151730904198673" "0.143491024078776" "0.169500453764739"
# [,11]                [,12]               [,13]              [,14]               [,15]              
# [1,] "8"                  "4"                 "8"                "9"                 "9"                
# [2,] "9"                  "9"                 "9"                "10"                "10"               
# [3,] "10"                 "10"                "25"               "29"                "30"               
# [4,] "DBP"                "CEBPA"             "DBP"              "ESR1"              "ESR1"             
# [5,] "ESR1"               "ESR1"              "ESR1"             "FOXA1"             "FOXA1"            
# [6,] "FOXA1"              "FOXA1"             "NR1H4"            "NR2F2"             "NR3C1"            
# [7,] "306.838284284331"   "314.658275327986"  "305.741041079001" "345.051304216669"  "338.19424981285"  
# [8,] "282.71650760362"    "282.878852958964"  "282.966610398561" "284.357684522839"  "285.798343899988" 
# [9,] "0.0786139732757666" "0.100996620336446" "0.07448928217182" "0.175897378019236" "0.154928435187342"

gene.tb <- gene.tb.tmp
log.gene.tb <- log.gene.tb.tmp
## ----------------------------------------- break peaks --------------------------
##------------------------------- table partition according to the two peaks -----------------------------
index.low <- which(Residual.Deviance<=400)
index.high <- which(Residual.Deviance >=600)

par(mfrow=c(1,2))
hist(Residual.Deviance[index.low], breaks=100, ylim=c(0,120)); 
hist(Residual.Deviance[index.high], breaks=100, ylim=c(0,120));

length(index.low)/(length(Residual.Deviance[index.high]))

length(index.low)/(length(Residual.Deviance[index.high])) # [1] 0.3641318
pg.tb <- table(gene.names.m[,index.low])
npg.tb <- table(gene.names.m[,index.high])

dim(pg.tb)
dim(npg.tb)

pg.tb[order(pg.tb)]
npg.tb[order(npg.tb)]

dim(pg.tb) # [1] 43
dim(npg.tb) # [1] 37

pg.tb[order(pg.tb)]
# NFE2L2   FOXA3    AHRR    ARNT   NR2F1 ONECUT1   CEBPG   HNF4G   NR0B2   PPARG   NCOA2   PPARA    THRA     VDR 
# 50      51      52      52      52      52      53      53      55      55      56      56      56      56 
# AHR     YY1   PPARD   HNF4A   NCOR1   NR5A2    RXRB   FOXA2   NR2F2   CEBPB   FOXA1   NCOA1   NR1H3   NCOA3 
# 57      57      58      59      59      60      61      64      64      65      65      65      66      67 
# NR1H4   NCOR2  PGRMC1    THRB   NR1D2    USF1   NR3C1   CEBPA   NR1I2    RXRA   CEBPD   NR1H2   NR1I3     DBP 
# 72      74      85      85      93      94     115     121     126     130     135     182     292     329 
# ESR1 
# 861 

npg.tb[order(npg.tb)]
# RXRA   CEBPA   NR1H3   CEBPD    THRB  PGRMC1   PPARD   NR1D2   NCOA3   CEBPG   NR2F2     VDR   NR1H4   FOXA1 
# 28      49     105     111     135     148     152     219     258     261     270     283     292     304 
# YY1   NCOR2   NR0B2   NCOA1    RXRB  NFE2L2   PPARG   NCOA2     AHR   HNF4A   NR2F1    USF1   FOXA2    ARNT 
# 310     311     311     342     389     394     407     410     412     413     414     415     417     429 
# ONECUT1   FOXA3   CEBPB   NR5A2   PPARA    THRA   NCOR1   HNF4G    AHRR 
# 437     442     450     453     458     458     459     480     485 

##-------------------------- gene 10 : CYP2C18 -----------------------------------------------
which(y.m[,10]==0)

log.CYP2C18 <- log(y.m[,10])
hist(log.CYP2C18, breaks=30) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP2C18, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2C18 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11.38   17.26   19.09   19.10   20.95   26.08 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002999 0.0332900 0.0559500 0.0637700 0.0859200 0.2652000

goodfit.index <- which((Residual.Deviance<=12.45) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]               
# [1,] "26"                "21"                "12"                "21"                "13"               
# [2,] "35"                "26"                "21"                "30"                "21"               
# [3,] "41"                "35"                "35"                "35"                "35"               
# [4,] "NR1I2"             "NR0B2"             "FOXA3"             "NR0B2"             "HNF4A"            
# [5,] "PPARD"             "NR1I2"             "NR0B2"             "NR3C1"             "NR0B2"            
# [6,] "USF1"              "PPARD"             "PPARD"             "PPARD"             "PPARD"            
# [7,] "12.6713032833281"  "13.7090355028478"  "14.6220392726295"  "13.7240035657325"  "14.4170454455044" 
# [8,] "11.380876059293"   "11.6897218749451"  "11.6953981681691"  "11.7552796083365"  "12.0206765582189" 
# [9,] "0.101838555607213" "0.147298008491061" "0.200152731769686" "0.143451140038448" "0.166217752197819"
# [,6]                [,7]                [,8]                 [,9]                [,10]              
# [1,] "21"                "13"                "21"                 "21"                "15"               
# [2,] "25"                "26"                "35"                 "35"                "21"               
# [3,] "35"                "35"                "40"                 "42"                "35"               
# [4,] "NR0B2"             "HNF4A"             "NR0B2"              "NR0B2"             "NCOA1"            
# [5,] "NR1H4"             "NR1I2"             "PPARD"              "PPARD"             "NR0B2"            
# [6,] "PPARD"             "PPARD"             "THRB"               "VDR"               "PPARD"            
# [7,] "14.344550618374"   "14.3301771064387"  "13.5445857455156"   "14.1841414307895"  "14.3807793605774" 
# [8,] "12.1655140346084"  "12.2615804530814"  "12.2650391403447"   "12.2667501267629"  "12.3353507229353" 
# [9,] "0.151906925614977" "0.144352483433569" "0.0944692314118598" "0.135178524084969" "0.142233503926034"
# [,11]               [,12]               [,13]               [,14]                [,15]              
# [1,] "12"                "21"                "21"                "21"                 "6"                
# [2,] "26"                "35"                "35"                "23"                 "21"               
# [3,] "35"                "36"                "41"                "35"                 "35"               
# [4,] "FOXA3"             "NR0B2"             "NR0B2"             "NR0B2"              "CEBPD"            
# [5,] "NR1I2"             "PPARD"             "PPARD"             "NR1H2"              "NR0B2"            
# [6,] "PPARD"             "PPARG"             "USF1"              "PPARD"              "PPARD"            
# [7,] "14.5361301686131"  "14.5395936175654"  "13.1185871883654"  "13.3383126364267"   "14.0807951878556" 
# [8,] "12.3532845533337"  "12.364116260854"   "12.4105092575225"  "12.4409639214366"   "12.4419997409121" 
# [9,] "0.150166900678469" "0.149624357731926" "0.053975166736768" "0.0672760295436078" "0.116385149068635"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP2C18, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2C18 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10.96   16.61   18.34   18.60   20.46   26.18 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003567 0.0362600 0.0618500 0.0694900 0.0942500 0.2889000 

goodfit.index <- which((Residual.Deviance<=12.01) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]               [,4]                [,5]               
# [1,] "19"                "21"                "11"               "12"                "24"               
# [2,] "21"                "35"                "21"               "26"                "26"               
# [3,] "40"                "40"                "35"               "35"                "35"               
# [4,] "NCOR2"             "NR0B2"             "FOXA2"            "FOXA3"             "NR1H3"            
# [5,] "NR0B2"             "PPARD"             "NR0B2"            "NR1I2"             "NR1I2"            
# [6,] "THRB"              "THRB"              "PPARD"            "PPARD"             "PPARD"            
# [7,] "13.4030723710224"  "12.5716264999904"  "13.3700943761402" "14.8702052993363"  "15.1761708705753" 
# [8,] "10.9550102043101"  "11.2366766371059"  "11.5008129391153" "11.5384000937075"  "11.5993334890278" 
# [9,] "0.182649328373774" "0.106187521788488" "0.13981063891073" "0.224059126189567" "0.235687737839228"
# [,6]                [,7]                 [,8]                [,9]                [,10]              
# [1,] "21"                "21"                 "21"                "21"                "26"               
# [2,] "32"                "35"                 "26"                "35"                "35"               
# [3,] "35"                "41"                 "35"                "43"                "41"               
# [4,] "NR0B2"             "NR0B2"              "NR0B2"             "NR0B2"             "NR1I2"            
# [5,] "ONECUT1"           "PPARD"              "NR1I2"             "PPARD"             "PPARD"            
# [6,] "PPARD"             "USF1"               "PPARD"             "YY1"               "USF1"             
# [7,] "13.7055059525197"  "12.7578961160667"   "13.2126403056188"  "13.8252196589121"  "13.4894115967311" 
# [8,] "11.6567312772791"  "11.6865178462719"   "11.7026004753085"  "11.7240886042845"  "11.8166915828706" 
# [9,] "0.149485519348085" "0.0839776605835167" "0.114287515241603" "0.151978131737903" "0.124002444574074"
# [,11]               [,12]               [,13]               [,14]               [,15]              
# [1,] "14"                "21"                "21"                "26"                "19"               
# [2,] "21"                "30"                "35"                "40"                "21"               
# [3,] "35"                "35"                "36"                "41"                "26"               
# [4,] "HNF4G"             "NR0B2"             "NR0B2"             "NR1I2"             "NCOR2"            
# [5,] "NR0B2"             "NR3C1"             "PPARD"             "THRB"              "NR0B2"            
# [6,] "PPARD"             "PPARD"             "PPARG"             "USF1"              "NR1I2"            
# [7,] "12.9865894051044"  "13.9411167437807"  "14.0627216531111"  "13.5919234783636"  "14.6057217021403" 
# [8,] "11.8300128570642"  "11.8430851516948"  "11.8752741709673"  "12.0054483848523"  "12.0099981769443" 
# [9,] "0.089059298939994" "0.150492362315366" "0.155549369183445" "0.116721897091076" "0.177719634683688"

##-------------------------- gene 11 : CYP2C19 -----------------------------------------------
which(y.m[,11]==0)

log.CYP2C19 <- log(y.m[,11])
hist(log.CYP2C19, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP2C19, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2C19 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 88.36  135.60  151.80  149.50  164.70  199.60 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001292 0.0195200 0.0341200 0.0391600 0.0534500 0.1921000 

goodfit.index <- which((Residual.Deviance<=95.45) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                 [,4]                 [,5]               
# [1,] "20"                "9"                 "9"                  "11"                 "20"               
# [2,] "27"                "11"                "20"                 "27"                 "27"               
# [3,] "28"                "27"                "27"                 "37"                 "29"               
# [4,] "NFE2L2"            "ESR1"              "ESR1"               "FOXA2"              "NFE2L2"           
# [5,] "NR1I3"             "FOXA2"             "NFE2L2"             "NR1I3"              "NR1I3"            
# [6,] "NR2F1"             "NR1I3"             "NR1I3"              "RXRA"               "NR2F2"            
# [7,] "99.6048339974748"  "99.6392625498832"  "97.7150901219613"   "98.4868149141387"   "105.905885570326" 
# [8,] "88.3638783799347"  "88.7969266666594"  "89.4536329673922"   "91.0484183907683"   "91.2782396975022" 
# [9,] "0.112855522833612" "0.108815898529916" "0.0845463801369654" "0.0755268258990326" "0.138119291426073"
# [,6]                 [,7]                 [,8]                 [,9]                 [,10]              
# [1,] "9"                  "6"                  "27"                 "9"                  "9"                
# [2,] "27"                 "9"                  "28"                 "27"                 "27"               
# [3,] "28"                 "27"                 "37"                 "31"                 "29"               
# [4,] "ESR1"               "CEBPD"              "NR1I3"              "ESR1"               "ESR1"             
# [5,] "NR1I3"              "ESR1"               "NR2F1"              "NR1I3"              "NR1I3"            
# [6,] "NR2F1"              "NR1I3"              "RXRA"               "NR5A2"              "NR2F2"            
# [7,] "101.581006101765"   "98.3278037303915"   "101.131406188198"   "96.7079010284402"   "104.243311772946" 
# [8,] "91.5393401591822"   "91.6558613588746"   "92.8362484500233"   "93.4509979751373"   "93.6146052083304" 
# [9,] "0.0988537752079644" "0.0678540770605526" "0.0820235577733177" "0.0336777348972254" "0.101960561151075"
# [,11]               [,12]                [,13]                [,14]                [,15]               
# [1,] "11"                "9"                  "11"                 "6"                  "4"                 
# [2,] "24"                "27"                 "27"                 "27"                 "9"                 
# [3,] "27"                "36"                 "31"                 "41"                 "27"                
# [4,] "FOXA2"             "ESR1"               "FOXA2"              "CEBPD"              "CEBPA"             
# [5,] "NR1H3"             "NR1I3"              "NR1I3"              "NR1I3"              "ESR1"              
# [6,] "NR1I3"             "PPARG"              "NR5A2"              "USF1"               "NR1I3"             
# [7,] "108.064780045286"  "104.578748487077"   "103.364798578116"   "100.933718868685"   "104.042370597691"  
# [8,] "93.7104714711732"  "94.8174646302701"   "94.8380119364999"   "95.126204423828"    "95.1834586900759"  
# [9,] "0.132830590763217" "0.0933390769924296" "0.0824921710186682" "0.0575379021990883" "0.0851471554975506"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP2C19, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2C19 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 87.48  136.60  151.60  150.00  165.40  198.40 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002492 0.0234400 0.0402300 0.0469900 0.0637700 0.2245000 

goodfit.index <- which((Residual.Deviance<=99.65) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]               [,2]                [,3]                [,4]                [,5]                
# [1,] "20"               "20"                "9"                 "9"                 "9"                 
# [2,] "27"               "27"                "11"                "27"                "20"                
# [3,] "28"               "29"                "27"                "28"                "27"                
# [4,] "NFE2L2"           "NFE2L2"            "ESR1"              "ESR1"              "ESR1"              
# [5,] "NR1I3"            "NR1I3"             "FOXA2"             "NR1I3"             "NFE2L2"            
# [6,] "NR2F1"            "NR2F2"             "NR1I3"             "NR2F1"             "NR1I3"             
# [7,] "103.620455919643" "109.515032164032"  "103.894109824168"  "104.647501458686"  "102.334632223923"  
# [8,] "87.4792731627393" "90.4276723448286"  "91.3819531767113"  "92.9552252078358"  "94.7911425324947"  
# [9,] "0.15577216500013" "0.174289861784581" "0.120431819172733" "0.111730104282194" "0.0737139473460213"
# [,6]                 [,7]                [,8]                [,9]                [,10]              
# [1,] "9"                  "9"                 "9"                 "4"                 "8"                
# [2,] "27"                 "24"                "27"                "9"                 "20"               
# [3,] "31"                 "27"                "29"                "27"                "27"               
# [4,] "ESR1"               "ESR1"              "ESR1"              "CEBPA"             "DBP"              
# [5,] "NR1I3"              "NR1H3"             "NR1I3"             "ESR1"              "NFE2L2"           
# [6,] "NR5A2"              "NR1I3"             "NR2F2"             "NR1I3"             "NR1I3"            
# [7,] "103.100279460875"   "109.537088157852"  "108.656931853377"  "108.624989072742"  "111.915999034286" 
# [8,] "94.9682229190261"   "96.7338880434332"  "97.0760047861897"  "97.4103251406038"  "97.4379945472735" 
# [9,] "0.0788752133783948" "0.116884612597773" "0.106582496575684" "0.103242025871492" "0.129364922012418"
# [,11]                [,12]                [,13]                [,14]                [,15]              
# [1,] "27"                 "9"                  "9"                  "9"                  "4"                
# [2,] "28"                 "21"                 "27"                 "14"                 "27"               
# [3,] "37"                 "27"                 "34"                 "27"                 "28"               
# [4,] "NR1I3"              "ESR1"               "ESR1"               "ESR1"               "CEBPA"            
# [5,] "NR2F1"              "NR0B2"              "NR1I3"              "HNF4G"              "NR1I3"            
# [6,] "RXRA"               "NR1I3"              "PPARA"              "NR1I3"              "NR2F1"            
# [7,] "105.703054694029"   "109.526702006865"   "109.569562201971"   "107.862036393993"   "116.469266507748" 
# [8,] "98.7156792776577"   "99.1378842248172"   "99.1837706074686"   "99.2850288892569"   "99.5923495917675" 
# [9,] "0.0661038173078048" "0.0948519182235287" "0.0947871962412149" "0.0795183160960029" "0.144904466405804"

##------------------------------- table partition according to the two peaks -----------------------------
index.low <- which(Residual.Deviance<=125)
par(mfrow=c(1,2))
hist(Residual.Deviance[index.low], breaks=52, ylim=c(0,250),xlim=c(87, 160)); 
hist(Residual.Deviance[-index.low], breaks=100, ylim=c(0,250));

length(index.low)/(length(Residual.Deviance[-index.low])) # [1] 0.1819749

pg.tb <- table(gene.names.m[,index.low])
npg.tb <- table(gene.names.m[,-index.low])

dim(pg.tb)
dim(npg.tb)

pg.tb[order(pg.tb)]
npg.tb[order(npg.tb)]

pg.tb[order(pg.tb)]
# THRA    AHRR    ARNT   NR1D2   NCOA1 ONECUT1   NR0B2   PPARG   CEBPB   FOXA1   NCOA2     VDR   CEBPD   HNF4G 
# 74      77      78      78      79      81      82      82      83      83      83      83      85      86 
# NR2F2  PGRMC1   CEBPG   FOXA3   CEBPA    USF1     AHR    RXRB   NCOA3   HNF4A   NCOR1   NCOR2  NFE2L2   NR1H2 
# 86      86      87      87      88      88      89      92      93      95      95      95      95      97 
# YY1   NR1H4   NR2F1   FOXA2   NR1H3   PPARD   NR5A2    RXRA     DBP    THRB   PPARA   NR3C1   NR1I2    ESR1 
# 97     101     101     102     106     107     110     113     116     116     120     192     219     832 
# NR1I3 
# 861 

npg.tb[order(npg.tb)]
# ESR1   NR1I2   NR3C1   PPARA     DBP    THRB    RXRA   NR5A2   PPARD   NR1H3   FOXA2   NR1H4   NR2F1   NR1H2 
# 29     642     669     741     745     745     748     751     754     755     759     760     760     764 
# YY1   HNF4A   NCOR1   NCOR2  NFE2L2   NCOA3    RXRB     AHR   CEBPA    USF1   CEBPG   FOXA3   HNF4G   NR2F2 
# 764     766     766     766     766     768     769     772     773     773     774     774     775     775 
# PGRMC1   CEBPD   CEBPB   FOXA1   NCOA2     VDR   NR0B2   PPARG ONECUT1   NCOA1    ARNT   NR1D2    AHRR    THRA 
# 775     776     778     778     778     778     779     779     780     782     783     783     784     787 

gene.tb <- gene.tb.tmp

# For CYP2C19,
# clear predictor: ESR1  NR1I3  
# potential predictor: NR3C1   NR1I2

##-------------------------- gene 12 : CYP2C9 -----------------------------------------------
which(y.m[,12]==0)

dim(gene.tb)
dim(log.gene.tb)
dim(gene.tb.tmp)
dim(log.gene.tb.tmp)

gene.tb <- gene.tb.tmp
log.gene.tb <- log.gene.tb.tmp

#gene.tb <- gene.tb[-which(y.m[,9]==0),]
#log.gene.tb <- log.gene.tb[-which(y.m[,9]==0),]
log.CYP2C9 <- log(y.m[,12])
hist(log.CYP2C9, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP2C9, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2C9 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 28.51   77.51   94.80   96.86  118.00  157.40 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001014 0.0364800 0.0659600 0.0820200 0.1115000 0.4092000

goodfit.index <- which((Residual.Deviance<=37.5) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]               [,2]                [,3]                [,4]                [,5]               
# [1,] "27"               "10"                "25"                "22"                "7"                
# [2,] "33"               "27"                "27"                "27"                "27"               
# [3,] "37"               "33"                "33"                "33"                "33"               
# [4,] "NR1I3"            "FOXA1"             "NR1H4"             "NR1D2"             "CEBPG"            
# [5,] "PGRMC1"           "NR1I3"             "NR1I3"             "NR1I3"             "NR1I3"            
# [6,] "RXRA"             "PGRMC1"            "PGRMC1"            "PGRMC1"            "PGRMC1"           
# [7,] "48.2488846865023" "48.8795705071804"  "51.3973602809268"  "56.5156461074226"  "51.0365553326472" 
# [8,] "28.5076746597596" "29.2503645769412"  "33.3566687688686"  "33.5791881273189"  "34.1491779886894" 
# [9,] "0.40915370697191" "0.401583027971893" "0.351004242502956" "0.405842621643341" "0.330887875051303"
# [,6]                [,7]                [,8]                [,9]                [,10]              
# [1,] "24"                "26"                "26"                "9"                 "27"               
# [2,] "27"                "33"                "27"                "27"                "33"               
# [3,] "33"                "37"                "33"                "33"                "35"               
# [4,] "NR1H3"             "NR1I2"             "NR1I2"             "ESR1"              "NR1I3"            
# [5,] "NR1I3"             "PGRMC1"            "NR1I3"             "NR1I3"             "PGRMC1"           
# [6,] "PGRMC1"            "RXRA"              "PGRMC1"            "PGRMC1"            "PPARD"            
# [7,] "56.7813733040797"  "45.9519634417571"  "52.2834010175321"  "56.4114404352313"  "50.2120174621509" 
# [8,] "34.2875835297251"  "34.3237555082979"  "34.9533443161965"  "36.3356708624746"  "36.4198604653445" 
# [9,] "0.396147336097248" "0.253051383717209" "0.331463836783003" "0.355881172646294" "0.274678407558564"
# [,11]               [,12]               [,13]               [,14]               [,15]              
# [1,] "27"                "21"                "27"                "23"                "26"               
# [2,] "32"                "27"                "30"                "27"                "27"               
# [3,] "33"                "33"                "33"                "33"                "37"               
# [4,] "NR1I3"             "NR0B2"             "NR1I3"             "NR1H2"             "NR1I2"            
# [5,] "ONECUT1"           "NR1I3"             "NR3C1"             "NR1I3"             "NR1I3"            
# [6,] "PGRMC1"            "PGRMC1"            "PGRMC1"            "PGRMC1"            "RXRA"             
# [7,] "50.6209315622547"  "55.8612320639316"  "57.0009359730306"  "50.7549463528293"  "60.7208356729699" 
# [8,] "36.4986558753266"  "36.5303966517435"  "37.0024423644264"  "37.3013724017959"  "37.3427196763561" 
# [9,] "0.278980952169168" "0.346051003494955" "0.350845003984958" "0.265069218229672" "0.385009786797462"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP2C9, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2C9 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 24.95   69.93   92.59   92.70  116.70  158.00 
summary(Residual.Deviance.Deduction)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000195 0.037190 0.067560 0.079920 0.110300 0.363800

goodfit.index <- which((Residual.Deviance<=30.5) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]               
# [1,] "27"                "22"                "13"                "4"                 "10"               
# [2,] "33"                "27"                "27"                "27"                "27"               
# [3,] "37"                "33"                "33"                "33"                "33"               
# [4,] "NR1I3"             "NR1D2"             "HNF4A"             "CEBPA"             "FOXA1"            
# [5,] "PGRMC1"            "NR1I3"             "NR1I3"             "NR1I3"             "NR1I3"            
# [6,] "RXRA"              "PGRMC1"            "PGRMC1"            "PGRMC1"            "PGRMC1"           
# [7,] "34.0233312667165"  "38.3394561461807"  "37.312739506758"   "38.2023048254285"  "34.614034118612"  
# [8,] "24.9534370701286"  "25.730483509195"   "26.4045417498707"  "26.7462931862437"  "27.2466985150505" 
# [9,] "0.266578664078685" "0.328877190873814" "0.292345132013468" "0.299877499316726" "0.212842443568291"
# [,6]                [,7]                [,8]                [,9]                [,10]              
# [1,] "27"                "25"                "6"                 "24"                "22"               
# [2,] "33"                "27"                "27"                "27"                "23"               
# [3,] "34"                "33"                "33"                "33"                "27"               
# [4,] "NR1I3"             "NR1H4"             "CEBPD"             "NR1H3"             "NR1D2"            
# [5,] "PGRMC1"            "NR1I3"             "NR1I3"             "NR1I3"             "NR1H2"            
# [6,] "PPARA"             "PGRMC1"            "PGRMC1"            "PGRMC1"            "NR1I3"            
# [7,] "38.3891437176878"  "35.7923925593426"  "38.8930853002228"  "39.2445564894459"  "42.3177814218746" 
# [8,] "27.9581597298638"  "28.4283615660025"  "28.5148795439074"  "29.0163152429666"  "29.1260691948824" 
# [9,] "0.271717026681632" "0.205742909785391" "0.266839353993237" "0.260628279726642" "0.311729768994298"
# [,11]               [,12]               [,13]               [,14]               [,15]              
# [1,] "11"                "27"                "27"                "27"                "20"               
# [2,] "27"                "31"                "33"                "32"                "27"               
# [3,] "33"                "33"                "35"                "33"                "33"               
# [4,] "FOXA2"             "NR1I3"             "NR1I3"             "NR1I3"             "NFE2L2"           
# [5,] "NR1I3"             "NR5A2"             "PGRMC1"            "ONECUT1"           "NR1I3"            
# [6,] "PGRMC1"            "PGRMC1"            "PPARD"             "PGRMC1"            "PGRMC1"           
# [7,] "38.3565645166958"  "39.2708329401909"  "37.5957987074233"  "37.2359737726973"  "38.8256540436075" 
# [8,] "29.3086260449606"  "30.1064313675487"  "30.2123388767039"  "30.3216414440288"  "30.4652209231493" 
# [9,] "0.235890220767734" "0.233364074212521" "0.196390556513474" "0.185689579944282" "0.215332705305312"

##----------------------------------------------- split peaks ---------------------------------------------

index.low <- which(Residual.Deviance<=35)
index.low.1 <- which((Residual.Deviance>35) & ( Residual.Deviance<=58))
index.high.1 <- which((Residual.Deviance>58) & ( Residual.Deviance<=81))
index.high.2 <- which((Residual.Deviance>81) & ( Residual.Deviance<115))
index.high <- which(Residual.Deviance>=115)


pg.tb <- table(gene.names.m[,index.low])
pg.tb.1 <- table(gene.names.m[,index.low.1])
npg.tb.1 <- table(gene.names.m[,index.high.1])
npg.tb.2 <- table(gene.names.m[,index.high.2])
npg.tb <- table(gene.names.m[,index.high])

pg.tb[order(pg.tb)]
pg.tb.1[order(pg.tb.1)]
npg.tb.1[order(npg.tb.1)]
npg.tb.2[order(npg.tb.2)]
npg.tb[order(npg.tb)]

pg.tb[order(pg.tb)]
# AHR    AHRR    ARNT   CEBPA   CEBPD   CEBPG     DBP   FOXA3   HNF4G   NCOA1   NCOA2   NCOA3   NCOR1  NFE2L2 
# 1       1       1       1       1       1       1       1       1       1       1       1       1       1 
# NR1H4   NR2F1   NR2F2   NR3C1   NR5A2 ONECUT1   PPARA   PPARG    RXRB    THRB     VDR     YY1   FOXA2    THRA 
# 1       1       1       1       1       1       1       1       1       1       1       1       2       2 
# CEBPB   HNF4A   NR0B2   NR1H2   NR1H3   NR1I2   PPARD   FOXA1    USF1    ESR1   NCOR2    RXRA   NR1D2  PGRMC1 
# 3       3       3       3       3       3       3       4       4       5       6       7       8      41 
# NR1I3 
# 60 

pg.tb.1[order(pg.tb.1)]
# VDR    THRA   HNF4G   NCOA1   NCOR1   NR1H4   CEBPD   PPARG    ARNT ONECUT1   FOXA1   NCOA2   CEBPB     YY1 
# 53      58      60      60      63      63      64      65      67      67      68      69      70      73 
# CEBPA   FOXA3   NR2F1   NR2F2   NCOA3  NFE2L2    RXRB    USF1   PPARA    THRB     AHR    AHRR   NR5A2   HNF4A 
# 74      74      74      75      76      77      80      80      81      81      83      83      87      88 
# NR3C1   FOXA2   NR1D2    RXRA   NR1H2     DBP   NCOR2   CEBPG   NR0B2  PGRMC1   NR1H3   PPARD    ESR1   NR1I2 
# 88      91      92      93      96     106     122     128     143     175     183     190     226     571 
# NR1I3 
# 801 

npg.tb.1[order(npg.tb.1)]
# AHRR   FOXA2   NCOA2   NR5A2  NFE2L2   PPARG   CEBPA   HNF4G    ARNT   NR2F2    THRB   HNF4A     VDR   FOXA3 
# 84      87      97      99     100     102     103     108     110     111     117     118     118     119 
# NR2F1    THRA    USF1   CEBPD    RXRB   NCOA1   NR1H4     AHR   NCOR1   PPARA   CEBPB ONECUT1     YY1   FOXA1 
# 123     125     126     127     135     136     140     147     150     152     153     159     164     169 
# NR0B2   NCOA3   NR1D2    RXRA   NR3C1   CEBPG     DBP   NCOR2   NR1H3   NR1I2  PGRMC1   NR1H2   PPARD    ESR1 
# 181     193     202     214     238     241     248     272     281     287     375     387     479     630 

npg.tb.2[order(npg.tb.2)]
# PPARD  PGRMC1    RXRB   FOXA2  NFE2L2    AHRR   PPARG     VDR ONECUT1   HNF4G   CEBPA   NCOA2    ARNT   NR1H4 
# 189     270     277     281     281     289     289     294     296     299     301     307     310     313 
# NCOA1   NR5A2    THRA   NCOR1   HNF4A   NR2F1    USF1   CEBPB   FOXA1   FOXA3   CEBPD   NR1H2   NR2F2     AHR 
# 314     321     321     327     330     333     333     343     346     350     351     375     375     380 
# PPARA   NR1H3   NCOA3     YY1   NR0B2    THRB   NR1D2   NCOR2    RXRA   CEBPG     DBP   NR3C1 
# 389     394     398     400     422     429     436     461     472     489     506     530 

npg.tb[order(npg.tb)]
# CEBPG   NR3C1    RXRA   NR0B2   NR1D2   NCOA3     YY1    THRB   PPARA     AHR   FOXA1   CEBPB   NR2F2   FOXA3 
# 2       4      75     112     123     193     223     233     238     250     274     292     299     317 
# CEBPD    USF1   NCOR1   HNF4A   NR2F1 ONECUT1   NR1H4   NCOA1   NR5A2    THRA    RXRB    ARNT   CEBPA   NCOA2 
# 318     318     320     322     330     338     344     350     353     355     368     373     382     387 
# HNF4G     VDR   FOXA2  NFE2L2    AHRR   PPARG 
# 393     395     400     402     404     404

# clear predictor: NR1I3
# potential predictor: PGRMC1   NR1H3   PPARD    ESR1   NR1I2

##-------------------------- gene 13 : CYP2C8 -----------------------------------------------
which(y.m[,13]==0)

dim(gene.tb)
dim(log.gene.tb)
dim(gene.tb.tmp)
dim(log.gene.tb.tmp)

gene.tb <- gene.tb.tmp
log.gene.tb <- log.gene.tb.tmp

#gene.tb <- gene.tb[-which(y.m[,9]==0),]
#log.gene.tb <- log.gene.tb[-which(y.m[,9]==0),]
log.CYP2C8 <- log(y.m[,13])
hist(log.CYP2C8, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP2C8, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2C8 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 49.29  136.70  165.50  167.40  200.90  253.50 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0005116 0.0315100 0.0577900 0.0736400 0.0999800 0.4206000 

goodfit.index <- which((Residual.Deviance<=53.5) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]               
# [1,] "10"                "27"                "7"                 "7"                 "17"               
# [2,] "27"                "33"                "27"                "27"                "27"               
# [3,] "33"                "37"                "33"                "37"                "37"               
# [4,] "FOXA1"             "NR1I3"             "CEBPG"             "CEBPG"             "NCOA3"            
# [5,] "NR1I3"             "PGRMC1"            "NR1I3"             "NR1I3"             "NR1I3"            
# [6,] "PGRMC1"            "RXRA"              "PGRMC1"            "RXRA"              "RXRA"             
# [7,] "85.0602017308255"  "85.045215285818"   "88.3874428425108"  "100.875525946747"  "100.504958444259" 
# [8,] "49.2871649742532"  "50.1822807796179"  "58.3419458978132"  "62.4062431439144"  "64.3548742482965" 
# [9,] "0.420561390975496" "0.409934108450823" "0.339929473898604" "0.381353975028003" "0.359684584278615"
# [,6]                [,7]                [,8]               [,9]               [,10]              
# [1,] "9"                 "26"                "22"               "23"               "7"                
# [2,] "27"                "27"                "27"               "27"               "21"               
# [3,] "37"                "37"                "33"               "37"               "27"               
# [4,] "ESR1"              "NR1I2"             "NR1D2"            "NR1H2"            "CEBPG"            
# [5,] "NR1I3"             "NR1I3"             "NR1I3"            "NR1I3"            "NR0B2"            
# [6,] "RXRA"              "RXRA"              "PGRMC1"           "RXRA"             "NR1I3"            
# [7,] "107.715802646697"  "104.973141047671"  "101.948696776833" "105.864341448754" "109.584504355096" 
# [8,] "64.3994236461341"  "64.5449636526568"  "66.4424439376127" "67.118597983653"  "67.2238351731534" 
# [9,] "0.402135786358467" "0.385128776671119" "0.34827569122285" "0.36599428036736" "0.386557108883547"
# [,11]               [,12]               [,13]               [,14]               [,15]              
# [1,] "27"                "9"                 "26"                "25"                "27"               
# [2,] "37"                "10"                "33"                "27"                "28"               
# [3,] "42"                "27"                "37"                "33"                "37"               
# [4,] "NR1I3"             "ESR1"              "NR1I2"             "NR1H4"             "NR1I3"            
# [5,] "RXRA"              "FOXA1"             "PGRMC1"            "NR1I3"             "NR2F1"            
# [6,] "VDR"               "NR1I3"             "RXRA"              "PGRMC1"            "RXRA"             
# [7,] "103.676427097772"  "112.892118682982"  "87.469087944816"   "94.1329689275675"  "100.847646595235" 
# [8,] "67.649974613922"   "67.8593789152391"  "68.0340068666025"  "68.3346901701457"  "68.5817321363051" 
# [9,] "0.347489332843955" "0.398900652172197" "0.222193709056107" "0.274062095898333" "0.319947123688801"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP2C8, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2C8 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 44.07  126.10  162.90  161.60  199.50  256.30 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000607 0.0328700 0.0590900 0.0744700 0.1016000 0.3668000 

goodfit.index <- which((Residual.Deviance<=50) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]               [,2]                [,3]                [,4]                [,5]               
# [1,] "27"               "10"                "27"                "27"                "22"               
# [2,] "33"               "27"                "33"                "28"                "27"               
# [3,] "37"               "33"                "34"                "37"                "33"               
# [4,] "NR1I3"            "FOXA1"             "NR1I3"             "NR1I3"             "NR1D2"            
# [5,] "PGRMC1"           "NR1I3"             "PGRMC1"            "NR2F1"             "NR1I3"            
# [6,] "RXRA"             "PGRMC1"            "PPARA"             "RXRA"              "PGRMC1"           
# [7,] "63.7195924988873" "62.679722391727"   "67.9830914707746"  "72.2594770060685"  "74.1829654261618" 
# [8,] "44.0686742952599" "46.2857823087704"  "48.2231486539109"  "48.9643337336456"  "49.7064909959189" 
# [9,] "0.30839679654214" "0.261550936369821" "0.290659668299409" "0.322381841629806" "0.329947371200813"
# [,6]                [,7]                [,8]                [,9]                [,10]              
# [1,] "4"                 "13"                "21"                "27"                "27"               
# [2,] "27"                "27"                "27"                "33"                "37"               
# [3,] "33"                "33"                "37"                "40"                "38"               
# [4,] "CEBPA"             "HNF4A"             "NR0B2"             "NR1I3"             "NR1I3"            
# [5,] "NR1I3"             "NR1I3"             "NR1I3"             "PGRMC1"            "RXRA"             
# [6,] "PGRMC1"            "PGRMC1"            "RXRA"              "THRB"              "RXRB"             
# [7,] "71.2471509275337"  "67.5782696146511"  "79.0297758661351"  "67.0010627127371"  "72.799056629331"  
# [8,] "49.7342797729809"  "49.7414677307086"  "52.0859698240462"  "53.2315284170025"  "53.4606127453323" 
# [9,] "0.301947107701664" "0.263942861894106" "0.340932335272313" "0.205512177542176" "0.265641407724055"
# [,11]               [,12]               [,13]               [,14]               [,15]              
# [1,] "6"                 "9"                 "7"                 "14"                "11"               
# [2,] "27"                "27"                "27"                "27"                "27"               
# [3,] "33"                "37"                "33"                "33"                "33"               
# [4,] "CEBPD"             "ESR1"              "CEBPG"             "HNF4G"             "FOXA2"            
# [5,] "NR1I3"             "NR1I3"             "NR1I3"             "NR1I3"             "NR1I3"            
# [6,] "PGRMC1"            "RXRA"              "PGRMC1"            "PGRMC1"            "PGRMC1"           
# [7,] "74.4219953255699"  "77.4143065200301"  "72.0330740255716"  "69.3569860706415"  "73.9085862161223" 
# [8,] "53.787205564374"   "54.0873401659616"  "54.5467311744772"  "54.6131919171035"  "54.6876054235074" 
# [9,] "0.277267354508919" "0.301326297459409" "0.242754360932684" "0.212578357117785" "0.260064246614179"

##----------------------------------------------- split peaks ---------------------------------------------
index.low <- which(Residual.Deviance<=92)
index.low.1 <- which((Residual.Deviance>92) & ( Residual.Deviance<=140))
index.high.1 <- which((Residual.Deviance>140) & ( Residual.Deviance<=185))
index.high <- which(Residual.Deviance>=185)

pg.tb <- table(gene.names.m[,index.low])
pg.tb.1 <- table(gene.names.m[,index.low.1])
npg.tb.1 <- table(gene.names.m[,index.high.1])
npg.tb <- table(gene.names.m[,index.high])

pg.tb[order(pg.tb)]
pg.tb.1[order(pg.tb.1)]
npg.tb.1[order(npg.tb.1)]
npg.tb[order(npg.tb)]

pg.tb[order(pg.tb)]
# NCOA2   NR1H4     VDR    ARNT   CEBPB   CEBPD   NR1H2   PPARG    AHRR   HNF4G   NCOA1   NCOR1 ONECUT1    USF1 
# 41      41      42      43      43      43      43      43      44      44      44      44      44      44 
# CEBPA   NR2F1    THRA   FOXA3   HNF4A  NFE2L2   NCOA3   NR2F2   NR5A2   PPARA     AHR    RXRB    THRB     DBP 
# 45      45      45      47      47      47      48      48      48      48      51      52      52      53 
# NR1D2   NR3C1     YY1   FOXA1   NCOR2   NR1H3   PPARD   FOXA2    RXRA   NR0B2   CEBPG  PGRMC1    ESR1   NR1I2 
# 53      53      53      54      58      58      62      65      77      79      88      94     121     188 
# NR1I3 
# 858 

pg.tb.1[order(pg.tb.1)]
# NR1I3   FOXA2    AHRR   HNF4G   NR2F2   NR5A2   PPARG   NCOA2     VDR   CEBPA    ARNT   CEBPD    THRB    THRA 
# 3     100     109     110     112     112     114     115     116     117     119     120     120     121 
# AHR  NFE2L2   FOXA3 ONECUT1    RXRB   NCOA1   HNF4A   NR1H4   NR2F1   CEBPB   PPARA    USF1   NCOR1     YY1 
# 122     122     123     123     125     127     129     130     130     133     139     140     141     150 
# FOXA1   NR1D2   NCOA3    RXRA   NR0B2     DBP   NR3C1   NCOR2   NR1H2   CEBPG   NR1H3   PPARD  PGRMC1   NR1I2 
# 183     187     227     242     247     254     271     283     312     337     359     394     465     673 
# ESR1 
# 740 

npg.tb.1[order(npg.tb.1)]
# RXRB   NR1H4  NFE2L2   PPARG   HNF4G   NCOA2    AHRR   FOXA2    ARNT   CEBPA   NCOA1   FOXA1     VDR    USF1 
# 252     258     261     262     264     266     267     268     271     275     275     278     278     279 
# NR5A2 ONECUT1    THRA   HNF4A   NR2F2   CEBPD   NR2F1  PGRMC1   CEBPB   NR1D2   FOXA3   NCOR1   PPARA   NCOA3 
# 281     286     287     292     294     302     302     302     305     315     317     324     329     342 
# AHR    THRB     YY1   NR0B2   PPARD    RXRA   CEBPG   NR1H3   NR1H2   NCOR2   NR3C1     DBP 
# 350     353     353     369     405     425     436     444     504     508     511     549 

npg.tb[order(npg.tb)]
# NR1H2     DBP   NCOR2   NR3C1    RXRA   NR0B2   NCOA3     YY1   NR1D2    THRB     AHR   PPARA   FOXA1   NCOR1 
# 2       5      12      26     117     166     244     305     306     336     338     345     346     352 
# FOXA3   CEBPB   NR2F1   HNF4A   CEBPD    USF1   NR2F2 ONECUT1    THRA   NCOA1   NR5A2   CEBPA     VDR    ARNT 
# 374     380     384     393     396     398     407     408     408     415     420     424     425     428 
# FOXA2  NFE2L2   NR1H4    RXRB   NCOA2    AHRR   PPARG   HNF4G 
# 428     431     432     432     439     441     442     443

##-------------------------- gene 14 : CYP2E1 -----------------------------------------------
which(y.m[,14]==0)

dim(gene.tb)
dim(log.gene.tb)
dim(gene.tb.tmp)
dim(log.gene.tb.tmp)

gene.tb <- gene.tb.tmp
log.gene.tb <- log.gene.tb.tmp

#gene.tb <- gene.tb[-which(y.m[,9]==0),]
#log.gene.tb <- log.gene.tb[-which(y.m[,9]==0),]
log.CYP2E1 <- log(y.m[,14])
hist(log.CYP2E1 , breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP2E1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2E1 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 38.97   54.98   58.70   58.43   62.47   71.53 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 9.082e-05 2.216e-02 3.874e-02 4.462e-02 6.080e-02 2.600e-01 

goodfit.index <- which((Residual.Deviance<=45) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                 [,5]               
# [1,] "9"                 "9"                 "29"                "27"                 "9"                
# [2,] "28"                "29"                "33"                "28"                 "33"               
# [3,] "33"                "33"                "40"                "33"                 "39"               
# [4,] "ESR1"              "ESR1"              "NR2F2"             "NR1I3"              "ESR1"             
# [5,] "NR2F1"             "NR2F2"             "PGRMC1"            "NR2F1"              "PGRMC1"           
# [6,] "PGRMC1"            "PGRMC1"            "THRB"              "PGRMC1"             "THRA"             
# [7,] "43.7954415798674"  "45.6917882934327"  "45.3542945796603"  "43.0924855671082"   "47.4715975760725" 
# [8,] "38.9720477111186"  "40.1430082826897"  "40.3512928942994"  "40.4719912372396"   "40.4869727175785" 
# [9,] "0.110134609784733" "0.121439326802197" "0.110309326420535" "0.0608109347925115" "0.147132711244891"
# [,6]                 [,7]                 [,8]                 [,9]                [,10]               
# [1,] "24"                 "28"                 "28"                 "6"                 "17"                
# [2,] "28"                 "33"                 "33"                 "9"                 "28"                
# [3,] "33"                 "40"                 "35"                 "33"                "33"                
# [4,] "NR1H3"              "NR2F1"              "NR2F1"              "CEBPD"             "NCOA3"             
# [5,] "NR2F1"              "PGRMC1"             "PGRMC1"             "ESR1"              "NR2F1"             
# [6,] "PGRMC1"             "THRB"               "PPARD"              "PGRMC1"            "PGRMC1"            
# [7,] "42.6538486652"      "43.7446960676813"   "43.3447346065298"   "52.4725454215626"  "43.2142438107475"  
# [8,] "40.6110616875423"   "40.6702308560098"   "40.9059986426877"   "40.9981090081242"  "41.0791214323838"  
# [9,] "0.0478922076573207" "0.0702820110331718" "0.0562637189033496" "0.218675048470646" "0.0494078384829385"
# [,11]               [,12]                [,13]               [,14]                [,15]               
# [1,] "9"                 "23"                 "26"                "28"                 "28"                
# [2,] "23"                "28"                 "28"                "33"                 "31"                
# [3,] "33"                "33"                 "33"                "43"                 "33"                
# [4,] "ESR1"              "NR1H2"              "NR1I2"             "NR2F1"              "NR2F1"             
# [5,] "NR1H2"             "NR2F1"              "NR2F1"             "PGRMC1"             "NR5A2"             
# [6,] "PGRMC1"            "PGRMC1"             "PGRMC1"            "YY1"                "PGRMC1"            
# [7,] "50.1409230690871"  "42.8127683771899"   "43.7130362507418"  "43.6138628361414"   "44.5685428678058"  
# [8,] "41.1045719149363"  "41.1081509690411"   "41.1836885850611"  "41.2895140388069"   "41.3535925056438"  
# [9,] "0.180219082558572" "0.0398156314754218" "0.057862548169204" "0.0532938072022473" "0.0721349668464093"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP2E1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2E1 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 37.47   53.88   58.38   57.72   62.32   72.69 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001081 0.0288100 0.0498900 0.0550900 0.0764800 0.2481000 

goodfit.index <- which((Residual.Deviance<=45) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                 [,2]                [,3]                [,4]                [,5]               
# [1,] "9"                  "3"                 "7"                 "14"                "5"                
# [2,] "28"                 "28"                "28"                "28"                "28"               
# [3,] "33"                 "33"                "33"                "33"                "33"               
# [4,] "ESR1"               "ARNT"              "CEBPG"             "HNF4G"             "CEBPB"            
# [5,] "NR2F1"              "NR2F1"             "NR2F1"             "NR2F1"             "NR2F1"            
# [6,] "PGRMC1"             "PGRMC1"            "PGRMC1"            "PGRMC1"            "PGRMC1"           
# [7,] "41.2464156117183"   "41.9806659513701"  "42.9491071138657"  "42.7440105475673"  "43.5606485694733" 
# [8,] "37.4728345166534"   "37.661144119197"   "38.2152350460556"  "38.2667110549801"  "38.3196824324934" 
# [9,] "0.0914887036630831" "0.102893123162382" "0.110220500166855" "0.104746827338645" "0.120314235648288"
# [,6]                [,7]                 [,8]                 [,9]                 [,10]               
# [1,] "11"                "24"                 "27"                 "28"                 "17"                
# [2,] "28"                "28"                 "28"                 "33"                 "28"                
# [3,] "33"                "33"                 "33"                 "35"                 "33"                
# [4,] "FOXA2"             "NR1H3"              "NR1I3"              "NR2F1"              "NCOA3"             
# [5,] "NR2F1"             "NR2F1"              "NR2F1"              "PGRMC1"             "NR2F1"             
# [6,] "PGRMC1"            "PGRMC1"             "PGRMC1"             "PPARD"              "PGRMC1"            
# [7,] "43.5535636800331"  "42.1518697155977"   "42.4097662296929"   "41.5641093803835"   "42.247132937435"   
# [8,] "38.3690645443961"  "38.4556124248304"   "38.5326400126495"   "38.6479063528028"   "38.6530013471548"  
# [9,] "0.119037311704847" "0.0876890471456265" "0.0914205986433579" "0.0701615665788087" "0.0850739763951049"
# [,11]                [,12]               [,13]               [,14]                [,15]              
# [1,] "28"                 "28"                "21"                "29"                 "28"               
# [2,] "33"                 "33"                "28"                "33"                 "33"               
# [3,] "43"                 "34"                "33"                "35"                 "36"               
# [4,] "NR2F1"              "NR2F1"             "NR0B2"             "NR2F2"              "NR2F1"            
# [5,] "PGRMC1"             "PGRMC1"            "NR2F1"             "PGRMC1"             "PGRMC1"           
# [6,] "YY1"                "PPARA"             "PGRMC1"            "PPARD"              "PPARG"            
# [7,] "41.9105734268809"   "43.4524816844983"  "43.521833794767"   "41.8957722087777"   "43.5258135492009" 
# [8,] "38.7332087131634"   "38.8580021492362"  "39.0026780203057"  "39.0119954881685"   "39.0659790669822" 
# [9,] "0.0758129620741375" "0.105735722268337" "0.103836520211258" "0.0688321653611862" "0.102464126883634"

##-------------------------- gene 15 : CYP2R1 -----------------------------------------------
which(y.m[,15]==0)

dim(gene.tb)
dim(log.gene.tb)
dim(gene.tb.tmp)
dim(log.gene.tb.tmp)

gene.tb <- gene.tb.tmp
log.gene.tb <- log.gene.tb.tmp

#gene.tb <- gene.tb[-which(y.m[,9]==0),]
#log.gene.tb <- log.gene.tb[-which(y.m[,9]==0),]
log.CYP2R1 <- log(y.m[,1])
hist(log.CYP2R1, breaks=20) 


index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP2R1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2R1 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.199  17.850  20.680  20.580  23.480  31.010 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0005935 0.0401800 0.0699500 0.0828900 0.1145000 0.3746000 

goodfit.index <- which((Residual.Deviance<=10.5) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]               
# [1,] "33"                "19"                "26"                "33"                "19"               
# [2,] "34"                "33"                "33"                "34"                "33"               
# [3,] "37"                "34"                "35"                "35"                "38"               
# [4,] "PGRMC1"            "NCOR2"             "NR1I2"             "PGRMC1"            "NCOR2"            
# [5,] "PPARA"             "PGRMC1"            "PGRMC1"            "PPARA"             "PGRMC1"           
# [6,] "RXRA"              "PPARA"             "PPARD"             "PPARD"             "RXRB"             
# [7,] "11.2600128821749"  "11.4689041703849"  "11.6351063469732"  "11.6413525221764"  "10.7641348449992" 
# [8,] "8.19931950865824"  "8.88589702693438"  "8.90553187550164"  "8.91705666101246"  "9.07783691492821" 
# [9,] "0.271819704430526" "0.225218303778347" "0.234598154075455" "0.234018844114053" "0.156658937699432"
# [,6]               [,7]                [,8]                [,9]               [,10]              
# [1,] "13"               "21"                "23"                "19"               "24"               
# [2,] "33"               "33"                "26"                "26"               "26"               
# [3,] "35"               "35"                "33"                "33"               "33"               
# [4,] "HNF4A"            "NR0B2"             "NR1H2"             "NCOR2"            "NR1H3"            
# [5,] "PGRMC1"           "PGRMC1"            "NR1I2"             "NR1I2"            "NR1I2"            
# [6,] "PPARD"            "PPARD"             "PGRMC1"            "PGRMC1"           "PGRMC1"           
# [7,] "12.3269792577911" "11.9773845489999"  "13.2478526365686"  "12.2183166360337" "12.857310289428"  
# [8,] "9.1786311980431"  "9.27879172330136"  "9.28577332785505"  "9.34071040650534" "9.34689566391791" 
# [9,] "0.25540304675682" "0.225307354427716" "0.299073322855123" "0.23551576827219" "0.273028693131604"
# [,11]               [,12]               [,13]               [,14]               [,15]              
# [1,] "17"                "32"                "12"                "10"                "19"               
# [2,] "33"                "33"                "19"                "33"                "24"               
# [3,] "35"                "34"                "33"                "35"                "33"               
# [4,] "NCOA3"             "ONECUT1"           "FOXA3"             "FOXA1"             "NCOR2"            
# [5,] "PGRMC1"            "PGRMC1"            "NCOR2"             "PGRMC1"            "NR1H3"            
# [6,] "PPARD"             "PPARA"             "PGRMC1"            "PPARD"             "PGRMC1"           
# [7,] "11.9859985931635"  "12.6220470331962"  "11.8407671148905"  "12.0949507871581"  "11.1255484500652" 
# [8,] "9.40687888038402"  "9.4837273070951"   "9.49459539774369"  "9.49480037165693"  "9.55633006346838" 
# [9,] "0.215177708618332" "0.248637936290939" "0.198143557286621" "0.214978172400824" "0.141046384691953"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP2R1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2R1 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.175  16.640  19.510  19.640  22.730  31.270 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003764 0.0403800 0.0712900 0.0829600 0.1128000 0.4199000 

goodfit.index <- which((Residual.Deviance<=9.75) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]               
# [1,] "26"                "19"                "26"                "27"                "33"               
# [2,] "33"                "33"                "31"                "30"                "34"               
# [3,] "34"                "38"                "33"                "33"                "37"               
# [4,] "NR1I2"             "NCOR2"             "NR1I2"             "NR1I3"             "PGRMC1"           
# [5,] "PGRMC1"            "PGRMC1"            "NR5A2"             "NR3C1"             "PPARA"            
# [6,] "PPARA"             "RXRB"              "PGRMC1"            "PGRMC1"            "RXRA"             
# [7,] "10.4529145185069"  "9.8967996423686"   "11.2647525460398"  "10.2503347006783"  "10.1725430126267" 
# [8,] "8.17540328406413"  "8.39387221507176"  "8.41014105874738"  "8.51932344138356"  "8.55161395420898" 
# [9,] "0.217882890978436" "0.151859942770059" "0.253410936070318" "0.168873632895147" "0.159343544323748"
# [,6]                [,7]                [,8]                [,9]                [,10]              
# [1,] "21"                "21"                "19"                "21"                "19"               
# [2,] "27"                "33"                "33"                "24"                "24"               
# [3,] "33"                "35"                "34"                "33"                "33"               
# [4,] "NR0B2"             "NR0B2"             "NCOR2"             "NR0B2"             "NCOR2"            
# [5,] "NR1I3"             "PGRMC1"            "PGRMC1"            "NR1H3"             "NR1H3"            
# [6,] "PGRMC1"            "PPARD"             "PPARA"             "PGRMC1"            "PGRMC1"           
# [7,] "10.5407155527472"  "10.9388030916448"  "10.681125947483"   "10.8642735701907"  "9.87752123892421" 
# [8,] "8.55851099816877"  "8.61828246637009"  "8.67163939858836"  "8.7305857388473"   "8.75310120257948" 
# [9,] "0.188052181529723" "0.212136611824303" "0.188134337032899" "0.196394891711658" "0.113836255994443"
# [,11]                [,12]               [,13]               [,14]               [,15]              
# [1,] "27"                 "26"                "19"                "24"                "27"               
# [2,] "33"                 "33"                "21"                "30"                "33"               
# [3,] "36"                 "38"                "33"                "33"                "34"               
# [4,] "NR1I3"              "NR1I2"             "NCOR2"             "NR1H3"             "NR1I3"            
# [5,] "PGRMC1"             "PGRMC1"            "NR0B2"             "NR3C1"             "PGRMC1"           
# [6,] "PPARG"              "RXRB"              "PGRMC1"            "PGRMC1"            "PPARA"            
# [7,] "9.44886621293677"   "10.4868066247862"  "11.6093087484901"  "10.8177761612624"  "9.99788008313004" 
# [8,] "8.79683000404175"   "8.8347647068918"   "8.85896006499287"  "8.88504196272295"  "8.88913583137568" 
# [9,] "0.0690068198872681" "0.157535270459715" "0.236908910175632" "0.178662801829865" "0.110897934615679"

##----------------------------------------------- split peaks ---------------------------------------------
index.low <- which(Residual.Deviance<=11.5)
index.low.1 <- which((Residual.Deviance>11.5) & ( Residual.Deviance<=14))
index.high <- which(Residual.Deviance>14)

pg.tb <- table(gene.names.m[,index.low])
pg.tb.1 <- table(gene.names.m[,index.low.1])
npg.tb <- table(gene.names.m[,index.high])

pg.tb[order(pg.tb)]
pg.tb.1[order(pg.tb.1)]
npg.tb[order(npg.tb)]

pg.tb[order(pg.tb)]
# CEBPA   NCOA1  NFE2L2   NR5A2    AHRR    ARNT   CEBPD   NCOA2   NR1H4   NR2F1   NR2F2   NR3C1    THRB     VDR 
# 9       9       9       9      10      10      10      10      10      10      10      10      10      10 
# FOXA2   NCOA3   NCOR1   PPARG    USF1     AHR   CEBPG    ESR1   HNF4G    THRA   NR1D2   HNF4A ONECUT1    RXRA 
# 11      11      11      11      11      12      12      12      12      13      14      16      16      16 
# FOXA1     YY1   FOXA3   NR1H2   CEBPB     DBP   NCOR2   NR0B2   NR1H3   NR1I2   PPARA   PPARD    RXRB   NR1I3 
# 17      18      20      26      28      41      41      41      41      41      41      41      41      42 
# PGRMC1 
# 395 

pg.tb.1[order(pg.tb.1)]
# DBP   FOXA1   NR1H3    THRA     VDR   CEBPG   NCOA1   NCOA2   NR3C1   CEBPA   NR2F2   NCOR1   NR1H2     YY1 
# 7      27      29      30      31      32      32      32      32      33      33      34      34      34 
# CEBPD   FOXA3   NR0B2   NR2F1    RXRA   NCOA3  NFE2L2   NR5A2   FOXA2    THRB    AHRR   HNF4A   PPARG    ARNT 
# 35      35      35      35      35      36      36      36      37      37      38      38      39      40 
# PPARA   HNF4G ONECUT1    USF1     AHR   NR1H4   NR1D2   NR1I2    ESR1   PPARD   NR1I3   CEBPB   NCOR2  PGRMC1 
# 40      43      43      44      45      51      56      60      66      70      71      79      79     466 

npg.tb[order(npg.tb)]
# NCOR2   NR1I3   PPARD   CEBPB   NR1I2   PPARA    ESR1   NR0B2   NR1D2   NR1H3   NR1H4   NR1H2 ONECUT1     AHR 
# 741     748     750     754     760     780     783     785     791     791     800     801     802     804 
# FOXA3   HNF4G    USF1   HNF4A     YY1    RXRA    ARNT   PPARG    AHRR     DBP   FOXA2   NCOA3    THRB   CEBPD 
# 806     806     806     807     809     810     811     811     813     813     813     814     814     816 
# NCOR1  NFE2L2   NR2F1   NR5A2   CEBPG   FOXA1   NR2F2    THRA   CEBPA   NCOA2   NR3C1   NCOA1    RXRB     VDR 
# 816     816     816     816     817     817     818     818     819     819     819     820     820     820 

##-------------------------- gene 16 : CYP1A1 -----------------------------------------------
which(y.m[,16]==0)

dim(gene.tb)
dim(log.gene.tb)
dim(gene.tb.tmp)
dim(log.gene.tb.tmp)

gene.tb <- gene.tb.tmp
log.gene.tb <- log.gene.tb.tmp

#gene.tb <- gene.tb[-which(y.m[,9]==0),]
#log.gene.tb <- log.gene.tb[-which(y.m[,9]==0),]
log.CYP1A1 <- log(y.m[,16])
hist(log.CYP1A1, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP1A1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP1A1 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 197.4   280.8   301.2   301.5   327.7   368.5 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002219 0.0168600 0.0292600 0.0351700 0.0479600 0.2157000

goodfit.index <- which((Residual.Deviance<=20) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                 [,4]                [,5]               
# [1,] "3"                 "3"                 "3"                  "3"                 "3"                
# [2,] "9"                 "24"                "7"                  "17"                "27"               
# [3,] "27"                "27"                "27"                 "27"                "43"               
# [4,] "ARNT"              "ARNT"              "ARNT"               "ARNT"              "ARNT"             
# [5,] "ESR1"              "NR1H3"             "CEBPG"              "NCOA3"             "NR1I3"            
# [6,] "NR1I3"             "NR1I3"             "NR1I3"              "NR1I3"             "YY1"              
# [7,] "231.059084811556"  "232.244419616611"  "223.435161011069"   "230.092635039346"  "233.9323152185"   
# [8,] "197.393004229136"  "201.21944686683"   "202.366218581992"   "202.980785383919"  "204.011918601327" 
# [9,] "0.145703340813787" "0.133587591904243" "0.0942955546196868" "0.117830149803756" "0.127901938598034"
# [,6]                [,7]                 [,8]                 [,9]                 [,10]               
# [1,] "6"                 "9"                  "3"                  "3"                  "6"                 
# [2,] "9"                 "14"                 "19"                 "27"                 "23"                
# [3,] "27"                "27"                 "27"                 "32"                 "27"                
# [4,] "CEBPD"             "ESR1"               "ARNT"               "ARNT"               "CEBPD"             
# [5,] "ESR1"              "HNF4G"              "NCOR2"              "NR1I3"              "NR1H2"             
# [6,] "NR1I3"             "NR1I3"              "NR1I3"              "ONECUT1"            "NR1I3"             
# [7,] "246.494325225673"  "228.213521072576"   "225.154101270635"   "231.944530031842"   "228.542063173189"  
# [8,] "205.753110433938"  "206.45854758501"    "206.60400642043"    "210.12384393801"    "211.663493056796"  
# [9,] "0.165282566868163" "0.0953272767771162" "0.0823884386094648" "0.0940771747919068" "0.0738532324511405"
# [,11]                [,12]                [,13]                [,14]                [,15]              
# [1,] "3"                  "7"                  "3"                  "14"                 "16"               
# [2,] "27"                 "14"                 "27"                 "27"                 "27"               
# [3,] "36"                 "27"                 "41"                 "32"                 "32"               
# [4,] "ARNT"               "CEBPG"              "ARNT"               "HNF4G"              "NCOA2"            
# [5,] "NR1I3"              "HNF4G"              "NR1I3"              "NR1I3"              "NR1I3"            
# [6,] "PPARG"              "NR1I3"              "USF1"               "ONECUT1"            "ONECUT1"          
# [7,] "225.608767350848"   "225.595001275707"   "228.623771622557"   "228.012925284145"   "238.85058045031"  
# [8,] "211.920249147446"   "212.336268371331"   "212.63099788601"    "212.709471445688"   "212.872744474187" 
# [9,] "0.0606736979423963" "0.0587722814308848" "0.0699523659462262" "0.0671166067422975" "0.108761870819599"

##------------------------------------------- log.gene.tb ------------------------

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP1A1, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP1A1 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 200.2   276.1   296.6   299.0   328.5   368.5 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003279 0.0181900 0.0307200 0.0351200 0.0474600 0.1821000

goodfit.index <- which((Residual.Deviance<=205) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                 [,3]                [,4]                 [,5]                
# [1,] "3"                 "3"                  "3"                 "3"                  "6"                 
# [2,] "9"                 "27"                 "24"                "27"                 "14"                
# [3,] "27"                "36"                 "27"                "29"                 "27"                
# [4,] "ARNT"              "ARNT"               "ARNT"              "ARNT"               "CEBPD"             
# [5,] "ESR1"              "NR1I3"              "NR1H3"             "NR1I3"              "HNF4G"             
# [6,] "NR1I3"             "PPARG"              "NR1I3"             "NR2F2"              "NR1I3"             
# [7,] "223.104204050725"  "222.767145648376"   "226.583825085131"  "227.641670049178"   "222.60873882047"   
# [8,] "200.207303091317"  "201.002611783145"   "202.68769703117"   "205.350950261599"   "205.384433156206"  
# [9,] "0.102628729282939" "0.0977008247867259" "0.105462638584123" "0.0979202084695842" "0.0773747955966596"
# [,6]                [,7]                 [,8]                 [,9]                 [,10]               
# [1,] "3"                 "3"                  "3"                  "3"                  "3"                 
# [2,] "27"                "14"                 "13"                 "20"                 "6"                 
# [3,] "43"                "27"                 "27"                 "27"                 "27"                
# [4,] "ARNT"              "ARNT"               "ARNT"               "ARNT"               "ARNT"              
# [5,] "NR1I3"             "HNF4G"              "HNF4A"              "NFE2L2"             "CEBPD"             
# [6,] "YY1"               "NR1I3"              "NR1I3"              "NR1I3"              "NR1I3"             
# [7,] "229.537661932704"  "217.952972927102"   "229.537081063701"   "222.639481110433"   "229.53184402284"   
# [8,] "205.85304308581"   "206.525086897473"   "206.668113234189"   "206.754590736461"   "206.820033315982"  
# [9,] "0.103184020641624" "0.0524328063809059" "0.0996308209703427" "0.0713480389675064" "0.0989484086774413"
# [,11]                [,12]                [,13]                [,14]                [,15]             
# [1,] "14"                 "14"                 "3"                  "1"                  "3"               
# [2,] "27"                 "27"                 "27"                 "14"                 "27"              
# [3,] "36"                 "33"                 "34"                 "27"                 "32"              
# [4,] "HNF4G"              "HNF4G"              "ARNT"               "AHR"                "ARNT"            
# [5,] "NR1I3"              "NR1I3"              "NR1I3"              "HNF4G"              "NR1I3"           
# [6,] "PPARG"              "PGRMC1"             "PPARA"              "NR1I3"              "ONECUT1"         
# [7,] "217.766384385754"   "221.540873139369"   "224.236600638626"   "220.672313469865"   "229.390327943791"
# [8,] "206.9405204198"     "208.203758993092"   "208.342950064166"   "208.4666237553"     "208.946627116529"
# [9,] "0.0497132006691044" "0.0602015960182984" "0.0708789311343272" "0.0553113778645015" "0.08912189546314"

##----------------------------------------------- split peaks ---------------------------------------------
index.low <- which(Residual.Deviance<=242)
index.high <- which(Residual.Deviance>242)

pg.tb <- table(gene.names.m[,index.low])
npg.tb <- table(gene.names.m[,index.high])

pg.tb[order(pg.tb)]
npg.tb[order(npg.tb)]

pg.tb[order(pg.tb)]
# AHR   CEBPA   CEBPB   CEBPD   FOXA1   FOXA3   NCOA1   NCOA3   NCOR1   NCOR2   NR1H4   NR2F1   NR2F2   NR5A2 
# 41      41      41      41      41      41      41      41      41      41      41      41      41      41 
# PGRMC1   PPARA   PPARG    RXRA    RXRB    THRA    THRB    USF1     VDR   FOXA2  NFE2L2   NR0B2   NR1D2 ONECUT1 
# 41      41      41      41      41      41      41      41      41      42      42      42      42      42 
# AHRR   HNF4A   NCOA2   NR1H2   NR3C1   PPARD   CEBPG   NR1I2     DBP     YY1   NR1H3    ARNT   HNF4G    ESR1 
# 43      43      43      43      43      43      44      44      45      45      47      50      50      64 
# NR1I3 
# 861 

npg.tb[order(npg.tb)]
# ESR1    ARNT   HNF4G   NR1H3     DBP     YY1   CEBPG   NR1I2    AHRR   HNF4A   NCOA2   NR1H2   NR3C1   PPARD 
# 797     811     811     814     816     816     817     817     818     818     818     818     818     818 
# FOXA2  NFE2L2   NR0B2   NR1D2 ONECUT1     AHR   CEBPA   CEBPB   CEBPD   FOXA1   FOXA3   NCOA1   NCOA3   NCOR1 
# 819     819     819     819     819     820     820     820     820     820     820     820     820     820 
# NCOR2   NR1H4   NR2F1   NR2F2   NR5A2  PGRMC1   PPARA   PPARG    RXRA    RXRB    THRA    THRB    USF1     VDR 
# 820     820     820     820     820     820     820     820     820     820     820     820     820     820 

##-------------------------- gene 17 : CYP1A2 -----------------------------------------------
which(y.m[,17]==0)

dim(gene.tb)
dim(log.gene.tb)
dim(gene.tb.tmp)
dim(log.gene.tb.tmp)

gene.tb.tmp <- gene.tb
log.gene.tb.tmp <- log.gene.tb

gene.tb <- gene.tb[-which(y.m[,17]==0),]
log.gene.tb <- log.gene.tb[-which(y.m[,17]==0),]

log.CYP1A2 <- log(y.m[,17])
log.CYP1A2 <- log.CYP1A2[-which(y.m[,17]==0)]
hist(log.CYP1A2, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP1A2, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP1A2 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 208.7   330.6   374.2   373.4   416.3   498.2 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 8.409e-05 2.429e-02 4.091e-02 4.855e-02 6.418e-02 2.589e-01 

goodfit.index <- which((Residual.Deviance<=210) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]               [,3]                [,4]               [,5]               
# [1,] "9"                 "9"                "9"                 "9"                "9"                
# [2,] "27"                "14"               "27"                "17"               "27"               
# [3,] "34"                "27"               "38"                "27"               "40"               
# [4,] "ESR1"              "ESR1"             "ESR1"              "ESR1"             "ESR1"             
# [5,] "NR1I3"             "HNF4G"            "NR1I3"             "NCOA3"            "NR1I3"            
# [6,] "PPARA"             "NR1I3"            "RXRB"              "NR1I3"            "THRB"             
# [7,] "281.634556077437"  "265.482524707936" "253.19261887052"   "251.468729567784" "293.34586170725"  
# [8,] "208.718650455503"  "214.660236487717" "215.048887129739"  "219.153855291455" "220.101613260252" 
# [9,] "0.258902553143676" "0.19143364813232" "0.150651041530906" "0.12850454341528" "0.249685637358993"
# [,6]                [,7]                [,8]                [,9]                [,10]               
# [1,] "9"                 "14"                "17"                "9"                 "17"                
# [2,] "27"                "27"                "27"                "27"                "27"                
# [3,] "31"                "33"                "41"                "41"                "33"                
# [4,] "ESR1"              "HNF4G"             "NCOA3"             "ESR1"              "NCOA3"             
# [5,] "NR1I3"             "NR1I3"             "NR1I3"             "NR1I3"             "NR1I3"             
# [6,] "NR5A2"             "PGRMC1"            "USF1"              "USF1"              "PGRMC1"            
# [7,] "288.215700543178"  "242.124618645714"  "260.922298590099"  "279.180861649138"  "237.986396293697"  
# [8,] "221.081445790729"  "222.313516481988"  "223.878926000947"  "225.202262835053"  "225.229109053791"  
# [9,] "0.232930595473898" "0.081821924075858" "0.141970896275698" "0.193346343639857" "0.0536051112104826"
# [,11]               [,12]               [,13]               [,14]              [,15]               
# [1,] "6"                 "1"                 "9"                 "5"                "7"                 
# [2,] "9"                 "9"                 "27"                "9"                "9"                 
# [3,] "27"                "27"                "29"                "27"               "29"                
# [4,] "CEBPD"             "AHR"               "ESR1"              "CEBPB"            "CEBPG"             
# [5,] "ESR1"              "ESR1"              "NR1I3"             "ESR1"             "ESR1"              
# [6,] "NR1I3"             "NR1I3"             "NR2F2"             "NR1I3"            "NR2F2"             
# [7,] "293.598267047767"  "295.226226272714"  "258.274952201286"  "294.901353700752" "253.773639145574"  
# [8,] "225.27129487657"   "226.128657029051"  "226.146605172692"  "228.797230877773" "228.855407999657"  
# [9,] "0.232722668489323" "0.234049562994565" "0.124395907364468" "0.22415672899914" "0.0981907783243899"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP1A2, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP1A2 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 191.2   327.3   370.1   367.7   420.5   504.2 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0004658 0.0285100 0.0481100 0.0539300 0.0737000 0.2378000 

goodfit.index <- which((Residual.Deviance<=205) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                 [,2]                 [,3]                 [,4]                 [,5]               
# [1,] "9"                  "9"                  "9"                  "8"                  "9"                
# [2,] "27"                 "14"                 "24"                 "9"                  "27"               
# [3,] "29"                 "27"                 "29"                 "29"                 "40"               
# [4,] "ESR1"               "ESR1"               "ESR1"               "DBP"                "ESR1"             
# [5,] "NR1I3"              "HNF4G"              "NR1H3"              "ESR1"               "NR1I3"            
# [6,] "NR2F2"              "NR1I3"              "NR2F2"              "NR2F2"              "THRB"             
# [7,] "195.27625737435"    "207.978002198323"   "219.128123375882"   "216.485641865648"   "231.981372175442" 
# [8,] "191.229831008407"   "199.913311183601"   "202.061720617941"   "203.505086471196"   "204.013678008119" 
# [9,] "0.0207215481305837" "0.0387766539224255" "0.0778832150570927" "0.0599603524861361" "0.120560085945915"
# [,6]                 [,7]                 [,8]                [,9]                 [,10]               
# [1,] "9"                  "7"                  "27"                "9"                  "9"                 
# [2,] "27"                 "9"                  "34"                "14"                 "29"                
# [3,] "34"                 "29"                 "41"                "28"                 "35"                
# [4,] "ESR1"               "CEBPG"              "NR1I3"             "ESR1"               "ESR1"              
# [5,] "NR1I3"              "ESR1"               "PPARA"             "HNF4G"              "NR2F2"             
# [6,] "PPARA"              "NR2F2"              "USF1"              "NR2F1"              "PPARD"             
# [7,] "214.27576558769"    "209.648710312305"   "229.971378262222"  "210.545964373537"   "213.570380125012"  
# [8,] "204.433923822552"   "205.869574852387"   "206.025557799206"  "207.88650392054"    "208.37201327326"   
# [9,] "0.0459307273416814" "0.0180260372424345" "0.104125220468579" "0.0126312582666213" "0.0243402987282675"
# [,11]                [,12]               [,13]                [,14]                [,15]               
# [1,] "9"                  "1"                 "9"                  "9"                  "9"                 
# [2,] "27"                 "9"                 "16"                 "23"                 "26"                
# [3,] "38"                 "27"                "27"                 "29"                 "29"                
# [4,] "ESR1"               "AHR"               "ESR1"               "ESR1"               "ESR1"              
# [5,] "NR1I3"              "ESR1"              "NCOA2"              "NR1H2"              "NR1I2"             
# [6,] "RXRB"               "NR1I3"             "NR1I3"              "NR2F2"              "NR2F2"             
# [7,] "215.151315301037"   "236.581987193084"  "230.761148627192"   "214.12153313282"    "217.548362675525"  
# [8,] "209.144690952679"   "209.489356796037"  "209.53827872273"    "209.962153822025"   "210.386593138907"  
# [9,] "0.0279181391010999" "0.114516877292671" "0.0919689905805974" "0.0194253200504358" "0.0329203559545959"

##----------------------------------------------- split peaks ---------------------------------------------
index.low <- which(Residual.Deviance<=260)
index.high <- which(Residual.Deviance>260)

pg.tb <- table(gene.names.m[,index.low])
npg.tb <- table(gene.names.m[,index.high])

pg.tb[order(pg.tb)]
npg.tb[order(npg.tb)]

pg.tb[order(pg.tb)]
# NR1H4 ONECUT1   FOXA1   NCOR1   NR1I2    THRA     VDR    AHRR   CEBPA   CEBPG     DBP   FOXA3   NCOA2  NFE2L2 
# 44      44      45      45      45      45      45      46      46      46      46      46      46      46 
# NR0B2   NR1D2     YY1   HNF4A   NCOA1   PPARG    ARNT   FOXA2    RXRA   CEBPB   PPARD    THRB   CEBPD   NR2F1 
# 46      46      46      47      47      47      48      48      48      49      49      49      50      51 
# NR1H3     AHR   HNF4G   NCOR2   NR3C1   NR5A2   NCOA3  PGRMC1   NR1H2   PPARA    RXRB    USF1   NR2F2   NR1I3 
# 52      53      53      54      63      63      64      66      73      73      80      80      84     261 
# ESR1 
# 860 

npg.tb[order(npg.tb)]
# ESR1   NR1I3   NR2F2    RXRB    USF1   NR1H2   PPARA  PGRMC1   NCOA3   NR3C1   NR5A2   NCOR2     AHR   HNF4G 
# 1     600     777     781     781     788     788     795     797     798     798     807     808     808 
# NR1H3   NR2F1   CEBPD   CEBPB   PPARD    THRB    ARNT   FOXA2    RXRA   HNF4A   NCOA1   PPARG    AHRR   CEBPA 
# 809     810     811     812     812     812     813     813     813     814     814     814     815     815 
# CEBPG     DBP   FOXA3   NCOA2  NFE2L2   NR0B2   NR1D2     YY1   FOXA1   NCOR1   NR1I2    THRA     VDR   NR1H4 
# 815     815     815     815     815     815     815     815     816     816     816     816     816     817 
# ONECUT1 
# 817 
gene.tb <- gene.tb.tmp
log.gene.tb <- log.gene.tb.tmp

##-------------------------- gene 18 : CYP4F12 -----------------------------------------------
which(y.m[,18]==0)

dim(gene.tb)
dim(log.gene.tb)
dim(gene.tb.tmp)
dim(log.gene.tb.tmp)

gene.tb <- gene.tb.tmp
log.gene.tb <- log.gene.tb.tmp

#gene.tb <- gene.tb[-which(y.m[,9]==0),]
#log.gene.tb <- log.gene.tb[-which(y.m[,9]==0),]
log.CYP4F12 <- log(y.m[,18])
hist(log.CYP4F12, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP4F12, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP4F12 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 22.12   54.61   71.19   68.66   83.27  107.70 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0005759 0.0258100 0.0453400 0.0592200 0.0770100 0.3824000 

goodfit.index <- which((Residual.Deviance<=28) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]               
# [1,] "8"                 "8"                 "8"                 "8"                 "8"                
# [2,] "22"                "27"                "27"                "9"                 "27"               
# [3,] "27"                "35"                "33"                "27"                "39"               
# [4,] "DBP"               "DBP"               "DBP"               "DBP"               "DBP"              
# [5,] "NR1D2"             "NR1I3"             "NR1I3"             "ESR1"              "NR1I3"            
# [6,] "NR1I3"             "PPARD"             "PGRMC1"            "NR1I3"             "THRA"             
# [7,] "27.5173738054892"  "27.0001032861992"  "27.8187321475984"  "27.4474413023594"  "28.2318908508024" 
# [8,] "22.1234441823789"  "22.3991136341072"  "22.495160889529"   "22.9256271032147"  "23.0554755876976" 
# [9,] "0.196019055497013" "0.170406372276501" "0.191366422805468" "0.164744471054063" "0.183353473929913"
# [,6]                [,7]                [,8]                [,9]                [,10]              
# [1,] "8"                 "5"                 "8"                 "4"                 "8"                
# [2,] "12"                "8"                 "27"                "8"                 "27"               
# [3,] "27"                "27"                "42"                "27"                "32"               
# [4,] "DBP"               "CEBPB"             "DBP"               "CEBPA"             "DBP"              
# [5,] "FOXA3"             "DBP"               "NR1I3"             "DBP"               "NR1I3"            
# [6,] "NR1I3"             "NR1I3"             "VDR"               "NR1I3"             "ONECUT1"          
# [7,] "28.4807396863168"  "27.6649427810754"  "28.4782858951332"  "28.3067476107312"  "27.1124332813734" 
# [8,] "23.156628847021"   "23.2230810155061"  "23.3514726833337"  "23.464306821493"   "23.4643102746228" 
# [9,] "0.186937238917769" "0.160559224745906" "0.180025343894579" "0.171070193433399" "0.134555352110608"
# [,11]               [,12]               [,13]               [,14]               [,15]              
# [1,] "8"                 "8"                 "8"                 "8"                 "8"                
# [2,] "26"                "27"                "24"                "27"                "13"               
# [3,] "27"                "30"                "27"                "36"                "27"               
# [4,] "DBP"               "DBP"               "DBP"               "DBP"               "DBP"              
# [5,] "NR1I2"             "NR1I3"             "NR1H3"             "NR1I3"             "HNF4A"            
# [6,] "NR1I3"             "NR3C1"             "NR1I3"             "PPARG"             "NR1I3"            
# [7,] "28.0439397003598"  "28.3811396901761"  "28.1448309064371"  "28.1406364718764"  "28.4208098348791" 
# [8,] "23.4962680567456"  "23.601039127342"   "23.6299535214775"  "23.693003263636"   "23.9594021544825" 
# [9,] "0.162162367064135" "0.168425250536671" "0.160415864638467" "0.158050199493011" "0.156976796450091"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP4F12, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP4F12 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 21.27   51.02   70.33   66.94   83.86  110.60 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001021 0.0255100 0.0456800 0.0540000 0.0726000 0.3235000

goodfit.index <- which((Residual.Deviance<=29.875) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                 [,4]                 [,5]                
# [1,] "8"                 "5"                 "6"                  "8"                  "8"                 
# [2,] "27"                "8"                 "8"                  "27"                 "20"                
# [3,] "32"                "27"                "27"                 "34"                 "27"                
# [4,] "DBP"               "CEBPB"             "CEBPD"              "DBP"                "DBP"               
# [5,] "NR1I3"             "DBP"               "DBP"                "NR1I3"              "NFE2L2"            
# [6,] "ONECUT1"           "NR1I3"             "NR1I3"              "PPARA"              "NR1I3"             
# [7,] "23.6671526828331"  "23.8508651187586"  "22.9455561840926"   "23.5770463599259"   "23.930845079095"   
# [8,] "21.2729917177448"  "21.3432087363527"  "21.484584491207"    "21.508131584329"    "21.7160599623502"  
# [9,] "0.101159653515266" "0.105139011516764" "0.0636712259735266" "0.0877512282078481" "0.0925493901040509"
# [,6]                 [,7]                 [,8]                 [,9]                 [,10]               
# [1,] "8"                  "8"                  "8"                  "8"                  "8"                 
# [2,] "27"                 "27"                 "23"                 "27"                 "10"                
# [3,] "31"                 "33"                 "27"                 "40"                 "27"                
# [4,] "DBP"                "DBP"                "DBP"                "DBP"                "DBP"               
# [5,] "NR1I3"              "NR1I3"              "NR1H2"              "NR1I3"              "FOXA1"             
# [6,] "NR5A2"              "PGRMC1"             "NR1I3"              "THRB"               "NR1I3"             
# [7,] "23.5753975572001"   "23.9888515922582"   "23.454233920235"    "23.5699914169639"   "23.2026127801039"  
# [8,] "21.8707451117642"   "21.9709934358256"   "22.0204191806165"   "22.0204698875502"   "22.0240089275901"  
# [9,] "0.0723064135525144" "0.0841164967264975" "0.0611324481752283" "0.0657412852640462" "0.0507961695384745"
# [,11]                [,12]                [,13]                [,14]                [,15]               
# [1,] "8"                  "8"                  "8"                  "8"                  "1"                 
# [2,] "27"                 "27"                 "27"                 "12"                 "8"                 
# [3,] "39"                 "41"                 "43"                 "27"                 "27"                
# [4,] "DBP"                "DBP"                "DBP"                "DBP"                "AHR"               
# [5,] "NR1I3"              "NR1I3"              "NR1I3"              "FOXA3"              "DBP"               
# [6,] "THRA"               "USF1"               "YY1"                "NR1I3"              "NR1I3"             
# [7,] "23.8143065974593"   "23.9801800992793"   "23.7013935373648"   "23.7850703730992"   "23.9731150037612"  
# [8,] "22.1301117790994"   "22.1967400923551"   "22.227595392743"    "22.2380123031264"   "22.2401013680517"  
# [9,] "0.0707219759461552" "0.0743714183772033" "0.0621819194849621" "0.0650432412309579" "0.0722898812039076"

##----------------------------------------------- split peaks ---------------------------------------------
index.low <- which(Residual.Deviance<=24)
index.low.1 <- which((Residual.Deviance>24) & ( Residual.Deviance<=34))
index.high.1 <- which((Residual.Deviance>34) & ( Residual.Deviance<=42))
index.high <- which(Residual.Deviance>42)

pg.tb <- table(gene.names.m[,index.low])
pg.tb.1 <- table(gene.names.m[,index.low.1])
npg.tb.1 <- table(gene.names.m[,index.high.1])
npg.tb <- table(gene.names.m[,index.high])

pg.tb[order(pg.tb)]
pg.tb.1[order(pg.tb.1)]
npg.tb.1[order(npg.tb.1)]
npg.tb[order(npg.tb)]

pg.tb[order(pg.tb)]
# AHR    AHRR    ARNT   CEBPA   CEBPB   CEBPD   CEBPG    ESR1   FOXA1   FOXA2   FOXA3   HNF4A   HNF4G   NCOA1 
# 1       1       1       1       1       1       1       1       1       1       1       1       1       1 
# NCOA2   NCOA3   NCOR1   NCOR2  NFE2L2   NR0B2   NR1H2   NR1H3   NR1I2   NR2F1   NR2F2   NR3C1   NR5A2 ONECUT1 
# 1       1       1       1       1       1       1       1       1       1       1       1       1       1 
# PGRMC1   PPARA   PPARD   PPARG    RXRA    RXRB    THRA    THRB    USF1     VDR     YY1   NR1D2   NR1H4     DBP 
# 1       1       1       1       1       1       1       1       1       1       1       2       2      41 
# NR1I3 
# 42 

pg.tb.1[order(pg.tb.1)]
# NR1H4    ARNT   CEBPB   FOXA3   HNF4G   NCOA2   NCOR1   NR5A2   PPARG     VDR    AHRR  NFE2L2   NR0B2   NR2F1 
# 42      43      43      43      43      43      43      43      43      43      44      44      44      44 
# ONECUT1    RXRA    THRB     AHR   CEBPA   CEBPD   CEBPG   FOXA2   NCOR2   NR3C1   FOXA1   NCOA1   NCOA3   NR1D2 
# 44      44      44      45      45      45      45      45      45      45      46      46      46      46 
# HNF4A   NR1H2   NR2F2   PPARA    RXRB    USF1    THRA     YY1   PPARD  PGRMC1   NR1I2    ESR1   NR1H3     DBP 
# 47      47      47      47      47      47      48      54      56      62      84      86      90     165 
# NR1I3 
# 819 

npg.tb.1[order(npg.tb.1)]
# PGRMC1   PPARA    THRA     AHR   FOXA2     YY1   NCOR1 ONECUT1     VDR    AHRR    ARNT   HNF4G  NFE2L2   PPARG 
# 32      40      40      42      42      42      43      43      43      44      44      44      44      44 
# FOXA1   FOXA3   HNF4A   NCOA2    THRB   NCOA3   NR0B2   NR2F1    USF1   NCOA1   CEBPA   CEBPB   NR1H4   CEBPD 
# 45      45      45      45      45      46      46      46      46      47      48      48      48      49 
# NR5A2   NR2F2    RXRB   CEBPG   NR1D2   NCOR2   PPARD    ESR1   NR1H2    RXRA   NR3C1   NR1I2   NR1H3     DBP 
# 49      50      51      60      61      69      69      73      77      77      92     130     245     655 

npg.tb[order(npg.tb)]
# NR1H3   NR1I2    ESR1   NR3C1   PPARD   NR1H2    RXRA   NCOR2   NR1D2   CEBPG    RXRB   NR2F2     YY1   CEBPD 
# 525     646     701     723     735     736     739     746     752     755     762     763     764     766 
# PGRMC1   CEBPA   NCOA1    USF1   HNF4A   NCOA3   NR5A2   CEBPB   FOXA1   NR1H4   NR0B2   NR2F1    THRB    AHRR 
# 766     767     767     767     768     768     768     769     769     769     770     770     771     772 
# FOXA3   NCOA2  NFE2L2    THRA     AHR    ARNT   FOXA2   HNF4G ONECUT1   PPARA   PPARG   NCOR1     VDR 
# 772     772     772     772     773     773     773     773     773     773     773     774     774 

##-------------------------- gene 19 : CYP4F11 -----------------------------------------------
which(y.m[,19]==0)

dim(gene.tb)
dim(log.gene.tb)
dim(gene.tb.tmp)
dim(log.gene.tb.tmp)

gene.tb <- gene.tb.tmp
log.gene.tb <- log.gene.tb.tmp

#gene.tb <- gene.tb[-which(y.m[,9]==0),]
#log.gene.tb <- log.gene.tb[-which(y.m[,9]==0),]
log.CYP4F11 <- log(y.m[,19])
hist(log.CYP4F11, breaks=20) 


index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP4F11, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP4F11 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.623  10.300  11.210  11.090  11.940  13.800 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0004723 0.0303600 0.0558000 0.0694100 0.0936200 0.3218000

large.interaction.index <- which((Residual.Deviance<=10) & (Residual.Deviance.Deduction>=0.30) )

goodfit.index <- which((Residual.Deviance<=8) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                 [,3]                 [,4]                 [,5]                
# [1,] "10"                "1"                  "1"                  "13"                 "13"                
# [2,] "13"                "13"                 "13"                 "28"                 "28"                
# [3,] "28"                "29"                 "28"                 "43"                 "36"                
# [4,] "FOXA1"             "AHR"                "AHR"                "HNF4A"              "HNF4A"             
# [5,] "HNF4A"             "HNF4A"              "HNF4A"              "NR2F1"              "NR2F1"             
# [6,] "NR2F1"             "NR2F2"              "NR2F1"              "YY1"                "PPARG"             
# [7,] "7.44915123490226"  "7.04823624752557"   "7.26500279820809"   "7.82580892870382"   "7.90898149438695"  
# [8,] "6.62264942945582"  "6.7633351843612"    "6.85135118039803"   "7.15497521543672"   "7.16381172003627"  
# [9,] "0.110952480273718" "0.0404216109050524" "0.0569375717118907" "0.0857206864336535" "0.0942181714395874"
# [,6]                 [,7]                [,8]                [,9]                [,10]              
# [1,] "10"                 "1"                 "3"                 "6"                 "1"                
# [2,] "13"                 "28"                "13"                "13"                "13"               
# [3,] "29"                 "37"                "43"                "15"                "17"               
# [4,] "FOXA1"              "AHR"               "ARNT"              "CEBPD"             "AHR"              
# [5,] "HNF4A"              "NR2F1"             "HNF4A"             "HNF4A"             "HNF4A"            
# [6,] "NR2F2"              "RXRA"              "YY1"               "NCOA1"             "NCOA3"            
# [7,] "7.31919054078508"   "8.09353625051143"  "9.73531207730462"  "9.65733046812982"  "9.21483259187958" 
# [8,] "7.17963175030524"   "7.19264532177862"  "7.23337324988172"  "7.31052016322634"  "7.31911103029098" 
# [9,] "0.0190675170569984" "0.111309926940264" "0.256996263453693" "0.243008180433319" "0.205725013741343"
# [,11]               [,12]                [,13]                [,14]               [,15]               
# [1,] "6"                 "10"                 "1"                  "13"                "1"                 
# [2,] "10"                "13"                 "13"                 "36"                "10"                
# [3,] "13"                "39"                 "39"                 "39"                "28"                
# [4,] "CEBPD"             "FOXA1"              "AHR"                "HNF4A"             "AHR"               
# [5,] "FOXA1"             "HNF4A"              "HNF4A"              "PPARG"             "FOXA1"             
# [6,] "HNF4A"             "THRA"               "THRA"               "THRA"              "NR2F1"             
# [7,] "8.42558322939061"  "7.96204693526417"   "7.63014987525941"   "8.38110951529726"  "8.15618916057845"  
# [8,] "7.35300993965621"  "7.36552502362121"   "7.37273554877933"   "7.40291999668429"  "7.41354610435245"  
# [9,] "0.127299589895806" "0.0749206725975133" "0.0337364705396863" "0.116713606572921" "0.0910527014031795"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP4F11, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP4F11 Residual Deviance")

summary(Residual.Deviance)
summary(Residual.Deviance.Deduction)

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.545  10.120  10.910  10.820  11.610  13.770 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0005497 0.0431600 0.0784900 0.0901000 0.1261000 0.3359000 

goodfit.index <- which((Residual.Deviance<=6.715) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]               
# [1,] "10"                "13"                "13"                "13"                "6"                
# [2,] "13"                "28"                "25"                "28"                "13"               
# [3,] "28"                "43"                "28"                "36"                "15"               
# [4,] "FOXA1"             "HNF4A"             "HNF4A"             "HNF4A"             "CEBPD"            
# [5,] "HNF4A"             "NR2F1"             "NR1H4"             "NR2F1"             "HNF4A"            
# [6,] "NR2F1"             "YY1"               "NR2F1"             "PPARG"             "NCOA1"            
# [7,] "7.71098267308445"  "7.81024168029515"  "8.54598965508894"  "8.08313803259145"  "9.56340209069743" 
# [8,] "6.54516896360122"  "6.61432629464921"  "6.71021324092201"  "6.71557331303644"  "6.83667899428607" 
# [9,] "0.151188734161284" "0.153121431397337" "0.214811448206443" "0.169187351996335" "0.285120616131337"
# [,6]                [,7]                [,8]                 [,9]                [,10]              
# [1,] "5"                 "1"                 "1"                  "6"                 "13"               
# [2,] "10"                "13"                "13"                 "10"                "15"               
# [3,] "13"                "23"                "29"                 "13"                "28"               
# [4,] "CEBPB"             "AHR"               "AHR"                "CEBPD"             "HNF4A"            
# [5,] "FOXA1"             "HNF4A"             "HNF4A"              "FOXA1"             "NCOA1"            
# [6,] "HNF4A"             "NR1H2"             "NR2F2"              "HNF4A"             "NR2F1"            
# [7,] "8.689829492572"    "9.22295742875153"  "7.16162654531275"   "8.74425202852187"  "8.23720518492227" 
# [8,] "6.85369689873303"  "6.87518570360854"  "6.96665392508502"   "6.99251301639784"  "7.01424255326787" 
# [9,] "0.211296734350022" "0.254557363327306" "0.0272246282313259" "0.200330343454218" "0.148468152024787"
# [,11]               [,12]               [,13]               [,14]               [,15]             
# [1,] "5"                 "13"                "13"                "1"                 "10"              
# [2,] "13"                "17"                "15"                "4"                 "13"              
# [3,] "28"                "28"                "35"                "28"                "37"              
# [4,] "CEBPB"             "HNF4A"             "HNF4A"             "AHR"               "FOXA1"           
# [5,] "HNF4A"             "NCOA3"             "NCOA1"             "CEBPA"             "HNF4A"           
# [6,] "NR2F1"             "NR2F1"             "PPARD"             "NR2F1"             "RXRA"            
# [7,] "8.38221980624216"  "8.02165637879419"  "9.6951510701561"   "8.46524677423484"  "8.76622554592706"
# [8,] "7.03682357895766"  "7.05065711800319"  "7.12731739206495"  "7.14969619027574"  "7.17977996233459"
# [9,] "0.160505958849062" "0.121047227023823" "0.264857520992688" "0.155406052421668" "0.18097248071943"

##-------------------------- gene 20 : CYP2A6 -----------------------------------------------
which(y.m[,20]==0)

dim(gene.tb)
dim(log.gene.tb)
dim(gene.tb.tmp)
dim(log.gene.tb.tmp)

gene.tb <- gene.tb.tmp
log.gene.tb <- log.gene.tb.tmp

#gene.tb <- gene.tb[-which(y.m[,9]==0),]
#log.gene.tb <- log.gene.tb[-which(y.m[,9]==0),]
log.CYP2A6 <- log(y.m[,20])
hist(log.CYP2A6, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP2A6, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2A6 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 153.0   255.9   295.9   292.4   328.8   405.9 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001628 0.0233500 0.0416000 0.0504900 0.0688000 0.2413000

goodfit.index <- which((Residual.Deviance<=159) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]                [,2]                [,3]                [,4]                [,5]                
# [1,] "13"                "13"                "9"                 "10"                "10"                
# [2,] "27"                "26"                "10"                "13"                "27"                
# [3,] "37"                "37"                "27"                "27"                "33"                
# [4,] "HNF4A"             "HNF4A"             "ESR1"              "FOXA1"             "FOXA1"             
# [5,] "NR1I3"             "NR1I2"             "FOXA1"             "HNF4A"             "NR1I3"             
# [6,] "RXRA"              "RXRA"              "NR1I3"             "NR1I3"             "PGRMC1"            
# [7,] "191.642845244155"  "187.069151237446"  "188.423837358925"  "192.849724874572"  "180.366309316741"  
# [8,] "153.030629622701"  "153.468863513829"  "154.57866343872"   "160.898408474686"  "163.567796174748"  
# [9,] "0.201480079114158" "0.179614262968288" "0.179622570024057" "0.165679865090119" "0.0931355373718557"
# [,6]                [,7]                [,8]                [,9]                [,10]              
# [1,] "10"                "20"                "10"                "21"                "25"               
# [2,] "27"                "26"                "27"                "27"                "27"               
# [3,] "28"                "30"                "41"                "37"                "37"               
# [4,] "FOXA1"             "NFE2L2"            "FOXA1"             "NR0B2"             "NR1H4"            
# [5,] "NR1I3"             "NR1I2"             "NR1I3"             "NR1I3"             "NR1I3"            
# [6,] "NR2F1"             "NR3C1"             "USF1"              "RXRA"              "RXRA"             
# [7,] "194.607815096124"  "220.039894920051"  "195.47869265692"   "204.483038967817"  "203.474358597952" 
# [8,] "166.805133790005"  "166.948291009899"  "167.843587739864"  "167.862924777829"  "167.93302153059"  
# [9,] "0.142865183972112" "0.241281718160437" "0.141371443308956" "0.179086316277564" "0.174672314055987"
# [,11]               [,12]               [,13]               [,14]               [,15]              
# [1,] "10"                "10"                "11"                "26"                "21"               
# [2,] "26"                "24"                "27"                "27"                "26"               
# [3,] "27"                "27"                "37"                "37"                "30"               
# [4,] "FOXA1"             "FOXA1"             "FOXA2"             "NR1I2"             "NR0B2"            
# [5,] "NR1I2"             "NR1H3"             "NR1I3"             "NR1I3"             "NR1I2"            
# [6,] "NR1I3"             "NR1I3"             "RXRA"              "RXRA"              "NR3C1"            
# [7,] "197.711411004375"  "205.315407341222"  "196.409387505313"  "200.617824948476"  "217.341581125626" 
# [8,] "168.276111700731"  "168.581822372791"  "168.772175602278"  "169.373974041332"  "169.667204889062" 
# [9,] "0.148880123580692" "0.178912948833801" "0.140712275793271" "0.155738159932541" "0.219352302443258"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP2A6, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2A6 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 136.1   243.7   293.3   284.1   330.7   401.4 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003957 0.0239000 0.0404100 0.0471100 0.0636700 0.2411000

goodfit.index <- which((Residual.Deviance<=145.65) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]               [,2]               [,3]                 [,4]                 [,5]               
# [1,] "21"               "21"               "26"                 "21"                 "1"                
# [2,] "26"               "27"               "30"                 "26"                 "26"               
# [3,] "37"               "37"               "40"                 "30"                 "30"               
# [4,] "NR0B2"            "NR0B2"            "NR1I2"              "NR0B2"              "AHR"              
# [5,] "NR1I2"            "NR1I3"            "NR3C1"              "NR1I2"              "NR1I2"            
# [6,] "RXRA"             "RXRA"             "THRB"               "NR3C1"              "NR3C1"            
# [7,] "155.827518496768" "162.965893835662" "155.202582512263"   "153.779085460486"   "162.933664200381" 
# [8,] "136.084979253601" "136.860863811745" "140.669170685845"   "140.784612555857"   "140.968139033329" 
# [9,] "0.12669481894866" "0.1601870760163"  "0.0936415592522063" "0.0845009115883152" "0.134812687573506"
# [,6]                [,7]                [,8]                [,9]               [,10]              
# [1,] "26"                "1"                 "11"                "21"               "10"               
# [2,] "30"                "26"                "26"                "27"               "27"               
# [3,] "38"                "37"                "27"                "30"               "28"               
# [4,] "NR1I2"             "AHR"               "FOXA2"             "NR0B2"            "FOXA1"            
# [5,] "NR3C1"             "NR1I2"             "NR1I2"             "NR1I3"            "NR1I3"            
# [6,] "RXRB"              "RXRA"              "NR1I3"             "NR3C1"            "NR2F1"            
# [7,] "162.79459964323"   "170.646653894847"  "177.599716202708"  "165.005564967764" "166.856750430484" 
# [8,] "141.628216367581"  "141.734801114946"  "141.797605534823"  "144.844152694017" "145.316231538494" 
# [9,] "0.130018952238193" "0.169425254583171" "0.201588783098175" "0.12218625643133" "0.129095879168305"
# [,11]               [,12]               [,13]                [,14]               [,15]               
# [1,] "11"                "22"                "26"                 "26"                "15"                
# [2,] "26"                "25"                "28"                 "37"                "26"                
# [3,] "30"                "27"                "37"                 "40"                "30"                
# [4,] "FOXA2"             "NR1D2"             "NR1I2"              "NR1I2"             "NCOA1"             
# [5,] "NR1I2"             "NR1H4"             "NR2F1"              "RXRA"              "NR1I2"             
# [6,] "NR3C1"             "NR1I3"             "RXRA"               "THRB"              "NR3C1"             
# [7,] "163.673293856503"  "170.760628336725"  "148.23303827085"    "169.414388191207"  "152.008192005212"  
# [8,] "145.453079734965"  "145.455670264837"  "145.50956063914"    "145.527817970288"  "145.537109359924"  
# [9,] "0.111320629604438" "0.148189651902599" "0.0183729461628746" "0.140994932460872" "0.0425706178063439"

##----------------------------------------------- split peaks ---------------------------------------------
index.low <- which(Residual.Deviance<=207)
index.high <- which(Residual.Deviance>207)

pg.tb <- table(gene.names.m[,index.low])
npg.tb <- table(gene.names.m[,index.high])

pg.tb[order(pg.tb)]
npg.tb[order(npg.tb)]

pg.tb[order(pg.tb)]
# ONECUT1    THRA    ARNT   HNF4G   NCOR1  NFE2L2    RXRB    AHRR   CEBPB   CEBPD   NR2F2     VDR     AHR   CEBPA 
# 82      82      83      83      83      83      83      84      84      84      84      84      85      85 
# NCOA2   PPARG    THRB    USF1   FOXA3   NCOA3   PPARD   NR1D2   NR1H2   PPARA     YY1   NCOR2   NR1H3   NR2F1 
# 85      85      85      85      86      86      86      87      87      87      87      88      88      88 
# NR1H4   CEBPG   FOXA2   NR5A2  PGRMC1     DBP   NCOA1   FOXA1    RXRA   HNF4A    ESR1   NR0B2   NR3C1   NR1I2 
# 90      91      92      92     101     103     104     107     110     112     134     169     186     861 
# NR1I3 
# 861 

npg.tb[order(npg.tb)]
# NR3C1   NR0B2    ESR1   HNF4A    RXRA   FOXA1   NCOA1     DBP  PGRMC1   FOXA2   NR5A2   CEBPG   NR1H4   NCOR2 
# 675     692     727     749     751     754     757     758     760     769     769     770     771     773 
# NR1H3   NR2F1   NR1D2   NR1H2   PPARA     YY1   FOXA3   NCOA3   PPARD     AHR   CEBPA   NCOA2   PPARG    THRB 
# 773     773     774     774     774     774     775     775     775     776     776     776     776     776 
# USF1    AHRR   CEBPB   CEBPD   NR2F2     VDR    ARNT   HNF4G   NCOR1  NFE2L2    RXRB ONECUT1    THRA 
# 776     777     777     777     777     777     778     778     778     778     778     779     779

##-------------------------- gene 21 :  CYP2B6             
which(y.m[,21]==0)

dim(gene.tb)
dim(log.gene.tb)
dim(gene.tb.tmp)
dim(log.gene.tb.tmp)

gene.tb <- gene.tb.tmp
log.gene.tb <- log.gene.tb.tmp

#gene.tb <- gene.tb[-which(y.m[,9]==0),]
#log.gene.tb <- log.gene.tb[-which(y.m[,9]==0),]
log.CYP2B6 <- log(y.m[,21])
hist(log.CYP2B6, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP2B6, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2B6 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 165.3   248.8   273.8   272.7   298.1   353.4 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002072 0.0219400 0.0392800 0.0478000 0.0653900 0.2378000 

goodfit.index <- which((Residual.Deviance<=180.195) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
# [,1]               [,2]                [,3]                [,4]                [,5]               
# [1,] "27"               "27"                "17"                "17"                "10"               
# [2,] "28"               "37"                "27"                "37"                "37"               
# [3,] "37"               "42"                "37"                "42"                "42"               
# [4,] "NR1I3"            "NR1I3"             "NCOA3"             "NCOA3"             "FOXA1"            
# [5,] "NR2F1"            "RXRA"              "NR1I3"             "RXRA"              "RXRA"             
# [6,] "RXRA"             "VDR"               "RXRA"              "VDR"               "VDR"              
# [7,] "200.564373399497" "208.88767687379"   "205.537774879034"  "208.652200688204"  "217.520764888093" 
# [8,] "165.322933796463" "166.054587874715"  "168.463210265777"  "170.737779848753"  "171.107449204452" 
# [9,] "0.17571136391625" "0.205053211563818" "0.180378349600584" "0.181711099688362" "0.213374183873979"
# [,6]                [,7]                [,8]                [,9]                [,10]              
# [1,] "13"                "27"                "10"                "13"                "9"                
# [2,] "17"                "29"                "27"                "27"                "27"               
# [3,] "37"                "37"                "33"                "37"                "37"               
# [4,] "HNF4A"             "NR1I3"             "FOXA1"             "HNF4A"             "ESR1"             
# [5,] "NCOA3"             "NR2F2"             "NR1I3"             "NR1I3"             "NR1I3"            
# [6,] "RXRA"              "RXRA"              "PGRMC1"            "RXRA"              "RXRA"             
# [7,] "194.082910158942"  "212.797015311734"  "207.052541405357"  "216.546392436664"  "228.414278446578" 
# [8,] "172.093513141319"  "172.313427596772"  "175.650798271503"  "175.960138956537"  "176.574304479653" 
# [9,] "0.113298986498168" "0.190245091810412" "0.151660747174203" "0.187425211860768" "0.226955925520433"
# [,11]               [,12]               [,13]              [,14]               [,15]              
# [1,] "10"                "11"                "27"               "3"                 "7"                
# [2,] "17"                "27"                "37"               "17"                "27"               
# [3,] "27"                "37"                "39"               "27"                "37"               
# [4,] "FOXA1"             "FOXA2"             "NR1I3"            "ARNT"              "CEBPG"            
# [5,] "NCOA3"             "NR1I3"             "RXRA"             "NCOA3"             "NR1I3"            
# [6,] "NR1I3"             "RXRA"              "THRA"             "NR1I3"             "RXRA"             
# [7,] "212.623541614979"  "213.299213749314"  "224.738495430846" "218.225954828351"  "228.431290144205" 
# [8,] "179.210379138448"  "179.696671123869"  "179.847535346962" "180.688617174604"  "181.003894412904" 
# [9,] "0.157147050711046" "0.157537114341816" "0.19974753322889" "0.172011334230488" "0.207622150631645"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP2B6, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2B6 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 144.2   244.9   273.5   268.4   298.4   351.0 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 9.013e-05 2.385e-02 4.117e-02 4.947e-02 6.656e-02 2.358e-01 

goodfit.index <- which((Residual.Deviance<=150) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]

goodfitr[,order(goodfitr[8,])]
#     [,1]                [,2]                [,3]                [,4]                [,5]                
# [1,] "27"                "10"                "11"                "10"                "26"                
# [2,] "28"                "11"                "27"                "27"                "28"                
# [3,] "37"                "27"                "37"                "28"                "37"                
# [4,] "NR1I3"             "FOXA1"             "FOXA2"             "FOXA1"             "NR1I2"             
# [5,] "NR2F1"             "FOXA2"             "NR1I3"             "NR1I3"             "NR2F1"             
# [6,] "RXRA"              "NR1I3"             "RXRA"              "NR2F1"             "RXRA"              
# [7,] "167.174535352162"  "188.645380620184"  "184.629188565003"  "174.88141479091"   "168.090798980595"  
# [8,] "144.230811671928"  "154.130519129068"  "154.80554489741"   "155.231641753563"  "155.285122998822"  
# [9,] "0.137244130105709" "0.182961604348041" "0.161532658510782" "0.112360556213708" "0.0761830871138393"
#     [,6]                [,7]                [,8]               [,9]                [,10]             
# [1,] "21"                "27"                "21"               "11"                "11"              
# [2,] "27"                "37"                "27"               "21"                "26"              
# [3,] "28"                "42"                "37"               "27"                "27"              
# [4,] "NR0B2"             "NR1I3"             "NR0B2"            "FOXA2"             "FOXA2"           
# [5,] "NR1I3"             "RXRA"              "NR1I3"            "NR0B2"             "NR1I2"           
# [6,] "NR2F1"             "VDR"               "RXRA"             "NR1I3"             "NR1I3"           
# [7,] "183.263919230378"  "193.986024248373"  "190.43487109319"  "193.613989373708"  "193.665536909358"
# [8,] "155.979582692263"  "157.43680730985"   "158.207076895816" "158.691349315408"  "159.396769410897"
# [9,] "0.148880023152927" "0.188411598619737" "0.16923263062259" "0.180372503925287" "0.1769481965937" 
#     [,11]               [,12]              [,13]               [,14]               [,15]               
# [1,] "27"                "11"               "3"                 "7"                 "22"                
# [2,] "29"                "17"               "17"                "11"                "27"                
# [3,] "37"                "27"               "27"                "27"                "28"                
# [4,] "NR1I3"             "FOXA2"            "ARNT"              "CEBPG"             "NR1D2"             
# [5,] "NR2F2"             "NCOA3"            "NCOA3"             "FOXA2"             "NR1I3"             
# [6,] "RXRA"              "NR1I3"            "NR1I3"             "NR1I3"             "NR2F1"             
# [7,] "179.812238324527"  "192.137416127799" "195.904005431616"  "197.208417017667"  "179.350833829611"  
# [8,] "159.518125439144"  "160.832676929645" "161.079842492461"  "161.131732705005"  "161.563989162579"  
# [9,] "0.112862801077847" "0.16292890697214" "0.177761362573625" "0.182936838387732" "0.0991734707178998"

##----------------------------------------------- split peaks ---------------------------------------------
index.low <- which(Residual.Deviance<=215)
index.high <- which(Residual.Deviance>215)

pg.tb <- table(gene.names.m[,index.low])
npg.tb <- table(gene.names.m[,index.high])

pg.tb[order(pg.tb)]
# HNF4G    ARNT   CEBPD   CEBPB   NCOA1   NR1H4 ONECUT1     YY1   CEBPA     VDR   NCOR1   NR1H3   PPARG   NCOR2 
# 66      67      68      69      69      69      69      69      70      70      71      75      75      77 
# USF1   HNF4A   NR1H2    THRA    AHRR   FOXA3   NCOA2    RXRB     AHR   NR1D2   NR2F2   FOXA1   PPARA    THRB 
# 78      79      80      80      81      82      82      82      83      83      83      85      85      85 
# NR5A2     DBP   FOXA2  NFE2L2   CEBPG   PPARD   NCOA3    RXRA   NR2F1  PGRMC1   NR3C1    ESR1   NR0B2   NR1I2 
# 86      87      87      88      92      92      96      96     103     109     116     126     140     753 
# NR1I3 
# 861 

npg.tb[order(npg.tb)]
# NR1I2   NR0B2    ESR1   NR3C1  PGRMC1   NR2F1   NCOA3    RXRA   CEBPG   PPARD  NFE2L2     DBP   FOXA2   NR5A2 
# 108     721     735     745     752     758     765     765     769     769     773     774     774     775 
# FOXA1   PPARA    THRB     AHR   NR1D2   NR2F2   FOXA3   NCOA2    RXRB    AHRR   NR1H2    THRA   HNF4A    USF1 
# 776     776     776     778     778     778     779     779     779     780     781     781     782     783 
# NCOR2   NR1H3   PPARG   NCOR1   CEBPA     VDR   CEBPB   NCOA1   NR1H4 ONECUT1     YY1   CEBPD    ARNT   HNF4G 
# 784     786     786     790     791     791     792     792     792     792     792     793     794     795 

# For CYP2B6:
#   clear predictor: NR1I2   NR1I3
# potential predictor:  none

##-------------------------- gene 22 : CYP2D6 -----------------------------------------------
which(y.m[,22]==0)
dim(gene.tb)
dim(log.gene.tb)
dim(gene.tb.tmp)
dim(log.gene.tb.tmp)

gene.tb <- gene.tb.tmp
log.gene.tb <- log.gene.tb.tmp

#gene.tb <- gene.tb[-which(y.m[,9]==0),]
#log.gene.tb <- log.gene.tb[-which(y.m[,9]==0),]
log.CYP2D6 <- log(y.m[,22])
hist(log.CYP2D6, breaks=20) 

index.m <- combn(c(1:43),3)
dim(index.m)
gene.names <- names(gene.tb)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})

dim(gene.tb)

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP2D6, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2D6 Residual Deviance")

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 27.44   46.94   54.91   54.47   61.81   77.95 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0006494 0.0326300 0.0558000 0.0665600 0.0895100 0.3257000

goodfit.index <- which((Residual.Deviance<=30.7095) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]
#     [,1]                [,2]                [,3]               [,4]                [,5]              
# [1,] "5"                 "27"                "24"               "24"                "23"              
# [2,] "23"                "35"                "27"               "35"                "27"              
# [3,] "27"                "38"                "35"               "43"                "35"              
# [4,] "CEBPB"             "NR1I3"             "NR1H3"            "NR1H3"             "NR1H2"           
# [5,] "NR1H2"             "PPARD"             "NR1I3"            "PPARD"             "NR1I3"           
# [6,] "NR1I3"             "RXRB"              "PPARD"            "YY1"               "PPARD"           
# [7,] "40.6883850942331"  "36.4721871170084"  "36.068548137857"  "34.227374996379"   "36.4340503161018"
# [8,] "27.4380697228873"  "28.2867918875582"  "28.3721279650266" "28.7885982477765"  "28.898518508496" 
# [9,] "0.325653508750925" "0.224428417281097" "0.21338314321425" "0.158901369128596" "0.20682662899753"
#     [,6]                [,7]                [,8]                [,9]                [,10]              
# [1,] "23"                "23"                "11"                "22"                "23"               
# [2,] "24"                "27"                "23"                "23"                "27"               
# [3,] "35"                "33"                "27"                "27"                "38"               
# [4,] "NR1H2"             "NR1H2"             "FOXA2"             "NR1D2"             "NR1H2"            
# [5,] "NR1H3"             "NR1I3"             "NR1H2"             "NR1H2"             "NR1I3"            
# [6,] "PPARD"             "PGRMC1"            "NR1I3"             "NR1I3"             "RXRB"             
# [7,] "35.8321848389473"  "41.136658664565"   "42.5971674643623"  "42.6221402031968"  "42.6499851706981" 
# [8,] "28.9558453326234"  "29.5086720327769"  "29.9413956519363"  "29.970201909823"   "30.0695762266466" 
# [9,] "0.191903997404305" "0.282667261009325" "0.297103600210368" "0.296839582270082" "0.294968659278566"
#     [,11]               [,12]               [,13]               [,14]               [,15]             
# [1,] "23"                "27"                "23"                "5"                 "18"              
# [2,] "27"                "35"                "24"                "23"                "23"              
# [3,] "30"                "36"                "27"                "24"                "27"              
# [4,] "NR1H2"             "NR1I3"             "NR1H2"             "CEBPB"             "NCOR1"           
# [5,] "NR1I3"             "PPARD"             "NR1H3"             "NR1H2"             "NR1H2"           
# [6,] "NR3C1"             "PPARG"             "NR1I3"             "NR1H3"             "NR1I3"           
# [7,] "42.5030138588911"  "36.4553478035918"  "38.8848781081479"  "38.0074329405789"  "41.1619473641541"
# [8,] "30.1002250396887"  "30.2100824839673"  "30.2897226516113"  "30.6153874869959"  "30.6928551807896"
# [9,] "0.291809631674104" "0.171312734506655" "0.221041080098836" "0.194489469076742" "0.25433908874004"

##------------------------------------------- log.gene.tb ------------------------
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- log.gene.tb[,i]; return(interaction3gene(log.CYP2D6, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

dev.new()
Residual.Deviance <- results[2,]
hist(Residual.Deviance,breaks=100,main="log.CYP2D6 Residual Deviance")

summary(Residual.Deviance)
summary(Residual.Deviance.Deduction)

summary(Residual.Deviance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 27.12   46.58   54.66   54.17   62.54   77.99 
summary(Residual.Deviance.Deduction)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001347 0.0304900 0.0524600 0.0636700 0.0848900 0.2927000

goodfit.index <- which((Residual.Deviance<=29.275) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])

goodfitr[,order(goodfitr[8,])]
#     [,1]               [,2]                [,3]                [,4]                [,5]               
# [1,] "5"                "27"                "23"                "11"                "23"               
# [2,] "27"               "35"                "27"                "23"                "24"               
# [3,] "38"               "38"                "38"                "27"                "27"               
# [4,] "CEBPB"            "NR1I3"             "NR1H2"             "FOXA2"             "NR1H2"            
# [5,] "NR1I3"            "PPARD"             "NR1I3"             "NR1H2"             "NR1H3"            
# [6,] "RXRB"             "RXRB"              "RXRB"              "NR1I3"             "NR1I3"            
# [7,] "33.9774223458261" "33.6650897857361"  "35.1629455969535"  "34.9396826710584"  "34.7526975279639" 
# [8,] "27.1166314886214" "27.4660859629165"  "27.5621630812716"  "27.65051169708"    "27.7840976672473" 
# [9,] "0.20192205245515" "0.184137451059051" "0.216158867997122" "0.208621556257471" "0.200519682108399"
#     [,6]                [,7]                [,8]                [,9]                [,10]              
# [1,] "23"                "5"                 "5"                 "24"                "19"               
# [2,] "27"                "23"                "19"                "27"                "24"               
# [3,] "33"                "27"                "27"                "35"                "27"               
# [4,] "NR1H2"             "CEBPB"             "CEBPB"             "NR1H3"             "NCOR2"            
# [5,] "NR1I3"             "NR1H2"             "NCOR2"             "NR1I3"             "NR1H3"            
# [6,] "PGRMC1"            "NR1I3"             "NR1I3"             "PPARD"             "NR1I3"            
# [7,] "34.3476790205327"  "34.1345246237898"  "34.0325393636875"  "33.5575359545533"  "34.5766323573562" 
# [8,] "27.9352258077626"  "28.0366687537005"  "28.5157518814609"  "28.5213805485548"  "28.8010229104436" 
# [9,] "0.186692475172391" "0.178641886397899" "0.162103315984497" "0.150075244285484" "0.167037940167815"
#     [,11]               [,12]               [,13]               [,14]              [,15]              
# [1,] "19"                "19"                "23"                "11"               "11"               
# [2,] "27"                "27"                "27"                "27"               "27"               
# [3,] "38"                "33"                "41"                "32"               "43"               
# [4,] "NCOR2"             "NCOR2"             "NR1H2"             "FOXA2"            "FOXA2"            
# [5,] "NR1I3"             "NR1I3"             "NR1I3"             "NR1I3"            "NR1I3"            
# [6,] "RXRB"              "PGRMC1"            "USF1"              "ONECUT1"          "YY1"              
# [7,] "35.1816133570608"  "34.6084274503969"  "35.1429648700297"  "34.5932330861535" "34.9223569326372" 
# [8,] "28.8356643570794"  "28.9932879956229"  "29.002755303474"   "29.2230246296701" "29.2570030242699" 
# [9,] "0.180376861503648" "0.162247749130525" "0.174720874839795" "0.15523869778546" "0.162227134877963"

##----------------------------------------------- split peaks ---------------------------------------------

index.low <- which(Residual.Deviance<=35.75)
index.high <- which(Residual.Deviance>35.75)
pg.tb <- table(gene.names.m[,index.low])
npg.tb <- table(gene.names.m[,index.high])

pg.tb[order(pg.tb)]
# AHR    AHRR    ARNT   CEBPA   CEBPD   FOXA1   FOXA2   FOXA3   HNF4A   HNF4G   NCOA2   NCOA3   NCOR1   NR1D2 
# 41      41      41      41      41      41      41      41      41      41      41      41      41      41 
# NR1H4   NR2F1   NR2F2   NR3C1   NR5A2 ONECUT1   PPARA   PPARG    RXRA    THRA   NCOA1    THRB     VDR     DBP 
# 41      41      41      41      41      41      41      41      41      41      42      42      42      43 
# USF1    ESR1  NFE2L2   CEBPG   NR0B2   CEBPB   NR1H2   NR1I2  PGRMC1     YY1   NCOR2    RXRB   PPARD   NR1H3 
# 44      45      45      47      47      49      51      52      53      53      54      55      63      88 
# NR1I3 
# 861 

npg.tb[order(npg.tb)]
# NR1H3   PPARD    RXRB   NCOR2  PGRMC1     YY1   NR1I2   NR1H2   CEBPB   CEBPG   NR0B2    ESR1  NFE2L2    USF1 
# 773     798     806     807     808     808     809     810     812     814     814     816     816     817 
# DBP   NCOA1    THRB     VDR     AHR    AHRR    ARNT   CEBPA   CEBPD   FOXA1   FOXA2   FOXA3   HNF4A   HNF4G 
# 818     819     819     819     820     820     820     820     820     820     820     820     820     820 
# NCOA2   NCOA3   NCOR1   NR1D2   NR1H4   NR2F1   NR2F2   NR3C1   NR5A2 ONECUT1   PPARA   PPARG    RXRA    THRA 
# 820     820     820     820     820     820     820     820     820     820     820     820     820     820 

# For CYP2D6:
#   clear predictor: NR1I3
# potential predictor: none




## -----------------  Quadruplets Visualization via igraph
setwd("/Users/ronglu/Downloads/Danxin_Figure_Request_12242018/Danxin_Data")
Danxin.node <- read.csv("Gephi_gene_node.csv", header=T, stringsAsFactors = F)
Danxin.edge <- read.csv("CYP3A4_gene_4comb_Gephi_5 (Edges_filtered_by_SI).csv", header=T, stringsAsFactors = F)

library(igraph)
head(Danxin.node)
head(Danxin.edge)

# > head(Danxin.node)
# Id  Label          SI Residual.Dev.Perc Null.Dev      Weight
# 1    AHR    AHR 0.068721840         0.8439045 184.8586 0.081433195
# 2  AHR.1  AHR.1 0.073308624         0.8338740 184.8586 0.087913312
# 3   AHRR   AHRR 0.008624813         0.9842717 184.8586 0.008762634
# 4   ARNT   ARNT 0.067964971         0.8464045 184.8586 0.080298450
# 5 ARNT.1 ARNT.1 0.038734057         0.9134514 184.8586 0.042404070
# 6  CEBPA  CEBPA 0.007946078         0.9819607 184.8586 0.008092053

# > head(Danxin.edge)
# Source  Target     Type   id                             label timeset    weight        si
# 1   ESR1 HNF4A.2 Directed    5 ESR1 - HNF4A.2 - NFE2L2 - PPARA.1      NA 0.6928397 0.3095288
# 2   ESR1  NFE2L2 Directed 1039                                        NA 0.6886607 0.3045233
# 3   ESR1  NFE2L2 Directed 1045                                        NA 0.6903179 0.3038208
# 4   ESR1  NFE2L2 Directed 1068                                        NA 0.6836823 0.3011575
# 5   ESR1  NFE2L2 Directed 1069                                        NA 0.6837338 0.3013028
# 6   ESR1  NFE2L2 Directed 1082 ESR1 - NFE2L2 - HNF4A.3 - PPARA.1      NA 0.6897985 0.3033020

# residual_dev_perc nul_.deviance
# 1         0.3071603      184.8586
# 2         0.3113393      184.8586
# 3         0.3096821      184.8586
# 4         0.3163177      184.8586
# 5         0.3162662      184.8586
# 6         0.3102015      184.8586

dim(Danxin.node) # [1] 78  6
length(unique(Danxin.node$Id))  ## [1] 78

which(Danxin.node$Id=="DBP")  # [1] 11
which(Danxin.node$Id=="YY1")  # [1] 78
Danxin.net <- graph_from_data_frame(d=Danxin.edge, vertices=Danxin.node[-c(11, 78),], directed=F) 
Danxin.net

plot(Danxin.net)

# Generate colors based on media type:
names(Danxin.node)
# [1] "Id"                "Label"             "SI"                "Residual.Dev.Perc" "Null.Dev"         
# [6] "Weight" 
unique.label <- unlist(lapply(strsplit(V(Danxin.net)$Label, "[.]"), function(x){return(x[1])}))
length(unique(unique.label))  ## [1] 43
library("RColorBrewer")
Danxin.colrs <- rainbow(43, alpha=.5) 
V(Danxin.net)$color <- Danxin.colrs[as.integer(as.factor(unique.label))]

## use univariate main-effect SI to set node size:
V(Danxin.net)$size <- V(Danxin.net)$SI*40
V(Danxin.net)$label.cex <- V(Danxin.net)$SI*5
E(Danxin.net)$label.cex <- rep(0.001, length(E(Danxin.net)))
E(Danxin.net)$size <- E(Danxin.net)$si
plot(Danxin.net)  ## NetVisual_6

V(Danxin.net)$label.cex <- V(Danxin.net)$SI*4  ## NetVisual_7

V(Danxin.net)$frame.color <- "white"  ## NetVisual_8

V(Danxin.net)$label.cex <- V(Danxin.net)$SI*6  ## NetVisual_9
plot(Danxin.net)

V(Danxin.net)$label.cex <- V(Danxin.net)$SI*4  ## NetVisual_10
l <- layout_with_fr(Danxin.net)
plot(Danxin.net, layout=l)

l <- layout_with_kk(Danxin.net)  ## ## NetVisual_11
plot(Danxin.net, layout=l)

l <- layout_in_circle(Danxin.net)
plot(Danxin.net, layout=l)  


## in order to investigate tiny circles in NetVisual_7, we need to fix node label size:
V(Danxin.net)$label.cex <- rep(0.5, length(V(Danxin.net)$SI))
l <- layout_with_fr(Danxin.net)
plot(Danxin.net, layout=l)   ## NetVisual_12


##  -------------------------- work of 12/26/2018 starts from this line: 

## Prof. Wolfgang commented: the labeling is quite small. Maybe ok online, but not in print:
## so ignore all nodes that is not connected:
connected.node <- unique(c(Danxin.edge$Source, Danxin.edge$Target))
connected.node
# [1] "ESR1"    "HNF4A.2" "NFE2L2"  "HNF4A.3" "PPARA.1" "NR1I3"   "NR1I3.1" "PPARA.4" "THRB"    "NR1I2"  
# [11] "ESR1.2"  "NR1I2.1" "VDR.2"   "PGRMC1"  "FOXA2"   "USF1"    "RXRA.1"  "NR2F2.1" "NR2F2"   "NCOR2.2"
# [21] "VDR"     "PPARA.3" "RXRB"    "FOXA1.1" "PPARA"   "PPARG.2" "NR1D2.1" "PPARA.2" "ARNT" 
connected.node.index <- sapply(connected.node, function(x){which(Danxin.node$Id==x)})
Danxin.net.2 <- graph_from_data_frame(d=Danxin.edge, vertices=Danxin.node[connected.node.index,], directed=F) 

unique.label <- unlist(lapply(strsplit(V(Danxin.net.2)$Label, "[.]"), function(x){return(x[1])}))
length(unique(unique.label))  ## [1] 19
library("RColorBrewer")
Danxin.colrs <- rainbow(19, alpha=.5) 
V(Danxin.net.2)$color <- Danxin.colrs[as.integer(as.factor(unique.label))]

## use univariate main-effect SI to set node size:
V(Danxin.net.2)$size <- V(Danxin.net.2)$SI*120
#V(Danxin.net.2)$label.cex <- V(Danxin.net.2)$SI*5
V(Danxin.net.2)$label.cex <- rep(1, length(V(Danxin.net.2)$SI))
V(Danxin.net.2)$label.degree <- 0
E(Danxin.net.2)$label.cex <- rep(0.001, length(E(Danxin.net.2)))
E(Danxin.net.2)$size <- E(Danxin.net.2)$si
V(Danxin.net.2)$frame.color <- "white"  
plot(Danxin.net.2)  ## NetVisual_6

save.image("Danxin_Quadruplets_Visualization_image_12262018.RData")


V(Danxin.net.2)$color <- Danxin.colrs[7]   ##   "lightskyblue"   ##   grep("blue", colors(), value=T)[4]
V(Danxin.net.2)$size <- V(Danxin.net.2)$SI*50
V(Danxin.net.2)$label.cex <- rep(1.25, length(V(Danxin.net.2)$SI))
V(Danxin.net.2)$label.cex <- rep(1.75, length(V(Danxin.net.2)$SI))
plot(Danxin.net.2)
l <- layout_with_fr(Danxin.net.2)
plot(Danxin.net.2, layout=l)

V(Danxin.net.2)$size <- V(Danxin.net.2)$SI*60
V(Danxin.net.2)$label.cex <- rep(3.5, length(V(Danxin.net.2)$SI))
pdf("Rplot25.pdf", 30, 30)
plot(Danxin.net.2)
dev.off()
