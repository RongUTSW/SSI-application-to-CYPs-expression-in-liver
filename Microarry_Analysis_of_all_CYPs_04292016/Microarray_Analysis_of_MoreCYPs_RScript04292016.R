## reference script: "CYP3A4_Microarray_RScript_withMoreCYPs 04 29 2016"

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

