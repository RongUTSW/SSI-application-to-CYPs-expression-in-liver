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



##------------------------------------------------------- 02 23 2016

##y1 <- CYPs.sorted.m[,30] ## CYP3A4

index.2m.all <- combn(c(1:78), 2)
dim(index.2m.all)

gene.names.2m.all <- apply(index.2m.all, 2, function(i){gene.names[i]})

all.gene.pairs.3 <- apply(index.2m.all, 2, function(i){SI.comb(y1, input.genes.m[,i], k=3)})

order.2m.3 <- order(-all.gene.pairs.3[2,], all.gene.pairs.3[1,],decreasing=TRUE)

gene.names.2m.all[,order.2m.3][,1:30]
all.gene.pairs.3[,order.2m.3][,1:30]

##------------------------------------------------------- 02 23 2016
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


> save.image("U:\\Data\\Danxin\\CYP3A4_Microarray\\gene_triplets_picks_02_23_2016_fullPoly3.Rdata")

dev.new();hist(all.gene.triplets.3[2,], breaks=1200)
dev.new();hist(all.gene.triplets.3[1,order.3m.3][1:200], breaks=50)

##-------------------------------------------------------------------------



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

SI.comb <- function(y, input, k=3){
  main.new.m <- apply(input, 2, function(x){ apply(matrix(1:k),1, function(i){x^i}) })
  main.SI <- apply(main.new.m, 2, function(x){SI(y,matrix(x, nr=dim(input)[1]))})
  
  comb.new.m <- cbind(matrix(main.new.m, nr=dim(input)[1]), comb.m(input, k))
  comb.SI <- SI(y, comb.new.m)
  
  return(cbind(comb.SI, main.SI))
}

y <- CYPs.sorted.m[,30]
input <- apply(input.genes.m, 2, function(x){ifelse(is.na(x), mean(x, na.rm=TRUE), x)})

##------------------------------------------------------- 02 23 2016

##y1 <- CYPs.sorted.m[,30] ## CYP3A4

index.2m.all <- combn(c(1:78), 2)
dim(index.2m.all)

gene.names.2m.all <- apply(index.2m.all, 2, function(i){gene.names[i]})

all.gene.pairs.3 <- apply(index.2m.all, 2, function(i){SI.comb(y1, input[,i], k=3)})

order.2m.3 <- order(-all.gene.pairs.3[2,], all.gene.pairs.3[1,],decreasing=TRUE)

gene.names.2m.all[,order.2m.3][,1:30]
all.gene.pairs.3[,order.2m.3][,1:30]

all.gene.pairs.3.stepwise <- all.gene.pairs.3
order.2m.3.stepwise <- order.2m.3

gene.names.2m.all[,order.2m.3.stepwise][,1:30]
all.gene.pairs.3.stepwise[,order.2m.3.stepwise][,1:30]

> dim(index.2m.all)
[1]    2 3003

##------------------------------------------------------- 02 23 2016
index.3m.all <- combn(c(1:78), 3)
dim(index.3m.all)

gene.names.3m.all <- apply(index.3m.all, 2, function(i){gene.names[i]})

all.gene.triplets.3.stepwise <- apply(index.3m.all, 2, function(i){SI.comb(y1, input[,i], k=3)})

order.3m.3.stepwise <- order(-all.gene.triplets.3.stepwise[2,], all.gene.triplets.3.stepwise[1,],decreasing=TRUE)

gene.names.3m.all[,order.3m.3.stepwise][,1:30]
all.gene.triplets.3.stepwise[,order.3m.3.stepwise][,1:30]

table(gene.names.2m.all[,order.2m.3.stepwise][,1:100])
table(gene.names.3m.all[,order.3m.3.stepwise][,1:100])

> save.image("U:\\Data\\Danxin\\CYP3A4_Microarray\\gene_triplets_picks_02_23_2016_Poly3stepwise.Rdata")

> gn.3 <- gene.names.3m.all[,order.3m.3.stepwise][,1:30]
> si.3 <- all.gene.triplets.3.stepwise[,order.3m.3.stepwise][,1:30]
> df.3 <- data.frame(gn.3[1,], rep("#",30), gn.3[2,], rep("#",30), gn.3[3,], rep("&", 30), si.3[1,]*10, rep("@",30))
> df.3[1:5,]

write.csv(df.3, "U:/Data/Danxin/CYP3A4_Microarray/graphviz_Poly3step_triplets.csv")

gn.2 <- gene.names.2m.all[,order.2m.3.stepwise][,1:30]
si.2 <- all.gene.pairs.3.stepwise[,order.2m.3.stepwise][,1:30]
df.2 <- data.frame(gn.2[1,], rep("#",30), gn.2[2,], rep("&", 30), si.2[1,]*10, rep("@",30))
df.2[1:5,]

write.csv(df.2, "U:/Data/Danxin/CYP3A4_Microarray/graphviz_Poly3step_pairs.csv")

gn.1 <- matrix(c(gn.3[1,], gn.3[2,], gn.3[3,], gn.2[1,], gn.2[2,]), nc=1)
si.1 <- matrix(c(si.3[4,], si.3[7,], si.3[10,], si.2[4,], si.2[7,]), nc=1)

df.1 <- data.frame(gn.1, rep("&", 30), si.1, rep("@",30))
df.1[1:5,]
df.1 <- df.1[order(gn.1),]

write.csv(df.1[!duplicated(df.1[,1]),], "U:/Data/Danxin/CYP3A4_Microarray/graphviz_Poly3step_single.csv")

##------------ plot 2nd-order analysis results only

gn.2 <- gene.names.2m.all[,order.2m.3.stepwise][,1:100]
si.2 <- all.gene.pairs.3.stepwise[,order.2m.3.stepwise][,1:100]
df.2 <- data.frame(gn.2[1,], rep("#",100), gn.2[2,], rep("&", 100), si.2[1,]*10, rep("@",100))
df.2[1:5,]

write.csv(df.2, "U:/Data/Danxin/CYP3A4_Microarray/graphviz_Poly3step_pairs_100.csv")

gn.1 <- matrix(c(gn.2[1,], gn.2[2,]), nc=1)
si.1 <- matrix(c(si.2[4,], si.2[7,]), nc=1)

df.1 <- data.frame(gn.1, rep("&", 100), si.1, rep("@",100))
df.1[1:5,]
df.1 <- df.1[order(gn.1),]

write.csv(df.1[!duplicated(df.1[,1]),], "U:/Data/Danxin/CYP3A4_Microarray/graphviz_Poly3step_single_100.csv")


##------------ plot 3rd-order analysis results only

gn.3 <- gene.names.3m.all[,order.3m.3.stepwise][,1:100]
si.3 <- all.gene.triplets.3.stepwise[,order.3m.3.stepwise][,1:100]
df.3 <- data.frame(gn.3[1,], rep("#",100), gn.3[2,], rep("#",100), gn.3[3,], rep("&", 100), si.3[1,]*10, rep("@",100))
df.3[1:5,]

write.csv(df.3, "U:/Data/Danxin/CYP3A4_Microarray/graphviz_Poly3step_triplets_100.csv")

gn.1 <- matrix(c(gn.3[1,], gn.3[2,], gn.3[3,]), nc=1)
si.1 <- matrix(c(si.3[4,], si.3[7,], si.3[10,]), nc=1)

df.1 <- data.frame(gn.1, rep("&", 100), si.1, rep("@",100))
df.1[1:5,]
df.1 <- df.1[order(gn.1),]

write.csv(df.1[!duplicated(df.1[,1]),], "U:/Data/Danxin/CYP3A4_Microarray/graphviz_Poly3step_tri_single_100.csv")


##------------------------------------------------------- 02 24 2016

all.gene.pairs.10.stepwise <- apply(index.2m.all, 2, function(i){SI.comb(y1, input[,i], k=10)})

order.2m.10.stepwise <- order(-all.gene.pairs.10.stepwise[2,], all.gene.pairs.10.stepwise[1,],decreasing=TRUE)

gene.names.2m.all[,order.2m.10.stepwise][,1:30]
all.gene.pairs.10.stepwise[,order.2m.10.stepwise][,1:30]


all.gene.triplets.10.stepwise <- apply(index.3m.all, 2, function(i){SI.comb(y1, input[,i], k=10)})

order.3m.10.stepwise <- order(-all.gene.triplets.10.stepwise[2,], all.gene.triplets.10.stepwise[1,],decreasing=TRUE)

gene.names.3m.all[,order.3m.10.stepwise][,1:30]
all.gene.triplets.10.stepwise[,order.3m.10.stepwise][,1:30]

save.image("U:\\Data\\Danxin\\CYP3A4_Microarray\\gene_triplets_picks_02_24_2016_Poly10stepwise.Rdata")


table(gene.names.2m.all[,order.2m.10.stepwise][,1:100])
table(gene.names.2m.all[,order.2m.10.stepwise][,1:100])

table(gene.names.3m.all[,order.3m.10.stepwise][,1:200])
table(gene.names.3m.all[,order.3m.10.stepwise][,1:200])


##------------------------------------------------------------------------------------------------
DBP.m <- apply(matrix(1:10, nc=1),1, function(i){input.genes[,11]^i}) 
summary(y1~., data=data.frame(DBP.m))

all.gene.pairs.2 <- apply(index.2m.all, 2, function(i){SI.comb(y1, input.genes.m[,i], k=2)})
order.2m.2 <- order(all.gene.pairs.2[1,],decreasing=TRUE)

gene.names.2m.all[,order.2m.2][,1:30]
all.gene.pairs.2[,order.2m.2][,1:30]





##-------------------------------------------------------

y1 <- CYPs.sorted.m[,30]
index.2m <- combn(sig.uni.index, 2) 
all.gene.pairs <- apply(index.2m, 2, function(i){SI.comb(y1, input.genes.m[,i])})
gene.names.2m <- apply(index.2m, 2, function(i){gene.names[i]})
dim(all.gene.pairs)
order.2m <- order(all.gene.pairs[1,],decreasing=TRUE)


> gene.pairs.pick <- gene.names.2m[,order.2m][,1:30]
all.gene.pairs[,order.2m][,1:30]

> gene.triplet.pick <- gene.names.3m[,order.3m][,1:30]

> xxxx_2 <- cbind(gene.pairs.pick[1,], rep("@",30), gene.pairs.pick[2,], rep("#",30), pair.ef, rep("$",30), round(pair.ef/100, 3), rep("&",30))
> xxxx_2[1:3,]
pair.ef                  
[1,] "NCOR2.2" "@" "DBP" "#" "319.865" "$" "3.199" "&"
[2,] "NR1D2.1" "@" "DBP" "#" "297.004" "$" "2.97"  "&"
[3,] "NCOR2"   "@" "DBP" "#" "275.45"  "$" "2.754" "&"
> xxxx_3 <- cbind(gene.triplet.pick[1,], rep("@",30), gene.triplet.pick[2,], rep("@",30),gene.triplet.pick[3,], rep("#",30), round(all.gene.triplets[5,order.3m][1:30]/100, 3), rep("$",30))
> xxxx_3[1:3,]
[,1]      [,2] [,3]      [,4] [,5]  [,6] [,7]     [,8]
[1,] "NCOR2"   "@"  "NR1D2.1" "@"  "DBP" "#"  "13.554" "$" 
[2,] "NCOR2"   "@"  "FOXA3"   "@"  "DBP" "#"  "8.632"  "$" 
[3,] "NCOR2.2" "@"  "NR1D2.1" "@"  "DBP" "#"  "8.526"  "$" 
> write.table(xxxx_3, "U:/graph_3.txt",sep="\t")
> write.table(xxxx_2, "U:/graph_2.txt",sep="\t")
> write.table(xxxx_1, "U:/graph_1.txt",sep="\t")




> save.image("U:\\Data\\Danxin\\CYP3A4_Microarray\\gene_triplets_picks_02 14 2016")
> quantile(all.gene.pairs[1,], c(1:99)/100)
1%          2%          3%          4%          5%          6%          7%          8% 
  0.03079295  0.03406651  0.03705681  0.03880997  0.04100650  0.04222772  0.04388659  0.04630706 
9%         10%         11%         12%         13%         14%         15%         16% 
  0.04918596  0.05125508  0.05271500  0.05458275  0.05649468  0.05795480  0.06034759  0.06290588 
17%         18%         19%         20%         21%         22%         23%         24% 
  0.06535383  0.06827969  0.07001938  0.07118430  0.07353976  0.07487934  0.07642496  0.07790125 
25%         26%         27%         28%         29%         30%         31%         32% 
  0.07968092  0.08185640  0.08276226  0.08387796  0.08524269  0.08694749  0.08848416  0.09048108 
33%         34%         35%         36%         37%         38%         39%         40% 
  0.09183981  0.09373995  0.09604224  0.09758099  0.10011823  0.10219187  0.10376550  0.10741549 
41%         42%         43%         44%         45%         46%         47%         48% 
  0.11053127  0.11338611  0.11567073  0.11818779  0.11974271  0.12134061  0.12306466  0.12490719 
49%         50%         51%         52%         53%         54%         55%         56% 
  0.12637441  0.12814393  0.13006095  0.13158420  0.13334655  0.13461691  0.13685847  0.13900013 
57%         58%         59%         60%         61%         62%         63%         64% 
  0.14012451  0.14197704  0.14309229  0.14454435  0.14592572  0.14695614  0.14798705  0.14951681 
65%         66%         67%         68%         69%         70%         71%         72% 
  0.15100778  0.15307506  0.15461402  0.15597065  0.15867126  0.16091208  0.16321747  0.16739373 
73%         74%         75%         76%         77%         78%         79%         80% 
  0.17075783  0.17444451  0.17944197  0.18529143  0.18794819  0.19072887  0.19347271  0.19497602 
81%         82%         83%         84%         85%         86%         87%         88% 
  0.19663776  0.19898225  0.20035869  0.20403333  0.20681200  0.21063764  0.21381570  0.22273597 
89%         90%         91%         92%         93%         94%         95%         96% 
  0.22668457  0.22861898  0.23031025  0.23217361  0.23547608  0.23829550  0.24247710  0.25104566 
97%         98%         99% 
0.26196990 11.35879668 32.56751552 
> length(which(all.gene.triplets[1,]>5))
[1] 716
> length(which(all.gene.pairs[1,]>11))
[1] 28


gene.names.2m <- apply(index.2m, 2, function(i){gene.names[i]})

> y1 <- CYPs.sorted.m[,30]
> index.2m <- combn(sig.uni.index, 2) 
> all.gene.pairs <- apply(index.2m, 2, function(i){SI.comb(y1, input.genes.m[,i])})
> 
  > dim(all.gene.pairs)
[1]    4 1378
> order.2m <- order(all.gene.pairs[1,],decreasing=TRUE)
> save.image("U:\\Data\\Danxin\\CYP3A4_Microarray\\gene_triplets_picks_02 14 2016")


##------
index.3m <- combn(c(1:78), 3) 

sig.uni.index <- order(cor.p.fdr, decreasing=F)[1:53]
index.3m <- combn(sig.uni.index, 3) 

names(CYPs.sorted)[30:32]

y1 <- CYPs.sorted.m[,30]
all.gene.triplets <- apply(index.3m, 2, function(i){SI.comb(y1, input.genes.m[,i])})
all.gene.triplets[,order.3m][,1:30]


> summary(all.gene.triplets[1,])
Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0334    0.1312    0.1713    3.4680    0.2183 1362.0000 
> dev.new(); hist(all.gene.triplets, breaks=200, xlim=c(-1,20))
> dev.new(); hist(all.gene.triplets[1,which(all.gene.triplets[1,]<20)], breaks=200)
> order(all.gene.triplets[1,],decreasing=TRUE)[1:10]
[1] 14074 14278  8033  8237 13788 19937  7715  4921 20288  4789
> all.gene.triplets[1,order(all.gene.triplets[1,],decreasing=TRUE)][1:10]
[1] 1361.7909  869.5875  859.0480  855.9750  788.0847  767.9188  673.3102  603.1572  602.7381
[10]  602.7303
> length(which(all.gene.triplets[1,]>5))
[1] 716

gene.names.3m <- apply(index.3m, 2, function(i){gene.names[i]})

> order.3m <- order(all.gene.triplets[1,],decreasing=TRUE)
> gene.names.3m[,order.3m[1:10]]
[,1]      [,2]    [,3]      [,4]      [,5]    [,6]      [,7]      [,8]      [,9]     
[1,] "NCOR2"   "NCOR2" "NCOR2.2" "NCOR2.2" "NCOR2" "NR5A2"   "NCOR2.2" "NR1I3.1" "PPARD"  
[2,] "NR1D2.1" "FOXA3" "NR1D2.1" "FOXA3"   "NR1H2" "NR1D2.1" "THRB"    "NR1D2"   "NR1D2.1"
[3,] "DBP"     "DBP"   "DBP"     "DBP"     "DBP"   "DBP"     "DBP"     "DBP"     "DBP"    
[,10]    
[1,] "NR1I3.1"
[2,] "NR1D2.1"
[3,] "DBP"    
> table(gene.names.3m[,order.3m[1:716]])

ARNT  ARNT.1   CEBPD   CEBPG CEBPG.1     DBP    ESR1  ESR1.2 FOXA1.1   FOXA2   FOXA3   HNF4A 
34      35      31      27      29     716      28      27      32      22      35      32 
HNF4A.2   NCOA1   NCOR1   NCOR2 NCOR2.1 NCOR2.2  NFE2L2   NR0B1   NR0B2   NR1D2 NR1D2.1   NR1H2 
30      32      26      42      37      42      31      29      35      32      42      39 
NR1H2.1   NR1H3   NR1I3 NR1I3.1   NR5A2  PGRMC1 PPARA.1 PPARA.5   PPARD PPARD.1 PPARD.2    RXRA 
33      34      40      42      34      36      33      31      39      31      36      34 
RXRA.1    RXRB    THRA    THRB    USF1     VDR   VDR.1   VDR.2 
34      30      31      36      32      31      33      33 
> tb.3m <- table(gene.names.3m[,order.3m[1:716]])
> tb.3m[order(tb.3m)]

FOXA2   NCOR1   CEBPG  ESR1.2    ESR1 CEBPG.1   NR0B1 HNF4A.2    RXRB   CEBPD  NFE2L2 PPARA.5 
22      26      27      27      28      29      29      30      30      31      31      31 
PPARD.1    THRA     VDR FOXA1.1   HNF4A   NCOA1   NR1D2    USF1 NR1H2.1 PPARA.1   VDR.1   VDR.2 
31      31      31      32      32      32      32      32      33      33      33      33 
ARNT   NR1H3   NR5A2    RXRA  RXRA.1  ARNT.1   FOXA3   NR0B2  PGRMC1 PPARD.2    THRB NCOR2.1 
34      34      34      34      34      35      35      35      36      36      36      37 
NR1H2   PPARD   NR1I3   NCOR2 NCOR2.2 NR1D2.1 NR1I3.1     DBP 
39      39      40      42      42      42      42     716 
> gene.names.3m[,order.3m][,1:10]
[,1]      [,2]    [,3]      [,4]      [,5]    [,6]      [,7]      [,8]      [,9]     
[1,] "NCOR2"   "NCOR2" "NCOR2.2" "NCOR2.2" "NCOR2" "NR5A2"   "NCOR2.2" "NR1I3.1" "PPARD"  
[2,] "NR1D2.1" "FOXA3" "NR1D2.1" "FOXA3"   "NR1H2" "NR1D2.1" "THRB"    "NR1D2"   "NR1D2.1"
[3,] "DBP"     "DBP"   "DBP"     "DBP"     "DBP"   "DBP"     "DBP"     "DBP"     "DBP"    
[,10]    
[1,] "NR1I3.1"
[2,] "NR1D2.1"
[3,] "DBP"    
> gene.names.3m[,order.3m][,1:10]
[,1]      [,2]    [,3]      [,4]      [,5]    [,6]      [,7]      [,8]      [,9]      [,10]    
[1,] "NCOR2"   "NCOR2" "NCOR2.2" "NCOR2.2" "NCOR2" "NR5A2"   "NCOR2.2" "NR1I3.1" "PPARD"   "NR1I3.1"
[2,] "NR1D2.1" "FOXA3" "NR1D2.1" "FOXA3"   "NR1H2" "NR1D2.1" "THRB"    "NR1D2"   "NR1D2.1" "NR1D2.1"
[3,] "DBP"     "DBP"   "DBP"     "DBP"     "DBP"   "DBP"     "DBP"     "DBP"     "DBP"     "DBP"    
> all.gene.triplets[,order.3m][,1:10]
[,1]         [,2]         [,3]         [,4]         [,5]         [,6]         [,7]
[1,] 1.361791e+03 869.58746377 859.04802232 855.97498964 788.08471035 767.91878319 673.31015309
[2,] 1.080161e-01   0.10801612   0.13547415   0.13547415   0.10801612   0.06264377   0.13547415
[3,] 3.416964e-02   0.01602659   0.03416964   0.01602659   0.06966843   0.03416964   0.09952252
[4,] 6.242475e+00   6.24247535   6.24247535   6.24247535   6.24247535   6.24247535   6.24247535
[5,] 1.355406e+03 863.22094572 852.63590318 849.58101356 781.66455046 761.57949444 666.83268107
[,8]         [,9]        [,10]
[1,] 603.15722911 602.73813436 602.73026663
[2,]   0.18976590   0.06109898   0.18976590
[3,]   0.02301699   0.03416964   0.03416964
[4,]   6.24247535   6.24247535   6.24247535
[5,] 596.70197087 596.40039039 596.26385573







##-------------------------------------------------- before 02 14 2016 ------------------------------------------------

load("U:/Data/Danxin/CYP3A4_Microarray/CYP3A4_Microarray_input_CYPs_separated.Rdata")
## CYPs.sorted, input.genes

CYPs.sorted.m <- as.matrix(CYPs.sorted)
input.genes.m <- as.matrix(input.genes)

library(MASS)
library(glmnet)
library(multcomp) ## multiT, multinormal
library(sgof)	## Beta-binomial SGoF procedure, for correcting positive correlation

## caculation of the numerator of the sobol main index

t.gs <- function(i, coef, mse, rhom){
  rho <- rhom[i,-i]
  return ( ( coef[i]+ matrix(coef[-i],nr=1) %*% ( matrix(mse[-i],nc=1) * matrix(rho,nc=1) ) / mse[i] )^2 * mse[i]^2  )
}

## y: response vector
## input: input variable matrix

glm.gs.full <- function(y, input){
  ## fit the sample with GLM, using iteratively reweighted least squares (IRLS)
  fit.1 <- glm(y~., data=data.frame(input))
  glm.coef.1 <- fit.1$coef
  glm.p.1 <- summary(fit.1)$coef[,4]; 
  p.fdr.irls <- p.adjust(glm.p.1[-1],"fdr")
  
  p.hochberg.irls <- p.adjust(glm.p.1[-1],"hochberg")
  
  p.hommel.irls <- p.adjust(glm.p.1[-1],"hommel")
  
  
  K <- diag( (dim(input)[2]+1) )[-1,]
  glnt_lmod <- glht(fit.1, linfct = K)
  p.multT.irls <- summary(glnt_lmod)$test[["pvalues"]][1:(dim(input)[2])] 
  
  ## res <- BBSGoF(glm.p.1[-1], adjusted.pvalues=T, blocks=2)
  ## p.BBSGoF.irls <- res$Adjusted.pvalues
  
  ## all types of p values have no value for intercept
  
  ## calculate emprical estimate of sobol ranking
  coef.est.1 <- glm.coef.1[-1]
  mse.est.1 <- apply(input, 2, function(x){var(x)^0.5})
  rhom.est.1 <- cor(input)
  
  est.gs.1 <- apply(matrix(1:(dim(input)[2]), nc=1), 1, function(i){t.gs(i, coef.est.1, mse.est.1, rhom.est.1)})  
  
  ##est.gs.order.1 <- order(est.gs.1, decreasing=TRUE, na.last=TRUE)
  
  return( list(p.irls=glm.p.1[-1], p.fdr.irls=p.fdr.irls, p.hochberg.irls=p.hochberg.irls, p.hommel.irls=p.hommel.irls, p.multT.irls=p.multT.irls, gs.irls=est.gs.1) )
}


glm.gs <- function(y, input){
  ## fit the sample with GLM, using iteratively reweighted least squares (IRLS)
  fit.1 <- glm(y~., data=data.frame(input))
  glm.coef.1 <- fit.1$coef
  
  ## all types of p values have no value for intercept
  
  ## calculate emprical estimate of sobol ranking
  coef.est.1 <- glm.coef.1[-1]
  mse.est.1 <- apply(input, 2, function(x){var(x)^0.5})
  rhom.est.1 <- cor(input)
  
  gs.irls <- apply(matrix(1:(dim(input)[2]), nc=1), 1, function(i){t.gs(i, coef.est.1, mse.est.1, rhom.est.1)})  
  
  return( gs.irls )
}

CYP.scan <- function(y, input){
  
  ##---------------- rank inputs without artificial noises when the correct model is unknown -------------------------
  fit.glm.gs.1 <- glm.gs.full(y, input)
  
  d <- dim(input)[2]
  n <- dim(input)[1]
  sn <- 5000
  
  noise.index <- sample(1:d, sn, replace=T)
  noise <- apply(matrix(noise.index, nc=1), 1, function(i){sample(input[,i]) })
  ##col.is.a.match <- unlist(apply(matrix(noise.index, nc=1), 1, function(i){ which(apply(noise, 2, identical, input[,i])) }) )
  ##if(length(col.is.a.match)>0) noise <- noise[,-col.is.a.match]
  
  ## generate artificial noise with similar mean and variance
  ##noise.mu.est <- quantile(apply(input,2,function(x){mean(x,na.rm=T)}),c(1:10)/10)[c(1,10)]
  ##noise.se.est <- quantile(apply(input,2,function(x){var(x,na.rm=T)^0.5}),c(1:10)/10)[c(1,10)]
  ##noise.rho.est <- quantile( cor(input, use="pairwise.complete.obs"),c(1:10)/10)[c(1,10)]
  
  ##noise.mu <- matrix(runif(1000, min=noise.mu.est[1], max=noise.se.est[2]),nc=1)
  ##noise.sd <- matrix(runif(1000, min=noise.se.est[1], max=noise.se.est[2]),nc=1)
  ##noise <- apply(cbind(noise.mu, noise.sd), 1, function(x){rnorm(n,x[1],x[2])})
  
  fit.glm.gs.noise <- apply(noise, 2, function(x){glm.gs(y,cbind(input, x))})
  
  gs.irls.m <- t(fit.glm.gs.noise)
  
  gs.irls.table <- table( unlist(apply(gs.irls.m, 1, function(x){which(x>x[d+1])}) ) )
  
  gs.irls.percentage <- gs.irls.table/sn
  
  gs.percentage.list <- as.integer(names(which(gs.irls.percentage>0.95)))
  
  my.list <- list(wrong.full.model=fit.glm.gs.1, gs.irls.percentage=gs.irls.percentage, gs.percentage.list=gs.percentage.list)
  
  return(my.list)
  
}

##---------------------- analysis code
ptm <- proc.time()
res <- CYP.scan(y, input)
proc.time()-ptm


## --------- above are standard function code ---- the next step: customize the analysis code part
y <- CYPs.sorted.m[,30]
input <- apply(input.genes.m, 2, function(x){ifelse(is.na(x), mean(x, na.rm=TRUE), x)})



> names(CYPs.sorted)[30:32]
[1] "CYP3A4"   "CYP3A4.1" "CYP3A4.2"

> names(fit.glm.gs.1)
[1] "p.irls"          "p.fdr.irls"      "p.hochberg.irls"
[4] "p.hommel.irls"   "p.multT.irls"    "gs.irls" 


dev.new();par(mfrow=c(2,3));lapply(fit.glm.gs.1, function(x){plot(x);abline(h=0.05,col="green");});

lapply(fit.glm.gs.1, function(x){gene.names[order(x,decreasing=F)]})

gene.names[order(fit.glm.gs.1[[1]],decreasing=F)]

lapply(fit.glm.gs.1, function(x){gene.names[which(x<0.05)]})
> lapply(fit.glm.gs.1, function(x){gene.names[which(x<0.05)]})
$p.irls
[1] "CEBPG.1" "ESR1.1"  "HNF4A.1" "NCOR2"   "NCOR2.2" "NFE2L2" 
[7] "NR1D2.1" "NR1H2.1" "NR1H3"   "PPARA.1"

$p.fdr.irls
[1] "NFE2L2"  "NR1D2.1"

$p.hochberg.irls
[1] "NFE2L2"  "NR1D2.1"

$p.hommel.irls
[1] "NFE2L2"  "NR1D2.1"

$p.multT.irls
[1] "NFE2L2"  "NR1D2.1"

$gs.irls
[1] "AHRR"    "ARNT.1"  "CEBPA"   "CEBPB"   "CEBPD"   "CEBPG.1"
[7] "DBP"     "ESR1.1"  "FOXA1"   "FOXA1.1" "FOXA3"   "HNF4A"  
[13] "HNF4A.1" "HNF4G"   "HNF4G.1" "NCOA1"   "NCOA2"   "NCOA3"  
[19] "NCOA3.1" "NCOR2.1" "NR0B1"   "NR0B2"   "NR1D2"   "NR1D2.1"
[25] "NR1H2.1" "NR1H4"   "NR2F1"   "NR2F1.1" "NR2F2"   "NR2F2.1"
[31] "NR3C1"   "ONECUT1" "PPARA"   "PPARA.2" "PPARA.3" "PPARA.4"
[37] "PPARA.5" "PPARD.1" "PPARD.2" "PPARG"   "PPARG.1" "PPARG.2"
[43] "RXRA"    "RXRB"    "RXRB.1"  "RXRG"    "RXRG.1"  "THRA"   
[49] "VDR"     "VDR.1"   "VDR.2"   "YY1"    


> Sobol.p[,order(Sobol.p, decreasing=F)]
AHR   AHR.1    ARNT  ARNT.1   CEBPD   CEBPG CEBPG.1    ESR1 
0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 
ESR1.2   FOXA1   FOXA2 HNF4A.2 HNF4A.3   NCOA1   NCOR1   NCOR2 
0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 
NCOR2.1 NCOR2.2  NFE2L2   NR0B2   NR1D2 NR1D2.1   NR1H2 NR1H2.1 
0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 
NR1H3   NR1I2 NR1I2.1   NR1I3 NR1I3.1   NR5A2  PGRMC1 PPARA.1 
0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 
PPARD PPARD.1 PPARD.2  RXRA.1    THRA    THRB    USF1     VDR 
0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 
VDR.1   VDR.2    RXRA   NR0B1   NR2F2 NR2F2.1    RXRB ONECUT1 
0.0000  0.0000  0.0002  0.0006  0.0016  0.0016  0.0020  0.0172 
RXRB.1 PPARA.5   CEBPB FOXA1.1   FOXA3 PPARG.2   CEBPA   NR3C1 
0.0174  0.0190  0.0234  0.0234  0.0294  0.0390  0.0406  0.0530 
PPARG.1   HNF4A NR2F1.1    AHRR   NR2F1   NCOA2 PPARA.2 NCOA3.1 
0.0596  0.0738  0.0852  0.0906  0.1122  0.1138  0.1182  0.1432 
HNF4G  ESR1.1     DBP HNF4A.1   NR1H4 PPARA.3   NCOA3 HNF4G.1 
0.1498  0.1994  0.2022  0.2072  0.2194  0.2202  0.2908  0.3732 
PPARA.4   PPARG   PPARA    RXRG  RXRG.1     YY1 
0.3962  0.4116  0.4850  0.5578  0.5650  0.7338 

Sobol.p.fdr.irls <- p.adjust(Sobol.p,"fdr")
Sobol.p.hochberg.irls <- p.adjust(Sobol.p,"hochberg")
Sobol.p.hommel.irls <- p.adjust(Sobol.p,"hommel")

> Sobol.p.hommel.irls <- matrix(p.adjust(Sobol.p,"hommel"),nr=1)
> colnames(Sobol.p.hommel.irls) <- gene.names
> Sobol.p.hommel.irls[,order(Sobol.p.hommel.irls, decreasing=F)]
AHR     AHR.1      ARNT    ARNT.1     CEBPD     CEBPG   CEBPG.1      ESR1    ESR1.2 
0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
FOXA1     FOXA2   HNF4A.2   HNF4A.3     NCOA1     NCOR1     NCOR2   NCOR2.1   NCOR2.2 
0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
NFE2L2     NR0B2     NR1D2   NR1D2.1     NR1H2   NR1H2.1     NR1H3     NR1I2   NR1I2.1 
0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
NR1I3   NR1I3.1     NR5A2    PGRMC1   PPARA.1     PPARD   PPARD.1   PPARD.2    RXRA.1 
0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
THRA      THRB      USF1       VDR     VDR.1     VDR.2      RXRA     NR0B1     NR2F2 
0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0072000 0.0204000 0.0512000 
NR2F2.1      RXRB   ONECUT1    RXRB.1   PPARA.5     CEBPB   FOXA1.1     FOXA3   PPARG.2 
0.0512000 0.0640000 0.3557077 0.3557077 0.3670000 0.3978000 0.3978000 0.4410000 0.5138000 
CEBPA     NR3C1   PPARG.1     HNF4A      AHRR   NR2F1.1     NR2F1     NCOA2   PPARA.2 
0.5278000 0.6163636 0.6215000 0.6356250 0.6457143 0.6457143 0.6732000 0.6780000 0.6780000 
HNF4G   NCOA3.1       DBP    ESR1.1   HNF4A.1   HNF4G.1     NCOA3     NR1H4     PPARA 
0.7062500 0.7062500 0.7338000 0.7338000 0.7338000 0.7338000 0.7338000 0.7338000 0.7338000 
PPARA.3   PPARA.4     PPARG      RXRG    RXRG.1       YY1 
0.7338000 0.7338000 0.7338000 0.7338000 0.7338000 0.7338000 
> length(which(Sobol.p.hommel.irls<0.05))
[1] 44


> Sobol <- matrix(fit.glm.gs.1[[6]],nr=1)
> colnames(Sobol) <- gene.names
> Sobol[,order(Sobol, decreasing=T)]
ESR1       ESR1.2       PGRMC1      NR1I3.1        NR1I3      NR1I2.1      NCOR2.2 
0.2186828684 0.2095962578 0.1790353067 0.1770680375 0.1722913065 0.1270753990 0.1270075726 
NR1I2       NFE2L2        FOXA2      HNF4A.2      PPARA.1      HNF4A.3        NCOR2 
0.1256777320 0.1252722587 0.1241786123 0.1059270793 0.1024012324 0.1001404843 0.0936928112 
THRB        NCOR1        CEBPG         ARNT         USF1        AHR.1        PPARD 
0.0697162466 0.0687030855 0.0663473583 0.0653010368 0.0641301779 0.0607185795 0.0574439865 
NR1H2          AHR       RXRA.1        NR1H3        NR5A2      PPARD.2      CEBPG.1 
0.0572046245 0.0566589096 0.0565223006 0.0535760587 0.0530384671 0.0441075186 0.0418390915 
CEBPD      NR1D2.1        NR0B2        FOXA1      PPARD.1       ARNT.1         THRA 
0.0407870222 0.0406842269 0.0389767959 0.0327755756 0.0325325862 0.0287878837 0.0257725703 
VDR.2      NR1H2.1          VDR        VDR.1      NCOR2.1        NCOA1        NR1D2 
0.0256724293 0.0252655979 0.0214634644 0.0196396742 0.0195187584 0.0176999493 0.0171598489 
RXRA        NR0B1      NR2F2.1        NR2F2         RXRB      ONECUT1       RXRB.1 
0.0153217717 0.0139055039 0.0109051029 0.0104827985 0.0098048473 0.0059915216 0.0059822316 
PPARA.5      FOXA1.1        CEBPB        FOXA3      PPARG.2        CEBPA        NR3C1 
0.0058766421 0.0054161452 0.0053824193 0.0049505500 0.0044584825 0.0043396184 0.0038989111 
PPARG.1        HNF4A      NR2F1.1         AHRR        NR2F1        NCOA2      PPARA.2 
0.0036977779 0.0032998085 0.0030697855 0.0029923142 0.0026094573 0.0025836219 0.0025223156 
NCOA3.1        HNF4G       ESR1.1          DBP      HNF4A.1        NR1H4      PPARA.3 
0.0022247522 0.0021563032 0.0016650235 0.0016315703 0.0016074651 0.0015204903 0.0015077354 
NCOA3      HNF4G.1      PPARA.4        PPARG        PPARA         RXRG       RXRG.1 
0.0011133556 0.0008100389 0.0007243664 0.0006711795 0.0004986936 0.0003503314 0.0003339842 
YY1 
0.0001171454 
> Sobol.p.hommel.irls[,order(Sobol.p.hommel.irls, decreasing=F)]
AHR     AHR.1      ARNT    ARNT.1     CEBPD     CEBPG   CEBPG.1      ESR1    ESR1.2 
0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
FOXA1     FOXA2   HNF4A.2   HNF4A.3     NCOA1     NCOR1     NCOR2   NCOR2.1   NCOR2.2 
0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
NFE2L2     NR0B2     NR1D2   NR1D2.1     NR1H2   NR1H2.1     NR1H3     NR1I2   NR1I2.1 
0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
NR1I3   NR1I3.1     NR5A2    PGRMC1   PPARA.1     PPARD   PPARD.1   PPARD.2    RXRA.1 
0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
THRA      THRB      USF1       VDR     VDR.1     VDR.2      RXRA     NR0B1     NR2F2 
0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0072000 0.0204000 0.0512000 
NR2F2.1      RXRB   ONECUT1    RXRB.1   PPARA.5     CEBPB   FOXA1.1     FOXA3   PPARG.2 
0.0512000 0.0640000 0.3557077 0.3557077 0.3670000 0.3978000 0.3978000 0.4410000 0.5138000 
CEBPA     NR3C1   PPARG.1     HNF4A      AHRR   NR2F1.1     NR2F1     NCOA2   PPARA.2 
0.5278000 0.6163636 0.6215000 0.6356250 0.6457143 0.6457143 0.6732000 0.6780000 0.6780000 
HNF4G   NCOA3.1       DBP    ESR1.1   HNF4A.1   HNF4G.1     NCOA3     NR1H4     PPARA 
0.7062500 0.7062500 0.7338000 0.7338000 0.7338000 0.7338000 0.7338000 0.7338000 0.7338000 
PPARA.3   PPARA.4     PPARG      RXRG    RXRG.1       YY1 
0.7338000 0.7338000 0.7338000 0.7338000 0.7338000 0.7338000 
> length(which(Sobol.p.hommel.irls<0.05))
[1] 44
> 
  
  
  ##------------- univariate linear regression ----------
> cor.p.fdr[,order(cor.p.fdr, decreasing=F)]
ESR1       ESR1.2       PGRMC1      NR1I3.1        NR1I3       NFE2L2      NCOR2.2 
2.093570e-63 2.302903e-61 2.216627e-48 3.843621e-47 3.968614e-47 3.553676e-34 3.025474e-32 
NR1I2.1        NR1I2        FOXA2      HNF4A.2      HNF4A.3      PPARA.1        NCOR2 
3.025474e-32 6.960304e-32 1.102170e-31 3.305313e-26 8.259728e-25 8.800963e-25 6.423584e-24 
NCOR1         USF1        NR1H3         ARNT        AHR.1         THRB        NR1H2 
2.354888e-18 1.083937e-17 1.698608e-15 6.119432e-15 9.767796e-15 9.767796e-15 1.625866e-14 
CEBPG       RXRA.1          AHR        NR5A2        PPARD      PPARD.2        NR0B2 
2.330884e-14 2.812705e-14 8.592025e-14 1.358710e-13 2.283816e-13 2.241094e-10 1.252932e-09 
CEBPD      CEBPG.1        FOXA1      NR1D2.1       ARNT.1      PPARD.1      NR1H2.1 
1.532065e-09 1.846863e-09 1.874560e-08 1.903248e-08 2.765289e-08 5.972647e-08 3.912570e-07 
VDR.2         THRA      NCOR2.1          VDR        NR1D2         RXRA        NCOA1 
8.767821e-07 2.050905e-06 2.502879e-06 4.474560e-06 3.661493e-05 6.376162e-05 1.235173e-04 
VDR.1        NR0B1      NR2F2.1        NR2F2         RXRB      PPARA.5        FOXA3 
5.719362e-04 7.508924e-04 1.739402e-03 2.141086e-03 5.348846e-03 2.180723e-02 2.682554e-02 
DBP        CEBPB        HNF4A      FOXA1.1        NR3C1        CEBPA      NR2F1.1 
3.305405e-02 3.742513e-02 4.070217e-02 4.723801e-02 6.470068e-02 1.046799e-01 1.070875e-01 
AHRR      PPARG.2      ONECUT1      NCOA3.1       RXRB.1        NCOA2      HNF4A.1 
1.080745e-01 1.113004e-01 1.253026e-01 1.269909e-01 1.303449e-01 1.326055e-01 1.342122e-01 
NR2F1      PPARG.1        HNF4G      PPARA.2        NCOA3        NR1H4      PPARA.3 
1.482715e-01 1.661189e-01 2.022964e-01 2.293953e-01 2.476173e-01 2.586312e-01 2.586312e-01 
ESR1.1       RXRG.1      HNF4G.1        PPARA      PPARA.4          YY1        PPARG 
2.806175e-01 3.682209e-01 3.811925e-01 4.655869e-01 4.655869e-01 4.954409e-01 6.327069e-01 
RXRG 
8.571622e-01 

> length(which(p.fdr<0.05))
[1] 53




##------------- Multivariate linear regression --------

> lapply(fit.glm.gs.1, function(x){gene.names[order(x,decreasing=F)]})
$p.irls
[1] "NFE2L2"  "NR1D2.1" "NCOR2.2" "ESR1.1"  "NCOR2"   "HNF4A.1"
[7] "PPARA.1" "NR1H3"   "CEBPG.1" "NR1H2.1" "ESR1"    "ARNT"   
[13] "NCOA3.1" "NR1H2"   "PPARA.5" "HNF4G.1" "VDR.2"   "CEBPA"  
[19] "ONECUT1" "NCOR2.1" "CEBPG"   "THRB"    "RXRA"    "ARNT.1" 
[25] "FOXA2"   "RXRG"    "NR0B1"   "PPARA.3" "FOXA1"   "DBP"    
[31] "HNF4G"   "NR1I2.1" "HNF4A"   "NR1I3.1" "NR1I3"   "AHR"    
[37] "NR3C1"   "RXRB.1"  "PPARA"   "FOXA3"   "PPARG"   "NR0B2"  
[43] "THRA"    "NCOA1"   "NR1H4"   "RXRB"    "ESR1.2"  "PPARG.2"
[49] "NCOR1"   "NR1I2"   "CEBPB"   "VDR"     "VDR.1"   "YY1"    
[55] "CEBPD"   "AHR.1"   "RXRA.1"  "PPARA.4" "HNF4A.2" "PGRMC1" 
[61] "NCOA2"   "FOXA1.1" "PPARG.1" "NR2F1"   "NR2F1.1" "NR2F2.1"
[67] "NCOA3"   "USF1"    "PPARD.2" "NR1D2"   "RXRG.1"  "PPARA.2"
[73] "PPARD"   "NR2F2"   "PPARD.1" "AHRR"    "NR5A2"   "HNF4A.3"

$p.fdr.irls
[1] "NFE2L2"  "NR1D2.1" "NCOR2.2" "ESR1.1"  "NCOR2"   "CEBPG.1"
[7] "HNF4A.1" "NR1H3"   "PPARA.1" "NR1H2.1" "ARNT"    "CEBPA"  
[13] "CEBPG"   "ESR1"    "HNF4G.1" "NCOA3.1" "NCOR2.1" "NR1H2"  
[19] "ONECUT1" "PPARA.5" "VDR.2"   "RXRA"    "THRB"    "ARNT.1" 
[25] "FOXA2"   "NR0B1"   "PPARA.3" "RXRG"    "FOXA1"   "DBP"    
[31] "HNF4G"   "HNF4A"   "NR1I2.1" "NR1I3.1" "NR1I3"   "AHR"    
[37] "AHR.1"   "CEBPB"   "CEBPD"   "ESR1.2"  "FOXA3"   "NCOA1"  
[43] "NCOR1"   "NR0B2"   "NR1H4"   "NR1I2"   "NR3C1"   "PPARA"  
[49] "PPARG"   "PPARG.2" "RXRA.1"  "RXRB"    "RXRB.1"  "THRA"   
[55] "VDR"     "VDR.1"   "YY1"     "PPARA.4" "HNF4A.2" "FOXA1.1"
[61] "NCOA2"   "PGRMC1"  "PPARG.1" "NR2F1"   "NR2F1.1" "NR2F2.1"
[67] "NCOA3"   "USF1"    "NR1D2"   "PPARD.2" "RXRG.1"  "AHRR"   
[73] "NR2F2"   "PPARA.2" "PPARD"   "PPARD.1" "NR5A2"   "HNF4A.3"

$p.hochberg.irls
[1] "NFE2L2"  "NR1D2.1" "NCOR2.2" "ESR1.1"  "NCOR2"   "AHR"    
[7] "AHR.1"   "AHRR"    "ARNT"    "ARNT.1"  "CEBPA"   "CEBPB"  
[13] "CEBPD"   "CEBPG"   "CEBPG.1" "DBP"     "ESR1"    "ESR1.2" 
[19] "FOXA1"   "FOXA1.1" "FOXA2"   "FOXA3"   "HNF4A"   "HNF4A.1"
[25] "HNF4A.2" "HNF4A.3" "HNF4G"   "HNF4G.1" "NCOA1"   "NCOA2"  
[31] "NCOA3"   "NCOA3.1" "NCOR1"   "NCOR2.1" "NR0B1"   "NR0B2"  
[37] "NR1D2"   "NR1H2"   "NR1H2.1" "NR1H3"   "NR1H4"   "NR1I2"  
[43] "NR1I2.1" "NR1I3"   "NR1I3.1" "NR2F1"   "NR2F1.1" "NR2F2"  
[49] "NR2F2.1" "NR3C1"   "NR5A2"   "ONECUT1" "PGRMC1"  "PPARA"  
[55] "PPARA.1" "PPARA.2" "PPARA.3" "PPARA.4" "PPARA.5" "PPARD"  
[61] "PPARD.1" "PPARD.2" "PPARG"   "PPARG.1" "PPARG.2" "RXRA"   
[67] "RXRA.1"  "RXRB"    "RXRB.1"  "RXRG"    "RXRG.1"  "THRA"   
[73] "THRB"    "USF1"    "VDR"     "VDR.1"   "VDR.2"   "YY1"    

$p.hommel.irls
[1] "NFE2L2"  "NR1D2.1" "NCOR2.2" "ESR1.1"  "NCOR2"   "HNF4A.1"
[7] "NR1H3"   "PPARA.1" "CEBPG.1" "NR1H2.1" "AHR"     "AHR.1"  
[13] "AHRR"    "ARNT"    "ARNT.1"  "CEBPA"   "CEBPB"   "CEBPD"  
[19] "CEBPG"   "DBP"     "ESR1"    "ESR1.2"  "FOXA1"   "FOXA1.1"
[25] "FOXA2"   "FOXA3"   "HNF4A"   "HNF4A.2" "HNF4A.3" "HNF4G"  
[31] "HNF4G.1" "NCOA1"   "NCOA2"   "NCOA3"   "NCOA3.1" "NCOR1"  
[37] "NCOR2.1" "NR0B1"   "NR0B2"   "NR1D2"   "NR1H2"   "NR1H4"  
[43] "NR1I2"   "NR1I2.1" "NR1I3"   "NR1I3.1" "NR2F1"   "NR2F1.1"
[49] "NR2F2"   "NR2F2.1" "NR3C1"   "NR5A2"   "ONECUT1" "PGRMC1" 
[55] "PPARA"   "PPARA.2" "PPARA.3" "PPARA.4" "PPARA.5" "PPARD"  
[61] "PPARD.1" "PPARD.2" "PPARG"   "PPARG.1" "PPARG.2" "RXRA"   
[67] "RXRA.1"  "RXRB"    "RXRB.1"  "RXRG"    "RXRG.1"  "THRA"   
[73] "THRB"    "USF1"    "VDR"     "VDR.1"   "VDR.2"   "YY1"    

$p.multT.irls
[1] "NFE2L2"  "NR1D2.1" "NCOR2.2" "ESR1.1"  "NCOR2"   "HNF4A.1"
[7] "PPARA.1" "NR1H3"   "CEBPG.1" "NR1H2.1" "ESR1"    "ARNT"   
[13] "NCOA3.1" "NR1H2"   "PPARA.5" "HNF4G.1" "VDR.2"   "CEBPA"  
[19] "ONECUT1" "NCOR2.1" "CEBPG"   "THRB"    "RXRA"    "ARNT.1" 
[25] "FOXA2"   "RXRG"    "NR0B1"   "PPARA.3" "FOXA1"   "DBP"    
[31] "HNF4G"   "NR1I2.1" "HNF4A"   "NR1I3.1" "NR1I3"   "AHR"    
[37] "NR3C1"   "RXRB.1"  "PPARA"   "FOXA3"   "PPARG"   "NR0B2"  
[43] "THRA"    "NCOA1"   "NR1H4"   "AHR.1"   "AHRR"    "CEBPB"  
[49] "CEBPD"   "ESR1.2"  "FOXA1.1" "HNF4A.2" "HNF4A.3" "NCOA2"  
[55] "NCOA3"   "NCOR1"   "NR1D2"   "NR1I2"   "NR2F1"   "NR2F1.1"
[61] "NR2F2"   "NR2F2.1" "NR5A2"   "PGRMC1"  "PPARA.2" "PPARA.4"
[67] "PPARD"   "PPARD.1" "PPARD.2" "PPARG.1" "PPARG.2" "RXRA.1" 
[73] "RXRB"    "RXRG.1"  "USF1"    "VDR"     "VDR.1"   "YY1"    

$gs.irls
[1] "YY1"     "RXRG.1"  "RXRG"    "PPARA"   "PPARG"   "PPARA.4"
[7] "HNF4G.1" "NCOA3"   "PPARA.3" "NR1H4"   "HNF4A.1" "DBP"    
[13] "ESR1.1"  "HNF4G"   "NCOA3.1" "PPARA.2" "NCOA2"   "NR2F1"  
[19] "AHRR"    "NR2F1.1" "HNF4A"   "PPARG.1" "NR3C1"   "CEBPA"  
[25] "PPARG.2" "FOXA3"   "CEBPB"   "FOXA1.1" "PPARA.5" "RXRB.1" 
[31] "ONECUT1" "RXRB"    "NR2F2"   "NR2F2.1" "NR0B1"   "RXRA"   
[37] "NR1D2"   "NCOA1"   "NCOR2.1" "VDR.1"   "VDR"     "NR1H2.1"
[43] "VDR.2"   "THRA"    "ARNT.1"  "PPARD.1" "FOXA1"   "NR0B2"  
[49] "NR1D2.1" "CEBPD"   "CEBPG.1" "PPARD.2" "NR5A2"   "NR1H3"  
[55] "RXRA.1"  "AHR"     "NR1H2"   "PPARD"   "AHR.1"   "USF1"   
[61] "ARNT"    "CEBPG"   "NCOR1"   "THRB"    "NCOR2"   "HNF4A.3"
[67] "PPARA.1" "HNF4A.2" "FOXA2"   "NFE2L2"  "NR1I2"   "NCOR2.2"
[73] "NR1I2.1" "NR1I3"   "NR1I3.1" "PGRMC1"  "ESR1.2"  "ESR1"   

> 
  
  
  ##----------

> gene.cor <- cor(input.genes.m, use="complete.obs")
> gene.cor[1:5, 1:5]
AHR      AHR.1        AHRR        ARNT     ARNT.1
AHR     1.0000000  0.9973423 -0.26824279  0.12603165 0.19934134
AHR.1   0.9973423  1.0000000 -0.27381404  0.13340858 0.19900164
AHRR   -0.2682428 -0.2738140  1.00000000 -0.06755174 0.09522031
ARNT    0.1260316  0.1334086 -0.06755174  1.00000000 0.31393083
ARNT.1  0.1993413  0.1990016  0.09522031  0.31393083 1.00000000

> summary(gene.cor[lower.tri(gene.cor)])
Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
-0.752000 -0.153500  0.009497  0.024450  0.192400  0.997500

dim( which(gene.cor>0.2 & gene.cor<0.5, arr.ind=TRUE) )
> dim( which(gene.cor>0.2 & gene.cor<0.5, arr.ind=TRUE) )
[1] 1224    2

y <- CYPs.sorted.m[,30]

## correlation test:

cor.test <- function(y, x){
  fit.1 <- glm(y~x)
  glm.coef.1 <- fit.1$coef
  glm.p.1 <- summary(fit.1)$coef[,4];
  return(glm.p.1[2])
}

cor.test.res <- apply(input, 2, function(x){cor.test(y,x)})

p <- cor.test.res
p.fdr <- p.adjust(cor.test.res,"fdr")
p.hochberg <- p.adjust(cor.test.res,"hochberg")
p.hommel <- p.adjust(cor.test.res,"hommel")

dev.new();par(mfrow=c(2,2));apply(rbind(p, p.fdr, p.hochberg, p.hommel), 1, function(x){plot(x);abline(h=0.05,col="green");});

> gene.names <- names(input.genes)

> gene.names[order(p, decreasing=F)]
[1] "ESR1"    "ESR1.2"  "PGRMC1"  "NR1I3.1" "NR1I3"   "NFE2L2"  "NCOR2.2" "NR1I2.1"
[9] "NR1I2"   "FOXA2"   "HNF4A.2" "HNF4A.3" "PPARA.1" "NCOR2"   "NCOR1"   "USF1"   
[17] "NR1H3"   "ARNT"    "THRB"    "AHR.1"   "NR1H2"   "CEBPG"   "RXRA.1"  "AHR"    
[25] "NR5A2"   "PPARD"   "PPARD.2" "NR0B2"   "CEBPD"   "CEBPG.1" "FOXA1"   "NR1D2.1"
[33] "ARNT.1"  "PPARD.1" "NR1H2.1" "VDR.2"   "THRA"    "NCOR2.1" "VDR"     "NR1D2"  
[41] "RXRA"    "NCOA1"   "VDR.1"   "NR0B1"   "NR2F2.1" "NR2F2"   "RXRB"    "PPARA.5"
[49] "FOXA3"   "DBP"     "CEBPB"   "HNF4A"   "FOXA1.1" "NR3C1"   "CEBPA"   "NR2F1.1"
[57] "AHRR"    "PPARG.2" "ONECUT1" "NCOA3.1" "RXRB.1"  "NCOA2"   "HNF4A.1" "NR2F1"  
[65] "PPARG.1" "HNF4G"   "PPARA.2" "NCOA3"   "NR1H4"   "PPARA.3" "ESR1.1"  "RXRG.1" 
[73] "HNF4G.1" "PPARA"   "PPARA.4" "YY1"     "PPARG"   "RXRG"   
> length(which(p<0.05))
[1] 54
> gene.names[order(p.fdr, decreasing=F)]
[1] "ESR1"    "ESR1.2"  "PGRMC1"  "NR1I3.1" "NR1I3"   "NFE2L2"  "NCOR2.2" "NR1I2.1"
[9] "NR1I2"   "FOXA2"   "HNF4A.2" "HNF4A.3" "PPARA.1" "NCOR2"   "NCOR1"   "USF1"   
[17] "NR1H3"   "ARNT"    "AHR.1"   "THRB"    "NR1H2"   "CEBPG"   "RXRA.1"  "AHR"    
[25] "NR5A2"   "PPARD"   "PPARD.2" "NR0B2"   "CEBPD"   "CEBPG.1" "FOXA1"   "NR1D2.1"
[33] "ARNT.1"  "PPARD.1" "NR1H2.1" "VDR.2"   "THRA"    "NCOR2.1" "VDR"     "NR1D2"  
[41] "RXRA"    "NCOA1"   "VDR.1"   "NR0B1"   "NR2F2.1" "NR2F2"   "RXRB"    "PPARA.5"
[49] "FOXA3"   "DBP"     "CEBPB"   "HNF4A"   "FOXA1.1" "NR3C1"   "CEBPA"   "NR2F1.1"
[57] "AHRR"    "PPARG.2" "ONECUT1" "NCOA3.1" "RXRB.1"  "NCOA2"   "HNF4A.1" "NR2F1"  
[65] "PPARG.1" "HNF4G"   "PPARA.2" "NCOA3"   "NR1H4"   "PPARA.3" "ESR1.1"  "RXRG.1" 
[73] "HNF4G.1" "PPARA"   "PPARA.4" "YY1"     "PPARG"   "RXRG"   
> length(which(p.fdr<0.05))
[1] 53

CYP3A4 <- y
dev.new();par(mfrow=c(3,3));apply(input.genes[,order(p, decreasing=F)[1:9]], 2, function(x){plot(x, CYP3A4)});

> dev.new();par(mfrow=c(3,3));apply(input.genes[,order(p, decreasing=F)[1:9]], 2, function(x){plot(x, CYP3A4)});
NULL

"ESR1"    "ESR1.2"  "PGRMC1"  "NR1I3.1" "NR1I3"   "NFE2L2"  "NCOR2.2" "NR1I2.1" "NR1I2"

> dev.new();par(mfrow=c(3,3));apply(input.genes[,order(p, decreasing=F)[45:53]], 2, function(x){plot(x, CYP3A4)});
NULL
> 
  "NR2F2.1" "NR2F2"   "RXRB"    "PPARA.5" "FOXA3"   "DBP"     "CEBPB"   "HNF4A"   "FOXA1.1"

> p.fdr[order(p.fdr, decreasing=F)[1:9]]
ESR1       ESR1.2       PGRMC1      NR1I3.1        NR1I3 
2.093570e-63 2.302903e-61 2.216627e-48 3.843621e-47 3.968614e-47 
NFE2L2      NCOR2.2      NR1I2.1        NR1I2 
3.553676e-34 3.025474e-32 3.025474e-32 6.960304e-32 
> p.fdr[order(p.fdr, decreasing=F)[45:53]]
NR2F2.1       NR2F2        RXRB     PPARA.5       FOXA3 
0.001739402 0.002141086 0.005348846 0.021807231 0.026825540 
DBP       CEBPB       HNF4A     FOXA1.1 
0.033054046 0.037425127 0.040702172 0.047238006



##-------------------------------- all input.genes summary statistics : --------
> apply(input.genes, 2, summary)
$AHR
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.87670 -0.24990 -0.04960 -0.07159  0.10920  0.55890       30 

$AHR.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.91470 -0.26090 -0.06310 -0.07854  0.10030  0.55020       12 

$AHRR
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.78060 -0.11030 -0.02645  0.01099  0.08807  2.00000        1 

$ARNT
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.317600 -0.084450  0.001400  0.003516  0.080150  0.335300         8 

$ARNT.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.51540 -0.16340 -0.07615 -0.07214  0.01638  0.39920        3 

$CEBPA
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.536200 -0.193000 -0.094050 -0.087650  0.002075  0.391500        13 

$CEBPB
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.32950 -0.06693  0.04220  0.05988  0.16430  1.14800        3 

$CEBPD
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.47210 -0.12470  0.01380  0.02394  0.17520  0.66000       10 

$CEBPG
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.53900 -0.15120 -0.04980 -0.03085  0.08350  0.42330        7 

$CEBPG.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.46630 -0.12980 -0.01480 -0.01131  0.09040  0.46460        6 

$DBP
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.252200 -0.076620 -0.038250 -0.032860  0.001425  2.000000         3 

$ESR1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-1.18700 -0.41600 -0.06875 -0.13830  0.16790  0.53670        1 

$ESR1.1
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
-0.4541  0.0465  0.1481  0.1724  0.2695  0.8518 

$ESR1.2
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
-2.0000 -0.4798 -0.0994 -0.1861  0.1582  0.5132       8 

$FOXA1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.70880 -0.15250 -0.04370 -0.04693  0.06440  0.39490       11 

$FOXA1.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.26470 -0.06000  0.01650  0.02795  0.09690  2.00000        2 

$FOXA2
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.44720 -0.13780 -0.00815 -0.02716  0.08592  0.36210        1 

$FOXA3
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.26410 -0.04665 -0.01760 -0.01516  0.01200  0.20950       16 

$HNF4A
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.22670 -0.02635  0.04010  0.05001  0.10980  0.98610        4 

$HNF4A.1
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
-0.5522 -0.1634 -0.1040 -0.1049 -0.0454  0.4451       4 

$HNF4A.2
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.71030 -0.18160 -0.03595 -0.05413  0.08855  0.36440        3 

$HNF4A.3
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.68740 -0.19270 -0.03610 -0.05182  0.09430  0.36560       32 

$HNF4G
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-2.00000 -0.10510 -0.04575 -0.09214  0.00335  0.22700        3 

$HNF4G.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-1.23400 -0.08885  0.06545  0.04040  0.18200  0.58000        9 

$NCOA1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.21660 -0.05210  0.00830  0.01111  0.06360  0.45260        2 

$NCOA2
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-2.000000 -0.140900 -0.074600 -0.088110 -0.007025  2.000000         3 

$NCOA3
Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
-0.313800 -0.052000 -0.002900  0.006585  0.054200  0.353900 

$NCOA3.1
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.369100 -0.049400  0.001400 -0.000031  0.063100  0.415800        26 

$NCOR1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.45610 -0.15880 -0.06280 -0.02993  0.08110  2.00000        6 

$NCOR2
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.457100 -0.109400 -0.013600 -0.001085  0.099800  0.390900         2 

$NCOR2.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.10580  0.00225  0.04110  0.04539  0.07570  0.39450        5 

$NCOR2.2
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.39570 -0.08420  0.01875  0.03320  0.13850  0.61610        1 

$NFE2L2
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.56020 -0.16690 -0.02100 -0.03212  0.10740  0.42330        5 

$NR0B1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-2.00000 -0.05425 -0.00970 -0.01067  0.04945  2.00000        4 

$NR0B2
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
-0.8689 -0.2818 -0.0607 -0.0771  0.1488  0.9219      25 

$NR1D2
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-2.00000 -0.18430 -0.03320 -0.04904  0.09840  0.71950        2 

$NR1D2.1
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.305800 -0.048100 -0.010200 -0.007857  0.028200  0.245300         2 

$NR1H2
Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
-0.16150 -0.05920 -0.02440 -0.01068  0.01600  0.28820 

$NR1H2.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.13730 -0.02650  0.01540  0.02296  0.06100  0.35700       14 

$NR1H3
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.56310 -0.13850 -0.04860 -0.04435  0.04095  0.40480        4 

$NR1H4
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-1.05200 -0.10890 -0.01520 -0.02797  0.05740  0.27280       58 

$NR1I2
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
-1.0640 -0.2941 -0.0305 -0.1067  0.1177  0.4092       4 

$NR1I2.1
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
-1.2450 -0.3146 -0.0266 -0.1118  0.1360  0.4762       3 

$NR1I3
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-2.00000 -0.28130 -0.08520 -0.11210  0.06013  0.66100        3 

$NR1I3.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-1.29300 -0.26740 -0.06800 -0.09552  0.07730  0.65000       10 

$NR2F1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.49430 -0.09995  0.00960  0.01875  0.12700  0.52950       36 

$NR2F1.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.47280 -0.09390  0.01280  0.02306  0.12900  0.52910        4 

$NR2F2
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.66170 -0.14550 -0.01855 -0.03459  0.08300  0.43250        5 

$NR2F2.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.66940 -0.14620 -0.01910 -0.03408  0.08260  0.43640       10 

$NR3C1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.55670 -0.11460 -0.02530 -0.02408  0.06682  0.29280        1 

$NR5A2
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.60750 -0.12920 -0.02625 -0.02438  0.07095  0.50850       31 

$ONECUT1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
-0.50800 -0.06455 -0.01500 -0.01003  0.04275  0.40480 

$PGRMC1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.59240 -0.16540 -0.01610 -0.03501  0.11520  0.40450        8 

$PPARA
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.363700 -0.081450  0.000750 -0.003853  0.081750  0.382600         1 

$PPARA.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.56470 -0.12570 -0.01785 -0.01280  0.09740  0.57050        1 

$PPARA.2
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.507500 -0.091000 -0.017100  0.001722  0.062800  0.980800         2 

$PPARA.3
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.66720 -0.15020 -0.02920 -0.03927  0.07965  0.51850        1 

$PPARA.4
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.80270 -0.13330 -0.00730 -0.04307  0.08235  0.42660       56 

$PPARA.5
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.61820 -0.11720 -0.02305 -0.03065  0.05790  0.51490       27 

$PPARD
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.474900 -0.154700 -0.053050 -0.002164  0.141100  0.912100         9 

$PPARD.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
-1.30700 -0.10420 -0.01320  0.02606  0.12860  0.83660 

$PPARD.2
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.347500 -0.131400 -0.033650  0.004966  0.104200  0.981500         9 

$PPARG
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-1.45900 -0.19830 -0.05580 -0.05456  0.09510  2.00000        2 

$PPARG.1
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.469900 -0.130900 -0.009000 -0.006644  0.103400  0.530400         4 

$PPARG.2
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.482400 -0.123800 -0.004300 -0.004518  0.109200  0.514700         2 

$RXRA
Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
-0.328500 -0.063200 -0.014300 -0.008356  0.045300  0.310500 

$RXRA.1
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.958000 -0.140300  0.040000  0.008409  0.183100  0.515400        10 

$RXRB
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.273400 -0.058300 -0.003450 -0.005627  0.047100  0.426600        23 

$RXRB.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.54510 -0.10630 -0.05780 -0.05653 -0.01100  0.36100        2 

$RXRG
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.58510 -0.12040 -0.06360 -0.06718 -0.01095  0.23510        4 

$RXRG.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-2.00000 -0.02835  0.02300  0.03713  0.07090  2.00000        4 

$THRA
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
-0.239600 -0.049250 -0.000400  0.002138  0.048780  0.262500         1 

$THRB
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.42530 -0.09000 -0.03375 -0.02945  0.02852  0.20980        1 

$USF1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-2.00000 -0.15910 -0.04020 -0.02054  0.11040  1.35000        2 

$VDR
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.85480 -0.25180 -0.06925 -0.06363  0.10000  1.42900        3 

$VDR.1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.28950 -0.06760 -0.02100 -0.02574  0.01760  0.75100        2 

$VDR.2
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.99500 -0.24710 -0.06130 -0.05981  0.11120  1.42100        2 

$YY1
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
-0.44870 -0.07450 -0.00100 -0.01244  0.05400  0.32970       10 

> 
  
  ##--------------------------------------- before 01 20 2016 ---------------------------------





data <- read.table("U:/Data/Danxin/CYP3A4_Microarray/CYP3A4_Microarray.csv", sep=",", header=T)
data[1:5,1:10]

data.1 <- read.table("U:/Data/Danxin/CYP3A4_Microarray/moreCYPs.csv", sep=",", header=T)
data.1[1:5,1:10]

> dim(data.2)
[1] 427  43
> dim(data)
[1] 427  82
> data[1:5,1:10]
X     AHR   AHR.1   AHRR    ARNT  ARNT.1   CEBPA  CEBPB  CEBPD   CEBPG
1 GSM242213 -0.6493 -0.6612 0.2826 -0.1619 -0.0474 -0.2233 0.4964 0.1738  0.0059
2 GSM242214      NA -0.3682 0.0336 -0.0940  0.0750 -0.0048 0.2513 0.1068 -0.2858
3 GSM242215 -0.2992 -0.3068 0.3418 -0.0732 -0.0996  0.1988 0.5854 0.2192  0.0251
4 GSM242216 -0.3256 -0.3363 0.2179 -0.0157  0.1440  0.0129 0.5688 0.0287  0.0670
5 GSM242217 -0.4974 -0.5120 0.1424 -0.0809  0.1111  0.1000 0.2296 0.2814 -0.1105
> data.2[1:5,1:10]
CYP4F11  CYP1A1  CYP3A7  CYP2U1 CYP51A1 CYP2D6  CYP2C9  CYP2B6 CYP2C19  CYP1B1
GSM242213  0.1840  0.2021 -0.6217 -0.0769 -0.0534 0.4206  0.1318  0.1976  0.3618  0.4384
GSM242214  0.0527 -0.2326 -1.4490  0.0717  0.2282 0.3601  0.0240 -0.2499  0.2895  0.2078
GSM242215  0.3110  0.3020  0.1535 -0.2346  0.3859 0.6635  0.2906 -0.1127  0.4624 -0.3988
GSM242216  0.0970 -0.0447 -0.9564  0.1200 -0.1218 0.5671 -0.1834 -0.8589  0.0688  0.1541
GSM242217  0.1087 -0.1591 -0.0728 -0.1518 -0.2899 0.4622  0.3110  0.3884  0.4838 -0.1837
> table(data[,1]==rownames(data.2))

TRUE 
427 
> data.with.more.CYPs <- cbind(data, data.2)
> data.with.more.CYPs[1:5,1:10]
X     AHR   AHR.1   AHRR    ARNT  ARNT.1   CEBPA  CEBPB  CEBPD   CEBPG
GSM242213 GSM242213 -0.6493 -0.6612 0.2826 -0.1619 -0.0474 -0.2233 0.4964 0.1738  0.0059
GSM242214 GSM242214      NA -0.3682 0.0336 -0.0940  0.0750 -0.0048 0.2513 0.1068 -0.2858
GSM242215 GSM242215 -0.2992 -0.3068 0.3418 -0.0732 -0.0996  0.1988 0.5854 0.2192  0.0251
GSM242216 GSM242216 -0.3256 -0.3363 0.2179 -0.0157  0.1440  0.0129 0.5688 0.0287  0.0670
GSM242217 GSM242217 -0.4974 -0.5120 0.1424 -0.0809  0.1111  0.1000 0.2296 0.2814 -0.1105
> dim(data.with.more.CYPs)
[1] 427 125

write.csv(data.with.more.CYPs, "U:/Data/Danxin/CYP3A4_Microarray/CYP3A4_Microarray_with_more_CYPs.csv")


save(CYPs.sorted, input.genes, file="U:/Data/Danxin/CYP3A4_Microarray/CYP3A4_Microarray_input_CYPs_separated.Rdata")

## --------------------- combining the previous microarray data with more CYPs data is done ------------

library(MASS)
##library(glmnet)

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

##------------------------------------------------------------------
##---------- select the best variable combination ------------------
##------------------------------------------------------------------
## cs: combination size 
## coef: the vector glm coef estimates, minus the intercept estimate
##covm: empirical estimate of the inputs covariance matrix

multiple.variable.selection.sub <- function(index, coef, covm){
  eta <- matrix(coef[index],nc=1) + solve(covm[index, index]) %*% covm[index,-index] %*% matrix(coef[-index],nc=1)
  return( t(eta) %*% covm[index, index] %*% eta )
}
## FYI: if solve() return error can try ginv() after loading MASS package.

multiple.variable.selection <- function(cs=5, coef, covm){
  index.m <- combn(c(1:length(coef)), cs)   ## index.m is a matrix with cs rows.
  mv.gs <- apply(index.m, 2, function(x){multiple.variable.selection.sub(x, coef,covm)})
  cc <- rbind(index.m, mv.gs)
  ccc <- cc[,order(mv.gs, decreasing=T)]
  return(list(ccc[-(cs+1),], ccc[(cs+1),]))
}


##---------------------- analysis code
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


save(CYPs.sorted, input.genes, CYPs.sorted.m, input.genes.m, scan.all.CYPs, scan.all.CYPs.against.noise, sig.est.gs, sig.est.gs.with.noise, file="U:/Data/Danxin/CYP3A4_Microarray/Sobol_selection_result1.Rdata")

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

save(CYPs.sorted, input.genes, CYPs.sorted.m, input.genes.m, scan.all.CYPs, scan.all.CYPs.against.noise, sig.est.gs, sig.est.gs.with.noise, scan.all.CYPs.against.noise.insig.genes, insig.est.gs.with.noise, file="U:/Data/Danxin/CYP3A4_Microarray/Sobol_selection_result2.Rdata")








##-------------- input normality test
input.normality.test <- unlist( apply(input.genes.m, 2, function(x){t <- shapiro.test(x); return(t$p.value)}) )
length(which(input.normality.test<=0.0001))
hist(input.genes[,3],breaks=30)








##--------------- output

> ls()
[1] "CYPs.sorted" "input.genes"
> class(CYPs.sorted)
[1] "data.frame"
> class(input.genes)
[1] "data.frame"
> dim(CYPs.sorted)
[1] 427  46
> dim(input.genes)
[1] 427  78


> CYPs.sorted.m <- as.matrix(CYPs.sorted)
> dim(CYPs.sorted.m)
[1] 427  46
> input.genes.m <- as.matrix(input.genes)
> dim(input.genes.m)
[1] 427  78
> CYPs.sorted.m[1:3,1:5]
CYP1A1 CYP1A1.1  CYP1B1 CYP1B1.1 CYP1B1.2
GSM242213  0.2021   0.2235  0.4384  -0.0338   0.4809
GSM242214 -0.2326  -0.9000  0.2078  -0.1180   0.2592
GSM242215  0.3020  -0.1735 -0.3988   0.0397  -0.3320
> CYPs.sorted[1:3,1:5]
CYP1A1 CYP1A1.1  CYP1B1 CYP1B1.1 CYP1B1.2
GSM242213  0.2021   0.2235  0.4384  -0.0338   0.4809
GSM242214 -0.2326  -0.9000  0.2078  -0.1180   0.2592
GSM242215  0.3020  -0.1735 -0.3988   0.0397  -0.3320

> input.genes[1:3,1:5]
AHR   AHR.1   AHRR    ARNT  ARNT.1
GSM242213 -0.6493 -0.6612 0.2826 -0.1619 -0.0474
GSM242214      NA -0.3682 0.0336 -0.0940  0.0750
GSM242215 -0.2992 -0.3068 0.3418 -0.0732 -0.0996
> input.genes.m[1:3,1:5]
AHR   AHR.1   AHRR    ARNT  ARNT.1
GSM242213 -0.6493 -0.6612 0.2826 -0.1619 -0.0474
GSM242214      NA -0.3682 0.0336 -0.0940  0.0750
GSM242215 -0.2992 -0.3068 0.3418 -0.0732 -0.0996
> 
  > ptm <- proc.time()
> scan.all.CYPs <- apply(CYPs.sorted.m, 2, function(y){scan.one.CYP(y,input.genes.m)})
> proc.time()-ptm
user  system elapsed 
3.96    0.07    4.02 
> 
  > names(scan.all.CYPs[[1]])
[1] "glm.coef.with.intercept" "glm.p.with.intercept"    "input.se.est"           
[4] "input.rho.est"           "est.gs"                  "est.gs.order"           
[7] "glm.p.order"            
> 
  > ##---------
> glm.p.c <- do.call(rbind, lapply(scan.all.CYPs, function(x){x[["glm.p.with.intercept"]]}))
> glm.p.c.fdr <- t(apply(glm.p.c, 1, function(x){p.adjust(x,"fdr")}))
> ## number of variables selected by controlling fdr at 0.05
  > table(apply(glm.p.c.fdr, 1, function(x){length(which(x<=0.05))}))

0  1  2  4 
38  6  1  1 
> 
  > ##----------- fdr adjusted p-value selection (fdr=0.05) ------------
> sig.glm.p.fdr <- apply(glm.p.c.fdr, 1, function(x){which(x<=0.05)})
> ##sig.glm.p.fdr.v <- unlist(sig.glm.p.fdr)
  > ##length(sig.glm.p.fdr.v)
  > sig.glm.p.fdr
$CYP1A1
named integer(0)

$CYP1A1.1
named integer(0)

$CYP1B1
named integer(0)

$CYP1B1.1
named integer(0)

$CYP1B1.2
named integer(0)

$CYP2A6
named integer(0)

$CYP2A6.1
named integer(0)

$CYP2B6
named integer(0)

$CYP2B6.1
named integer(0)

$CYP2C18
NR1H3 
41 

$CYP2C18.1
NR1H3 
41 

$CYP2C19
named integer(0)

$CYP2C19.1
named integer(0)

$CYP2C8
named integer(0)

$CYP2C8.1
named integer(0)

$CYP2C9
named integer(0)

$CYP2C9.1
named integer(0)

$CYP2D6
named integer(0)

$CYP2D6.1
named integer(0)

$CYP2E1
named integer(0)

$CYP2E1.1
named integer(0)

$CYP2F1
(Intercept)       HNF4G       NR0B1        USF1 
1          24          35          75 

$CYP2F1.1
named integer(0)

$CYP2J2
named integer(0)

$CYP2J2.1
named integer(0)

$CYP2R1
HNF4G 
24 

$CYP2R1.1
named integer(0)

$CYP2U1
named integer(0)

$CYP2U1.1
NCOR2.1    THRB 
32      74 

$CYP3A4
named integer(0)

$CYP3A4.1
named integer(0)

$CYP3A4.2
named integer(0)

$CYP3A43
named integer(0)

$CYP3A5
named integer(0)

$CYP3A5.1
named integer(0)

$CYP3A7
named integer(0)

$CYP3A7.1
named integer(0)

$CYP4F11
named integer(0)

$CYP4F11.1
NR1H2 
39 

$CYP4F12
named integer(0)

$CYP4F12.1
named integer(0)

$CYP51A1
named integer(0)

$CYP51A1.1
named integer(0)

$CYP51A1.2
named integer(0)

$CYP7A1
(Intercept) 
1 

$CYP7A1.1
(Intercept) 
1 

> 
  > glm.p.c.fdr.order <- t(apply(glm.p.c.fdr,1,function(x){order(x,decreasing=FALSE, na.last=TRUE)}))
> dim(glm.p.c.fdr.order)
[1] 46 79
> input.gene.names <- names(input.genes)
> names <- c("intercept",input.gene.names)
> apply(glm.p.c.fdr.order[,1:10],1,function(x){names[x]})
CYP1A1    CYP1A1.1  CYP1B1      CYP1B1.1    CYP1B1.2    CYP2A6      CYP2A6.1    CYP2B6     
[1,] "ONECUT1" "ONECUT1" "intercept" "NR0B1"     "intercept" "intercept" "intercept" "intercept"
[2,] "NR1I3.1" "RXRA"    "RXRA.1"    "intercept" "RXRA.1"    "AHR.1"     "AHR"       "NFE2L2"   
[3,] "AHRR"    "ARNT"    "AHR"       "PPARD.1"   "AHR"       "AHRR"      "AHR.1"     "AHR.1"    
[4,] "HNF4G.1" "CEBPG"   "AHR.1"     "RXRG"      "AHR.1"     "CEBPD"     "AHRR"      "NCOR2.1"  
[5,] "NCOR2"   "DBP"     "ARNT"      "VDR.1"     "ARNT"      "CEBPG"     "ARNT.1"    "AHRR"     
[6,] "NCOR2.2" "HNF4G.1" "ARNT.1"    "CEBPD"     "ARNT.1"    "CEBPG.1"   "CEBPB"     "AHR"      
[7,] "NR0B1"   "NCOA3.1" "CEBPB"     "FOXA1"     "CEBPB"     "ESR1"      "CEBPD"     "HNF4A.3"  
[8,] "NR1I2"   "NCOR2"   "CEBPG.1"   "FOXA1.1"   "CEBPG.1"   "ESR1.2"    "CEBPG"     "HNF4A.2"  
[9,] "NR1I2.1" "NR1H2"   "ESR1"      "NCOA1"     "DBP"       "FOXA1"     "CEBPG.1"   "NCOA1"    
[10,] "NR1I3"   "PPARA"   "ESR1.2"    "NCOR2.1"   "ESR1"      "HNF4A"     "DBP"       "DBP"      

CYP2B6.1    CYP2C18   CYP2C18.1 CYP2C19     CYP2C19.1   CYP2C8      CYP2C8.1    CYP2C9     
[1,] "NFE2L2"    "NR1H3"   "NR1H3"   "intercept" "intercept" "intercept" "intercept" "intercept"
[2,] "intercept" "CEBPA"   "NR1I3.1" "NR3C1"     "AHR"       "DBP"       "DBP"       "AHRR"     
[3,] "AHR.1"     "NR1I3.1" "CEBPA"   "CEBPA"     "AHR.1"     "HNF4A.2"   "FOXA1"     "CEBPA"    
[4,] "NCOR2.1"   "HNF4A.1" "HNF4A.1" "NCOA3.1"   "AHRR"      "HNF4A.3"   "NCOA3.1"   "CEBPD"    
[5,] "AHRR"      "NCOR2.1" "NR1H4"   "NCOR2"     "ARNT"      "NCOA1"     "NR1H2"     "HNF4A.2"  
[6,] "DBP"       "NR0B1"   "NCOR2.1" "NCOR2.2"   "ARNT.1"    "AHR.1"     "NR3C1"     "HNF4A.3"  
[7,] "NR2F1"     "NR1H4"   "NR0B1"   "NR1H2"     "CEBPA"     "NFE2L2"    "PGRMC1"    "NCOA3.1"  
[8,] "AHR"       "NCOR2"   "NCOR2"   "PGRMC1"    "CEBPB"     "NR1H2"     "PPARG"     "NCOR2"    
[9,] "HNF4A.2"   "DBP"     "NR2F1.1" "NCOA1"     "CEBPD"     "NR2F1"     "HNF4A.2"   "NCOR2.1"  
[10,] "HNF4A.3"   "FOXA3"   "FOXA3"   "PPARG.2"   "CEBPG"     "NR2F1.1"   "HNF4A.3"   "NCOR2.2"  

CYP2C9.1  CYP2D6    CYP2D6.1    CYP2E1    CYP2E1.1  CYP2F1      CYP2F1.1    CYP2J2    CYP2J2.1 
[1,] "NR3C1"   "ARNT"    "intercept" "CEBPG"   "CEBPG.1" "HNF4G"     "intercept" "NR1I3.1" "NCOR1"  
[2,] "AHRR"    "CEBPA"   "AHR"       "CEBPG.1" "FOXA1"   "NR0B1"     "AHR"       "NCOR1"   "NR1I3.1"
[3,] "CEBPD"   "ESR1"    "AHR.1"     "FOXA1"   "FOXA3"   "USF1"      "AHR.1"     "NR1I3"   "NR1I3"  
[4,] "NR1H2"   "ESR1.1"  "AHRR"      "FOXA1.1" "HNF4A.1" "intercept" "AHRR"      "PPARG.1" "CEBPD"  
[5,] "NR1I2"   "ESR1.2"  "ARNT"      "FOXA3"   "NCOR1"   "NCOR2.1"   "ARNT"      "PPARG.2" "PPARA.2"
[6,] "PGRMC1"  "FOXA1.1" "ARNT.1"    "HNF4A.1" "NR3C1"   "RXRB"      "ARNT.1"    "NCOA3"   "PPARG.1"
[7,] "PPARG.2" "FOXA2"   "CEBPA"     "HNF4G.1" "PPARA.5" "ONECUT1"   "CEBPA"     "PPARA.2" "PPARG.2"
[8,] "RXRA"    "HNF4G"   "CEBPB"     "NCOR1"   "PPARD"   "NR2F1"     "CEBPB"     "AHRR"    "ARNT"   
[9,] "CEBPA"   "HNF4G.1" "CEBPD"     "NCOR2.1" "PPARD.2" "NFE2L2"    "CEBPD"     "ARNT"    "NR1I2.1"
[10,] "RXRG.1"  "NCOA1"   "CEBPG"     "NR0B1"   "NR1I3.1" "NR1I3.1"   "CEBPG"     "CEBPD"   "AHRR"   

CYP2R1    CYP2R1.1  CYP2U1    CYP2U1.1  CYP3A4      CYP3A4.1    CYP3A4.2    CYP3A43    
[1,] "HNF4G"   "NR0B1"   "HNF4G"   "NCOR2.1" "intercept" "intercept" "intercept" "intercept"
[2,] "AHRR"    "VDR"     "THRB"    "THRB"    "NFE2L2"    "NFE2L2"    "RXRB.1"    "ESR1"     
[3,] "ARNT"    "VDR.2"   "FOXA3"   "HNF4G"   "NR0B2"     "RXRB.1"    "NFE2L2"    "ESR1.2"   
[4,] "CEBPA"   "CEBPA"   "NCOR1"   "YY1"     "NR1D2.1"   "AHRR"      "NR2F1"     "FOXA2"    
[5,] "CEBPD"   "HNF4A"   "RXRG"    "NR0B1"   "NR2F1"     "ARNT"      "NR2F1.1"   "HNF4G"    
[6,] "CEBPG.1" "HNF4G"   "ESR1.2"  "PPARG"   "NR2F1.1"   "DBP"       "DBP"       "NCOA1"    
[7,] "DBP"     "HNF4G.1" "NR1I3"   "AHR"     "PPARD.2"   "ESR1"      "NCOA3.1"   "NR1I3"    
[8,] "ESR1.1"  "NR1H3"   "CEBPG.1" "AHR.1"   "RXRB.1"    "ESR1.2"    "PPARD.2"   "PPARG"    
[9,] "FOXA1"   "ONECUT1" "ESR1"    "NR3C1"   "AHRR"      "FOXA2"     "FOXA2"     "VDR.2"    
[10,] "HNF4A"   "PPARA.3" "NR1H2"   "RXRB"    "ARNT"      "NCOA3.1"   "HNF4A.2"   "RXRA.1"   

CYP3A5    CYP3A5.1  CYP3A7      CYP3A7.1    CYP4F11   CYP4F11.1 CYP4F12   CYP4F12.1 CYP51A1    
[1,] "NR1D2.1" "NR1H2.1" "intercept" "intercept" "NR1H2"   "NR1H2"   "NR3C1"   "CEBPA"   "NR2F1"    
[2,] "NR1H2.1" "NCOA3"   "RXRB.1"    "AHRR"      "NR1I3.1" "RXRB.1"  "ONECUT1" "NR1I2.1" "NR2F1.1"  
[3,] "RXRA"    "RXRA"    "AHRR"      "ARNT"      "NR0B1"   "THRA"    "PPARA.4" "NR3C1"   "intercept"
[4,] "NCOA1"   "PPARD.1" "ARNT"      "ESR1"      "CEBPG.1" "CEBPD"   "THRB"    "NR1I2"   "ESR1"     
[5,] "PPARD.1" "NCOA1"   "DBP"       "ESR1.2"    "FOXA1"   "HNF4A.2" "NR1I2"   "NR2F2"   "ESR1.2"   
[6,] "FOXA3"   "NR1D2.1" "ESR1"      "FOXA1"     "FOXA3"   "HNF4A.3" "NR1I2.1" "NR2F2.1" "FOXA2"    
[7,] "NCOA3"   "USF1"    "ESR1.2"    "FOXA2"     "HNF4G.1" "RXRB"    "USF1"    "ONECUT1" "NCOA1"    
[8,] "NFE2L2"  "NCOR1"   "FOXA1"     "NFE2L2"    "NCOR1"   "NR1H4"   "NR1D2.1" "CEBPD"   "NCOR2.1"  
[9,] "ARNT"    "RXRG.1"  "FOXA2"     "NR0B2"     "NCOR2.1" "FOXA3"   "NR2F2.1" "DBP"     "NCOR2.2"  
[10,] "ESR1.1"  "PPARA.4" "HNF4A.2"   "RXRB.1"    "NR1I3"   "HNF4G"   "VDR.1"   "THRA"    "NR0B1"    

CYP51A1.1 CYP51A1.2 CYP7A1      CYP7A1.1   
[1,] "CEBPB"   "FOXA2"   "intercept" "intercept"
[2,] "CEBPG.1" "NCOA1"   "NR1D2"     "AHR"      
[3,] "NCOR1"   "NCOR2.1" "PPARA.1"   "NR1D2"    
[4,] "NR1D2"   "NR0B1"   "VDR.2"     "NR5A2"    
[5,] "PPARA.2" "NR1H2.1" "FOXA3"     "PPARA.1"  
[6,] "PPARA.5" "NR2F1"   "NCOR2.2"   "PPARA.5"  
[7,] "RXRG.1"  "NR2F1.1" "PPARA.5"   "VDR.2"    
[8,] "YY1"     "PPARD.1" "NCOR2"     "AHR.1"    
[9,] "CEBPG"   "PPARG.1" "VDR"       "FOXA3"    
[10,] "NR5A2"   "PPARG.2" "AHR"       "HNF4G"    
> 
  
  
  > sig.est.gs[1:5]
$CYP1A1
[1]  1  2  5  9 10 11 12 14 21 22 28 29 32 33 34 37 40 42 43 44 45 53 57 58 59 60 61 62 66 71 76

$CYP1A1.1
[1]  5  9 10 11 12 14 21 22 25 28 42 43 44 45 53 57 58 59 60 61 62 66 72 77

$CYP1B1
[1]  4  6 11 15 19 21 22 25 35 41 42 43 44 45 54 56 57 58 59 61 68 73 75 77 78

$CYP1B1.1
[1]  1  2  4  8  9 10 12 14 17 18 22 30 31 32 34 38 42 43 44 45 52 53 55 56 60 61 62 69 72 73 74 76

$CYP1B1.2
[1]  4  6 11 15 19 20 21 22 25 28 35 41 42 43 44 45 48 49 56 57 58 59 61 68 73 75 77 78

> unlist(lapply(sig.est.gs, function(x){length(x)}))
CYP1A1  CYP1A1.1    CYP1B1  CYP1B1.1  CYP1B1.2    CYP2A6  CYP2A6.1    CYP2B6  CYP2B6.1   CYP2C18 CYP2C18.1 
31        24        25        32        28        25        24        28        29        23        23 
CYP2C19 CYP2C19.1    CYP2C8  CYP2C8.1    CYP2C9  CYP2C9.1    CYP2D6  CYP2D6.1    CYP2E1  CYP2E1.1    CYP2F1 
25        29        29        27        25        27        32        29        26        28        10 
CYP2F1.1    CYP2J2  CYP2J2.1    CYP2R1  CYP2R1.1    CYP2U1  CYP2U1.1    CYP3A4  CYP3A4.1  CYP3A4.2   CYP3A43 
27        24        25        24         9        32        29        25        26        26        28 
CYP3A5  CYP3A5.1    CYP3A7  CYP3A7.1   CYP4F11 CYP4F11.1   CYP4F12 CYP4F12.1   CYP51A1 CYP51A1.1 CYP51A1.2 
25        23        27        23        29        29        25        27        33        19        33 
CYP7A1  CYP7A1.1 
29        28 
> 
  > dim(input.genes)
[1] 427  78

##---------------------------------
> ptm <- proc.time()
> scan.all.CYPs.against.noise <- apply(matrix(1:46,nr=1), 2, function(i){test.sig.against.noise(CYPs.sorted.m[,i],input.genes.m,sig.est.gs[[i]])})
> proc.time()-ptm
user  system elapsed 
3.01    0.15    3.18 
>  
  > names(scan.all.CYPs.against.noise[[1]])
NULL

> scan.all.CYPs.against.noise[[1]][[2]]
[1] 31
> names(scan.all.CYPs.against.noise[[1]][[1]])
[1] "glm.coef.with.intercept" "glm.p.with.intercept"    "input.se.est"            "input.rho.est"          
[5] "est.gs"                  "est.gs.order"            "glm.p.order"            
> 
  > sig.est.gs.with.noise[[1]]
[1]  1  2  3  4  5  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
> sig.est.gs[[1]]
[1]  1  2  5  9 10 11 12 14 21 22 28 29 32 33 34 37 40 42 43 44 45 53 57 58 59 60 61 62 66 71 76

> sig.est.gs.with.noise[[2]]
[1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
> sig.est.gs[[2]]
[1]  5  9 10 11 12 14 21 22 25 28 42 43 44 45 53 57 58 59 60 61 62 66 72 77

> sig.est.gs.with.noise[[3]]
[1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
> sig.est.gs[[3]]
[1]  4  6 11 15 19 21 22 25 35 41 42 43 44 45 54 56 57 58 59 61 68 73 75 77 78

> sig.est.gs.with.noise[[4]]
[1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
> sig.est.gs[[4]]
[1]  1  2  4  8  9 10 12 14 17 18 22 30 31 32 34 38 42 43 44 45 52 53 55 56 60 61 62 69 72 73 74 76

> sig.est.gs.with.noise[[5]]
[1]  1  2  3  4  5  7  8  9 11 12 13 14 15 16 19 20 21 22 23 24 25 26 27 28
> sig.est.gs[[5]]
[1]  4  6 11 15 19 20 21 22 25 28 35 41 42 43 44 45 48 49 56 57 58 59 61 68 73 75 77 78

> sig.est.gs.with.noise[[6]]
[1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
> sig.est.gs[[6]]
[1]  4  8  9 10 12 13 14 15 17 21 22 30 33 38 42 43 44 45 53 55 61 63 72 73 74

> sig.est.gs.with.noise[[7]]
[1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
> sig.est.gs[[7]]
[1]  4  8  9 10 12 14 15 17 21 22 30 33 38 42 43 44 45 53 55 61 63 72 73 74

> sig.est.gs.with.noise[[8]]
[1]  2  3  4  6  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 27 28
> sig.est.gs[[8]]
[1]  3  4  9 10 11 12 13 14 15 17 21 22 30 33 40 42 43 44 45 53 55 60 61 62 63 71 72 74
> 
  
  > ##------------- confirm no signal genes --------------
> index <- c(1:78)
> ptm <- proc.time()
> scan.all.CYPs.against.noise.insig.genes <- apply(matrix(1:46,nr=1), 2, function(i){test.sig.against.noise(CYPs.sorted.m[,i],input.genes,index[-sig.est.gs[[i]] ] )})
> proc.time()-ptm
user  system elapsed 
6.02    0.17    6.19 
> 
  > insig.est.gs.with.noise <- lapply(scan.all.CYPs.against.noise.insig.genes, function(x){gs <- x[[1]][["est.gs"]];return( which(gs[c(1:x[[2]])] > max(gs[-c(1:x[[2]])] ) )); })
> 
  > 
  > 
  > dev.new()
> 
  > plot(scan.all.CYPs.against.noise.insig.genes[[1]][[1]][["est.gs"]])
> 
  > insig.est.gs.with.noise[[1]]
[1]  1  3  6  9 13 15 16 18 19 20 21 22 23 24 26 30 31 32 39 45 46
> index[-sig.est.gs[[1]] ]
[1]  3  4  6  7  8 13 15 16 17 18 19 20 23 24 25 26 27 30 31 35 36 38 39 41 46 47 48 49 50 51 52 54 55 56 63 64
[37] 65 67 68 69 70 72 73 74 75 77 78
> insig.est.gs.with.noise[[2]]
[1] 20 21 22 23 25 28 30 31 38 49 52
> index[-sig.est.gs[[2]] ]
[1]  1  2  3  4  6  7  8 13 15 16 17 18 19 20 23 24 26 27 29 30 31 32 33 34 35 36 37 38 39 40 41 46 47 48 49 50
[37] 51 52 54 55 56 63 64 65 67 68 69 70 71 73 74 75 76 78
> insig.est.gs.with.noise[[3]]
[1]  8 14 16 18 23 30 43 51
> index[-sig.est.gs[[3]] ]
[1]  1  2  3  5  7  8  9 10 12 13 14 16 17 18 20 23 24 26 27 28 29 30 31 32 33 34 36 37 38 39 40 46 47 48 49 50
[37] 51 52 53 55 60 62 63 64 65 66 67 69 70 71 72 74 76
> insig.est.gs.with.noise[[4]]
[1] 10 23 36
> index[-sig.est.gs[[4]] ]
[1]  3  5  6  7 11 13 15 16 19 20 21 23 24 25 26 27 28 29 33 35 36 37 39 40 41 46 47 48 49 50 51 54 57 58 59 63
[37] 64 65 66 67 68 70 71 75 77 78
> insig.est.gs.with.noise[[5]]
[1]  8 12 13 14 15 17 19 21 28 34 40 48 49 50
> index[-sig.est.gs[[5]] ]
[1]  1  2  3  5  7  8  9 10 12 13 14 16 17 18 23 24 26 27 29 30 31 32 33 34 36 37 38 39 40 46 47 50 51 52 53 54
[37] 55 60 62 63 64 65 66 67 69 70 71 72 74 76
> 
  > save(CYPs.sorted, input.genes, CYPs.sorted.m, input.genes.m, scan.all.CYPs, scan.all.CYPs.against.noise, sig.est.gs, sig.est.gs.with.noise, scan.all.CYPs.against.noise.insig.genes, insig.est.gs.with.noise, file="U:/Data/Danxin/CYP3A4_Microarray/Sobol_selection_result2.Rdata")
> 
  
  ##------------------------------------- 01 07 2016 post-analysis summary
  
  > load("U:/Data/Danxin/CYP3A4_Microarray/Sobol_selection_result2.Rdata")
> ls()
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
> 
  
  > names(scan.all.CYPs[[1]])
[1] "glm.coef.with.intercept" "glm.p.with.intercept"   
[3] "input.se.est"            "input.rho.est"          
[5] "est.gs"                  "est.gs.order"           
[7] "glm.p.order"      

## ---------fdr picks, significance level 0.05:##----------- fdr adjusted p-value selection (fdr=0.05) ------------
glm.p.c <- do.call(rbind, lapply(scan.all.CYPs, function(x){x[["glm.p.with.intercept"]]}))
glm.p.c.fdr <- t(apply(glm.p.c, 1, function(x){p.adjust(x,"fdr")}))
## number of variables selected by controlling fdr at 0.05
##table(apply(glm.p.c.fdr, 1, function(x){length(which(x<=0.05))}))
sig.glm.p.fdr <- apply(glm.p.c.fdr, 1, function(x){which(x<=0.05)})
sig.glm.p.fdr

table(apply(glm.p.c.fdr, 1, function(x){length(which(x<=0.05))}))

##---------- sobol ranking:
input.gene.names <- names(input.genes)
> input.gene.names
[1] "AHR"     "AHR.1"   "AHRR"    "ARNT"    "ARNT.1"  "CEBPA"   "CEBPB"   "CEBPD"   "CEBPG"  
[10] "CEBPG.1" "DBP"     "ESR1"    "ESR1.1"  "ESR1.2"  "FOXA1"   "FOXA1.1" "FOXA2"   "FOXA3"  
[19] "HNF4A"   "HNF4A.1" "HNF4A.2" "HNF4A.3" "HNF4G"   "HNF4G.1" "NCOA1"   "NCOA2"   "NCOA3"  
[28] "NCOA3.1" "NCOR1"   "NCOR2"   "NCOR2.1" "NCOR2.2" "NFE2L2"  "NR0B1"   "NR0B2"   "NR1D2"  
[37] "NR1D2.1" "NR1H2"   "NR1H2.1" "NR1H3"   "NR1H4"   "NR1I2"   "NR1I2.1" "NR1I3"   "NR1I3.1"
[46] "NR2F1"   "NR2F1.1" "NR2F2"   "NR2F2.1" "NR3C1"   "NR5A2"   "ONECUT1" "PGRMC1"  "PPARA"  
[55] "PPARA.1" "PPARA.2" "PPARA.3" "PPARA.4" "PPARA.5" "PPARD"   "PPARD.1" "PPARD.2" "PPARG"  
[64] "PPARG.1" "PPARG.2" "RXRA"    "RXRA.1"  "RXRB"    "RXRB.1"  "RXRG"    "RXRG.1"  "THRA"   
[73] "THRB"    "USF1"    "VDR"     "VDR.1"   "VDR.2"   "YY1"    
> 
  
  Sobol.picks <- lapply(scan.all.CYPs, function(x){index <- order(x[["est.gs"]],decreasing=T); gi <- input.gene.names[index]; return(gi[1:length(which(x[["est.gs"]]>mean(x[["est.gs"]]) ) )])})

> Sobol.picks <- lapply(scan.all.CYPs, function(x){index <- order(x[["est.gs"]],decreasing=T); gi <- input.gene.names[index]; return(gi[1:length(which(x[["est.gs"]]>mean(x[["est.gs"]]) ) )])})
> 
  > Sobol.picks[1:5]
$CYP1A1
[1] "NR1I3.1" "NR1I3"   "VDR.1"   "PPARD"   "NCOR1"   "PGRMC1"  "CEBPG"   "PPARD.2" "ESR1"   
[10] "PPARA.4" "ARNT.1"  "DBP"     "ESR1.2"  "HNF4A.3" "PPARA.5" "PPARD.1" "PPARA.3" "CEBPG.1"
[19] "NR1D2.1" "NR1I2"   "RXRA"    "AHR.1"   "AHR"     "HNF4A.2" "NCOA3.1" "NR1H3"   "NR1I2.1"
[28] "NCOR2.2" "RXRG.1"  "NR0B1"   "NFE2L2" 

$CYP1A1.1
[1] "DBP"     "NR1I3.1" "NR1I3"   "CEBPG"   "PPARD"   "PPARD.2" "PPARA.4" "RXRA"    "PPARA.5"
[10] "CEBPG.1" "PPARD.1" "PGRMC1"  "HNF4A.3" "NCOA3.1" "HNF4A.2" "NR1I2"   "ESR1"    "PPARA.3"
[19] "ESR1.2"  "NR1I2.1" "ARNT.1"  "NCOA1"   "VDR.2"   "THRA"   

$CYP1B1
[1] "HNF4A.2" "HNF4A.3" "NR1I2"   "NR1I2.1" "NR0B2"   "ARNT"    "PPARA.4" "PPARA.5" "FOXA1"  
[10] "RXRB"    "YY1"     "HNF4A"   "VDR.2"   "DBP"     "VDR"     "THRB"    "NR1I3.1" "CEBPA"  
[19] "PPARA.3" "NR1H4"   "NR1I3"   "PPARD.1" "PPARA.2" "NCOA1"   "PPARA"  

$CYP1B1.1
[1] "NR0B1"   "VDR.1"   "THRA"    "PGRMC1"  "ARNT"    "ESR1"    "NCOR2"   "CEBPG.1" "ESR1.2" 
[10] "USF1"    "AHR"     "NR1I3.1" "AHR.1"   "NR1I2"   "RXRB.1"  "NR1I3"   "CEBPG"   "NR1I2.1"
[19] "NCOR2.2" "PPARA.1" "NCOR2.1" "FOXA2"   "THRB"    "PPARD.1" "PPARD"   "PPARA.2" "FOXA3"  
[28] "PPARD.2" "HNF4A.3" "NR1H2"   "CEBPD"   "ONECUT1"

$CYP1B1.2
[1] "HNF4A.2" "HNF4A.3" "NR1I2"   "NR1I2.1" "NR0B2"   "ARNT"    "PPARA.5" "PPARA.4" "FOXA1"  
[10] "RXRB"    "VDR"     "VDR.2"   "CEBPA"   "YY1"     "DBP"     "HNF4A.1" "NR1I3.1" "HNF4A"  
[19] "THRB"    "NR1I3"   "PPARA.3" "NR2F2.1" "NR2F2"   "NR1H4"   "NCOA1"   "PPARD.1" "PPARA.2"
[28] "NCOA3.1"

> length(Sobol.picks)
[1] 46
> names(CYPs.sorted)
[1] "CYP1A1"    "CYP1A1.1"  "CYP1B1"    "CYP1B1.1"  "CYP1B1.2"  "CYP2A6"    "CYP2A6.1" 
[8] "CYP2B6"    "CYP2B6.1"  "CYP2C18"   "CYP2C18.1" "CYP2C19"   "CYP2C19.1" "CYP2C8"   
[15] "CYP2C8.1"  "CYP2C9"    "CYP2C9.1"  "CYP2D6"    "CYP2D6.1"  "CYP2E1"    "CYP2E1.1" 
[22] "CYP2F1"    "CYP2F1.1"  "CYP2J2"    "CYP2J2.1"  "CYP2R1"    "CYP2R1.1"  "CYP2U1"   
[29] "CYP2U1.1"  "CYP3A4"    "CYP3A4.1"  "CYP3A4.2"  "CYP3A43"   "CYP3A5"    "CYP3A5.1" 
[36] "CYP3A7"    "CYP3A7.1"  "CYP4F11"   "CYP4F11.1" "CYP4F12"   "CYP4F12.1" "CYP51A1"  
[43] "CYP51A1.1" "CYP51A1.2" "CYP7A1"    "CYP7A1.1" 
> Sobol.picks[30:35]
$CYP3A4
[1] "ESR1.2"  "ESR1"    "PGRMC1"  "NR1I3"   "NR1I2.1" "NR1I2"   "USF1"    "FOXA2"   "NFE2L2" 
[10] "PPARA.1" "NR1I3.1" "HNF4A.2" "HNF4A.3" "NCOR2"   "NCOR2.2" "NR1H3"   "NR1H2"   "ARNT"   
[19] "AHR"     "AHR.1"   "CEBPG"   "CEBPG.1" "NCOR1"   "THRA"    "NR5A2"  

$CYP3A4.1
[1] "ESR1.2"  "ESR1"    "USF1"    "NR1I2"   "NR1I2.1" "PGRMC1"  "FOXA2"   "NR1I3"   "PPARA.1"
[10] "NCOR2"   "HNF4A.2" "NFE2L2"  "ARNT"    "HNF4A.3" "NR1I3.1" "THRA"    "NR1H2"   "CEBPG.1"
[19] "NCOR2.2" "AHR"     "AHR.1"   "CEBPG"   "NCOR1"   "NR1H3"   "PPARD.1" "DBP"    

$CYP3A4.2
[1] "USF1"    "ESR1.2"  "ESR1"    "NR1I2"   "NR1I2.1" "PGRMC1"  "PPARA.1" "FOXA2"   "NR1I3"  
[10] "HNF4A.2" "NFE2L2"  "NCOR2"   "HNF4A.3" "ARNT"    "CEBPG.1" "THRA"    "NR1I3.1" "NR1H2"  
[19] "AHR"     "NCOR2.2" "AHR.1"   "CEBPG"   "NR1H3"   "DBP"     "NCOR1"   "PPARD.1"

$CYP3A43
[1] "ESR1"    "NR1I3"   "ESR1.2"  "NCOR1"   "FOXA2"   "PGRMC1"  "NR1I3.1" "NR1I2"   "NR1I2.1"
[10] "HNF4A.2" "NCOR2"   "USF1"    "HNF4A.3" "NCOR2.2" "CEBPG"   "CEBPG.1" "NFE2L2"  "ARNT"   
[19] "PPARA.1" "AHR"     "AHR.1"   "PPARD"   "PPARD.1" "CEBPD"   "PPARD.2" "NR1H3"   "RXRA.1" 
[28] "THRA"   

$CYP3A5
[1] "ESR1.2"  "NR1I3"   "NR1I3.1" "ESR1"    "NCOR1"   "NR1I2"   "PGRMC1"  "NR1I2.1" "HNF4A.3"
[10] "FOXA2"   "HNF4A.2" "ARNT"    "NFE2L2"  "THRA"    "NR1H3"   "RXRA"    "CEBPG"   "PPARD.1"
[19] "NCOR2.2" "PPARD"   "NR5A2"   "CEBPG.1" "PPARA.1" "AHR"     "NR0B2"  

$CYP3A5.1
[1] "NCOR1"   "NR1I3.1" "NR1I3"   "ESR1.2"  "ESR1"    "HNF4A.3" "NR1I2"   "HNF4A.2" "NR1I2.1"
[10] "PGRMC1"  "RXRA"    "FOXA2"   "ARNT"    "THRA"    "PPARD.1" "PPARD"   "NR5A2"   "PPARA.4"
[19] "NR1H3"   "CEBPG"   "CEBPG.1" "NR0B2"   "PPARD.2"

> 
  
  
  ##---------- sobol ranking tested against artificial noises:
  ## the basic idea is: although the amount of information collected in experiments is fixed, 
  ## we can simulated as much as independent artificial noise as we want, and use them as the control group
  ## for testing for significance of gene effect. 
  
  
  dev.new()
par(mfrow=c(1,2))
plot(scan.all.CYPs.against.noise[[1]][[1]][["est.gs"]])
abline(v=length(sig.est.gs[[1]]), col="red")
plot(scan.all.CYPs.against.noise.insig.genes[[1]][[1]][["est.gs"]])
abline(v= (78-length(sig.est.gs[[1]])) , col="red")

dev.new()
plot(scan.all.CYPs[[1]][["est.gs"]])
abline(h= mean(scan.all.CYPs[[1]][["est.gs"]]) , col="red")


dev.new()
par(mfrow=c(1,2))
plot(scan.all.CYPs.against.noise[[2]][[1]][["est.gs"]])
abline(v=length(sig.est.gs[[2]]), col="red")
plot(scan.all.CYPs.against.noise.insig.genes[[2]][[1]][["est.gs"]])
abline(v= (78-length(sig.est.gs[[2]])) , col="red")

dev.new()
plot(scan.all.CYPs[[2]][["est.gs"]])
abline(h= mean(scan.all.CYPs[[2]][["est.gs"]]) , col="red")


dev.new()
par(mfrow=c(1,2))
plot(scan.all.CYPs.against.noise[[3]][[1]][["est.gs"]])
abline(v=length(sig.est.gs[[3]]), col="red")
plot(scan.all.CYPs.against.noise.insig.genes[[3]][[1]][["est.gs"]])
abline(v= (78-length(sig.est.gs[[3]])) , col="red")

dev.new()
plot(scan.all.CYPs[[3]][["est.gs"]])
abline(h= mean(scan.all.CYPs[[3]][["est.gs"]]) , col="red")

























## ------------------ original code for evaluation of Sobol indices ---------------

library(MASS)
library(glmnet)


## caculation of the numerator of the sobol main index

t.gs <- function(i, coef, mse, rhom){
  rho <- rhom[i,-i]
  return ( ( coef[i]+ matrix(coef[-i],nr=1) %*% ( matrix(mse[-i],nc=1) * matrix(rho,nc=1) ) / mse[i] )^2 * mse[i]^2  )
}

## d: dimension of multinormal, the total number of inputs 
## rho: correlation between inputs
## n : sample size
## r : the first r columns in the input matrix are the true predictors, the remainings are fake predictors

one.sample.fit <- function(d=55, rho.t=0.1, rho.f=0.8, n=1000, r=5){
  ## multinormal mean
  mu <- runif(d, min=-50, max=50)
  ## multinormal marginal S.E.
  mse <- runif(d, min=0, max=10)
  
  ## multinormal covariance matrix, true predictors and fake predictors are correlated
  varm <- matrix(mse[1:r], nc=1) %*% matrix(mse[1:r], nr=1) * rho.t
  varm[row(varm)==col(varm)] <- varm[row(varm)==col(varm)] / rho.t
  varm1 <- varm
  if(d>r){
    varm <- matrix(mse[(r+1):d], nc=1) %*% matrix(mse[(r+1):d], nr=1) * rho.f
    varm[row(varm)==col(varm)] <- varm[row(varm)==col(varm)] / rho.f
    varm2 <- matrix(0,nr=r,nc=(d-r))
    varm1 <- cbind(varm1, varm2)
    varm <- cbind(t(varm2), varm)
    varm <- rbind(varm1, varm)
  }
  
  ## input is n (sample size) by d (multinormal dimension) matrix
  input <- mvrnorm(n, mu, varm)
  
  # generate ture coefficients for generating mean of response, except for beta0, all coefficients are positive. 
  coef <- runif(r, min=-1, max=1)
  beta0 <- runif(1, min=-1, max=1)
  
  # first r columns for generating mean of response y
  y.mean <- t(coef) %*% t(input[,1:r ]) + beta0
  y.se <- runif(1,min=3.5, max=5)
  
  ## generate one sample of Gaussian GLM response
  y <- apply(matrix(y.mean, nc=1), 1, function(x){rnorm(1, x, y.se)})
  
  ## fit the sample with GLM via Iteratively Reweighted Least Square (IRLS)
  fit <- glm(y~., data=data.frame(input))
  glm.coef <- fit$coef
  glm.p <- summary(fit)$coef[,4]; glm.p <- glm.p[-1];
  glm.p.order <- order(glm.p, na.last=TRUE)
  
  ## calculate emprical estimate of sobol ranking using glm.coef
  coef.est <- glm.coef[-1]
  mse.est <- apply(input, 2, function(x){var(x)^0.5})
  rhom.est <- cor(input)
  
  est.gs <- apply(matrix(1:d, nc=1), 1, function(i){t.gs(i, coef.est, mse.est, rhom.est)})  
  est.gs.order <- order(est.gs, decreasing=TRUE, na.last=TRUE)
  
  ## calculate emprical estimate of sobol ranking using glmnet.coef ??????????????????????????
  ## lasso penalty (alpha=1, default, indifference to very correlated predictors)???
  ## ridge-regression penalty (alpha=0, good for correlated inputs) ???
  
  my.list <- list(input.mu=mu, input.se=mse, input.cov=varm, true.beta0=beta0, true.coef=coef, true.y.mean=y.mean, true.y.se=y.se, glm.coef=glm.coef, glm.p=glm.p, input.se.est=mse.est, input.rho.est=rhom.est, est.gs=est.gs, est.gs.order=est.gs.order, glm.p.order=glm.p.order)
  
  return(my.list)
}


##------------------------------------------------------------------
##---------- select the best variable combination ------------------
##------------------------------------------------------------------
## cs: combination size 
## coef: the vector glm coef estimates, minus the intercept estimate
##covm: empirical estimate of the inputs covariance matrix

multiple.variable.selection.sub <- function(index, coef, covm){
  eta <- matrix(coef[index],nc=1) + solve(covm[index, index]) %*% covm[index,-index] %*% matrix(coef[-index],nc=1)
  return( t(eta) %*% covm[index, index] %*% eta )
}
## FYI: if solve() return error can try ginv() after loading MASS package.

multiple.variable.selection <- function(cs=5, coef, covm){
  index.m <- combn(c(1:length(coef)), cs)   ## index.m is a matrix with cs rows.
  mv.gs <- apply(index.m, 2, function(x){multiple.variable.selection.sub(x, coef,covm)})
  cc <- rbind(index.m, mv.gs)
  ccc <- cc[,order(mv.gs, decreasing=T)]
  return(list(ccc[-(cs+1),], ccc[(cs+1),]))
}


##---------
glm.p.c <- do.call(rbind, lapply(multi.sample.fit, function(x){x[["glm.p"]]}))
glm.p.c.fdr <- t(apply(glm.p.c, 1, function(x){p.adjust(x,"fdr")}))
## number of variables selected by controlling fdr at 0.05
table(apply(glm.p.c.fdr, 1, function(x){length(which(x<=0.05))}))

##----------- fdr adjusted p-value selection (fdr=0.05) ------------
sig.glm.p.fdr <- apply(glm.p.c.fdr[,-1], 1, function(x){which(x<=0.05)})
sig.glm.p.fdr.v <- unlist(sig.glm.p.fdr)
length(which(sig.glm.p.fdr.v<=5))
length(sig.glm.p.fdr.v)

glm.p.c.fdr.order <- t(apply(glm.p.c.fdr[,-1],1,function(x){order(x,decreasing=FALSE, na.last=TRUE)}))
dim(glm.p.c.fdr.order)

length(which(est.gs.order.c[,1:5]<=5))
length(which(glm.p.order.c[,1:5]<=5))
length(which(glm.p.c.fdr.order[,1:5]<=5))

length(which(est.gs.order.c[,1:10]<=5))
length(which(glm.p.order.c[,1:10]<=5))
length(which(glm.p.c.fdr.order[,1:10]<=5))

##---------- simple average cutoff ------------------------
est.gs.c <- do.call(rbind, lapply(multi.sample.fit, function(x){x[["est.gs"]]}))

sig.est.gs.hclust <- unlist(apply(est.gs.c, 1, function(x){which( x > mean(x) )}))
length(which(sig.est.gs.hclust<=5))
length(sig.est.gs.hclust)

##-----------
table(apply(est.gs.c, 1, function(x){length(which( x > mean(x) ))}))
dev.new()
plot(est.gs.c[1,])



##------------------------------------------------------ code before adding more CYPs ---------------------------------

data <- read.table("U:/Data/Danxin/CYP3A4_Microarray/CYP3A4_Microarray.csv", sep=",", header=T)
data[1:5,1:10]

gene.names <- data$ORF
length(gene.names)
length(unique(gene.names))
tb <- table(gene.names)

unique.gene.names <- names(tb)

> unique.gene.names
[1] "AHR"     "AHRR"    "ARNT"    "CEBPA"   "CEBPB"   "CEBPD"   "CEBPG"  
[8] "CYP3A4"  "DBP"     "ESR1"    "FOXA1"   "FOXA2"   "FOXA3"   "HNF4A"  
[15] "HNF4G"   "NCOA1"   "NCOA2"   "NCOA3"   "NCOR1"   "NCOR2"   "NFE2L2" 
[22] "NR0B1"   "NR0B2"   "NR1D2"   "NR1H2"   "NR1H3"   "NR1H4"   "NR1I2"  
[29] "NR1I3"   "NR2F1"   "NR2F2"   "NR3C1"   "NR5A2"   "ONECUT1" "PGRMC1" 
[36] "PPARA"   "PPARD"   "PPARG"   "RXRA"    "RXRB"    "RXRG"    "THRA"   
[43] "THRB"    "USF1"    "VDR"     "YY1" 

csv.table <- read.table("U:/Data/Danxin/consider_gender_effects/CYP3A4_GTEx_Liver_subject_info_matched_final.csv", sep=",", header=T)
dim(csv.table)
names(csv.table)
original.gene.list <- names(csv.table[,9:54])
original.gene.list <- c(original.gene.list[1:7],"CYP3A4",original.gene.list[8:45])
original.gene.list
unique.gene.names



index.m <- combn(c(1:81),3)
dim(index.m)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
dim(gene.names.m)
filter <- apply(matrix(c(1:85320),nc=1),1,function(i){comb <- gene.names.m[,i];if(comb[1]==comb[2]|comb[1]==comb[3]|comb[3]==comb[2]){return(i)}else{return(0)}; })
length(filter)
index.m <- index.m[,which(filter==0)]
dim(index.m) ##85320 -> 81114
gene.names.m <- gene.names.m[,which(filter==0)]
gene.names.m[,sample(c(1:81114),10)]
filter <- apply(matrix(c(1:81114),nc=1),1,function(i){comb <- gene.names.m[,i];if(comb[1]=="CYP3A4"|comb[2]=="CYP3A4"|comb[3]=="CYP3A4"){return(i)}else{return(0)}; })
index.m <- index.m[,which(filter==0)]
dim(index.m) ##81114 -> 71870
gene.names.m <- gene.names.m[,which(filter==0)]
which(gene.names.m=="CYP3A4")



## this type of data is even better: buy treating each coloumn as a sample, we can actually identify genes of which the changes between one sample vs the pooled sample predicts the changes in CYP3A4 changes the best
## we should do the same thing using the RNA-seq data. In that case, we can actually get a lot more sample than the original 137 "model samples". Just make sure "independence" across "samples" when doing so. 
## check outlier of each gene (across samples) to see if the outlier corresponds to a outlier in CYP3A4: yes -> keep the sample; no -> delete the sample.
## check gene outliers which have very low expression measure -- which cannot be done since all measure are relative measure. 


write.csv(data.C, "U:/Data/Danxin/CYP3A4_Microarray/CYP3A4_Microarray.csv")

save(input.genes, CYPs.sorted, file="U:/Data/Danxin/CYP3A4_Microarray/CYP3A4_Microarray_input_CYPs_separated.Rdata")

##------------------------------- pre-processing is done

library(MASS)


## caculation of the numerator of the sobol main index

t.gs <- function(i, coef, mse, rhom){
  rho <- rhom[i,-i]
  return ( ( coef[i]+ matrix(coef[-i],nr=1) %*% ( matrix(mse[-i],nc=1) * matrix(rho,nc=1) ) / mse[i] )^2 * mse[i]^2  )
}

## d: dimension of multinormal, the total number of inputs 
## rho: correlation between inputs
## n : sample size
## r : the first r columns in the input matrix are the true predictors, the remainings are fake predictors

one.y.fit <- function(y, input){
  d <- dim(input)[2]
  
  ## fit the sample with GLM
  fit <- glm(y~., data=data.frame(input))
  glm.coef <- fit$coef
  glm.p <- summary(fit)$coef[,4]; glm.p <- glm.p[-1];
  glm.p.order <- order(glm.p, na.last=TRUE)
  
  ## calculate emprical estimate of sobol ranking
  coef.est <- glm.coef[-1]
  mse.est <- apply(input, 2, function(x){var(x, na.rm=T)^0.5})
  rhom.est <- cor(input, use="pairwise.complete.obs")
  
  est.gs <- apply(matrix(1:d, nc=1), 1, function(i){t.gs(i, coef.est, mse.est, rhom.est)})  
  est.gs.order <- order(est.gs, decreasing=TRUE, na.last=TRUE)
  
  #my.list <- list(input.mu=mu, input.se=mse, input.cov=varm, true.beta0=beta0, true.coef=coef, true.y.mean=y.mean, true.y.se=y.se, glm.coef=glm.coef, glm.p=glm.p, input.se.est=mse.est, input.rho.est=rhom.est, est.gs=est.gs, est.gs.order=est.gs.order, glm.p.order=glm.p.order)
  my.list <- list(glm.coef=glm.coef, glm.p=glm.p, input.se.est=mse.est, input.rho.est=rhom.est, est.gs=est.gs, est.gs.order=est.gs.order, glm.p.order=glm.p.order)
  
  return(my.list)
}

##--------------------- step 1: initial selection (top 10) without artificial fake inputs control ---------------------------
result.1  <- apply(CYPs.sorted, 2, function(y){one.y.fit(y, input.genes)})

[1] "glm.coef"      "glm.p"         "input.se.est"  "input.rho.est" "est.gs"        "est.gs.order" 
[7] "glm.p.order"  

names(result.1[[1]])
result.1[[1]][["est.gs.order"]]
result.1[[1]][["glm.p.order"]]

plot(result.1[[1]][["est.gs"]])
plot(result.1[[1]][["glm.p"]])

##---------
glm.p.c <- do.call(rbind, lapply(result.1 , function(x){x[["glm.p"]]}))
glm.p.c.fdr <- t(apply(glm.p.c, 1, function(x){p.adjust(x,"fdr")}))
## number of variables selected by controlling fdr at 0.05
table(apply(glm.p.c.fdr, 1, function(x){length(which(x<=0.05))}))

##----------- fdr adjusted p-value selection (fdr=0.05) ------------
sig.glm.p.fdr <- apply(glm.p.c.fdr[,-1], 1, function(x){which(x<=0.05)})
sig.glm.p.fdr ## here shows genes picked up by fdr-adjusted glm p-value

##---------
est.gs.order.c <- do.call(rbind, lapply(result.1, function(x){x[["est.gs.order"]]}))  ## est.gs.order.c[,1:10] are the top 10 picks
glm.p.c.fdr.order <- t(apply(glm.p.c.fdr[,-1],1,function(x){order(x,decreasing=FALSE, na.last=TRUE)}))

est.gs.order.c[1,1:10]
glm.p.c.fdr.order[1,1:10]
table(c(est.gs.order.c[1,1:10], glm.p.c.fdr.order[1,1:10]))


> dim(input.genes)
[1] 427  78
> choose(78,1)
[1] 78
> choose(78,2)
[1] 3003
> choose(78,3)
[1] 76076


##--------------------- step 2: add artificial fake inputs to top 10 picks for comparison ------------------------

gene.mean <- apply(input.genes,2, function(x){mean(x, na.rm=T)})
gene.mean.pool <- quantile(gene.mean, c(1:20)/20)[-20]

gene.se <- apply(input.genes,2, function(x){var(x, na.rm=T)^0.5})
gene.se.pool <- quantile(gene.se,c(1:20)/20)[-20]

gene.corr <- cor(input.genes, use="pairwise.complete.obs")
gene.corr.offdiag <- gene.corr[row(gene.corr)!=col(gene.corr)]
gene.corr.pool <- quantile(gene.corr.offdiag, c(1:20)/20)[-20]

n.fake <- 30
mu.fake <- sample(gene.mean.pool, n.fake, replace=T)
se.fake <- sample(gene.se.pool, n.fake, replace=T)
rho.fake <- matrix(sample(gene.corr.pool, n.fake^2 , replace=T), nr=n.fake)
diag(rho.fake) <- 1

> matrix(1:9,nr=3)
[,1] [,2] [,3]
[1,]    1    4    7
[2,]    2    5    8
[3,]    3    6    9
> a <- matrix(1:9,nr=3)
> lower.tri(a)
[,1]  [,2]  [,3]
[1,] FALSE FALSE FALSE
[2,]  TRUE FALSE FALSE
[3,]  TRUE  TRUE FALSE
> a[lower.tri(a)] <- a[upper.tri(a)]
> a
[,1] [,2] [,3]
[1,]    1    4    7
[2,]    4    5    8
[3,]    7    8    9



## multinormal marginal S.E.
mse <- runif(d, min=0, max=10)

## multinormal covariance matrix, true predictors and fake predictors are correlated
varm <- matrix(mse[1:r], nc=1) %*% matrix(mse[1:r], nr=1) * rho.t
varm[row(varm)==col(varm)] <- varm[row(varm)==col(varm)] / rho.t
varm1 <- varm














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

index.m <- combn(c(1:21),3)
dim(index.m)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})

gene <- rpkm[,-7]

ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP3A4, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
hist(Residual.Deviance.Deduction,breaks=100)

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
summary(Residual.Deviance.Deduction)

> summary(Residual.Deviance)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
28.51   77.51   94.80   96.86  118.00  157.40 
> summary(Residual.Deviance.Deduction)
Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0001014 0.0364800 0.0659600 0.0820200 0.1115000 0.4092000

large.interaction.index <- which((Residual.Deviance<=) & (Residual.Deviance.Deduction>=) )

gene.names.m[,large.interaction.index]
results[,large.interaction.index]
index.m[,large.interaction.index]



goodfit.index <- which((Residual.Deviance<=37.5) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]

> goodfitr[,order(goodfitr[8,])]
[,1]               [,2]                [,3]                [,4]                [,5]               
[1,] "27"               "10"                "25"                "22"                "7"                
[2,] "33"               "27"                "27"                "27"                "27"               
[3,] "37"               "33"                "33"                "33"                "33"               
[4,] "NR1I3"            "FOXA1"             "NR1H4"             "NR1D2"             "CEBPG"            
[5,] "PGRMC1"           "NR1I3"             "NR1I3"             "NR1I3"             "NR1I3"            
[6,] "RXRA"             "PGRMC1"            "PGRMC1"            "PGRMC1"            "PGRMC1"           
[7,] "48.2488846865023" "48.8795705071804"  "51.3973602809268"  "56.5156461074226"  "51.0365553326472" 
[8,] "28.5076746597596" "29.2503645769412"  "33.3566687688686"  "33.5791881273189"  "34.1491779886894" 
[9,] "0.40915370697191" "0.401583027971893" "0.351004242502956" "0.405842621643341" "0.330887875051303"
[,6]                [,7]                [,8]                [,9]                [,10]              
[1,] "24"                "26"                "26"                "9"                 "27"               
[2,] "27"                "33"                "27"                "27"                "33"               
[3,] "33"                "37"                "33"                "33"                "35"               
[4,] "NR1H3"             "NR1I2"             "NR1I2"             "ESR1"              "NR1I3"            
[5,] "NR1I3"             "PGRMC1"            "NR1I3"             "NR1I3"             "PGRMC1"           
[6,] "PGRMC1"            "RXRA"              "PGRMC1"            "PGRMC1"            "PPARD"            
[7,] "56.7813733040797"  "45.9519634417571"  "52.2834010175321"  "56.4114404352313"  "50.2120174621509" 
[8,] "34.2875835297251"  "34.3237555082979"  "34.9533443161965"  "36.3356708624746"  "36.4198604653445" 
[9,] "0.396147336097248" "0.253051383717209" "0.331463836783003" "0.355881172646294" "0.274678407558564"
[,11]               [,12]               [,13]               [,14]               [,15]              
[1,] "27"                "21"                "27"                "23"                "26"               
[2,] "32"                "27"                "30"                "27"                "27"               
[3,] "33"                "33"                "33"                "33"                "37"               
[4,] "NR1I3"             "NR0B2"             "NR1I3"             "NR1H2"             "NR1I2"            
[5,] "ONECUT1"           "NR1I3"             "NR3C1"             "NR1I3"             "NR1I3"            
[6,] "PGRMC1"            "PGRMC1"            "PGRMC1"            "PGRMC1"            "RXRA"             
[7,] "50.6209315622547"  "55.8612320639316"  "57.0009359730306"  "50.7549463528293"  "60.7208356729699" 
[8,] "36.4986558753266"  "36.5303966517435"  "37.0024423644264"  "37.3013724017959"  "37.3427196763561" 
[9,] "0.278980952169168" "0.346051003494955" "0.350845003984958" "0.265069218229672" "0.385009786797462"




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
summary(Residual.Deviance.Deduction)

> summary(Residual.Deviance)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
24.95   69.93   92.59   92.70  116.70  158.00 
> summary(Residual.Deviance.Deduction)
Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.000195 0.037190 0.067560 0.079920 0.110300 0.363800

large.interaction.index <- which((Residual.Deviance<=) & (Residual.Deviance.Deduction>=) )

gene.names.m[,large.interaction.index]
results[,large.interaction.index]
index.m[,large.interaction.index]


goodfit.index <- which((Residual.Deviance<=30.5) )
length(goodfit.index)

goodfitr <- rbind(index.m[,goodfit.index], gene.names.m[,goodfit.index], results[,goodfit.index])
goodfitr[,order(goodfitr[8,])]

> goodfitr[,order(goodfitr[8,])]
[,1]                [,2]                [,3]                [,4]                [,5]               
[1,] "27"                "22"                "13"                "4"                 "10"               
[2,] "33"                "27"                "27"                "27"                "27"               
[3,] "37"                "33"                "33"                "33"                "33"               
[4,] "NR1I3"             "NR1D2"             "HNF4A"             "CEBPA"             "FOXA1"            
[5,] "PGRMC1"            "NR1I3"             "NR1I3"             "NR1I3"             "NR1I3"            
[6,] "RXRA"              "PGRMC1"            "PGRMC1"            "PGRMC1"            "PGRMC1"           
[7,] "34.0233312667165"  "38.3394561461807"  "37.312739506758"   "38.2023048254285"  "34.614034118612"  
[8,] "24.9534370701286"  "25.730483509195"   "26.4045417498707"  "26.7462931862437"  "27.2466985150505" 
[9,] "0.266578664078685" "0.328877190873814" "0.292345132013468" "0.299877499316726" "0.212842443568291"
[,6]                [,7]                [,8]                [,9]                [,10]              
[1,] "27"                "25"                "6"                 "24"                "22"               
[2,] "33"                "27"                "27"                "27"                "23"               
[3,] "34"                "33"                "33"                "33"                "27"               
[4,] "NR1I3"             "NR1H4"             "CEBPD"             "NR1H3"             "NR1D2"            
[5,] "PGRMC1"            "NR1I3"             "NR1I3"             "NR1I3"             "NR1H2"            
[6,] "PPARA"             "PGRMC1"            "PGRMC1"            "PGRMC1"            "NR1I3"            
[7,] "38.3891437176878"  "35.7923925593426"  "38.8930853002228"  "39.2445564894459"  "42.3177814218746" 
[8,] "27.9581597298638"  "28.4283615660025"  "28.5148795439074"  "29.0163152429666"  "29.1260691948824" 
[9,] "0.271717026681632" "0.205742909785391" "0.266839353993237" "0.260628279726642" "0.311729768994298"
[,11]               [,12]               [,13]               [,14]               [,15]              
[1,] "11"                "27"                "27"                "27"                "20"               
[2,] "27"                "31"                "33"                "32"                "27"               
[3,] "33"                "33"                "35"                "33"                "33"               
[4,] "FOXA2"             "NR1I3"             "NR1I3"             "NR1I3"             "NFE2L2"           
[5,] "NR1I3"             "NR5A2"             "PGRMC1"            "ONECUT1"           "NR1I3"            
[6,] "PGRMC1"            "PGRMC1"            "PPARD"             "PGRMC1"            "PGRMC1"           
[7,] "38.3565645166958"  "39.2708329401909"  "37.5957987074233"  "37.2359737726973"  "38.8256540436075" 
[8,] "29.3086260449606"  "30.1064313675487"  "30.2123388767039"  "30.3216414440288"  "30.4652209231493" 
[9,] "0.235890220767734" "0.233364074212521" "0.196390556513474" "0.185689579944282" "0.215332705305312"

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

> pg.tb[order(pg.tb)]

AHR    AHRR    ARNT   CEBPA   CEBPD   CEBPG     DBP   FOXA3   HNF4G   NCOA1   NCOA2   NCOA3   NCOR1  NFE2L2 
1       1       1       1       1       1       1       1       1       1       1       1       1       1 
NR1H4   NR2F1   NR2F2   NR3C1   NR5A2 ONECUT1   PPARA   PPARG    RXRB    THRB     VDR     YY1   FOXA2    THRA 
1       1       1       1       1       1       1       1       1       1       1       1       2       2 
CEBPB   HNF4A   NR0B2   NR1H2   NR1H3   NR1I2   PPARD   FOXA1    USF1    ESR1   NCOR2    RXRA   NR1D2  PGRMC1 
3       3       3       3       3       3       3       4       4       5       6       7       8      41 
NR1I3 
60 
> pg.tb.1[order(pg.tb.1)]

VDR    THRA   HNF4G   NCOA1   NCOR1   NR1H4   CEBPD   PPARG    ARNT ONECUT1   FOXA1   NCOA2   CEBPB     YY1 
53      58      60      60      63      63      64      65      67      67      68      69      70      73 
CEBPA   FOXA3   NR2F1   NR2F2   NCOA3  NFE2L2    RXRB    USF1   PPARA    THRB     AHR    AHRR   NR5A2   HNF4A 
74      74      74      75      76      77      80      80      81      81      83      83      87      88 
NR3C1   FOXA2   NR1D2    RXRA   NR1H2     DBP   NCOR2   CEBPG   NR0B2  PGRMC1   NR1H3   PPARD    ESR1   NR1I2 
88      91      92      93      96     106     122     128     143     175     183     190     226     571 
NR1I3 
801 
> npg.tb.1[order(npg.tb.1)]

AHRR   FOXA2   NCOA2   NR5A2  NFE2L2   PPARG   CEBPA   HNF4G    ARNT   NR2F2    THRB   HNF4A     VDR   FOXA3 
84      87      97      99     100     102     103     108     110     111     117     118     118     119 
NR2F1    THRA    USF1   CEBPD    RXRB   NCOA1   NR1H4     AHR   NCOR1   PPARA   CEBPB ONECUT1     YY1   FOXA1 
123     125     126     127     135     136     140     147     150     152     153     159     164     169 
NR0B2   NCOA3   NR1D2    RXRA   NR3C1   CEBPG     DBP   NCOR2   NR1H3   NR1I2  PGRMC1   NR1H2   PPARD    ESR1 
181     193     202     214     238     241     248     272     281     287     375     387     479     630 
> npg.tb.2[order(npg.tb.2)]

PPARD  PGRMC1    RXRB   FOXA2  NFE2L2    AHRR   PPARG     VDR ONECUT1   HNF4G   CEBPA   NCOA2    ARNT   NR1H4 
189     270     277     281     281     289     289     294     296     299     301     307     310     313 
NCOA1   NR5A2    THRA   NCOR1   HNF4A   NR2F1    USF1   CEBPB   FOXA1   FOXA3   CEBPD   NR1H2   NR2F2     AHR 
314     321     321     327     330     333     333     343     346     350     351     375     375     380 
PPARA   NR1H3   NCOA3     YY1   NR0B2    THRB   NR1D2   NCOR2    RXRA   CEBPG     DBP   NR3C1 
389     394     398     400     422     429     436     461     472     489     506     530 
> npg.tb[order(npg.tb)]

CEBPG   NR3C1    RXRA   NR0B2   NR1D2   NCOA3     YY1    THRB   PPARA     AHR   FOXA1   CEBPB   NR2F2   FOXA3 
2       4      75     112     123     193     223     233     238     250     274     292     299     317 
CEBPD    USF1   NCOR1   HNF4A   NR2F1 ONECUT1   NR1H4   NCOA1   NR5A2    THRA    RXRB    ARNT   CEBPA   NCOA2 
318     318     320     322     330     338     344     350     353     355     368     373     382     387 
HNF4G     VDR   FOXA2  NFE2L2    AHRR   PPARG 
393     395     400     402     404     404







##------------------------------ code for the previous RNA-seq data
## Null Deviance = 2(LL(Saturated Model) - LL(Null Model)) on df = df_Sat - df_Null
## Residual Deviance = 2(LL(Saturated Model) - LL(Proposed Model)) df = df_Sat - df_Res

interaction3gene <- function(y,gender,gene){
  xm <- data.frame(gender,gene)
  m1 <- glm(y ~ .,data=xm); ##print(summary(m1));
  m2 <- glm(y ~ .^4, data=xm);##print(summary(m2));
  return(c(m1$deviance, m2$deviance, (m1$deviance-m2$deviance)/m1$deviance))
  
}
index.m <- combn(c(1:43),3)
dim(index.m)
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP3A4, gender,gene))})
proc.time()-ptm
##------------------------------ code for the previous RNA-seq data



data[1:5,1:8]
gene.tb <- t(data[,5:431]); colnames(gene.tb) <- gene.names;
gene.tb[1:5,1:8]

which(gene.names=="CYP3A4") ## 11 12 13
log.CYP3A4.D <- gene.tb[,11]

interaction3gene <- function(y,gene){
  xm <- data.frame(gene)
  m1 <- glm(y ~ .,data=xm); ##print(summary(m1));
  m2 <- glm(y ~ .^4, data=xm);##print(summary(m2));
  return(c(m1$deviance, m2$deviance, (m1$deviance-m2$deviance)/m1$deviance))
  
}
ptm <- proc.time()
results <- apply(index.m, 2, function(i){gene <- gene.tb[,i]; return(interaction3gene(log.CYP3A4.D, gene))})
proc.time()-ptm

Residual.Deviance.Deduction <- results[3,]
dev.new()
hist(Residual.Deviance.Deduction,breaks=100)

Residual.Deviance <- results[2,]
dev.new()
hist(Residual.Deviance,breaks=100)


index.m[,which(Residual.Deviance.Deduction>=0.05 & Residual.Deviance<=65)]
gene.names.m[,which(Residual.Deviance.Deduction>=0.05 & Residual.Deviance<=65)]
results[,which(Residual.Deviance.Deduction>=0.05 & Residual.Deviance<=65)]


> index.m[,which(Residual.Deviance.Deduction>=0.05 & Residual.Deviance<=65)]
[,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
[1,]    8   15   17   17   17   17   24   25
[2,]   15   48   24   25   47   48   36   36
[3,]   25   61   36   36   61   61   58   58
> gene.names.m[,which(Residual.Deviance.Deduction>=0.05 & Residual.Deviance<=65)]
[,1]    [,2]    [,3]     [,4]     [,5]    [,6]    [,7]     [,8]    
[1,] "CEBPD" "ESR1"  "ESR1"   "ESR1"   "ESR1"  "ESR1"  "HNF4A"  "HNF4A" 
[2,] "ESR1"  "NR1I3" "HNF4A"  "HNF4A"  "NR1I3" "NR1I3" "NFE2L2" "NFE2L2"
[3,] "HNF4A" "PPARA" "NFE2L2" "NFE2L2" "PPARA" "PPARA" "PPARA"  "PPARA" 
> results[,which(Residual.Deviance.Deduction>=0.05 & Residual.Deviance<=65)]
[,1]        [,2]        [,3]        [,4]       [,5]        [,6]
[1,] 70.3262649 65.96760793 67.12186713 64.42281558 67.9827806 66.05411055
[2,] 64.7198199 61.78228876 63.48765629 60.97463284 64.0331782 61.26384951
[3,]  0.0797205  0.06344506  0.05414347  0.05352425  0.0580971  0.07252026
[,7]        [,8]
[1,] 68.02441991 64.92921385
[2,] 63.95360642 61.09553918
[3,]  0.05984341  0.05904391


