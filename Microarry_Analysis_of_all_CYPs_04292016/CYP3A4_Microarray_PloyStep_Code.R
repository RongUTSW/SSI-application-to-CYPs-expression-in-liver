
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

##------------------------------------------------------- example 1:

##y1 <- CYPs.sorted.m[,30] ## CYP3A4

y <- CYPs.sorted.m[,30]
input <- apply(input.genes.m, 2, function(x){ifelse(is.na(x), mean(x, na.rm=TRUE), x)})

index.2m.all <- combn(c(1:78), 2)
dim(index.2m.all)

gene.names.2m.all <- apply(index.2m.all, 2, function(i){gene.names[i]})

all.gene.pairs.3 <- apply(index.2m.all, 2, function(i){SI.comb(y1, input[,i], k=3)})

order.2m.3 <- order(-all.gene.pairs.3[2,], all.gene.pairs.3[1,],decreasing=TRUE)

all.gene.pairs.3.stepwise <- all.gene.pairs.3
order.2m.3.stepwise <- order.2m.3

gene.names.2m.all[,order.2m.3.stepwise][,1:30]
all.gene.pairs.3.stepwise[,order.2m.3.stepwise][,1:30]



##------------------------------------------------------- example 2:
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


## ------------------------------------------------------ log file:
## U:\\Data\\Danxin\\CYP3A4_Microarray\\CYP3A4_Microarray_RScript_withMoreCYPs 02 25 2016.txt