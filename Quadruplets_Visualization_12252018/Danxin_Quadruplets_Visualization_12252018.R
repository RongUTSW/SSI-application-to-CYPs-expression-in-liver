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
