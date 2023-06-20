library(circlize)
library(networkD3)
library(igraph)
library(reshape2)
library(sysfonts)
library(showtext)
library(caret)
library(pliable)

# this script generates the chord diagrams (Figure 11 and 12)

plot.circ<-function(beta,theta,X,Z,adj = c(.15, .8),cex = .5){
  
  data_mat = merge(as.data.frame(beta), as.data.frame(theta), by='row.names')
  rownames(data_mat) <- data_mat$Row.names
  data_mat = subset(data_mat, select = -c(Row.names) )
  colnames(data_mat) <- c("main", colnames(Z) )
  new_data<-data_mat[order( data_mat$main,decreasing = T),]
  new_data <- new_data[new_data$main!=0,]
  new_data<-as.matrix(new_data)
  new_data <- new_data[,-1]
  
  df = data.frame(from = rep(rownames(new_data), times = ncol(new_data)),
                  to = rep(colnames(new_data), each = nrow(new_data)),
                  value = as.vector(new_data),
                  stringsAsFactors = FALSE)
  df<-df[df$value!=0,]
  
  
  chordDiagram(df,big.gap = 10,transparency = 0,scale = F, annotationTrack = "grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df))))))

  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = adj)#,cex = cex)
  }, bg.border = NA) # here set bg.border to NA is important
  
  
}





# get X and Z

GDSC<-readRDS("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/GDSC_data_matteo.rds")

name_drug <- c("Axitinib")

# extract the drugs' pharmacological profiling and tissue dummy
col_filter <- colnames(GDSC$y) %in% paste("IC50.", name_drug,sep="")
YX0 <- cbind(
  GDSC$y[, col_filter],
  GDSC$z[,c(1:12)]
)

X23 <- GDSC$x1
X1 <- log2(X23)

X2=GDSC$x2;Z=GDSC$z[,-c(6)]

rownames(X1)<-GDSC$names.cell.line
rownames(X2)<-GDSC$names.cell.line
rownames(Z)<-GDSC$names.cell.line

X = cbind(X1,X2)





# plotting

# beta_values and theta_values contain the average value
# of beta and theta coefficients from the GDSC analysis

# beta_counts and theta_counts containt the number of times
# a certain coefficients was included in the model

beta = beta_values * (beta_counts>=4)
theta = theta_values* (theta_counts>=4)

plot.circ(beta,theta,X,Z)
