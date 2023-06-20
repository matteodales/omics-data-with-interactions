









# script to create figures 9 and 10 regarding IL.1.R4 and PLXB2













library(RColorBrewer)
library(latex2exp)
library(ggplot2)
library(caret)
library(gridExtra)

load("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/labor_onset_data.rda")

# Function to extract the second character from a string
extract_second_char <- function(string) {
    return(substr(string, 2, 2))
}

# Apply the function to each element of the vector
Timepoint <- sapply(Timepoint, extract_second_char)

y = DOS
X1 = Proteomics
X2 = Metabolomics
Z = cbind(as.integer(Timepoint=='G1'),as.integer(Timepoint=='G3'))

X1_raw = X1

preprocess_X1 = preProcess(X1_raw, method = c("center", "scale"))
X1 = predict(preprocess_X1, X1_raw)

myvar = X1[,'IL.1.R4']

mydata1 = data.frame(myvar,y,Timepoint)

f = function(x) {
  if(x<=1.5){8.76}
  else if(x>1.5 && x<=2.5){7.49}
  else{7.1}
}

ps = seq(0.5, 3.49, length.out = 1001)

mydata2 = data.frame(x = ps, y = sapply(ps,f))

plot_name = paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/IL1R4_plots.pdf")
pdf(file=plot_name, width=18, height=8)
#par(mfrow=c(1,2), las=2, cex.main=2.33, cex.axis = 2.15, cex.lab=2, mar=c(10.3,6.6,6.1,0.9)) #, srt=45)

p1 = ggplot(mydata1, aes(x=myvar, y=y, col=Timepoint)) + xlab('IL.1.R4') + ylab('Days to labor') +
  theme(plot.title = element_text(face='bold', size='22'),legend.position = "none", axis.text=element_text(size=26), axis.title=element_text(size=24), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA),
        plot.margin = margin(1,1,1.5,1.2, "cm")) + geom_point(size=3) +
  geom_smooth(method="lm", se=FALSE, size=3, show.legend=FALSE)

p2 = ggplot(mydata2, aes(x = x, y = y, col=as.character(floor(ps-0.5)+1))) + xlab('Timepoint') + ylab('Coefficient') +
  theme(plot.title = element_text(face='bold', size='22'),legend.text = element_text(size=22), legend.title = element_text(size=24),
        legend.key.size = unit(1, 'cm'), axis.text=element_text(size=26), axis.title=element_text(size=26), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA),
        plot.margin = margin(1,1,1.5,1.2, "cm")) + geom_point(size=5) + labs(colour="Timepoint")

# p2 = ggplot(mydata2, aes(x=tp,y=Coefficient, fill = tp)) + geom_bar(stat='identity', width=0.7) + #geom_line(c(1,2,3),mydata2['Coefficient']) +
#   theme(axis.title=element_text(size=22), axis.text=element_text(size=18),
#         legend.text = element_text(size=16), legend.title = element_text(size=18), legend.key.size = unit(1, 'cm'),
#         plot.margin = margin(1,1,1.5,1.2, "cm")) +
#   xlab('Timepoint') + labs(fill='Timepoint') + geom_line(aes(x = as.numeric(tp), y = Coefficient))

grid.arrange(p1, p2, nrow = 1)

dev.off()
