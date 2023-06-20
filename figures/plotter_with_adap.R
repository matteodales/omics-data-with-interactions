library(RColorBrewer)
library(latex2exp)







## Script to generate the results for low and high dimensional simulation studies

## Figures 6 and 7










filename = 'C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/100_ntrain200_ntest9800_pimp30_px1500_px2500_tx10_tx20_sigma24_factstr8_SNR3_with_adap.csv'
data_df <- read.csv(filename)[,-1]

num_sims = 4
filenames = c(
  "C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/100_ntrain200_ntest9800_pimp30_px1500_px2500_tx12_tx22_sigma33_factstr4_SNR1_with_adap.csv",
  "C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/100_ntrain200_ntest9800_pimp30_px1500_px2500_tx14_tx24_sigma40_factstr4_SNR1_with_adap.csv",
  "C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/100_ntrain200_ntest9800_pimp30_px1500_px2500_tx10_tx20_sigma22_factstr4_SNR2_with_adap.csv",
  "C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/100_ntrain200_ntest9800_pimp30_px1500_px2500_tx10_tx20_sigma24_factstr8_SNR3_with_adap.csv"
)



simN = 10
alphalist = c(0,0.2,0.4,0.6,0.8,1,3,5,9)
alpha_axis = sapply(as.character(round(alphalist,1)),
                    function(x) paste0("(alpha=",x,")"))
color = c(rep("white",4), "#B6D7A8","#83B9E3")

plot_name = paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/color_plot_lowdim.pdf")
pdf(file=plot_name, width=18, height=20)
par(mfrow=c(num_sims,4), las=2, cex.main=2.33, cex.axis = 2.15, cex.lab=1.95, mar=c(10.3,6.6,6.1,0.9)) #, srt=45)

for(i in seq(num_sims)){
  
  filename = filenames[i]
  data_df <- read.csv(filename)[,-1]

  p1 = boxplot(data_df[,c(1:6)],
               col=color,
               names=c("Separate X1","Separate X2", "Early fusion", "Late fusion", "     Coop", "Adap Coop"),
               xlab="",ylab="", xaxt = "n")#, main="Test MSE")
  tick = seq_along(p1$names)
  axis(1, at = tick, labels = F)
  text(tick-0.6, par("usr")[3]-0.2*(par("usr")[4]-par("usr")[3]), p1$names, srt = 45, xpd = T, cex=2.0)
  title(main="Test MSE", line=1.8, cex.lab=2.75)

  boxplot(data_df[7:12],
          col=color,
          names=c("Separate X","Separate Z", "Early fusion", "Late fusion", "Coop", "Adap Coop"),
          xlab="", ylab="", main="", xaxt = "n")
  tick = seq_along(p1$names)
  axis(1, at = tick, labels = F)
  text(tick-0.6, par("usr")[3]-0.2*(par("usr")[4]-par("usr")[3]), p1$names, srt = 45, xpd = T, cex=2.0)
  title(main="Number of Features Selected", line=1.8, cex.lab=2.75)

  boxplot(data_df[13:18],
          col=color,
          names=c("Separate X","Separate Z", "Early fusion", "Late fusion", "Coop", "Adap Coop"),
          xlab="", ylab="", main="", xaxt = "n")
  tick = seq_along(p1$names)
  axis(1, at = tick, labels = F)
  text(tick-0.6, par("usr")[3]-0.2*(par("usr")[4]-par("usr")[3]), p1$names, srt = 45, xpd = T, cex=2.0)
  title(main="Number of Interaction Terms", line=1.8, cex.lab=2.75)

  counts = sapply(alphalist, function(x) sum(data_df$chosen_rhos == x))
  counts_adap = sapply(alphalist, function(x) sum(data_df$chosen_rhos_adap == x))
  barplot(rbind(counts,counts_adap),
          beside=TRUE,
          main="",
          names=c(as.character(round(alphalist,1))),
          ylim=c(0,10), cex.lab=2.0,
          density=c(66,66),
          angle=c(36,36),
          col=c("#B6D7A8", "#83B9E3"))
  legend("topleft", c("Coop","Adap Coop"),
         cex=2.33,
         pch=15,
         bty = "n",
         col=c("#B6D7A8", "#83B9E3"))
  title(main=expression(bold(paste("Selected ", rho, " by CV (Counts)"))), line=2.2, cex.lab=2.75)
  mtext(TeX(sprintf("$rho$")), side=1, las=1, line=6.3, cex=1.72)
}

dev.off()
