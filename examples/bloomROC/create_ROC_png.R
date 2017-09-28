# Generates ROC curves.
tags <-  c( "0p0075")
for (FPR_text in tags) {
   bloom <- read.csv(paste0("example_ROC_",FPR_text,"_bloom.csv"))
   FPR = bloom[,4]
   TPR = bloom[,2]
   pdf(paste0("ROC_",FPR_text,"_bloom.png"))
   plot(FPR,TPR, main="ROC curves with option -p 0.0075", 
     xlab="False positive rate",ylab="sensitivity", 
     xlim = c(min(FPR),max(FPR)),
     ylim = c(min(TPR),max(TPR)),
     type="o", col="blue")
   dev.off()
}
