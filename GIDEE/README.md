# GIDEE
GIDEE requires that you have results from S-PrediXcan with elastic net v8 models ran on all the 49 tissues. Once you have the results, you can run enrichment tests as follows

# Brown's method
```
library(EmpiricalBrownsMethod)
ExpressionFileNames <-  list.files(path = ".", pattern = "*.normalized_expression.csv") #Expression matrices downloaded from GTEx  
fileNames <- list.files(path = ".", pattern = "^results") #S-PrediXcan results for all 49 tissues  
medianR2 <- data.frame()
out <- data.frame()
for(i in 1:49){
  Expression <- read.table(ExpressionFileNames[i], header = TRUE)
  file <- read.csv(fileNames[i], header = TRUE)
  me <- median(file$pred_perf_r2)
  mpredperfR2 <- file[which(file$pred_perf_r2 >= me),]
  pval = mpredperfR2$pvalue
  colnames(Expression)[4] <- "gene"
  a <- match_df(Expression,mpredperfR2, on = "gene")
  BrownInputMatrix <- subset(a, select = -c(chr, start, end, gene))
  rownames(BrownInputMatrix) <- a$gene
  BrownVal <-  EmpiricalBrownsMethod::empiricalBrownsMethod(data_matrix = BrownInputMatrix, p_values = pval, extra_info = TRUE)
  #Getting tissue names
  name <- unlist(strsplit(ExpressionFileNames[i], "/"))[5]
  TissueName <- unlist(strsplit(name, ".", fixed = TRUE))[1]
  b <- cbind(TissueName, nrow(Expression), nrow(file), nrow(BrownInputMatrix), BrownVal$P_test, BrownVal$P_Fisher, BrownVal$Scale_Factor_C, BrownVal$DF)
  out <- rbind(out, b)
  medianR2 <- rbind(medianR2, me)
  print(i)
}

colnames(out) <- c("TissueName", "nrow(Expression)", "nrow(metaXcanTissue)", "nrow(BrownInputMatrix)", "BrownVal$P_test", "BrownVal$P_Fisher", "BrownVal$Scale_Factor_C", "BrownVal$DF")
write.csv(out, "Brown-R2.csv", row.names = FALSE)
```
# Mean squared z-score 
```
fileNames <- list.files(pattern = "^results") #S-PrediXcan results for all 49 tissues 
out1 <- data.frame()

for(i in 1:49){
  f <- read.csv(fileNames[i], header = TRUE)
  me <- median(f$pred_perf_r2)
  fR2 <- f[which(f$pred_perf_r2 > me),]
  fR2$zscSqr <- fR2$zscore*fR2$zscore
  meanZsq <- mean(fR2$zscSqr)
  #Getting tissue names
  tn <- unlist(strsplit(fileNames[i], "results-en_"))[2]
  TissueName <- unlist(strsplit(tn, ".", fixed = TRUE))[1]
  out <- cbind(TissueName, meanZsq)
  out1 <- rbind(out1, out)
  print(i)
}
colnames(out1) <- c("TissueName", "meanZsq")
write.csv(out1, "meanZsq-R2.csv", row.names = FALSE)
```
# Getting Residuals for Brown's method and mean squared z-score adjusted with GTEx sample size
```
out1 <- data.frame()
gtxSampleSize <- read.csv("GTEx-Sample-Size.csv", header = TRUE) #file contains the number of genotyped and RNA-seq samples for each of the 49 tissues
out1$SampleSize <- gtxSampleSize$RNASeq.and.Genotyped.samples

brownPval <- read.csv("Brown-R2.csv", header = TRUE)
Zsq <- read.csv("meanZsq-R2.csv", header = TRUE)

out1$BrownPv <- brownPval$BrownVal.P_test
out1$MeanZsq <- Zsq$meanZsq

#calculating residuals for -log10 brown p-value
out1$log10BrownPv <- -log10(out1$BrownPv)
b1 <- lm(out1$log10BrownPv~out1$SampleSize)
out1$brownSS <- b1$residuals

#calculating residuals for mean Z square
out1$meanZsq <- as.numeric(as.character(out1$meanZsq))
t1 <- lm(out1$meanZsq~out1$SampleSize)
out1$ZsqSS <- resid(t1)

write.csv(out1, "GIDEE.csv", row.names = FALSE)
```
