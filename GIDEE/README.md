# GIDEE
GIDEE requires that you have results from S-PrediXcan with elastic net v8 models ran on all the 49 tissues. Once you have the results, you can run enrichment tests as follows

# Brown's method
```
library(EmpiricalBrownsMethod)
ExpressionFileNames <-  list.files(path = ".", pattern = "*.normalized_expression.csv") #Expression matrices downloaded from GTEx  
MetaXcanFileNames <- list.files(path = ".", pattern = "^results") #S-PrediXcan results for all 49 tissues  
out <- data.frame()  
for(i in 1:49)  
{        
Expression <- read.csv(ExpressionFileNames[i])  
metaXcanTissue <- read.csv(MetaXcanFileNames[i])
colnames(Expression)[4] <- "gene"
a <- match_df(Expression, metaXcanTissue, on = "gene")
pval <- metaXcanTissue$pvalue
BrownInputMatrix <- subset(a, select = -c(X.chr, start, end, gene))
rownames(BrownInputMatrix) <- a$gene
BrownVal <- EmpiricalBrownsMethod::empiricalBrownsMethod(data_matrix = BrownInputMatrix, p_values = pval, extra_info = TRUE)
b <- cbind(ExpressionFileNames[i], nrow(Expression), nrow(metaXcanTissue), nrow(BrownInputMatrix), BrownVal$P_test, BrownVal$P_Fisher, BrownVal$Scale_Factor_C, BrownVal$DF)
out <- rbind(out, b)
print(i)
}
write.csv(out, "BrownsMethod.csv",row.names = F)
```
# Mean squared z-score 
```


