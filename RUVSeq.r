library(RUVSeq)

Step 1: #load the data.
filtered <- read.delim("filt_raw_count_matrix.txt",header=T,row.names=1)

#sample data

#Gene_Symbol     HI102   HI103   HI113   HI115   HI118   HI121   HI122   HI123   HI124   HI127   HI97    HI99    HI102low        HI103low        HI113low        HI115low        HI118low        HI121low        HI122low        HI123low        HI124low        HI127low        HI97low HI99low
#5S_rRNA 3       1       5       1       1       5       9       4       3       9       6       7       9       1       4       7       6       2       1       0       1       5       2       1
#7SK     101     32      84      32      81      129     64      53      27      58      55      95      88      46      45      42      99      43      78      78      67      59      28      54

Step 2: #create design
treat=c(rep("High",12), rep("Low",12))
subjects=factor(c(rep(1:12), rep(1:12)))
design <- model.matrix(~subjects+treat)

Step 3: #Normalise 
set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(design, row.names=colnames(filtered)))
set <- betweenLaneNormalization(set, which="upper")

Step 4: #create empirical data set
y <- DGEList(counts=filtered)
y <- calcNormFactors(y,method="upperquartile")
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
top <- topTags(lrt, n=nrow(y))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

Step 5: #normalise using empirical data set
set2 <- RUVg(set, empirical, k=1)

Step 6: #DE analysis
design <- model.matrix(~subjects+treat+ W_1, data=pData(set2)) 
y <- DGEList(counts=counts(set), group=treat)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design) 
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=3)
topTags(lrt)


<code data-gist-id="5457662" data-gist-line="3,18"></code>
