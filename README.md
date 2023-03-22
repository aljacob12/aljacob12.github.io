# Find a gene project

```{r}
library(bio3d)


#ide.mat <- seqidentity(fasta)

# Plot identity matrix
#plot.dmat(ide.mat, color.palette=mono.colors,
         # main="Sequence Identity", xlab="Structure No.",
          #ylab="Structure No.")
```


```{r}
library(bio3d)

fast <- read.fasta("multiple_alignment.fa")
```
```{r}
ide.mat <- seqidentity(fast)

# Plot identity matrix
heatmap(ide.mat , margins = c(8,8),cexRow = 0.7, cexCol = 0.7, png("S7.png", width=800, height=800))
          #resnum.1 = fast[ide.mat],
          #main="Sequence Identity", xlab="Structure No.",
          #ylab="Structure No.")

```
```{r}

```
```{r}
"install.packages("ggmsa")

library(Biostrings)
library(ggmsa)

protein_sequences <- system.file("extdata", "multiple_alignment", package = "ggmsa")

aln = readAAMultipleAlignment(protein_sequences)

ggmsa(protein_sequences, start = 265, end = 300)


aln = unmasked(aln)
names(aln)[1]


ref = aln[1]

bm = sapply(1:length(aln),function(i){
  as.numeric(as.matrix(aln[i])==as.matrix(ref))
})

bm = t(bm)
rownames(bm) = names(aln)

library(pheatmap)
pheatmap(bm[nrow(bm):1,265:300],cluster_rows=FALSE,cluster_cols=FALSE)

```

# Halloween project

```{r}
candy_file <- "candy-data.csv"

```

```{r}
candy = read.csv("candy-data.csv", row.names=1)
head(candy)
```

#Q1. How many different candy types are in this dataset?
```{r}
dim(candy)
```
12 types of candy

#Q2. How many fruity candy types are in the dataset?
```{r}
sum(candy$fruity)
```
38 fruit types of candy


##2. What is your favorate candy?
Twix

#Q3. What is your favorite candy in the dataset and what is it’s winpercent value?
```{r}
candy["Twix", ]$winpercent
```

#Q4. What is the winpercent value for “Kit Kat”?
```{r}
candy["Kit Kat", ]$winpercent
```

#Q5. What is the winpercent value for “Tootsie Roll Snack Bars”?
```{r}
candy["Tootsie Roll Snack Bars", ]$winpercent
```

```{r}
library("skimr")
skim(candy)
```




#Q6. Is there any variable/column that looks to be on a different scale to the majority of the other columns in the dataset?

The hist column is a graph rather than a numerical readout
    
#Q7. What do you think a zero and one represent for the candy$chocolate column?
Zeros would mean that the candy in question is not chocolate while a 1 means that it is chocolate 


##Graph

#Q8. Plot a histogram of winpercent values
```{r}
hist(candy$winpercent)
```

#Q9. Is the distribution of winpercent values symmetrical?
no it is not
#Q10. Is the center of the distribution above or below 50%?
it is below 50%
#Q11. On average is chocolate candy higher or lower ranked than fruit candy?
```{r}
choco_win = candy$winpercent[as.logical(candy$chocolate)]
mean(choco_win)
```

```{r}
candy_win = candy$winpercent[as.logical(candy$fruity)]
mean(candy_win)
```
On avg chocolate is higher ranked than fruit candy
#Q12. Is this difference statistically significant?
```{r}
t.test(choco_win,candy_win)
```
yes it is statistically significant

##Overall Candy Rankings

#Q13. What are the five least liked candy types in this set?
```{r}
head(candy[order(candy$winpercent),], n=5)
```

#Q14. What are the top 5 all time favorite candy types out of this set?

```{r}
head(candy[order(candy$winpercent, na.last = TRUE, decreasing = TRUE),], n=5)
```
#Graph:
```{r}
library(ggplot2)

ggplot(candy) + 
  aes(x=winpercent, y=rownames(candy)) +
  geom_col()
```
#organized graph:
```{r}
library(ggplot2)

my_cols=rep("black", nrow(candy))
my_cols[as.logical(candy$chocolate)] = "chocolate"
my_cols[as.logical(candy$bar)] = "brown"
my_cols[as.logical(candy$fruity)] = "pink"

ggplot(candy) +
  aes(winpercent, reorder(rownames(candy),winpercent))+
  geom_col(fill=my_cols)
```
#Q17. What is the worst ranked chocolate candy?
The worst ranked chocolate candy is the charleston chew

#Q18. What is the best ranked fruity candy?
The best ranked fruity candy is starburt


##Taking a look at pricepercent

```{r}
library(ggrepel)

# How about a plot of price vs win
ggplot(candy) +
  aes(winpercent, pricepercent, label=rownames(candy)) +
  geom_point(col=my_cols) + 
  geom_text_repel(col=my_cols, size=3.3, max.overlaps = 5)

```

#Q19. Which candy type is the highest ranked in terms of winpercent for the least money - i.e. offers the most bang for your buck?
Reese's Miniatures

#Q20. What are the top 5 most expensive candy types in the dataset and of these which is the least popular?
```{r}
ord <- order(candy$pricepercent, decreasing = TRUE)
head( candy[ord,c(11,12)], n=5 )
```
The top 5 most expensive candies are Nik L Nip , Nestle Smarties, Ring pop, Hershey's Krackel, Hershey's Milk Chocolate. The least popular is Nik L Nip.

##Exploring the correlation structure

```{r}
cij <- cor(candy)

```


#Q22. Examining this plot what two variables are anti-correlated (i.e. have minus values)?
chocolate and fruity are strongly anti-correlated

#Q23. Similarly, what two variables are most positively correlated?
besides the same variables being correlated with themselves, winpercent and chocolate and strongly positively correlated, and bar and chocolate are similarly strongly correlated.


##6. Principal Component Analysis

```{r}
pca <- prcomp(candy, scale=TRUE)
summary(pca)

plot(pca$x[,1:2], col=my_cols, pch=16)
my_data <- cbind(candy, pca$x[,1:3])
p <- ggplot(my_data) + 
        aes(x=PC1, y=PC2, 
            size=winpercent/100,  
            text=rownames(my_data),
            label=rownames(my_data)) +
        geom_point(col=my_cols)

p
```
```{r}
library(ggrepel)

p + geom_text_repel(size=3.3, col=my_cols, max.overlaps = 7)  + 
  theme(legend.position = "none") +
  labs(title="Halloween Candy PCA Space",
       subtitle="Colored by type: chocolate bar (dark brown), chocolate other (light brown), fruity (red), other (black)",
       caption="Data from 538")
```
```{r}
par(mar=c(8,4,2,2))
barplot(pca$rotation[,1], las=2, ylab="PC1 Contribution")
```

#Q24. What original variables are picked up strongly by PC1 in the positive direction? Do these make sense to you?

Fruity, hard, and pluribus variables are picked up strongly by PC1 in the positive direction, which makes sense because these types of candies are less popular than the average and would be "overestimated" in PC1 plot. 

# Class 12
##Bioconductor and DESeq2 setup
```{r}
library(BiocManager)
library(DESeq2)
```

##Import countData and colData
```{r}
counts <- read.csv("airway_scaledcounts.csv", row.names=1)
metadata <-  read.csv("airway_metadata.csv")
head(counts)
head(metadata)
```

>Q1. How many genes are in this dataset?
>Q2. How many ‘control’ cell lines do we have?


##Toy differential gene expression
```{r}
control <- metadata[metadata[,"dex"]=="control",]
control.counts <- counts[ ,control$id]
control.mean <- rowSums( control.counts )/4 
head(control.mean)
```

>Q3. How would you make the above code in either approach more robust?
>Q4. Follow the same procedure for the treated samples (i.e. calculate the mean per gene across drug treated samples and assign to a labeled vector called treated.mean)


```{r}
library(dplyr)
control <- metadata %>% filter(dex=="control")
control.counts <- counts %>% select(control$id) 
control.mean <- rowSums(control.counts)/4
head(control.mean)
```

```{r}
treated <- metadata[metadata[,"dex"]=="treated",]
treated.mean <- rowSums( counts[ ,treated$id] )/4 
names(treated.mean) <- counts$ensgene
```

```{r}
meancounts <- data.frame(control.mean, treated.mean)
```

>Q5 (a). Create a scatter plot showing the mean of the treated samples against the mean of the control samples. Your plot should look something like the following
>Q5 (b).You could also use the ggplot2 package to make this figure producing the plot below. What geom_?() function would you use for this plot

```{r}
plot(meancounts[,1],meancounts[,2], xlab="Control", ylab="Treated")
```

```{r}
meancounts$log2fc <- log2(meancounts[,"treated.mean"]/meancounts[,"control.mean"])
head(meancounts)
```


```{r}
zero.vals <- which(meancounts[,1:2]==0, arr.ind=TRUE)

to.rm <- unique(zero.vals[,1])
mycounts <- meancounts[-to.rm,]
head(mycounts)
```


```{r}
up.ind <- mycounts$log2fc > 2
down.ind <- mycounts$log2fc < (-2)
```


##DESeq2 analysis


## Gene annotation

Use one of bioconductors's main annotation packages to help with mapping between various ID schemes.
```{r}
head(res)
row.names(res)
```


```{r}
library(AnnotationDbi)
library(org.Hs.eg.db)
```


Look at what types of IDs I can translate between from the 'org.Hs.eg.b'
```{r}
columns(org.Hs.eg.db)
```


```{r}
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),  # Our genenames
                     keytype="ENSEMBL",        # The format of our genenames
                     column="SYMBOL",          # The new format we want to add
                     multiVals="first")
```


```{r}
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

res$uniprot <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="UNIPROT",
                     keytype="ENSEMBL",
                     multiVals="first")

res$genename <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

head(res)
```

## Pathway analysis

Here we play with just one, the GAGE package (which stands for Generally Applicable Gene set Enrichment), to do KEGG pathway enrichment analysis on our RNA-seq based differential expression results.


```{r}
library(pathview)
library(gage)
library(gageData)

data(kegg.sets.hs)
head(kegg.sets.hs, 2)
```


The main gage() function requires a named vector of fold changes, where the names of the values are the Entrez gene IDs.

```{r}
foldchanges <- res$log2FoldChange
names(foldchanges) <- res$entrez
head(foldchanges)
```

Now, let’s run the gage pathway analysis.

```{r}
keggres = gage(foldchanges, gsets=kegg.sets.hs)
```

Now lets look at the object returned from gage(). our results here:

```{r}
attributes(keggres)
```


```{r}
# Look at the first three down (less) pathways
head(keggres$less, 3)
```


Lets pull up the highlighted pathways and show our differentially expressed genes on the pathway. 

```{r}
pathview(gene.data=foldchanges, pathway.id="hsa05310")
```


put this into my document
![The asthma pathway with my highlighted differentially expressed genes in color](hsa05310.pathview.png)
