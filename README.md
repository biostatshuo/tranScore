# tranScore

# tranScore：A summary statistics-based gene recognition method, by integrating shared genetic information obtained from the auxiliary population under the transfer learning framework


# Backgroubd
tranScore is a R procedure for borrowing the idea of transfer learning to integrate useful genetic information available from auxiliary populationfor association analysis of complex diseases and traits. 
In tranScore, the gene-based association analysis consists of two components: the first component examines the auxiliary population influence, while the second component assesses the target population impact.

With the increase in GWAS sample size, the number of statistically significant loci detected is also increased. However, the gene-based association studies in underrepresented populations are imminent. We developed efficient statistical methods to leverage existing data by incorporating widespread cross-ethnic genetic similarly for the same trait in auxiliary populations with larger sample sizes. However, there are fewer trait studies of underrepresented populations on how to use this shared information more effectively in association analyses.

In this study, we refer to major populations (e.g., EUR) as the auxiliary population and underrepresented populations (e.g., EAS) as the target population. In our application context, transfer learning is employed to enhance the statistical power of association in underrepresented populations. By borrowing this idea, we propose a new statistical method for cross-population association studies and refer to the proposed method as tranScore. We hope to discover loci with a higher probability by exploiting cross-population genetic similarities between the non-EUR and EUR populations.

Unlike previous studies analyzing individual SNPs, tranScore focuses on a set of SNPs within a given gene and evaluates their joint effects on the phenotype of interest. Particularly, the construction of tranScore is based on summary statistics, avoiding the reliance on individual-level data, thus addressing privacy concerns and data accessibility challenges, making it more widely applicable in the post-GWAS era.

# Example
```ruby
library(data.table)
library(harmonicmeanp)
setwd("/public/home/shuozhang/fuct")
source("tranScore.R")
source("ACAT_function.R")

a<-fread("data.txt") 
b<-fread("reference.txt")
estimate=as.numeric(a[,1])
var=as.numeric(a[,2])
weight=as.matrix(as.numeric(a[,3]))
ZAFR=as.numeric(a[,4])
R=1
p=dim(a)[1]

cov.G=cov(b)
cov.G=apply(cov.G,2,as.numeric)
cov.G=apply(cov.G,1,as.numeric)

mis=tranScore(estimate=estimate,var =var,cov=cov.G,weight=weight,R=1,p=p,regularization=FALSE)
po=mis$pvalue[3]
pa=mis$pvalue[4]
pf=mis$pvalue[5]
ph<-cbind(po,pa,pf)
phmp<-as.vector(c(p.hmp(ph,L=length(ph))))
pACAT<-ACAT(ph)

$pvalue

  pvalue.tranScore-optim    0.936118 

  pvalue.tranScore-adapt    0.9149943

  pvalue.tranScore-fisher   0.8964409

  pvalue.tranScore-HMP      0.7904413

  pvalue.tranScore-ACAT     0.919011                        
```
  
# Cite
Shuo Zhang$, Jike Qi$, Yuchen Jiang, Hua Lin, Xinyi Wang, Ting Wang, Hongyan Cao and Ping Zeng# (2025). An integrative association analysis for complex diseases in underrepresented by leveraging the trans-ethnic genetic similarity.

# Contact
We are very grateful to any questions, comments, or bugs reports; and please contact Ping Zeng via zpstat@xzhmu.edu.cn.
