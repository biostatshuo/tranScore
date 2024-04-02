# tranScore

# tranScore：A summary statistics-based gene recognition method, by integrating shared genetic information obtained from the auxiliary population under the transfer learning framework


# Introduction
TLegene is a R procedure for borrowing the idea of transfer learning to integrate useful genetic information available from auxiliary populationfor association analysis of complex diseases and traits. In TLegene, the identification of eGene consists of two components: the first component represents the indirect influence of the auxiliary study 
after transfer learning, and the second component represents the direct effect of the target study.

Specifically, let e be a n by 1 vector of gene expression level on n individuals in the target study, X is a n by p matrix for covariates, G is a n by m matrix for 
genotypes of cis-SNPs for for a given gene in the target study, θ quantifies the association between the gene expression level and the weighted genetic score Gγ which 
is the indirect effect of auxiliary study, b quantifies the association between gene expression level and genotypes G which is the direct effect not completely 
interpreted by auxiliary data and α quantifies a p-vector of fixed effect sizes for clinical covariates. We relate e, Gγ, X and G by a linear mixed model:
<p align="center">
e= Xα + (Gγ) × θ + Gb,  b ~ N(0, τ)
</p>
Above, τ is the genetic variance which is the direct effect not completely interpreted by auxiliary data.

TLegene examines the association of G and Gγ with e (while controlling for X) by testing for:
<p align="center">
H0: θ = 0 and b = 0 <==> H0: θ = 0 and τ = 0
</p>
This is a joint test which requires simultaneously assessing the significance of both fixed effects and random effects: the first part of H0 evaluates the indirect 
influence of auxiliary samples, whereas the second part assesses the direct impact of target samples. Briefly, we derive the test statistic for θ under H0: θ = 0 and τ 
= 0 as usual, while we derive the score statistic for τ under τ = 0 but without the constraint of θ = 0. By doing this, we ensure that these two statistics are 
independent. This strategy substantially eases the development of test statistics for the joint test. In conclusion, under this framework two asymptotically 
independent statistics can be derived. Finally, in order to aggregate the two independent test statistics, we propose three p-value combination approaches (i.e. 
TLegene-oScore, TLegene-aScore, and TLegene-fScore).

# Example
```ruby
library(data.table)
library(harmonicmeanp)
source("TLegene.R")
source("null_model_fit.R")
source("numerical_approximation.R")
data = read.table("data.txt",head=T)
weights = read.table("weights.txt",head=T)
data<-as.matrix(data)
weights<-as.matrix(weights)
p=dim(data)[2]-3
result=TLegene(data,d=2,p,R=1,
               outcome_type="Continuous",
               weight_method= "User",
               user_weight = weights)

po=result$pvalue[3]
pa=result$pvalue[4]
pf=result$pvalue[5]
ph<-cbind(po,pa,pf)
phmp<-as.vector(c(p.hmp(ph,L=length(ph))))
$pvalue

  pvalue.TLegene-oScore   4.93183049954382e-10 

  pvalue.TLegene-aScore   2.86890693215321e-16

  pvalue.TLegene-fScore  2.94312008307861e-15

  pvalue.TLegene-HMP     7.8422646713986e-16
                             
```
  
# Cite
Shuo Zhang<sup>$</sup>, Zhou Jiang<sup>$</sup> and Ping Zeng<sup>#</sup> (2022). Incorporating genetic similarity of auxiliary samples into eGene identification under the transfer learning framework.

# Contact
We are very grateful to any questions, comments, or bugs reports; and please contact Ping Zeng via zpstat@xzhmu.edu.cn.
