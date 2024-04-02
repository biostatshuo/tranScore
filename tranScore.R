#This function performs using summary statistics
#@param estimate, marginal point estimate from target population GWAS
#@param var, marginal variance estimate from target population GWAS
#@param cov, covariance of genotype from target population
#@param weight, marginal point estimate from auxiliary population GWAS
#@param R, Input: Number of burden scores
#@param p, Input: Number of genetic variants in the set.
#@param regularization=TRUE, A logical value (default: TRUE) specifying if regularization on covariance needs to be carried out
#@param n, sample size used to estimate covariance matrix
#@param burden_scale=TRUE, Input: A logical value (default: TRUE) specifying if standardization on burden components need to be carried out.
#@param burden_name=paste0(rep("burden",R),seq(1:R)), The names of the burden components
#@param exclusive_burden=rep(0, R), 1_r=exclude the rth burden component component for joint burden testing, 0 otherwise.Used for adjusting for known loci. In normal case, set= rep(0, R)
#@param exclusive_variance=rep(0, p), 1_p=exclude the pth variant for variance component testing, 0 otherwise. Used for adjusting for known loci. In normal case, set= rep(0, p)
#@param chisq_app = "3M", Input: moment matching (Liu's) method in quantile approximation in optimal linear combination method.
#@param combinations_return=TRUE, Input: indicating if the combination methods, including the optimal linear combination, data-adpative weighted combination and Fisher's combination method, are implemented. By assigning FALSE, only the burden and variance component tests are conducted.
#@param combination_preference = "All", Input:other options: "OptMin", "AdptWt", "Fisher"
#@param acc = 5e-10, Input:A numerical value indicating the precision of Davies method for p-value calculation. Default is 5e-10.
#@param acc_auto = TRUE, Input: if "TRUE": acc is determined by the magnitute of min_pvalue in OptMin
#@param accurate_app_threshold = -log10(0.05), Input: if the min_pvalue<0.05, the accurate approximation on min_quantile in OptMin combination is used instead of approximation with Liu's method.
#@param max_core = 4 ), Input: An integer specifying the maximum number of cores can be recruited in parallel package. Default is 4 cores.)
#' @return A list of pvalues and statistics
#' @export
tranScore=function(estimate, # marginal point estimate from target population GWAS
                        var, # marginal variance estimate from target population GWAS
                        cov, # covariance of genotype from target population
                        weight, # marginal point estimate from auxiliary population GWAS
                        R, # Input: Number of burden scores
                        p,# Input: Number of genetic variants in the set.
                        regularization=TRUE, #A logical value (default: FALSE) specifying if regularization on covariance needs to be carried out
                        n, #sample size used to estimate covariance matrix
                        burden_scale=TRUE,  # Input: A logical value (default: TRUE) specifying if standardization on burden components need to be carried out.
                        burden_name=paste0(rep("burden",R),seq(1:R)), # The names of the burden components
                        exclusive_burden=rep(0, R), # 1_r=exclude the rth burden component component for joint burden testing, 0 otherwise.Used for adjusting for known loci. In normal case, set= rep(0, R)
                        exclusive_variance=rep(0, p), # 1_p=exclude the pth variant for variance component testing, 0 otherwise. Used for adjusting for known loci. In normal case, set= rep(0, p)
                        chisq_app = "3M",# Input: moment matching (Liu's) method in quantile approximation in optimal linear combination method.
                        combinations_return=TRUE, # Input: indicating if the combination methods, including the optimal linear combination, data-adpative weighted combination and Fisher's combination method, are implemented. By assigning FALSE, only the burden and variance component tests are conducted.
                        combination_preference = "All", # Input:other options: "OptMin", "AdptWt", "Fisher"
                        acc = 5e-10, # Input:A numerical value indicating the precision of Davies method for p-value calculation. Default is 5e-10.
                        acc_auto = TRUE, # Input: if "TRUE": acc is determined by the magnitute of min_pvalue in OptMin
                        accurate_app_threshold = -log10(0.05), # Input: if the min_pvalue<0.05, the accurate approximation on min_quantile in OptMin combination is used instead of approximation with Liu's method.
                        max_core = 4 )# Input: An integer specifying the maximum number of cores can be recruited in parallel package. Default is 4 cores.)
{

  if(regularization==TRUE){
   cov_name=rownames(cov)
   cor=cov2cor(as.matrix(cov))
   tune=1/sqrt(n)
  cor[row(cor)!=col(cor)]=cor[row(cor)!=col(cor)]/(1+tune)
  cov=diag((diag(cov))^0.5)%*%cor%*%diag((diag(cov))^0.5)
  rownames(cov)=cov_name}else{
    cov=cov
  }

  if (!is.matrix(weight)) {weight <-as.matrix(weight)}
  if(burden_scale==TRUE){
    for(i in 1:ncol(weight)){
    scale_factor=(t(weight[,i])%*%cov%*%weight[,i])^0.5
    weight[,i]=weight[,i]/rep(scale_factor, length(weight[,i]))}

  }else{
    weight=weight
  }



 #numerical approximation functions
  WeightSumChisq = function(weight,df=rep(1,length(weight))){
    # Calculate the df and the noncentral parameter of an approximated noncentral chisq distribution for the weighted sum of independent chisq distributions.
    # Option weight: a vector of weights for each chisq.
    # Option df: a vector of df for each chisq.
    # Return the df(l) and the noncentrality (\delta) of the approximation noncentral chisq.
    # Reference: Liu (2009). Match the skewness.

    # Sum of weights to different powers and terms related to criterions
    S1 = sum(weight)    # equal to c1[1]
    S2 = sum(weight^2)  # equal to c1[2]
    S3 = sum(weight^3)
    S4 = sum(weight^4)
    tmp_1 = S4/(S2^2)
    tmp_2 = (S3^2)/(S2^3)

    # Calculate a=sqrt(l+2\delta), \delta, and l
    if(tmp_1<tmp_2){
      a = 1/(S3/(S2^(3/2))-sqrt(tmp_2-tmp_1))
      delta = a^3*S3/(S2^(3/2))-a^2
      l = a^2-2*delta
    }
    else{
      a = (S2^(3/2))/S3
      delta = 0
      l = a^2
    }
    #print(c(tmp_1,tmp_2,a,delta,l))
    return(c(mean.Q=S1,SD.Q=sqrt(2*S2),df=l,noncentral=delta,mean.X=l+delta,SD.X=sqrt(2*(l+2*delta))))
  }

  quan_rho_calculator = function(grid=seq(0,1,by=0.05),df_f,lambda_r,min_pvalue,acc=5e-10,lim=1e4,max_core=4){
    grid_inner = grid[2:(length(grid)-1)]
    numCores = parallel::detectCores()

    quan_rho_grid_inner = mclapply(grid_inner,FUN=function(grid_point){
      diff.prob.davies = function(q,tailprob,lambda,df=rep(1,length(lambda)),delta=rep(0,length(lambda)),lim=1e4,acc=5e-10){
        tailprob-CompQuadForm::davies(q,lambda=lambda,lim=lim,acc=acc)$Qq
      }
      quan_tmp = uniroot(diff.prob.davies,lower=0,upper=1e7,tailprob=min_pvalue,lambda=c(grid_point*rep(1,df_f),(1-grid_point)*lambda_r),acc=acc,lim=lim)$root
    },
    mc.cores=min(numCores,max_core))
    quan_rho_grid_inner = as.numeric(unlist(quan_rho_grid_inner))
    names(quan_rho_grid_inner) = as.character(grid_inner)
    return(list(grid=grid,quan_rho_grid_inner=quan_rho_grid_inner))
  }

  pvalue_minimizePValue_davies = function(u_f,u_r,df_f,lambda_r,min_pvalue=NULL,quan_grid=NULL,acc=5e-10,lim=1e4,ChisqApp="3M",max_core=4){
    if(is.null(min_pvalue)){
      p_rho = function(rho,df_f,lambda_r,u_f,u_r,acc=5e-10){
        pvalue = CompQuadForm::davies(q=rho*u_f+(1-rho)*u_r,
                                      lambda=c(rho*rep(1,df_f),(1-rho)*lambda_r),acc=acc,lim=1/acc)$Qq
        if(pvalue<=acc) pvalue = pvalue.liu(Q=rho*u_f+(1-rho)*u_r,
                                            weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),
                                            method=ChisqApp)$pvalue
        return(pvalue)
      }
      min_pvalue_inside01 = optimize(p_rho,interval=c(0,1),maximum=FALSE,df_f=df_f,lambda_r=lambda_r,u_f=u_f,u_r=u_r,acc=acc)$objective
      min_pvalue = min(c(min_pvalue_inside01,
                         p_rho(0,df_f=df_f,lambda_r=lambda_r,u_f=u_f,u_r=u_r,acc=acc),
                         p_rho(1,df_f=df_f,lambda_r=lambda_r,u_f=u_f,u_r=u_r,acc=acc)))
    }


    integrand = function(a,min_pvalue,df_f,lambda_r,acc=5e-10,grid=NULL,quan_rho_grid_inner=NULL,max_core=4){
      quan_rho = function(rho,df_f,lambda_r,p_rho_min,u_f){
        if(rho==1) stop("rho is not allowed to be 1.")
        else{
          diff.prob.davies = function(q,tailprob,lambda,df=rep(1,length(lambda)),delta=rep(0,length(lambda)),lim=1e4,acc=5e-10){
            tailprob-CompQuadForm::davies(q,lambda=lambda,lim=lim,acc=acc)$Qq
          }
          quan_tmp = uniroot(diff.prob.davies,lower=0,upper=1e7,tailprob=p_rho_min,lambda=c(rho*rep(1,df_f),(1-rho)*lambda_r),acc=acc,lim=lim)$root
          names(quan_tmp) = c()
          quan_tmp_u_r = (quan_tmp - rho*u_f)/(1-rho)
          return(quan_tmp_u_r)
        }
      }

      optimize_interval = c(0,1)

      numCores = parallel::detectCores()
      res = mclapply(a,FUN=function(a.tmp){
        if(!is.null(grid) & !is.null(quan_rho_grid_inner)){
          grid_inner = as.numeric(names(quan_rho_grid_inner))
          ### refne interval for optimization by running the following over a grid (no 0 and 1)
          grid_inner_index = order((quan_rho_grid_inner - grid_inner*a.tmp)/(1-grid_inner))[1]
          optimize_interval = c(grid[grid_inner_index],grid[grid_inner_index+2])
        }
        min_quantile = optimize(quan_rho,interval=optimize_interval,df_f=df_f,lambda_r=lambda_r,p_rho_min=min_pvalue,u_f=a.tmp)$objective
        tail.prob = 0
        if(class(try(CompQuadForm::davies(q=min_quantile,lambda=lambda_r,acc=acc,lim=1/acc)))!="try-error") tail.prob = CompQuadForm::davies(q=min_quantile,lambda=lambda_r,acc=acc,lim=1/acc)$Qq
        if(tail.prob<=acc) tail.prob = pvalue.liu(Q=min_quantile,weight=lambda_r,method=ChisqApp)$pvalue
        (1-tail.prob) * dchisq(a.tmp,df=df_f)
      },
      mc.cores=min(numCores,max_core))
      res = as.numeric(unlist(res))
      return(res)
    }

    pvalue = 1 - integrate(integrand,lower=0,upper=qchisq(min_pvalue,df=df_f,lower.tail=FALSE),min_pvalue=min_pvalue,df_f=df_f,lambda_r=lambda_r,acc=acc,grid=quan_grid[[1]],quan_rho_grid_inner=quan_grid[[2]],max_core=max_core)$value
    if(pvalue<0) pvalue = abs(pvalue)
    return(pvalue)
  }

  pvalue_minimizePValue_liu = function(u_f,u_r,df_f,lambda_r,min_pvalue=NULL,acc=5e-10,lim=1e4,ChisqApp="4M"){
    if(is.null(min_pvalue)){
      p_rho = function(rho,df_f,lambda_r,u_f,u_r,acc=5e-10){
        pvalue = CompQuadForm::davies(q=rho*u_f+(1-rho)*u_r,
                                      lambda=c(rho*rep(1,df_f),(1-rho)*lambda_r),acc=acc,lim=1/acc)$Qq
        if(pvalue<=acc) pvalue = pvalue.liu(Q=rho*u_f+(1-rho)*u_r,
                                            weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),
                                            method=ChisqApp)$pvalue
      }
      min_pvalue_inside01 = optimize(p_rho,interval=c(0,1),maximum=FALSE,df_f=df_f,lambda_r=lambda_r,u_f=u_f,u_r=u_r,acc=acc)$objective
      min_pvalue = min(c(min_pvalue_inside01,pvalue.r,pvalue.f))
    }

    integrand = function(a,min_pvalue,df_f,lambda_r,acc=5e-10){
      quan_rho = function(rho,df_f,lambda_r,p_rho_min,u_f){
        if(rho==1) stop("rho is not allowed to be 1.")
        else{
          if(ChisqApp=="4M") parameters = WeightSumChisq_mod(weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),df=rep(1,length(weight)))
          if(ChisqApp=="3M") parameters = WeightSumChisq(weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),df=rep(1,length(weight)))
          # mean.Q, SD.Q, df, noncentral, mean.X, SD.X
          quan_tmp = (qchisq(p_rho_min,df=parameters["df"],ncp=parameters["noncentral"],lower.tail=FALSE) - parameters["mean.X"]) / parameters["SD.X"] * parameters["SD.Q"] + parameters["mean.Q"] # Note: use upper tail for quantile calcualtion to avoid Inf when p_rho_min is extremely small.
          names(quan_tmp) = c()
          quan_tmp_u_r = (quan_tmp - rho*u_f)/(1-rho)
          return(quan_tmp_u_r)
        }
      }
      res = sapply(a,function(a.tmp){
        # min_quantile = optimize(quan_rho,interval=c(0,1-1e-3),df_f=df_f,lambda_r=lambda_r,p_rho_min=min_pvalue,u_f=a.tmp)$objective
        min_quantile_tmp = optimize(quan_rho,interval=c(0,1),df_f=df_f,lambda_r=lambda_r,p_rho_min=min_pvalue,u_f=a.tmp)$objective # try 2017-04-06
        min_quantile = min(quan_rho(rho=0,df_f=df_f,lambda_r=lambda_r,p_rho_min=min_pvalue,u_f=a.tmp),
                           min_quantile_tmp) # try 2017-04-06
        tail.prob = 0
        if(class(try(CompQuadForm::davies(q=min_quantile,lambda=lambda_r,acc=acc,lim=1/acc)))!="try-error") tail.prob = CompQuadForm::davies(q=min_quantile,lambda=lambda_r,acc=acc,lim=1/acc)$Qq
        if(tail.prob<=acc) tail.prob = pvalue.liu(Q=min_quantile,weight=lambda_r,method=ChisqApp)$pvalue
        (1-tail.prob) * dchisq(a.tmp,df=df_f)
      })
      return(res)
    }
    pvalue = 1 - integrate(integrand,lower=0,upper=qchisq(min_pvalue,df=df_f,lower.tail=FALSE),min_pvalue=min_pvalue,df_f=df_f,lambda_r=lambda_r,acc=acc)$value
    return(pvalue)
  }

  pvalue.liu = function(Q,weight,method="3M",df=rep(1,length(weight))){
    if(method=="4M") para = WeightSumChisq_mod(weight=weight,df=df)
    else para = WeightSumChisq(weight=weight,df=df)
    mean.Q = para[1]
    SD.Q = para[2]
    mean.X = para[5]
    SD.X = para[6]
    df.app = para[3]
    noncentral.app = para[4]
    Q.new = (Q-mean.Q)/SD.Q * SD.X + mean.X
    pvalue = pchisq(Q.new,df=df.app,ncp=noncentral.app,lower.tail=FALSE)
    return(list(pvalue=unname(pvalue),df=unname(df.app),ncp=unname(noncentral.app)))
  }


  weight=as.matrix(weight)

  #point estimate of burden
  dd = diag(cov)
  Sigma1 = solve(t(weight)%*%cov%*%weight)%*%t(weight)%*%diag(dd)
  betaB = Sigma1%*%estimate
  betaB
  point=betaB
  cov2=cov2cor(cov)
  d=diag(var^0.5)
  covvv=d%*%cov2%*%d

  #variance of burden estimate
  var2=solve(t(weight)%*%cov%*%weight)%*%t(weight)%*%diag(dd)%*%covvv%*%t(solve(t(weight)%*%cov%*%weight)%*%t(weight)%*%diag(dd) )

# glm(y~as.matrix(design)%*%as.matrix(weight)+as.matrix(confounder), family="binomial")



  #pvalues of burden
  #indicator_burden=c(1,0)
  indicator_burden=exclusive_burden
  idx=which(indicator_burden==0)
  burden_stat=t(betaB[idx])%*%solve(var2[idx, idx])%*%betaB[idx]
  burden_stat_ind=betaB/(diag(var2)^0.5)
  z_stats=burden_stat
  pvalue.f.sum_ind=2 * pnorm(abs(burden_stat_ind), lower.tail = FALSE)
  pvalue.f.sum=pchisq(z_stats, df=length(idx), lower.tail = FALSE)
  SE=(diag(var2)^0.5)
  score=z_stats/(((diag(var2))^0.5)[1])
    #z_stats/((diag(var2))[1])^0.5
 # z=estimate/var^0.5
  #u=t(weight[,1])%*%(z/var^0.5)
  #score=u
  #variance component
  cov1 = matrix(0,R+1,R+1)
  estimate.adj.summary = vector()
  exclusive_indicator=exclusive_variance

  for (i in 1:p){
    #weight_combined=cbind(weight, cov[i,])
    for(j in 1:(R)){for( k in 1:(R)){
      cov1[j,k]=t(weight[,j])%*%cov%*%weight[,k]
    }}
    cov1[R+1,]=c(cov[i,]%*%weight, cov[i,i])
    cov1[,R+1]=c(cov[i,]%*%weight, cov[i,i])
    vec=rep(0, p)
    vec[i]=1
    tt = try(solve(cov1),silent = TRUE)
    if (class(tt)=="try-error"){
      estimate.adj.summary[i]<-NA
      exclusive_indicator[i]<-1
    } else{
    estimate.adj.summary[i] = t(c(rep(0, R), 1))%*%solve(cov1)%*%rbind(t(weight),vec)%*%diag(dd)%*%estimate
  }}
  cov1 = matrix(0,R+1,R+1)
  var.alpha=vector()
  for (i in 1:p){
   # for(j in 1:(R)){for( k in 1:(R)){
   #   cov1[j,k]=t(weight[,j])%*%cov%*%weight[,k]
   # }}
   cov1[1:R,1:R]<-t(weight)%*%cov%*%weight
    cov1[R+1,]=c(cov[i,]%*%weight, cov[i,i])
    cov1[,R+1]=c(cov[i,]%*%weight, cov[i,i])
    vec=rep(0, p)
    vec[i]=1
    tt = try(solve(cov1),silent = TRUE)
    if (class(tt)=="try-error"){
      var.alpha[i]<-NA
      exclusive_indicator[i]<-1
    } else{
    var.alpha[i] =  t(c(rep(0, R), 1))%*%solve(cov1)%*%rbind(t(weight),vec)%*%diag(dd)%*%covvv%*%t(t(c(rep(0, R), 1))%*%solve(cov1)%*%rbind(t(weight),vec)%*%diag(dd))
  }}

  cov1 = matrix(0,R+1,R+1)
  cov_adj = matrix(0, nrow=p, ncol=p)
  for (i in 1:p){

    cov1[1:R,1:R]<-t(weight)%*%cov%*%weight

#    for(j in 1:(R)){for( k in 1:(R)){
 #     cov1[j,k]=t(weight[,j])%*%cov%*%weight[,k]
   # }}
    cov1[R+1,]=c(cov[i,]%*%weight, cov[i,i])
    cov1[,R+1]=c(cov[i,]%*%weight, cov[i,i])







    vec=rep(0, p)
    vec[i]=1
    tt = try(solve(cov1),silent = TRUE)
    if (class(tt)=="try-error"){
      cov_adj[i,]<-rep(NA,p)
      exclusive_indicator[i]<-1
    } else{
    cov_adj[i,] =t(c(rep(0, R), 1))%*%solve(cov1)%*%rbind(t(weight),vec)%*%diag(dd)
  }}
  idx=which(is.na(estimate.adj.summary)==TRUE)
  #cov_adj=cov_adj[-idx,-idx]
 # covvv2=covvv[-idx, -idx]
 cor.adj.summary=cov_adj%*%covvv%*%t(cov_adj)
  #cor.adj.summary=cov2cor(cor.adj.summary)

  exclusive_indicator[idx]=1
  ind=exclusive_indicator
  num=which(ind==0)
  estimate.adj.summary=estimate.adj.summary[num]
  var.alpha=var.alpha[num]
  U=(estimate.adj.summary/var.alpha)
  stat = sum((estimate.adj.summary/var.alpha)^2)
  random_stat=stat
  cor.adj.summary=cor.adj.summary[num, num]
  cor.adj.summary=cov2cor(cor.adj.summary)
  V = diag(1/var.alpha^0.5)%*%cor.adj.summary%*%diag(1/var.alpha^0.5)

  lambda=eigen(V)$values

  #dist5=davies(q=stat,lambda=lambda,acc=acc, lim=1/acc)$Qq
  dist4=CompQuadForm::davies(q=stat,lambda=lambda,acc=acc, lim=1/acc)$Qq
  if(dist4<=acc){dist4= pvalue.liu(Q=stat,weight=lambda,method=chisq_app)$pvalue
  }
  result=list(z_stats, pvalue.f.sum, U, V, dist4)
  names(result)=c("z score of burden","burden pvalue", "U", "V","variance pvalue")
#  if(combinations_return==FALSE){
  #  return(c(pvalue.f.sum, pvalue.f.sum_ind,dist4))
 # }else{

  #combined test
  pvalue = c(pvalue.f=NA,pvalue.r=NA,pvalue.otranScore=NA,pvalue.atranScore=NA,pvalue.ftranScore=NA)
  stat = c(stat.f=NA,stat.r=NA,stat.otranScore=NA,stat.atranScore=NA,stat.ftranScore=NA)
  pvalue.f.ind = NA
  stat.f.ind = NA
  rho = c(rho.otranScore.f=NA,rho.otranScore.r=NA,rho.atranScore.f=NA,rho.atranScore.r = NA)

  #Q.f=burden_stat_ind[1]
  Q.f=z_stats
  Q.r=random_stat
  EigenvalueD.r=lambda
  #S.f = as.matrix(burden)
  pvalue.f=pvalue.f.sum
  pvalue.r=dist4
  pvalue["pvalue.r"] = pvalue.r
  stat["stat.r"] = Q.r
  pvalue["pvalue.f"] = pvalue.f
  stat["stat.f"] = Q.f

  pvalue.f.ind = pvalue.f.sum_ind
  stat.f.ind=burden_stat_ind^2


  if(combinations_return & (("All" %in% combination_preference) | ("OptMin" %in% combination_preference))){
    # OptMin method - Calculate the optimal weight for linear combinations of Q.f and Q.r by minimizing the pvalues
    p_rho = function(rho,df_f,lambda_r,u_f,u_r,acc=5e-10){
      pvalue = CompQuadForm::davies(q=rho*u_f+(1-rho)*u_r,
                                    lambda=c(rho*rep(1,df_f),(1-rho)*lambda_r),acc=acc,lim=1/acc)$Qq
      # Note on 2017/2/6: previous version doesn't include acc and lim input here and it causes the returned pvalue >1 in some simulation runs. Caused error in the following calculation of quantiles based on returned pvalues. acc and lim are added to control this.
      if(pvalue<=acc) pvalue = pvalue.liu(Q=rho*u_f+(1-rho)*u_r,weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),method=chisq_app)$pvalue
      return(pvalue)
    }
    rho_min_pvalue_inside01 = optimize(p_rho,interval=c(0,1),maximum=FALSE,df_f=R,lambda_r=EigenvalueD.r,u_f=Q.f,u_r=Q.r,acc=acc)$minimum
    min_pvalue_inside01 = optimize(p_rho,interval=c(0,1),maximum=FALSE,df_f=R,lambda_r=EigenvalueD.r,u_f=Q.f,u_r=Q.r,acc=acc)$objective
    ### Note: optimize function does NOT include the comparison at the endpoints of the interval. Add comparison with p_values at 0 and 1 manually below
    if(min_pvalue_inside01>pvalue.r | min_pvalue_inside01>pvalue.f){
      rho_min_pvalue = ifelse(pvalue.r<pvalue.f,0,1)
      min_pvalue = ifelse(pvalue.r<pvalue.f,pvalue.r,pvalue.f)
    }
    else{
      rho_min_pvalue = rho_min_pvalue_inside01
      min_pvalue = min_pvalue_inside01
    }

    Q.OptMin = Q.f*rho_min_pvalue + Q.r*(1-rho_min_pvalue)

    if(acc_auto==TRUE){
      acc_supp = max(min(min_pvalue,0.01)*0.1,acc)
      lim_supp = as.integer(1/acc_supp)
    }
    else{
      acc_supp = acc
      lim_supp = as.integer(1/acc)
    }

    if(-log10(min_pvalue)>accurate_app_threshold){
      quan_grid = quan_rho_calculator(grid=seq(0,1,by=0.05),df_f=R,lambda_r=EigenvalueD.r,min_pvalue=min_pvalue,acc=acc_supp,lim=lim_supp,max_core=max_core)
      pvalue.OptMin = pvalue_minimizePValue_davies(u_f=Q.f,u_r=Q.r,df_f=R,lambda_r=EigenvalueD.r,min_pvalue=min_pvalue,acc=acc,ChisqApp=chisq_app,max_core=max_core)
    }
    else pvalue.OptMin = pvalue_minimizePValue_liu(u_f=Q.f,u_r=Q.r,df_f=R,lambda_r=EigenvalueD.r,min_pvalue=min_pvalue,acc=acc,ChisqApp=chisq_app)
    pvalue["pvalue.otranScore"] = pvalue.OptMin
    stat["stat.otranScore"] = Q.OptMin
    rho["rho.otranScore.f"] = rho_min_pvalue
    rho["rho.otranScore.r"] = 1-rho_min_pvalue
  }

  if(combinations_return & (("All" %in% combination_preference) | ("AdptWt" %in% combination_preference))){
    # Adaptive weighted combination
    Z.f = -2*log(pvalue.f)
    if(is.null(log(pvalue.r))) Z.r = -2*log(1e-16)
    else Z.r = -2*log(pvalue.r)
    esti.weight.f = Z.f/sqrt(Z.f^2+Z.r^2)
    esti.weight.r = Z.r/sqrt(Z.f^2+Z.r^2)
    Q.AdptWt = sqrt(Z.f^2+Z.r^2)[1] # [1,] is added to avoid operation problems. It guarantees it's a value.
    integrand = function(y){
      pchisq(sqrt(Q.AdptWt^2-y^2),df=2,lower.tail=FALSE)*dchisq(y,df=2)
    }
    pvalue.AdptWt = pchisq(Q.AdptWt,df=2,lower.tail=FALSE)+integrate(integrand,lower=0,upper=Q.AdptWt)$value
    pvalue["pvalue.atranScore"] = pvalue.AdptWt
    stat["stat.atranScore"] = Q.AdptWt
    rho["rho.atranScore.f"] = esti.weight.f
    rho["rho.atranScore.r"] = esti.weight.r
  }

  if(combinations_return & (("All" %in% combination_preference) | ("Fisher" %in% combination_preference))){
    # Fisher's combination
    Q.Fisher = -2*log(pvalue.f)-2*log(pvalue.r)
    pvalue.Fisher = pchisq(as.numeric(Q.Fisher),df=4,lower.tail=FALSE)
    pvalue["pvalue.ftranScore"] = pvalue.Fisher
    stat["stat.ftranScore"] = Q.Fisher

  }

  ci=matrix(NA, ncol=3, nrow=length(point))
  for(i in 1:length(point)){
  ci[i,]=c(exp(point[i]), exp(c(point[i]-1.96*SE[i])),  exp(point[i]+1.96*SE[i]))
  }
  colnames(ci)=c("Odds Ratio",  "95% CI lower", "95% CI higher")

 rownames(ci)=burden_name
 rownames(pvalue.f.ind)=burden_name
 rownames(stat.f.ind)=burden_name
  #output
   return(list(stat=stat,
              pvalue=pvalue,
              stat.f.ind=stat.f.ind,
              pvalue.f.ind=pvalue.f.ind,
              odds_ratio.f.ind=ci,
              rho=rho,
              data.info=c(p=p,R=R)
  ))
  }









