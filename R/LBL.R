
#' Logistic Bayesian Lasso for Rare (or Common) Haplotype Association
#' @param dat the non-SNP and SNP data as a data frame. If the twoBinaryPheno option is FALSE (default) and the complex.sampling option is FALSE (default), the first column of the non-SNP data is the affection status, others (optional) are environmental covariates; if the complex.sampling option is set to be TRUE, the non-SNP data should consists of affection status, sampling weights, stratifying variables and environmental covariates (optional). If the twoBinaryPheno option is set to be TRUE, then there should be no environmental covariate and the first two columns should be two binary phenotypes. SNP data should comprise the last 2*numSNPs columns (allelic format) or last numSNPs columns (genotypic format). Missing allelic data should be coded as NA or "" and missing genotypic data should be coded as, e.g., "A" if one allele is missing and "" if both alleles are missing. Covariates (including stratifying variables) should be coded as dummy variables, e.g., 0, 1, etc.
#' @param numSNPs number of SNPs per haplotype.
#' @param maxMissingGenos maximum number of single-locus genotypes with missing data to allow for each subject. (Subjects with more missing data, or with missing non-SNP data are removed.) The default is 1.
#' @param allelic TRUE if single-locus SNP genotypes are in allelic format and FALSE if in genotypic format; default is TRUE.
#' @param haplo.baseline haplotype to be used for baseline coding; default is the most frequent haplotype according to the initial haplotype frequency estimates returned by pre.hapassoc.
#' @param cov.baseline Needed only if the non-SNP data contains stratifying variables or environmental covariates. Indicates the baseline level(s) for the covariates (including stratifying variables). Note that they should be listed in the same order as in the actual data. The default is the level(s) that is coded as 0 for each covariate. This option is ignored if twoBinaryPheno = TRUE.
#' @param complex.sampling whether complex sampling with frequency matching will be used; default is FALSE. Specifically, when this option is set to be TRUE, G-E and/or G-S dependence is assumed, which needs to be further specified by the names.dep option. This option is ignored if twoBinaryPheno = TRUE.
#' @param n.stra Needed only if the complex.sampling option is set to be TRUE. Indicates number of stratifying variables.
#' @param interaction.stra Needed only if the complex.sampling option is set to be TRUE. Indicates whether or not to model interaction between haplotypes and stratifying variables in the model; default is TRUE. This option is ignored if twoBinaryPheno = TRUE.
#' @param interaction.env Needed only if the non-SNP data contains environmental covariates. Indicates whether or not to model interaction between haplotypes and environmental covariates in the model; default is TRUE. This option is ignored if twoBinaryPheno = TRUE.
#' @param interaction.model Needed only if the complex.sampling option is set to be FALSE and the interaction.cov option is set to be TRUE. Indicates whether G-E independence is assumed or not for fitting haplotype-environment interactions. "i" represents G-E independent model, "d" represents G-E dependent model, and "u" represents uncertainty about G-E independence, i.e., allows possibility of both models. The default is "i". This option is ignored if twoBinaryPheno = TRUE.
#' @param names.dep Needed only if the complex.sampling option is set to be TRUE or interaction.model option is set to be "d" or "u". Indicates the covariates that are believed to cause G-E dependence. The default is a vector consisting of all covariates, however, if the number of covariates is large, then this will lead to a very large and complicated G-E dependence model so a judicious choice of covariates for this model is recommended in that case.
#' @param a first hyperparameter of the prior for regression coefficients, beta. The prior variance of beta is 2/lambda^2 and lambda has Gamma(a,b) prior. The Gamma parameters a and b are such that the mean and variance of the Gamma distribution are a/b and a/b^2. The default is 20.
#' @param b b parameter of the Gamma(a,b) distribution described above; default is 20.
#' @param start.beta starting value of all regression coefficients, beta; default is 0.01.
#' @param gamma starting value of the gamma parameters (slopes), which are used to model G-E dependence through a multinomial logistic regression model; default is 0.01. This option is ignored if twoBinaryPheno = TRUE.
#' @param lambda starting value of the lambda parameter described above; default is 1.
#' @param D starting value of the D parameter, which is the within-population inbreeding coefficient; default is 0.
#' @param e a (small) number epsilon in the null hypothesis of no association, H0: |beta| <= epsilon. Changing e from default of 0.1 may need choosing a different threshold for Bayes Factor (one of the outputs) to infer association. The default is 0.1.
#' @param seed the seed to be used for the MCMC in Bayesian Lasso; default is a random seed. If exactly same results need to be reproduced, seed should be fixed to the same number.
#' @param burn.in burn-in period of the MCMC sampling scheme; default is 20000 for model with a single univariate phenotype and 50000 for model with two binary phenotypes.
#' @param num.it total number of MCMC iterations including burn-in. When the complex.sampling option is set to be FALSE, default is 50000 if there are no covariates or interaction.model = "i"; default values are 70000 and 100000, respectively, if interaction.model = "d" and "u". When the complex.sampling option is set to be TRUE, the default value of num.it is 120000. When the twoBinaryPheno option is set to be TRUE, the default value of num.it is 200000.
#' @param twoBinaryPheno whether two binary correlated phenotypes will be used (no environmental covariate allowed). the options of "complex.sampling", "interaction.stra", "interaction.env", and "interaction.model" will be ignored; default is FALSE.
#' @param start.u Needed only if twoBinaryPheno=TRUE. Starting value of u (subject-specific latent variables); ui induces correlation between the two phenotypes of i-th individual; default is 0.01.
#' @param sigma_sq_u Needed only if twoBinaryPheno=TRUE. Starting value of sigma_sq_u parameter, which is the variance of u elements; ui is assumed to follow N(0, sigma_sq_u) distribution; default is 1.
#' @param start.f00 Needed only if twoBinaryPheno=TRUE. Starting value of the f00 parameter vector, which consists of haplotype frequencies in the population of controls for both phenotypes; if it set to be None (default), the initFreq returned by pre.hapassoc when applied to the corresponding sample will be used.
#' @param start.f10 Needed only if twoBinaryPheno=TRUE. Starting value of the f10 parameter vector, which consists of haplotype frequencies in the population of cases for the first phenotype and controls for the second phenotype; if it set to be None (default), the initFreq returned by pre.hapassoc when applied to the corresponding sample will be used.
#' @param start.f01 Needed only if twoBinaryPheno=TRUE. Starting value of the f01 parameter vector, which consists of haplotype frequencies in the population of controls for the first phenotype and cases for the second phenotype; if it set to be None (default), the initFreq returned by pre.hapassoc when applied to the corresponding sample will be used.
#' @param e_allHap Needed only if twoBinaryPheno=TRUE. Epsilon in the null hypothesis for testing all haplotypes together in a block, H0: |beta| <= epsilon for all beta coefficients corresponding to all haplotypes in a block and both diseases. The default is 0.4.
#' @param print.freq.ci Needed only if twoBinaryPheno=TRUE. Whether the 95\% credible sets for f00, f10, and f01 are to be printed. The default is FALSE.
#' @param print.lambda.ci Needed only if twoBinaryPheno=TRUE. Whether the 95\% credible set for lambda is to be printed. The default is FALSE.
#' @param print.D.ci Needed only if twoBinaryPheno=TRUE. Whether the 95\% credible set for D is to be printed. The default is FALSE.
#' @param print.sigma_sq_u.ci Needed only if twoBinaryPheno=TRUE. Whether the 95\% credible set for sigma_sq_u is to be printed. The default is FALSE.
#'
#' @return
#' \item{BF }{For single phenotype. A vector of Bayes Factors for all regression coefficients. If BF exceeds a certain threshold (e.g., 2 or 3) association may be concluded.}
#' \item{OR }{For single phenotype. A vector of estimated odds ratios of the corresponding haplotype against the reference haplotype (haplo.baseline). This is the exponential of the posterior means of the regression coefficients.}
#' \item{CI.OR }{For single phenotype. 95\% credible sets for the ORs. If CI.OR excludes 1, association may be concluded.}
#' \item{freq }{For single phenotype. A vector of posterior means of the haplotype frequencies.}
#' \item{CI.freq }{In univariate model. 95\% credible sets for each haplotype frequency.}
#' \item{percentage.indep}{For single phenotype. Available only if the interaction.model option is set to be "u". Percentage of iterations in which independent model is choosen.}
#' \item{percentage.dep}{For single phenotype. Available only if the interaction.model option is set to be "u". Percentage of iterations in which dependent model is chosen.}
#' \item{CI.gamma }{For single phenotype. Available only if the interaction.model option is set to be "d" or "u". 95\% credible sets for the gamma parameters as described above.}
#' \item{CI.lambda }{95\% credible sets for the lambda parameter as described above. For two binary phenotypes, it is optional, only shown if print.lambda.ci=TRUE}
#' \item{CI.D }{95\% credible sets for D as described above. For two binary phenotypes, it is optional, only shown if print.D.ci=TRUE}
#' \item{BF_bivariate_hap }{For two binary phenotypes. A vector of Bayes Factors for testing association of each haplotype with both phenotypes jointly. If a BF exceeds a certain threshold, the corresponding haplotype may be associated with at least one of the two phenotypes.}
#' \item{BF_bivariate_allHap }{For two binary phenotypes. The joint Bayes Factor for testing association of all haplotypes in a block together with both phenotypes jointly. If joint BF exceeds a certain threshold, then at least one of the haplotypes may be associated with at least one of the two phenotypes.}
#' \item{beta1 }{For two binary phenotypes. A vector of estimated posterior means of the regression coefficients for the first phenotype.}
#' \item{beta2 }{For two binary phenotypes. A vector of estimated posterior means of the regression coefficients for the second  phenotype.}
#' \item{CI.beta1 }{For two binary phenotypes. 95\% credible sets for the beta1s. These are based on marginal distribution of each beta1 coefficient and should not be used for inference about association with the two phenotypes jointly.}
#' \item{CI.beta2 }{For two binary phenotypes. 95\% credible sets for the beta2s. These are based on marginal distribution of each beta2 coefficient and should not be used for inference about association with the two phenotypes jointly.}
#' \item{freq00 }{For two binary phenotypes. A vector of posterior means of the haplotype frequencies in the population of controls for both phenotypes.}
#' \item{freq10 }{For two binary phenotypes. A vector of posterior means of the haplotype frequencies in the population of cases for the first phenotype and controls for the second phenotype.}
#' \item{freq01 }{For two binary phenotypes. A vector of posterior means of the haplotype frequencies in the population of controls for the first phenotypes and cases for the second phenotype.}
#' \item{CI.freq00 }{For two binary phenotypes. 95\% credible sets for each haplotype frequency in the population of controls for both phenotypes. Optional, only shown if print.freq.ci=TRUE}
#' \item{CI.freq10 }{For two binary phenotypes. 95\% credible sets for each haplotype frequency in the population of cases for the first phenotype and controls for the second phenotype. Optional, only shown if print.freq.ci=TRUE}
#' \item{CI.freq01 }{For two binary phenotypes. 95\% credible sets for each haplotype frequency in the population of controls for the first phenotypes and cases for the second phenotype. Optional, only shown if print.freq.ci=TRUE}
#' \item{CI.sigma_sq_u}{For two binary phenotypes. 95\% credible sets for sigma_sq_u as described above. Optional, only shown if print.sigma_sq_u.ci=TRUE}
#' @export LBL
#' @import hapassoc dummies
#' @importFrom stats quantile
#' @importFrom hapassoc pre.hapassoc
#' @useDynLib LBLGXE, .registration = TRUE
#' @seealso \code{\link{pre.hapassoc}}
#' @description Bayesian LASSO is used to find the posterior distributions of logistic
#' regression coefficients, which are then used to calculate Bayes Factor and credible
#' set to test for association with haplotypes, environmental covariates, and interactions.
#' It can handle complex sampling data, in particular, frequency matched cases and controls
#' with controls obtained using stratified sampling. This version can also be applied to a
#' dataset with no environmental covariate and two correlated binary phenotypes. The function
#' first calls pre.hapassoc function from the hapassoc package, and some of the options such
#' as "dat", "numSNPs", "maxMissingGenos" and "allelic" are used by pre.hapassoc. It takes as
#' an argument a dataframe with non-SNP and SNP data. The rows of the input data frame should
#' correspond to subjects. Missing single-locus genotypes, up to a maximum of maxMissingGenos
#' (see below), are allowed, but subjects with missing data in more than maxMissingGenos,
#' or with missing non-SNP data, are removed.
#' @references
#' Yuan X and Biswas S (2019). Bivariate Logistic Bayesian LASSO for Detecting Rare Haplotype Association with Two Correlated Phenotypes. Genetic Epidemiology, in press.
#'
#' Zhang Y, Hofmann J, Purdue M, Lin S, and Biswas S.
#' Logistic Bayesian LASSO for Genetic Association Analysis of Data from Complex Sampling Designs. Journal of Human Genetics, 62:819-829.
#'
#' Zhang Y, Lin S, and Biswas S.
#' Detecting Rare and common Haplotype-Environment Interaction under Uncertainty of Gene-Environment Independence Assumption. Biometrics, 73:344-355.
#'
#' Zhang, Y. and Biswas, S (2015). An Improved Version of Logistic Bayesian LASSO for Detecting Rare Haplotype-Environment Interactions With Application to Lung Cancer, Cancer Informatics, 14(S2): 11-16.
#'
#' Biswas S, Xia S and Lin S (2014). Detecting Rare Haplotype-Environment Interaction with Logistic Bayesian LASSO. Genetic Epidemiology, 38: 31-41.
#'
#' Biswas S, Lin S (2012). Logistic Bayesian LASSO for Identifying Association with Rare Haplotypes and Application to Age-related Macular Degeneration. Biometrics, 68(2): 587-97.
#'
#' Burkett K, Graham J and McNeney B (2006). hapassoc: Software for Likelihood Inference of Trait Associations with SNP Haplotypes and Other Attributes. Journal of Statistical Software, 16(2): 1-19.
#' @author Xiaochen Yuan, Yuan Zhang, Shuang Xia, Swati Biswas, Shili Lin
#' @examples
#' # Load example datasets
#' # This dataset consists of affection status, a binary environmental covariate, and SNP data.
#' data(LBL.ex1)
#' # This dataset consists of affection status, complex sampling weights, a binary stratifying
#' # variable, a binary environmental covariate, and SNP data.
#' data(LBL.ex2)
#' # This dataset consists of two correlated affection statuses, no environmental covariate,
#' #and SNP data.
#' data(LBL.ex3)
#' # Install hapassoc and dummies package
#' library(hapassoc)
#' library(dummies)
#' # Run LBL to make inference on haplotype associations and interactions. Note the default
#' # setting for burn.in and num.it are larger in the LBL function. However, you may want to
#' # use smaller numbers for a quick check to make sure the package is loaded properly. With
#' # such shorts runs, the results may not be meaningful.
#' ## Analyzing LBL.ex1 under G-E independence assumption.
#' out.LBL<-LBL(LBL.ex1, numSNPs=5, burn.in=0, num.it=5)
#'
#' ## Analyzing LBL.ex1 under uncertainty of G-E independence assumption.
#' out.LBL<-LBL(LBL.ex1, numSNPs=5, interaction.model="u", burn.in=0, num.it=5)
#'
#' ## Analyzing LBL.ex2 which comes from complex sampling design with frequency matching.
#' out.LBL<-LBL(LBL.ex2, numSNPs=5, complex.sampling=TRUE, n.stra=1, names.dep="stra",
#' burn.in=0, num.it=5)
#'
#' ## Analyzing LBL.ex3 using the bivariate LBL method.
#' out.LBL<-LBL(LBL.ex3, numSNPs=5, twoBinaryPheno=TRUE, burn.in=0, num.it=5)

LBL<-function(dat, numSNPs, maxMissingGenos = 1, allelic = TRUE, haplo.baseline = "missing", cov.baseline = "missing", complex.sampling = FALSE, n.stra = NULL, interaction.stra = TRUE, interaction.env = TRUE, interaction.model = "i", names.dep = "missing", a = 20, b = 20, start.beta = 0.01, gamma = 0.01, lambda = 1, D = 0, e = 0.1, seed = NULL, burn.in = NULL, num.it = NULL, twoBinaryPheno=FALSE, start.u=0.01,sigma_sq_u=1, start.f00=NULL, start.f10=NULL, start.f01=NULL, e_allHap=0.4, print.freq.ci=FALSE, print.lambda.ci=FALSE, print.D.ci=FALSE, print.sigma_sq_u.ci=FALSE)
{

  haplos.list<-pre.hapassoc(dat, numSNPs=numSNPs, pooling.tol=0, allelic=allelic, maxMissingGenos=maxMissingGenos)
  ####################################################################################################
  # Bivariate LBL
  ####################################################################################################
  if (twoBinaryPheno == TRUE){
    n.Y <- dim(haplos.list$nonHaploDM)[2] # num of phenotypes
    temp=cbind(haplos.list$ID, haplos.list$wt, haplos.list$nonHaploDM, haplos.list$haploDM)
    temp=temp[order(temp[,1]),]
    colnames(temp)=c("ID","wt","Y1","Y2",colnames(temp[,5:ncol(temp)]))

    haplos.list$ID=temp[,1]
    haplos.list$wt=temp[,2]
    haplos.list$nonHaploDM=temp[,3:(3+n.Y-1)]
    haplos.list$haploDM=temp[(3+n.Y):ncol(temp)]

    haplo.names <- names(haplos.list$initFreq)
    freq <- haplos.list$initFreq
    if (haplo.baseline=="missing") haplo.baseline <- haplo.names[which.max(freq)]
    column.subset <- colnames(haplos.list$haploDM) != haplo.baseline
    haplo.mat=haplos.list$haploDM[,column.subset]
    haplo.names=colnames(haplo.mat)

    hdat <- cbind(haplos.list$nonHaploDM, haplos.list$haploDM[, column.subset])

    ID <- haplos.list$ID
    N <- sum(haplos.list$wt)
    N00 <- sum(temp[temp$Y1==0&temp$Y2==0,]$wt)
    N10 <- sum(temp[temp$Y1==1&temp$Y2==0,]$wt)
    N01 <- sum(temp[temp$Y1==0&temp$Y2==1,]$wt)
    N11 <- sum(temp[temp$Y1==1&temp$Y2==1,]$wt)


    y1 <- as.numeric(hdat[, 1])
    y2 <- as.numeric(hdat[, 2])
    freq.new<-freq[column.subset]
    freq.new[length(freq.new)+1]<-freq[haplo.baseline]
    names(freq.new) = c(haplo.names, haplo.baseline)

    # calculate initial values of f00, f10, and f01
    colnames(dat) <- c("Y1","Y2",colnames(dat[,3:ncol(dat)]))
    dat00 <- dat[dat$Y1==0 & dat$Y2==0,]
    dat10 <- dat[dat$Y1==1 & dat$Y2==0,]
    dat01 <- dat[dat$Y1==0 & dat$Y2==1,]

    haplos.list<-pre.hapassoc(dat00, numSNPs=numSNPs, pooling.tol=0, allelic=allelic, maxMissingGenos=maxMissingGenos)
    f00_new <- haplos.list$initFreq[names(freq.new)]
    num_na <- sum(is.na(f00_new))
    f00_new[is.na(f00_new)] <- 0.00001
    f00_new[which.max(f00_new)] <- f00_new[which.max(f00_new)]-0.00001*num_na
    f00_new <- as.vector(f00_new)

    haplos.list<-pre.hapassoc(dat10, numSNPs=numSNPs, pooling.tol=0, allelic=allelic, maxMissingGenos=maxMissingGenos)
    f10_new <- haplos.list$initFreq[names(freq.new)]
    num_na <- sum(is.na(f10_new))
    f10_new[is.na(f10_new)] <- 0.00001
    f10_new[which.max(f10_new)] <- f10_new[which.max(f10_new)]-0.00001*num_na
    f10_new <- as.vector(f10_new)

    haplos.list<-pre.hapassoc(dat01, numSNPs=numSNPs, pooling.tol=0, allelic=allelic, maxMissingGenos=maxMissingGenos)
    f01_new <- haplos.list$initFreq[names(freq.new)]
    num_na <- sum(is.na(f01_new))
    f01_new[is.na(f01_new)] <- 0.00001
    f01_new[which.max(f01_new)] <- f01_new[which.max(f01_new)]-0.00001*num_na
    f01_new <- as.vector(f01_new)

    num.haplo.id<-as.vector(table(ID))
    x.length<-as.integer(dim(haplo.mat)[2])

    haplo.map<-matrix(rep(NA,2*(dim(haplo.mat)[1])),ncol=2)
    for(i in 1:dim(haplo.mat)[1])
    {
      for(j in 1:dim(haplo.mat)[2])
      {
        if(haplo.mat[i,j]==2)
        {
          haplo.map[i,1]=j
          haplo.map[i,2]=j
          break
        }
        if(haplo.mat[i,j]==1)
        {
          if(is.na(haplo.map[i,1])==TRUE) haplo.map[i,1]=j
          else
          {
            haplo.map[i,2]=j
            break
          }
        }
      }
      if(is.na(haplo.map[i,1])==TRUE)
      {
        haplo.map[i,1]=dim(haplo.mat)[2]+1
        haplo.map[i,2]=dim(haplo.mat)[2]+1
      }
      if(is.na(haplo.map[i,2])==TRUE) haplo.map[i,2]=dim(haplo.mat)[2]+1
    }
    uniq.mat<-unique(haplo.mat)
    uniq.map<-unique(haplo.map)

    if (is.null(start.f00)==TRUE) f00=f00_new else f00=start.f00
    if (is.null(start.f10)==TRUE) f10=f10_new else f10=start.f10
    if (is.null(start.f01)==TRUE) f01=f01_new else f01=start.f01
    if (is.null(seed)==FALSE) set.seed(seed)
    if (is.null(burn.in)==TRUE) burn.in=50000
    if (is.null(num.it)==TRUE) num.it=200000
    beta1=rep(start.beta, x.length+1)
    beta2=rep(start.beta, x.length+1)
    u=rep(start.u, N)

    beta1.out<-numeric((num.it-burn.in)*(1+x.length))
    beta2.out<-numeric((num.it-burn.in)*(1+x.length))
    lambda.out<-numeric(num.it-burn.in)
    freq00.out<-numeric((num.it-burn.in)*(x.length+1))
    freq10.out<-numeric((num.it-burn.in)*(x.length+1))
    freq01.out<-numeric((num.it-burn.in)*(x.length+1))
    D.out<-numeric(num.it-burn.in)
    #u.out=numeric((num.it-burn.in)*N)
    sigma_sq_u.out<-numeric(num.it-burn.in)

    #dyn.load("parallel_LBL_2BinaryPheno.so")
    out=.C("mcmc_twoBinaryPheno", as.integer(haplo.map),as.integer(dim(haplo.mat)[1]), as.integer(y1),as.integer(y2), as.integer(N),as.integer(N00),as.integer(N10),as.integer(N01), as.integer(num.haplo.id), as.integer(x.length), as.double(f00),as.double(f10),as.double(f01), as.double(D), as.double(beta1),as.double(beta2), as.double(a), as.double(b), as.integer(t(uniq.map)), as.integer(dim(uniq.mat)[1]), as.double(lambda),as.double(u),as.double(sigma_sq_u), as.integer(num.it), as.integer(burn.in), beta1.out=as.double(beta1.out), beta2.out=as.double(beta2.out),lambda.out=as.double(lambda.out), freq00.out=as.double(freq00.out),freq10.out=as.double(freq10.out),freq01.out=as.double(freq01.out), D.out=as.double(D.out),sigma_sq_u.out=as.double(sigma_sq_u.out))

    beta1.out<-matrix(out$beta1.out,nrow=num.it-burn.in, byrow=TRUE)
    beta2.out<-matrix(out$beta2.out,nrow=num.it-burn.in, byrow=TRUE)
    freq00.out<-matrix(out$freq00.out,nrow=num.it-burn.in, byrow=TRUE)
    freq10.out<-matrix(out$freq10.out,nrow=num.it-burn.in, byrow=TRUE)
    freq01.out<-matrix(out$freq01.out,nrow=num.it-burn.in, byrow=TRUE)
    #u.out<-matrix(out$u.out,nrow=num.it-burn.in, byrow=TRUE)
    ############################################################################
    hap.name <- colnames(hdat[, -c(1,2)])
    n.hap <- length(hap.name)
    colnames(beta1.out) <- c("beta0",hap.name)
    colnames(beta2.out) <- c("beta0",hap.name)
    beta.out <- cbind(beta1.out[,hap.name],beta2.out[,hap.name])
    ###############################################################################
    ci.beta1.95 <- numeric((1+x.length)*2)
    ci.beta2.95 <- numeric((1+x.length)*2)
    post.mean.beta1 <- numeric(1+x.length)
    post.mean.beta2 <- numeric(1+x.length)
    post.mean.freq00<-numeric(x.length+1)
    post.mean.freq10<-numeric(x.length+1)
    post.mean.freq01<-numeric(x.length+1)

    k<-1
    for (i in 1:(1+x.length)){

      ci.beta1.95[k]<-quantile(beta1.out[,i], probs=0.025)
      ci.beta1.95[k+1]<-quantile(beta1.out[,i], probs=0.975)
      ci.beta2.95[k]<-quantile(beta2.out[,i], probs=0.025)
      ci.beta2.95[k+1]<-quantile(beta2.out[,i], probs=0.975)
      post.mean.beta1[i]<- mean(beta1.out[,i])
      post.mean.beta2[i]<- mean(beta2.out[,i])
      post.mean.freq00[i]<- mean(freq00.out[,i])
      post.mean.freq10[i]<- mean(freq10.out[,i])
      post.mean.freq01[i]<- mean(freq01.out[,i])
      k<-k+2
    }

    ci.beta1.95<-matrix(ci.beta1.95,nrow=1+x.length, ncol=2, byrow=TRUE)
    ci.beta2.95<-matrix(ci.beta2.95,nrow=1+x.length, ncol=2, byrow=TRUE)

    BF_bivariate_hap=rep(0,x.length)
    max_all_beta=apply(abs(beta.out),1,max)

    # BF_bivariate_allHap
    prior_H0=Hprob(x.length*2, 0, a, b, e_allHap)
    prior_Ha=1-prior_H0
    post_H0=length(max_all_beta[max_all_beta<=e_allHap])/length(max_all_beta)
    post_Ha=1-post_H0
    all_BF=(post_Ha/post_H0)/(prior_Ha/prior_H0)

    # BF_bivariate_hap
    for(j in 1:x.length){
      ebeta=beta.out[,c(j,j+x.length)]
      maxbeta=apply(abs(ebeta),1,max)
      prior_H0=Hprob(2, 0, a, b, e)
      prior_Ha=1-prior_H0
      post_H0=length(maxbeta[maxbeta<=e])/length(maxbeta)
      post_Ha=1-post_H0
      BF_bivariate_hap[j]=round((post_Ha/post_H0)/(prior_Ha/prior_H0), 4)
      if(BF_bivariate_hap[j]>100) BF_bivariate_hap[j] <- ">100"
    }
    names(BF_bivariate_hap) <- hap.name
    post.mean.beta1<-round(post.mean.beta1,3)
    post.mean.beta2<-round(post.mean.beta2,3)
    ci.beta1.95<-round(ci.beta1.95,2)
    ci.beta2.95<-round(ci.beta2.95,2)

    names(post.mean.beta1)<-c("beta0",hap.name)
    names(post.mean.beta2)<-c("beta0",hap.name)
    ci.beta1.95<-data.frame(c("beta0",hap.name),ci.beta1.95)
    colnames(ci.beta1.95)<-c("Hap", "Lower", "Upper")
    ci.beta2.95<-data.frame(c("beta0",hap.name),ci.beta2.95)
    colnames(ci.beta2.95)<-c("Hap", "Lower", "Upper")
    names(post.mean.freq00)<-c(hap.name, haplo.baseline)
    names(post.mean.freq10)<-c(hap.name, haplo.baseline)
    names(post.mean.freq01)<-c(hap.name, haplo.baseline)

    # CI for lambda (optional)
    if (print.lambda.ci==TRUE){
      ci.lambda<-numeric(2)
      ci.lambda<-c(quantile(out$lambda.out, probs=0.025),quantile(out$lambda.out, probs=0.975))
    } else ci.lambda <- ""

    # CI for D (optional)
    if (print.D.ci==TRUE){
      ci.D<-numeric(2)
      ci.D<-c(quantile(out$D.out, probs=0.025),quantile(out$D.out, probs=0.975))
    } else ci.D <- ""

    # CI for sigma_sq_u (optional)
    if (print.sigma_sq_u.ci==TRUE){
      ci.sigma_sq_u <- numeric(2)
      ci.sigma_sq_u<-c(quantile(out$sigma_sq_u.out, probs=0.025),quantile(out$sigma_sq_u.out, probs=0.975))
    } else ci.sigma_sq_u <- ""

    # CI for frequencies (optional)
    if (print.freq.ci==TRUE){
      ci.freq00<-numeric((x.length+1)*2)
      ci.freq10<-numeric((x.length+1)*2)
      ci.freq01<-numeric((x.length+1)*2)
      k<-1
      for (i in 1:(1+x.length)){

        ci.freq00[k]<-quantile(freq00.out[,i], probs=0.025)
        ci.freq00[k+1]<-quantile(freq00.out[,i], probs=0.975)
        ci.freq10[k]<-quantile(freq10.out[,i], probs=0.025)
        ci.freq10[k+1]<-quantile(freq10.out[,i], probs=0.975)
        ci.freq01[k]<-quantile(freq01.out[,i], probs=0.025)
        ci.freq01[k+1]<-quantile(freq01.out[,i], probs=0.975)
        k<-k+2
      }
      ci.freq00<-matrix(ci.freq00,nrow=(x.length+1), ncol=2, byrow=TRUE)
      ci.freq10<-matrix(ci.freq10,nrow=(x.length+1), ncol=2, byrow=TRUE)
      ci.freq01<-matrix(ci.freq01,nrow=(x.length+1), ncol=2, byrow=TRUE)
      ci.freq00<-data.frame(c(hap.name,haplo.baseline), ci.freq00)
      ci.freq10<-data.frame(c(hap.name,haplo.baseline), ci.freq10)
      ci.freq01<-data.frame(c(hap.name,haplo.baseline), ci.freq01)
      colnames(ci.freq00)<-c("Hap", "Lower", "Upper")
      colnames(ci.freq10)<-c("Hap", "Lower", "Upper")
      colnames(ci.freq01)<-c("Hap", "Lower", "Upper")

    }
    else{
      ci.freq00<-""
      ci.freq10<-""
      ci.freq01<-""
    }

    #################################################################################
    ans <- list(BF_bivariate_hap=BF_bivariate_hap, BF_bivariate_allHap=all_BF,
                beta1=post.mean.beta1,beta2=post.mean.beta2, CI.beta1=ci.beta1.95, CI.beta2=ci.beta2.95,
                freq00=post.mean.freq00, freq10=post.mean.freq10, freq01=post.mean.freq01, CI.freq00=ci.freq00, CI.freq10=ci.freq10, CI.freq01=ci.freq01,
                CI.lambda=ci.lambda,CI.D=ci.D, CI.sigma_sq_u=ci.sigma_sq_u)
  }

  ####################################################################################################
  #Non-complex sampling methods
  ####################################################################################################
  if((twoBinaryPheno == FALSE) & (complex.sampling==FALSE))
  {
    var.baseline = cov.baseline
    interaction.cov = interaction.env
    n.cov=dim(haplos.list$nonHaploDM)[2]-1#num of covariates
    cov.baseline=var.baseline

    ##################################################
    #LBL
    ##################################################
    if(n.cov==0){

      temp=cbind(haplos.list$ID, haplos.list$wt, haplos.list$nonHaploDM, haplos.list$haploDM)
      temp=temp[order(temp[,1]),]
      haplos.list$ID=temp[,1]
      haplos.list$wt=temp[,2]
      haplos.list$nonHaploDM=temp[,3]
      haplos.list$haploDM=temp[4:ncol(temp)]

      if (is.null(seed)==FALSE) set.seed(seed)
      if (is.null(burn.in)==TRUE) burn.in=20000
      if (is.null(num.it)==TRUE) num.it=50000
      haplo.names <- names(haplos.list$initFreq)
      freq <- haplos.list$initFreq
      freq = freq[names(haplos.list$haploDM)]
      if (haplo.baseline=="missing") haplo.baseline <- haplo.names[which.max(freq)]
      column.subset <- colnames(haplos.list$haploDM) != haplo.baseline
      hdat <- cbind(haplos.list$nonHaploDM, haplos.list$haploDM[, column.subset])
      colnames(hdat) <- c("affected", colnames(haplos.list$haploDM[, column.subset]))
      ID <- haplos.list$ID
      N <- sum(haplos.list$wt)
      y <- as.numeric(hdat[, 1])# y <- as.numeric(hdat[colnames(hdat)!="pheno"])
      x <- data.matrix(hdat[, -1])# x <- data.matrix(hdat[colnames(hdat)!="pheno"])
      colnames(x)<-NULL

      freq.new<-freq[names(haplos.list$haploDM[, column.subset])]
      freq.new[length(freq.new)+1]<-freq[haplo.baseline]
      freq.new<-as.vector(freq.new)
      num.haplo.id<-as.vector(table(ID))
      x.length<-as.integer(dim(x)[2])
      unique.x<-unique(x)

      beta=rep(start.beta, x.length)
      beta.out<-numeric((num.it-burn.in)*x.length)
      lambda.out<-numeric(num.it-burn.in)
      freq.out<-numeric((num.it-burn.in)*(x.length+1))
      D.out<-numeric(num.it-burn.in)

      out<-.C("mcmc_without_cov", as.integer(x), as.integer(dim(x)[1]), as.integer(y), as.integer(N), as.integer(num.haplo.id), as.integer(x.length), as.double(freq.new), as.double(D), as.double(beta), as.double(a), as.double(b), as.integer(t(unique.x)), as.integer(dim(unique.x)[1]), as.double(lambda), as.integer(num.it), as.integer(burn.in), beta.out=as.double(beta.out), lambda.out=as.double(lambda.out), freq.out=as.double(freq.out), D.out=as.double(D.out))

      beta.out<-matrix(out$beta.out,nrow=num.it-burn.in, byrow=TRUE)
      freq.out<-matrix(out$freq.out,nrow=num.it-burn.in, byrow=TRUE)

      post.mean.beta<-numeric(x.length)
      ci.beta<-numeric(x.length*2)
      post.mean.freq<-numeric(x.length+1)
      ci.freq<-numeric((x.length+1)*2)
      ci.lambda<-numeric(2)
      ci.D<-numeric(2)

      k<-1
      for (i in 1:x.length)
      {
        ci.beta[k]<-quantile(beta.out[,i], probs=0.025)
        ci.beta[k+1]<-quantile(beta.out[,i], probs=0.975)
        k<-k+2
        post.mean.beta[i]<- mean(beta.out[,i])
      }
      ci.beta<-matrix(ci.beta,nrow=x.length, ncol=2, byrow=TRUE)
      ci.OR<-data.frame(exp(ci.beta))
      OR<-exp(post.mean.beta)

      k<-1
      for (i in 1:(x.length+1))
      {
        ci.freq[k]<-quantile(freq.out[,i], probs=0.025)
        ci.freq[k+1]<-quantile(freq.out[,i], probs=0.975)
        k<-k+2
        post.mean.freq[i]<- mean(freq.out[,i])
      }
      ci.freq<-matrix(ci.freq,nrow=(x.length+1), ncol=2, byrow=TRUE)

      ci.lambda<-c(quantile(out$lambda.out, probs=0.025),quantile(out$lambda.out, probs=0.975))
      ci.D<-c(quantile(out$D.out, probs=0.025),quantile(out$D.out, probs=0.975))

      prob.alt<-numeric(x.length)
      BF<-numeric(x.length)

      for (i in 1:x.length)
      {
        prob.alt[i]<-length(beta.out[,i][abs(beta.out[,i]) > e])/(num.it-burn.in)
      }

      for (i in 1:x.length)
      {
        prior.prob<-(b/(e+b))^a
        prior.odds<-prior.prob/(1-prior.prob)
        if (prob.alt[i]<=(100*prior.odds)/(100*prior.odds+1))
        {
          BF[i]<-round((prob.alt[i]/(1-prob.alt[i]))/prior.odds,4)
        } else
          BF[i]<-">100"
      }

      OR<-round(OR,4)
      ci.OR<-round(ci.OR,4)
      post.mean.freq<-round(post.mean.freq,4)
      ci.freq<-round(ci.freq,4)
      ci.lambda<-round(ci.lambda,4)
      ci.D<-round(ci.D,4)

      names(BF)<-colnames(hdat[, -1])
      names(OR)<-colnames(hdat[, -1])

      ci.OR<-data.frame(colnames(hdat[, -1]),ci.OR)
      colnames(ci.OR)<-c("Hap", "Lower", "Upper")

      ci.freq<-data.frame(c(colnames(hdat[, -1]),haplo.baseline), ci.freq)
      colnames(ci.freq)<-c("Hap", "Lower", "Upper")
      names(post.mean.freq)<-c(colnames(hdat[, -1]), haplo.baseline)

      names(ci.lambda)<-c("Lower", "Upper")
      names(ci.D)<-c("Lower", "Upper")

      ans <- list(BF=BF, OR=OR, CI.OR=ci.OR, freq=post.mean.freq, CI.freq=ci.freq, CI.lambda=ci.lambda, CI.D=ci.D)
    }#end here

    ##################################################
    #LBL-GXE-I
    ##################################################
    if(n.cov>=1&&interaction.model=="i"){

      temp=cbind(haplos.list$ID, haplos.list$wt, haplos.list$nonHaploDM, haplos.list$haploDM)
      temp=temp[order(temp[,1]),]
      haplos.list$ID=temp[,1]
      haplos.list$wt=temp[,2]
      haplos.list$nonHaploDM=temp[,3:(3+n.cov)]
      haplos.list$haploDM=temp[(3+n.cov+1):ncol(temp)]

      if (is.null(seed)==FALSE) set.seed(seed)
      if (is.null(burn.in)==TRUE) burn.in=20000
      if (is.null(num.it)==TRUE) num.it=50000
      freq <- haplos.list$initFreq
      freq = freq[names(haplos.list$haploDM)]
      freq.names <- names(freq)
      if (haplo.baseline=="missing") haplo.baseline <- freq.names[which.max(freq)]
      column.subset <- colnames(haplos.list$haploDM) != haplo.baseline
      haplo.mat=haplos.list$haploDM[,column.subset]
      haplo.names=colnames(haplo.mat)

      freq.new<-freq[column.subset]
      freq.new[length(freq.new)+1]<-freq[haplo.baseline]

      if (cov.baseline=="missing") cov.baseline = rep(0, n.cov)
      haplos.list$nonHaploDM=as.matrix(haplos.list$nonHaploDM)
      for(i in 1:n.cov)
      {
        cov.temp=dummy(haplos.list$nonHaploDM[,i+1])#transform the ith covariate into dummy variables
        cov.levels=sort(unique(haplos.list$nonHaploDM[,i+1]))
        cov.subset=(cov.levels!=cov.baseline[i])
        cov.temp=cov.temp[,cov.subset]#remove the baseline from dummy variables
        colnames.cov.temp=paste(colnames(haplos.list$nonHaploDM)[i+1], cov.levels[cov.subset], sep="")
        if(i==1)
        {
          cov.mat=cov.temp
          cov.names=colnames.cov.temp
        }else{
          cov.mat=cbind(cov.mat, cov.temp)#bind all the covariates into one matrix
          cov.names=c(cov.names, colnames.cov.temp)#bind all the covariate names into one matrix
        }
      }

      y <- haplos.list$nonHaploDM[,1]
      N <- sum(haplos.list$wt)
      n=dim(haplo.mat)[1]
      num.haplo.id=as.vector(table(haplos.list$ID))

      h.length=ncol(haplo.mat)+1
      if(interaction.cov==TRUE){
        if (is.matrix(cov.mat)==TRUE){x.length=ncol(haplo.mat)+ncol(cov.mat)+ncol(haplo.mat)*ncol(cov.mat)
        }else{x.length=ncol(haplo.mat)+1+ncol(haplo.mat)}
      }else{
        if (is.matrix(cov.mat)==TRUE){x.length=ncol(haplo.mat)+ncol(cov.mat)
        }else{x.length=ncol(haplo.mat)+1}
      }

      haplo.map<-matrix(rep(NA,2*(dim(haplo.mat)[1])),ncol=2)
      for(i in 1:dim(haplo.mat)[1])
      {
        for(j in 1:dim(haplo.mat)[2])
        {
          if(haplo.mat[i,j]==2)
          {
            haplo.map[i,1]=j
            haplo.map[i,2]=j
            break
          }
          if(haplo.mat[i,j]==1)
          {
            if(is.na(haplo.map[i,1])==TRUE) haplo.map[i,1]=j
            else
            {
              haplo.map[i,2]=j
              break
            }
          }
        }
        if(is.na(haplo.map[i,1])==TRUE)
        {
          haplo.map[i,1]=dim(haplo.mat)[2]+1
          haplo.map[i,2]=dim(haplo.mat)[2]+1
        }
        if(is.na(haplo.map[i,2])==TRUE) haplo.map[i,2]=dim(haplo.mat)[2]+1
      }

      uniq.mat<-unique(haplo.mat)
      uniq.map<-matrix(rep(NA,2*(dim(uniq.mat)[1])),ncol=2)
      for(i in 1:dim(uniq.mat)[1])
      {
        for(j in 1:dim(uniq.mat)[2])
        {
          if(uniq.mat[i,j]==2)
          {
            uniq.map[i,1]=j
            uniq.map[i,2]=j
            break
          }
          if(uniq.mat[i,j]==1)
          {
            if(is.na(uniq.map[i,1])==TRUE) uniq.map[i,1]=j
            else
            {
              uniq.map[i,2]=j
              break
            }
          }
        }
        if(is.na(uniq.map[i,1])==TRUE)
        {
          uniq.map[i,1]=dim(uniq.mat)[2]+1
          uniq.map[i,2]=dim(uniq.mat)[2]+1
        }
        if(is.na(uniq.map[i,2])==TRUE) uniq.map[i,2]=dim(uniq.mat)[2]+1
      }

      if (is.matrix(cov.mat)==TRUE){cov.vec=apply(cov.mat, 1, function(x) paste(x, collapse=""))#concatenate each row of cov.mat
      }else{cov.vec=cov.mat}
      cov.vec.new=unique(cbind(haplos.list$ID, cov.vec))[,-1]
      y.new=unique(cbind(haplos.list$ID, y))[,-1]
      num.E=as.vector(table(cov.vec.new))#distribution of E's

      value.E=unique(cov.mat)#the dimention of value.E is len.E X len.dummy
      if (is.matrix(value.E)==TRUE){
        value.E=value.E[order(unique(cov.vec)),]#reorder value.E to be consistant with num.E
        len.E=nrow(value.E)#num of different E's
        len.dummy=ncol(value.E)#num of dummy variables
      }else{
        value.E=value.E[order(unique(cov.vec))]#reorder value.E to be consistant with num.E
        len.E=length(value.E)#num of different E's
        len.dummy=1#num of dummy variables
      }

      if (is.matrix(value.E)==TRUE){value.E.vec=apply(value.E, 1, function(x) paste(x, collapse=""))#concatenate each row of value.E
      }else{value.E.vec=value.E}
      index.E=match(cov.vec.new, value.E.vec)

      beta=rep(start.beta, x.length+1)
      beta.out<-numeric((num.it-burn.in)*(x.length+1))
      lambda.out<-numeric(num.it-burn.in)
      freq.out<-numeric((num.it-burn.in)*(h.length))
      D.out<-numeric(num.it-burn.in)

      if(interaction.cov==FALSE)no.interaction=1
      else no.interaction=0

      out<-.C("mcmc_indep", as.integer(y.new), as.integer(N), as.integer(n), as.integer(num.haplo.id), as.integer(x.length), as.integer(h.length), as.integer(haplo.map), as.integer(dim(uniq.mat)[1]), as.integer(t(uniq.map)), as.integer(index.E), as.integer(t(value.E)), as.integer(num.E), as.integer(len.E), as.integer(len.dummy), as.double(beta), as.double(lambda), as.double(freq.new), as.double(D), as.double(a), as.double(b), beta.out=as.double(beta.out), lambda.out=as.double(lambda.out), freq.out=as.double(freq.out), D.out=as.double(D.out), as.integer(num.it), as.integer(burn.in), as.integer(no.interaction))

      beta.out<-matrix(out$beta.out,nrow=num.it-burn.in, byrow=TRUE)
      freq.out<-matrix(out$freq.out,nrow=num.it-burn.in, byrow=TRUE)

      post.mean.beta<-numeric(x.length+1)
      ci.beta<-numeric((x.length+1)*2)
      post.mean.freq<-numeric(h.length)
      ci.freq<-numeric((h.length)*2)
      ci.lambda<-numeric(2)
      ci.D<-numeric(2)

      k<-1
      for (i in 1:(x.length+1))
      {
        ci.beta[k]<-quantile(beta.out[,i], probs=0.025)
        ci.beta[k+1]<-quantile(beta.out[,i], probs=0.975)
        k<-k+2
        post.mean.beta[i]<-mean(beta.out[,i])
      }
      ci.beta<-matrix(ci.beta,nrow=x.length+1, ncol=2, byrow=TRUE)
      ci.OR<-data.frame(exp(ci.beta))
      OR<-exp(post.mean.beta)

      k<-1
      for (i in 1:(h.length))
      {
        ci.freq[k]<-quantile(freq.out[,i], probs=0.025)
        ci.freq[k+1]<-quantile(freq.out[,i], probs=0.975)
        k<-k+2
        post.mean.freq[i]<- mean(freq.out[,i])
      }

      ci.lambda<-c(quantile(out$lambda.out, probs=0.025),quantile(out$lambda.out, probs=0.975))
      ci.D<-c(quantile(out$D.out, probs=0.025),quantile(out$D.out, probs=0.975))

      prob.alt<-numeric(x.length+1)
      BF<-numeric(x.length+1)

      for (i in 1:(x.length+1))
      {
        prob.alt[i]<-length(beta.out[,i][abs(beta.out[,i]) > e])/(num.it-burn.in)
      }

      for (i in 1:(x.length+1))
      {
        prior.prob<-(b/(e+b))^a
        prior.odds<-prior.prob/(1-prior.prob)
        if (prob.alt[i]<=(100*prior.odds)/(100*prior.odds+1))
        {
          BF[i]<-round((prob.alt[i]/(1-prob.alt[i]))/prior.odds,4)
        } else
          BF[i]<-">100"
      }

      OR<-round(OR,4)
      ci.OR<-round(ci.OR,4)
      post.mean.freq<-round(post.mean.freq,4)
      ci.freq<-round(ci.freq,4)
      ci.lambda<-round(ci.lambda,4)
      ci.D<-round(ci.D,4)

      if(interaction.cov==TRUE){
        for(i in 1:len.dummy)
        {
          name.paste.temp<-paste(cov.names[i],haplo.names,sep=".")
          if(i==1){name.paste=name.paste.temp
          }else{name.paste=c(name.paste,name.paste.temp)}
        }
        pname<-c("beta0",haplo.names,name.paste,cov.names)
      }else{
        pname<-c("beta0",haplo.names,cov.names)
      }
      names(BF)<-pname
      names(OR)<-pname

      ci.OR<-data.frame(pname,ci.OR)
      colnames(ci.OR)<-c("Name","Lower","Upper")

      ci.freq<-matrix(ci.freq,nrow=(h.length),ncol=2,byrow=TRUE)
      ci.freq<-data.frame(c(haplo.names,haplo.baseline),ci.freq)
      colnames(ci.freq)<-c("Name","Lower","Upper")
      names(post.mean.freq)<-c(haplo.names,haplo.baseline)

      names(ci.lambda)<-c("Lower","Upper")
      names(ci.D)<-c("Lower","Upper")

      ans <- list(BF=BF, OR=OR, CI.OR=ci.OR, freq=post.mean.freq, CI.freq=ci.freq, CI.lambda=ci.lambda, CI.D=ci.D)
    }#ends here

    ##################################################
    #LBL-GXE-D
    ##################################################
    if(n.cov>=1&&interaction.model=="d"){

      names.cov=colnames(haplos.list$nonHaploDM)[2:(1+n.cov)]
      if (names.dep=="missing") names.dep=names.cov
      temp=cbind(haplos.list$ID, haplos.list$wt, haplos.list$nonHaploDM, haplos.list$haploDM)
      temp=temp[order(temp[,1]),]
      haplos.list$ID=temp[,1]
      haplos.list$wt=temp[,2]
      haplos.list$nonHaploDM=temp[,3:(3+n.cov)]
      haplos.list$haploDM=temp[(3+n.cov+1):ncol(temp)]

      if (is.null(seed)==FALSE) set.seed(seed)
      if (is.null(burn.in)==TRUE) burn.in=20000
      if (is.null(num.it)==TRUE) num.it=70000
      freq <- haplos.list$initFreq
      freq = freq[names(haplos.list$haploDM)]
      freq.names <- names(freq)
      if (haplo.baseline=="missing") haplo.baseline <- freq.names[which.max(freq)]
      column.subset <- colnames(haplos.list$haploDM) != haplo.baseline
      haplo.mat=haplos.list$haploDM[,column.subset]
      haplo.names=colnames(haplo.mat)

      freq.new<-freq[column.subset]
      freq.new[length(freq.new)+1]<-freq[haplo.baseline]

      if (cov.baseline=="missing") cov.baseline = rep(0, n.cov)
      haplos.list$nonHaploDM=as.matrix(haplos.list$nonHaploDM)
      num.dummies.each.cov=rep(0, n.cov)
      for(i in 1:n.cov)
      {
        cov.temp=dummy(haplos.list$nonHaploDM[,i+1])
        cov.levels=sort(unique(haplos.list$nonHaploDM[,i+1]))
        cov.subset=(cov.levels!=cov.baseline[i])
        cov.temp=cov.temp[,cov.subset]
        colnames.cov.temp=paste(colnames(haplos.list$nonHaploDM)[i+1], cov.levels[cov.subset], sep="")
        num.dummies.each.cov[i]=length(cov.levels)-1
        if(i==1)
        {
          cov.mat=cov.temp
          cov.names=colnames.cov.temp
        }else{
          cov.mat=cbind(cov.mat, cov.temp)
          cov.names=c(cov.names, colnames.cov.temp)
        }
      }

      y <- haplos.list$nonHaploDM[,1]
      N <- sum(haplos.list$wt)
      n=dim(haplo.mat)[1]
      num.haplo.id=as.vector(table(haplos.list$ID))

      h.length=ncol(haplo.mat)+1
      if (is.matrix(cov.mat)==TRUE){x.length=ncol(haplo.mat)+ncol(cov.mat)+ncol(haplo.mat)*ncol(cov.mat)
      }else{x.length=ncol(haplo.mat)+1+ncol(haplo.mat)}

      haplo.map<-matrix(rep(NA,2*(dim(haplo.mat)[1])),ncol=2)
      for(i in 1:dim(haplo.mat)[1])
      {
        for(j in 1:dim(haplo.mat)[2])
        {
          if(haplo.mat[i,j]==2)
          {
            haplo.map[i,1]=j
            haplo.map[i,2]=j
            break
          }
          if(haplo.mat[i,j]==1)
          {
            if(is.na(haplo.map[i,1])==TRUE) haplo.map[i,1]=j
            else
            {
              haplo.map[i,2]=j
              break
            }
          }
        }
        if(is.na(haplo.map[i,1])==TRUE)
        {
          haplo.map[i,1]=dim(haplo.mat)[2]+1
          haplo.map[i,2]=dim(haplo.mat)[2]+1
        }
        if(is.na(haplo.map[i,2])==TRUE) haplo.map[i,2]=dim(haplo.mat)[2]+1
      }

      uniq.mat<-unique(haplo.mat)
      uniq.map<-matrix(rep(NA,2*(dim(uniq.mat)[1])),ncol=2)
      for(i in 1:dim(uniq.mat)[1])
      {
        for(j in 1:dim(uniq.mat)[2])
        {
          if(uniq.mat[i,j]==2)
          {
            uniq.map[i,1]=j
            uniq.map[i,2]=j
            break
          }
          if(uniq.mat[i,j]==1)
          {
            if(is.na(uniq.map[i,1])==TRUE) uniq.map[i,1]=j
            else
            {
              uniq.map[i,2]=j
              break
            }
          }
        }
        if(is.na(uniq.map[i,1])==TRUE)
        {
          uniq.map[i,1]=dim(uniq.mat)[2]+1
          uniq.map[i,2]=dim(uniq.mat)[2]+1
        }
        if(is.na(uniq.map[i,2])==TRUE) uniq.map[i,2]=dim(uniq.mat)[2]+1
      }

      if (is.matrix(cov.mat)==TRUE){cov.vec=apply(cov.mat, 1, function(x) paste(x, collapse=""))
      }else{cov.vec=cov.mat}
      cov.vec.new=unique(cbind(haplos.list$ID, cov.vec))[,-1]
      y.new=unique(cbind(haplos.list$ID, y))[,-1]
      num.E=as.vector(table(cov.vec.new))

      value.E=unique(cov.mat)
      if (is.matrix(value.E)==TRUE){
        value.E=value.E[order(unique(cov.vec)),]
        len.E=nrow(value.E)
        len.dummy=ncol(value.E)
      }else{
        value.E=value.E[order(unique(cov.vec))]
        len.E=length(value.E)
        len.dummy=1
      }

      if (is.matrix(value.E)==TRUE){value.E.vec=apply(value.E, 1, function(x) paste(x, collapse=""))
      }else{value.E.vec=value.E}
      index.E=match(cov.vec.new, value.E.vec)

      n.dep=length(names.dep)
      which.dep=rep(0, n.dep)
      which.dep=match(names.dep, names.cov)
      which.dep=sort(which.dep)

      which.dummy.dep=rep(0, sum(num.dummies.each.cov[which.dep]))
      len.dummy.dep=length(which.dummy.dep)
      l=0
      for (i in 1:n.dep)
      {
        step=num.dummies.each.cov[which.dep[i]]
        end.dummy.dep=sum(num.dummies.each.cov[1:which.dep[i]])
        which.dummy.dep[(l+1):(l+step)]=(end.dummy.dep+1-step):end.dummy.dep
        l=l+step
      }

      if (is.matrix(value.E)==TRUE){value.dep=value.E[,which.dummy.dep]
      }else{value.dep=value.E}
      if (is.matrix(value.dep)==TRUE){value.dep.vec=apply(value.dep, 1, function(x) paste(x, collapse=""))
      }else{value.dep.vec=value.dep}
      o=order(value.dep.vec)
      t=table(value.dep.vec)
      if (is.matrix(value.dep)==TRUE){value.dep=unique(value.dep[o,])
      }else{value.dep=unique(value.dep[o])}
      len.dep=length(t)
      index.Etodep=rep(0, len.E)
      l=0
      for (i in 1:len.dep)
      {
        index.Etodep[o[(l+1):(l+t[i])]]=i
        l=l+t[i]
      }

      if (is.matrix(cov.mat)==TRUE){dep.mat=cov.mat[,which.dummy.dep]
      }else{dep.mat=cov.mat}
      if (is.matrix(dep.mat)==TRUE){dep.vec=apply(dep.mat, 1, function(x) paste(x, collapse=""))
      }else{dep.vec=dep.mat}

      j=1
      for(i in which.dep)
      {
        dep.temp0=dummy(haplos.list$nonHaploDM[,i+1])
        dep.levels0=sort(unique(haplos.list$nonHaploDM[,i+1]))
        colnames.dep.temp0=paste(colnames(haplos.list$nonHaploDM)[i+1], dep.levels0, sep="")
        if(j==1)
        {
          dep.mat0=dep.temp0
          dep.names0=colnames.dep.temp0
        }else{
          dep.mat0=cbind(dep.mat0, dep.temp0)
          dep.names0=c(dep.names0, colnames.dep.temp0)
        }
        j=j+1
      }

      value.dep0=unique(dep.mat0)
      value.dep0=value.dep0[order(unique(dep.vec)),]

      beta=rep(start.beta, x.length+1)
      beta.out<-numeric((num.it-burn.in)*(x.length+1))
      lambda.out<-numeric(num.it-burn.in)
      b.out<-numeric((num.it-burn.in)*(h.length-1)*(len.dummy.dep+1))
      freq.out<-numeric((num.it-burn.in)*len.dep*h.length)
      D.out<-numeric(num.it-burn.in)

      out<-.C("mcmc_dep", as.integer(y.new), as.integer(N), as.integer(n), as.integer(num.haplo.id), as.integer(x.length), as.integer(h.length), as.integer(haplo.map), as.integer(dim(uniq.mat)[1]), as.integer(t(uniq.map)), as.integer(index.E), as.integer(t(value.E)), as.integer(num.E), as.integer(len.E), as.integer(len.dummy), as.integer(index.Etodep), as.integer(t(value.dep)), as.integer(len.dep), as.integer(len.dummy.dep), as.double(beta), as.double(lambda), as.double(freq.new), as.double(D), as.double(a), as.double(b), beta.out=as.double(beta.out), lambda.out=as.double(lambda.out), b.out=as.double(b.out), freq.out=as.double(freq.out), D.out=as.double(D.out), as.integer(num.it), as.integer(burn.in), as.double(gamma))

      beta.out<-matrix(out$beta.out,nrow=num.it-burn.in, byrow=TRUE)
      b.out<-matrix(out$b.out,nrow=num.it-burn.in, byrow=TRUE)
      freq.out<-matrix(out$freq.out,nrow=num.it-burn.in, byrow=TRUE)

      post.mean.beta<-numeric(x.length+1)
      ci.beta<-numeric((x.length+1)*2)
      post.mean.freq<-numeric(len.dep*h.length)
      ci.freq<-numeric(len.dep*h.length*2)
      ci.b<-numeric((h.length-1)*(len.dummy.dep+1)*2)
      ci.lambda<-numeric(2)
      ci.D<-numeric(2)

      k<-1
      for (i in 1:(x.length+1))
      {
        ci.beta[k]<-quantile(beta.out[,i], probs=0.025)
        ci.beta[k+1]<-quantile(beta.out[,i], probs=0.975)
        k<-k+2
        post.mean.beta[i]<-mean(beta.out[,i])
      }
      ci.beta<-matrix(ci.beta,nrow=x.length+1, ncol=2, byrow=TRUE)
      ci.OR<-data.frame(exp(ci.beta))
      OR<-exp(post.mean.beta)

      k<-1
      for (i in 1:(len.dep*h.length))
      {
        ci.freq[k]<-quantile(freq.out[,i], probs=0.025)
        ci.freq[k+1]<-quantile(freq.out[,i], probs=0.975)
        k<-k+2
        post.mean.freq[i]<-mean(freq.out[,i])
      }
      post.mean.freq<-matrix(post.mean.freq, nrow=h.length, byrow=FALSE)
      ci.freq=round(ci.freq,4)
      ci.freq<-matrix(ci.freq, ncol=2, byrow=TRUE)
      ci.freq<-apply(ci.freq, 1, function(x) paste(x, collapse=", "))
      ci.freq<-paste("(", ci.freq, sep="")
      ci.freq<-paste(ci.freq, ")", sep="")
      ci.freq<-matrix(ci.freq, nrow=h.length, byrow=FALSE)

      k<-1
      for (i in 1:((h.length-1)*(len.dummy.dep+1)))
      {
        ci.b[k]<-quantile(b.out[,i], probs=0.025)
        ci.b[k+1]<-quantile(b.out[,i], probs=0.975)
        k<-k+2
      }
      ci.b=round(ci.b,4)
      ci.b<-matrix(ci.b, ncol=2, byrow=TRUE)
      ci.b<-apply(ci.b, 1, function(x) paste(x, collapse=", "))
      ci.b<-paste("(", ci.b, sep="")
      ci.b<-paste(ci.b, ")", sep="")
      ci.b<-matrix(ci.b, nrow=(h.length-1), byrow=TRUE)

      ci.lambda<-c(quantile(out$lambda.out, probs=0.025),quantile(out$lambda.out, probs=0.975))
      ci.D<-c(quantile(out$D.out, probs=0.025),quantile(out$D.out, probs=0.975))

      prob.alt<-numeric(x.length+1)
      BF<-numeric(x.length+1)

      for (i in 1:(x.length+1))
      {
        prob.alt[i]<-length(beta.out[,i][abs(beta.out[,i]) > e])/(num.it-burn.in)
      }

      for (i in 1:(x.length+1))
      {
        prior.prob<-(b/(e+b))^a
        prior.odds<-prior.prob/(1-prior.prob)
        if (prob.alt[i]<=(100*prior.odds)/(100*prior.odds+1))
        {
          BF[i]<-round((prob.alt[i]/(1-prob.alt[i]))/prior.odds,4)
        } else
          BF[i]<-">100"
      }

      OR<-round(OR,4)
      ci.OR<-round(ci.OR,4)
      post.mean.freq<-round(post.mean.freq,4)
      ci.lambda<-round(ci.lambda,4)
      ci.D<-round(ci.D,4)

      for(i in 1:len.dummy)
      {
        name.paste.temp<-paste(cov.names[i],haplo.names,sep=".")
        if(i==1){name.paste=name.paste.temp
        }else{name.paste=c(name.paste,name.paste.temp)}
      }
      pname<-c("beta0",haplo.names,name.paste,cov.names)
      names(BF)<-pname
      names(OR)<-pname

      ci.OR<-data.frame(pname,ci.OR)
      colnames(ci.OR)<-c("Name","Lower","Upper")

      rownames(post.mean.freq)=c(haplo.names,haplo.baseline)
      rownames(ci.freq)=c(haplo.names,haplo.baseline)
      names.step1=rep(dep.names0,len.dep)[as.logical(as.vector(t(value.dep0)))]
      names.step2=matrix(names.step1,nrow=len.dep,byrow=TRUE)
      names.step3=apply(names.step2, 1, function(x) paste(x, collapse=":"))
      colnames(post.mean.freq)=names.step3
      colnames(ci.freq)=names.step3

      rownames(ci.b)=haplo.names
      colnames(ci.b)=paste("gamma",0:len.dummy.dep,sep="")

      names(ci.lambda)<-c("Lower","Upper")
      names(ci.D)<-c("Lower","Upper")

      ans <- list(BF=BF, OR=OR, CI.OR=ci.OR, freq=post.mean.freq, CI.freq=ci.freq, CI.gamma=ci.b, CI.lambda=ci.lambda, CI.D=ci.D)
    }#ends here

    ##################################################
    #LBL-GXE
    ##################################################
    if(n.cov>=1&&interaction.model=="u"){

      names.cov=colnames(haplos.list$nonHaploDM)[2:(1+n.cov)]
      if (names.dep=="missing") names.dep=names.cov
      temp=cbind(haplos.list$ID, haplos.list$wt, haplos.list$nonHaploDM, haplos.list$haploDM)
      temp=temp[order(temp[,1]),]
      haplos.list$ID=temp[,1]
      haplos.list$wt=temp[,2]
      haplos.list$nonHaploDM=temp[,3:(3+n.cov)]
      haplos.list$haploDM=temp[(3+n.cov+1):ncol(temp)]

      if (is.null(seed)==FALSE) set.seed(seed)
      if (is.null(burn.in)==TRUE) burn.in=20000
      if (is.null(num.it)==TRUE) num.it=100000
      freq <- haplos.list$initFreq
      freq = freq[names(haplos.list$haploDM)]
      freq.names <- names(freq)
      if (haplo.baseline=="missing") haplo.baseline <- freq.names[which.max(freq)]
      column.subset <- colnames(haplos.list$haploDM) != haplo.baseline
      haplo.mat=haplos.list$haploDM[,column.subset]
      haplo.names=colnames(haplo.mat)

      freq.new<-freq[column.subset]
      freq.new[length(freq.new)+1]<-freq[haplo.baseline]

      if (cov.baseline=="missing") cov.baseline = rep(0, n.cov)
      haplos.list$nonHaploDM=as.matrix(haplos.list$nonHaploDM)
      num.dummies.each.cov=rep(0, n.cov)
      for(i in 1:n.cov)
      {
        cov.temp=dummy(haplos.list$nonHaploDM[,i+1])
        cov.levels=sort(unique(haplos.list$nonHaploDM[,i+1]))
        cov.subset=(cov.levels!=cov.baseline[i])
        cov.temp=cov.temp[,cov.subset]
        colnames.cov.temp=paste(colnames(haplos.list$nonHaploDM)[i+1], cov.levels[cov.subset], sep="")
        num.dummies.each.cov[i]=length(cov.levels)-1
        if(i==1)
        {
          cov.mat=cov.temp
          cov.names=colnames.cov.temp
        }else{
          cov.mat=cbind(cov.mat, cov.temp)
          cov.names=c(cov.names, colnames.cov.temp)
        }
      }

      y <- haplos.list$nonHaploDM[,1]
      N <- sum(haplos.list$wt)
      n=dim(haplo.mat)[1]
      num.haplo.id=as.vector(table(haplos.list$ID))

      h.length=ncol(haplo.mat)+1
      if (is.matrix(cov.mat)==TRUE){x.length=ncol(haplo.mat)+ncol(cov.mat)+ncol(haplo.mat)*ncol(cov.mat)
      }else{x.length=ncol(haplo.mat)+1+ncol(haplo.mat)}

      haplo.map<-matrix(rep(NA,2*(dim(haplo.mat)[1])),ncol=2)
      for(i in 1:dim(haplo.mat)[1])
      {
        for(j in 1:dim(haplo.mat)[2])
        {
          if(haplo.mat[i,j]==2)
          {
            haplo.map[i,1]=j
            haplo.map[i,2]=j
            break
          }
          if(haplo.mat[i,j]==1)
          {
            if(is.na(haplo.map[i,1])==TRUE) haplo.map[i,1]=j
            else
            {
              haplo.map[i,2]=j
              break
            }
          }
        }
        if(is.na(haplo.map[i,1])==TRUE)
        {
          haplo.map[i,1]=dim(haplo.mat)[2]+1
          haplo.map[i,2]=dim(haplo.mat)[2]+1
        }
        if(is.na(haplo.map[i,2])==TRUE) haplo.map[i,2]=dim(haplo.mat)[2]+1
      }

      uniq.mat<-unique(haplo.mat)
      uniq.map<-matrix(rep(NA,2*(dim(uniq.mat)[1])),ncol=2)
      for(i in 1:dim(uniq.mat)[1])
      {
        for(j in 1:dim(uniq.mat)[2])
        {
          if(uniq.mat[i,j]==2)
          {
            uniq.map[i,1]=j
            uniq.map[i,2]=j
            break
          }
          if(uniq.mat[i,j]==1)
          {
            if(is.na(uniq.map[i,1])==TRUE) uniq.map[i,1]=j
            else
            {
              uniq.map[i,2]=j
              break
            }
          }
        }
        if(is.na(uniq.map[i,1])==TRUE)
        {
          uniq.map[i,1]=dim(uniq.mat)[2]+1
          uniq.map[i,2]=dim(uniq.mat)[2]+1
        }
        if(is.na(uniq.map[i,2])==TRUE) uniq.map[i,2]=dim(uniq.mat)[2]+1
      }

      if (is.matrix(cov.mat)==TRUE){cov.vec=apply(cov.mat, 1, function(x) paste(x, collapse=""))
      }else{cov.vec=cov.mat}
      cov.vec.new=unique(cbind(haplos.list$ID, cov.vec))[,-1]
      y.new=unique(cbind(haplos.list$ID, y))[,-1]
      num.E=as.vector(table(cov.vec.new))

      value.E=unique(cov.mat)
      if (is.matrix(value.E)==TRUE){
        value.E=value.E[order(unique(cov.vec)),]
        len.E=nrow(value.E)
        len.dummy=ncol(value.E)
      }else{
        value.E=value.E[order(unique(cov.vec))]
        len.E=length(value.E)
        len.dummy=1
      }

      if (is.matrix(value.E)==TRUE){value.E.vec=apply(value.E, 1, function(x) paste(x, collapse=""))
      }else{value.E.vec=value.E}
      index.E=match(cov.vec.new, value.E.vec)

      n.dep=length(names.dep)
      which.dep=rep(0, n.dep)
      which.dep=match(names.dep, names.cov)
      which.dep=sort(which.dep)

      which.dummy.dep=rep(0, sum(num.dummies.each.cov[which.dep]))
      len.dummy.dep=length(which.dummy.dep)
      l=0
      for (i in 1:n.dep)
      {
        step=num.dummies.each.cov[which.dep[i]]
        end.dummy.dep=sum(num.dummies.each.cov[1:which.dep[i]])
        which.dummy.dep[(l+1):(l+step)]=(end.dummy.dep+1-step):end.dummy.dep
        l=l+step
      }

      if (is.matrix(value.E)==TRUE){value.dep=value.E[,which.dummy.dep]
      }else{value.dep=value.E}
      if (is.matrix(value.dep)==TRUE){value.dep.vec=apply(value.dep, 1, function(x) paste(x, collapse=""))
      }else{value.dep.vec=value.dep}
      o=order(value.dep.vec)
      t=table(value.dep.vec)
      if (is.matrix(value.dep)==TRUE){value.dep=unique(value.dep[o,])
      }else{value.dep=unique(value.dep[o])}
      len.dep=length(t)
      index.Etodep=rep(0, len.E)
      l=0
      for (i in 1:len.dep)
      {
        index.Etodep[o[(l+1):(l+t[i])]]=i
        l=l+t[i]
      }

      if (is.matrix(cov.mat)==TRUE){dep.mat=cov.mat[,which.dummy.dep]
      }else{dep.mat=cov.mat}
      if (is.matrix(dep.mat)==TRUE){dep.vec=apply(dep.mat, 1, function(x) paste(x, collapse=""))
      }else{dep.vec=dep.mat}

      j=1
      for(i in which.dep)
      {
        dep.temp0=dummy(haplos.list$nonHaploDM[,i+1])
        dep.levels0=sort(unique(haplos.list$nonHaploDM[,i+1]))
        colnames.dep.temp0=paste(colnames(haplos.list$nonHaploDM)[i+1], dep.levels0, sep="")
        if(j==1)
        {
          dep.mat0=dep.temp0
          dep.names0=colnames.dep.temp0
        }else{
          dep.mat0=cbind(dep.mat0, dep.temp0)
          dep.names0=c(dep.names0, colnames.dep.temp0)
        }
        j=j+1
      }

      value.dep0=unique(dep.mat0)
      value.dep0=value.dep0[order(unique(dep.vec)),]

      beta=rep(start.beta, x.length+1)
      beta.out<-numeric((num.it-burn.in)*(x.length+1))
      lambda.out<-numeric(num.it-burn.in)
      freq.out<-numeric((num.it-burn.in)*h.length)
      freq.out.indep<-numeric((num.it-burn.in)*h.length)
      freq.out.dep<-numeric((num.it-burn.in)*len.dep*h.length)
      b.out.dep<-numeric((num.it-burn.in)*(h.length-1)*(len.dummy.dep+1))
      D.out<-numeric(num.it-burn.in)
      count.model.1<-0
      count.model.2<-0

      out<-.C("mcmc", as.integer(y.new), as.integer(N), as.integer(n), as.integer(num.haplo.id), as.integer(x.length), as.integer(h.length), as.integer(haplo.map), as.integer(dim(uniq.mat)[1]), as.integer(t(uniq.map)), as.integer(index.E), as.integer(t(value.E)), as.integer(num.E), as.integer(len.E), as.integer(len.dummy), as.integer(index.Etodep), as.integer(t(value.dep)), as.integer(len.dep), as.integer(len.dummy.dep), as.double(beta), as.double(lambda), as.double(freq.new), as.double(D), as.double(a), as.double(b), beta.out=as.double(beta.out), lambda.out=as.double(lambda.out), freq.out=as.double(freq.out), freq.out.indep=as.double(freq.out.indep), freq.out.dep=as.double(freq.out.dep), b.out.dep=as.double(b.out.dep), D.out=as.double(D.out), as.integer(num.it), as.integer(burn.in), count.model.1=as.integer(count.model.1), count.model.2=as.integer(count.model.2))

      beta.out<-matrix(out$beta.out,nrow=num.it-burn.in, byrow=TRUE)
      freq.out<-matrix(out$freq.out,nrow=num.it-burn.in, byrow=TRUE)
      freq.out.indep<-matrix(out$freq.out.indep,nrow=num.it-burn.in, byrow=TRUE)
      freq.out.dep<-matrix(out$freq.out.dep,nrow=num.it-burn.in, byrow=TRUE)
      b.out.dep<-matrix(out$b.out.dep,nrow=num.it-burn.in, byrow=TRUE)

      post.mean.beta<-numeric(x.length+1)
      ci.beta<-numeric((x.length+1)*2)
      post.mean.freq<-numeric(h.length)
      post.mean.freq.indep<-numeric(h.length)
      post.mean.freq.dep<-numeric(len.dep*h.length)
      ci.freq<-numeric(h.length*2)
      ci.freq.indep<-numeric(h.length*2)
      ci.freq.dep<-numeric(len.dep*h.length*2)
      ci.b.dep<-numeric((h.length-1)*(len.dummy.dep+1)*2)
      ci.lambda<-numeric(2)
      ci.D<-numeric(2)

      k<-1
      for (i in 1:(x.length+1))
      {
        ci.beta[k]<-quantile(beta.out[,i], probs=0.025)
        ci.beta[k+1]<-quantile(beta.out[,i], probs=0.975)
        k<-k+2
        post.mean.beta[i]<-mean(beta.out[,i])
      }
      ci.beta<-matrix(ci.beta,nrow=x.length+1, ncol=2, byrow=TRUE)
      ci.OR<-data.frame(exp(ci.beta))
      OR<-exp(post.mean.beta)

      k<-1
      sub=apply(freq.out.indep, 1, function(x) any(x!=0))
      freq.out.indep=freq.out.indep[sub,]
      for (i in 1:(h.length))
      {
        ci.freq[k]<-quantile(freq.out[,i], probs=0.025)
        ci.freq[k+1]<-quantile(freq.out[,i], probs=0.975)
        ci.freq.indep[k]<-quantile(freq.out.indep[,i], probs=0.025)
        ci.freq.indep[k+1]<-quantile(freq.out.indep[,i], probs=0.975)
        k<-k+2
        post.mean.freq[i]<-mean(freq.out[,i])
        post.mean.freq.indep[i]<-mean(freq.out.indep[,i])
      }
      k<-1
      sub=apply(freq.out.dep, 1, function(x) any(x!=0))
      freq.out.dep=freq.out.dep[sub,]
      for (i in 1:(len.dep*h.length))
      {
        ci.freq.dep[k]<-quantile(freq.out.dep[,i], probs=0.025)
        ci.freq.dep[k+1]<-quantile(freq.out.dep[,i], probs=0.975)
        k<-k+2
        post.mean.freq.dep[i]<-mean(freq.out.dep[,i])
      }
      post.mean.freq.dep<-matrix(post.mean.freq.dep, nrow=h.length, byrow=FALSE)
      ci.freq.dep=round(ci.freq.dep,4)
      ci.freq.dep<-matrix(ci.freq.dep, ncol=2, byrow=TRUE)
      ci.freq.dep<-apply(ci.freq.dep, 1, function(x) paste(x, collapse=", "))
      ci.freq.dep<-paste("(", ci.freq.dep, sep="")
      ci.freq.dep<-paste(ci.freq.dep, ")", sep="")
      ci.freq.dep<-matrix(ci.freq.dep, nrow=h.length, byrow=FALSE)

      k<-1
      sub=apply(b.out.dep, 1, function(x) any(x!=0))
      b.out.dep=b.out.dep[sub,]
      for (i in 1:((h.length-1)*(len.dummy.dep+1)))
      {
        ci.b.dep[k]<-quantile(b.out.dep[,i], probs=0.025)
        ci.b.dep[k+1]<-quantile(b.out.dep[,i], probs=0.975)
        k<-k+2
      }
      ci.b.dep=round(ci.b.dep,4)
      ci.b.dep<-matrix(ci.b.dep, ncol=2, byrow=TRUE)
      ci.b.dep<-apply(ci.b.dep, 1, function(x) paste(x, collapse=", "))
      ci.b.dep<-paste("(", ci.b.dep, sep="")
      ci.b.dep<-paste(ci.b.dep, ")", sep="")
      ci.b.dep<-matrix(ci.b.dep, nrow=(h.length-1), byrow=TRUE)

      ci.lambda<-c(quantile(out$lambda.out, probs=0.025),quantile(out$lambda.out, probs=0.975))
      ci.D<-c(quantile(out$D.out, probs=0.025),quantile(out$D.out, probs=0.975))

      prob.alt<-numeric(x.length+1)
      BF<-numeric(x.length+1)

      for (i in 1:(x.length+1))
      {
        prob.alt[i]<-length(beta.out[,i][abs(beta.out[,i]) > e])/(num.it-burn.in)
      }

      for (i in 1:(x.length+1))
      {
        prior.prob<-(b/(e+b))^a
        prior.odds<-prior.prob/(1-prior.prob)
        if (prob.alt[i]<=(100*prior.odds)/(100*prior.odds+1))
        {
          BF[i]<-round((prob.alt[i]/(1-prob.alt[i]))/prior.odds,4)
        } else
          BF[i]<-">100"
      }

      OR<-round(OR,4)
      ci.OR<-round(ci.OR,4)
      post.mean.freq<-round(post.mean.freq,4)
      post.mean.freq.indep<-round(post.mean.freq.indep,4)
      post.mean.freq.dep<-round(post.mean.freq.dep,4)
      ci.freq<-round(ci.freq,4)
      ci.freq.indep<-round(ci.freq.indep,4)
      ci.lambda<-round(ci.lambda,4)
      ci.D<-round(ci.D,4)

      for(i in 1:len.dummy)
      {
        name.paste.temp<-paste(cov.names[i],haplo.names,sep=".")
        if(i==1){name.paste=name.paste.temp
        }else{name.paste=c(name.paste,name.paste.temp)}
      }
      pname<-c("beta0",haplo.names,name.paste,cov.names)
      names(BF)<-pname
      names(OR)<-pname

      ci.OR<-data.frame(pname,ci.OR)
      colnames(ci.OR)<-c("Name","Lower","Upper")

      ci.freq<-matrix(ci.freq,nrow=(h.length),ncol=2,byrow=TRUE)
      ci.freq.indep<-matrix(ci.freq.indep,nrow=(h.length),ncol=2,byrow=TRUE)
      ci.freq<-data.frame(c(haplo.names,haplo.baseline),ci.freq)
      ci.freq.indep<-data.frame(c(haplo.names,haplo.baseline),ci.freq.indep)
      colnames(ci.freq)<-c("Name","Lower","Upper")
      colnames(ci.freq.indep)<-c("Name","Lower","Upper")
      names(post.mean.freq)<-c(haplo.names,haplo.baseline)
      names(post.mean.freq.indep)<-c(haplo.names,haplo.baseline)

      rownames(post.mean.freq.dep)=c(haplo.names,haplo.baseline)
      rownames(ci.freq.dep)=c(haplo.names,haplo.baseline)
      names.step1=rep(dep.names0,len.dep)[as.logical(as.vector(t(value.dep0)))]
      names.step2=matrix(names.step1,nrow=len.dep,byrow=TRUE)
      names.step3=apply(names.step2, 1, function(x) paste(x, collapse=":"))
      colnames(post.mean.freq.dep)=names.step3
      colnames(ci.freq.dep)=names.step3

      rownames(ci.b.dep)=haplo.names
      colnames(ci.b.dep)=paste("gamma",0:len.dummy.dep,sep="")

      names(ci.lambda)<-c("Lower","Upper")
      names(ci.D)<-c("Lower","Upper")

      percentage.indep=out$count.model.1/(num.it-burn.in)
      percentage.dep=out$count.model.2/(num.it-burn.in)

      ans <- list(BF=BF, OR=OR, CI.OR=ci.OR, freq=post.mean.freq, CI.freq=ci.freq, CI.gamma=ci.b.dep, CI.lambda=ci.lambda, CI.D=ci.D, percentage.indep=percentage.indep, percentage.dep=percentage.dep)
    }#ends here
  }

  ####################################################################################################
  #Complex sampling methods
  ####################################################################################################
  if((twoBinaryPheno == FALSE) & (complex.sampling==TRUE))
  {
    var.baseline = cov.baseline
    interaction.cov = interaction.env
    if (is.null(n.stra)==TRUE) print("Warning: n.stra is missing!")

    weight=haplos.list$nonHaploDM[,2]
    haplos.list$nonHaploDM=haplos.list$nonHaploDM[,-2]
    n.var=dim(haplos.list$nonHaploDM)[2]-1
    names.var=colnames(haplos.list$nonHaploDM)[2:(1+n.var)]

    if (n.var==n.stra) interaction.cov=F

    if (names.dep=="missing") names.dep=names.var
    temp=cbind(haplos.list$ID, haplos.list$wt, weight, haplos.list$nonHaploDM, haplos.list$haploDM)
    temp=temp[order(temp[,1]),]
    haplos.list$ID=temp[,1]
    haplos.list$wt=temp[,2]
    weight=temp[,3]
    haplos.list$nonHaploDM=temp[,4:(4+n.var)]
    haplos.list$haploDM=temp[(4+n.var+1):ncol(temp)]

    if (is.null(seed)==FALSE) set.seed(seed)
    if (is.null(burn.in)==TRUE) burn.in=20000
    if (is.null(num.it)==TRUE) num.it=120000
    freq <- haplos.list$initFreq
    freq = freq[names(haplos.list$haploDM)]
    freq.names <- names(freq)
    if (haplo.baseline=="missing") haplo.baseline <- freq.names[which.max(freq)]
    column.subset <- colnames(haplos.list$haploDM) != haplo.baseline
    haplo.mat=haplos.list$haploDM[,column.subset]
    haplo.names=colnames(haplo.mat)

    freq.new<-freq[column.subset]
    freq.new[length(freq.new)+1]<-freq[haplo.baseline]

    if(var.baseline=="missing") var.baseline = rep(0, n.var)
    haplos.list$nonHaploDM=as.matrix(haplos.list$nonHaploDM)
    num.dummies.each.var=rep(0, n.var)
    for(i in 1:n.var)
    {
      var.temp=dummy(haplos.list$nonHaploDM[,i+1])
      var.levels=sort(unique(haplos.list$nonHaploDM[,i+1]))
      var.subset=(var.levels!=var.baseline[i])
      var.temp=var.temp[,var.subset]
      colnames.var.temp=paste(colnames(haplos.list$nonHaploDM)[i+1], var.levels[var.subset], sep=".")
      num.dummies.each.var[i]=length(var.levels)-1
      if(i==1)
      {
        var.mat=var.temp
        var.names=colnames.var.temp
      }else{
        var.mat=cbind(var.mat, var.temp)
        var.names=c(var.names, colnames.var.temp)
      }
    }

    y <- haplos.list$nonHaploDM[,1]
    N <- sum(haplos.list$wt)
    n=dim(haplo.mat)[1]
    num.haplo.id=as.vector(table(haplos.list$ID))

    h.length=ncol(haplo.mat)+1
    len.dummy.stra=sum(num.dummies.each.var[1:n.stra])
    len.dummy.cov=sum(num.dummies.each.var[-c(1:n.stra)])
    if (is.matrix(var.mat)==TRUE){
      if (interaction.stra==TRUE & interaction.cov==TRUE){
        x.length=ncol(haplo.mat)+ncol(var.mat)+ncol(haplo.mat)*ncol(var.mat)
      } else if (interaction.stra==TRUE & interaction.cov==FALSE){
        x.length=ncol(haplo.mat)+ncol(var.mat)+ncol(haplo.mat)*len.dummy.stra
      } else if (interaction.stra==FALSE & interaction.cov==TRUE){
        x.length=ncol(haplo.mat)+ncol(var.mat)+ncol(haplo.mat)*len.dummy.cov
      } else {x.length=ncol(haplo.mat)+ncol(var.mat)}
    }else{
      if (interaction.stra==TRUE){
        x.length=ncol(haplo.mat)+1+ncol(haplo.mat)
      } else {
        x.length=ncol(haplo.mat)+1}}

    haplo.map<-matrix(rep(NA,2*(dim(haplo.mat)[1])),ncol=2)
    for(i in 1:dim(haplo.mat)[1])
    {
      for(j in 1:dim(haplo.mat)[2])
      {
        if(haplo.mat[i,j]==2)
        {
          haplo.map[i,1]=j
          haplo.map[i,2]=j
          break
        }
        if(haplo.mat[i,j]==1)
        {
          if(is.na(haplo.map[i,1])==TRUE) haplo.map[i,1]=j
          else
          {
            haplo.map[i,2]=j
            break
          }
        }
      }
      if(is.na(haplo.map[i,1])==TRUE)
      {
        haplo.map[i,1]=dim(haplo.mat)[2]+1
        haplo.map[i,2]=dim(haplo.mat)[2]+1
      }
      if(is.na(haplo.map[i,2])==TRUE) haplo.map[i,2]=dim(haplo.mat)[2]+1
    }

    uniq.mat<-unique(haplo.mat)
    uniq.map<-matrix(rep(NA,2*(dim(uniq.mat)[1])),ncol=2)
    for(i in 1:dim(uniq.mat)[1])
    {
      for(j in 1:dim(uniq.mat)[2])
      {
        if(uniq.mat[i,j]==2)
        {
          uniq.map[i,1]=j
          uniq.map[i,2]=j
          break
        }
        if(uniq.mat[i,j]==1)
        {
          if(is.na(uniq.map[i,1])==TRUE) uniq.map[i,1]=j
          else
          {
            uniq.map[i,2]=j
            break
          }
        }
      }
      if(is.na(uniq.map[i,1])==TRUE)
      {
        uniq.map[i,1]=dim(uniq.mat)[2]+1
        uniq.map[i,2]=dim(uniq.mat)[2]+1
      }
      if(is.na(uniq.map[i,2])==TRUE) uniq.map[i,2]=dim(uniq.mat)[2]+1
    }

    if (is.matrix(var.mat)==TRUE){var.vec=apply(var.mat, 1, function(x) paste(x, collapse=""))
    }else{var.vec=var.mat}
    var.vec.new=unique(cbind(haplos.list$ID, var.vec))[,-1]
    y.new=unique(cbind(haplos.list$ID, y))[,-1]
    weight=unique(cbind(haplos.list$ID, weight))[,-1]
    weight0=weight[order(var.vec.new)]
    num.E=as.vector(table(var.vec.new))
    num.E.wt=rep(0,length(num.E))
    l=0
    for (i in 1:length(num.E.wt))
    {
      num.E.wt[i]=sum(weight0[(l+1):(l+num.E[i])])
      l=l+num.E[i]
    }

    value.E=unique(var.mat)
    if (is.matrix(value.E)==TRUE){
      value.E=value.E[order(unique(var.vec)),]
      len.E=nrow(value.E)
      len.dummy=ncol(value.E)
    }else{
      value.E=value.E[order(unique(var.vec))]
      len.E=length(value.E)
      len.dummy=1
    }

    if (is.matrix(value.E)==TRUE){value.E.vec=apply(value.E, 1, function(x) paste(x, collapse=""))
    }else{value.E.vec=value.E}
    index.E=match(var.vec.new, value.E.vec)

    n.dep=length(names.dep)
    which.dep=rep(0, n.dep)
    which.dep=match(names.dep, names.var)
    which.dep=sort(which.dep)

    which.dummy.dep=rep(0, sum(num.dummies.each.var[which.dep]))
    len.dummy.dep=length(which.dummy.dep)
    l=0
    for (i in 1:n.dep)
    {
      step=num.dummies.each.var[which.dep[i]]
      end.dummy.dep=sum(num.dummies.each.var[1:which.dep[i]])
      which.dummy.dep[(l+1):(l+step)]=(end.dummy.dep+1-step):end.dummy.dep
      l=l+step
    }

    if (is.matrix(value.E)==TRUE){value.dep=value.E[,which.dummy.dep]
    }else{value.dep=value.E}
    if (is.matrix(value.dep)==TRUE){value.dep.vec=apply(value.dep, 1, function(x) paste(x, collapse=""))
    }else{value.dep.vec=value.dep}
    o=order(value.dep.vec)
    t=table(value.dep.vec)
    if (is.matrix(value.dep)==TRUE){value.dep=unique(value.dep[o,])
    }else{value.dep=unique(value.dep[o])}
    len.dep=length(t)
    index.Etodep=rep(0, len.E)
    l=0
    for (i in 1:len.dep)
    {
      index.Etodep[o[(l+1):(l+t[i])]]=i
      l=l+t[i]
    }

    if (is.matrix(var.mat)==TRUE){dep.mat=var.mat[,which.dummy.dep]
    }else{dep.mat=var.mat}
    if (is.matrix(dep.mat)==TRUE){dep.vec=apply(dep.mat, 1, function(x) paste(x, collapse=""))
    }else{dep.vec=dep.mat}

    j=1
    for(i in which.dep)
    {
      dep.temp0=dummy(haplos.list$nonHaploDM[,i+1])
      dep.levels0=sort(unique(haplos.list$nonHaploDM[,i+1]))
      colnames.dep.temp0=paste(colnames(haplos.list$nonHaploDM)[i+1], dep.levels0, sep=".")
      if(j==1)
      {
        dep.mat0=dep.temp0
        dep.names0=colnames.dep.temp0
      }else{
        dep.mat0=cbind(dep.mat0, dep.temp0)
        dep.names0=c(dep.names0, colnames.dep.temp0)
      }
      j=j+1
    }

    value.dep0=unique(dep.mat0)
    value.dep0=value.dep0[order(unique(dep.vec)),]

    beta=rep(start.beta, x.length+1)
    beta.out<-numeric((num.it-burn.in)*(x.length+1))
    lambda.out<-numeric(num.it-burn.in)
    b.out<-numeric((num.it-burn.in)*(h.length-1)*(len.dummy.dep+1))
    freq.out<-numeric((num.it-burn.in)*len.dep*h.length)
    D.out<-numeric(num.it-burn.in)

    out<-.C("mcmc_complex", as.integer(y.new), as.double(weight), as.integer(N), as.integer(n), as.integer(num.haplo.id), as.integer(x.length), as.integer(h.length), as.integer(haplo.map), as.integer(dim(uniq.mat)[1]), as.integer(t(uniq.map)), as.integer(index.E), as.integer(t(value.E)), as.integer(num.E.wt), as.integer(len.E), as.integer(len.dummy), as.integer(index.Etodep), as.integer(t(value.dep)), as.integer(len.dep), as.integer(len.dummy.dep), as.double(beta), as.double(lambda), as.double(freq.new), as.double(D), as.double(a), as.double(b), beta.out=as.double(beta.out), lambda.out=as.double(lambda.out), b.out=as.double(b.out), freq.out=as.double(freq.out), D.out=as.double(D.out), as.integer(num.it), as.integer(burn.in), as.integer(interaction.stra), as.integer(interaction.cov), as.integer(len.dummy.stra), as.integer(len.dummy.cov))

    beta.out<-matrix(out$beta.out,nrow=num.it-burn.in, byrow=TRUE)
    b.out<-matrix(out$b.out,nrow=num.it-burn.in, byrow=TRUE)
    freq.out<-matrix(out$freq.out,nrow=num.it-burn.in, byrow=TRUE)

    ci.beta<-numeric((x.length+1)*2)
    post.mean.beta<-numeric(x.length+1)
    ci.beta<-numeric((x.length+1)*2)
    post.mean.freq<-numeric(len.dep*h.length)
    ci.freq<-numeric(len.dep*h.length*2)
    ci.b<-numeric((h.length-1)*(len.dummy.dep+1)*2)
    ci.lambda<-numeric(2)
    ci.D<-numeric(2)

    k<-1
    for (i in 1:(x.length+1))
    {
      ci.beta[k]<-quantile(beta.out[,i], probs=0.025)
      ci.beta[k+1]<-quantile(beta.out[,i], probs=0.975)
      k<-k+2
      post.mean.beta[i]<-mean(beta.out[,i])
    }
    ci.beta<-matrix(ci.beta,nrow=x.length+1, ncol=2, byrow=TRUE)
    ci.OR<-data.frame(exp(ci.beta))
    OR<-exp(post.mean.beta)

    k<-1
    for (i in 1:(len.dep*h.length))
    {
      ci.freq[k]<-quantile(freq.out[,i], probs=0.025)
      ci.freq[k+1]<-quantile(freq.out[,i], probs=0.975)
      k<-k+2
      post.mean.freq[i]<-mean(freq.out[,i])
    }
    post.mean.freq<-matrix(post.mean.freq, nrow=h.length, byrow=FALSE)
    ci.freq=round(ci.freq,4)
    ci.freq<-matrix(ci.freq, ncol=2, byrow=TRUE)
    ci.freq<-apply(ci.freq, 1, function(x) paste(x, collapse=", "))
    ci.freq<-paste("(", ci.freq, sep="")
    ci.freq<-paste(ci.freq, ")", sep="")
    ci.freq<-matrix(ci.freq, nrow=h.length, byrow=FALSE)

    k<-1
    for (i in 1:((h.length-1)*(len.dummy.dep+1)))
    {
      ci.b[k]<-quantile(b.out[,i], probs=0.025)
      ci.b[k+1]<-quantile(b.out[,i], probs=0.975)
      k<-k+2
    }
    ci.b=round(ci.b,4)
    ci.b<-matrix(ci.b, ncol=2, byrow=TRUE)
    ci.b<-apply(ci.b, 1, function(x) paste(x, collapse=", "))
    ci.b<-paste("(", ci.b, sep="")
    ci.b<-paste(ci.b, ")", sep="")
    ci.b<-matrix(ci.b, nrow=(h.length-1), byrow=TRUE)

    ci.lambda<-c(quantile(out$lambda.out, probs=0.025),quantile(out$lambda.out, probs=0.975))
    ci.D<-c(quantile(out$D.out, probs=0.025),quantile(out$D.out, probs=0.975))

    prob.alt<-numeric(x.length+1)
    BF<-numeric(x.length+1)

    for (i in 1:(x.length+1))
    {
      prob.alt[i]<-length(beta.out[,i][abs(beta.out[,i]) > e])/(num.it-burn.in)
    }

    for (i in 1:(x.length+1))
    {
      prior.prob<-(b/(e+b))^a
      prior.odds<-prior.prob/(1-prior.prob)
      if (prob.alt[i]<=(100*prior.odds)/(100*prior.odds+1))
      {
        BF[i]<-round((prob.alt[i]/(1-prob.alt[i]))/prior.odds,4)
      } else
        BF[i]<-">100"
    }

    OR<-round(OR,4)
    ci.OR<-round(ci.OR,4)
    post.mean.freq<-round(post.mean.freq,4)
    ci.lambda<-round(ci.lambda,4)
    ci.D<-round(ci.D,4)

    if (interaction.stra==TRUE & interaction.cov==TRUE){
      for(i in 1:len.dummy)
      {
        name.paste.temp<-paste(var.names[i],haplo.names,sep=" x ")
        if(i==1){name.paste=name.paste.temp
        }else{name.paste=c(name.paste,name.paste.temp)}
      }
      pname<-c("beta0",haplo.names,name.paste,var.names)
    } else if (interaction.stra==TRUE & interaction.cov==FALSE){
      for(i in 1:len.dummy.stra)
      {
        name.paste.temp<-paste(var.names[i],haplo.names,sep=" x ")
        if(i==1){name.paste=name.paste.temp
        }else{name.paste=c(name.paste,name.paste.temp)}
      }
      pname<-c("beta0",haplo.names,name.paste,var.names)
    } else if (interaction.stra==FALSE & interaction.cov==TRUE){
      for(i in 1:len.dummy.cov)
      {
        name.paste.temp<-paste(var.names[i+len.dummy.stra],haplo.names,sep=" x ")
        if(i==1){name.paste=name.paste.temp
        }else{name.paste=c(name.paste,name.paste.temp)}
      }
      pname<-c("beta0",haplo.names,name.paste,var.names)
    } else {pname<-c("beta0",haplo.names,var.names)}
    names(BF)<-pname
    names(OR)<-pname

    ci.OR<-data.frame(pname,ci.OR)
    colnames(ci.OR)<-c("Name","Lower","Upper")

    rownames(post.mean.freq)=c(haplo.names,haplo.baseline)
    rownames(ci.freq)=c(haplo.names,haplo.baseline)
    names.step1=rep(dep.names0,len.dep)[as.logical(as.vector(t(value.dep0)))]
    names.step2=matrix(names.step1,nrow=len.dep,byrow=TRUE)
    names.step3=apply(names.step2, 1, function(x) paste(x, collapse=":"))
    colnames(post.mean.freq)=names.step3
    colnames(ci.freq)=names.step3

    rownames(ci.b)=haplo.names
    colnames(ci.b)=paste("gamma",0:len.dummy.dep,sep="")

    names(ci.lambda)<-c("Lower","Upper")
    names(ci.D)<-c("Lower","Upper")

    ans <- list(BF=BF, OR=OR, CI.OR=ci.OR, freq=post.mean.freq, CI.freq=ci.freq, CI.gamma=ci.b, CI.lambda=ci.lambda, CI.D=ci.D)
  }#ends here

  return (ans)
}


