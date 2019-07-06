# Probability of the null hypothesis H0: z betas with |beta|<=e and nz betas with |beta|>e.

Hprob<-function(z,nz,a,b,e){

 coeff=1

 if(z>0){
   for(i in 1:z){
     c1=c(coeff,0)
     c2=c(0,coeff)
     coeff=c1+(-1)*c2
   }
 }

 prob=0
 for(i in 1:(z+1)){
   prob=prob+coeff[i]/(b+(nz+i-1)*e)^a
 }
 prob=prob*b^a
 return(prob)
}
