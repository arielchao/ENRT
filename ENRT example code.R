library(dplyr)
library(gmodels)
library(cfdecomp)
library(geepack)

# Simulate dataset --------------------------------------------------------

set.seed(10594)

# Parameters set to similar values observed in HPTN 037
P_M=0.6
P_R=0.5
P_y0=0.3
RR=0.3
ICC=0.15

# 150 networks with 3 participants each (1 ego + 2 network members)
k = 150
n_k = 3
N = k*n_k

data <- data.frame(matrix(nrow=N,ncol=2))
colnames(data) <- c('i','k')
data$i <- rep(1:n_k,k)
data$k <- rep(1:k, each=n_k)

# clustering parameters
var.b <- 3.29*ICC/(1-ICC)
data$b <- rep(rnorm(k,0,sqrt(var.b)),each=n_k)

data$b <- ifelse(data$b>min(-log(P_y0),(-log(P_y0)-log(RR))),min(-log(P_y0),(-log(P_y0)-log(RR))),data$b)

# if(ASpE=='RD'){
#   data$b <- ifelse(RR>1 & data$b >(1-P_y0-(P_y0*RR-P_y0)),(1-P_y0-(P_y0*RR-P_y0)),data$b)
#   data$b <- ifelse(RR>1 & data$b <(-P_y0),(-P_y0),data$b)
#   data$b <- ifelse(RR<1 & data$b >(1-P_y0),(1-P_y0),data$b)
#   data$b <- ifelse(RR<1 & data$b <(-P_y0-(P_y0*RR-P_y0)),(-P_y0-(P_y0*RR-P_y0)),data$b)
# }else{
#   data$b <- ifelse(data$b>min(-log(P_y0),(-log(P_y0)-log(RR))),min(-log(P_y0),(-log(P_y0)-log(RR))),data$b)
# }

# intervention assignment 
assign <- rbinom(k,1,P_R)

# subset data to network members only
main <- subset(data,data$i !=1)

# true spillover exposure (note: used for outcome construction but only available only in internal validation study) 
for(i in 1:nrow(main)){
  main$G[i] <- assign[main$k[i]]
}

# mismeasured spillover exposure
main$u <- runif(nrow(main))

for(i in 1:nrow(main)){
  if(main$u[i]<=P_M){
    main$k.prime[i] <- main$k[i]
  }else{
    main$k.prime[i] <- sample(c(1:k)[-main$k[i]],1)
  }
  main$G_star[i] <- assign[main$k.prime[i]]
}

# outcome
main$Y <- rbinom(nrow(main),1,(P_y0*(RR^main$G)*exp(main$b))/(1+P_y0*(RR^main$G)*exp(main$b)))

# internal validation study
main$V <- sample(c(rep(0,nrow(main)*0.9),rep(1,nrow(main)*0.1)),nrow(main))
val <- subset(main,V==1)

# true spillover should not observed in main study, unless participants are part of the internal validation study
main$G<- ifelse(main$V==1,main$G,NA)


# Obtain sensitivity (theta) and specificity (phi) from validation study --------------------------------------------------------

# misclassification model
theta <- sum(val$G==1 & val$G_star==1)/sum(val$G==1)
phi <- sum(val$G==0 & val$G_star==0)/sum(val$G==0)

# Bias correction ---------------------------------------------------------

A <- sum(main$G_star==1 & main$Y==1)
B <- sum(main$G_star==0 & main$Y==1)
C <- sum(main$G_star==1 & main$Y==0)
D <- sum(main$G_star==0 & main$Y==0)
N1 <- A+C
N0 <- B+D
m1 <- A+B
m0 <- C+D
n <- nrow(main)

# RD under misclassified exposure
rd_mis <- (A/N1)-(B/N0)
var_rd_mis <- ((A/N1)*(1-(A/N1))/N1)+((B/N0)*(1-(B/N0))/N0)
lcl_rd_mis <- rd_mis - 1.96*sqrt(var_rd_mis)
ucl_rd_mis <- rd_mis + 1.96*sqrt(var_rd_mis)
rd_mis_CI <- paste(format(round(rd_mis,2),nsmall=2), ' (', format(round(lcl_rd_mis,2),nsmall=2), ',', format(round(ucl_rd_mis,2),nsmall=2), ')')

# Corrected RD (matrix method)
rd_corrected <- (B-(phi*m1))/(N0-(phi*n)) - (A-(theta*m1))/(N1-(theta*n))

VarA <- A*(1-(A/N1))
VarB <- B*(1-(B/N0))
VarTheta <- theta*(1-theta)/(sum(val$G==1))
VarPhi <- phi*(1-phi)/(sum(val$G==0))
dfdA <- (-(1-theta)/(N1-theta*n))-(phi/(N0-phi*n))
dfdB <- ((theta)/(N1-theta*n))+((1-phi)/(N0-phi*n))
dfdTheta <- (m1*N1-A*n)/((N1-theta*n)^2)
dfdPhi <- (B*n-m1*N0)/((N0-phi*n)^2)
var_rd_cor <- VarA*(dfdA^2) + VarB*(dfdB^2) + VarTheta*(dfdTheta^2) + VarPhi*(dfdPhi^2)
lcl_rd_cor <- rd_corrected - 1.96*sqrt(var_rd_cor)
ucl_rd_cor <- rd_corrected+ 1.96*sqrt(var_rd_cor)
rd_cor_CI <- paste(format(round(rd_corrected,2),nsmall=2), ' (', format(round(lcl_rd_cor,2),nsmall=2), ',', format(round(ucl_rd_cor,2),nsmall=2), ')')

# Corrected RD (inverse matrix method)
rr <- (sum(val$Y==1&val$G==1,na.rm=T)/sum(val$G==1))/(sum(val$Y==1&val$G==0,na.rm=T)/sum(val$G==0))
rr.tilde <- (sum(val$Y==0&val$G==1,na.rm=T)/sum(val$G==1))/(sum(val$Y==0&val$G==0,na.rm=T)/sum(val$G==0))

pz <- sum(main$G_star==1)/nrow(main)
pr_g1 <- (pz-1+phi)/(theta+phi-1)
pr_g0 <- 1-pr_g1

ppv1 <- (theta*rr*pr_g1)/(((1-phi)*pr_g0)+(theta*rr*pr_g1))
ppv0 <- (theta*rr.tilde*pr_g1)/(((1-phi)*pr_g0)+(theta*rr.tilde*pr_g1))
npv1 <- (phi*(1/rr)*pr_g0)/(((1-theta)*pr_g1)+(phi*(1/rr)*pr_g0))
npv0 <- (phi*(1/rr.tilde)*pr_g0)/(((1-theta)*pr_g1)+(phi*(1/rr.tilde)*pr_g0))

a.hat.invmatrix <- (ppv1*A)+((1-npv1)*B)
b.hat.invmatrix <- ((1-ppv1)*A)+(npv1*B)
c.hat.invmatrix <- (ppv0*C)+((1-npv0)*D)
d.hat.invmatrix <- ((1-ppv0)*C)+(npv0*D)

rd_corrected_inv <- (a.hat.invmatrix/(a.hat.invmatrix+c.hat.invmatrix))-(b.hat.invmatrix/(b.hat.invmatrix+d.hat.invmatrix))

num1 <- ppv1*A+(1-npv1)*B
num2 <- (1-ppv1)*A+npv1*B
den1 <- ppv1*A+(1-npv1)*B+ppv0*C+(1-npv0)*D
den2 <- (1-ppv1)*A+npv1*B+(1-ppv0*C)+npv0*D
dfdA <- ((ppv1/den1)-(ppv1*num1/den1^2))-(((1-ppv1)/den2)-((1-ppv1)*num2)/den2^2)
dfdB <- (((1-npv1)/den1)-((1-npv1)*num1/den1^2))-((npv1/den2)-(npv1*num2)/den2^2)
dfdC <- -(ppv0*num1/den1^2)+((1-ppv0)*num2/den2^2)
dfdD <- -((1-npv0)*num1/den1^2)+(npv0*num2/den2^2)
dfdppv1 <- ((A/den1)-(A*num1/den1^2))-((-A/den2)-(-A*num2/den2^2))
dfdppv0 <- -(C/num1/den1^2)+(-C*num2/den2^2)
dfdnpv1 <- ((-B/den1)-(-B*num1/den1^2))-((B/den2)-(B*num2/den2^2))
dfdnpv0 <- (D/num1/den1^2)+(D*num2/den2^2)
Varppv1 <- ppv1*(1-ppv1)/sum(val$Y==1&val$G_star==1)
Varppv0 <- ppv0*(1-ppv0)/sum(val$Y==0&val$G_star==1)
Varnpv1 <- npv1*(1-npv1)/sum(val$Y==1&val$G_star==0)
Varnpv0 <- npv0*(1-npv0)/sum(val$Y==0&val$G_star==0)
VarC <- C*(1-(C/N1))
VarD <- D*(1-(D/N0)) 
var_rd_inv <- VarA*(dfdA^2) + VarB*(dfdB^2) + VarC*(dfdC^2) + VarD*(dfdD^2) + Varppv1*(dfdppv1^2) + Varppv0*(dfdppv0^2) + Varnpv1*(dfdnpv1^2) + Varnpv0*(dfdnpv0^2)
lcl_rd_inv <- rd_corrected_inv - 1.96*sqrt(var_rd_inv)
ucl_rd_inv <- rd_corrected_inv + 1.96*sqrt(var_rd_inv)
rd_inv_CI <- paste(format(round(rd_corrected_inv,2),nsmall=2), ' (', format(round(lcl_rd_inv,2),nsmall=2), ',', format(round(ucl_rd_inv,2),nsmall=2), ')')

# RR under misclassified exposure
rr_mis <- (A/N1) / (B/N0)
log.rr_mis <- log(rr_mis)
var.log.rr.mis <- ((C/A)/N1)+(((N0-B)/B)/N0)
lcl_rr_mis <- exp(log.rr_mis  - 1.96*sqrt(var.log.rr.mis))
ucl_rr_mis <- exp(log.rr_mis  + 1.96*sqrt(var.log.rr.mis))
rr_mis_CI <- paste(format(round(rr_mis,2),nsmall=2), ' (', format(round(lcl_rr_mis,2),nsmall=2), ',', format(round(ucl_rr_mis,2),nsmall=2), ')')

# Corrected RR (matrix method)
rr_corrected <- ((B-(phi*m1))/(N0-(phi*n))) / ((A-(theta*m1))/(N1-(theta*n)))

dgdA <- -((1-theta-phi)*B)/((A-theta*m1)*(B-phi*m1))
dgdB <- ((1-theta-phi)*A)/((A-theta*m1)*(B-phi*m1))
dgdTheta <- ((A*nrow(main))-m1*N1)/((N1-theta*n)*(A-theta*m1))
dgdPhi <- ((B*n)-m1*N0)/((N0-phi*n)*(B-phi*m1))
var_logrr <- VarA*(dgdA^2) + VarB*(dgdB^2) + VarTheta*(dgdTheta^2) + VarPhi*(dgdPhi^2)
corrected_logrr <- log(rr_corrected)
lcl_rr <- exp(corrected_logrr  - 1.96*sqrt(var_logrr))
ucl_rr <- exp(corrected_logrr  + 1.96*sqrt(var_logrr))
rr_cor_CI <- paste(format(round(rr_corrected,2),nsmall=2), ' (', format(round(lcl_rr,2),nsmall=2), ',', format(round(ucl_rr,2),nsmall=2), ')')


# Corrected RR (inverse matrix method)
rr_corrected_invmatrix <- ((a.hat.invmatrix/(a.hat.invmatrix+c.hat.invmatrix)))/((b.hat.invmatrix/(b.hat.invmatrix+d.hat.invmatrix)))

one <- ppv1*A+(1-npv1)*B
two <- (1-ppv1)*A+npv1*B+(1-ppv0*C)+npv0*D
three<- (1-ppv1)*A+npv1*B
four <- ppv1*A+(1-npv1)*B+ppv0*C+(1-npv0)*D

dgdA <- (ppv1/one) + ((1-ppv1)/two) - ((1-ppv1)/three) - (ppv1/four)
dgdB <- ((1-npv1)/one) + (npv1/two) - (npv1/three) - ((1-npv1)/four)
dgdC <- ((1-ppv0)/two) - (ppv0/four)
dgdD <- (npv0/two) - ((1-npv0)/four)
dgdppv1 <- (A/one) + (-A/two) - (-A/three) - (A/four)
dgdppv0 <- (-C/two) - (C/four)
dgdnpv1 <- (-B/one) + (B/two) - (B/three) - (-B/four)
dgdnpv0 <- (D/two) - (-D/four)

corrected_logrr_inv <- log(rr_corrected_invmatrix)

var_rr_inv <- VarA*(dgdA^2) + VarB*(dgdB^2) + VarC*(dgdC^2) + VarD*(dgdD^2) + Varppv1*(dgdppv1^2) + Varppv0*(dgdppv0^2) + Varnpv1*(dgdnpv1^2) + Varnpv0*(dgdnpv0^2)
lcl_rr_inv <- exp(corrected_logrr_inv - 1.96*sqrt(var_rr_inv))
ucl_rr_inv <- exp(corrected_logrr_inv + 1.96*sqrt(var_rr_inv))
rr_inv_CI <- paste(format(round(corrected_logrr_inv,2),nsmall=2), ' (', format(round(lcl_rr_inv,2),nsmall=2), ',', format(round(ucl_rr_inv,2),nsmall=2), ')')

# Cluster variance - inflate matrix methods corrected RD/RR by design effect ---------------------------------------------------

main <- main[order(main$k), ]
reg.star <- geeglm(Y~G_star, data = main , id=factor (k), family = binomial(link='log'), corstr="exchangeable")
icc <- reg.star$geese$alpha
m.bar <- mean(table(main$k))*(1+(sd(table(main$k))/mean(table(main$k)))^2)
vif <- 1+(m.bar-1)*icc

var_rd_vif <- var_rd_cor*vif
lcl_rd_vif <- rd_corrected - 1.96*sqrt(var_rd_vif)
ucl_rd_vif <- rd_corrected + 1.96*sqrt(var_rd_vif)
rd_vif_CI <- paste(format(round(rd_corrected,2),nsmall=2), ' (', format(round(lcl_rd_vif,2),nsmall=2), ',', format(round(ucl_rd_vif,2),nsmall=2), ')')

var_rr_vif <- var_logrr*vif
lcl_rr_vif <- exp(corrected_logrr  - 1.96*sqrt(var_rr_vif))
ucl_rr_vif <- exp(corrected_logrr  + 1.96*sqrt(var_rr_vif))
rr_vif_CI <- paste(format(round(rr_corrected,2),nsmall=2), ' (', format(round(lcl_rr_vif,2),nsmall=2), ',', format(round(ucl_rr_vif,2),nsmall=2), ')')


# Cluster variance - bootstrap --------------------------------------------

boot.rd <- c()
boot.rr <- c()

for(j in 1:1000){
  
  boot <- cluster.resample(main,'k')
  boot.val <- cluster.resample(val,'k')

  theta <- sum(boot.val$G==1 & boot.val$G_star==1)/sum(boot.val$G==1)
  phi <- sum(boot.val$G==0 & boot.val$G_star==0)/sum(boot.val$G==0)

  #MAIN
  A <- sum(boot$G_star==1 & boot$Y==1)
  B <- sum(boot$G_star==0 & boot$Y==1)
  C <- sum(boot$G_star==1 & boot$Y==0)
  D <- sum(boot$G_star==0 & boot$Y==0)
  N1 <- A+C
  N0 <- B+D
  m1 <- A+B
  m0 <- C+D
  n <- nrow(boot)
  
  boot.rd[j] <- (B-(phi*m1))/(N0-(phi*n)) - (A-(theta*m1))/(N1-(theta*n))
  boot.rr[j] <- ((B-(phi*m1))/(N0-(phi*n))) / ((A-(theta*m1))/(N1-(theta*n)))
  

  if(theta > max(A/m1,C/m0) & phi > max(B/m1,D/m0)){
    boot.rd[j] <- (B-(phi*m1))/(N0-(phi*n)) - (A-(theta*m1))/(N1-(theta*n))
    boot.rr[j] <- ((B-(phi*m1))/(N0-(phi*n))) / ((A-(theta*m1))/(N1-(theta*n)))
  }else{
    boot.rd[j] <- NA
    boot.rr[j] <- NA
  }
}

rd_boot_CI <- paste(format(round(rd_corrected,2),nsmall=2), ' (', format(round(rd_corrected-1.96*sd(boot.rd,na.rm=T),2),nsmall=2), ',', format(round(rd_corrected+1.96*sd(boot.rd,na.rm=T),2),nsmall=2), ')') 
rr_boot_CI <- paste(format(round(rr_corrected,2),nsmall=2), ' (', format(round(exp(log(rr_corrected)-1.96*sd(boot.rr,na.rm=T)),2),nsmall=2), ',', format(round(exp(log(rr_corrected)+1.96*sd(boot.rr,na.rm=T)),2),nsmall=2), ')') 


# Cluster variance - likelihood -------------------------------------------

# in each network, a unit can be included in either 1) only main study, 2) both main and validation study
# likelihood is summed over all cluster in main study

# data: 
count <- matrix(nrow=length(unique(main$k)),ncol=12)
colnames(count) <- c('y1G_star1','y1G_star0','y0G_star1','y0G_star0','y1g1','y1g0','y0g1','y0g0','g0G_star0','g1G_star0','g0G_star1','g1G_star1')

for(i in 1:nrow(count)){
  count[i,1] <- sum(main$Y==1&main$G_star==1&main$V==0&main$k==unique(main$k)[i])
  count[i,2] <- sum(main$Y==1&main$G_star==0&main$V==0&main$k==unique(main$k)[i])
  count[i,3] <- sum(main$Y==0&main$G_star==1&main$V==0&main$k==unique(main$k)[i])
  count[i,4] <- sum(main$Y==0&main$G_star==0&main$V==0&main$k==unique(main$k)[i])
  count[i,5] <- sum(main$Y==1&main$G==1&main$V==1&main$k==unique(main$k)[i])
  count[i,6] <- sum(main$Y==1&main$G==0&main$V==1&main$k==unique(main$k)[i])
  count[i,7] <- sum(main$Y==0&main$G==1&main$V==1&main$k==unique(main$k)[i])
  count[i,8] <- sum(main$Y==0&main$G==0&main$V==1&main$k==unique(main$k)[i])
  count[i,9] <- sum(main$G==0&main$G_star==0&main$V==1&main$k==unique(main$k)[i])
  count[i,10] <- sum(main$G==1&main$G_star==0&main$V==1&main$k==unique(main$k)[i])
  count[i,11] <- sum(main$G==0&main$G_star==1&main$V==1&main$k==unique(main$k)[i])
  count[i,12] <- sum(main$G==1&main$G_star==1&main$V==1&main$k==unique(main$k)[i])
}

count <- ifelse(is.na(count)==T,0,count)

# expressions in integrals 

# Risk difference
integrand.G <- function(b){
  one <- d0+d1*ppv+b
  two <- d0+d1*(1-npv)+b
  three <- 1-d0-d1*ppv-b
  four <- 1-d0-d1*(1-npv)-b

  (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4])*exp((-(b)^2)/(2*sig))
}

integrand.F <- function(b){
  one <- d0+d1*ppv+b
  two <- d0+d1*(1-npv)+b
  three <- 1-d0-d1*ppv-b
  four <- 1-d0-d1*(1-npv)-b
  five <- d0+d1+b
  six <- d0+b
  seven <- 1-d0-d1-b
  eight <- 1-d0-b

  (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4]*five^count[k,5]*six^count[k,6]*seven^count[k,7]*eight^count[k,8])*exp((-(b)^2)/(2*sig))
}

integrand.d0.F <- function(b){
  one <- d0+d1*ppv+b
  two <- d0+d1*(1-npv)+b
  three <- 1-d0-d1*ppv-b
  four <- 1-d0-d1*(1-npv)-b
  five <- d0+d1+b
  six <- d0+b
  seven <- 1-d0-d1-b
  eight <- 1-d0-b

  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4]*five^count[k,5]*six^count[k,6]*seven^count[k,7]*eight^count[k,8])*exp((-(b)^2)/(2*sig))

  Br <- (count[k,1]/one)+(count[k,2]/two)-(count[k,3]/three)-(count[k,4]/four)+(count[k,5]/five)+(count[k,6]/six)-(count[k,7]/seven)-(count[k,8]/eight)

  Fun*Br
}

integrand.d0.G <- function(b){
  one <- d0+d1*ppv+b
  two <- d0+d1*(1-npv)+b
  three <- 1-d0-d1*ppv-b
  four <- 1-d0-d1*(1-npv)-b

  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4])*exp((-(b)^2)/(2*sig))

  Br <- (count[k,1]/one)+(count[k,2]/two)-(count[k,3]/three)-(count[k,4]/four)

  Fun*Br
}

integrand.d1.F <- function(b){
  one <- d0+d1*ppv+b
  two <- d0+d1*(1-npv)+b
  three <- 1-d0-d1*ppv-b
  four <- 1-d0-d1*(1-npv)-b
  five <- d0+d1+b
  six <- d0+b
  seven <- 1-d0-d1-b
  eight <- 1-d0-b

  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4]*five^count[k,5]*six^count[k,6]*seven^count[k,7]*eight^count[k,8])*exp((-(b)^2)/(2*sig))

  Br <- (ppv*count[k,1]/one)+((1-npv)*count[k,2]/two)-(ppv*count[k,3]/three)-((1-npv)*count[k,4]/four)+(count[k,5]/five)-(count[k,7]/seven)

  Fun*Br
}

integrand.d1.G <- function(b){
  one <- d0+d1*ppv+b
  two <- d0+d1*(1-npv)+b
  three <- 1-d0-d1*ppv-b
  four <- 1-d0-d1*(1-npv)-b

  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4])*exp((-(b)^2)/(2*sig))

  Br <- (ppv*count[k,1]/one)+((1-npv)*count[k,2]/two)-(ppv*count[k,3]/three)-((1-npv)*count[k,4]/four)

  Fun*Br
}

integrand.sig.F <- function(b){
  one <- d0+d1*ppv+b
  two <- d0+d1*(1-npv)+b
  three <- 1-d0-d1*ppv-b
  four <- 1-d0-d1*(1-npv)-b
  five <- d0+d1+b
  six <- d0+b
  seven <- 1-d0-d1-b
  eight <- 1-d0-b

  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4]*five^count[k,5]*six^count[k,6]*seven^count[k,7]*eight^count[k,8])*exp((-(b)^2)/(2*sig))

  (1/(2*sig^2))*Fun*(b^2)
}

integrand.sig.G <- function(b){
  one <- d0+d1*ppv+b
  two <- d0+d1*(1-npv)+b
  three <- 1-d0-d1*ppv-b
  four <- 1-d0-d1*(1-npv)-b

  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4])*exp((-(b)^2)/(2*sig))

  (1/(2*sig^2))*(b^2)*Fun
}

integrand.npv.F <- function(b){
  one <- d0+d1*ppv+b
  two <- d0+d1*(1-npv)+b
  three <- 1-d0-d1*ppv-b
  four <- 1-d0-d1*(1-npv)-b
  five <- d0+d1+b
  six <- d0+b
  seven <- 1-d0-d1-b
  eight <- 1-d0-b

  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4]*five^count[k,5]*six^count[k,6]*seven^count[k,7]*eight^count[k,8])*exp((-(b)^2)/(2*sig))

  Br <- (-d1*count[k,2]/two)+(d1*count[k,4]/four)

  Fun*Br
}

integrand.npv.G <- function(b){
  one <- d0+d1*ppv+b
  two <- d0+d1*(1-npv)+b
  three <- 1-d0-d1*ppv-b
  four <- 1-d0-d1*(1-npv)-b

  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4])*exp((-(b)^2)/(2*sig))

  Br <- (-d1*count[k,2]/two)+(d1*count[k,4]/four)

  Fun*Br
}

integrand.ppv.F <- function(b){
  one <- d0+d1*ppv+b
  two <- d0+d1*(1-npv)+b
  three <- 1-d0-d1*ppv-b
  four <- 1-d0-d1*(1-npv)-b
  five <- d0+d1+b
  six <- d0+b
  seven <- 1-d0-d1-b
  eight <- 1-d0-b

  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4]*five^count[k,5]*six^count[k,6]*seven^count[k,7]*eight^count[k,8])*exp((-(b)^2)/(2*sig))

  Br <- (d1*count[k,1]/one)-(d1*count[k,3]/three)

  Fun*Br
}

integrand.ppv.G <- function(b){
  one <- d0+d1*ppv+b
  two <- d0+d1*(1-npv)+b
  three <- 1-d0-d1*ppv-b
  four <- 1-d0-d1*(1-npv)-b

  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4])*exp((-(b)^2)/(2*sig))

  Br <- (d1*count[k,1]/one)-(d1*count[k,3]/three)

  Fun*Br
}


fun <- function(param){
  d0 <<- param[1]
  d1 <<- param[2]
  sig <<- param[3]
  npv <<- param[4]
  ppv <<- param[5]
  
  if(d1>0){
    lower <<- -d0
    upper <<- 1-d0-d1
  }else{
    lower  <<- -d0-d1
    upper <<- 1-d0
  }
  
  l <- c()
  k <<- 1
  while(k <= length(unique(main$k))){
    l[k] <- log((1/integrate(integrand.G,lower=lower,upper=upper)$value)*integrate(integrand.F,lower=lower,upper=upper)$value)+(count[k,9]*log(npv))+(count[k,10]*log(1-npv))+(count[k,11]*log(1-ppv))+(count[k,12]*log(ppv))
    k=k+1
  }
  sum(-l)
}

# initial values set at around estimated quantities
opt <- optim(c((A-(theta*m1))/(N1-(theta*n)),rd_corrected,var.b,npv1,ppv1),fun,method='L-BFGS-B',lower=c((A-(theta*m1))/(N1-(theta*n))-0.05,rd_corrected-0.05,0.01,0.01,0.01),upper=c((A-(theta*m1))/(N1-(theta*n))+0.05,rd_corrected+0.05,0.99,0.99,0.99)) 

d0 <- opt$par[1]
d1 <- opt$par[2]
sig <- opt$par[3]
npv <- opt$par[4]
ppv <- opt$par[5]


delta=10e-8
lower <- -d0-d1
upper <- 1-d0

network.score <- matrix(nrow=length(unique(main$k)),ncol=5)
colnames(network.score) <- c('dldd0','dldd1','dldsig','dldnpv','dldppv')

if(d1>0){
  for(k in 1:length(unique(main$k))){
    network.score[k,1] <- ((-integrand.F(upper)+integrand.F(lower)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrand.G(lower)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score[k,2] <- ((-integrand.F(upper)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)+((integrand.G(upper)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
    network.score[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
    network.score[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
    
  }
}else{
  for(k in 1:length(unique(main$k))){
    network.score[k,1] <- ((-integrand.F(upper)+integrand.F(lower)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrand.G(lower)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score[k,2] <- ((integrand.F(lower)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)+((integrand.G(lower)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
    network.score[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
    network.score[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
    
  }
}

score <- matrix(c(sum(network.score[,1]),sum(network.score[,2]),sum(network.score[,3]),sum(network.score[,4]),sum(network.score[,5])),nrow=5,ncol=1)

network.score.delta1 <- matrix(nrow=length(unique(main$k)),ncol=5)
network.score.delta2 <- matrix(nrow=length(unique(main$k)),ncol=5)
network.score.delta3 <- matrix(nrow=length(unique(main$k)),ncol=5)
network.score.delta4 <- matrix(nrow=length(unique(main$k)),ncol=5)
network.score.delta5 <- matrix(nrow=length(unique(main$k)),ncol=5)

if(d1 >0){
  for(k in 1:length(unique(main$k))){
    
    d0 <- d0+delta*d0
    
    network.score.delta1[k,1] <- ((-integrand.F(upper)+integrand.F(lower)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrand.G(lower)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta1[k,2] <- ((-integrand.F(upper)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)+((integrand.G(upper)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta1[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
    network.score.delta1[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
    network.score.delta1[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
    
    d0 <- d0/(1+delta)
    d1 <- d1+delta*d1
    network.score.delta2[k,1] <- ((-integrand.F(upper)+integrand.F(lower)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrand.G(lower)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta2[k,2] <- ((-integrand.F(upper)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)+((integrand.G(upper)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta2[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
    network.score.delta2[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
    network.score.delta2[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
    
    d1 <- d1/(1+delta)
    sig <- sig+delta*sig
    network.score.delta3[k,1] <- ((-integrand.F(upper)+integrand.F(lower)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrand.G(lower)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta3[k,2] <- ((-integrand.F(upper)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)+((integrand.G(upper)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta3[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
    network.score.delta3[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
    network.score.delta3[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
    
    sig <- sig/(1+delta)
    npv <- npv+delta*npv
    network.score.delta4[k,1] <- ((-integrand.F(upper)+integrand.F(lower)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrand.G(lower)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta4[k,2] <- ((-integrand.F(upper)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)+((integrand.G(upper)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta4[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
    network.score.delta4[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
    network.score.delta4[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
    
    npv <- npv/(1+delta)
    ppv <- ppv+delta*ppv
    network.score.delta5[k,1] <- ((-integrand.F(upper)+integrand.F(lower)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrand.G(lower)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta5[k,2] <- ((-integrand.F(upper)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)+((integrand.G(upper)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta5[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
    network.score.delta5[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
    network.score.delta5[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
    
    ppv <- ppv/(1+delta)
  }
}else{
  for(k in 1:length(unique(main$k))){
    d0 <- d0+delta*d0
    
    network.score.delta1[k,1] <- ((-integrand.F(upper)+integrand.F(lower)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrand.G(lower)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta1[k,2] <- ((integrand.F(lower)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)+((integrand.G(lower)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta1[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
    network.score.delta1[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
    network.score.delta1[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
    
    
    d0 <- d0/(1+delta)
    d1 <- d1+delta*d1
    network.score.delta2[k,1] <- ((-integrand.F(upper)+integrand.F(lower)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrand.G(lower)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta2[k,2] <- ((integrand.F(lower)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)+((integrand.G(lower)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta2[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
    network.score.delta2[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
    network.score.delta2[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
    
    d1 <- d1/(1+delta)
    sig <- sig+delta*sig
    network.score.delta3[k,1] <- ((-integrand.F(upper)+integrand.F(lower)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrand.G(lower)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta3[k,2] <- ((integrand.F(lower)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)+((integrand.G(lower)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta3[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
    network.score.delta3[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
    network.score.delta3[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
    
    sig <- sig/(1+delta)
    npv <- npv+delta*npv
    network.score.delta4[k,1] <- ((-integrand.F(upper)+integrand.F(lower)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrand.G(lower)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta4[k,2] <- ((integrand.F(lower)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)+((integrand.G(lower)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta4[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
    network.score.delta4[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
    network.score.delta4[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
    
    npv <- npv/(1+delta)
    ppv <- ppv+delta*ppv
    network.score.delta5[k,1] <- ((-integrand.F(upper)+integrand.F(lower)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrand.G(lower)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta5[k,2] <- ((integrand.F(lower)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)+((integrand.G(lower)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
    network.score.delta5[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
    network.score.delta5[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
    network.score.delta5[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
    
    ppv <- ppv/(1+delta)
  }
}

score.delta1 <- matrix(c(sum(network.score.delta1[,1]),sum(network.score.delta1[,2]),sum(network.score.delta1[,3]),sum(network.score.delta1[,4]),sum(network.score.delta1[,5])),nrow=5,ncol=1)
score.delta2 <- matrix(c(sum(network.score.delta2[,1]),sum(network.score.delta2[,2]),sum(network.score.delta2[,3]),sum(network.score.delta2[,4]),sum(network.score.delta2[,5])),nrow=5,ncol=1)
score.delta3 <- matrix(c(sum(network.score.delta3[,1]),sum(network.score.delta3[,2]),sum(network.score.delta3[,3]),sum(network.score.delta3[,4]),sum(network.score.delta3[,5])),nrow=5,ncol=1)
score.delta4 <- matrix(c(sum(network.score.delta4[,1]),sum(network.score.delta4[,2]),sum(network.score.delta4[,3]),sum(network.score.delta4[,4]),sum(network.score.delta4[,5])),nrow=5,ncol=1)
score.delta5 <- matrix(c(sum(network.score.delta5[,1]),sum(network.score.delta5[,2]),sum(network.score.delta5[,3]),sum(network.score.delta5[,4]),sum(network.score.delta5[,5])),nrow=5,ncol=1)

hessian <- matrix(nrow=5,ncol=5)

for(i in 1:5){
  hessian[1,i] <- sum((network.score.delta1[,i]-network.score[,i])/(delta*d0))
  hessian[2,i] <- sum((network.score.delta2[,i]-network.score[,i])/(delta*d1))
  hessian[3,i] <- sum((network.score.delta3[,i]-network.score[,i])/(delta*sig))
  hessian[4,i] <- sum((network.score.delta4[,i]-network.score[,i])/(delta*npv))
  hessian[5,i] <- sum((network.score.delta5[,i]-network.score[,i])/(delta*ppv))
}

var_rd_lik <- solve(-hessian)[2,2] 
lcl_rd_lik <- d1 - 1.96*sqrt(var_rd_lik)
ucl_rd_lik <- d1 + 1.96*sqrt(var_rd_lik)
rd_lik_CI <- paste(format(round(d1,2),nsmall=2), ' (', format(round(lcl_rd_lik,2),nsmall=2), ',', format(round(ucl_rd_lik,2),nsmall=2), ')') 


# Risk ratio
integrand.G <- function(b){
  one <- exp(d0+b)*(1-ppv)+exp(d0+d1+b)*ppv
  two <- exp(d0+b)*npv+exp(d0+d1+b)*(1-npv)
  three <- (1-exp(d0+b))*(1-ppv)+(1-exp(d0+d1+b))*ppv
  four <- (1-exp(d0+b))*npv+(1-exp(d0+d1+b))*(1-npv)
  
  (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4])*exp((-(b)^2)/(2*sig))
}

integrand.F <- function(b){
  one <- exp(d0+b)*(1-ppv)+exp(d0+d1+b)*ppv
  two <- exp(d0+b)*npv+exp(d0+d1+b)*(1-npv)
  three <- (1-exp(d0+b))*(1-ppv)+(1-exp(d0+d1+b))*ppv
  four <- (1-exp(d0+b))*npv+(1-exp(d0+d1+b))*(1-npv)
  five <- exp(d0+d1+b)
  six <- exp(d0+b)
  seven <- 1-exp(d0+d1+b)
  eight <- 1-exp(d0+b)
  
  (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4]*five^count[k,5]*six^count[k,6]*seven^count[k,7]*eight^count[k,8])*exp((-(b)^2)/(2*sig))
}

integrand.d0.F <- function(b){
  one <- exp(d0+b)*(1-ppv)+exp(d0+d1+b)*ppv
  two <- exp(d0+b)*npv+exp(d0+d1+b)*(1-npv)
  three <- (1-exp(d0+b))*(1-ppv)+(1-exp(d0+d1+b))*ppv+0.000001
  four <- (1-exp(d0+b))*npv+(1-exp(d0+d1+b))*(1-npv)+0.000001
  five <- exp(d0+d1+b)
  six <- exp(d0+b)
  seven <- 1-exp(d0+d1+b)+0.000001
  eight <- 1-exp(d0+b)+0.000001
  
  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4]*five^count[k,5]*six^count[k,6]*seven^count[k,7]*eight^count[k,8])*exp((-(b)^2)/(2*sig))
  
  Br <- count[k,1]+count[k,2]+(count[k,3]*(-exp(d0+b)*(1-ppv)-exp(d0+d1+b)*ppv)/three)+(count[k,4]*(-exp(d0+b)*npv-exp(d0+d1+b)*(1-npv))/four)+count[k,5]+count[k,6]-(count[k,7]*exp(d0+d1+b)/seven)-(count[k,8]*exp(d0+b)/eight)
  
  Fun*Br
}

integrand.d0.G <- function(b){
  one <- exp(d0+b)*(1-ppv)+exp(d0+d1+b)*ppv
  two <- exp(d0+b)*npv+exp(d0+d1+b)*(1-npv)
  three <- (1-exp(d0+b))*(1-ppv)+(1-exp(d0+d1+b))*ppv+0.000001
  four <- (1-exp(d0+b))*npv+(1-exp(d0+d1+b))*(1-npv)+0.000001
  
  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4])*exp((-(b)^2)/(2*sig))
  
  Br <- count[k,1]+count[k,2]+(count[k,3]*(-exp(d0+b)*(1-ppv)-exp(d0+d1+b)*ppv)/three)+(count[k,4]*(-exp(d0+b)*npv-exp(d0+d1+b)*(1-npv))/four)
  
  Fun*Br
}

integrand.d1.F <- function(b){
  one <- exp(d0+b)*(1-ppv)+exp(d0+d1+b)*ppv+0.000001
  two <- exp(d0+b)*npv+exp(d0+d1+b)*(1-npv)+0.000001
  three <- (1-exp(d0+b))*(1-ppv)+(1-exp(d0+d1+b))*ppv+0.000001
  four <- (1-exp(d0+b))*npv+(1-exp(d0+d1+b))*(1-npv)+0.000001
  five <- exp(d0+d1+b)
  six <- exp(d0+b)
  seven <- 1-exp(d0+d1+b)+0.00001
  eight <- 1-exp(d0+b)
  
  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4]*five^count[k,5]*six^count[k,6]*seven^count[k,7]*eight^count[k,8])*exp((-(b)^2)/(2*sig))
  
  Br <- (count[k,1]*ppv*exp(d0+d1+b)/one)+(count[k,2]*(1-npv)*exp(d0+d1+b)/two)-(count[k,3]*ppv*exp(d0+d1+b)/three)-(count[k,4]*(1-npv)*exp(d0+d1+b)/four)+count[k,5]-(count[k,7]*exp(d0+d1+b)/seven)

  Fun*Br
}

integrand.d1.G <- function(b){
  one <- exp(d0+b)*(1-ppv)+exp(d0+d1+b)*ppv+0.000001
  two <- exp(d0+b)*npv+exp(d0+d1+b)*(1-npv)+0.000001
  three <- (1-exp(d0+b))*(1-ppv)+(1-exp(d0+d1+b))*ppv+0.000001
  four <- (1-exp(d0+b))*npv+(1-exp(d0+d1+b))*(1-npv)+0.000001
  
  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4])*exp((-(b)^2)/(2*sig))
  
  Br <- (count[k,1]*ppv*exp(d0+d1+b)/one)+(count[k,2]*(1-npv)*exp(d0+d1+b)/two)-(count[k,3]*ppv*exp(d0+d1+b)/three)-(count[k,4]*(1-npv)*exp(d0+d1+b)/four)

  Fun*Br
}


integrand.sig.F <- function(b){
  one <- exp(d0+b)*(1-ppv)+exp(d0+d1+b)*ppv
  two <- exp(d0+b)*npv+exp(d0+d1+b)*(1-npv)
  three <- (1-exp(d0+b))*(1-ppv)+(1-exp(d0+d1+b))*ppv
  four <- (1-exp(d0+b))*npv+(1-exp(d0+d1+b))*(1-npv)
  five <- exp(d0+d1+b)
  six <- exp(d0+b)
  seven <- 1-exp(d0+d1+b)
  eight <- 1-exp(d0+b)
  
  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4]*five^count[k,5]*six^count[k,6]*seven^count[k,7]*eight^count[k,8])*exp((-(b)^2)/(2*sig))
  
  (1/(2*sig^2))*Fun*(b^2)
}

integrand.sig.G <- function(b){
  one <- exp(d0+b)*(1-ppv)+exp(d0+d1+b)*ppv
  two <- exp(d0+b)*npv+exp(d0+d1+b)*(1-npv)
  three <- (1-exp(d0+b))*(1-ppv)+(1-exp(d0+d1+b))*ppv
  four <- (1-exp(d0+b))*npv+(1-exp(d0+d1+b))*(1-npv)

  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4])*exp((-(b)^2)/(2*sig))
  
  (1/(2*sig^2))*Fun*(b^2)
}

integrand.npv.F <- function(b){
  one <- exp(d0+b)*(1-ppv)+exp(d0+d1+b)*ppv
  two <- exp(d0+b)*npv+exp(d0+d1+b)*(1-npv)+0.000001
  three <- (1-exp(d0+b))*(1-ppv)+(1-exp(d0+d1+b))*ppv
  four <- (1-exp(d0+b))*npv+(1-exp(d0+d1+b))*(1-npv)+0.000001
  five <- exp(d0+d1+b)
  six <- exp(d0+b)
  seven <- 1-exp(d0+d1+b)
  eight <- 1-exp(d0+b)
  
  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4]*five^count[k,5]*six^count[k,6]*seven^count[k,7]*eight^count[k,8])*exp((-(b)^2)/(2*sig))
  
  Br <- (count[k,2]*(exp(d0+b)-exp(d0+d1+b))/two)+(count[k,4]*(-exp(d0+b)+exp(d0+d1+b))/four)

  Fun*Br
}

integrand.npv.G <- function(b){
  one <- exp(d0+b)*(1-ppv)+exp(d0+d1+b)*ppv
  two <- exp(d0+b)*npv+exp(d0+d1+b)*(1-npv)+0.000001
  three <- (1-exp(d0+b))*(1-ppv)+(1-exp(d0+d1+b))*ppv
  four <- (1-exp(d0+b))*npv+(1-exp(d0+d1+b))*(1-npv)+0.000001
  
  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4])*exp((-(b)^2)/(2*sig))
  
  Br <- (count[k,2]*(exp(d0+b)-exp(d0+d1+b))/two)+(count[k,4]*(-exp(d0+b)+exp(d0+d1+b))/four)
  
  Fun*Br
}


integrand.ppv.F <- function(b){
  one <- exp(d0+b)*(1-ppv)+exp(d0+d1+b)*ppv+0.000001
  two <- exp(d0+b)*npv+exp(d0+d1+b)*(1-npv)
  three <- (1-exp(d0+b))*(1-ppv)+(1-exp(d0+d1+b))*ppv+0.000001
  four <- (1-exp(d0+b))*npv+(1-exp(d0+d1+b))*(1-npv)
  five <- exp(d0+d1+b)
  six <- exp(d0+b)
  seven <- 1-exp(d0+d1+b)
  eight <- 1-exp(d0+b)
  
  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4]*five^count[k,5]*six^count[k,6]*seven^count[k,7]*eight^count[k,8])*exp((-(b)^2)/(2*sig))
  
  Br <- (count[k,1]*(-exp(d0+b)+exp(d0+d1+b))/one)+(count[k,3]*(-exp(d0+d1+b)+exp(d0+b))/three)
  
  Fun*Br
}

integrand.ppv.G <- function(b){
  one <- exp(d0+b)*(1-ppv)+exp(d0+d1+b)*ppv+0.000001
  two <- exp(d0+b)*npv+exp(d0+d1+b)*(1-npv)
  three <- (1-exp(d0+b))*(1-ppv)+(1-exp(d0+d1+b))*ppv+0.000001
  four <- (1-exp(d0+b))*npv+(1-exp(d0+d1+b))*(1-npv)
  
  Fun <- (one^count[k,1]*two^count[k,2]*three^count[k,3]*four^count[k,4])*exp((-(b)^2)/(2*sig))
  
  Br <- (count[k,1]*(-exp(d0+b)+exp(d0+d1+b))/one)+(count[k,3]*(-exp(d0+d1+b)+exp(d0+b))/three)
  
  Fun*Br
}

fun <- function(param){
  d0 <<- param[1]
  d1 <<- param[2]
  sig <<- param[3]
  npv <<- param[4]
  ppv <<- param[5]
  
  lower=-d0-d1
  upper <<- 1-d0
  
  l <- c()
  k <<- 1
  while(k <= length(unique(main$k))){
    l[k] <- log((1/integrate(integrand.G,lower=lower,upper=upper)$value)*integrate(integrand.F,lower=lower,upper=upper)$value)+(count[k,9]*log(npv))+(count[k,10]*log(1-npv))+(count[k,11]*log(1-ppv))+(count[k,12]*log(ppv))
    k=k+1
  }
  sum(-l)
}

opt <- optim(c(log((A-(theta*m1))/(N1-(theta*n))),log(rr_corrected),var.b,npv1,ppv1),fun,method='L-BFGS-B',lower=c(log((A-(theta*m1))/(N1-(theta*n)))-0.05,log(rr_corrected)-0.05,0.01,0.01,0.01),upper=c(log((A-(theta*m1))/(N1-(theta*n)))+0.05,log(rr_corrected)+0.05,0.99,0.99,0.99)) #initial values at estimated quantities

d0 <- opt$par[1]
d1 <- opt$par[2]
sig <- opt$par[3]
npv <- opt$par[4]
ppv <- opt$par[5]

delta=10e-8
lower <- -d0-d1
upper <- 1-d0

network.score <- matrix(nrow=length(unique(main$k)),ncol=5)
colnames(network.score) <- c('dldd0','dldd1','dldsig','dldnpv','dldppv')

for(k in 1:length(unique(main$k))){
  network.score[k,1] <- ((-integrand.F(upper)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
  network.score[k,2] <- ((-integrand.F(upper)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
  network.score[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
  network.score[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
  network.score[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
}

score <- matrix(c(sum(network.score[,1]),sum(network.score[,2]),sum(network.score[,3]),sum(network.score[,4]),sum(network.score[,5])),nrow=5,ncol=1)

network.score.delta1 <- matrix(nrow=length(unique(main$k)),ncol=5)
network.score.delta2 <- matrix(nrow=length(unique(main$k)),ncol=5)
network.score.delta3 <- matrix(nrow=length(unique(main$k)),ncol=5)
network.score.delta4 <- matrix(nrow=length(unique(main$k)),ncol=5)
network.score.delta5 <- matrix(nrow=length(unique(main$k)),ncol=5)

for(k in 1:length(unique(main$k))){
  d0 <- d0+delta*d0
  network.score.delta1[k,1] <- ((-integrand.F(upper)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
  network.score.delta1[k,2] <- ((-integrand.F(upper)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
  network.score.delta1[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
  network.score.delta1[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
  network.score.delta1[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
  
  d0 <- d0/(1+delta)
  d1 <- d1+delta*d1
  network.score.delta2[k,1] <- ((-integrand.F(upper)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
  network.score.delta2[k,2] <- ((-integrand.F(upper)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
  network.score.delta2[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
  network.score.delta2[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
  network.score.delta2[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
  
  d1 <- d1/(1+delta)
  sig <- sig+delta*sig
  network.score.delta3[k,1] <- ((-integrand.F(upper)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
  network.score.delta3[k,2] <- ((-integrand.F(upper)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
  network.score.delta3[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
  network.score.delta3[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
  network.score.delta3[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
  
  sig <- sig/(1+delta)
  npv <- npv+delta*npv
  network.score.delta4[k,1] <- ((-integrand.F(upper)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
  network.score.delta4[k,2] <- ((-integrand.F(upper)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
  network.score.delta4[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
  network.score.delta4[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
  network.score.delta4[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
  
  npv <- npv/(1+delta)
  ppv <- ppv+delta*ppv
  network.score.delta5[k,1] <- ((-integrand.F(upper)+integrate(integrand.d0.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrate(integrand.d0.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
  network.score.delta5[k,2] <- ((-integrand.F(upper)+integrate(integrand.d1.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((-integrand.G(upper)+integrate(integrand.d1.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value)
  network.score.delta5[k,3] <- ((integrate(integrand.sig.F,lower=lower,upper=upper)$value)/(integrate(integrand.F,lower=lower,upper=upper)$value))-((integrate(integrand.sig.G,lower=lower,upper=upper)$value)/(integrate(integrand.G,lower=lower,upper=upper)$value))
  network.score.delta5[k,4] <- (((integrate(integrand.npv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.npv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))+(count[k,9]/npv)-(count[k,10]/(1-npv))
  network.score.delta5[k,5] <- (((integrate(integrand.ppv.F,lower=lower,upper=upper)$value)/integrate(integrand.F,lower=lower,upper=upper)$value)-((integrate(integrand.ppv.G,lower=lower,upper=upper)$value)/integrate(integrand.G,lower=lower,upper=upper)$value))-(count[k,11]/(1-ppv))+(count[k,12]/ppv)
  
  ppv <- ppv/(1+delta)
}

hessian <- matrix(nrow=5,ncol=5)

for(i in 1:5){
  hessian[1,i] <- sum((network.score.delta1[,i]-network.score[,i])/(delta*d0))
  hessian[2,i] <- sum((network.score.delta2[,i]-network.score[,i])/(delta*d1))
  hessian[3,i] <- sum((network.score.delta3[,i]-network.score[,i])/(delta*sig))
  hessian[4,i] <- sum((network.score.delta4[,i]-network.score[,i])/(delta*npv))
  hessian[5,i] <- sum((network.score.delta5[,i]-network.score[,i])/(delta*ppv))
}

var_rr_lik <- solve(-hessian)[2,2] 
lcl_rr_lik <- exp(d1 - 1.96*sqrt(var_rr_lik))
ucl_rr_lik <- exp(d1 + 1.96*sqrt(var_rr_lik))
rr_lik_CI <- paste(format(round(exp(d1),2),nsmall=2), ' (', format(round(lcl_rr_lik,2),nsmall=2), ',', format(round(ucl_rr_lik,2),nsmall=2), ')') 

# All results
all_results <- matrix(c(rd_mis_CI,rd_cor_CI,rd_inv_CI,rd_vif_CI,rd_boot_CI,rd_lik_CI,
                           rr_mis_CI,rr_cor_CI,rr_inv_CI,rr_vif_CI,rr_boot_CI,rr_lik_CI),ncol=2)
rownames(all_results) <- c('ITT','Matrix method','Inverse matrix method','Design effect','Bootstrap','MLE')
colnames(all_results) <- c('RD','RR')
all_results

# Note: as noted in the manuscript, inverse matrix method results may not be accurate with small validation study size.

