library("glmnet")
library("doParallel")
library("foreach")
library("spdep")
library("ROI")
library("ROI.plugin.qpoases")
library("plotrix")
library("copula")
library("VineCopula")
library("optimParallel")
## raw data

load("~/spatialGARCH/Real_estate_data.Rdata")
load("~/spatialGARCH/Real_estate_data_2.Rdata")



plz <- c()
for (i in 1:length(list_object[[3]])) {
    plz[i] <- slot(slot(list_object[[3]][i], "polygons")[[1]], 
        "ID")
}

sum_prices <- list_object_2[[1]]
areas <- list_object_2[[2]]

T <- ncol(areas)

sum_prices[is.na(areas)] <- mean(sum_prices, na.rm = TRUE)
areas[is.na(areas)] <- mean(areas, na.rm = TRUE)

areas <- areas + 1e-09  

group_ind <- match(substr(plz, 1, 3), unique(substr(plz, 1, 3)))
n <- max(group_ind)
rownames(sum_prices) <- group_ind
rownames(areas) <- group_ind

prices <- rowsum(sum_prices, rownames(sum_prices))/rowsum(areas, 
    rownames(areas))
prices <- t(prices)

save.colnames <- colnames(prices)

prices.norm <- matrix(, nrow(prices), ncol(prices))

for (i.norm in 1:ncol(prices.norm)) {
    prices.norm[, i.norm] <- qnorm(ecdf(prices[, i.norm])(prices[, i.norm]), 0, 1)  
## get the empirical distribution function and the corresponding quantile for each price, then pretend they are normal quantiles and get the normalized value with qnorm, see ziel et al. 2017
    prices.norm[which(prices.norm[, i.norm] == Inf), i.norm] <- max(prices.norm[-which(prices.norm[, i.norm] == Inf), i.norm]) + 1e-04  
## ATTENTION if this isnt overwritten we get INF
}

prices <- prices.norm
colnames(prices) <- save.colnames


group_ind <- group_ind[-which(group_ind == as.numeric(colnames(prices)[15]))]
plzs <- unique(substr(plz, 1, 3))[-as.numeric(colnames(prices)[15])]
prices <- prices[, -15]  ## lots of 0s 
n <- n - 1


Y.mat <- prices
Y <- as.numeric(prices)

Y_quantiles=pnorm(Y.mat)

Y_quantiles==0










#now vine copula


#C-Vine structure matrix for dim 23
structure_matt_23=t(matrix(c(rep(1,23),
              0,rep(2,22),
              rep(0,2),rep(3,21),
              rep(0,3),rep(4,20),
              rep(0,4),rep(5,19),
              rep(0,5),rep(6,18),
              rep(0,6),rep(7,17),
              rep(0,7),rep(8,16),
              rep(0,8),rep(9,15),
              rep(0,9),rep(10,14),
              rep(0,10),rep(11,13),
              rep(0,11),rep(12,12),
              rep(0,12),rep(13,11),
              rep(0,13),rep(14,10),
              rep(0,14),rep(15,9),
              rep(0,15),rep(16,8),
              rep(0,16),rep(17,7),
              rep(0,17),rep(18,6),
              rep(0,18),rep(19,5),
              rep(0,19),rep(20,4),
              rep(0,20),rep(21,3),
              rep(0,21),rep(22,2),
              rep(0,22),rep(23,1))
    
            ,ncol=23,nrow=23))



#define the types of each pair copula (must be upper triagular by definition)
#here, we just take the clayton copula throughout
vine_matrix_23=-3*diag(23)+3
vine_matrix_23[lower.tri(vine_matrix_23)] <- 0

#matrix of copula parameters
#must be strictly lower triangular
lower_par_mat_23=diag(23)*(-1)+1
lower_par_mat_23[upper.tri(lower_par_mat_23)] <- 0

#set up a copula object for comprising of the structure matrix, copula type matrix and copula parameter matrix
copula_obj_23=RVineMatrix(Matrix=structure_matt_23,family=vine_matrix_23,par=lower_par_mat_23)

#define the log-likelihood
loglike=function(par){return(-RVineLogLik(data=Y_quantiles,RVM=copula_obj_23,par=matrix(par,23,23))$loglik)}
nc=4
cl<-makeCluster(nc)
setDefaultCluster(cl=cl)
clusterExport(cl,"Y_quantiles")
clusterExport(cl,"loglike")
clusterExport(cl,"copula_obj_23")
clusterEvalQ(cl, library(VineCopula))


result1=optimParallel(par = as.vector(lower_par_mat_23),fn=loglike,lower=as.vector(matrix(0.001,23,23)),upper=as.vector(matrix(28,23,23)))


#now get the correlation matrix 

cor_mat_23=
  
                     
#takes veery long (hence, slice matrix down to four dimensions)
    
    #C-Vine structure matrix for dim 4
    structure_matt_3=t(matrix(c(rep(1,3),
                                 0,rep(2,2),
                                 rep(0,2),rep(3,1)
                                 )           ,ncol=3,nrow=3))
#structure_matt_3 <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)


#define the types of each pair copula (must be upper triagular by definition)
#here, we just take the clayton copula throughout
vine_matrix_3=-4*diag(3)+4
vine_matrix_3[lower.tri(vine_matrix_3)] <- 0

#matrix of copula parameters
#must be strictly lower triangular
lower_par_mat_3=diag(3)*(-17)+17
lower_par_mat_3[upper.tri(lower_par_mat_3)] <- 0

#set up a copula object for comprising of the structure matrix, copula type matrix and copula parameter matrix
copula_obj_3=RVineMatrix(Matrix=structure_matt_3,family=vine_matrix_3,par=lower_par_mat_3)
#function to convert a vector into an upper triagular matrix

RVinePDF(c(1,1,1),copula_obj_3)



upper_matrix=function(vec){
    n=length(vec)
    m=matrix()
    return()}
#define the log-likelihood
loglike=function(par){return(-RVineLogLik(data=Y_quantiles[,c(1,2,3,4)],RVM=copula_obj_4,par=matrix(par,4,4))$loglik)}
nc=4
cl<-makeCluster(nc)
setDefaultCluster(cl=cl)
clusterExport(cl,"Y_quantiles")
clusterExport(cl,"loglike")
clusterExport(cl,"copula_obj_4")
clusterEvalQ(cl, library(VineCopula))


result=optimParallel(par = as.vector(lower_par_mat_4),fn=loglike,lower=as.vector(matrix(0.001,4,4)),upper=as.vector(matrix(28,4,4)))
#setting the estimated copula parameters
copula_obj_4_re=RVineMatrix(Matrix=structure_matt_4,family=vine_matrix_4,par=matrix(result$par,4))

#estimating the correlation matrix 
correlation_mat_4=cor(RVineSim(100000,copula_obj_4_re))
                       
save(correlation_mat_4,file="correlation_mat_4.RData")              
