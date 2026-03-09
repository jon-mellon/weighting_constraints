set.seed(20260309)
suppressPackageStartupMessages(library(survey))
N_pop<-90000; p_vars<-10; target_rr<-0.25
z<-rnorm(N_pop)
X<-matrix(0L,N_pop,p_vars)
for(j in 1:p_vars){ strength<-runif(1,0.35,0.85); shift<-runif(1,-0.5,0.5); eta<-shift+strength*z+rnorm(N_pop,sd=0.8); X[,j]<-rbinom(N_pop,1,plogis(eta)) }
colnames(X)<-paste0('X',1:p_vars)
beta<-seq(0.75,0.2,length.out=p_vars)
lin<-as.vector(X%*%beta+0.2*z)
int<-uniroot(function(a) mean(plogis(a+lin))-target_rr, c(-8,3))$root
rho<-plogis(int+lin)
means<-colMeans(X)
cat('mean rr',mean(rho),'\n')
sample_n<-3200
for(r in 1:20){
  s<-sample.int(N_pop,sample_n)
  resp<-rbinom(sample_n,1,rho[s])==1
  rr<-mean(resp); ridx<-s[resp]; df<-as.data.frame(X[ridx,,drop=FALSE]); des<-svydesign(ids=~1,data=df,weights=~1)
  vars<-paste0('X',1:3); form<-as.formula(paste('~',paste(vars,collapse=' + ')))
  pop_totals<-c('(Intercept)'=nrow(df), setNames(nrow(df)*means[1:3], vars))
  std_msg<-NULL; con_msg<-NULL
  std <- withCallingHandlers(tryCatch(calibrate(des,formula=form,population=pop_totals,calfun='raking',maxit=300,epsilon=1e-8),error=function(e)e),warning=function(w){std_msg<<-conditionMessage(w); invokeRestart('muffleWarning')})
  con <- withCallingHandlers(tryCatch(calibrate(des,formula=form,population=pop_totals,calfun='raking',bounds=c(rr,Inf),maxit=300,epsilon=1e-8),error=function(e)e),warning=function(w){con_msg<<-conditionMessage(w); invokeRestart('muffleWarning')})
  if(!is.null(std_msg)||!is.null(con_msg)||inherits(std,'error')||inherits(con,'error')) cat('rep',r,'rr',round(rr,3),'std',std_msg,'con',con_msg,'errstd',inherits(std,'error'),'errcon',inherits(con,'error'),'\n')
}
