library(ROI)
library(ROI.plugin.ecos)
library(cvAUC)

logbin <- function(y_var,x_var,d1_) {
  
  logbin_table <- function(x_i) {
    p<-sprintf("%3.7f", out$coefficients[x_i+1,4])
    RR<-sprintf("%3.2f", exp(out$coefficients[x_i+1,1]))
    if (exp(out$coefficients[x_i+1,1])>999) {RR<-">999"}
    l<-sprintf("%3.2f", exp(out$coefficients[x_i+1,2]))
    if (exp(out$coefficients[x_i+1,2])>999) {l<-"NA"}
    if (exp(out$coefficients[x_i+1,2])<0.001) {l<-"NA"}
    u<-sprintf("%3.2f", exp(out$coefficients[x_i+1,3]))
    if (exp(out$coefficients[x_i+1,3])>999) {u<-"NA"}
    if (exp(out$coefficients[x_i+1,3])<0.001) {u<-"NA"}
    return(paste(p," | ",RR," (",l,"-",u,")",sep=""))
  }
  
  logbin_test <- function(x_i) {
    bounds_oo <- bounds_o
    bounds_oo$lower$val[x_i+1] <- 0
    bounds_oo$upper$val[x_i] <- 0
    bounds(o) <- bounds_oo
    return(1-pchisq(-2*ROI_solve(o, solver = "ecos", control=control)$objval-out$deviance,1))
  }
  
  logbin_0 <- function(x_i) {
    bounds_oo <- bounds_o
    bounds_oo$lower$val[x_i+1] <- out$coefficients[x_i+1]
    bounds(o) <- bounds_oo
    t <- 1
    while (abs(t)>1e-08) {
      while(pchisq(-2*ROI_solve(o, solver = "ecos", control=control)$objval-out$deviance,1)<.95 & bounds_oo$lower$val[x_i+1]<12 ){
        bounds_oo$lower$val[x_i+1] <- bounds_oo$lower$val[x_i+1]+t
        bounds(o) <- bounds_oo
      }
      bounds_oo$lower$val[x_i+1] <- bounds_oo$lower$val[x_i+1]-t
      bounds(o) <- bounds_oo
      t<-t/2
    }
    u_lim <- bounds_oo$lower$val[x_i+1]
    bounds_oo <- bounds_o
    bounds_oo$upper$val[x_i] <- out$coefficients[x_i+1]
    bounds(o) <- bounds_oo
    t <- -1
    while (abs(t)>1e-08) {
      while(pchisq(-2*ROI_solve(o, solver = "ecos", control=control)$objval-out$deviance,1)<.95 & bounds_oo$upper$val[x_i]>-12 ){
        bounds_oo$upper$val[x_i] <- bounds_oo$upper$val[x_i]+t
        bounds(o) <- bounds_oo
      }
      bounds_oo$upper$val[x_i] <- bounds_oo$upper$val[x_i]-t
      bounds(o) <- bounds_oo
      t<-t/2
    }
    l_lim <- bounds_oo$upper$val[x_i]
    return(c(l_lim,u_lim))
  }
  
  out <- list()
  d1<-d1_[complete.cases(d1_[ ,c(y_var,x_var)]),]
  control <- list(tol <- 1e-08)
  y<-t(as.matrix(d1[,y_var]))
  x<-as.matrix(cbind(1,d1[,x_var]))
  y_is_0 <- y == 0L
  n_y_is_0 <- sum(y_is_0)
  o <- OP(c(y %*% x, double(n_y_is_0), rep(1, n_y_is_0)), maximum = TRUE)
  L1 <- cbind(x, matrix(0, nrow(x), 2 * n_y_is_0))
  log1exp <- function(xi, j, n_y_is_0) {
    M <- matrix(0, nrow = 6, ncol = length(xi) + 2 * n_y_is_0)
    M[1, seq_along(xi)] <- -xi
    M[3, length(xi) + j] <- -1
    M[4, length(xi) + n_y_is_0 + j] <- -1
    M[6, length(xi) + j] <- 1
    M
  }
  L2 <- mapply(log1exp, split(x[y_is_0, ], seq_len(n_y_is_0)), 
               seq_len(n_y_is_0), MoreArgs = list(n_y_is_0 = n_y_is_0), 
               SIMPLIFY = FALSE)
  rhs <- c(c(0, 1, 0), c(0, 1, 1))
  rhs <- c(rep(-1e-07, nrow(x)), rep(rhs, n_y_is_0))
  cones <- c(K_lin(nrow(x)), K_expp(2 * n_y_is_0))
  L <- do.call(rbind, c(list(L1), L2))
  constraints(o) <- C_constraint(L, cones, rhs)
  bounds_o <- V_bound(ld = -Inf, ui = c(2:ncol(x)), ub = c(rep(0, ncol(x)-1)), li = c(2:ncol(x)), lb = c(rep(0, ncol(x)-1)), nobj = length(objective(o)))
  bounds(o) <- bounds_o
  so <- ROI_solve(o, solver = "ecos", control=control)
  out$deviance.null <- -2*so$objval
  bounds_o <- V_bound(ld = -Inf, ui = c(2:ncol(x)), ub = c(rep(12, ncol(x)-1)), li = c(2:ncol(x)), lb = c(rep(-12, ncol(x)-1)), nobj = length(objective(o)))
  bounds(o) <- bounds_o
  so <- ROI_solve(o, solver = "ecos", control=control)
  out$deviance <- -2*so$objval
  out$coefficients <- head(solution(so, force = TRUE), NCOL(x))
  fit<-glm(formula = as.formula(paste(paste("`",y_var,"`"," ~ ", sep = ""),paste("`",paste(x_var, collapse = "` + `"),"`", sep = ""), sep = "")), family = binomial(link = "log"), data = d1,start=out$coefficients)
  out$glm<-fit
  names(out$coefficients)<-names(fit$coefficients)
  limits <- c(NA,NA)
  for (i in(1:length(x_var))) {limits<-rbind(limits,logbin_0(i))}
  pval <- NA
  for (i in(1:length(x_var))) {pval<-rbind(pval,logbin_test(i))}
  out$coefficients<-cbind(out$coefficients,limits,pval)
  table<-NULL
  for (i in(1:length(x_var))) {table<-rbind(table,logbin_table(i))}
  out$table<-table
  rownames(out$table)<-x_var
  out$predicted<-exp(x %*% out$coefficients[,1])
  out$AUC<-AUC(out$predicted,as.vector(y))
  return(out)
}