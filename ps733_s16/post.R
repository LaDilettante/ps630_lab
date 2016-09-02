###### post interface ######

post <- function(model,x1name,x1vals,x2name=NULL,x2vals=NULL,holds=NULL,
                 n.sims=1000,cut=NULL,quantiles=c(.025,.975),weights=NULL){
  
  sims <- sim(model, n.sims=n.sims)
  
  if (class(model)[1] == "glm"){post.glm(model,x1name,x1vals,sims,quantiles=quantiles,
                                         x2name=x2name,x2vals=x2vals,holds=holds,weights=weights)}
  
  
  else if (class(model)[1] == "lm"){post.lm(model,x1name,x1vals,sims,quantiles=quantiles,
                                            x2name=x2name,x2vals=x2vals,holds=holds,weights=weights)}
  
  
  else if (class(model)[1] == "polr"){post.polr(model,x1name,x1vals,sims,cut=cut,quantiles=quantiles,
                                                x2name=x2name,x2vals=x2vals,holds=holds,weights=weights)}
  
  else{return("Model class not supported")}
}

####### post.glm #######

post.glm <- function(model,x1name,x1vals,sims,quantiles=c(.025,.975),
                     x2name=NULL,x2vals=NULL,holds=NULL,weights=NULL){
  
  if (is.null(weights)){wi <- c(rep(1, length(model$model[,1])))} else{wi <- weights}
  
  if (model$family[2]=="probit"){link <- pnorm}
    else if (model$family[2]=="logit"){link <- plogis}
    else if (model$family[2]=="log"){link <- exp}
    else if (model$family[2]=="cloglog"){link <- function(x){1-exp(-exp(x))}}
    else if (model$family[2]=="identity"){link <- identity}
    else {stop("Link function is missing or not supported")}
  
  n.sims <- length(coef(sims)[ ,1])
  n.x1 <- length(x1vals)
  n.obs <- length(model.matrix(lm(model$formula, model$model))[,1])
  k <- length(model.matrix(lm(model$formula, model$model))[1,])
  n.q <- length(quantiles)
  
  if (is.null(x2name)){

    X <- array(NA, c(n.obs,k,n.x1))
    
    for (i in 1:(n.x1)){
      newdata <- data.frame(model$model)
      if (!is.null(holds)){
        for (j in 1:length(holds)){
          newdata[ ,names(holds)[j]] <- as.numeric(holds[j])
        }
      }
      newdata[ ,x1name] <- x1vals[i]
      X[ , ,i] <- model.matrix(lm(model$formula, data=newdata))
    }
    
    X <- aperm(X, c(2,1,3))
    l1 <- apply(apply(X, c(2,3), function(x) link(coef(sims) %*% x)), c(1,3), function(x) weighted.mean(x, wi))
    l2 <- array(NA, c(n.x1+1,n.q+1))
    for (i in 1:n.x1){
      l2[i,1] <- mean(l1[,i])
      l2[i,2:(n.q+1)] <- quantile(l1[,i], probs=quantiles)
    }
    l2[nrow(l2),1] <- mean(l1[ ,n.x1] - l1[ ,1])
    l2[nrow(l2),2:(n.q+1)] <- quantile(l1[ ,n.x1] - l1[ ,1], probs=quantiles)
    rownames(l2) <- c(paste(c(rep(paste(x1name," ="),n.x1),
                              paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),
                            c(x1vals,"")))  
    colnames(l2) <- c("mean",quantiles)
    out <- list("est"=l2, "sims"=l1)
    invisible(out)
  } 
  
  else{

    n.x2 <- length(x2vals)
    X <- array(NA, c(n.obs,k,n.x1,n.x2))
    
    for (j in 1:n.x2){
      for (i in 1:n.x1){
        newdata <- data.frame(model$model)
        if (!is.null(holds)){
          for (k in 1:length(holds)){
            newdata[ ,names(holds)[k]] <- as.numeric(holds[k])
          }
        }
        newdata[ ,x1name] <- x1vals[i]
        newdata[ ,x2name] <- x2vals[j]
        X[ , ,i,j] <- model.matrix(lm(model$formula, data=newdata))
      }
    }
    
    X <- aperm(X, c(2,1,3,4))
    l1 <- apply(apply(X, c(2,3,4), function(x) link(coef(sims) %*% x)), c(1,3,4), function(x) weighted.mean(x, wi))
    l2 <- array(NA, c(n.x1+1,n.q+1,n.x2))
    for (j in 1:n.x2){
      for (i in 1:n.x1){
        l2[i,1,j] <- mean(l1[,i,j])
        l2[i,2:(n.q+1),j] <- quantile(l1[,i,j], probs=quantiles)
      }
      l2[nrow(l2),1,j] <- mean(l1[ ,n.x1,j] - l1[ ,1,j])
      l2[nrow(l2),2:(n.q+1),j] <- quantile(l1[ ,n.x1,j] - l1[ ,1,j], probs=quantiles)
    }
    dimnames(l2) <- list(paste(c(rep(paste(x1name," ="),n.x1),paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),c(x1vals,"")),
                         c("mean",quantiles),
                         paste(c(rep(paste(x2name," ="),n.x2)),
                               c(x2vals)))
    out <- list("est"=l2, "sims"=l1)
    invisible(out)
  }
}

####### post.lm #######

post.lm <- function(model,x1name,x1vals,sims,quantiles=c(.025,.975),
                    x2name=NULL,x2vals=NULL,holds=NULL,weights=NULL){
  
  if (is.null(weights)){wi <- c(rep(1, length(model$model[,1])))} else{wi <- weights}
  
  n.sims <- length(coef(sims)[ ,1])
  n.x1 <- length(x1vals)
  n.obs <- length(model.matrix(lm(formula(model), model$model))[,1])
  k <- length(model.matrix(lm(formula(model), model$model))[1,])
  n.q <- length(quantiles)
  
  if (is.null(x2name)){

    X <- array(NA, c(n.obs,k,n.x1))
    
    for (i in 1:(n.x1)){
      newdata <- data.frame(model$model)
      if (!is.null(holds)){
        for (j in 1:length(holds)){
          newdata[ ,names(holds)[j]] <- as.numeric(holds[j])
        }
      }
      newdata[ ,x1name] <- x1vals[i]
      X[ , ,i] <- model.matrix(lm(formula(model), data=newdata))
    }
    
    X <- aperm(X, c(2,1,3))
    l1 <- apply(apply(X, c(2,3), function(x) coef(sims) %*% x), c(1,3), function(x) weighted.mean(x, wi))
    l2 <- array(NA, c(n.x1+1,n.q+1))
    for (i in 1:n.x1){
      l2[i,1] <- mean(l1[,i])
      l2[i,2:(n.q+1)] <- quantile(l1[,i], probs=quantiles)
    }
    l2[nrow(l2),1] <- mean(l1[ ,ncol(l1)] - l1[ ,1])
    l2[nrow(l2),2:(n.q+1)] <- quantile(l1[ ,ncol(l1)] - l1[ ,1], probs=quantiles)
    rownames(l2) <- c(paste(c(rep(paste(x1name," ="),n.x1),
                              paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),
                            c(x1vals,"")))  
    colnames(l2) <- c("mean",quantiles)
    
    out <- list("est"=l2, "sims"=l1)
    invisible(out)
  } 
  
  else{
    n.sims <- length(coef(sims)[ ,1])
    n.x1 <- length(x1vals)
    n.x2 <- length(x2vals)
    n.obs <- length(model.matrix(lm(formula(model), model$model))[,1])
    k <- length(model.matrix(lm(formula(model), model$model))[1,])
    n.q <- length(quantiles)
    
    X <- array(NA, c(n.obs,k,n.x1,n.x2))
    
    for (j in 1:n.x2){
      for (i in 1:n.x1){
        newdata <- data.frame(model$model)
        if (!is.null(holds)){
          for (k in 1:length(holds)){
            newdata[ ,names(holds)[k]] <- as.numeric(holds[k])
          }
        }
        newdata[ ,x1name] <- x1vals[i]
        newdata[ ,x2name] <- x2vals[j]
        X[ , ,i,j] <- model.matrix(lm(formula(model), data=newdata))
      }
    }
    
    X <- aperm(X, c(2,1,3,4))
    l1 <- apply(apply(X, c(2,3,4), function(x) coef(sims) %*% x), c(1,3,4), function(x) weighted.mean(x, wi))
    l2 <- array(NA, c(n.x1+1,n.q+1,n.x2))
    for (j in 1:n.x2){
      for (i in 1:n.x1){
        l2[i,1,j] <- mean(l1[,i,j])
        l2[i,2:(n.q+1),j] <- quantile(l1[,i,j], probs=quantiles)
      }
      l2[nrow(l2),1,j] <- mean(l1[ ,n.x1,j] - l1[ ,1,j])
      l2[nrow(l2),2:(n.q+1),j] <- quantile(l1[ ,n.x1,j] - l1[ ,1,j], probs=quantiles)
    }
    dimnames(l2) <- list(paste(c(rep(paste(x1name," ="),n.x1),paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),c(x1vals,"")),
                         c("mean",quantiles),
                         paste(c(rep(paste(x2name," ="),n.x2)),
                               c(x2vals)))
    out <- list("est"=l2, "sims"=l1)
    invisible(out)
  }
}

####### post.polr #######

post.polr <- function(model,x1name,x1vals,sims,cut=NULL,quantiles=c(.025,.975),
                      x2name=NULL,x2vals=NULL,holds=NULL,weights=NULL){
  
  if (is.null(weights)){wi <- c(rep(1, length(model$model[,1])))} else{wi <- weights}
  
  if (model$method=="probit"){link <- pnorm}
    else if (model$method=="logistic"){link <- plogis}
    else if (model$method=="cloglog"){link <- function(x){1-exp(-exp(x))}}
  
  n.sims <- length(coef(sims)[ ,1])
  n.x1 <- length(x1vals)
  n.obs <- length(model$model[,1])
  k <- length(model.matrix(polr(getCall(model)$formula, model$model))[1,])
  n.q <- length(quantiles)
  n.y <- length(levels(model$model[,1]))
  n.z <- length(model$zeta)
  tau <- array(NA, c(n.sims,n.z+2))
  tau[,1] <- -Inf
  tau[,2:(ncol(tau)-1)] <- coef(sims)[,1:n.z]
  tau[,ncol(tau)] <- Inf
  beta <- coef(sims)[,(n.z+1):ncol(coef(sims))]
  
  if (is.null(cut)){
    if (is.null(x2name)){
      
      X_temp <- array(NA, c(n.obs,k,n.x1))
      X <- array(NA, c(n.obs,k-1,n.x1))
      
      for (i in 1:(n.x1)){
        newdata <- data.frame(model$model)
        if (!is.null(holds)){
          for (j in 1:length(holds)){
            newdata[ ,names(holds)[j]] <- as.numeric(holds[j])
          }
        }
        newdata[ ,x1name] <- x1vals[i]
        X_temp[ , ,i] <- suppressWarnings(model.matrix(polr(getCall(model)$formula, data=newdata)))
        X[ , ,i] <- X_temp[,-1,i]
      }
      
      l1 <- array(NA, c(n.sims, n.obs, n.x1, n.y))
      X <- aperm(X, c(2,1,3))
      for (z in 1:n.y){
        l1[,,,z] <- apply(X, c(2,3), function(x) (link(tau[,z+1] - beta %*% x) - link(tau[,z] - beta %*% x)))
      }
      
      l2 <- apply(l1, c(1,3,4), function(x) weighted.mean(x, wi))
      l3 <- array(NA, c(n.x1+1, n.q+1, n.y))
      for (j in 1:n.y){
        for (i in 1:n.x1){
          l3[i,1,j] <- mean(l2[,i,j])
          l3[i,2:(n.q+1),j] <- quantile(l2[,i,j], probs=quantiles)
        }
        l3[nrow(l3),1,j] <- mean(l2[ ,n.x1,j] - l2[ ,1,j])
        l3[nrow(l3),2:(n.q+1),j] <- quantile(l2[ ,n.x1,j] - l2[ ,1,j], probs=quantiles)
      }
      dimnames(l3) <- list(paste(c(rep(paste(x1name," ="),n.x1),paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),c(x1vals,"")),
                           c("mean",quantiles),
                           paste(c(rep("Y =",length(levels(model$model[,1])))),
                                 c(1:length(levels(model$model[,1]))))) 
      out <- list("est"=l3, "sims"=l2)
      invisible(out)
    }
    else{
      
      n.x2 <- length(x2vals)
      
      X_temp <- array(NA, c(n.obs,k,n.x1,n.x2))
      X <- array(NA, c(n.obs,k-1,n.x1,n.x2))
      
      for (j in 1:n.x2){
        for (i in 1:(n.x1)){
          newdata <- data.frame(model$model)
          if (!is.null(holds)){
            for (k in 1:length(holds)){
              newdata[ ,names(holds)[k]] <- as.numeric(holds[k])
            }
          }
          newdata[ ,x1name] <- x1vals[i]
          newdata[ ,x2name] <- x2vals[j]
          X_temp[ , ,i,j] <- suppressWarnings(model.matrix(polr(getCall(model)$formula, data=newdata)))
          X[ , ,i,j] <- X_temp[,-1,i,j]
        }
      }
      
      X <- aperm(X, c(2,1,3,4))
  
      l1 <- array(NA, c(n.sims, n.obs, n.x1, n.x2, n.y))
      for (z in 1:n.y){
        l1[,,,,z] <- apply(X, c(2,3,4), function(x) (link(tau[,z+1] - beta %*% x) - link(tau[,z] - beta %*% x)))
      }
      
      l2 <- apply(l1, c(1,3,4,5), function(x) weighted.mean(x, wi))
      l3 <- array(NA, c(n.x1+1, n.q+1, n.x2, n.y))
      for (k in 1:n.y){
        for (j in 1:n.x2){
          for (i in 1:n.x1){
            l3[i,1,j,k] <- mean(l2[,i,j,k])
            l3[i,2:(n.q+1),j,k] <- quantile(l2[,i,j,k], probs=quantiles)
          }
          l3[n.x1+1,1,j,k] <- mean(l2[,n.x1,j,k] - l2[,1,j,k])
          l3[n.x1+1,2:(n.q+1),j,k] <- quantile(l2[,n.x1,j,k] - l2[,1,j,k], probs=quantiles)
        }
      }
      dimnames(l3) <- list(paste(c(rep(paste(x1name," ="),n.x1),paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),c(x1vals,"")),
                           c("mean",quantiles),
                           paste(c(rep(paste(x2name," ="),n.x2)),x2vals),
                           paste(c(rep("Y =",n.y)), c(1:n.y))) 
      out <- list("est"=l3, "sims"=l2)
      invisible(out)
    }
  }
  else{
    if (is.null(x2name)){
      
      X_temp <- array(NA, c(n.obs,k,n.x1))
      X <- array(NA, c(n.obs,k-1,n.x1))
      
      for (i in 1:(n.x1)){
        newdata <- data.frame(model$model)
        if (!is.null(holds)){
          for (j in 1:length(holds)){
            newdata[ ,names(holds)[j]] <- as.numeric(holds[j])
          }
        }
        newdata[ ,x1name] <- x1vals[i]
        X_temp[ , ,i] <- suppressWarnings(model.matrix(polr(getCall(model)$formula, data=newdata)))
        X[ , ,i] <- X_temp[,-1,i]
      }
      
      X <- aperm(X, c(2,1,3))
      l1 <- apply(apply(X, c(2,3), function(x) link(-tau[,cut+1] + beta %*% x)), 
                  c(1,3), function(x) weighted.mean(x, wi))
      l2 <- array(NA, c(n.x1+1,n.q+1))
      for (i in 1:n.x1){
        l2[i,1] <- mean(l1[,i])
        l2[i,2:(n.q+1)] <- quantile(l1[,i], probs=quantiles)
      }
      l2[nrow(l2),1] <- mean(l1[ ,ncol(l1)] - l1[ ,1])
      l2[nrow(l2),2:(n.q+1)] <- quantile(l1[ ,ncol(l1)] - l1[ ,1], probs=quantiles)
      rownames(l2) <- c(paste(c(rep(paste(x1name," ="),n.x1),
                                paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),
                              c(x1vals,"")))  
      colnames(l2) <- c("mean",quantiles)
      out <- list("est"=l2, "sims"=l1)
      invisible(out)
    } 
    
    else{
      
      n.x2 <- length(x2vals)
      
      X_temp <- array(NA, c(n.obs,k,n.x1,n.x2))
      X <- array(NA, c(n.obs,k-1,n.x1,n.x2))
      
      for (j in 1:n.x2){
        for (i in 1:n.x1){
          newdata <- data.frame(model$model)
          if (!is.null(holds)){
            for (k in 1:length(holds)){
              newdata[ ,names(holds)[k]] <- as.numeric(holds[k])
            }
          }
          newdata[ ,x1name] <- x1vals[i]
          newdata[ ,x2name] <- x2vals[j]
          X_temp[ , ,i,j] <- suppressWarnings(model.matrix(polr(getCall(model)$formula, data=newdata)))
          X[ , ,i,j] <- X_temp[,-1,i,j]
        }
      }
      
      X <- aperm(X, c(2,1,3,4))
      l1 <- apply(apply(X, c(2,3,4), function(x) link(-tau[,cut+1] + beta %*% x)), c(1,3,4), function(x) weighted.mean(x, wi))
      l2 <- array(NA, c(n.x1+1,n.q+1,n.x2))
      for (j in 1:n.x2){
        for (i in 1:n.x1){
          l2[i,1,j] <- mean(l1[,i,j])
          l2[i,2:(n.q+1),j] <- quantile(l1[,i,j], probs=quantiles)
        }
        l2[nrow(l2),1,j] <- mean(l1[ ,n.x1,j] - l1[ ,1,j])
        l2[nrow(l2),2:(n.q+1),j] <- quantile(l1[ ,n.x1,j] - l1[ ,1,j], probs=quantiles)
      }
      dimnames(l2) <- list(paste(c(rep(paste(x1name," ="),n.x1),paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),c(x1vals,"")),
                           c("mean",quantiles),
                           paste(c(rep(paste(x2name," ="),n.x2)),
                                 c(x2vals)))
      out <- list("est"=l2, "sims"=l1)
      invisible(out)
    }
  }
}