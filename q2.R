# X: T number of points
# M: M number of densities, here M == 3
require("mvtnorm");
require("ellipse");

vsqr <- function(v) { matrix(v) %*% t(v) }

gaussian.density.sum <- function(x.t, weight.vec, mean.vec, sigma.vec) {
    t1 <- apply(array(1:length(weight.vec)), 1, function(i) weight.vec[[i]] * dmvnorm(x.t, mean.vec[[i]], sigma.vec[[i]]));
    if (is.nan(sum(t1))) {
      stop("error");
    }
    sum(t1);
}

aposteriori <- function(i, x.t, weight.vec, mean.vec, sigma.vec) {
    t1 <- weight.vec[[i]] * dmvnorm(x.t, mean.vec[[i]], sigma.vec[[i]]);
    t2 <- gaussian.density.sum(x.t, weight.vec, mean.vec, sigma.vec);
    if(t2==0) {
      browser();
    }
    t1 / t2;
}

weight.next <- function(i, X, weight.vec, mean.vec, sigma.vec) {
    t1 <- apply(X, 1, function(x.t) aposteriori(i, x.t, weight.vec, mean.vec, sigma.vec));
    sum(t1) / nrow(X);
}

mean.next <- function(i, X, weight.vec, mean.vec, sigma.vec) {
    t1 <- apply(X, 1, function(x.t) aposteriori(i, x.t, weight.vec, mean.vec, sigma.vec) * x.t);
    t2 <- apply(X, 1, function(x.t) aposteriori(i, x.t, weight.vec, mean.vec, sigma.vec));
    if (sum(t2) == 0) {
      stop("error");
    }
    rowSums(t1) / sum(t2);
}

sigma.next <- function(i, X, weight.vec, mean.vec, sigma.vec, mean.i) {
    t1 <- apply(X, 1, function(x.t) aposteriori(i, x.t, weight.vec, mean.vec, sigma.vec) * vsqr(x.t - mean.i));
    t2 <- apply(X, 1, function(x.t) aposteriori(i, x.t, weight.vec, mean.vec, sigma.vec));
    if(sum(t2)==0) {
        browser();
        stop("error");
    }
    matrix(rowSums(t1) / sum(t2), ncol=2);
}

log.likelihood <- function(X, weight.vec, mean.vec, sigma.vec) {
  t1 <- apply(X, 1, function(x.t) gaussian.density.sum(x.t, weight.vec, mean.vec, sigma.vec));
  if (is.nan(t1[1])) {
    stop("error");
  }
  t2 <- log(t1);
  t3 <- sum(t2);
  
}

p <- function(sigma) { sigma[1,2] / sqrt(sigma[1,1]) / sqrt(sigma[2,2]); }

main <- function() {
    weight.vec <- list(1/3, 1/3, 1/3);
    mean.vec <- list(c(0.5,2), c(-2,-1), c(6.5,0));
    sigma.vec <- rep(list(diag(2)), 3);
    M <- length(weight.vec);
    threshold <- 0.1;
    
    data <- read.table("/Users/leizhao/Desktop/course/hw5/dataming/mixtureData.txt");
    old.ll <- -Inf;
    new.ll <- log.likelihood(data, weight.vec, mean.vec, sigma.vec);

    plot(data);
    while(abs(new.ll - old.ll) > threshold) {
        print(new.ll);
        print(old.ll);
        lines(ellipse(x=p(sigma.vec[[1]]), centre=mean.vec[[1]]), type="l", col="red");
        lines(ellipse(x=p(sigma.vec[[2]]), centre=mean.vec[[2]]), type="l", col="green");
        lines(ellipse(x=p(sigma.vec[[3]]), centre=mean.vec[[3]]), type="l", col="blue");
        weight.vec <- lapply(as.list(1:M), function(i) weight.next(i, data, weight.vec, mean.vec, sigma.vec));
        sigma.vec <- lapply(as.list(1:M), function(i) sigma.next(i, data, weight.vec, mean.vec, sigma.vec, mean.vec[[i]]));
        mean.vec <- lapply(as.list(1:M), function(i) mean.next(i, data, weight.vec, mean.vec, sigma.vec));
        old.ll <- new.ll;
        new.ll <- log.likelihood(data, weight.vec, mean.vec, sigma.vec);
    }    
    
}

main();
