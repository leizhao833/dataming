# (0,1),(1,2),(2,1),(2,2),(2,3),(3,2),(3,5),(4,3),(4,5),(5,4),(5,5),(6,3),(6,4),(6,5)
data <- c(0,1,1,2,2,1,2,2,2,3,3,2,3,5,4,3,4,5,5,4,5,5,6,3,6,4,6,5)
#data <- c(0,1,1,2,2,1,2,2,2,3,3,2,3,5,4,3,4,5,5,4,5,5,6,3,6,4,6,5,99,99)
#data <- c(0,1,1,2,2,1,99,99)
data <- t(matrix(data=data, nrow=2))
K <- 2;

print.cluster <- function(cluster) {
    for(i in 1:nrow(cluster)) {
        cat("("); cat(cluster[i,1]); cat(","); cat(cluster[i,2]); cat(")");
    }
}

print.cluster.list <- function(cluster.list) {
    for(i in 1:length(cluster.list)) {
        print.cluster(cluster.list[[i]]);
        cat("\n");
    }
}

distance <- function(p1, p2) {
    dist(rbind(p1,p2), method="euclidean")
}


k.mean <- function(points, k) {
    num <- nrow(points);
    group.size <- function(n, m) {
        v <- floor(n/m);
        mod <- n%%m;
        c(rep(v+1, mod), rep(v, m-mod));
    }
    group.size.vec <- group.size(num, k);
    init.clusters <- function(points, groups) {
        ret <- list();
        for(size in groups) {
            idx <- sort(sample(1:nrow(points), size));
            ret <- c(ret, list(matrix(points[idx,], ncol=2)));
            points <- matrix(points[-idx,], ncol=2);
        }
        return(ret);
    }
    cluster.list <- init.clusters(points, group.size.vec);
    changed <- TRUE;
    while(changed) {
        changed <- FALSE;
        mean2d <- function(points) {
            return (c(mean(points[,1]), mean(points[,2])));
        }
        centroids <- lapply(cluster.list, mean2d);
#        print.cluster.list(cluster.list);
#        print(centroids);
        create.buffer <- function(size) {
            buffer <- list();
            for(i in 1:size) {
                buffer <- c(buffer, list(matrix(numeric(0),0,2)))
            }
            return(buffer);
        }
        cluster.buffer <- create.buffer(length(cluster.list));
        for(cur.cluster in cluster.list) {
            if(nrow(cur.cluster) == 0) next;
            for(i in 1:nrow(cur.cluster)) {
                pt <- cur.cluster[i,];
                dists <- lapply(centroids, function(cntr) { distance(pt, cntr); });
                idx.min <- which.min(dists);
                dest.cluster <- cluster.list[[idx.min]];
                if(!identical(dest.cluster, cur.cluster)) {
                    changed <- TRUE;
                }
                dest.cluster.buf <- cluster.buffer[[idx.min]];
                dest.cluster.buf <- rbind(dest.cluster.buf, pt);
                cluster.buffer[[idx.min]] <- dest.cluster.buf;
            }
        }
        cluster.list <- cluster.buffer;
    }
    return(cluster.list);
}

comp.cost <- function(medoids, others, points) {
    cost.matrix <- numeric(0);
    for(oi in others) {
        row <- as.numeric(lapply(medoids, function(mi) { distance(points[mi,], points[oi,]); }));
        cost.matrix <- rbind(cost.matrix, row);
    }
    cost.vec <- apply(cost.matrix, 1, min);
    index.vec <- apply(cost.matrix, 1, which.min);
    list(cost=sum(cost.vec), index=index.vec);
}

k.medoids <- function(points, k) {
    num <- nrow(points);
    medoids <- sample(1:nrow(points), k);
    others <- (1:nrow(points))[-medoids];
    cost.info <- comp.cost(medoids, others, points);
    cost <- cost.info$'cost';
    final.medoids <- medoids;
    final.others <- others;
    final.assignment <- cost.info$'index';
    changed <- TRUE;
    while(changed) {
        changed <- FALSE;
        cat("\nmedoids:\n");
        print.cluster(points[medoids, ]);
        cat("\nothers:\n");
        print.cluster(points[others, ]);
        cat("\ncost:\n");
        print(cost);
        for(mi in medoids) {
            for(oi in others) {
                new.medoids <- union(setdiff(medoids, mi), oi);
                new.others <- union(setdiff(others, oi), mi);
                cost.info <- comp.cost(new.medoids, new.others, points);
                new.cost <- cost.info$'cost';
                if(new.cost < cost) {
                    cat("swap:\n");
                    print.cluster(matrix(points[mi, ], 1, 2));
                    print.cluster(matrix(points[oi, ], 1, 2));
                    changed <- TRUE;
                    final.medoids <- new.medoids;
                    final.others <- new.others;
                    final.assignment <- cost.info$'index';
                    medoids <- new.medoids;
                    others <- new.others;
                    cost <- new.cost;
                    break;
                }
            }
            if(changed) break;
        }
    }
    cluster.index.list <- as.list(final.medoids);
    for(i in 1:length(final.others)) {
        oi <- final.others[i];
        idx <- final.assignment[i];
        cl <- cluster.index.list[[idx]];
        cluster.index.list[[idx]] <- union(cl, oi);
    }
    index2matrix <- function(index.vec) {
        points[sort(index.vec),];
    }
    lapply(cluster.index.list, index2matrix);
}

agnes <- function(points, k) {
    num.clusters <- nrow(points);
    dissimilarity <- matrix(Inf, num.clusters, num.clusters);
    cluster.list <- list();
    for(i in 1:num.clusters) {
        cluster.list[[i]] <- matrix(points[i,], ncol=2);
    }
    for(i in 1:(num.clusters-1)) {
        for(j in (i+1):num.clusters) {
            dissimilarity[i,j] <- distance(points[i,], points[j,]);
        }
    }
    which.min.2d <- function(m) {
        r <- nrow(m);
        idx <- which.min(m);
        c(idx%%r, ceiling(idx/r));
    }
    while(num.clusters > k) {
        #print("--- dissimilarity ---");
        #print(dissimilarity);
        #print("--- clusters ---");
        print.cluster.list(cluster.list);
        closest.pair <- which.min.2d(dissimilarity);
        #print("--- closest pair ---");
        print(closest.pair);
        k1 <- closest.pair[1];
        k2 <- closest.pair[2];
        first.row <- rep(Inf, num.clusters-1);
        j <- 1;
        for(i in 1:num.clusters) {
            if(i==k1 || i==k2) next;
            j <- j+1;
            dist.min <- min(
                dissimilarity[k1,i],
                dissimilarity[i,k1],
                dissimilarity[k2,i],
                dissimilarity[i,k2]);
            first.row[j] <- dist.min;
        }
        dissimilarity <- dissimilarity[c(-k1,-k2), c(-k1,-k2)];
        dissimilarity <- cbind(rep(Inf, num.clusters-2), dissimilarity);
        dissimilarity <- rbind(first.row, dissimilarity);
        new.cluster <- rbind(cluster.list[[k1]], cluster.list[[k2]]);
        cluster.list <- c(list(new.cluster), cluster.list[c(-k1,-k2)]);
        num.clusters <- num.clusters - 1;
    }
    return(cluster.list);
}

#print("-----   K mean  -----");
#l <- k.mean(data, K);
#print.cluster.list(l);
print("----- K medoids -----");
l <- k.medoids(data, K);
print.cluster.list(l);
#print("-----   AGNES  -----");
#l <- agnes(data, K);
#print.cluster.list(l);
