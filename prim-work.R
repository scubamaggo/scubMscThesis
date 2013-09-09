source("prim-source.R")

################################################################
# Test performance
############
sim <- function(x1.ord, x1.ab, x2, times, alpha.peel=0.2, alpha.paste=0.1) {

    exp.p <- double(times)
    true.p <- double(times)
    results <- list()
    for (i in 1:times) {
        # training data
        x11.d <- runif(x1.ord, 0, 10)
        x11.d <- c(x11.d, runif(x1.ab, 0, 2))
        x12.d <- runif(x1.ord+x1.ab, 0, 10)

        x21.d <- runif(x2, 0, 10)
        x22.d <- runif(x2, 0, 10)

        x1.d <- cbind(x11.d, x12.d)
        x2.d <- cbind(x21.d, x22.d)

        y1.d <- c(sample(c("other", "suspicious"), x1.ord, replace=TRUE, prob=c(3/4, 1/4)),
                  sample(c("other", "suspicious"), x1.ab, replace=TRUE, prob=c(3/4, 1/4)))

        y2.d <- sample(c("other", "suspicious"), x2, TRUE, c(3/4, 1/4))

        # test data
        x11.t <- runif(x1.ord, 0, 10)
        x11.t <- c(x11.t, runif(x1.ab, 0, 2))
        x12.t <- runif(x1.ord+x1.ab, 0, 10)

        x21.t <- runif(x2, 0, 10)
        x22.t <- runif(x2, 0, 10)

        x1.t <- cbind(x11.t, x12.t)
        x2.t <- cbind(x21.t, x22.t)

        y1.t <- c(sample(c("other", "suspicious"), x1.ord, replace=TRUE, prob=c(3/4, 1/4)),
                  sample(c("other", "suspicious"), x1.ab, replace=TRUE, prob=c(3/4, 1/4)))

        y2.t <- sample(c("other", "suspicious"), x2, TRUE, c(3/4, 1/4))

        # estimator
        res.t <- prim(x1.d, x2.d, y1.d, y2.d, alpha.peel, alpha.paste, 50, 5)

        # p-value correct box
        t.p <- lrt.pois(count.fires(y1.t[x1.t[, 1] <= 2], y2.t[x2.t[, 1] <= 2], 1, 1))
        # p-value estimated box
        e.p <- lrt.pois(count.fires(y1.t[x1.t[ ,1] <= res.t$box[2,1] &
                                            x1.t[, 1] >= res.t$box[1,1] &
                                            x1.t[, 2] <= res.t$box[2,2] &
                                            x1.t[, 2] >= res.t$box[1,2]],
                                       y2.t[x2.t[, 1] <= res.t$box[2,1] &
                                            x2.t[, 1] >= res.t$box[1,1] &
                                            x2.t[, 2] <= res.t$box[2,2] &
                                            x2.t[, 2] >= res.t$box[1,2]], 1, 1))


        ## box.size <- (res.t$box[2,1] - res.t$box[1,1]) * (res.t$box[2,2] -
        ##                                                  res.t$box[1,2])
        ## box.size.correct <- (ifelse(res.t$box[2,1]<2,res.t$box[2,1],2) -
        ##                      res.t$box[1,1]) * (res.t$box[2,2] - res.t$box[1,2])
        ## if (box.size.correct < 0)
        ##     box.size.correct <- 0

        ## box.size.inccorect <- box.size - box.size.correct

        ## eo.1 <- x1.ord * 3/4 * box.size/100 + x1.ab * 3/4 * box.size.correct/20
        ## es.1 <- x1.ord * 1/4 * box.size/100 + x1.ab * 1/4 * box.size.correct/20
        ## eo.2 <- x2 * 3/4 * box.size/100
        ## es.2 <- x2 * 1/4 * box.size/100

        ##e.p <- lrt.pois(cbind(c(eo.1, es.1), c(eo.2, es.2)))

        exp.p[i] <- e.p
        true.p[i] <- t.p
        results[[i]] <- res.t
        #if (e.p > 1e-15)
        if (res.t$box[2,1] > 6)
            browser()
            #print(res.t$box)

        ##true.p[i] <- lrt.pois(count.fires(y1.d[x11.d < 2], y2.d[x21.d < 2], 1, 1))
    }
    list(exp.p, true.p, results)
}

plot.box = function(x1, y1, x2, y2, box, add=FALSE, ...) {
  if(!add) {
    plot(x1[,1], x1[,2], col="blue")
    points(x2[,1], x2[,2], col="black", pch=3, ...)
  }
  rect(box[1,1], box[1,2], box[2,1], box[2,2], border="red")
}

sim.res <- sim(400, 50, 400, 1000)

es.dens <- density(log(sim.res[[1]]))
tr.dens <- density(log(sim.res[[2]]))

plot(es.dens, xlim=c(xmi <- min(c(es.dens$x, tr.dens$x)),
              xma <- max(c(es.dens$x, tr.dens$x))),
              ylim=ymi <- c(min(c(es.dens$y, tr.dens$y)),
              yma <- max(c(es.dens$y, tr.dens$y))),
              main="P-value of resulting boxes", xlab="P (log scale)", xaxt="n")
lines(tr.dens, lty=2)
axis(1, at=xax <- seq(xmi, xma, by=5), label=signif(exp(xax), 2))
legend(xmi, yma, c("Estimated Box", "True Box"), lty=c(1,2))

hist(log(sim.res[[1]]), ,xlim=c(xmi, xma),  breaks=30, xaxt="n",
     main="P-values of estimated box", xlab="P (log scale)")
axis(1, at=xax <- seq(xmi, xma, by=5), label=signif(exp(xax), 2))

hist(log(sim.res[[2]]), ,xlim=c(xmi, xma),  breaks=30, xaxt="n",
     main="P-values of true box", xlab="P (log scale)")
axis(1, at=xax <- seq(xmi, xma, by=5), label=signif(exp(xax), 2))





sim.res2 <- sim(400, 200, 400, 50, alpha.paste=0.10)
median(sim.res2[[1]])
median(sim.res2[[2]])

################################################################
# Working Progress
############

#######
# generate some dummy data


dum <- function(x1.ord, x1.ab, x2, seed) {
    set.seed(seed+1)
    x11.d <- runif(x1.ord, 0, 10)
    x11.d <- c(x11.d, runif(x1.ab, 0, 2))
    set.seed(seed+2)
    x12.d <- runif(x1.ord+x1.ab, 0, 10)

    set.seed(seed+3)
    x21.d <- runif(x2, 0, 10)
    set.seed(seed+4)
    x22.d <- runif(x2, 0, 10)

    x1.d <- cbind(x11.d, x12.d)
    x2.d <- cbind(x21.d, x22.d)

    set.seed(seed+5)
    y1.d <- c(sample(c("other", "suspicious"), x1.ord, replace=TRUE, prob=c(3/4, 1/4)),
              sample(c("other", "suspicious"), x1.ab, replace=TRUE, prob=c(3/4, 1/4)))

    set.seed(seed+6)
    y2.d <- sample(c("other", "suspicious"), x2, TRUE, c(3/4, 1/4))

    init.peel = peeling(x1.d, x2.d, y1.d, y2.d, 0.05, 50)
    ## pasting(x1.d, x2.d, y1.d, y2.d,
    ##         init.peel[[length(init.peel)]]$box,
    ##         init.peel[[length(init.peel)]]$wts,
    ##         init.peel[[length(init.peel)]]$pvalue,
    ##         0.05)
    res.t <- prim(x1.d, x2.d, y1.d, y2.d, 0.05, 50, 5)
    res.t$x1 <- x1.d
    res.t$x2 <- x2.d
    res.t$y1 <- y1.d
    res.t$y2 <- y2.d
    res.t
}

my.lrt <- function(x1, x2, y1, y2) {
    w <- 2 * (dpois(x1, x1, TRUE) + dpois(y1, y1, TRUE) - dpois(x1, mean(c(x1, y1)), TRUE) - dpois(y1, mean(c(x1, y1)), TRUE)) +
        2 * (dpois(x2, x2, TRUE) + dpois(y2, y2, TRUE) - dpois(x2, mean(c(x2, y2)), TRUE) - dpois(y2, mean(c(x2, y2)), TRUE))
    p <- pchisq(w, 2, lower.tail=FALSE)
    attr(p, 'stat') <- w
    p
}


for (i in 1:10){
    tester <- dum(400, 200, 400, 30+i)
}
tester <- dum(400, 200, 400, 30+4)
#tester <- dum(200, 200, 400)
tester$box
tester$cv.p
sum(tester$wts)

hist(tester$x1)
hist(tester$x2)

p.t <- peeling(tester$x1, tester$x2, tester$y1, tester$y2, 0.05, 50)
p.t2 <- peeling(tester$x1, tester$x2, tester$y1, tester$y2, 0.15, 50)
p.t2[[length(p.t2)]]$box

unlist(sapply(p.t, function(m) m$pvalue[1]))
sapply(p.t, function(m) m$box)
sapply(p.t, function(m) m$pvalues)

## r = dum(200, 200, 400)

## sum(tester$wts > 0)


## dum2 <- function(x1.ord, x1.ab, x2, nx) {
##     x11.d <- runif(x1.ord, 0, 10)
##     x11.d <- c(x11.d, runif(x1.ab, 0, 2))
##     x12.d <- matrix(runif((nx-1)*(x1.ord+x1.ab), 0, 10), ncol=nx-1)
##     x2.d <- matrix(runif(nx*x2, 0, 10), ncol=nx)

##     x1.d <- cbind(x11.d, x12.d)

##     y1.d <- c(sample(c("other", "suspicious"), x1.ord, replace=TRUE, prob=c(3/4, 1/4)),
##               sample(c("other", "suspicious"), x1.ab, replace=TRUE, prob=c(2/4, 2/4)))

##     y2.d <- sample(c("other", "suspicious"), x2, TRUE, c(1/4, 3/4))


##     res.t <- prim(x1.d, x2.d, y1.d, y2.d, 0.05, 50, 5)
##     res.t
## }

## tester <- dum2(400, 200, 200, 5)
## tester$box
## tester$cv.

################################################
# Do the peeling manually
################################################

### pre-analysis
par(mfrow=c(2,2))
hist(tester$x1[, 1], main="Year 1, Variable 1")
hist(tester$x2[, 1], main="Year 2, Variable 1")
hist(tester$x1[, 2], main="Year 2, Variable 1")
hist(tester$x2[, 2], main="Year 2, Variable 2")

par(mfrow=c(2,1))
hist(c(tester$x1[, 1], tester$x2[, 1]), main="Variable 1")
hist(c(tester$x1[, 2], tester$x2[, 2]), main="Variable 2")

### first step
#
wts <- rep(1, length(tester$y1) + length(tester$y2))
tot.n <- length(tester$y1) + length(tester$y2)
n.rem <- floor(tot.n * 0.05)

# get the 2 variables over both years as single vectors
x.var1 <- c(tester$x1[, 1], tester$x2[, 1])
x.var2 <- c(tester$x1[, 2], tester$x2[, 2])

# sort them and return the index for this
x.v1.sort.ix <- sort(x.var1, index.return=TRUE)$ix
x.v2.sort.ix <- sort(x.var2, index.return=TRUE)$ix

# fires in each year
table(tester$y1)
table(tester$y2)

# p-value using simple lrt
my.lrt(table(tester$y1)[1], table(tester$y1)[2], table(tester$y2)[1], table(tester$y2)[2])

### lower boundary
# variable 1
wts.temp <- wts
wts.temp[x.v1.sort.ix[1:n.rem]] <- 0
p1.lw <- lrt.pois(count.fires(tester$y1, tester$y2, wts.temp[1:400], wts.temp[401:800]))
# analysis
table(tester$y1[wts.temp[1:400] > 0])
table(tester$y2[wts.temp[401:800] > 0])
# box after removal
my.lrt(table(tester$y1[wts.temp[1:400] > 0])[1], table(tester$y1[wts.temp[1:400] > 0])[2],
       table(tester$y2[wts.temp[401:800] > 0])[1], table(tester$y2[wts.temp[401:800] > 0])[2])
# removed box
my.lrt(table(tester$y1[wts.temp[1:400] == 0])[1], table(tester$y1[wts.temp[1:400] == 0])[2],
       table(tester$y2[wts.temp[401:800] == 0])[1], table(tester$y2[wts.temp[401:800] == 0])[2])


# variable 2
wts.temp <- wts
wts.temp[x.v2.sort.ix[1:n.rem]] <- 0
p2.lw <- lrt.pois(count.fires(tester$y1, tester$y2, wts.temp[1:400], wts.temp[401:800]))
# analysis
table(tester$y1[wts.temp[1:400] > 0])
table(tester$y2[wts.temp[401:800] > 0])
# box after removal
my.lrt(table(tester$y1[wts.temp[1:400] > 0])[1], table(tester$y1[wts.temp[1:400] > 0])[2],
       table(tester$y2[wts.temp[401:800] > 0])[1], table(tester$y2[wts.temp[401:800] > 0])[2])
# removed box
my.lrt(table(tester$y1[wts.temp[1:400] == 0])[1], table(tester$y1[wts.temp[1:400] == 0])[2],
       table(tester$y2[wts.temp[401:800] == 0])[1], table(tester$y2[wts.temp[401:800] == 0])[2])

### upper boundary
# variable 1
wts.temp <- wts
wts.temp[x.v1.sort.ix[(801-n.rem):800]] <- 0
p1.up <- lrt.pois(count.fires(tester$y1, tester$y2, wts.temp[1:400], wts.temp[401:800]))
# analysis
table(tester$y1[wts.temp[1:400] > 0])
table(tester$y2[wts.temp[401:800] > 0])
# box after removal
my.lrt(table(tester$y1[wts.temp[1:400] > 0])[1], table(tester$y1[wts.temp[1:400] > 0])[2],
       table(tester$y2[wts.temp[401:800] > 0])[1], table(tester$y2[wts.temp[401:800] > 0])[2])
# removed box
my.lrt(table(tester$y1[wts.temp[1:400] == 0])[1], table(tester$y1[wts.temp[1:400] == 0])[2],
       table(tester$y2[wts.temp[401:800] == 0])[1], table(tester$y2[wts.temp[401:800] == 0])[2])

# variable 2
wts.temp <- wts
wts.temp[x.v2.sort.ix[(801-n.rem):800]] <- 0
p2.up <- lrt.pois(count.fires(tester$y1, tester$y2, wts.temp[1:400], wts.temp[401:800]))
# analysis
table(tester$y1[wts.temp[1:400] > 0])
table(tester$y2[wts.temp[401:800] > 0])
# box after removal
my.lrt(table(tester$y1[wts.temp[1:400] > 0])[1], table(tester$y1[wts.temp[1:400] > 0])[2],
       table(tester$y2[wts.temp[401:800] > 0])[1], table(tester$y2[wts.temp[401:800] > 0])[2])
# removed box
my.lrt(table(tester$y1[wts.temp[1:400] == 0])[1], table(tester$y1[wts.temp[1:400] == 0])[2],
       table(tester$y2[wts.temp[401:800] == 0])[1], table(tester$y2[wts.temp[401:800] == 0])[2])


# p value matrix
pvalues <- matrix(c(p1.lw, p1.up, p2.lw, p2.up), nrow=2)
pvalues



### post-analysis
par(mfrow=c(2,2))
hist(subset(tester$x1[, 1] * wts.temp[1:400], tester$x1[, 1]*wts.temp[1:400]>0), main="Year 1, Variable 1", xlab="x")
hist(subset(tester$x1[, 2] * wts.temp[1:400], tester$x1[, 2]*wts.temp[1:400]>0), main="Year 1, Variable 2", xlab="x")
hist(subset(tester$x2[, 1] * wts.temp[401:800], tester$x2[, 1]*wts.temp[401:800]>0), main="Year 2, Variable 1", xlab="x")
hist(subset(tester$x2[, 2] * wts.temp[401:800], tester$x2[, 2]*wts.temp[401:800]>0), main="Year 2, Variable 2", xlab="x")

par(mfrow=c(2,1))
hist((x.var1 * wts.temp)[x.var1*wts.temp>0], main="Variable 1")
hist((x.var2 * wts.temp)[x.var2*wts.temp>0], main="Variable 2")

table(tester$y1[wts.temp[1:400]>0])
table(tester$y2[wts.temp[401:800]>0])


### second step
#

wts[x.v1.sort.ix[1:n.rem]] <- 0
n.rem <- floor(sum(wts) * 0.05)

# lower boundary
# variable 1
wts.temp <- wts
wts.temp[x.v1.sort.ix[41:(41+n.rem)]] <- 0
p1.lw <- lrt.pois(count.fires(tester$y1, tester$y2, wts.temp[1:400], wts.temp[401:800]))
# variable 2
wts.temp <- wts
wts.temp[x.v2.sort.ix[41:(41+n.rem)]] <- 0
p2.lw <- lrt.pois(count.fires(tester$y1, tester$y2, wts.temp[1:400], wts.temp[401:800]))

# upper boundary
# variable 1
wts.temp <- wts
wts.temp[x.v1.sort.ix[(801-n.rem):800]] <- 0
p1.up <- lrt.pois(count.fires(tester$y1, tester$y2, wts.temp[1:400], wts.temp[401:800]))
# variable 2
wts.temp <- wts
wts.temp[x.v2.sort.ix[(801-n.rem):800]] <- 0
p2.up <- lrt.pois(count.fires(tester$y1, tester$y2, wts.temp[1:400], wts.temp[401:800]))

# p value matrix
pvalues









####################################
# Extra Functions

