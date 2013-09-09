####### Master Thesis #########
####
# ToDo:
###
# Change pasting. Use average w values per observation. Add on remaining obs if it is <n.add. Use greedy procedure. Still end up with same sized boxes. DONE
# Test performance for longer peeling/pasting sequences
# Use w and different pvalues for simulation
##################################3
# Searching for significant differential rules #

# --- Read in data ---


lrt.pois <- function(O) {
# Calculate the p-value for the difference
# log-likelihood ratio test using Possion distributions for count data
#
# input:
# O: matrix containing the number of other and suspicious fires
#    in the 2 datasets, 1 col for each dataset
#
# output:
#  The p-value and the value of the likelihood ratio test statistic

  if(length(dim(O)) == 2) dim(O) <- c(1,dim(O))
  O2 <- apply(O, 1:2, mean)
  w <- 2 * ( apply(dpois.cont(O, O, log=TRUE), 1, sum) -
    apply(dpois.cont(O, O2, log=TRUE), 1, sum) )
  p <- pchisq(w, dim(O)[2], lower=FALSE)
  attr(p, "stat") <- w
  p
}

## dpois.cont <- function(x, lambda, log = FALSE) {
##     if (class(x) == "array") {
##         if (log)
##             apply(x, c(1,2), function(x){x * log(lambda) - lambda - lgamma(x + 1)})
##         else
##             lambda^x * exp(-lambda) / gamma(x + 1)
##     }else{
##         if (log){
##             print("pups")
##             x * log(lambda) - lambda - lgamma(x + 1)
##         }else
##             lambda^x * exp(-lambda) / gamma(x + 1)
##     }
## }

dpois.cont <- function(x, lambda, log = FALSE) {
    x2 <- as.vector(x)
    l2 <- as.vector(lambda)
    cons <- rep(1, length(l2))
    cons <- sapply(l2, function(x) {ifelse(x <= 10, 1/my.int(x)$value, 1)})
    if (log)
        res <- log(cons) + x2 * log(l2) - l2 - lgamma(x2 + 1)
    else
        res <- cons * l2^x2 * exp(-l2) / gamma(x2 + 1)
    res <- sapply(res, function(x) {ifelse(is.nan(x), 0, x)})
    if (length(res) > 2)
        res <- matrix(res, nrow=2)
    dim(res) <- c(1, dim(res))
    res
}

norm.const <- function(x, lambda = 130) {
    exp(-lambda + x * log(lambda) - lgamma(x + 1))
}

my.int <- function(lambda) {
    integrate(norm.const, 0, Inf, lambda = lambda)
}

count.fires <- function(data1, data2, w1, w2){
# Count the number of 'suspicious' and 'other' fires in 2 datasets
#
# input:
#   two data sets with 'label' column for fire data
#   'label' is a factor with levels 'other' and 'suspicious'
#
# output:
#   matrix with counts for 'other' and 'suspicious' fires
#   Output matrix can be used for lrt.pois()

    o1 <- sum((data1 == 'other') * w1)
    o2 <- sum((data2 == 'other') * w2)

    s1 <- sum((data1 == 'suspicious') * w1)
    s2 <- sum((data2 == 'suspicious') * w2)

    matrix(c(o1,s1,o2,s2), nr=2)
}


####
# PRIM
####
##########################################################
##########################################################


peel.step <- function(x.ordered.ix, x.ordered, wts, y1, y2, box, alpha, beta.0, d, forced=FALSE){

  ######################################################################
  # Peel of a small number of the current wts. The wts chosen are the ones
  # that result in the smallest pvalue
  #
  # Arguments:
  #   x.ordered: data.frame of initial columns with cont. variables,
  #                           ordered with decr.=TRUE
  #   x.ordered.ix: data.frame of indicies that link ordered x to the
  #                           original x
  #   y1, y2: vectors with the dependent variables of current box
  #   box: the current box boundaries
  #   alpha: percentage of obs of current box that are added in the step
  #   categorical: index of variables that are treated as categorical
  #   d: number of variables(columns)
  ######################################################################
  ####
  # Build nec. matrices and calc. new boundaries
  ####


  # nr observations to remove
  n.rem = alpha * sum(wts)

  # will store new boundaries of box
  new.boundaries <- matrix(NA, nrow=2, ncol=d)
  # will store p-values for all the new boundares
  new.pvalue <- matrix(1, nrow=2, ncol=d)
  # store the new box boundaries (update later)
  new.box <- box


  ####
  # Do one step of the peeling procedure
  ####

  for (i in 1:d) {

      ## upper boundary
      # weights index for weights that need to be altered
      #ix <- wts[x.ordered.ix[, i]] > 0 & cumsum(wts[x.ordered.ix[, i]]) <= n.rem
      ix <- cumsum(wts[x.ordered.ix[, i]]) <= n.rem & cumsum(wts[x.ordered.ix[, i]]) > 0
      # index of first weight that is positive
      min.ix <- which(ix)[1]
      # index of last weight that needs to be altered
      max.ix <- which(ix)[length(which(ix))]

      # if no ix meets the condition (i.e. first weight is larger than n.rem)
      # set max.ix to the index before the first weight
      if (all(ix == FALSE)) {
          to.remove.ix <- ix
          max.ix <- min(which(cumsum(wts[x.ordered.ix[, i]]) > 0)) - 1
          wt.to.remove <- n.rem
      } else {
          # if the next lowest variable is unequal to the one for max.ix,
          # remove all weights
          if (x.ordered[max.ix, i] != x.ordered[max.ix + 1, i]) {
              # ix of weights that will be set to zero
              to.remove.ix <- ix
              # total weight removed
              wt.removed <- sum(wts[x.ordered.ix[to.remove.ix, i]])
              # weight left to remove
              wt.to.remove <- n.rem - wt.removed
          }else{
              # check if any weights need to be set to 0
              if (x.ordered[min.ix, i] != x.ordered[max.ix, i]) {
                  # ix of weights that will be set to zero
                  to.remove.ix <- ix & (x.ordered[, i] != x.ordered[max.ix, i])
                  # total weight removed
                  wt.removed <- sum(wts[x.ordered.ix[to.remove.ix, i]])
                  # weight left to remove
                  wt.to.remove <- n.rem - wt.removed
              }else {
                  wt.to.remove <- n.rem
                  to.remove.ix <- rep(FALSE, length(wts))
              }
          }
      }

    # get temp weights and set some to 0
    wts.temp <- wts
    wts.temp[x.ordered.ix[to.remove.ix, i]] <- 0

    # number of duplicates at the boundary
    # nr.dupl <- sum(x.ordered[, i] == x.ordered[max.ix + 1, i])
    if (wt.to.remove > 0) {
        # total weight of all the duplicates at the boundary
        wts.bound.ttl <- sum(wts.temp[x.ordered.ix[x.ordered[, i] == x.ordered[max.ix + 1, i], i]])
        wts.temp[x.ordered.ix[x.ordered[, i] == x.ordered[max.ix + 1, i], i]] <- wts.temp[x.ordered.ix[x.ordered[, i] == x.ordered[max.ix + 1, i], i]] * (wts.bound.ttl - wt.to.remove)/wts.bound.ttl
    }

    count.temp <- count.fires(y1, y2, wts.temp[1:length(y1)], wts.temp[(length(y1) + 1):length(wts.temp)])
    p.temp <- lrt.pois(count.temp)

    # check if this is the lowest p-value so far
    if (p.temp < min(new.pvalue)) {
        boundary <- 2
        new.wts <- wts.temp
        var.select <- i
        new.p <- p.temp
    }

    new.pvalue[2, i] <- p.temp
    new.boundaries[2, i] <- x.ordered[max.ix + 1, i]



    ## lower boundary
    # weights index for weights that need to be altered
    ix <- rev(cumsum(rev(wts[x.ordered.ix[, i]]))) <= n.rem & rev(cumsum(rev(wts[x.ordered.ix[, i]]))) > 0
    # index of first weight that is positive
    min.ix <- which(ix)[1]
    # index of last weight that needs to be altered
    max.ix <- which(ix)[length(which(ix))]

    if (all(ix == FALSE)){
        to.remove.ix <- ix
        min.ix <- max(which(rev(cumsum(rev(wts[x.ordered.ix[, i]]))) > 0)) + 1
        wt.to.remove <- n.rem
    } else {
        # if the next lowest variable is unequal to the one for min.ix,
        # remove all weights
        if (x.ordered[min.ix, i] != x.ordered[min.ix - 1, i]) {
            # ix of weights that will be set to zero
            to.remove.ix <- ix
            # total weight removed
            wt.removed <- sum(wts[x.ordered.ix[to.remove.ix, i]])
            # weight left to remove
            wt.to.remove <- n.rem - wt.removed
        }else{
            # check if any weights need to be set to 0
            if (x.ordered[min.ix, i] != x.ordered[max.ix, i]) {
                # ix of weights that will be set to zero
                to.remove.ix <- ix & (x.ordered[, i] != x.ordered[min.ix, i])
                # total weight removed
                wt.removed <- sum(wts[x.ordered.ix[to.remove.ix, i]])
                # weight left to remove
                wt.to.remove <- n.rem - wt.removed
            }else {
                wt.to.remove <- n.rem
                to.remove.ix <- rep(FALSE, length(wts))
            }
        }
    }

    # get temp weights and set some to 0
    wts.temp <- wts
    wts.temp[x.ordered.ix[to.remove.ix, i]] <- 0

    # number of duplicates at the boundary
    #nr.dupl <- sum(x.ordered[, i] == x.ordered[min.ix - 1, i])
    if (wt.to.remove > 0) {
        # total weight of all the duplicates at the boundary
        wts.bound.ttl <- sum(wts.temp[x.ordered.ix[x.ordered[, i] == x.ordered[min.ix - 1, i], i]])
        wts.temp[x.ordered.ix[x.ordered[, i] == x.ordered[min.ix - 1, i], i]] <- wts.temp[x.ordered.ix[x.ordered[, i] == x.ordered[min.ix - 1, i], i]] * (wts.bound.ttl - wt.to.remove)/wts.bound.ttl
    }

    ## # remove some weight of the boundary weights
    ##   if (any(wts.temp[x.ordered.ix[x.ordered[, i] == x.ordered[min.ix - 1, i], i]] < wt.to.remove/nr.dupl)) {
    ##       wts.temp2 <- wts.temp[x.ordered.ix[x.ordered[, i] == x.ordered[min.ix - 1, i], i]]
    ##       to.remove.ix2 <- which(wts.temp2 < wt.to.remove/nr.dupl)
    ##       wt.to.remove <- wt.to.remove - sum(wts.temp2[to.remove.ix2])
    ##       wts.temp2[to.remove.ix2] <- 0
    ##       wts.temp[x.ordered.ix[x.ordered[, i] == x.ordered[min.ix - 1, i], i]] <- wts.temp2
    ##   }
    ## wts.temp[x.ordered.ix[x.ordered[, i] == x.ordered[min.ix - 1, i] & wts.temp[x.ordered.ix[ ,i]] > 0, i]] <- wts.temp[x.ordered.ix[x.ordered[, i] == x.ordered[min.ix - 1, i] & wts.temp[x.ordered.ix[ ,i]] > 0, i]] - wt.to.remove/nr.dupl

    count.temp <- count.fires(y1, y2, wts.temp[1:length(y1)], wts.temp[(length(y1) + 1):length(wts.temp)])
    p.temp <- lrt.pois(count.temp)

    # check if this is the lowest p-value so far
    if (p.temp < min(new.pvalue)) {
        boundary <- 1
        new.wts <- wts.temp
        var.select <- i
        new.p <- p.temp
    }

      new.pvalue[1, i] <- p.temp
      new.boundaries[1, i] <- x.ordered[min.ix - 1, i]
  }

  new.box[boundary, var.select] <- new.boundaries[boundary, var.select]

  if (sum(new.wts) >= beta.0 | forced) {
      return(list(wts=new.wts, box=new.box, pvalue=new.p, pvalues=new.pvalue,
                  updated=c(boundary, var.select)))
  }
}

paste.step <- function(x.ordered, x.ordered.ix, x.init, wts, y1, y2, box, alpha,
                       d, forced=FALSE, reduced.wt=NULL){
  ######################################################################
  # Paste a small box onto the current supplied box. Box chosen is the one
  # with smallest p-value
  #
  # Arguments:
  #   x1, x2: datasets containing independent variables of the current box
  #   x1.init, x2.init: datatasets containing original independent variables
  #                     (i.e. of box that contains all obs)
  #   x.ordered: list of initial columns with cont. variables,
  #                           ordered with decr.=TRUE
  #   y1, y2: vectors with the dependent variables of current box
  #   y1.init, y2.init: vectors with original dependent variables
  #                     (i.e. of box that contains all obs)
  #   box: the current box boundaries
  #   alpha: percentage of obs of current box that are added in the step
  #   categorical: index of variables that are treated as categorical
  #   d: number of variables(columns)
  #   forced: if true, the step is always done regardless of p-value
  #   reduced.wt: if weight is not determined by box size but set fixed
  ######################################################################

  ####
  # Build nec. matrices and calc. new boundaries
  ####

  # get the number of variables in the box
  n <- sum(wts)

  # store new upper and lower boundaries of all new boxes
  new.boundaries <- matrix(0, nrow=2, ncol=d)
  # will store p-values for all the new boundares
  new.pvalue <- matrix(1, nrow=2, ncol=d)
  # will store the stat value per observation
  new.wvalue <- matrix(0, nrow=2, ncol=d)

  # get the amount of weight added through the new box
  if (is.null(reduced.wt))
      wt.to.add <- n * alpha
  else
      wt.to.add <- reduced.wt

  var.select <- FALSE

  if (!forced)
      curr.p <- lrt.pois(count.fires(y1, y2, wts[1:length(y1)], wts[(length(y1) + 1):length(wts)]))
  else
      curr.p <- 1

  # each column holds the index of variables that are eligible for pasting
  # for the variable corresponding to that column
  lower.add <- matrix(TRUE, nrow = dim(x.init)[1], ncol = dim(x.init)[2])
  upper.add <- matrix(TRUE, nrow = dim(x.init)[1], ncol = dim(x.init)[2])

  for (i in 1:dim(x.init)[2]) {
      for (j in 1:dim(x.init)[2]) {
          if (j == i) {
              lower.add[, i] <- lower.add[, i] & x.init[, j] <= box[1, j]
              upper.add[, i] <- upper.add[, i] & x.init[, j] >= box[2, j]
              next
          }
          lower.add[, i] <- lower.add[, i] & x.init[, j] >= box[1, j] & x.init[, j] <= box[2, j]
          upper.add[, i] <- upper.add[, i] & x.init[, j] >= box[1, j] & x.init[, j] <= box[2, j]
      }
  }

  # Lower Boundary
  for (i in 1:d) {
      # reduced will be set to TRUE, if less than the required weight is left
      reduced.temp <- FALSE
      wts.temp <- wts
      # total weight outside the box eligible for pasting
      wt.left <- sum((1 - wts)[lower.add[, i]])

      # if no weight is left, skip this variable
      if (wt.left == 0)
          next

      # if not enough weight is left to add, add the remaining weight
      if (wt.to.add > wt.left) {
          reduced.temp <- TRUE
          wts.eli <- lower.add[x.ordered.ix[, i], i] * (1 - wts)[x.ordered.ix[, i]]
          to.one <- which(wts.eli > 0)
          wts.temp[x.ordered.ix[to.one, i]] <- 1
      } else {
          # weights of eligible variables ordered by according to var i
          wts.eli <- lower.add[x.ordered.ix[, i], i] * (1 - wts)[x.ordered.ix[, i]]
          # ix of wts that shld be added by paste step
          ix.in <- which(cumsum(wts.eli) < wt.to.add & wts.eli > 0)
          # ix of wts that are not added (except boundary)
          ix.out <- which(cumsum(wts.eli) >= wt.to.add & wts.eli > 0)

          # if there are no obs that are completely added to the box, just do boundary
          # handling
          if (length(ix.in) == 0) {
              to.one <- FALSE
              wt.left.to.add <- wt.to.add
          } else {
              # if the last value of variable to modify weight for is different from the one before,
              # set every weight but the last one to 1
              if (x.ordered[tail(ix.in, 1), i] != x.ordered[head(ix.out, 1), i]) {
                  to.one <- ix.in
                  wt.left.to.add <- wt.to.add - cumsum(wts.eli)[tail(ix.in, 1)]
              } else {
                  # set all weights to one where the observation is unequal to the on at the boundary
                  to.one <- ix.in[x.ordered[ix.in, i] != x.ordered[tail(ix.in, 1), i]]
                  if (length(to.one) == 0)
                      to.one <- FALSE
                  wt.left.to.add <- wt.to.add - sum(wts.eli[to.one])
              }
          }
#print(paste("lower:", i))
          # modify the weight at boundary
          wts.temp[x.ordered.ix[to.one, i]] <- 1
          bord.dupl <- c(ix.in, ix.out)[x.ordered[c(ix.in, ix.out), i] == x.ordered[head(ix.out, 1), i]]
          # check if any wts will reach 1 with the new addition
          while (any(wts.eli[bord.dupl] < wt.left.to.add/length(bord.dupl))) {
              #wts.temp[x.ordered.ix[wts.eli < wt.left.to.add/length(bord.dupl), i]] <- 1
              # ix of the weights that will reach 1 with the new addition
              ix.temp <- wts.eli[bord.dupl] < wt.left.to.add/length(bord.dupl)
              # set those weights to 1
              wts.temp[x.ordered.ix[bord.dupl[ix.temp], i]] <- 1
              # reduce the amount of weight left to add
              wt.left.to.add <- wt.left.to.add - sum(wts.eli[bord.dupl[ix.temp]])
              # deduct the modified observations from the observations at the boundary
              # that are yet to be modified
              bord.dupl <- bord.dupl[!ix.temp]
          }

          wts.temp[x.ordered.ix[bord.dupl, i]] <- wts.temp[x.ordered.ix[bord.dupl, i]] +
              wt.left.to.add/length(bord.dupl)
      }


      count.temp <- count.fires(y1, y2, wts.temp[1:length(y1)], wts.temp[(length(y1) + 1):length(wts.temp)])
      p.temp <- lrt.pois(count.temp)

      # calculate w value per observation
      if (reduced.temp)
          w.temp <- attr(p.temp, 'stat') / wt.left
      else
          w.temp <- attr(p.temp, 'stat') / wt.to.add

      if (w.temp > max(new.wvalue) & p.temp < curr.p) {
          boundary <- 1
          reduced <- reduced.temp
          new.wts <- wts.temp
          var.select <- i
          new.p <- p.temp
          if (reduced)
              new.boundary <- min(x.init[, i])
          else
              new.boundary <- x.ordered[head(ix.out, 1), i]
      }

      new.pvalue[1, i] <- p.temp
      new.wvalue[1, i] <- w.temp
  }

  # Upper boundary
  for (i in 1:d) {
      # reduced will be set to TRUE, if less than the required weight is left
      reduced.temp <- FALSE
      wts.temp <- wts
      # total weight outside the box eligible for pasting
      wt.left <- sum((1 - wts)[upper.add[, i]])

      # if no weight is left, skip this variable
      if (wt.left == 0)
          next

      # if not enough weight is left to add, add remaining weight
      if (wt.to.add > wt.left) {
          reduced.temp <- TRUE
          wts.eli <- upper.add[x.ordered.ix[, i], i] * (1 - wts)[x.ordered.ix[, i]]
          to.one <- which(wts.eli > 0)
          wts.temp[x.ordered.ix[to.one, i]] <- 1
      } else {
          # weights of eligible variables ordered by according to var i
          wts.eli <- upper.add[x.ordered.ix[, i], i] * (1 - wts)[x.ordered.ix[, i]]
          # ix of wts that shld be added by paste step
          ix.in <- which(rev(cumsum(rev(wts.eli))) < wt.to.add & wts.eli > 0)
          # ix of wts that are not added (except boundary)
          ix.out <- which(rev(cumsum(rev(wts.eli))) >= wt.to.add & wts.eli > 0)

          # if there are no obs that are completely added to the box, just do boundary
          # handling
          if (length(ix.in) == 0) {
              to.one <- FALSE
              wt.left.to.add <- wt.to.add
          } else {
              # if the last value of variable to modify weight for is different from the one before,
              # set every weight but the last one to 1
              if (x.ordered[head(ix.in, 1), i] != x.ordered[tail(ix.out, 1), i]) {
                  to.one <- ix.in
                  #wt.left.to.add <- wt.to.add - rev(cumsum(rev(wts.eli)))[head(ix.in, 1)]
                  wt.left.to.add <- wt.to.add - sum(wts.eli[ix.in])
              } else {
                  # set all weights to one where the observation is unequal to the on at the boundary
                  to.one <- ix.in[x.ordered[ix.in, i] != x.ordered[head(ix.in, 1), i]]
                  if (length(to.one) == 0)
                      to.one <- FALSE
                  wt.left.to.add <- wt.to.add - sum(wts.eli[to.one])
              }
          }
#print(paste("upper:", i))
          # modify the weight at boundary
          wts.temp[x.ordered.ix[to.one, i]] <- 1
          bord.dupl <- c(ix.in, ix.out)[x.ordered[c(ix.in, ix.out), i] == x.ordered[tail(ix.out, 1), i]]
          # check if any wts will reach 1 with the new addition
          while (any(wts.eli[bord.dupl] < wt.left.to.add/length(bord.dupl))) {
              #wts.temp[x.ordered.ix[wts.eli < wt.left.to.add/length(bord.dupl), i]] <- 1
              # ix of the weights that will reach 1 with the new addition
              ix.temp <- wts.eli[bord.dupl] < wt.left.to.add/length(bord.dupl)
              # set those weights to 1
              wts.temp[x.ordered.ix[bord.dupl[ix.temp], i]] <- 1
              # reduce the amount of weight left to add
              wt.left.to.add <- wt.left.to.add - sum(wts.eli[bord.dupl[ix.temp]])
              # deduct the modified observations from the observations at the boundary
              # that are yet to be modified
              bord.dupl <- bord.dupl[!ix.temp]
          }

          wts.temp[x.ordered.ix[bord.dupl, i]] <- wts.temp[x.ordered.ix[bord.dupl, i]] +
              wt.left.to.add/length(bord.dupl)
      }

      count.temp <- count.fires(y1, y2, wts.temp[1:length(y1)], wts.temp[(length(y1) + 1):length(wts.temp)])
      p.temp <- lrt.pois(count.temp)

      # calculate w value per observation
      if (reduced.temp)
          w.temp <- attr(p.temp, 'stat') / wt.left
      else
          w.temp <- attr(p.temp, 'stat') / wt.to.add

      if (w.temp > max(new.wvalue) & p.temp < curr.p) {
          boundary <- 2
          reduced <- reduced.temp
          new.wts <- wts.temp
          var.select <- i
          new.p <- p.temp
          if (reduced)
              new.boundary <- max(x.init[, i])
          else
              new.boundary <- x.ordered[tail(ix.out, 1), i]
      }

      new.pvalue[2, i] <- p.temp
      new.wvalue[2, i] <- w.temp
  }

  if (var.select) {

      new.box <- box
      new.box[boundary, var.select] <- new.boundary
      result <- list(wts = new.wts, box = new.box, pvalue = new.p, updated = c(boundary, var.select))
      if (reduced & sum(new.wts) < length(x.init[,1]))
          result <- paste.step(x.ordered, x.ordered.ix, x.init, new.wts, y1, y2, new.box, alpha,
                               d, forced=TRUE, reduced.wt=wt.to.add - (sum(new.wts) - sum(wts)))
      result
  } #else {
    #  browser()
    #  print(new.pvalue)
    #  print(sum(wts))
  #}
}

peeling <- function(x1, x2, y1, y2, alpha, beta.0, maxiter=1000, forced.steps=FALSE){

  # get number of variables
  if (is.null(dim(x1)))
    d <- 1
  else
    d <- dim(x1)[2]

  x.ordered <- apply(rbind(x1, x2), 2, function(x) {sort(x, decreasing=TRUE)})
  x.ordered.ix <- apply(rbind(x1, x2), 2, function(x) {sort(x, decreasing=TRUE, index.return=TRUE)$ix})

  box.org <- apply(rbind(x1, x2), 2, range)
  #updated <- matrix(c(0, 0), nrow=1)
  # initial p-value
  p.org <- lrt.pois(count.fires(y1, y2, rep(1, length(y1)), rep(1, length(y2))))

  # initial box
  wts <- rep(1, length(y1) + length(y2))
  result.temp <- list(wts=wts, box=box.org, updated=c(0,0), pvalue=p.org)
  box.new <- list()

  if (identical(forced.steps, FALSE))
      forced <- FALSE
  else
      forced <- TRUE

  min.p <- p.org
  i = 1
  while (!is.null(result.temp)){
    #print(i)
    if (forced & i > forced.steps )
        break

    box.new[[i]] <- result.temp
    box.new[[i]]$step <- i

    ## if (result.temp$pvalue > 100000*min.p & i >=20 & !forced)
    ##     break
    ## if (result.temp$pvalue < min.p)
    ##     min.p <- result.temp$pvalue

    result.temp <- peel.step(x.ordered.ix, x.ordered, result.temp$wts,
                             y1, y2, result.temp$box, alpha, beta.0, d, forced)
    #if(i == 1) plot.box(x1, y1, x2, y2, result.temp$box)
    #else plot.box(x1, y1, x2, y2, result.temp$box, add=TRUE)
    if (i == maxiter){
      warning("maxiter reached")
      break
    }

    i = i+1

  }

  box.new
}

pasting <- function(x1, x2, y1, y2, box, wts, p.init, alpha, forced.steps=FALSE) {

  # get number of variables
  if (is.null(dim(x1)))
    d <- 1
  else
    d <- dim(x1)[2]

  x.ordered <- apply(rbind(x1, x2), 2, function(x) {sort(x, decreasing=TRUE)})
  x.ordered.ix <- apply(rbind(x1, x2), 2, function(x) {sort(x, decreasing=TRUE, index.return=TRUE)$ix})

  result.temp <- list(wts=wts, box=box, pvalue=p.init)
  box.new <- list()

  if (identical(forced.steps, FALSE))
      forced <- FALSE
  else
      forced <- TRUE

  i = 1
  while (!is.null(result.temp)) {
      #print(i)
      if (forced & i > forced.steps )
          break
      if (any(result.temp$wts > 1) | any(result.temp$wts < 0))
          print("wt error")

      box.new[[i]] <- result.temp
      box.new[[i]]$step = i
      result.temp <- paste.step(x.ordered, x.ordered.ix, rbind(x1, x2),
                                result.temp$wts, y1, y2, result.temp$box,
                                alpha, d, forced)

      #if (is.null(result.temp))
      #    print("was null")
      i = i+1
  }


  box.new
}

cv <- function(x1, x2, y1, y2, peel, alpha.paste, beta.0, K = 10, maxiter = 1000, paste.steps) {

    #set.seed(42)
    # build the data sets for the cv
    n <- length(y1) + length(y2)
    rs <- runif(n)
    id <- seq(n)[order(rs)]
    k <- as.integer(n*seq(1, K-1)/K)
    k <- matrix(c(0,rep(k,each=2),n),ncol=2,byrow=TRUE)
    k[,1] <- k[,1]+1
    l <- lapply(seq.int(K),function(x,k,d)
                list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
                     test=d[seq(k[x,1],k[x,2])]),k=k,d=id)

    x.init <- rbind(x1, x2)
    y.init <- factor(c(as.character(y1), as.character(y2)))

    size <- list()
    w <- matrix(nrow=paste.steps, ncol=0)

    # get the weights for the original box
    wts <- rep(0, n)
    wts[peel$wts>0] <- 1

    for (i in 1:K) {
        paste <- pasting(x.init[l[[i]]$train[l[[i]]$train <= length(y1)], ],
                        x.init[l[[i]]$train[l[[i]]$train > length(y1)], ],
                        y.init[l[[i]]$train[l[[i]]$train <= length(y1)]],
                        y.init[l[[i]]$train[l[[i]]$train > length(y1)]],
                        peel$box, wts[l[[i]]$train], peel$pvalue, alpha.paste,
                        paste.steps)

        # matrices with variables as rows and boundaries as columns
        lw.bounds <- sapply(paste, function(m) m$box[1, ])
        up.bounds <- sapply(paste, function(m) m$box[2, ])

        # list of matrices with the lower and upper boundaries repeated according to the size of the test set
        lw.mats <- lapply(1:length(paste), function(j) matrix(rep(lw.bounds[, j], rep(length(l[[i]]$test), length(lw.bounds[, j]))), ncol = length(lw.bounds[, j])))
        up.mats <- lapply(1:length(paste), function(j) matrix(rep(up.bounds[, j], rep(length(l[[i]]$test), length(up.bounds[, j]))), ncol = length(up.bounds[, j])))

        # list of logical matrices indicating which observations are within the box
        lw.mats.in <- lapply(lw.mats, function(m) x.init[l[[i]]$test, ] >= m)
        up.mats.in <- lapply(up.mats, function(m) x.init[l[[i]]$test, ] <= m)

        # list of logical vectors indicating which rows are in the box
        lw.vec.in <- lapply(lw.mats.in, function(m) apply(m, 1, all))
        up.vec.in <- lapply(up.mats.in, function(m) apply(m, 1, all))

        # list of logical vectors combining both upper and lower boundary restriction
        in.box <- lapply(mapply(cbind, lw.vec.in, up.vec.in, SIMPLIFY=FALSE), function(m) apply(m, 1, all))

        # list of w values for each box from the pasting procedure
        w <- cbind(w, unlist(sapply(in.box, function(m) attr(lrt.pois(count.fires(y.init[l[[i]]$test[m & l[[i]]$test <= length(y1)]],
                                                        y.init[l[[i]]$test[m & l[[i]]$test > length(y1)]],
                                                        1,1)), 'stat'))))

        #w <- cbind(w, unlist(sapply(paste, function(m) attr(m$pvalue, 'stat'))))
        size[[i]] <- unlist(sapply(paste, function(m) sum(m$wts)))
    }
    w
}

prim <- function(x1, x2, y1, y2, alpha.peel, alpha.paste, beta, cv.K) {
    init.peel <- peeling(x1, x2, y1, y2, alpha.peel, beta)

    # get all the pvalue from the peeling sequence
    pvalues <- sapply(init.peel, function(m) m$pvalue[1])
    ## paste.steps <- which(pvalues < min(pvalues)*100000)[1]
    ## if (is.na(paste.steps))
    ##     paste.steps <- 0

    ## # cant paste back to a full set
    ## if (paste.steps < 2)
    ##     paste.steps <- 2

    ## paste.steps <- length(pvalues) - paste.steps

    paste.steps <- Inf

    init.paste <- pasting(x1, x2, y1, y2,
                          init.peel[[length(init.peel)]]$box,
                          init.peel[[length(init.peel)]]$wts,
                          init.peel[[length(init.peel)]]$pvalue,
                          alpha.paste, forced.steps=paste.steps)

    cv.w <- cv(x1, x2, y1, y2, init.peel[[length(init.peel)]], alpha.paste,
               beta, cv.K, paste.steps = length(init.paste))

    # get the average w stat for each box size
    w.average <- apply(cv.w, 1, median)
    best.size <- which(w.average == max(w.average))[1]

    tryCatch(init.paste[[best.size]], error = function(e) browser())
    box.chosen <- init.paste[[best.size]]
    p <- pchisq(max(w.average), 2, lower=FALSE)
    attr(p, "stat") <- max(w.average)
    box.chosen$cv.p <- p
    box.chosen
}


######################################################################################
# Testing
######################################################################################

arson = read.csv("arson.csv")
i = with(arson, day <= 731)
arson1 = arson[i,]
arson2 = arson[!i,]
arson2[,"day"] = arson2[,"day"] - 731

# declare continuous and categorical variables
contin.variables <- c(1,2,8,9)
categ.variables <- c(3,4,5,6,7,10)

#arson1[is.na(arson1[, 6]), 6] <- 10
#rson2[is.na(arson2[, 6]), 6] <- 10
arson1[is.na(arson1)] <- 0
arson2[is.na(arson2)] <- 0

arson1[, c(3:7, 10)] <- lapply(arson1[, c(3:7, 10)], factor)
arson2[, c(3:7, 10)] <- lapply(arson2[, c(3:7, 10)], factor)

x1 <- arson1[, -11]
x2 <- arson2[, -11]
y1 <- arson1[, 11]
y2 <- arson2[, 11]

d <- dim(x1)[2]

# categorical
categorical <- NULL
if (is.null(categorical)){
    # get the columns that are factors
    categ <- which(sapply(x1, is.factor))
    categorical <- rep(FALSE, d)
    categorical[categ] <- TRUE
}

# turn factors in numeric variables
x1[, categorical] <- sapply(x1[, categorical], function(x) {as.numeric(levels(x))[x]})
x2[, categorical] <- sapply(x2[, categorical], function(x) {as.numeric(levels(x))[x]})
