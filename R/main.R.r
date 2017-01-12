#' Empirical Bayes parameter ranking for parallel estimation scenarios
#'
#' Empirical Bayes ranking applicable to parallel-estimation settings where the estimated parameters are asymptotically unbiased and normal, with known standard errors.  A mixture normal prior for each parameter is estimated using Empirical Bayes methods, subsequentially ranks for each parameter are simulated from the resulting joint posterior over all parameters (The marginal posterior densities for each parameter are assumed independent). Finally, experiments are ordered by expected posterior  rank, although computations minimizing other plausible rank-loss functions are also given.  
#' @param betahat estimated effect sizes for each experiment
#' @param sebeta standard error of estimated effect sizes
#' @param Jmin minimum number of non-null clusters fit
#' @param Jmax maximum number of non-null clusters fit
#' @param maxiter maximum number of iterations for EM algorithm
#' @param tol EM algorithm is considered to have converged if the sum of the squared Euclidean distances between the parameter estimates on 2 iterations is less than tol
#' @param nsim number of simulations from posterior distribution
#' @param cutoff controls which experiments are included for posterior rank simulation. If a numeric between 0 and 1, it specifies the minimum posterior probability for inclusion in posterior rank simulations. If equal to 'f' then experiements in posterior rank simulation had p-values that were significant according to a Benjamini Hochberg correction at BH_FDR, if equal to 'b' posterior simulations correspond to experiments with Bonferoni significant p-values at level alpha.
#' @param maxpar maximum number of experiments to simulate
#' @param multiplestart if TRUE, multiple start points are used for the EM-algorithm based fitting of the mixture normals (for a given number of clusters) 
#' @param sigmabig the standard deviation for the 1st non-null cluster component
#' @param fixedcluster2 TRUE if the standard deviation for the 1st non-null cluster of the marginal distribution is fixed at sigmabig and its mean is fixed at 0.  If set to FALSE, the estimated mean and standard deviation of cluster 2 are free to vary.
#' @param penfactor factor for dirichlet penalization for cluster probabilities at each step of the EM algorithm. The larger this is, the smaller the Dirichlet penalization
#' @param fudge small constant added to cluster probabilies at each EM step to ensure stability
#' @param alpha represents Bonferroni-corrected significance threshold when cutoff="b"
#' @param FDR_BH represents FDR-corrected significance threshold when cutoff="f"
#' @param topvec a vector representing values for K such that posterior probabilities that the parameter for each experiment is within the set of K parameters having the largest absolute values are given.
#' @return A list of the top ranked experiments
#' @importFrom stats dnorm kmeans p.adjust pnorm qnorm quantile rnorm
#' @importFrom utils flush.console head
#' @export
#' @examples
#' truetheta <- c(rep(0,900),rnorm(100))
#' setheta <- pmax(rexp(1000,1),.1)
#' esttheta <- rnorm(length(truetheta),mean=truetheta,sd=setheta)
#' # just rank experiments that are significant at 5% FDR
#' stuff <- rankEM(esttheta,setheta,cutoff='f',FDR_BH=.05)
#' # rank all experiments (slower)
#' # stuff <- rankEM(esttheta,setheta,cutoff='f',FDR_BH=1)
rankEM <- function(betahat, sebeta, Jmin = 1, Jmax = 4, maxiter = 200, 
    tol = 1e-05, nsim = 10000, cutoff = 0.5, maxpar = 40000, 
    multiplestart = FALSE, sigmabig = 10, fixedcluster2 = TRUE, 
    penfactor = 5000, fudge = 0.001, alpha = 0.05, FDR_BH = 0.05, 
    topvec = c(10, 100, 1000, 10000)) 
{
    Jmin=Jmin+1
    Jmax=Jmax+1
    tokeep <- (1:length(betahat))[!is.na(betahat) & !is.na(sebeta)]
    betahat <- betahat[tokeep]
    sebeta <- sebeta[tokeep]
    x <- betahat/sebeta
    stuff <- choose_model(x, Jmin = Jmin, Jmax = Jmax, maxiter = maxiter, 
        tol = tol, multiplestart = multiplestart, sigmabig = sigmabig, 
        fixedcluster2 = fixedcluster2, penfactor = penfactor, 
        fudge = fudge)
    mu <- stuff$mu
    sigma <- stuff$sigma
    pi <- stuff$pi
    J <- length(mu)
    scaledpost <- getpost(x = x, mu = mu, sigma = sigma, pi = pi)
    mixprobs <- scaledpost$mixprobs
    unscaled_mu_post <- matrix(rep(sebeta, J), ncol = J, byrow = FALSE) * 
        scaledpost$mu_post
    unscaled_mu_overall <- apply(unscaled_mu_post * mixprobs, 
        1, sum)
    unscaled_sigma_post <- sebeta %o% scaledpost$sigma_post
    print("simulating from posterior - please wait")
    if (cutoff == "b") {
        thres <- qnorm(1 - 0.5 * alpha/nrow(unscaled_mu_post))
        tokeep <- (1:nrow(unscaled_mu_post))[abs(betahat)/sebeta > 
            thres]
        unscaled_mu_post <- unscaled_mu_post[tokeep, ]
        unscaled_sigma_post <- unscaled_sigma_post[tokeep, ]
        unscaled_mu_overall <- unscaled_mu_overall[tokeep]
        mixprobs <- mixprobs[tokeep, ]
        betahat <- betahat[tokeep]
        sebeta <- sebeta[tokeep]
        the.subset <- mixprobs[, 1]
        outputs <- simpost_new(unscaled_mu_post, unscaled_sigma_post, 
            mixprobs, nsim = nsim, cutoff = -0.1, topvec = topvec)
        rankedlist <- cbind(indexes = (1:length(the.subset)), 
            `original beta` = betahat, `se beta` = sebeta, `posterior mean` = unscaled_mu_overall, 
            `prob alternative` = 1 - mixprobs[, 1], `Expected rank` = outputs$mean_ranks, 
            SE_ranks = outputs$se_ranks, mean_topK = outputs$mean_topK)
        colnames(rankedlist)[(ncol(rankedlist) - length(topvec) + 
            1):ncol(rankedlist)] <- paste("prob top", topvec, 
            sep = "")
        rankedlist <- rankedlist[order(outputs$mean_ranks, decreasing = TRUE), 
            ]
        print("top hits are")
        print(head(rankedlist))
        rankedlist
        return(rankedlist)
    }
    if (cutoff == "f") {
        pvec <- 2 * (1 - pnorm(abs(betahat)/sebeta))
        a_pvec <- p.adjust(pvec, method = "BH")
        tokeep <- (1:nrow(unscaled_mu_post))[a_pvec < FDR_BH]
        unscaled_mu_post <- unscaled_mu_post[tokeep, ]
        unscaled_sigma_post <- unscaled_sigma_post[tokeep, ]
        unscaled_mu_overall <- unscaled_mu_overall[tokeep]
        mixprobs <- mixprobs[tokeep, ]
        betahat <- betahat[tokeep]
        sebeta <- sebeta[tokeep]
        the.subset <- mixprobs[, 1]
        outputs <- simpost_new(unscaled_mu_post, unscaled_sigma_post, 
            mixprobs, nsim = nsim, cutoff = -0.1, topvec = topvec)
        rankedlist <- cbind(indexes = (1:length(the.subset)), 
            `original beta` = betahat, `se beta` = sebeta, `posterior mean` = unscaled_mu_overall, 
            `prob alternative` = 1 - mixprobs[, 1], `Expected rank` = outputs$mean_ranks, 
            SE_ranks = outputs$se_ranks, mean_topK = outputs$mean_topK)
        colnames(rankedlist)[(ncol(rankedlist) - length(topvec) + 
            1):ncol(rankedlist)] <- paste("prob top", topvec, 
            sep = "")
        rankedlist <- rankedlist[order(outputs$mean_ranks, decreasing = TRUE), 
            ]
        print("top hits are")
        print(head(rankedlist))
        return(rankedlist)
    }
    the.subset <- mixprobs[, 1]
    maxnullprob <- 1 - cutoff
    maxnullprob <- min(maxnullprob, quantile(mixprobs[, 1], min(maxpar, 
        length(mixprobs[, 1]))/length(mixprobs[, 1])))
    cutoff <- 1 - maxnullprob
    outputs <- simpost_new(unscaled_mu_post, unscaled_sigma_post, 
        mixprobs, nsim = nsim, cutoff = cutoff, topvec = topvec)
    rankedlist <- cbind(indexes = (1:length(the.subset))[the.subset <= 
        (maxnullprob)], `original beta` = betahat[the.subset <= 
        (maxnullprob)], `se beta` = sebeta[the.subset <= (maxnullprob)], 
        `posterior mean` = unscaled_mu_overall[the.subset <= 
            (maxnullprob)], `prob alternative` = 1 - mixprobs[the.subset <= 
            (maxnullprob), 1], `Expected rank` = outputs$mean_ranks, 
        SE_ranks = outputs$se_ranks, mean_topK = outputs$mean_topK)
    colnames(rankedlist)[(ncol(rankedlist) - length(topvec) + 
        1):ncol(rankedlist)] <- paste("prob top", topvec, sep = "")
    rankedlist <- rankedlist[order(outputs$mean_ranks, decreasing = TRUE), 
        ]
    print("top hits are")
    print(head(rankedlist))
    return(rankedlist)
}

choose_model <- 

function (x, Jmin = 4, Jmax = 6, maxiter = 200, tol = 1e-05, 
    multiplestart = FALSE, verbose = TRUE, sigmabig = 10, fixedcluster2 = TRUE, 
    penfactor = 5, fudge = 0.001) 
{
    for (k in Jmin:Jmax) {
        results <- fit(x, k, maxiter = maxiter, tol = tol, multiplestart = multiplestart, 
            sigmabig = sigmabig, fixedcluster2 = fixedcluster2, 
            penfactor = penfactor, fudge = fudge)
        if (k == Jmin) {
            BICvec <- results$BIC
            mucurrent <- results$mu
            sigmacurrent <- results$sigma
            picurrent <- results$pi
            clustercurrent <- Jmin
        }
        if (k >= (Jmin + 1)) {
            BICvec <- c(BICvec, results$BIC)
            names(BICvec) <- c(Jmin:k)
            if (BICvec[(k - (Jmin - 1))] <= min(BICvec)) {
                mucurrent <- results$mu
                sigmacurrent <- results$sigma
                picurrent <- results$pi
                clustercurrent <- k
            }
        }
        if (verbose) {
            print(k)
            print(results$mu)
            print(results$sigma)
            print(results$pi)
            print(BICvec[length(BICvec)])
        }
    }
    return(list(nclust = clustercurrent, BICvec = BICvec, mu = mucurrent, 
        sigma = sigmacurrent, pi = picurrent))
}

fit <- 

function (x, J = 4, maxiter = 100, tol = 1e-05, multiplestart = FALSE, 
    fixedcluster2 = TRUE, sigmabig = 10, penfactor = 5e+05, fudge = 0.001) 
{
    muInt <- matrix(0, J, 2)
    sigInt <- matrix(0, J, 2)
    if (J == 1) {
        print("J must be larger than 1")
        return("")
    }
    if (fixedcluster2) {
        muInt[1:2, 1] <- c(0, 0)
        muInt[1:2, 2] <- c(0, 0)
        sigInt[1:2, 1] <- c(1, sigmabig)
        sigInt[1:2, 2] <- c(1, sigmabig)
    }
    if (!fixedcluster2) {
        muInt[1:2, 1] <- c(0, -Inf)
        muInt[1:2, 2] <- c(0, Inf)
        sigInt[1:2, 1] <- c(1, 1)
        sigInt[1:2, 2] <- c(1, Inf)
    }
    if (J > 2) {
        muInt[3:J, 1] <- rep(-Inf, J - 2)
        muInt[3:J, 2] <- rep(Inf, J - 2)
        sigInt[3:J, 1] <- rep(1, J - 2)
        sigInt[3:J, 2] <- rep(Inf, J - 2)
    }
    if (J <= 3) 
        mustartmat <- cbind(rep(0, J), rep(0, J))
    if (J > 3) {
         theseq <- seq(from = 1, to = (2 * (J - 2) - 1), length = J - 
            2)
        theseq <- theseq/(2 * (J - 2))
        mustartmat <- cbind(c(0, 0, kmeans(x, centers = quantile(x, 
            theseq))$centers), c(0, 0, sample(x, J - 2)), c(0, 
            0, sample(x, J - 2)))
    }
    L <- 1
    if (multiplestart) 
        L <- ncol(mustartmat)
    for (k in 1:L) {
        if (J <= 2) 
            starts <- cbind(pi = c(0.99, 0.001), mu = c(0, 0), 
                sig = c(1, 10))
        if (J > 2) 
            starts <- cbind(pi = c(0.99, 0.001, rep(1 - 0.99 - 
                0.001, J - 2)/(J - 2)), mu = c(0, 0, mustartmat[3:J, 
                k]), sig = c(1, 10, rep(1, (J - 2))))
        pen <- length(x)/penfactor
        p <- c(1, rep(0, J - 1))
        EMstep <- function(x, pi, mu, sig, pen, p, muInt, sigInt, 
            fudge = 0.001) {
            n <- length(x)
            J <- length(mu)
            phiMat <- matDens(x, mu, sig)
            denominators <- (phiMat %*% pi)[, 1]
            denominators[denominators < 10^-10] <- 10^-10
            G <- phiMat * ((1/denominators) %o% pi)
            tG <- t(G)
            numeratorsMu = (tG %*% x)
            zplusj <- rowSums(tG)
            zplusj[zplusj < 10^-10] <- 10^-10
            newMu <- (numeratorsMu/zplusj)
            newSig <- tG %*% (x * x)/zplusj - newMu^2
            newMu <- pmax(pmin(muInt[, 2], newMu), muInt[, 1])
            newSig <- pmax(pmin(sigInt[, 2], sqrt(newSig)), sigInt[, 
                1])
            newPi <- (pen * p + zplusj)/(sum(zplusj) + pen)
            newPi <- (newPi + fudge)
            newPi <- pmax(0, pmin(newPi, 1))
            newPi <- newPi/sum(newPi)
            if (any(is.na(newPi), is.na(newMu), is.na(newSig))) {
                print(pi)
                print(mu)
                print(sig)
            }
            phiMat_new <- matDens(x, newMu, newSig)
            denominators <- (phiMat_new %*% pi)[, 1]
            denominators[denominators < 10^-10] <- 10^-10
            if (fixedcluster2) 
                BIC <- -2 * sum(log(denominators)) + log(length(x)) * 
                  (3 * (J - 2) + 1)
            if (!fixedcluster2) 
                BIC <- -2 * sum(log(denominators)) + log(length(x)) * 
                  (3 * (J - 1))
            return(list(pi = newPi, mu = newMu, sigma = newSig, 
                BIC = BIC))
        }
        pi <- starts[, 1]
        mu <- starts[, 2]
        sigma <- starts[, 3]
        iter = 0
        converged = FALSE
        while (!converged && iter < maxiter) {
            res <- EMstep(x, pi, mu, sigma, pen = pen, p, muInt, 
                sigInt, fudge = fudge)
            distance <- sum((pi - res$pi)^2) + sum((mu - res$mu)^2) + 
                sum((sigma - res$sigma)^2)
            if (is.na(distance)) {
                print(starts)
                print(muInt)
                print(sigInt)
            }
            if (distance < tol) 
                converged <- TRUE
            iter <- iter + 1
            pi <- res$pi
            mu <- res$mu
            sigma <- res$sigma
            BIC <- res$BIC
        }
        if (k == 1) {
            pimat <- pi
            mumat <- mu
            sigmamat <- sigma
            BICvec <- BIC
        }
        if (k > 1) {
            pimat <- cbind(pimat, pi)
            mumat <- cbind(mumat, mu)
            sigmamat <- cbind(sigmamat, sigma)
            BICvec <- c(BICvec, BIC)
        }
    }
    if (L == 1) 
        return(list(pi = pimat, mu = mumat, sigma = sigmamat, 
            BIC = BICvec, converged = converged, nIter = iter))
    if (L > 1) 
        return(list(pi = pimat[, which.min(BICvec)], mu = mumat[, 
            which.min(BICvec)], sigma = sigmamat[, which.min(BICvec)], 
            BIC = BICvec[which.min(BICvec)]))
}

getpost <- function (x, mu, sigma, pi) 
{
    phiMat = matDens(x, mu, sigma)
    sigma_b <- pmax(sigma - 1, 0)
    denominators = (phiMat %*% pi)
    mixprobs <- t(t(phiMat) * pi)
    mixprobs <- mixprobs/matrix(rep(denominators, length(mu)), 
        nrow = length(denominators), byrow = FALSE)
    mu_post <- (x %o% (sigma_b^2) + matrix(rep(mu, length(x)), 
        nrow = length(x), byrow = TRUE))/(1 + matrix(rep(sigma_b^2, 
        length(x)), nrow = length(x), byrow = TRUE))
    sigma_post <- sqrt(sigma_b^2/(sigma_b^2 + 1))
    return(list(mixprobs = mixprobs, mu_post = mu_post, sigma_post = sigma_post))
}

matDens <- function (x, mu, sig) 
{
    n <- length(x)
    J <- length(mu)
    res <- matrix(x, n, J) - matrix(mu, n, J, byrow = TRUE)
    ss <- matrix(1/sig, n, J, byrow = TRUE)
    res <- res * ss
    res <- dnorm(0) * exp(-0.5 * res^2)
    res <- res * ss
    return(res)
}

simpost_new <- 
function (mu_post, sigma_post, mixprobs, nsim = NA, cutoff, topvec = c(10, 
    100)) 
{
    mu_post <- mu_post[mixprobs[, 1] <= (1 - cutoff), ]
    sigma_post <- sigma_post[mixprobs[, 1] <= (1 - cutoff), ]
    mixprobs <- mixprobs[mixprobs[, 1] <= (1 - cutoff), ]
    J <- ncol(mu_post)
    if (is.na(nsim)) 
        nsim <- floor(10^8/nrow(mu_post))
    ranks <- matrix(0, ncol = 0, nrow = nrow(mu_post))
    sum_ranks <- numeric(nrow(mu_post))
    sum_ranks2 <- numeric(nrow(mu_post))
    theN <- nrow(mu_post)
    sum_topK <- matrix(0, nrow = nrow(mu_post), ncol = length(topvec))
    for (k in 1:ceiling(nsim/100)) {
        randompart <- matrix(rnorm(100 * nrow(mu_post)), nrow = nrow(mu_post))
        mysample <- function(x) {
            sample(1:J, size = 100, replace = TRUE, prob = x)
        }
        components <- t(apply(mixprobs, 1, mysample))
        meanmat <- matrix(0, nrow = nrow(randompart), ncol = ncol(randompart))
        sdmat <- matrix(0, nrow = nrow(randompart), ncol = ncol(randompart))
        for (i in 1:J) {
            meanmat[components == i] <- matrix(rep(mu_post[, 
                i], 100), byrow = FALSE, ncol = 100)[components == 
                i]
            sdmat[components == i] <- matrix(rep(sigma_post[, 
                i], 100), byrow = FALSE, ncol = 100)[components == 
                i]
        }
        output <- randompart * sdmat + meanmat
        ranks <- apply(output, 2, function(x) {
            rank(abs(x))
        })
        sum_ranks <- sum_ranks + apply(ranks, 1, sum)
        sum_ranks2 <- sum_ranks2 + apply(ranks, 1, function(x) {
            sum(x^2)
        })
        for (i in 1:length(topvec)) {
            sum_topK[, i] <- sum_topK[, i] + apply(ranks > theN - 
                topvec[i], 1, sum)
        }
        flush.console()
        print(k * 100)
    }
    mean_ranks <- sum_ranks/nsim
    mean_topK <- sum_topK/nsim
    sd_ranks <- sqrt((sum_ranks2 - nsim * (mean_ranks^2))/(nsim - 
        1))
    se_topK <- sqrt(mean_topK * (1 - mean_topK))/sqrt(nsim)
    se_ranks <- sd_ranks/sqrt(nsim)
    return(list(mean_ranks = mean_ranks, sd_ranks = sd_ranks, 
        se_ranks = se_ranks, mean_topK = mean_topK, se_topK = se_topK))
}