MCMCSubclass <- function(indices, niter, subvox){
    ind <- indices[t(subvox),]
    nmask <- nrow(subvox)
    p <- matrix(0, nrow=nmask,ncol=ncol(indices))
    for(i in 1:8){
        a <- seq(1,(nmask*8)-7,by=8) + i - 1
        p <- p + ind[a,]
    }
    p.vox <- p/8/niter

    ind <- indices[t(subvox),]
    mode <- matrix(max.col(ind),nrow=8)
    class.modemean=cbind(colSums(mode==1),colSums(mode==2), colSums(mode==3))/8

    list(fussy.vox = p.vox, class.modemean=class.modemean)
}

updateY <- function(yobs, indices, subvox, mu, sigma)
    .Call("updateY", yobs, subvox, indices, mu, sigma)

getDenSub <- function(y, mu, sigma)
    .Call("getDenSub", y, mu, sigma)

updateIndices <- function(neighbors, nneigh, blocks, nblocks, k, indices, check,den){
    .Call("updateIndices", blocks, neighbors, nneigh, k, indices, check, den)
}

updateMusB <- function(Nj, S1j, sigma, eta, xi, y, nvert, k){
    oldtau <- 1/(sigma^2)
    newpre <- eta + oldtau * Nj
    newmean <- (eta*xi + oldtau * S1j) / newpre
    mu <- rnorm(k, mean=newmean, sd=1/sqrt(newpre))
    mu
}

updateSdsB <- function(Nj, ybar, S2j, mu, lambda, phi, y, nvert,k){
    shape <- Nj/2 + lambda
    mid <- S2j + Nj * (ybar - mu)^2
    scale <- 2 * phi / (2 + phi*mid)
    newtau <- rgamma(k, shape=shape, scale=scale)
    sigma <- 1/sqrt(newtau)
    sigma
}

yIndicesSummaries <- function(y, indices, k)
    structure(.Call("yIndicesSummaries", y, indices, k), names = c("Nj", "ybar", "S2j"))


mritc.bayes <- function(y, neighbors, blocks, sub, subvox,
                        spatialMat=(if(sub) matrix(c(2,0,-1,0,2,0,-1,0,2), nrow=3)
                                    else diag(1,3)),
                        beta=ifelse(sub, 0.3, 0.7), mu, sigma, niter=100, verbose){
    #r <- range(yobs)
    r <- 256
    xi <- r/2
    m1 <- 0.3
    eta <- m1/r^2
    lambda <- m1
    phi <- 1/r^2

    checkErrors(mu=mu, sigma=sigma, err=NULL)

    k <- length(mu)

    nneigh <- ncol(neighbors)
    nblocks <- length(blocks)
    neighbors <- structure(as.integer(neighbors), dim = dim(neighbors))

    y <- as.double(y)

    check <- matrix(as.double(getCheck(nneigh, k, beta, spatialMat)), ncol=k)

    if(sub==F){
        yunique <- sort(unique(y))
        n.yunique <- length(yunique)
        nvert <- length(y)
        ymatch <- match(y, yunique)
        indices <- initialIndices(y, nvert, mu, sigma, k, sub=FALSE)
    }

    else{
        yobs <- y
        nobs <- length(yobs)
        nvert <- nobs*8
        mu <- mu/8
        sigma <- sqrt(sigma^2)/8
        y <- rep(0, nvert)
        y[subvox] <- rep(yobs/8, 8)
        xi <- xi/8
        eta <- eta*64
        lambda <- lambda
        phi <- phi*64
        indices <- initialIndices(yobs, nvert, mu, sigma, k, sub=TRUE, subvox)
    }


    indicesMCMC <- matrix(0L, nrow=nrow(indices), ncol=ncol(indices))
    musave <- matrix(0, ncol=k, nrow=niter)
    sigmasave <- matrix(0, ncol=k, nrow=niter)

    for (i in 1:niter){
        if(sub==F)
            den <- getDen(yunique, n.yunique, ymatch, mu, sigma)
        else den <- getDenSub(y, mu, sigma)

        indices <- updateIndices(neighbors, nneigh, blocks, nblocks, k, indices,
                           check, den)

        indicesMCMC <- .Call("updateCounts", indicesMCMC, indices)

        if(sub==T){
            y <- updateY(yobs, indices, subvox, mu, sigma)
        }

        yIndicesSums <- yIndicesSummaries(y, indices, k)
        Nj <- yIndicesSums$Nj
        ybar <- yIndicesSums$ybar
        S2j <- yIndicesSums$S2j
        S1j <- Nj * ybar
        mu <- updateMusB(Nj, S1j, sigma, eta, xi, y, nvert, k)
        sigma <- updateSdsB(Nj, ybar, S2j, mu, lambda, phi, y, nvert,k)

        musave[i,] <- mu
        sigmasave[i,] <- sigma

        if (verbose && i %% 10 ==0)
            cat(paste("Iteration ", i, " has finished", "\n", sep=""))

    }
    if(sub==FALSE)
        prob <- indicesMCMC[-(nvert+1),]/niter
    else
        prob <- class <- MCMCSubclass(indicesMCMC[-(nvert+1),], niter, subvox)$fussy.vox
    list(prob=prob, mu=colMeans(musave), sigma=colMeans(sigmasave))
}

