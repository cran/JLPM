
plotAsso <- function(x, xlim = NULL)
{
    sharedtype <- x$sharedtype
    nonlin <- x$nonlin
    nea <- sum(x$idea)
    
    if(sharedtype == 1)
    {
        nomRE <- x$Names$Xnames[x$idea==1]
        simRE <- rangeRE(x)$RE

        if(nonlin[1] == 0)
        {
            asso <- simRE * x$b[sum(x$N[1:3]) + 1:x$N[4]]
        }
        else
        {
            if(nonlin[1] == 1)
            {
                avt <- sum(x$N[1:3])
                asso <- matrix(0, nrow(simRE), ncol(simRE))
                for(q in 1:nea)
                {
                    asso[q, ] <- simRE[q, ] * x$b[avt + 1] + simRE[q, ]^2 * x$b[avt + 2] + simRE[q, ]^3 * x$b[avt + 3]

                    avt <- avt + 3
                }
            }
        }

        par(mfrow = n2mfrow(nea))
        for(q in 1:nea)
            plot(simRE[q, ], asso[q, ], type = "p", col = 1, xlab = nomRE[q], ylab = paste("asso(", nomRE[q], ")", sep = ""))

        res <- list(x = simRE, y = asso)
    }

    if(sharedtype %in% c(2, 4))
    {
        simCL <- rangeCL(x)$CL

        if(nonlin[2] == 0)
        {
            asso <- simCL * x$b[sum(x$N[1:3]) + 1]
        }
        else
        {
            if(nonlin[2] == 1)
            {
                avt <- sum(x$N[1:3])
                asso <- simCL * x$b[avt + 1] + simCL^2 * x$b[avt + 2] + simCL^3 * x$b[avt + 3]
            }
        }

            plot(simCL, asso, type = "p", col = 1, xlab = "CL", ylab = "asso(CL)")

        res <- list(x = simCL, y = asso)
    }

    if(sharedtype %in% c(3, 4))
    {
        simCS <- rangeCS(x)$CS

        avt <- sum(x$N[1:3])
        if(sharedtype == 4)
        {
            if(nonlin[2] == 0)
                avt <- avt + 1
            else if(nonlin[2] == 1)
                avt <- avt + 3         
        }

        if(nonlin[3] == 0)
        {
            asso <- simCS * x$b[avt + 1]
        }
        else
        {
            if(nonlin[3] == 1)
            {
                asso <- simCS * x$b[avt + 1] + simCS^2 * x$b[avt + 2] + simCS^3 * x$b[avt + 3]
            }
        }

            plot(simCS, asso, type = "p", col = 1, xlab = "CS", ylab = "asso(CS)")

        if(sharedtype == 4)
            res <- list(CL = res, CS = list(x = simCS, y = asso))
        else
            res <- list(x = simCS, y = asso)
    }


    
    return(invisible(res))
}



### fct rangeRE

rangeRE <- function(x)
{
    nea <- sum(x$idea)
    nMC <- x$nMC
    if((nea == 0) | (nMC == 0)) stop("nothing to compute in rangeRE")

    ## MC sequence
    sequnif <- spacefillr::generate_sobol_owen_set(nMC, nea)
    seqMC <- apply(sequnif, 2, qnorm)

    ## cholesky of the RE
    ch <- matrix(0, nea, nea)
    if(nea > 1)
    {
        ch[upper.tri(ch, diag = TRUE)] <- c(1, x$b[sum(x$N[1:6]) + 1:x$N[7]])
    }
    else
    {
        ch[1, 1] <- 1
    }
    ch <- t(ch)

    ## simulated RE used in computation
    simRE <- apply(seqMC, 1, function(u) ch %*% u)
    
    ## range
    range_simRE <- apply(simRE, 1, range)

    return(list(range = range_simRE, RE = simRE))
}

rangeCL <- function(x)
{
    nea <- sum(x$idea)
    nMC <- x$nMC
    var.time <- as.character(x$call$var.time)
    
    beta <- c(0, x$b[sum(x$N[1:4]) + 1:x$N[5]])
    nvar <- length(x$Names$Xvar)
    d <- as.data.frame(matrix(0, 100, nvar))
    d[, var.time] <- seq(x$rangeTlong[1], x$rangeTlong[2], length.out = 100)

    Xfixed <- model.matrix(x$form$fixed, data = d)
    mu <- Xfixed %*% beta

    zu <- 0
    if((nea > 0) & (nMC > 0))
    {
        u <- rangeRE(x)$RE
        Xrandom <- model.matrix(x$form$random, data = d)
        zu <- Xrandom %*% u
    }

    cl <- sweep(zu, 1, mu, "+")
    
    return(list(range = range(cl), CL = cl))
}


rangeCS <- function(x)
{
    nea <- sum(x$idea)
    nMC <- x$nMC
    var.time <- as.character(x$call$var.time)
    
    beta <- x$b[sum(x$N[1:4]) + 1:x$N[5]]
    nvar <- length(x$Names$Xvar)
    d <- as.data.frame(matrix(0, 100, nvar))
    d[, var.time] <- seq(x$rangeTlong[1], x$rangeTlong[2], length.out = 100)

    ddt <- derivMat(x$form$fixed, x$form$random, d, var.time)

    #Xfixed <- model.matrix(x$form$fixed, data = d)
    Xfixed <- ddt[, 1:length(beta), drop = FALSE]
    mu <- Xfixed %*% beta

    zu <- 0
    if((nea > 0) & (nMC > 0))
    {
        u <- rangeRE(x)$RE
        #Xrandom <- model.matrix(x$form$random, data = d)
        Xrandom <- ddt[, x$N[5] + 1:nea, drop = FALSE]
        zu <- Xrandom %*% u
    }

    cs <- sweep(zu, 1, mu, "+")
    
    return(list(range = range(cs), CS = cs))
}
