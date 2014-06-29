## @knitr inputPump

pumpCode <- modelCode({ 
  for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta);
      lambda[i] <- theta[i]*t[i];
      x[i] ~ dpois(lambda[i])
  }
  alpha ~ dexp(1.0);
  beta ~ dgamma(0.1,1.0);
})

pumpConsts <- list(N = 10,
               t = c(94.3, 15.7, 62.9, 126, 5.24,
                 31.4, 1.05, 1.05, 2.1, 10.5))

pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

pumpInits <- list(alpha = 1, beta = 1,
              theta = rep(0.1, pumpConsts$N))

## @knitr explorePump

pump <- nimbleModel(code = pumpCode, name = 'pump', constants = pumpConsts,
                    data = pumpData, inits = pumpInits)

pump$getNodeNames()
pump$x
pump$alpha
pump$theta

## @knitr plotPump

plot(pump$graph)

## @knitr manipPump

pump$getDependencies(c('alpha', 'beta'))
pump$getDependencies(c('alpha', 'beta'), determOnly = TRUE)

set.seed(0) ## This makes the simulations here reproducible
calculate(pump, pump$getDependencies(c('alpha', 'beta'), determOnly = TRUE))
simulate(pump, 'theta')
pump$theta   ## the new theta values
pump$lambda  ## lambda hasn't been calculated yet
calculate(pump, pump$getDependencies(c('theta')))
pump$lambda  ## now it has

## @knitr compilePump

Cpump <- compileNimble(pump)
Cpump$theta

## @knitr mcmcPump
pumpSpec <- MCMCspec(pump, print = TRUE)
pumpSpec$addMonitors(c('alpha', 'beta', 'theta'))

pumpMCMC <- buildMCMC(pumpSpec)
CpumpMCMC <- compileNimble(pumpMCMC, project = pump)

niter <- 1000
set.seed(0)
CpumpMCMC(niter)

samples <- as.matrix(nfVar(CpumpMCMC, 'mvSamples'))

par(mfrow = c(1, 4), mai = c(.5, .5, .1, .2))
plot(samples[ , 'alpha'], type = 'l', xlab = 'iteration',
     ylab = expression(alpha))
plot(samples[ , 'beta'], type = 'l', xlab = 'iteration',
     ylab = expression(beta))
plot(samples[ , 'alpha'], samples[ , 'beta'], xlab = expression(alpha),
     ylab = expression(beta))
plot(samples[ , 'theta[1]'], type = 'l', xlab = 'iteration',
     ylab = expression(theta[1]))

acf(samples[, 'alpha']) ## plot autocorrelation of alpha sample
acf(samples[, 'beta'])  ## plot autocorrelation of beta  sample

## @knitr mcmcPump2

pumpSpec$addSampler('RW_block', list(targetNodes = c('alpha', 'beta'),
                                     adaptInterval = 100))
pumpMCMC2 <- buildMCMC(pumpSpec)
# need to reset the nimbleFunctions in order to add the new MCMC
CpumpNewMCMC <- compileNimble(pumpMCMC2, project  = pump, resetFunctions = TRUE)

set.seed(0);
CpumpNewMCMC(niter)
samplesNew <- as.matrix(nfVar(CpumpNewMCMC, 'mvSamples'))

par(mfrow = c(1, 4), mai = c(.5, .5, .1, .2))
plot(samplesNew[ , 'alpha'], type = 'l', xlab = 'iteration',
     ylab = expression(alpha))
plot(samplesNew[ , 'beta'], type = 'l', xlab = 'iteration',
     ylab = expression(beta))
plot(samplesNew[ , 'alpha'], samplesNew[ , 'beta'], xlab = expression(alpha),
     ylab = expression(beta))
plot(samplesNew[ , 'theta[1]'], type = 'l', xlab = 'iteration',
     ylab = expression(theta[1]))

acf(samplesNew[, 'alpha']) ## plot autocorrelation of alpha sample
acf(samplesNew[, 'beta'])  ## plot autocorrelation of beta  sample

## @knitr mcemPump

pump2 <- pump$newModel()

nodes <- pump2$getNodeNames(stochOnly = TRUE)

box = list( list(c('alpha','beta'), c(0, Inf)))

pumpMCEM <- buildMCEM(model = pump2, latentNodes = 'theta[1:10]',
                       boxConstraints = box)

pumpMLE <- pumpMCEM()
pumpMLE

## @knitr nfPump

simNodesMany <- nimbleFunction(
    setup = function(model, nodes) {
        mv <- modelValues(model)
        deps <- model$getDependencies(nodes)
        allNodes <- model$getNodeNames()
    },
    run = function(n = integer()) {
        resize(mv, n)
        for(i in 1:n) {
            simulate(model, nodes)
            calculate(model, deps)
            copy(from = model, nodes = allNodes, to = mv, rowTo = i, logProb = TRUE)
        }
    })

simNodesTheta1to5 <- simNodesMany(pump, 'theta[1:5]')

## @knitr runPumpSimsR
set.seed(0)  ## make the calculation repeatable
pump$alpha <- pumpMLE[1]
pump$beta <- pumpMLE[2]
calculate(pump, pump$getDependencies(c('alpha','beta'), determOnly = TRUE))
saveTheta <- pump$theta
simNodesTheta1to5(10)
nfVar(simNodesTheta1to5, 'mv')[['theta']][1:2]
nfVar(simNodesTheta1to5, 'mv')[['logProb_x']][1:2]

## @knitr runPumpSimsC
CsimNodesTheta1to5 <- compileNimble(simNodesTheta1to5,
                                    project  = pump, resetFunctions = TRUE)
Cpump$alpha <- pumpMLE[1]
Cpump$beta <- pumpMLE[2]
calculate(Cpump, Cpump$getDependencies(c('alpha','beta'), determOnly = TRUE))
Cpump$theta <- saveTheta

set.seed(0)
CsimNodesTheta1to5(10)
nfVar(CsimNodesTheta1to5, 'mv')[['theta']][1:2]
nfVar(CsimNodesTheta1to5, 'mv')[['logProb_x']][1:2]
