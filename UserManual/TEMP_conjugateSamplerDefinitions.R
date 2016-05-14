
sampler_conjugate_dbeta <- nimbleFunction(contains = sampler_BASE, setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    calcNodesDeterm <- model$getDependencies(target, determOnly = TRUE)
    target_nodeFunctionList <- nimbleFunctionList(node_stoch_dbeta)
    target_nodeFunctionList[[1]] <- model$nodeFunctions[[target]]
    dep_dbern_nodeNames <- control$dep_dbern
    dep_dbern_nfs <- nimbleFunctionList(node_stoch_dbern)
    for (i in seq_along(dep_dbern_nodeNames)) {
        dep_dbern_nfs[[i]] <- model$nodeFunctions[[dep_dbern_nodeNames[i]]]
    }
    dep_dbin_nodeNames <- control$dep_dbin
    dep_dbin_nfs <- nimbleFunctionList(node_stoch_dbin)
    for (i in seq_along(dep_dbin_nodeNames)) {
        dep_dbin_nfs[[i]] <- model$nodeFunctions[[dep_dbin_nodeNames[i]]]
    }
    dep_dnegbin_nodeNames <- control$dep_dnegbin
    dep_dnegbin_nfs <- nimbleFunctionList(node_stoch_dnegbin)
    for (i in seq_along(dep_dnegbin_nodeNames)) {
        dep_dnegbin_nfs[[i]] <- model$nodeFunctions[[dep_dnegbin_nodeNames[i]]]
    }
}, run = function() {
    prior_shape1 <- nfMethod(target_nodeFunctionList[[1]], 'get_shape1')()
    prior_shape2 <- nfMethod(target_nodeFunctionList[[1]], 'get_shape2')()
    declare(dep_dbern_values, double(1, length(dep_dbern_nfs)))
    declare(dep_dbern_prob, double(1, length(dep_dbern_nfs)))
    for (i in seq_along(dep_dbern_nfs)) {
        dep_dbern_values[i] <- nfMethod(dep_dbern_nfs[[i]], 'get_value')()
        dep_dbern_prob[i] <- nfMethod(dep_dbern_nfs[[i]], 'get_prob')()
    }
    declare(dep_dbin_values, double(1, length(dep_dbin_nfs)))
    declare(dep_dbin_prob, double(1, length(dep_dbin_nfs)))
    declare(dep_dbin_size, double(1, length(dep_dbin_nfs)))
    for (i in seq_along(dep_dbin_nfs)) {
        dep_dbin_values[i] <- nfMethod(dep_dbin_nfs[[i]], 'get_value')()
        dep_dbin_prob[i] <- nfMethod(dep_dbin_nfs[[i]], 'get_prob')()
        dep_dbin_size[i] <- nfMethod(dep_dbin_nfs[[i]], 'get_size')()
    }
    declare(dep_dnegbin_values, double(1, length(dep_dnegbin_nfs)))
    declare(dep_dnegbin_prob, double(1, length(dep_dnegbin_nfs)))
    declare(dep_dnegbin_size, double(1, length(dep_dnegbin_nfs)))
    for (i in seq_along(dep_dnegbin_nfs)) {
        dep_dnegbin_values[i] <- nfMethod(dep_dnegbin_nfs[[i]], 'get_value')()
        dep_dnegbin_prob[i] <- nfMethod(dep_dnegbin_nfs[[i]], 'get_prob')()
        dep_dnegbin_size[i] <- nfMethod(dep_dnegbin_nfs[[i]], 'get_size')()
    }
    contribution_shape1 <- 0
    contribution_shape2 <- 0
    for (i in seq_along(dep_dbern_nfs)) {
        contribution_shape1 <- contribution_shape1 + dep_dbern_values[i]
        contribution_shape2 <- contribution_shape2 + (1 - dep_dbern_values[i])
    }
    for (i in seq_along(dep_dbin_nfs)) {
        contribution_shape1 <- contribution_shape1 + dep_dbin_values[i]
        contribution_shape2 <- contribution_shape2 + (dep_dbin_size[i] - dep_dbin_values[i])
    }
    for (i in seq_along(dep_dnegbin_nfs)) {
        contribution_shape1 <- contribution_shape1 + dep_dnegbin_size[i]
        contribution_shape2 <- contribution_shape2 + dep_dnegbin_values[i]
    }
    newValue <- rbeta(1, shape1 = prior_shape1 + contribution_shape1, shape2 = prior_shape2 + contribution_shape2)
    model[[target]] <<- newValue
    calculate(model, calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
}, methods = list(getPosteriorLogDensity = function() {
    prior_shape1 <- nfMethod(target_nodeFunctionList[[1]], 'get_shape1')()
    prior_shape2 <- nfMethod(target_nodeFunctionList[[1]], 'get_shape2')()
    declare(dep_dbern_values, double(1, length(dep_dbern_nfs)))
    declare(dep_dbern_prob, double(1, length(dep_dbern_nfs)))
    for (i in seq_along(dep_dbern_nfs)) {
        dep_dbern_values[i] <- nfMethod(dep_dbern_nfs[[i]], 'get_value')()
        dep_dbern_prob[i] <- nfMethod(dep_dbern_nfs[[i]], 'get_prob')()
    }
    declare(dep_dbin_values, double(1, length(dep_dbin_nfs)))
    declare(dep_dbin_prob, double(1, length(dep_dbin_nfs)))
    declare(dep_dbin_size, double(1, length(dep_dbin_nfs)))
    for (i in seq_along(dep_dbin_nfs)) {
        dep_dbin_values[i] <- nfMethod(dep_dbin_nfs[[i]], 'get_value')()
        dep_dbin_prob[i] <- nfMethod(dep_dbin_nfs[[i]], 'get_prob')()
        dep_dbin_size[i] <- nfMethod(dep_dbin_nfs[[i]], 'get_size')()
    }
    declare(dep_dnegbin_values, double(1, length(dep_dnegbin_nfs)))
    declare(dep_dnegbin_prob, double(1, length(dep_dnegbin_nfs)))
    declare(dep_dnegbin_size, double(1, length(dep_dnegbin_nfs)))
    for (i in seq_along(dep_dnegbin_nfs)) {
        dep_dnegbin_values[i] <- nfMethod(dep_dnegbin_nfs[[i]], 'get_value')()
        dep_dnegbin_prob[i] <- nfMethod(dep_dnegbin_nfs[[i]], 'get_prob')()
        dep_dnegbin_size[i] <- nfMethod(dep_dnegbin_nfs[[i]], 'get_size')()
    }
    contribution_shape1 <- 0
    contribution_shape2 <- 0
    for (i in seq_along(dep_dbern_nfs)) {
        contribution_shape1 <- contribution_shape1 + dep_dbern_values[i]
        contribution_shape2 <- contribution_shape2 + (1 - dep_dbern_values[i])
    }
    for (i in seq_along(dep_dbin_nfs)) {
        contribution_shape1 <- contribution_shape1 + dep_dbin_values[i]
        contribution_shape2 <- contribution_shape2 + (dep_dbin_size[i] - dep_dbin_values[i])
    }
    for (i in seq_along(dep_dnegbin_nfs)) {
        contribution_shape1 <- contribution_shape1 + dep_dnegbin_size[i]
        contribution_shape2 <- contribution_shape2 + dep_dnegbin_values[i]
    }
    targetValue <- model[[target]]
    posteriorLogDensity <- dbeta(targetValue, shape1 = prior_shape1 + contribution_shape1, shape2 = prior_shape2 + contribution_shape2, log = 1)
    returnType(double())
    return(posteriorLogDensity)
}, reset = function() {
}), where = getLoadingNamespace())




sampler_conjugate_ddirch <- nimbleFunction(contains = sampler_BASE, setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    calcNodesDeterm <- model$getDependencies(target, determOnly = TRUE)
    target_nodeFunctionList <- nimbleFunctionList(node_stoch_ddirch)
    target_nodeFunctionList[[1]] <- model$nodeFunctions[[target]]
    dep_dmulti_nodeNames <- control$dep_dmulti
    dep_dmulti_nfs <- nimbleFunctionList(node_stoch_dmulti)
    for (i in seq_along(dep_dmulti_nodeNames)) {
        dep_dmulti_nfs[[i]] <- model$nodeFunctions[[dep_dmulti_nodeNames[i]]]
    }
    d <- max(determineNodeIndexSizes(target))
}, run = function() {
    prior_alpha <- nfMethod(target_nodeFunctionList[[1]], 'get_alpha')()
    declare(dep_dmulti_values, double(2, c(length(dep_dmulti_nfs), d)))
    declare(dep_dmulti_prob, double(2, c(length(dep_dmulti_nfs), d)))
    for (i in seq_along(dep_dmulti_nfs)) {
        dep_dmulti_values[i, 1:d] <- nfMethod(dep_dmulti_nfs[[i]], 'get_value')()
        dep_dmulti_prob[i, 1:d] <- nfMethod(dep_dmulti_nfs[[i]], 'get_prob')()
    }
    contribution_alpha <- numeric(length = d)
    for (i in seq_along(dep_dmulti_nfs)) {
        contribution_alpha <- contribution_alpha + dep_dmulti_values[i, 1:d]
    }
    newValue <- rdirch(1, alpha = prior_alpha + contribution_alpha)
    model[[target]] <<- newValue
    calculate(model, calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
}, methods = list(getPosteriorLogDensity = function() {
    prior_alpha <- nfMethod(target_nodeFunctionList[[1]], 'get_alpha')()
    declare(dep_dmulti_values, double(2, c(length(dep_dmulti_nfs), d)))
    declare(dep_dmulti_prob, double(2, c(length(dep_dmulti_nfs), d)))
    for (i in seq_along(dep_dmulti_nfs)) {
        dep_dmulti_values[i, 1:d] <- nfMethod(dep_dmulti_nfs[[i]], 'get_value')()
        dep_dmulti_prob[i, 1:d] <- nfMethod(dep_dmulti_nfs[[i]], 'get_prob')()
    }
    contribution_alpha <- numeric(length = d)
    for (i in seq_along(dep_dmulti_nfs)) {
        contribution_alpha <- contribution_alpha + dep_dmulti_values[i, 1:d]
    }
    targetValue <- model[[target]]
    posteriorLogDensity <- ddirch(targetValue, alpha = prior_alpha + contribution_alpha, log = 1)
    returnType(double())
    return(posteriorLogDensity)
}, reset = function() {
}), where = getLoadingNamespace())




sampler_conjugate_dgamma <- nimbleFunction(contains = sampler_BASE, setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    calcNodesDeterm <- model$getDependencies(target, determOnly = TRUE)
    target_nodeFunctionList <- nimbleFunctionList(node_stoch_dgamma)
    target_nodeFunctionList[[1]] <- model$nodeFunctions[[target]]
    dep_dpois_nodeNames <- control$dep_dpois
    dep_dpois_nfs <- nimbleFunctionList(node_stoch_dpois)
    for (i in seq_along(dep_dpois_nodeNames)) {
        dep_dpois_nfs[[i]] <- model$nodeFunctions[[dep_dpois_nodeNames[i]]]
    }
    dep_dnorm_nodeNames <- control$dep_dnorm
    dep_dnorm_nfs <- nimbleFunctionList(node_stoch_dnorm)
    for (i in seq_along(dep_dnorm_nodeNames)) {
        dep_dnorm_nfs[[i]] <- model$nodeFunctions[[dep_dnorm_nodeNames[i]]]
    }
    dep_dlnorm_nodeNames <- control$dep_dlnorm
    dep_dlnorm_nfs <- nimbleFunctionList(node_stoch_dlnorm)
    for (i in seq_along(dep_dlnorm_nodeNames)) {
        dep_dlnorm_nfs[[i]] <- model$nodeFunctions[[dep_dlnorm_nodeNames[i]]]
    }
    dep_dgamma_nodeNames <- control$dep_dgamma
    dep_dgamma_nfs <- nimbleFunctionList(node_stoch_dgamma)
    for (i in seq_along(dep_dgamma_nodeNames)) {
        dep_dgamma_nfs[[i]] <- model$nodeFunctions[[dep_dgamma_nodeNames[i]]]
    }
    dep_dexp_nodeNames <- control$dep_dexp
    dep_dexp_nfs <- nimbleFunctionList(node_stoch_dexp)
    for (i in seq_along(dep_dexp_nodeNames)) {
        dep_dexp_nfs[[i]] <- model$nodeFunctions[[dep_dexp_nodeNames[i]]]
    }
}, run = function() {
    prior_shape <- nfMethod(target_nodeFunctionList[[1]], 'get_shape')()
    prior_rate <- nfMethod(target_nodeFunctionList[[1]], 'get_rate')()
    declare(dep_dpois_values, double(1, length(dep_dpois_nfs)))
    declare(dep_dpois_lambda, double(1, length(dep_dpois_nfs)))
    for (i in seq_along(dep_dpois_nfs)) {
        dep_dpois_values[i] <- nfMethod(dep_dpois_nfs[[i]], 'get_value')()
        dep_dpois_lambda[i] <- nfMethod(dep_dpois_nfs[[i]], 'get_lambda')()
    }
    declare(dep_dnorm_values, double(1, length(dep_dnorm_nfs)))
    declare(dep_dnorm_tau, double(1, length(dep_dnorm_nfs)))
    declare(dep_dnorm_mean, double(1, length(dep_dnorm_nfs)))
    for (i in seq_along(dep_dnorm_nfs)) {
        dep_dnorm_values[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_value')()
        dep_dnorm_tau[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_tau')()
        dep_dnorm_mean[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_mean')()
    }
    declare(dep_dlnorm_values, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dlnorm_taulog, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dlnorm_meanlog, double(1, length(dep_dlnorm_nfs)))
    for (i in seq_along(dep_dlnorm_nfs)) {
        dep_dlnorm_values[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_value')()
        dep_dlnorm_taulog[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_taulog')()
        dep_dlnorm_meanlog[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_meanlog')()
    }
    declare(dep_dgamma_values, double(1, length(dep_dgamma_nfs)))
    declare(dep_dgamma_rate, double(1, length(dep_dgamma_nfs)))
    declare(dep_dgamma_shape, double(1, length(dep_dgamma_nfs)))
    for (i in seq_along(dep_dgamma_nfs)) {
        dep_dgamma_values[i] <- nfMethod(dep_dgamma_nfs[[i]], 'get_value')()
        dep_dgamma_rate[i] <- nfMethod(dep_dgamma_nfs[[i]], 'get_rate')()
        dep_dgamma_shape[i] <- nfMethod(dep_dgamma_nfs[[i]], 'get_shape')()
    }
    declare(dep_dexp_values, double(1, length(dep_dexp_nfs)))
    declare(dep_dexp_rate, double(1, length(dep_dexp_nfs)))
    for (i in seq_along(dep_dexp_nfs)) {
        dep_dexp_values[i] <- nfMethod(dep_dexp_nfs[[i]], 'get_value')()
        dep_dexp_rate[i] <- nfMethod(dep_dexp_nfs[[i]], 'get_rate')()
    }
    declare(dep_dpois_offset, double(1, length(dep_dpois_nfs)))
    declare(dep_dpois_coeff, double(1, length(dep_dpois_nfs)))
    declare(dep_dnorm_offset, double(1, length(dep_dnorm_nfs)))
    declare(dep_dnorm_coeff, double(1, length(dep_dnorm_nfs)))
    declare(dep_dlnorm_offset, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dlnorm_coeff, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dgamma_offset, double(1, length(dep_dgamma_nfs)))
    declare(dep_dgamma_coeff, double(1, length(dep_dgamma_nfs)))
    declare(dep_dexp_offset, double(1, length(dep_dexp_nfs)))
    declare(dep_dexp_coeff, double(1, length(dep_dexp_nfs)))
    model[[target]] <<- 0
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dpois_nfs)) {
        dep_dpois_offset[i] <- nfMethod(dep_dpois_nfs[[i]], 'get_lambda')()
    }
    for (i in seq_along(dep_dnorm_nfs)) {
        dep_dnorm_offset[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_tau')()
    }
    for (i in seq_along(dep_dlnorm_nfs)) {
        dep_dlnorm_offset[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_taulog')()
    }
    for (i in seq_along(dep_dgamma_nfs)) {
        dep_dgamma_offset[i] <- nfMethod(dep_dgamma_nfs[[i]], 'get_rate')()
    }
    for (i in seq_along(dep_dexp_nfs)) {
        dep_dexp_offset[i] <- nfMethod(dep_dexp_nfs[[i]], 'get_rate')()
    }
    model[[target]] <<- 1
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dpois_nfs)) {
        dep_dpois_coeff[i] <- nfMethod(dep_dpois_nfs[[i]], 'get_lambda')() - dep_dpois_offset[i]
    }
    for (i in seq_along(dep_dnorm_nfs)) {
        dep_dnorm_coeff[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_tau')() - dep_dnorm_offset[i]
    }
    for (i in seq_along(dep_dlnorm_nfs)) {
        dep_dlnorm_coeff[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_taulog')() - dep_dlnorm_offset[i]
    }
    for (i in seq_along(dep_dgamma_nfs)) {
        dep_dgamma_coeff[i] <- nfMethod(dep_dgamma_nfs[[i]], 'get_rate')() - dep_dgamma_offset[i]
    }
    for (i in seq_along(dep_dexp_nfs)) {
        dep_dexp_coeff[i] <- nfMethod(dep_dexp_nfs[[i]], 'get_rate')() - dep_dexp_offset[i]
    }
    contribution_shape <- 0
    contribution_rate <- 0
    for (i in seq_along(dep_dpois_nfs)) {
        contribution_shape <- contribution_shape + dep_dpois_values[i]
        contribution_rate <- contribution_rate + dep_dpois_coeff[i]
    }
    for (i in seq_along(dep_dnorm_nfs)) {
        contribution_shape <- contribution_shape + 1/2
        contribution_rate <- contribution_rate + dep_dnorm_coeff[i]/2 * (dep_dnorm_values[i] - dep_dnorm_mean[i])^2
    }
    for (i in seq_along(dep_dlnorm_nfs)) {
        contribution_shape <- contribution_shape + 1/2
        contribution_rate <- contribution_rate + dep_dlnorm_coeff[i]/2 * (log(dep_dlnorm_values[i]) - dep_dlnorm_meanlog[i])^2
    }
    for (i in seq_along(dep_dgamma_nfs)) {
        contribution_shape <- contribution_shape + dep_dgamma_shape[i]
        contribution_rate <- contribution_rate + dep_dgamma_coeff[i] * dep_dgamma_values[i]
    }
    for (i in seq_along(dep_dexp_nfs)) {
        contribution_shape <- contribution_shape + 1
        contribution_rate <- contribution_rate + dep_dexp_coeff[i] * dep_dexp_values[i]
    }
    newValue <- rgamma(1, shape = prior_shape + contribution_shape, scale = 1/(prior_rate + contribution_rate))
    model[[target]] <<- newValue
    calculate(model, calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
}, methods = list(getPosteriorLogDensity = function() {
    prior_shape <- nfMethod(target_nodeFunctionList[[1]], 'get_shape')()
    prior_rate <- nfMethod(target_nodeFunctionList[[1]], 'get_rate')()
    declare(dep_dpois_values, double(1, length(dep_dpois_nfs)))
    declare(dep_dpois_lambda, double(1, length(dep_dpois_nfs)))
    for (i in seq_along(dep_dpois_nfs)) {
        dep_dpois_values[i] <- nfMethod(dep_dpois_nfs[[i]], 'get_value')()
        dep_dpois_lambda[i] <- nfMethod(dep_dpois_nfs[[i]], 'get_lambda')()
    }
    declare(dep_dnorm_values, double(1, length(dep_dnorm_nfs)))
    declare(dep_dnorm_tau, double(1, length(dep_dnorm_nfs)))
    declare(dep_dnorm_mean, double(1, length(dep_dnorm_nfs)))
    for (i in seq_along(dep_dnorm_nfs)) {
        dep_dnorm_values[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_value')()
        dep_dnorm_tau[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_tau')()
        dep_dnorm_mean[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_mean')()
    }
    declare(dep_dlnorm_values, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dlnorm_taulog, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dlnorm_meanlog, double(1, length(dep_dlnorm_nfs)))
    for (i in seq_along(dep_dlnorm_nfs)) {
        dep_dlnorm_values[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_value')()
        dep_dlnorm_taulog[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_taulog')()
        dep_dlnorm_meanlog[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_meanlog')()
    }
    declare(dep_dgamma_values, double(1, length(dep_dgamma_nfs)))
    declare(dep_dgamma_rate, double(1, length(dep_dgamma_nfs)))
    declare(dep_dgamma_shape, double(1, length(dep_dgamma_nfs)))
    for (i in seq_along(dep_dgamma_nfs)) {
        dep_dgamma_values[i] <- nfMethod(dep_dgamma_nfs[[i]], 'get_value')()
        dep_dgamma_rate[i] <- nfMethod(dep_dgamma_nfs[[i]], 'get_rate')()
        dep_dgamma_shape[i] <- nfMethod(dep_dgamma_nfs[[i]], 'get_shape')()
    }
    declare(dep_dexp_values, double(1, length(dep_dexp_nfs)))
    declare(dep_dexp_rate, double(1, length(dep_dexp_nfs)))
    for (i in seq_along(dep_dexp_nfs)) {
        dep_dexp_values[i] <- nfMethod(dep_dexp_nfs[[i]], 'get_value')()
        dep_dexp_rate[i] <- nfMethod(dep_dexp_nfs[[i]], 'get_rate')()
    }
    declare(dep_dpois_offset, double(1, length(dep_dpois_nfs)))
    declare(dep_dpois_coeff, double(1, length(dep_dpois_nfs)))
    declare(dep_dnorm_offset, double(1, length(dep_dnorm_nfs)))
    declare(dep_dnorm_coeff, double(1, length(dep_dnorm_nfs)))
    declare(dep_dlnorm_offset, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dlnorm_coeff, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dgamma_offset, double(1, length(dep_dgamma_nfs)))
    declare(dep_dgamma_coeff, double(1, length(dep_dgamma_nfs)))
    declare(dep_dexp_offset, double(1, length(dep_dexp_nfs)))
    declare(dep_dexp_coeff, double(1, length(dep_dexp_nfs)))
    model[[target]] <<- 0
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dpois_nfs)) {
        dep_dpois_offset[i] <- nfMethod(dep_dpois_nfs[[i]], 'get_lambda')()
    }
    for (i in seq_along(dep_dnorm_nfs)) {
        dep_dnorm_offset[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_tau')()
    }
    for (i in seq_along(dep_dlnorm_nfs)) {
        dep_dlnorm_offset[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_taulog')()
    }
    for (i in seq_along(dep_dgamma_nfs)) {
        dep_dgamma_offset[i] <- nfMethod(dep_dgamma_nfs[[i]], 'get_rate')()
    }
    for (i in seq_along(dep_dexp_nfs)) {
        dep_dexp_offset[i] <- nfMethod(dep_dexp_nfs[[i]], 'get_rate')()
    }
    model[[target]] <<- 1
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dpois_nfs)) {
        dep_dpois_coeff[i] <- nfMethod(dep_dpois_nfs[[i]], 'get_lambda')() - dep_dpois_offset[i]
    }
    for (i in seq_along(dep_dnorm_nfs)) {
        dep_dnorm_coeff[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_tau')() - dep_dnorm_offset[i]
    }
    for (i in seq_along(dep_dlnorm_nfs)) {
        dep_dlnorm_coeff[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_taulog')() - dep_dlnorm_offset[i]
    }
    for (i in seq_along(dep_dgamma_nfs)) {
        dep_dgamma_coeff[i] <- nfMethod(dep_dgamma_nfs[[i]], 'get_rate')() - dep_dgamma_offset[i]
    }
    for (i in seq_along(dep_dexp_nfs)) {
        dep_dexp_coeff[i] <- nfMethod(dep_dexp_nfs[[i]], 'get_rate')() - dep_dexp_offset[i]
    }
    contribution_shape <- 0
    contribution_rate <- 0
    for (i in seq_along(dep_dpois_nfs)) {
        contribution_shape <- contribution_shape + dep_dpois_values[i]
        contribution_rate <- contribution_rate + dep_dpois_coeff[i]
    }
    for (i in seq_along(dep_dnorm_nfs)) {
        contribution_shape <- contribution_shape + 1/2
        contribution_rate <- contribution_rate + dep_dnorm_coeff[i]/2 * (dep_dnorm_values[i] - dep_dnorm_mean[i])^2
    }
    for (i in seq_along(dep_dlnorm_nfs)) {
        contribution_shape <- contribution_shape + 1/2
        contribution_rate <- contribution_rate + dep_dlnorm_coeff[i]/2 * (log(dep_dlnorm_values[i]) - dep_dlnorm_meanlog[i])^2
    }
    for (i in seq_along(dep_dgamma_nfs)) {
        contribution_shape <- contribution_shape + dep_dgamma_shape[i]
        contribution_rate <- contribution_rate + dep_dgamma_coeff[i] * dep_dgamma_values[i]
    }
    for (i in seq_along(dep_dexp_nfs)) {
        contribution_shape <- contribution_shape + 1
        contribution_rate <- contribution_rate + dep_dexp_coeff[i] * dep_dexp_values[i]
    }
    targetValue <- model[[target]]
    posteriorLogDensity <- dgamma(targetValue, shape = prior_shape + contribution_shape, scale = 1/(prior_rate + contribution_rate), log = 1)
    returnType(double())
    return(posteriorLogDensity)
}, reset = function() {
}), where = getLoadingNamespace())




sampler_conjugate_dnorm <- nimbleFunction(contains = sampler_BASE, setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    calcNodesDeterm <- model$getDependencies(target, determOnly = TRUE)
    target_nodeFunctionList <- nimbleFunctionList(node_stoch_dnorm)
    target_nodeFunctionList[[1]] <- model$nodeFunctions[[target]]
    dep_dnorm_nodeNames <- control$dep_dnorm
    dep_dnorm_nfs <- nimbleFunctionList(node_stoch_dnorm)
    for (i in seq_along(dep_dnorm_nodeNames)) {
        dep_dnorm_nfs[[i]] <- model$nodeFunctions[[dep_dnorm_nodeNames[i]]]
    }
    dep_dlnorm_nodeNames <- control$dep_dlnorm
    dep_dlnorm_nfs <- nimbleFunctionList(node_stoch_dlnorm)
    for (i in seq_along(dep_dlnorm_nodeNames)) {
        dep_dlnorm_nfs[[i]] <- model$nodeFunctions[[dep_dlnorm_nodeNames[i]]]
    }
}, run = function() {
    prior_mean <- nfMethod(target_nodeFunctionList[[1]], 'get_mean')()
    prior_tau <- nfMethod(target_nodeFunctionList[[1]], 'get_tau')()
    declare(dep_dnorm_values, double(1, length(dep_dnorm_nfs)))
    declare(dep_dnorm_mean, double(1, length(dep_dnorm_nfs)))
    declare(dep_dnorm_tau, double(1, length(dep_dnorm_nfs)))
    for (i in seq_along(dep_dnorm_nfs)) {
        dep_dnorm_values[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_value')()
        dep_dnorm_mean[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_mean')()
        dep_dnorm_tau[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_tau')()
    }
    declare(dep_dlnorm_values, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dlnorm_meanlog, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dlnorm_taulog, double(1, length(dep_dlnorm_nfs)))
    for (i in seq_along(dep_dlnorm_nfs)) {
        dep_dlnorm_values[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_value')()
        dep_dlnorm_meanlog[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_meanlog')()
        dep_dlnorm_taulog[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_taulog')()
    }
    declare(dep_dnorm_offset, double(1, length(dep_dnorm_nfs)))
    declare(dep_dnorm_coeff, double(1, length(dep_dnorm_nfs)))
    declare(dep_dlnorm_offset, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dlnorm_coeff, double(1, length(dep_dlnorm_nfs)))
    model[[target]] <<- 0
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dnorm_nfs)) {
        dep_dnorm_offset[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_mean')()
    }
    for (i in seq_along(dep_dlnorm_nfs)) {
        dep_dlnorm_offset[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_meanlog')()
    }
    model[[target]] <<- 1
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dnorm_nfs)) {
        dep_dnorm_coeff[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_mean')() - dep_dnorm_offset[i]
    }
    for (i in seq_along(dep_dlnorm_nfs)) {
        dep_dlnorm_coeff[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_meanlog')() - dep_dlnorm_offset[i]
    }
    contribution_mean <- 0
    contribution_tau <- 0
    for (i in seq_along(dep_dnorm_nfs)) {
        contribution_mean <- contribution_mean + dep_dnorm_coeff[i] * (dep_dnorm_values[i] - dep_dnorm_offset[i]) * dep_dnorm_tau[i]
        contribution_tau <- contribution_tau + dep_dnorm_coeff[i]^2 * dep_dnorm_tau[i]
    }
    for (i in seq_along(dep_dlnorm_nfs)) {
        contribution_mean <- contribution_mean + dep_dlnorm_coeff[i] * (log(dep_dlnorm_values[i]) - dep_dlnorm_offset[i]) * dep_dlnorm_taulog[i]
        contribution_tau <- contribution_tau + dep_dlnorm_coeff[i]^2 * dep_dlnorm_taulog[i]
    }
    newValue <- rnorm(1, mean = (prior_mean * prior_tau + contribution_mean)/(prior_tau + contribution_tau), sd = (prior_tau + contribution_tau)^(-0.5))
    model[[target]] <<- newValue
    calculate(model, calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
}, methods = list(getPosteriorLogDensity = function() {
    prior_mean <- nfMethod(target_nodeFunctionList[[1]], 'get_mean')()
    prior_tau <- nfMethod(target_nodeFunctionList[[1]], 'get_tau')()
    declare(dep_dnorm_values, double(1, length(dep_dnorm_nfs)))
    declare(dep_dnorm_mean, double(1, length(dep_dnorm_nfs)))
    declare(dep_dnorm_tau, double(1, length(dep_dnorm_nfs)))
    for (i in seq_along(dep_dnorm_nfs)) {
        dep_dnorm_values[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_value')()
        dep_dnorm_mean[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_mean')()
        dep_dnorm_tau[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_tau')()
    }
    declare(dep_dlnorm_values, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dlnorm_meanlog, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dlnorm_taulog, double(1, length(dep_dlnorm_nfs)))
    for (i in seq_along(dep_dlnorm_nfs)) {
        dep_dlnorm_values[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_value')()
        dep_dlnorm_meanlog[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_meanlog')()
        dep_dlnorm_taulog[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_taulog')()
    }
    declare(dep_dnorm_offset, double(1, length(dep_dnorm_nfs)))
    declare(dep_dnorm_coeff, double(1, length(dep_dnorm_nfs)))
    declare(dep_dlnorm_offset, double(1, length(dep_dlnorm_nfs)))
    declare(dep_dlnorm_coeff, double(1, length(dep_dlnorm_nfs)))
    model[[target]] <<- 0
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dnorm_nfs)) {
        dep_dnorm_offset[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_mean')()
    }
    for (i in seq_along(dep_dlnorm_nfs)) {
        dep_dlnorm_offset[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_meanlog')()
    }
    model[[target]] <<- 1
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dnorm_nfs)) {
        dep_dnorm_coeff[i] <- nfMethod(dep_dnorm_nfs[[i]], 'get_mean')() - dep_dnorm_offset[i]
    }
    for (i in seq_along(dep_dlnorm_nfs)) {
        dep_dlnorm_coeff[i] <- nfMethod(dep_dlnorm_nfs[[i]], 'get_meanlog')() - dep_dlnorm_offset[i]
    }
    contribution_mean <- 0
    contribution_tau <- 0
    for (i in seq_along(dep_dnorm_nfs)) {
        contribution_mean <- contribution_mean + dep_dnorm_coeff[i] * (dep_dnorm_values[i] - dep_dnorm_offset[i]) * dep_dnorm_tau[i]
        contribution_tau <- contribution_tau + dep_dnorm_coeff[i]^2 * dep_dnorm_tau[i]
    }
    for (i in seq_along(dep_dlnorm_nfs)) {
        contribution_mean <- contribution_mean + dep_dlnorm_coeff[i] * (log(dep_dlnorm_values[i]) - dep_dlnorm_offset[i]) * dep_dlnorm_taulog[i]
        contribution_tau <- contribution_tau + dep_dlnorm_coeff[i]^2 * dep_dlnorm_taulog[i]
    }
    targetValue <- model[[target]]
    posteriorLogDensity <- dnorm(targetValue, mean = (prior_mean * prior_tau + contribution_mean)/(prior_tau + contribution_tau), sd = (prior_tau + contribution_tau)^(-0.5), log = 1)
    returnType(double())
    return(posteriorLogDensity)
}, reset = function() {
}), where = getLoadingNamespace())




sampler_conjugate_dmnorm <- nimbleFunction(contains = sampler_BASE, setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    calcNodesDeterm <- model$getDependencies(target, determOnly = TRUE)
    target_nodeFunctionList <- nimbleFunctionList(node_stoch_dmnorm)
    target_nodeFunctionList[[1]] <- model$nodeFunctions[[target]]
    dep_dmnorm_nodeNames <- control$dep_dmnorm
    dep_dmnorm_nfs <- nimbleFunctionList(node_stoch_dmnorm)
    for (i in seq_along(dep_dmnorm_nodeNames)) {
        dep_dmnorm_nfs[[i]] <- model$nodeFunctions[[dep_dmnorm_nodeNames[i]]]
    }
    d <- max(determineNodeIndexSizes(target))
}, run = function() {
    prior_prec <- nfMethod(target_nodeFunctionList[[1]], 'get_prec')()
    prior_mean <- nfMethod(target_nodeFunctionList[[1]], 'get_mean')()
    declare(dep_dmnorm_values, double(2, c(length(dep_dmnorm_nfs), d)))
    declare(dep_dmnorm_mean, double(2, c(length(dep_dmnorm_nfs), d)))
    declare(dep_dmnorm_prec, double(3, c(length(dep_dmnorm_nfs), d, d)))
    for (i in seq_along(dep_dmnorm_nfs)) {
        dep_dmnorm_values[i, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_value')()
        dep_dmnorm_mean[i, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_mean')()
        dep_dmnorm_prec[i, 1:d, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_prec')()
    }
    declare(dep_dmnorm_offset, double(2, c(length(dep_dmnorm_nfs), d)))
    declare(dep_dmnorm_coeff, double(3, c(length(dep_dmnorm_nfs), d, d)))
    model[[target]] <<- model[[target]] * 0
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dmnorm_nfs)) {
        dep_dmnorm_offset[i, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_mean')()
    }
    for (sizeIndex in 1:d) {
        unitVector <- model[[target]] * 0
        unitVector[sizeIndex] <- 1
        model[[target]] <<- unitVector
        calculate(model, calcNodesDeterm)
        for (i in seq_along(dep_dmnorm_nfs)) {
            dep_dmnorm_coeff[i, 1:d, sizeIndex] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_mean')() - dep_dmnorm_offset[i, 1:d]
        }
    }
    contribution_prec <- array(dim = c(d, d))
    contribution_mean <- numeric(length = d)
    for (i in seq_along(dep_dmnorm_nfs)) {
        contribution_prec <- contribution_prec + calc_dmnormConjugacyContributions(dep_dmnorm_coeff[i, 1:d, 1:d], dep_dmnorm_prec[i, 1:d, 1:d], 2)
        contribution_mean <- contribution_mean + (calc_dmnormConjugacyContributions(dep_dmnorm_coeff[i, 1:d, 1:d], dep_dmnorm_prec[i, 1:d, 1:d], 1) %*% asCol(dep_dmnorm_values[i, 1:d] - dep_dmnorm_offset[i, 1:d]))[, 1]
    }
    R <- chol(prior_prec + contribution_prec)
    A <- prior_prec %*% asCol(prior_mean) + asCol(contribution_mean)
    mu <- backsolve(R, forwardsolve(t(R), A))[, 1]
    newValue <- rmnorm_chol(1, mean = mu, cholesky = R, prec_param = 1)
    model[[target]] <<- newValue
    calculate(model, calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
}, methods = list(getPosteriorLogDensity = function() {
    prior_prec <- nfMethod(target_nodeFunctionList[[1]], 'get_prec')()
    prior_mean <- nfMethod(target_nodeFunctionList[[1]], 'get_mean')()
    declare(dep_dmnorm_values, double(2, c(length(dep_dmnorm_nfs), d)))
    declare(dep_dmnorm_mean, double(2, c(length(dep_dmnorm_nfs), d)))
    declare(dep_dmnorm_prec, double(3, c(length(dep_dmnorm_nfs), d, d)))
    for (i in seq_along(dep_dmnorm_nfs)) {
        dep_dmnorm_values[i, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_value')()
        dep_dmnorm_mean[i, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_mean')()
        dep_dmnorm_prec[i, 1:d, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_prec')()
    }
    declare(dep_dmnorm_offset, double(2, c(length(dep_dmnorm_nfs), d)))
    declare(dep_dmnorm_coeff, double(3, c(length(dep_dmnorm_nfs), d, d)))
    model[[target]] <<- model[[target]] * 0
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dmnorm_nfs)) {
        dep_dmnorm_offset[i, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_mean')()
    }
    for (sizeIndex in 1:d) {
        unitVector <- model[[target]] * 0
        unitVector[sizeIndex] <- 1
        model[[target]] <<- unitVector
        calculate(model, calcNodesDeterm)
        for (i in seq_along(dep_dmnorm_nfs)) {
            dep_dmnorm_coeff[i, 1:d, sizeIndex] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_mean')() - dep_dmnorm_offset[i, 1:d]
        }
    }
    contribution_prec <- array(dim = c(d, d))
    contribution_mean <- numeric(length = d)
    for (i in seq_along(dep_dmnorm_nfs)) {
        contribution_prec <- contribution_prec + calc_dmnormConjugacyContributions(dep_dmnorm_coeff[i, 1:d, 1:d], dep_dmnorm_prec[i, 1:d, 1:d], 2)
        contribution_mean <- contribution_mean + (calc_dmnormConjugacyContributions(dep_dmnorm_coeff[i, 1:d, 1:d], dep_dmnorm_prec[i, 1:d, 1:d], 1) %*% asCol(dep_dmnorm_values[i, 1:d] - dep_dmnorm_offset[i, 1:d]))[, 1]
    }
    R <- chol(prior_prec + contribution_prec)
    A <- prior_prec %*% asCol(prior_mean) + asCol(contribution_mean)
    mu <- backsolve(R, forwardsolve(t(R), A))[, 1]
    targetValue <- model[[target]]
    posteriorLogDensity <- dmnorm_chol(targetValue, mean = mu, cholesky = R, prec_param = 1, log = 1)
    returnType(double())
    return(posteriorLogDensity)
}, reset = function() {
}), where = getLoadingNamespace())




sampler_conjugate_dwish <- nimbleFunction(contains = sampler_BASE, setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    calcNodesDeterm <- model$getDependencies(target, determOnly = TRUE)
    target_nodeFunctionList <- nimbleFunctionList(node_stoch_dwish)
    target_nodeFunctionList[[1]] <- model$nodeFunctions[[target]]
    dep_dmnorm_nodeNames <- control$dep_dmnorm
    dep_dmnorm_nfs <- nimbleFunctionList(node_stoch_dmnorm)
    for (i in seq_along(dep_dmnorm_nodeNames)) {
        dep_dmnorm_nfs[[i]] <- model$nodeFunctions[[dep_dmnorm_nodeNames[i]]]
    }
    d <- max(determineNodeIndexSizes(target))
}, run = function() {
    prior_R <- nfMethod(target_nodeFunctionList[[1]], 'get_R')()
    prior_df <- nfMethod(target_nodeFunctionList[[1]], 'get_df')()
    declare(dep_dmnorm_values, double(2, c(length(dep_dmnorm_nfs), d)))
    declare(dep_dmnorm_prec, double(3, c(length(dep_dmnorm_nfs), d, d)))
    declare(dep_dmnorm_mean, double(2, c(length(dep_dmnorm_nfs), d)))
    for (i in seq_along(dep_dmnorm_nfs)) {
        dep_dmnorm_values[i, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_value')()
        dep_dmnorm_prec[i, 1:d, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_prec')()
        dep_dmnorm_mean[i, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_mean')()
    }
    declare(dep_dmnorm_offset, double(3, c(length(dep_dmnorm_nfs), d, d)))
    declare(dep_dmnorm_coeff, double(3, c(length(dep_dmnorm_nfs), d, d)))
    I <- identityMatrix(d)
    model[[target]] <<- I
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dmnorm_nfs)) {
        dep_dmnorm_offset[i, 1:d, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_prec')()
    }
    model[[target]] <<- I * 2
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dmnorm_nfs)) {
        dep_dmnorm_coeff[i, 1:d, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_prec')()
    }
    for (i in seq_along(dep_dmnorm_nfs)) {
        dep_dmnorm_coeff[i, 1:d, 1:d] <- dep_dmnorm_coeff[i, 1:d, 1:d] - dep_dmnorm_offset[i, 1:d, 1:d]
        dep_dmnorm_offset[i, 1:d, 1:d] <- dep_dmnorm_offset[i, 1:d, 1:d] - dep_dmnorm_coeff[i, 1:d, 1:d]
    }
    contribution_R <- array(dim = c(d, d))
    contribution_df <- 0
    for (i in seq_along(dep_dmnorm_nfs)) {
        contribution_R <- contribution_R + asCol(dep_dmnorm_values[i, 1:d] - dep_dmnorm_mean[i, 1:d]) %*% (asRow(dep_dmnorm_values[i, 1:d] - dep_dmnorm_mean[i, 1:d]) %*% dep_dmnorm_coeff[i, 1:d, 1:d])
        contribution_df <- contribution_df + 1
    }
    newValue <- rwish_chol(1, cholesky = chol(prior_R + contribution_R), df = prior_df + contribution_df, scale_param = 0)
    model[[target]] <<- newValue
    calculate(model, calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
}, methods = list(getPosteriorLogDensity = function() {
    prior_R <- nfMethod(target_nodeFunctionList[[1]], 'get_R')()
    prior_df <- nfMethod(target_nodeFunctionList[[1]], 'get_df')()
    declare(dep_dmnorm_values, double(2, c(length(dep_dmnorm_nfs), d)))
    declare(dep_dmnorm_prec, double(3, c(length(dep_dmnorm_nfs), d, d)))
    declare(dep_dmnorm_mean, double(2, c(length(dep_dmnorm_nfs), d)))
    for (i in seq_along(dep_dmnorm_nfs)) {
        dep_dmnorm_values[i, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_value')()
        dep_dmnorm_prec[i, 1:d, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_prec')()
        dep_dmnorm_mean[i, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_mean')()
    }
    declare(dep_dmnorm_offset, double(3, c(length(dep_dmnorm_nfs), d, d)))
    declare(dep_dmnorm_coeff, double(3, c(length(dep_dmnorm_nfs), d, d)))
    I <- identityMatrix(d)
    model[[target]] <<- I
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dmnorm_nfs)) {
        dep_dmnorm_offset[i, 1:d, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_prec')()
    }
    model[[target]] <<- I * 2
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dep_dmnorm_nfs)) {
        dep_dmnorm_coeff[i, 1:d, 1:d] <- nfMethod(dep_dmnorm_nfs[[i]], 'get_prec')()
    }
    for (i in seq_along(dep_dmnorm_nfs)) {
        dep_dmnorm_coeff[i, 1:d, 1:d] <- dep_dmnorm_coeff[i, 1:d, 1:d] - dep_dmnorm_offset[i, 1:d, 1:d]
        dep_dmnorm_offset[i, 1:d, 1:d] <- dep_dmnorm_offset[i, 1:d, 1:d] - dep_dmnorm_coeff[i, 1:d, 1:d]
    }
    contribution_R <- array(dim = c(d, d))
    contribution_df <- 0
    for (i in seq_along(dep_dmnorm_nfs)) {
        contribution_R <- contribution_R + asCol(dep_dmnorm_values[i, 1:d] - dep_dmnorm_mean[i, 1:d]) %*% (asRow(dep_dmnorm_values[i, 1:d] - dep_dmnorm_mean[i, 1:d]) %*% dep_dmnorm_coeff[i, 1:d, 1:d])
        contribution_df <- contribution_df + 1
    }
    targetValue <- model[[target]]
    posteriorLogDensity <- dwish_chol(targetValue, cholesky = chol(prior_R + contribution_R), df = prior_df + contribution_df, scale_param = 0, log = 1)
    returnType(double())
    return(posteriorLogDensity)
}, reset = function() {
}), where = getLoadingNamespace())




