library(assist)
library(dplyr)

logit <- function(p) {
    log(p/(1-p))
}

fit_splines <- function(values, groups, times) {
    # Fits a periodic spline with a random effects model to the `values`
    # where `times` is in hours, the period is fixed to 24 hours,
    # and where groupings of the random effects are specified by `groups`

    # returns a summary results object as well as the actual fit object itself

    df <- data.frame(value = values, time = times, group = groups)
    df$t <- (df$time %% 24) / 24 # Rescale time to [0,1)

    # Approximate a starting value for the fixed effects
    # This has to be put into a global variable since `snm` forces
    # evaluation in the global environment for some reason
    mesor_guess <<- mean(values)

    # Fit the spline `f` with random effects in amplitude, phase, and mesor
    # amplitude is constrained to be positive by using `exp(logAmp)`
    # phase is constrained to be within (-0.5,0.5) by using `alogit(phi)-0.5`
    # Random effects for phi and logAmp should be average 0 across all groups
    # so the overall amplitude and phase are encoded in the function `f`
    fit <- snm(value ~ exp(logAmp)*f(t - (alogit(phi) - 0.5)) + mesor,
               func=f(u)~list(periodic(u)),
               #func=f(u)~list(~sin(2*pi*u)+cos(2*pi*u)-1, lspline(u, type="sine0")),
               fixed=list(mesor~1),
               random=list(logAmp~1, phi~1, mesor~1),
               group=~group,
               start=c(mesor_guess),
               spar="m",
               control=list(
                    converg="PRSS",
                    rkpk.control=list(
                        limnla=c(-5,3) #Min and max values for lambda
                                       # In my experience, too low values are no good
                                       # (i.e. overfits). Default is c(-10,3)
                    )
                ),
               data=df)

    # Extract the shape of the overall fitted curve
    # compute it at 4 points per hour
    # as well as a 95% confidence interval
    curves <- data.frame(u=seq(from=0,to=1, length.out=24*4+1))
    curves.ci <- intervals(fit, newdata=curves)
    curves$fit <- curves.ci$fit[,1]
    curves$pstd <- curves.ci$pstd

    # Time of peak of overall fit in hours since time 0
    fit_peak <- curves$u[which.max(curves$fit)] * 24
    fit_trough <- curves$u[which.min(curves$fit)] * 24

    # Amplitude of overall fit (peak-to-trough height)
    fit_amplitude <- max(curves$fit) - min(curves$fit)

    # Compute phases from the 'phi' value
    phases <- fit_peak + logit(fit$coefficients$random$group[,'phi'])

    # A `t`-statistic for the curve being non-zero.
    # But it's not clear what the DoF for it are
    # and whether there's multiple-hypothesis testing
    # due to checking all the timepoints
    t_stat <- max( abs(curves$fit / curves$pstd) )

    # Random effects structure
    random_effs_structure <- as.matrix(fit$nlmeObj$modelStruct$reStruct$group)

    # Gather summary results
    fit_summary <- summary(fit)
    nlme_summary <- summary(fit$nlmeObj)
    results <- list(
        fit_peak_time = fit_peak,
        fit_trough_time = fit_trough,
        fit_amplitude = fit_amplitude,
        fit_mesor = fit$coefficients$fixed['mesor'],
        mesor_stderr = nlme_summary$tTable[, 'Std.Error'],
        re_logAmp_sd = sqrt(random_effs_structure['logAmp', 'logAmp'])*fit$sigma,
        re_phi_sd = sqrt(random_effs_structure['phi', 'phi'])*fit$sigma,
        re_mesor_sd = sqrt(random_effs_structure['mesor', 'mesor'])*fit$sigma,
        lambda = fit$lambda, # Smothing parameter
        sigma = fit$sigma, # Noise SD
        iteration_times = fit$iteration$times,
        iteration_converg = fit$iteration$converg,
        funcDf = fit$funcDf,
        log_likelihood = nlme_summary$logLik,
        AIC = nlme_summary$AIC,
        BIC = nlme_summary$BIC,
        t_stat = t_stat
    )

    return(list(
         summary = results,
         fit = fit,
         random_effects = fit$coefficients$random$group,
         random_effects_structure = random_effs_structure,
         curves = curves
    ))
}
