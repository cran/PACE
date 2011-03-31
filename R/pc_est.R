pc_est = function(y, tt, mu, xcov, sigma, noeig, error = 1, method = "CE",
                  shrink = 0, out1, out21, regular = 0, rho = "CV", verbose = "on"){

    r = getEigens(xcov, out1, out21, noeig,1)
    lambda = r$lambda
    phi = r$phi
    eigens = r$eigen
    noeig = r$noeig
    rm(r)

    #calculate the fitted covariance
    #based on lambda and eigens
    xcovfit = eigens%*%diag(lambda, length(lambda), length(lambda))%*%t(eigens)
    r = getScores(y, tt, mu, phi, lambda, sigma, noeig, error, method, shrink, out1, regular, rho, verbose)
    list(xi_est = r$xi_est, xi_var = r$xi_var, lambda = lambda, phi = phi, 
         eigen = eigens, noeig = noeig, xcovfit = xcovfit, y_predOrig = r$y_predOrig, 
         rho_opt = r$rho_opt, sig1 = r$sig1)
 
}