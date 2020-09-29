"""
armijosc(xt, xc, ft, fc, d, residm, ItRules, derivative_is_old)
Line search for Newton's method. Read the notebook or print book
for the explanation. This is an internal function and I did not
design it to be hackable by the novice.
"""
function armijosc(xt, xc, ft, fc, d, residm, ItRules, derivative_is_old)
    idid = true
    alpha = 1.e-4
    iarm = -1
    lambda = 1.0
    lam0 = 0.0
    lamc = lambda
    lamm = lamc
    f=ItRules.f
#    fp=ItRules.fp
#    dx=ItRules.dx
    ResidC = residm
    armmax=ItRules.armmax
    armfix=ItRules.armfix
    #
    if derivative_is_old
       armmax=0
    end
    #
    #   Take the full step and, if happy, go home.
    #
    (xt, ft, residt) = UpdateIteration(xt, xc, ft, lambda, d, ItRules)
    armfail = residt > (1 - alpha * lambda) * residm
    iarm+=1
    #
    #
    # At this point I've taken a full step. I'll enter the loop only if
    # that full step has failed.
    #
    ffc = residt^2
    ff0 = residm^2
    ffm = ffc
    while armfail && iarm < armmax
        #
        #   At this point the full step has failed. Now it's time to be 
        #   serious about the line search.
        #
        lambda = update_lambda(iarm, armfix, lambda, lamc, ff0, ffc, ffm)
        (xt, ft, residt) = UpdateIteration(xt, xc, ft, lambda, d, ItRules)
        ffm = ffc
        ffc = residt^2
        iarm += 1
        armfail = residt > (1 - alpha * lambda) * residm
    end
    if iarm >= armmax 
        idid = false
    end
    return (ax = xt, afc = ft, resnorm = residt, aiarm = iarm, idid = idid)
end

function update_lambda(iarm, armfix, lambda, lamc, ff0, ffc, ffm)
       if iarm == 0 || armfix == true
            lambda = lambda * 0.5
        else
            lamm=lamc
            lamc=lambda
            lambda = parab3p(lamc, lamm, ff0, ffc, ffm)
        end
return lambda
end

"""
parab3p(lambdac, lambdam, ff0, ffc, ffm)

Three point parabolic line search.

input:\n
       lambdac = current steplength
       lambdam = previous steplength
       ff0 = value of || F(x_c) ||^2
       ffc = value of || F(x_c + lambdac d) ||^2
       ffm = value of || F(x_c + lambdam d) ||^2

output:\n
       lambdap = new value of lambda

internal parameters:\n
       sigma0 = .1, sigma1=.5, safeguarding bounds for the linesearch

You get here if cutting the steplength in half doesn't get you
sufficient decrease. Now you have three points and can build a parabolic
model. I do not like cubic models because they either need four points
or a derivative. 

So let's think about how this works. I cheat a bit and check the model
for negative curvature, which I don't want to see.

 The polynomial is

 p(lambda) = ff0 + (c1 lambda + c2 lambda^2)/d1

 d1 = (lambdac - lambdam)*lambdac*lambdam < 0
 So if c2 > 0 we have negative curvature and default to
      lambdap = sigma0 * lambda
 The logic is that negative curvature is telling us that
 the polynomial model is not helping much, so it looks better
 to take the smallest possible step. This is not what I did in the
 matlab code because I did it wrong. I have sinced fixed it.

 So (Students, listen up!) if c2 < 0 then all we gotta do is minimize
 (c1 lambda + c2 lambda^2)/d1 over [.1* lambdac, .5*lambdac]
 This means to MAXIMIZE c1 lambda + c2 lambda^2 becase d1 < 0.
 So I find the zero of the derivative and check the endpoints.

"""
function parab3p(lambdac, lambdam, ff0, ffc, ffm)
    #
    # internal parameters
    #
    sigma0 = 0.1
    sigma1 = 0.5
    #
    c2 = lambdam * (ffc - ff0) - lambdac * (ffm - ff0)
    if c2 >= 0
        #
        # Sanity check for negative curvature
        #
        lambdap = sigma0 * lambdac
    else
        #
        # It's a convex parabola, so use calculus!
        #
        c1 = lambdac * lambdac * (ffm - ff0) - lambdam * lambdam * (ffc - ff0)
        lambdap = -c1 * 0.5 / c2
        #
        lambdaup = sigma1 * lambdac
        lambdadown = sigma0 * lambdac
        lambdap = max(lambdadown, min(lambdaup, lambdap))
    end
end
