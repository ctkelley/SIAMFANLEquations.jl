"""
armijosc(xt, xc, ft, fc, d, residm, ItRules, derivative_is_old)
Line search for scalar equations. Read the notebook or print book
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
    fp=ItRules.fp
    dx=ItRules.dx
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
