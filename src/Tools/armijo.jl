"""
armijosc(fc, d, xm, residm, ItRules, derivative_is_old)

Line search for scalar equations. Read the notebook or print book
for the explanation. This is an internal function and I did not
design it to be hackable by the novice.
"""
function armijosc(fc, d, xm, residm, ItRules, derivative_is_old)
    idid = true
    alpha = 1.e-4
    iarm = -1
    lambda = 1.0
    x = xm
    lam0 = 0.0
    lamc = lambda
    lamm = lamc
    newjac = 0
    f=ItRules.f
    fp=ItRules.fp
    dx=ItRules.dx
    armmax=ItRules.armmax
    #
    # If I have an old derivative I will not tolerate a failure in
    # the line search. 
    #
    if derivative_is_old
    armmax=0
    end
    armfix=ItRules.armfix
    #
    #   Take the full step and, if happy, go home.
    #
#    x = xm + lambda * d
#    fc = f(x)
#    residc=norm(fc)
    (x, residc, fc) = UpdateIteration(x, xm, lambda, d, ItRules)
    armfail = residc > (1 - alpha * lambda) * residm
    iarm+=1
    #
    #
    # At this point I've taken a full step. I'll enter the loop only if
    # that full step has failed.
    #
    ffc = residc^2
    ff0 = residm^2
    ffm = ffc
    while armfail && iarm < armmax
        #
        #   At this point fp = f'(xm) then it's time to be serious 
        #   about the line  search.
        #
        lambda = update_lambda(iarm, armfix, lambda, lamc, ff0, ffc, ffm)
        x = xm + lambda * d
        fc = f(x)
        residc=norm(fc)
        ffm = ffc
        ffc = residc^2
        iarm += 1
        armfail = residc > (1 - alpha * lambda) * residm
    end
    if iarm >= armmax 
        idid = false
    end
    return (ax = x, afc = fc, aiarm = iarm, newjac = newjac,
            adfo = derivative_is_old, ad = d, idid = idid)
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
