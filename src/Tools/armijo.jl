"""
armijosc(fc, d, xm, fm, f, h, fp, armmax, armfix, derivative_is_old)

Line search for scalar equations. Read the notebook or print book
for the explanation. This is an internal function and I did not
design it to be hackable by the novice.
"""
function armijosc(fc, d, xm, fm, f, h, fp, armmax, armfix,
               derivative_is_old)
    idid = true
    alpha = 1.e-4
    iarm = -1
    lambda = 1.0
    x = xm
    lam0 = 0.0
    lamc = lambda
    lamm = lamc
    #
    #   Take the full step and, if happy, go home.
    #
    x = xm + lambda * d
    fc = f(x)
    armfail = norm(fc) > (1 - alpha * lambda) * norm(fm)
    #
    # If I have an old derivative I will not tolerate a failure in
    # the line search. 
    #
    # jflag tells me that I had to refresh the derivative
    # liarm is the counter that = iarm unless I refresh the Jacobian
    #
    jflag = false
    if armfail && derivative_is_old
        df = fpeval_newton(xm, f, fm, fp, h)
        d = -fm / df
        x = xm + lambda * d
        fc = f(x)
        armfail = norm(fc) > (1 - alpha * lambda) * norm(fm)
        derivative_is_old = false
        jflag = true
    end
    liarm = 0
    iarm = 0
    #
    # At this point I've taken a full step. I'll enter the loop only if
    # that full step has failed.
    #
    ffc = norm(fc)^2
    ff0 = norm(fm)^2
    ffm = ffc
    while armfail && iarm < armmax
        #
        #   At this point fp = f'(xm) then it's time to be serious 
        #   about the line  search.
        #
        if iarm == 0 || armfix == true
            lambda = lambda * 0.5
        else
            lamm=lamc
            lamc=lambda
            lambda = parab3p(lamc, lamm, ff0, ffc, ffm)
        end
        x = xm + lambda * d
        fc = f(x)
        ffm = ffc
        ffc = norm(fc)^2
        iarm += 1
        liarm += 1
        armfail = norm(fc) > (1 - alpha * lambda) * norm(fm)
    end
    if iarm >= armmax 
        idid = false
    end
    return (ax = x, afc = fc, aiarm = iarm, 
            adfo = derivative_is_old, ad = d, idid = idid)
end
