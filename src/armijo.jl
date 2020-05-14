"""
armijo(fc, d, xm, fm, f, h, fp, armmax, derivative_is_old)

Line search for scalar equations. Read the notebook or print book
for the explanation. This is an internal function and I did not
design it to be hackable by the novice.
"""
function armijo(fc, d, xm, fm, f, h, fp, armmax, derivative_is_old)
    idid=true
    alpha = 1.e-4
    iarm = -1
    lambda = 1.0
    x = xm
    lam0=0.0
    lamc=lambda
    lamm=lamc
    armfail=abs(fc) > (1 - alpha * lambda) * abs(fm)
#
# If I have an old derivative I will not tolerate a failure in
# the line search. 
#
# jflag tells me that I had to refresh the derivative
# liarm is the counter that = iarm unless I refresh the Jacobian
#
    jflag=false
    liarm=-1
    while armfail && iarm < armmax
        #
        #   At this point fp = f'(xm) then it's time to be serious 
        #   about the line  search.
        #
        x = xm + lambda * d
        fc = f(x)
        iarm += 1
        liarm += 1
        armfail=abs(fc) > (1 - alpha * lambda) * abs(fm)
        if armfail && derivative_is_old
        df = fpeval_newton(xm, f, fm, fp, h)
        dfold = df
        d = -fm / df
        x = xm + lambda * d
        fc = f(x)
        armfail=abs(fc) > (1 - alpha * lambda) * abs(fm)
        derivative_is_old = false
        lambda=1.0
        else
        lambda = lambda * .5
        end
    end
    if iarm >= armmax
       idid=false
       println("Linesearch failure")
    end
    return (ax = x, afc = fc, aiarm = iarm, adfo = derivative_is_old, ad = d,
            idid=idid)
end

