function knl_continue(FFUN, qdata, pval, nval, x, x0, lambda)
#
# preallocated arrays from qdata
#
    FS = qdata.FS
    FPS = qdata.FPS
    n = length(x0)
    KDS = nkl_init(n, "gmres")
#
# precomputed data for the function, coming from qdata
#
    fdata = qdata.fdata
#
# storage for xdot and xold
#
    xdot= qdata.xdot
    xold= qdata.xold
#
# range for the parameter
#
    dlam = qdata.dlam
    lambdamax = qdata.lambdamax
    lambda = lambda + dlam
#
# functions to update the parameter and collect the data for
# the bifurcation diagram
#
    setlam=qdata.setlam
    bif_update=qdata.bif_update
#
# 
#
    dlamm1=1.0/dlam
    idid=true
    lambdaz=lambda
    while lambda <= lambdamax  && idid
#
# setlam informs FFUN about xold, xdot, and lambda by updating 
# qdata.fdata. This is the biggest, but not the only, part of this 
# deal that is not for general use.
#
        setlam(qdata, lambda, xdot, xold)
        debug && println("lambda = $lambda")
#
# When I send fdata to FFUN I tell it about xdot and xold so it can
# compute the normalization.
#
        nout = nsoli(FFUN, x0, FS, FPS; pdata = fdata, eta=.01, 
                     Krylov_Data=KDS, fixedeta=true, atol=1.e-8)
#
#  Stop the continuation if nsoli fails. This will usually be
#  because the change in x becomes too much for the predictor to track.
#
        idid = nout.idid
#
#  Record the solution and compute the derivative in lambda. 
#  xdot needs to go to the pseudo-arclength computation.  
#
        x = nout.solution
        xdot .= xold
        axpby!(dlamm1, x, -dlamm1, xdot)
#
#  Linear predictor
#
        x0 .= xold
        axpby!(2.0, x, -1.0, x0)
        xold .= x;
#
# update the data for the diagram
#
        bif_update(pval, nval, x, lambda)
#
# and continue to continue
#
        lambdaz=lambda
        lambda = lambda + dlam
        #
        if abs(lambda - lambdamax) < 1.e-12
            lambda = lambdamax
        end
        #
    end
    #
    return (pval, nval, x, lambdaz)
end
