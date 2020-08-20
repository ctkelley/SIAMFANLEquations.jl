"""
bvp_test()

Test nsold on the boundary value problem from Chapter 2 with
the LAPACK band solver.

Compare two small grids.
"""
function bvp_test(nsmall=101)
smallout=bvp_solve(nsmall)
hsmall=20.0/(nsmall-1)
nbig=2*nsmall
bigout=bvp_solve(nbig)
hbig=20.0/(nbig-1)
statss=smallout.bvpout.stats
statsb=bigout.bvpout.stats
armok=norm(statss.iarm-statsb.iarm) == 0
hs=smallout.bvpout.history./sqrt(hsmall)
bs=bigout.bvpout.history./sqrt(hbig)
outok = norm(hs-bs,Inf) < .05
bvpok = outok && armok
end

function bvp_solve(n = 801, T = Float64)
    bdata = bvpinit(n, T)
    U0 = zeros(2n)
    FV = zeros(2n)
    FPV = bdata.JacS
    tv = bdata.tv
    sv = -.1 * tv .* tv
    v = exp.(sv)
    vp = -.2 * tv .* v
    U0[1:2:2n-1] = v
    U0[2:2:2n] = vp
    bvpout = nsold(Fbvp!, U0, FV, FPV, Jbvp!; rtol = 1.e-10,
             pdata = bdata)
    return (bvpout=bvpout, tv=tv)
end


