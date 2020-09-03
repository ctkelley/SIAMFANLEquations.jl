"""
bvp_test()

Test nsold on the boundary value problem from Chapter 2 with
the LAPACK band solver.

Compare two small grids.
"""
function bvp_test(nsmall=101)
#
smallout=bvp_solve(nsmall; bfact=qr!)
smallout2=bvp_solve(nsmall)
hsmall=20.0/(nsmall-1)
statss=smallout.bvpout.stats
hs=smallout.bvpout.history./sqrt(hsmall)
hs2=smallout2.bvpout.history./sqrt(hsmall)
smok = norm(hs-hs2,Inf)<1.e-13
#
nbig=2*nsmall
bigout=bvp_solve(nbig;bfact=qr!)
hbig=20.0/(nbig-1)
statsb=bigout.bvpout.stats
bs=bigout.bvpout.history./sqrt(hbig)
#
armok=norm(statss.iarm-statsb.iarm) == 0
outok = norm(hs-bs,Inf) < .05
bvpok = outok && armok && smok
end

function bvp_solve(n = 801, T = Float64; bfact=qr)
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
    if bfact == qr
# Test to see if the default works.
    bvpout = nsold(Fbvp!, U0, FV, FPV, Jbvp!; rtol = 1.e-10,
             pdata = bdata)
    else
    bvpout = nsold(Fbvp!, U0, FV, FPV, Jbvp!; rtol = 1.e-10,
             pdata = bdata, jfact=bfact)
    end
    return (bvpout=bvpout, tv=tv)
end
