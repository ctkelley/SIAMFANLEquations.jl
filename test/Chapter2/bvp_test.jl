"""
bvp_test()

Test nsold on the boundary value problem from Chapter 2 with
the LAPACK band solver.

Compare two small grids.
"""
function bvp_test(nsmall=101)
#
smallout=BVP_solve(nsmall; bfact=qr!)
smallout2=BVP_solve(nsmall)
hsmall=20.0/(nsmall-1)
statss=smallout.bvpout.stats
hs=smallout.bvpout.history./sqrt(hsmall)
hs2=smallout2.bvpout.history./sqrt(hsmall)
smok = norm(hs-hs2,Inf)<1.e-13
#
nbig=2*nsmall
bigout=BVP_solve(nbig;bfact=qr!)
hbig=20.0/(nbig-1)
statsb=bigout.bvpout.stats
bs=bigout.bvpout.history./sqrt(hbig)
#
armok=norm(statss.iarm-statsb.iarm) == 0
outok = norm(hs-bs,Inf) < .05
bvpok = outok && armok && smok
end
