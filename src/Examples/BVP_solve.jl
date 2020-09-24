"""
BVP_solve(n = 801, T = Float64; bfact=qr)
Solve the BVP for the Chapter 2 figures and testing.
"""
function BVP_solve(n = 801, T = Float64; bfact=qr)
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
# Test to see if the default (qr instead of qr!) works.
    bvpout = nsol(Fbvp!, U0, FV, FPV, Jbvp!; rtol = 1.e-10,
             pdata = bdata)
    else
    bvpout = nsol(Fbvp!, U0, FV, FPV, Jbvp!; rtol = 1.e-10,
             pdata = bdata, jfact=bfact)
    end
    return (bvpout=bvpout, tv=tv)
end
