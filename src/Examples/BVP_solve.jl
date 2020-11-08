"""
BVP_solve(n = 801, T = Float64; bfact=qr!)
Solve the BVP for the Chapter 2 figures and testing.
"""
function BVP_solve(n = 801, T = Float64; bfact = qr!)
    # Set it up
    bdata = bvpinit(n, T)
    #
    U0 = zeros(2n)
    FV = zeros(2n)
    # Banded matrix with the correct number of bands
    FPV = BandedMatrix{T}(Zeros(2n, 2n), (2, 4))
    #
    # Build the initial iterate
    #
    BVP_U0!(U0, n, bdata)
    #
    if bfact == qr!
        bvpout = nsol(
            Fbvp!,
            U0,
            FV,
            FPV,
            Jbvp!;
            rtol = 1.e-10,
            sham = 1,
            pdata = bdata,
            jfact = bfact,
        )
    else
        # Test for default of qr. Used for CI only.
        bvpout = nsol(Fbvp!, U0, FV, FPV, Jbvp!; rtol = 1.e-10, sham = 1, pdata = bdata)
    end
    return (bvpout = bvpout, tv = bdata.tv)
end

function BVP_U0!(U0, n, bdata)
    tv = bdata.tv
    view(U0, 1:2:2n-1) .= exp.(-.1 .* tv .* tv)
    view(U0, 2:2:2n) .= -.2 .* view(U0, 1:2:2n-1) .* tv
end
