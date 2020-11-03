"""
Fbvp! and Jbvp! are the function and Jacobian evaluations for the
boundary value problem example in Chapter 2.

The Jacobian is banded and I've padded the storage so I can use lu! or qr!
for the linear solver. You need to be careful about this and RTFM (ie the
LAPACK or LINPACK) manual to get the padding right. The short story is
that if your upper/lower bandwidths are lu/ll, then you must store in 
in a matrix with bandwidts lu+2/ll with zeros in the unused bands.

That aside, there is not much here that I did not explain in the book.
"""

function Fbvp!(FV, U, bdata)
    n2 = length(U)
    n = bdata.n
    n2 == 2n || error("dimension error in Fbvp")
    force = bdata.force
    tv = bdata.tv
    tvdag = bdata.tvdag
    h = bdata.h
    FV[1] = U[2]
    FV[2n] = U[2n-1]
    v = view(U,1:2:2n-1)
    vp= view(U,2:2:2n)
    force .= Phi.(tv, tvdag, vp, v)
    h2 = 0.5 * h
    @inbounds @simd for ip = 1:n-1
        FV[2*ip+1] = v[ip+1] - v[ip] - h2 * (vp[ip] + vp[ip+1])
        FV[2*ip] = vp[ip+1] - vp[ip] + h2 * (force[ip] + force[ip+1])
    end
end

function Jbvp!(FVP, FV, x, bdata)
    n = bdata.n
    tv = bdata.tv
    tvdag = bdata.tvdag
    h = bdata.h
    h2 = h * 0.5
    zdat = bdata.zdat
    DiagFP = bdata.DiagFP
    jacinit!(FVP, DiagFP)
    #
    # Using @view to avoid allocations. Build the vector I'll
    # need to populate the Jacobian.
    #
#    @views zdat[1:n] .= (h .* tv[1:n] .* x[1:2:2n-1] .- h2)
    @views zdat .= (h .* tv[1:n] .* x[1:2:2n-1] .- h2)
    #
    # The diagnd function gets the the diagonals so I can populate
    # them without allocations.
    #
    nup = diagind(FVP, 1)
    FUP = view(FVP, nup)
    @views FUP[2:2:2n-2] .= zdat[2:n]
    ndown = diagind(FVP, -1)
    FDOWN = view(FVP, ndown)
    @views FDOWN[1:2:2n-3] .= zdat[1:n-1]
end


function Phi(t, tdag, vp, v)
    phi = 4.0 * tdag * vp + (t * v - 1.0) * v
    return phi
end

function bvpinit(n, T = Float64)
#
# Allocate space for the Jacobian, compute the parts of the Jacobian
# that don't depend on the iteration, and store a few vectors.
#
    h = 20.0 / (n - 1)
    h2 = h * 0.5
    tv = collect(0:h:20.0)
    tvdag = collect(0:h:20.0)
    @views tvdag[2:n] .= (1.0 ./ tv[2:n])
    force = zeros(n)
    D = ones(T,2n)
    D[1] = 0.0
    D[2n] = 0.0
    h4=4*h2
    @views D[2:2:2n-2] .= (-1 .+ h4 .* tvdag[1:n-1])
    D1 = zeros(T,2n - 1)
    D1[1] = 1.0
#    @views D1[3:2:2n-1] .= -h2
    view(D1,3:2:2n-1) .= -h2
    Dm1 = zeros(T,2n - 1)
    view(Dm1,2:2:2n-2).= -h2
    Dm1[2n-1] = 1.0
#    @views Dm1[2:2:2n-2] .= -h2
    Dm2 = zeros(T,2n - 2)
    view(Dm2,1:2:2n-3) .= -1.0
#    @views Dm2[1:2:2n-3] .= -1.0
    D2 = zeros(T,2n - 2)
    @views D2[2:2:2n-2] .= (1.0 .+  h4 .* tvdag[2:n])
#
# The bandwidths are lu=ll=2, so my padded matrix gets lu=4.
# Allocate the storage and precompute the bands that don't change. 
#
#    DiagFP = (Dm2 = Dm2, Dm1 = Dm1, D = D, D1 = D1, D2 = D2)
    DiagFP = [Dm2, Dm1, D, D1, D2]
    zdat = zeros(T,n)
    return ( h = h, tv = tv, force = force, tvdag = tvdag, zdat = zdat,
        n = n, DiagFP = DiagFP,)
end

function jacinit!(FVP, DiagFP)
#
# Fill the unused padding bands with zeros.
#
    for ip = 3:4
#        view(FVP,band(ip)) .= 0.0
        FVP[band(ip)] .= 0.0
    end
#
# Get the bands you've computed.
# Put the good bands in the right place.
#
for ip=-2:2
    ib=ip+3
    FVP[band(ip)] .= DiagFP[ib]
end
#     FVP[band(-2)] .= DiagFP.Dm2
#     FVP[band(-1)] .= DiagFP.Dm1
#     FVP[band(0)] .= DiagFP.D
#     FVP[band(1)] .= DiagFP.D1
#     FVP[band(2)] .= DiagFP.D2
#    view(FVP,band(-2)) .= DiagFP.Dm2
#    view(FVP,band(-1)) .= DiagFP.Dm1
#    view(FVP,band(0)) .= DiagFP.D
#    view(FVP,band(1)) .= DiagFP.D1
#    view(FVP,band(2)) .= DiagFP.D2
end
