"""
kl_gmres(x0, b, atv, V, eta, ptv=nothing; 
            orth = "cgs2", side="right", lmaxit=-1, pdata=nothing)

Gmres linear solver. Handles preconditioning and (coming soon) restarts. 
Uses gmres_base which is completely oblivious to these things.

The deal is

Input:\n
x0:  initial iterate, this is usually zero for nonlinear solvers

b: right hand side (duh!)

atv:  matrix-vector product which depends on precomputed data pdta
      I expect you to use pdata most or all of the time, so it is not
      an optional argument, even if it's nothing (at least for now). 
      If your mat-vec is just A*v, you have to write a function where 
      A is the precomputed data.
      API for atv is av=atv(v,pdata)

V:  Preallocated n x K array for the Krylov vectors. I store the initial
    normalized residual in column 1, so  you have at most K-1 iterations
    before gmres\\_base returns a failure. kl\\_gmres will handle the restarts.

eta: Termination happens when ||b - Ax|| <= eta || b ||

ptv:  preconditioner-vector product, which will also use pdata. The
      default is nothing, which is no preconditioning at all.
      API for ptv is px=ptv(x,pdata)

Keyword arguments

pdata: precomputed data. The default is nothing, but that ain't gonna
        work well for nonlinear equations.

orth: your choice of the wise default, classical Gram-Schmidt twice,
       or something slower and less stable. Those are classical once (really
       bad) or a couple variants of modified Gram-Schmidt. mgs2 is what I
       used in my old matlab codes. Not terrible, but far from great.

side: left or right preconditioning. The default is "right".

lmaxit: maximum number of linear iterations. The default is -1, which
        means that the maximum number of linear iterations is K-1, which
        is all V will allow without restarts.

Other parameters on the way.

Output:\n
A named tuple (sol, reshist, lits, idid)

where

sol= final result
reshist = residual norm history
lits = number of iterations
idid = status of the iteration
       true -> converged 
       false -> failed to converge
              
 
"""
function kl_gmres(
    x0,
    b,
    atv,
    V,
    eta,
    ptv = nothing;
    orth = "cgs2",
    side = "right",
    lmaxit = -1,
    pdata = nothing,
)
    #
    # Build some precomputed data to inform KL_atv about preconditioning ...
    #
    if side == "right" || ptv == nothing
        rhs = b
    else
        rhs = ptv(b, pdata)
    end
    (n,K) = size(V)
    lmaxit > 0 || (lmaxit=K-1)
#
#   Temporary fix before I get restarts in there
#
    lmaxit > K-1 && (lmaxit=K-1)
    Kpdata = (pdata = pdata, side = side, ptv = ptv, atv = atv, lmaxit=lmaxit)
    gout = gmres_base(x0, rhs, Katv, V, eta, Kpdata; orth = orth)
    #
    # Fixup the solution if preconditioning from the right.
    #
    if side == "left" || ptv == nothing
        return gout
    else
        sol = ptv(gout.sol, pdata)
        return (sol = sol, reshist = gout.reshist, lits = gout.lits, idid = gout.idid)
    end
end

"""
Katv(x,Kpdata)

Builds a matrix-vector product to hand to gmres_base. Puts the preconditioner
in there on the correct side.
"""
function Katv(x, Kpdata)
    pdata = Kpdata.pdata
    ptv = Kpdata.ptv
    atv = Kpdata.atv
    side = Kpdata.side
    sideok = (side == "left") || (side == "right")
    sideok || error(
        "Bad preconditioner side in kl_gmres, input side = ",
        side,
        ". Side must be \"left\" or \"right\" ",
    )
    if ptv == nothing
        y = atv(x, pdata)
        return y
    elseif side == "left"
        y = atv(x, pdata)
        return ptv(y, pdata)
    elseif side == "right"
        y = ptv(x, pdata)
        return atv(y, pdata)
    end
end

"""
gmres_base(x0, b, atv, V, eta, pdata; orth="mgs1")

Base GMRES solver. This is GMRES(m) with no restarts and no preconditioning.
The idea for the future is that it'll be called by kl_gmres (linear
solver) which
is the backend of klgmres.
"""
function gmres_base(x0, b, atv, V, eta, pdata; orth = "mgs1")
    (n, m) = size(V)
    #
    # Allocate for Givens
    #
#    kmax = m - 1
    kmax=pdata.lmaxit
    kmax > m-1 && error("lmaxit error in gmres_base")
    r = copy(b)
    T = eltype(V)
    h = zeros(T, kmax + 1, kmax + 1)
    c = zeros(kmax + 1)
    s = zeros(kmax + 1)
    #
    # Don't do the mat-vec if the intial iterate is zero
    #
    (norm(x0) == 0.0) || (r .-= atv(x0, pdata))
    #
    #
    rho = norm(r)
    #
    # Initial residual = 0? This can't be good.
    #
    rho == 0.0 && error("Initial resdiual in kl_gmres is zero. Why?")
    #
    g = zeros(size(c))
    g[1] = rho
    errtol = eta * rho
    reshist = []
    #
    # Initialize
    #
    idid = true
    push!(reshist, rho)
    k = 0
    #
    # Showtime!
    #
    @views V[:, 1] = r / rho
    beta = rho
    while (rho > errtol) && (k < kmax)
        k += 1
        @views V[:, k+1] .= atv(V[:, k], pdata)
        @views vv = vec(V[:, k+1])
        @views hv = vec(h[1:k+1, k])
        @views Vkm = V[:, 1:k]
        #
        # Don't mourn. Orthogonalize!
        #
        Orthogonalize!(Vkm, hv, vv, orth)
        #
        # Build information for new Givens rotations.
        #   
        if k > 1
            hv = @view h[1:k, k]
            giveapp!(c[1:k-1], s[1:k-1], hv, k - 1)
        end
        nu = norm(h[k:k+1, k])
        if nu != 0
            c[k] = conj(h[k, k] / nu)
            s[k] = -h[k+1, k] / nu
            h[k, k] = c[k] * h[k, k] - s[k] * h[k+1, k]
            h[k+1, k] = 0.0
            gv = @view g[k:k+1]
            giveapp!(c[k], s[k], gv, 1)
        end
        #
        # Update the residual norm.
        #
        rho = abs(g[k+1])
        push!(reshist, rho)
    end
    #
    # At this point either k = kmax or rho < errtol.
    # It's time to compute x and check out.
    #
    y = h[1:k, 1:k] \ g[1:k]
    qmf = view(V, 1:n, 1:k)
    #    mul!(r, qmf, y)
    #    r .= qmf*y    
    #    x .+= r
    r .= x0
    mul!(r, qmf, y, 1.0, 1.0)
    (rho <= errtol) || (idid = false)
    k > 0 || println("GMRES iteration terminates on entry.")
    return (sol = r, reshist = Float64.(reshist), lits = k, idid = idid)
end

function giveapp!(c, s, vin, k)
    nv = length(vin)
    for i = 1:k
        w1 = c[i] * vin[i] - s[i] * vin[i+1]
        w2 = s[i] * vin[i] + c[i] * vin[i+1]
        vin[i:i+1] .= [w1, w2]
    end
    return vin
end
