"""
kl_gmres(x0, b, atv, V, eta, ptv=nothing; kl_store=zeros(1,1); 
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

kl_store: You may at some point have the option of giving me some room
         for the vectors gmres needs to store copies of x0 and b,
         which I will not overwrite and a couple of vectors I use
         in the iteration. If you're only doing a linear solve, it
         does no harm to let me allocate those vectores in kl_gmres.
         If the solver is inside a loop, you should allocate this
         storage. nsoli and ptscoli allocate this without your having
         to do anything. Right now I'm not clear on an efficient way
         to do this.

pdata: precomputed data. The default is nothing, but that ain't gonna
        work well for nonlinear equations.

orth: your choice of the wise default, classical Gram-Schmidt twice,
       or something slower and less stable. Those are classical once (really
       bad) or a couple variants of modified Gram-Schmidt. mgs2 is what I
       used in my old matlab codes. Not terrible, but far from great.

side: left or right preconditioning. The default is "right".

lmaxit: maximum number of linear iterations. The default is -1, which
        means that the maximum number of linear iterations is K-1, which
        is all V will allow without restarts. If lmaxit > K-1, then the
        iteration will restart until you consume lmaxit iterations or
        terminate successfully.

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
              
# Examples: In these examples you have the matrix and use 
```
function atv(x, A)
    return A * x
end
to compute the matvec.
```

#### Three dimensional problem. Will converge in the correct three iterations
only if you orthogonalize with CGS twice. 

```jldoctest
julia> function atv(x, A)
           return A * x
       end
atv (generic function with 1 method)

julia> A = [0.001 0 0; 0 0.0011 0; 0 0 1.e4];

julia> V = zeros(3, 10); b = [1.0; 1.0; 1.0]; x0 = zeros(3);

julia> gout = kl_gmres(x0, b, atv, V, 1.e-10; pdata = A);

julia> gout.reshist
4-element Array{Float64,1}:
 1.73205e+00
 1.41421e+00
 6.72673e-02
 1.97712e-34

julia> norm(b - A*gout.sol,Inf)
1.28536e-10
```

#### Integral equation. Notice that pdata has the kernel of the 
operator and we do the matvec directly. Just like the previous example.
We put the grid information and, for this artifical example, the solution
in the precoputed data.

```jldoctest
julia> function integop(u, pdata)
           K = pdata.K
           return u - K * u
       end
integop (generic function with 1 method)

julia> function integopinit(n)
           h = 1 / n
           X = collect(0.5*h:h:1.0-0.5*h)
           K = [ker(x, y) for x in X, y in X]
           K .*= h
           sol = [usol(x) for x in X]
           f = sol - K * sol
           pdata = (K = K, xe = sol, f = f)
           return pdata
       end
integopinit (generic function with 1 method)

julia> function usol(x)
           return exp.(x) .* log.(2.0 * x .+ 1.0)
       end
usol (generic function with 1 method)

julia> function ker(x, y)
           ker = 0.1 * sin(x + exp(y))
       end
ker (generic function with 1 method)

julia> n=100; pdata = integopinit(n); ue = pdata.xe; f=pdata.f;

julia> u0 = zeros(size(f)); V = zeros(n, 20); V32=zeros(Float32,n,20);

julia> gout = kl_gmres(u0, f, integop, V, 1.e-10; pdata = pdata);

julia> gout32 = kl_gmres(u0, f, integop, V32, 1.e-10; pdata = pdata);

julia> [norm(gout.sol-ue,Inf) norm(gout32.sol-ue,Inf)]
1×2 Array{Float64,2}:
 4.44089e-16  2.93700e-07

julia> [gout.reshist gout32.reshist]
4×2 Array{Float64,2}:
 1.48252e+01  1.48252e+01
 5.52337e-01  5.52337e-01
 1.77741e-03  1.77742e-03
 1.29876e-19  8.73568e-11
```
 
"""
function kl_gmres(
    x0,
    b,
    atv,
    V,
    eta,
    ptv = nothing;
    kl_store = nothing,
    orth = "cgs2",
    side = "right",
    lmaxit = -1,
    pdata = nothing,
)
    
    # Build some precomputed data to inform KL_atv about 
    # preconditioning ...
    # Do not overwrite the initial iterate or the right hand side.
    n=length(x0)
    if kl_store != nothing
    y0=kl_store[1]
    y0.=x0
    rhs=kl_store[2]
    rhs .= b
    linsol=kl_store[3]
    restmp=kl_store[4]
    else
    y0=copy(x0)
    rhs=copy(b)
    # gmres_base needs two vectors
    linsol=zeros(size(b))
    restmp=zeros(size(b))
    end
    #
    if side == "right" || ptv == nothing
        itsleft = false
    else
        itsleft = true
        rhs .= ptv(rhs, pdata)
    end
    (n, K) = size(V)
    K > 1 || error("Must allocate for GMRES iterations. V must have 
                   at least two columns")
    klmaxit = lmaxit
    lmaxit > 0 || (lmaxit = K - 1)
    #
    itvec = maxitvec(K, lmaxit)
    ip=1
    idid=false
    Kpdata = (pdata = pdata, side = side, ptv = ptv, 
              atv = atv, linsol=linsol, restmp=restmp)
    gout=[]
#
# Restarted GMRES loop. 
#
    while ip <= length(itvec) && idid==false
    localout = gmres_base(y0, rhs, Katv, V, eta, Kpdata; 
                         lmaxit=itvec[ip], orth = orth)
    idid=localout.idid
    gout=outup(gout, localout, ip, klmaxit)
    reslen=length(localout.reshist)
#
# Update the termination criterion for the restart.
# gmres_base overwrites y0 with the solution
#
    idid || (eta = eta*localout.rho0/localout.reshist[reslen])
    ip += 1
    end
    #
    # Fixup the solution if preconditioning from the right.
    #
    sol=y0
    if side == "left" || ptv == nothing
        return (sol = sol, reshist = gout.reshist, lits = gout.lits, 
              idid = gout.idid)
    else
        sol .= ptv(sol, pdata)
        return (sol = sol, reshist = gout.reshist, lits = gout.lits, 
               idid = gout.idid)
    end
end

"""
Katv(x,Kpdata)

Builds a matrix-vector product to hand to gmres_base. Puts the preconditioner
in there on the correct side.
"""
function Katv(x, Kpdata)
#    y=copy(x)
    y=Kpdata.linsol
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
        y .= atv(x, pdata)
        return y
    elseif side == "left"
        y .= atv(x, pdata)
        return ptv(y, pdata)
    elseif side == "right"
        y .= ptv(x, pdata)
        return atv(y, pdata)
    end
end

"""
gmres_base(x0, b, atv, V, eta, pdata; orth="cgs2", lmaxit=-1)

Base GMRES solver. This is GMRES(m) with no restarts and no preconditioning.
The idea for the future is that it'll be called by kl_gmres (linear
solver) which is the backend of klgmres.

gmres_base overwrites x0 with the solution. This is one of many reasons
that you should not invoke it directly.
"""
function gmres_base(x0, b, atv, V, eta, pdata; orth = "cgs2", lmaxit=-1)
          
    (n, m) = size(V)
    #
    # Allocate for Givens
    #
    #    kmax = m - 1
    kmax = m 
    lmaxit == -1 || (kmax=lmaxit)
    kmax > m - 1 && error("lmaxit error in gmres_base")
    r = pdata.restmp
    r .= b
    T = eltype(V)
    h = zeros(T, kmax + 1, kmax + 1)
    c = zeros(kmax + 1)
    s = zeros(kmax + 1)
    #
    # Don't do the mat-vec if the intial iterate is zero
    #
    y=pdata.linsol
    (norm(x0) == 0.0) || (r .-= atv(x0, pdata))
#    (norm(x0) == 0.0) || (y .= atv(x0, pdata); r .-=y;)
    #
    #
    rho0 = norm(r)
    rho=rho0
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
    @views V[:, 1] .= r / rho
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
#    sol = x0
#    mul!(sol, qmf, y, 1.0, 1.0)
    mul!(x0, qmf, y, 1.0, 1.0)
    (rho <= errtol) || (idid = false)
    k > 0 || println("GMRES iteration terminates on entry.")
    return (rho0=rho0, reshist = Float64.(reshist), 
        lits = k, idid = idid)
end

function giveapp!(c, s, vin, k)
    for i = 1:k
        w1 = c[i] * vin[i] - s[i] * vin[i+1]
        w2 = s[i] * vin[i] + c[i] * vin[i+1]
        vin[i:i+1] .= [w1, w2]
    end
    return vin
end

#
# The functions maxitvec and outup manage the restarts.
# There is no reason to look at them or fiddle with them.
#

function maxitvec(K, lmaxit)
levels=Int.(ceil(lmaxit/(K-1)))
itvec=ones(Int,levels);
itvec[1:levels-1] .= K-1;
remainder=lmaxit-(levels-1)*(K-1) ;
itvec[levels]=remainder
return itvec
end

function outup(gout, localout, ip, klmaxit)
idid = localout.idid
#
# If I'm doing restarts I won't store the last residual
# unless the iteration is successful. The reason is that
# I will add that residual to the list when I restart.
#
if idid || klmaxit==-1
    lreshist=localout.reshist
else
   lk=length(localout.reshist)
   lreshist=localout.reshist[1:lk-1]
end
if ip == 1
  reshist=lreshist
  lits = localout.lits
else
   reshist = gout.reshist 
   append!(reshist,lreshist)
   lits = gout.lits + localout.lits
end
   gout = (reshist = reshist, lits = lits, 
                   idid = idid)
return gout
end

