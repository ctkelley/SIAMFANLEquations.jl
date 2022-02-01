"""
aasol(GFix!, x0, m, Vstore; maxit=20,
      rtol=1.e-10, atol=1.e-10, beta=1.0, pdata=nothing, keepsolhist = false)

C. T. Kelley, 2021

Julia code for Anderson acceleration. Nothing fancy.

Solvers fixed point problems x = G(x).

You must allocate storage for the function and fixed point map
history --> in the calling program <-- in the array Vstore.

For an n dimensional problem with Anderson(m), Vstore must have
at least 2m + 4 columns and 3m + 3 is better.  If m=0 (Picard) then
V must have at least 4 columns.

Inputs:\n

- GFix!: fixed-point map, the ! indicates that GFix! overwrites xout, your
    preallocated storage for the function value xout=G(xin).\n
    So xout=GFix!(xout,xin) or xout=GFix!(xout,xin,pdata) returns
    xout=G(xin)

- x0: Initial iterate. It is a vector of size N\n
  You should store it as (N,) and design G! to use vectors of size (N,).
  If you use (N,1) consistently instead, the solvers may work, but I make
  no guarantees.

- m: depth for Anderson acceleration. m=0 is Picard iteration

- Vstore: Working storage array. For an n dimensional problem Vstore
  should have at least 3m+3 columns unless you are storage bound. If storage
  is a problem, then you can allocate a minimum of 2m+4 columns. The smaller
  allocation exacts a performance penalty, especially for small problems
  and small values of m. So for Anderson(3), Vstore should be no smaller 
  than zeros(N,8) with zeros(N,11) a better choice. Vstore needs to
  allocate for the history of differences of the residuals and fixed
  point maps. The extra m-1 columns are for storing intermediate results
  in the downdating phase of the QR factorization for the coefficeint 
  matrix of the optimization problem. See the notebook or the print book 
  for the details of this mess. 

  If m=0, then Vstore needs 4 columns.

Keyword Arguments (kwargs):\n

maxit: default = 20\n
limit on nonlinear iterations\n

rtol and atol: default = 1.e-10\n
relative and absolute error tolerances\n

beta:\n
Anderson mixing parameter. Changes G(x) to (1-beta)x + beta G(x).
Equivalent to accelerating damped Picard iteration. The history
vector is the one for the damped fixed point map, not the original
one. Keep this in mind when comparing results.

pdata:\n
precomputed data for the fixed point map.
Things will go better if you use this rather than hide the data
in global variables within the module for your function.

keepsolhist: default = false\n
Set this to true to get the history of the iteration in the output
tuple. This is on by default for scalar equations and off for systems.
Only turn it on if you have use for the data, which can get REALLY LARGE.

Output:\n
- A named tuple (solution, functionval, history, stats, idid,
               errcode, solhist)
where

   -- solution = converged result

   -- functionval = G(solution)
      You might want to use functionval as your solution since it's
      a Picard iteration applied to the converged Anderson result. If G
      is a contraction it will be better than the solution.

   -- history = the vector of residual norms (||x-G(x)||) for the iteration

   -- stats = named tuple (condhist, alphanorm) of the history of the
              condition numbers of the optimization problem
              and l1 norm of the coefficients. 
This is only for diagosing
problems and research. Condihist[k] and alphanorm[k] are
the condition number and coefficient norm for the optimization
problem that computes iteration k+1 from iteration k. 

I record this for iterations k=1, ... until the final iteration 
K. So I do not record the stats for k=0 or the final iteration. 
We did record the data for the final iteration in Toth/Kelley 
2015 at the cost of an extra optimization problem solve. 
Since we've already terminated, there's not any point in 
collecting that data.\n
Bottom line: if history has length K+1 for iterations 
0 ... K, then condhist and alphanorm have length K-1.
 
   -- idid=true if the iteration succeeded and false if not.

   -- errcode = 0 if if the iteration succeeded

        = -1 if the initial iterate satisfies the termination criteria

        = -2 if || residual || > div_test || residual_0 ||
             I have fixed div_test = 1.e4 for now. I terminate 
             the iteration when this happens to avoid generating 
             Infs and/or NaNs.

        = 10 if no convergence after maxit iterations

   -- solhist:\n
This is the entire history of the iteration if you've set
keepsolhist=true\n

solhist is an N x K array where N is the length of x and K is the number of
iterations + 1. 

### Examples for aasol

#### Duplicate Table 1 from Toth-Kelley 2015.


The final entries in the condition number and coefficient norm statistics
are never used in the computation and we don't compute them in Julia.
See the docstrings, notebook, and the print book for the story on this.

```jldoctest
julia> function tothk!(G, u)
       G[1]=cos(.5*(u[1]+u[2]))
       G[2]=G[1]+ 1.e-8 * sin(u[1]*u[1])
       return G
       end
tothk! (generic function with 1 method)

julia> u0=ones(2,); m=2; vdim=3*m+3; Vstore = zeros(2, vdim);
julia> aout = aasol(tothk!, u0, m, Vstore; rtol = 1.e-10);
julia> aout.history
8-element Vector{Float64}:
 6.50111e-01
 4.48661e-01
 2.61480e-02
 7.25389e-02
 1.53107e-04
 1.18513e-05
 1.82466e-08
 1.04725e-13

julia> [aout.stats.condhist aout.stats.alphanorm]
6×2 Matrix{Float64}:
 1.00000e+00  1.00000e+00
 2.01556e+10  4.61720e+00
 1.37776e+09  2.15749e+00
 3.61348e+10  1.18377e+00
 2.54948e+11  1.00000e+00
 3.67694e+10  1.00171e+00
```

Now we put a mixing or damping parameter in there with beta = .5. This
example is nasty enough to make mixing do ok. Keep in mind
that the history is for the damped residual, not the original one.

```
julia> bout=aasol(tothk!, u0, m, Vstore; rtol = 1.e-10, beta=.5);

julia> bout.history
7-element Vector{Float64}:
 3.25055e-01
 3.70140e-02
 1.81111e-03
 9.55308e-04
 1.25936e-05
 1.40854e-09
 2.18196e-12
```


#### H-equation example with m=2. This takes more iterations than
Newton, which should surprise no one.

```jldoctest
julia> n=16; x0=ones(n,); Vstore=zeros(n,20); m=2;
julia> hdata=heqinit(x0,.99);
julia> hout=aasol(HeqFix!, x0, m, Vstore; pdata=hdata);
julia> hout.history
12-element Vector{Float64}:
 1.47613e+00
 7.47800e-01
 2.16609e-01
 4.32017e-02
 2.66867e-02
 6.82965e-03
 2.70779e-04
 6.51027e-05
 7.35581e-07
 1.85649e-09
 4.94803e-10
 5.18866e-12

```
"""
function aasol(
    GFix!,
    x0,
    m,
    Vstore;
    maxit = 20,
    rtol = 1.e-10,
    atol = 1.e-10,
    beta = 1.0,
    pdata = nothing,
    keepsolhist = false,
)
    #
    # Startup
    #
    # Set up the storage
    #
    (sol, gx, df, dg, res, DG, QP, Qd, solhist) =
           Anderson_Init(x0, Vstore, m, maxit, beta, keepsolhist)
    #
    #   Iteration 1
    #
    k = 0
    ~keepsolhist || (@views solhist[:, k+1] .= sol)
    gx = EvalF!(GFix!, gx, sol, pdata)
    (beta == 1.0) || (gx=betafix!(gx, sol, beta))
    copy!(res,gx)
    axpy!(-1.0,sol,res)
#    res .= gx - sol
    resnorm = norm(res)
    resnorm_up_bd = 1.e4*resnorm
    tol = rtol * resnorm + atol
    ItData = ItStatsA(resnorm)
    toosoon = (resnorm <= tol)
    if ~toosoon
    #
    #   If we need more iterations, get organized.
    #
#    sol .= gx
    copy!(sol, gx)
    alpha = zeros(m + 1)
    k = k + 1
    ~keepsolhist || (@views solhist[:, k+1] .= sol)
    (gx, dg, df, res, resnorm) =
        aa_point!(gx, GFix!, sol, res, dg, df, beta, pdata)
    updateHist!(ItData, resnorm)
    end
    n=length(x0)
    RF=zeros(m,m)
    RP=zeros(m,m)
    ThetA=zeros(m)
    TmPReS=zeros(m)
    while ((k < maxit) && (resnorm > tol) && ~toosoon 
                      && (resnorm < resnorm_up_bd))
        if m == 0
            alphanrm = 1.0
            condit = 1.0
            copy!(sol, gx)
 #           sol .= gx
        else
            BuildDG!(DG,m,k+1,dg)
            (QP, RP) = aa_qr_update!(QP, RP, df, m, k-1, Qd)
            mk = min(m, k)
            @views QA=QP[:,1:mk]
            @views RA=RP[1:mk,1:mk]
            @views theta=ThetA[1:mk]
            @views tres=TmPReS[1:mk]
            mul!(tres,QA',res)
            theta .= RA\tres
            condit = cond(RA)
            alphanrm = falpha(alpha, theta, min(m, k))
            copy!(sol, gx)
#            @views sol .-= DG[:, 1:mk] * theta
            @views mul!(sol, DG[:,1:mk], theta, -1.0, 1.0)
        end
        updateStats!(ItData, condit, alphanrm)
        k += 1
        ~keepsolhist || (@views solhist[:, k+1] .= sol)
        (gx, dg, df, res, resnorm) =
            aa_point!(gx, GFix!, sol, res, dg, df, beta, pdata)
        updateHist!(ItData, resnorm)
    end
    (idid, errcode) = AndersonOK(resnorm, tol, k, m, toosoon, resnorm_up_bd)
    aaout=CloseIteration(sol, gx, ItData, idid, errcode, keepsolhist, solhist)
    return aaout
end

"""
BuildDG!(DG,m,k,dg)

Keeps the history of the fixed point map differences
"""
function BuildDG!(DG,m,k,dg)
    if m == 1
        @views copy!(DG[:, 1], dg)
    elseif k > m + 1
        for ic = 1:m-1
#            @views DG[:, ic] .= DG[:, ic+1]
            @views copy!(DG[:,ic], DG[:, ic+1])
        end
        @views copy!(DG[:, m], dg)
    else
        @views copy!(DG[:, k-1], dg)
    end
end

"""
aa_point!(gx, gfix, sol, res, dg, df, pdata)

Evaluate the fixed point map at the new point. 
Keep the books to get ready to update the coefficient matrix
for the optimization problem.
"""
function aa_point!(gx, gfix, sol, res, dg, df, beta, pdata)
#    dg .= -gx
    copy!(dg, -gx)
    gx = EvalF!(gfix, gx, sol, pdata)
    (beta == 1.0) || (gx=betafix!(gx, sol, beta))
#    dg .+= gx
    axpy!(1.0, gx, dg)
#    df .= -res
    copy!(df,-res)
#    res .= gx - sol
#    res .= gx
#    res .-= sol
    copy!(res,gx)
    axpy!(-1.0,sol,res)
    axpy!(1.0, res, df)
#    df .+= res
    resnorm = norm(res)
    return (gx, dg, df, res, resnorm)
end

"""
betafix(gx, sol, dg, beta)

Put the mixing parameter beta in the right place.
"""
function betafix!(gx, sol, beta)
gx=axpby!((1.0-beta),sol,beta,gx)
return gx
end

"""
aa_qr_update(Q, R, vnew, m, k, Qd)

Update the QR factorization for the Anderson acceleration optimization
problem.

Still need to make the allocation for Qd go away.
"""
function aa_qr_update!(Q, R, vnew, m, k, Qd)
(n,m)=size(Q)
    aaqr_dim_check(Q, R, vnew, m, k)
    if k == 0
        R[1, 1] = norm(vnew)
        @views Q[:, 1] .= vnew / norm(vnew)
    else
        if k > m - 1
            downdate_aaqr!(Q, R, m, Qd)
        end # inner if block
        kq = min(k, m - 1)
        update_aaqr!(Q, R, vnew, m, kq)
    end # outer if block
    return (Q, R)
end

function update_aaqr!(Q, R, vnew, m, k)
    (nq, mq) = size(Q)
    (k > m - 1) && error("Dimension error in Anderson QR")
    @views Qkm=Q[:,1:k]
    @views hv = vec(R[1:k+1, k+1])
    Orthogonalize!(Qkm, hv, vnew, "cgs2")
    @views R[1:k+1, k+1] .= hv
    @views Q[:, k+1] .= vnew
#    return (Q = Q, R = R)
end

function downdate_aaqr!(Q, R, m, Qd)
(nq,mq)=size(Q)
(pd,md)=size(Qd)
(md == m-1) || @error("dimension error in downdate")
    @views Rp = R[:, 2:m]
    G = qr!(Rp)
    Rd = Matrix(G.R)
    Qx = Matrix(G.Q)
    @views R[1:m-1, 1:m-1] .= Rd
    @views R[:,m].=0.0
if (pd==nq)
    mul!(Qd,Q,Qx)
    @views Q[:,1:m-1] .= Qd
else
    blocksize=pd
    (dlow,dhigh)=blockdim(nq,blocksize)
    blen=length(dlow)
    for il=1:blen
    asize=dhigh[il]-dlow[il]+1
    @views QZ=Qd[1:asize,:]
    @views Qsec=Q[dlow[il]:dhigh[il],:]
    @views mul!(QZ,Qsec,Qx)
    @views Qsec[:,1:m-1] .= QZ
    end
end
    @views Q[:, m] .= 0.0
    return (Q, R)
end

function aaqr_dim_check(Q, R, vnew, m, k)
    (mq, nq) = size(Q)
    (mr, nr) = size(R)
    n = length(vnew)
    dimqok = ((mq == n) && (nq == m))
    dimrok = ((mr == m) && (nr == m))
    dimok = (dimqok && dimrok)
    dimok || error("array size error in AA update")
end

function blockdim(n, block)
    p = Int(floor(n / block))
    res = n - p * block
    ilow = Int64[]
    ihigh = Int64[]
    for jb = 1:p
        lowval = (jb - 1) * block + 1
        push!(ilow, lowval)
        highval = ilow[jb] + block - 1
        push!(ihigh, highval)
    end
    if res > 0
        lowval = p * block + 1
        push!(ilow, lowval)
        push!(ihigh, n)
    end
    return (ilow, ihigh)
end
