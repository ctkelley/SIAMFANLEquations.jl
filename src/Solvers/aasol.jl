"""
aasol(GFix!, x0, m, Vstore; maxit=20,
      rtol=1.e-10, atol=1.e-10, beta=1.0, pdata=nothing, keepsolhist = false)

C. T. Kelley, 2021

Julia code for Anderson acceleration. Nothing fancy.

Solvers fixed point problems x = G(x).

WARNING! This is an alpha version of the solver. The API may change
and the management of the optimization problem WILL change.

You must allocate storage for the function and fixed point map
history --> in the calling program <-- in the array Vstore.

For an n dimensional problem with Anderson(m), Vstore must have
at least 2m+2 columns. THIS MAY CHANGE!!!

Inputs:\n

- GFix!: fixed-point map, the ! indicates that GFix! overwrites xout, your
    preallocated storage for the function value xout=G(xin).\n
    So xout=GFix!(xout,xin) or xout=GFix!(xout,xin,pdata) returns
    xout=G(xin)

- sol: Initial iterate

- m: depth for Anderson acceleration. m=0 is Picard iteration

- Vstore: Working storage array. For an n dimensional problem Vstore
  should have at least 2m+2 columns. So for Anderson(3), Vstore should
  be no smaller than zeros(n,8)

Keyword Arguments (kwargs):\n

maxit: default = 20\n
limit on nonlinear iterations\n

rtol and atol: default = 1.e-10\n
relative and absolute error tolerances\n

beta:\n
Anderson mixing parameter. Changes G(x) to (1-beta)x + beta G(x).
Equivalent to accelerating damped Picard iteration.

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
2015 at the cost of an extra optimiztion problem solve. 
Since we've already terminated, there's not any point in 
collecting that data.\n
Bottom line: if history has length K+1 for iterations 
0 ... K, then condhist and alphanorm have length K-1.
 
   -- idid=true if the iteration succeeded and false if not.

   -- errcode = 0 if if the iteration succeeded

        = -1 if the initial iterate satisfies the termination criteria

        = 10 if no convergence after maxit iterations

   -- solhist:\n
This is the entire history of the iteration if you've set
keepsolhist=true\n

solhist is an N x K array where N is the length of x and K is the number of
iterations + 1. 

### Examples from the docstrings for aasol

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

julia> u0=ones(2,); m=2; vdim=2*(m+1); Vstore = zeros(2, vdim);
julia> aout = aasol(tothk!, u0, m, Vstore; rtol = 1.e-10);
julia> aout.history
8-element Vector{Float64}:
 6.50111e-01
 4.48661e-01
 2.61480e-02
 7.25389e-02
 1.53107e-04
 1.18512e-05
 1.82476e-08
 1.04804e-13

julia> [aout.stats.condhist aout.stats.alphanorm]
6×2 Matrix{Float64}:
 1.00000e+00  1.00000e+00
 2.01556e+10  4.61720e+00
 1.37776e+09  2.15749e+00
 3.61344e+10  1.18377e+00
 2.54947e+11  1.00000e+00
 3.67672e+10  1.00171e+00
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
    (sol, DG, DF, solhist) =
           Anderson_Init(x0, Vstore, m, maxit, beta, keepsolhist)
    gx = @views Vstore[:,2*m+1]
    df = copy(sol)
    dg = copy(sol)
    gold = copy(sol)
    res = copy(sol)
    resold = copy(sol)
    #
    #   Iteration 1
    #
    k = 0
    ~keepsolhist || (@views solhist[:, k+1] .= sol)
    gx = EvalF!(GFix!, gx, sol, pdata)
    (beta == 1.0) || (gx=betafix!(gx, sol, beta))
    gold .= gx
    res .= gx - sol
    resnorm = norm(res)
    tol = rtol * resnorm + atol
    ItData = ItStatsA(resnorm)
    toosoon = (resnorm <= tol)
#    if toosoon
#        println("aasol terminates on entry")
#    else     
    if ~toosoon
    #
    #   If we need more iterations, get organized.
    #
    sol .= gx
    alpha = zeros(m + 1)
    k = k + 1
    ~keepsolhist || (@views solhist[:, k+1] .= sol)
    (gx, dg, df, res, resold, resnorm) =
        aa_point!(gx, GFix!, gold, sol, res, resold, dg, df, beta, pdata)
    updateHist!(ItData, resnorm)
    end
    while (k < maxit) && resnorm > tol && ~toosoon
        if m == 0
            alphanrm = 1.0
            condit = 1.0
            sol .= gx
        else
            BuildDG!(DG,m,k+1,dg)
            (QA, RA) = BuildDF!(DF,m,k+1,df)
            mk = min(m, k)
            theta = RA\ (QA'*res)
            condit = cond(RA)
            alphanrm = falpha(alpha, theta, min(m, k))
            copy!(sol, gx)
            @views sol .-= DG[:, 1:mk] * theta
        end
        updateStats!(ItData, condit, alphanrm)
        k += 1
        ~keepsolhist || (@views solhist[:, k+1] .= sol)
        (gx, dg, df, res, resold, resnorm) =
            aa_point!(gx, GFix!, gold, sol, res, resold, dg, df, beta, pdata)
        updateHist!(ItData, resnorm)
    end
    (idid, errcode) = AndersonOK(resnorm, tol, k, toosoon)
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
            @views DG[:, ic] .= DG[:, ic+1]
        end
        @views copy!(DG[:, m], dg)
    else
        @views copy!(DG[:, k-1], dg)
    end
end



"""
BuildDF!(DG, DF, m, k, dg, df)

Builds the coefficient matrix for the optimization problem. 
This will get replaced by something that up/down dates the QR
factorization before v0.4.3 is released.
"""
function BuildDF!(DF, m, k, df)
    if m == 1
        @views copy!(DF[:, 1], df)
    elseif k > m + 1
        for ic = 1:m-1
            @views DF[:, ic] .= DF[:, ic+1]
        end
        @views copy!(DF[:, m], df)
    else
        @views copy!(DF[:, k-1], df)
    end
mk=min(k-1,m)
@views AC=DF[:,1:mk]
QZ=qr(AC)
QA=Matrix(QZ.Q)
RA=Matrix(QZ.R)
return (QA, RA)
end

"""
aa_point!(gx, gfix, gold, sol, res, resold, dg, df, pdata)

Evaluate the fixed point map at the new point. 
Keep the books to get ready to update the coefficient matrix
for the optimization problem.
"""
function aa_point!(gx, gfix, gold, sol, res, resold, dg, df, beta, pdata)
    gold .= gx
    gx = EvalF!(gfix, gx, sol, pdata)
    (beta == 1.0) || (gx=betafix!(gx, sol, beta))
    dg .= gx - gold
    resold .= res
    res .= gx - sol
    df .= res - resold
    resnorm = norm(res)
    return (gx, dg, df, res, resold, resnorm)
end

"""
betafix(gx, sol, dg, beta)

Put the mixing parameter beta in the right place.
"""
function betafix!(gx, sol, beta)
gx=axpby!((1.0-beta),sol,beta,gx)
return gx
end
