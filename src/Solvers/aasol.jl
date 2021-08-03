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
              and l1 norm of the coefficients. This is only for diagosing
              problems and research.
 
   -- idid=true if the iteration succeeded and false if not.

   -- errcode = 0 if if the iteration succeeded

        = -1 if the initial iterate satisfies the termination criteria

        = 10 if no convergence after maxit iterations

   -- solhist:\n
This is the entire history of the iteration if you've set
keepsolhist=true\n

solhist is an N x K array where N is the length of x and K is the number of
iterations + 1. 
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
            Abuild4!(DG, DF, m, k + 1, dg, df)
            mk = min(m, k)
            @views AC = DF[:, 1:mk]
            @views theta = AC \ res
            condit = cond(AC)
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
Abuild4!(DG, DF, m, k, dg, df)

Builds the coefficient matrix for the optimization problem. 
This will get replaced by something that up/down dates the QR
factorization before v0.4.3 is released.
"""
function Abuild4!(DG, DF, m, k, dg, df)
    if m == 1
        @views copy!(DF[:, 1], df)
        @views copy!(DG[:, 1], dg)
    elseif k > m + 1
        for ic = 1:m-1
            @views DF[:, ic] .= DF[:, ic+1]
            @views DG[:, ic] .= DG[:, ic+1]
        end
        @views copy!(DF[:, m], df)
        @views copy!(DG[:, m], dg)
    else
        @views copy!(DF[:, k-1], df)
        @views copy!(DG[:, k-1], dg)
    end
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
