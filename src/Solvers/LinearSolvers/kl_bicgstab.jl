"""
kl\\_bicgstab( x0, b, atv, V, eta, ptv = nothing; lmaxit = 10, 
      pdata = nothing, side = "right",)

BiCGSTAB linear solver. Deals with preconditioning. 
Uses bicgstab\\_base with is oblivious to that.

The code works and does what it needs to do, but ...\n
The user interface is unstable and, even worse,
the nonlinear solvers are not hooked up yet.

The way this works it 

Input:\n
x0:  initial iterate, this is usually zero for nonlinear solvers

b: right hand side (duh!)

atv:  matrix-vector product which depends on precomputed data pdta
      I expect you to use pdata most or all of the time, so it is not
      an optional argument, even if it's nothing (at least for now).
      If your mat-vec is just A*v, you have to write a function where
      A is the precomputed data.
      API for atv is av=atv(v,pdata)

V: a vector for me to store a Jacobian-vector product. It goes where 
   FPS would go in gmres

eta: Termination happens when ||b - Ax|| <= eta || b ||

ptv:  preconditioner-vector product, which will also use pdata. The
      default is nothing, which is no preconditioning at all.
      API for ptv is px=ptv(x,pdat) just like kl\\_gmres

Keyword arguments

pdata: precomputed data. The default is nothing, but that ain't gonna
        work well for nonlinear equations.

side: left or right preconditioning. The default is "right".

lmaxit: maximum number of linear iterations. The default is 10.

Output:\n
A named tuple (sol, reshist, ...)

This part is not finished. Watch this space.


"""
function kl_bicgstab( x0, b, atv, V, eta, ptv = nothing; lmaxit = 10, 
      pdata = nothing, side = "right",)
    rhs = V
    rhs .= b
    if side == "right" || ptv == nothing
        itsleft = false
    else
        itsleft = true
        rhs .= ptv(rhs, pdata)
    end
    linsol = copy(b)
    y0 = copy(x0)
    Kpdata = (pdata = pdata, side = side, ptv = ptv, atv = atv, linsol = linsol)
    bout = bicgstab_base(y0, rhs, Katv, eta; 
              lmaxit = lmaxit, pdata = Kpdata)
    if side == "left" || ptv == nothing
        return bout
    else
        sol = y0
        sol .= ptv(sol, pdata)
        return (sol = sol, reshist = bout.reshist)
    end
end

"""
bicgstab_base(x0, rhs, atv, eta; 
         lmaxit = 10, pdata = nothing)

Base BiCGSTAB. Overwrites initial iterate and right hand side.
"""
function bicgstab_base(x0, rhs, atv, eta; lmaxit = 10, pdata = nothing)
    r = rhs
    x = x0
    (norm(x0) == 0.0) || (r .-= atv(x0, pdata))
    k = 0
    rho = zeros(lmaxit + 2)
    rho[1] = 1.0
    rho[2] = r' * r
    alpha = 1.0
    omega = 1.0
    r0 = copy(r)
    rnorm = norm(r0)
    v = zeros(size(x0))
    p = zeros(size(x0))
    s = zeros(size(x0))
    t = zeros(size(x0))
    tol = eta * norm(rhs)
    k = 0
    reshist = []
    push!(reshist, rnorm)
    while rnorm > tol && k < lmaxit
        k += 1
        abs(omega) > 0 || error("Breakdown omega = 0")
        beta = (rho[k+1] / rho[k]) * (alpha / omega)
        axpy!(-omega, v, p)
        #        p .= r + beta * (p - omega * v)
        #        p .= r + beta * p
        axpby!(1.0, r, beta, p)
        v .= atv(p, pdata)
        tau=r0'*v
        abs(tau) > 0 || error("Breakdown r0'*v = 0")
        alpha = rho[k+1] / tau
        #        s .= r - alpha * v
        copy!(s, r)
        axpy!(-alpha, v, s)
        t .= atv(s, pdata)
        norm(t) > 0 || error("Breakdown t = 0")
        omega = (t' * s) / (t' * t)
        rho[k+2] = -omega * (r0' * t)
        #        r .= s - omega * t
        copy!(r, s)
        axpy!(-omega, t, r)
        #        x .= x + alpha * p + omega * s
        copy!(t, s)
        axpby!(alpha, p, omega, t)
        #        x .= x + t
        x .+= t
        rnorm = norm(r)
        push!(reshist, rnorm)
    end
    return (sol = x, reshist = reshist)
end
