function kl_bicgstab(
    x0,
    b,
    atv,
    eta,
    ptv = nothing;
    lmaxit = 10,
    pdata = nothing,
    side = "right",
)
    rhs = copy(b)
    if side == "right" || ptv == nothing
        itsleft = false
    else
        itsleft = true
        rhs .= ptv(rhs, pdata)
    end
    linsol = copy(b)
    y0 = copy(x0)
    Kpdata = (pdata = pdata, side = side, ptv = ptv, atv = atv, linsol = linsol)
    bout = kl_bicgstab_base(y0, rhs, Katv, eta; 
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
kl_bicgstab_base(x0, rhs, atv, eta; 
         lmaxit = 10, pdata = nothing)

Base BiCGSTAB. Overwrites initial iterate and right hand side.
"""
function kl_bicgstab_base(x0, rhs, atv, eta; lmaxit = 10, pdata = nothing)
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
        beta = (rho[k+1] / rho[k]) * (alpha / omega)
        axpy!(-omega, v, p)
        #        p .= r + beta * (p - omega * v)
        #        p .= r + beta * p
        axpby!(1.0, r, beta, p)
        v .= atv(p, pdata)
        alpha = rho[k+1] / (r0' * v)
        #        s .= r - alpha * v
        copy!(s, r)
        axpy!(-alpha, v, s)
        t .= atv(s, pdata)
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
