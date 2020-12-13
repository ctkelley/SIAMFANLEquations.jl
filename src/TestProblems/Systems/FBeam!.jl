"""
FBeam!(FV, U, bdata)
Function evaluation for PTC example.
F(u) = -u'' - lambda sin(u)
"""
function FBeam!(FV, U, bdata)
    D2 = bdata.D2
    lambda = bdata.lambda
    su = lambda * sin.(U)
    FV .= (D2 * U - su)
end

"""
BeamJ!(FP,FV,U,bdata)

Jacobian for the beam problem
F(u) = -u'' - lambda sin(u)
so
F'(u) w  = D2 w - lambda cos(u) w
"""
function BeamJ!(FP, FV, U, bdata)
    D2 = bdata.D2
    lambda = bdata.lambda
    cu = lambda * cos.(U)
    n = length(U)
    zr = zeros(n - 1)
    FP .= D2 - Tridiagonal(zr, cu, zr)
end

"""
BeamtdJ!(FP, FV, U, bdata)

Jacobian evaluation for the time-dependent beam problem.
If F^n(w) = w - u_n + dt F(w) = 0 then
F^n(w)' = I + dt F'(w)
"""
function BeamtdJ!(FP, FV, U, bdata)
    FP .= BeamJ!(FP, FV, U, bdata)
    dt = bdata.dt
    FP .= I + dt * FP
end

"""
FBeamtd!(FV, U, bdata)

Function evaluation for the time-dependent beam problem.
The implicit Euler step for u_t = - F(u) is
u_{n+1} = u_n - dt F(u_{n+1})
so the nonlinear equation is
F^n(w) = w - u_n + dt F(w) = 0
"""
function FBeamtd!(FV, U, bdata)
    un = bdata.UN
    dt = bdata.dt
    FV .= FBeam!(FV, U, bdata)
    dU = U - un
    FV .= dU + dt * FV
    #axpby!(1.0,dU,dt,FV)
end

"""
beaminit(n,dt,lambda=20.0)

Set up the beam problem with n interior grid points.
"""

function beaminit(n, dt, lambda = 20.0)
    D2 = Lap1d(n)
    dx = 1.0 / (n + 1)
    x = collect(dx:dx:1.0-dx)
    UN = zeros(size(x))
    bdata = (D2 = D2, x = x, dx = dx, dt = dt, lambda = lambda, UN = UN)
end

"""
Lap1d(n)

returns -d^2/dx^2 on [0,1] zero BC
"""
function Lap1d(n)
    dx = 1 / (n + 1)
    d = 2.0 * ones(n)
    sup = -ones(n - 1)
    slo = -ones(n - 1)
    D2 = Tridiagonal(slo, d, sup)
    D2 = D2 / (dx * dx)
    return D2
end
