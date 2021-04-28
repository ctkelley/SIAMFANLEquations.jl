"""
pdeF!.jl

This file contains everything you need to run the Ellptic PDE examples.  
This includes the version with an explict sparse matrix Jacobian and
the fixed point formulations using the fish2d.jl preconditioner.

I've also parked the exact solution in here so you can do the grid refinement
study.

Look at pdeinit for the construction of the precomputed data. There is
a lot of it.

If you only want to run the examples, you should not have to look
at the code.
"""
### And now for the functions ...
"""
pdeF!(FV, u, pdata)

Residual using sparse matrix-vector multiplication
"""
function pdeF!(FV, u, pdata)
    D2 = pdata.D2
    CV = pdata.CV
    rhs = pdata.RHS
    p1 = pdata.jvect1
    #    FV .= D2 * u + 20.0 * u .* (CV * u) - rhs
    #    FV .= D2*u
    #    p1 .= CV*u
    mul!(FV, D2, u)
    mul!(p1, CV, u)
    p1 .*= 20.0
    p1 .*= u
    FV .+= p1
    FV .-= rhs
end

"""
pdeJ!(FP, F, u, pdata)

Sparse matrix Jacobian. The package does not do its own sparse
differencing. The Jacobian for this problem is easy enough to 
compute analytically.
"""
function pdeJ!(FP, F, u, pdata)
    D2 = pdata.D2
    CV = pdata.CV
    CT = pdata.CT
    cu = CV * u
    #DC=spdiagm(0 => 20*cu); DU=spdiagm(0 => 20*u)
    DC = Diagonal(20 * cu)
    DU = Diagonal(20 * u)
    #
    # The easy way to compute the Jacobian is 
    #FP .= D2 + DU*CV + DC
    # but you allocate yourself silly with that one.
    # So we preallocate room for DU*CV in CT and sum the terms for FP
    # one at a time. I have to use Diagonal instead of spdiagm if I want
    # mul! to work fast.
    #
    FP .= D2
    FP .+= DC
    mul!(CT, DU, CV)
    #CT .= CV; lmul!(DU,CT); 
    FP .+= CT
    # I should be able to do mul!(FP,DU,CV), but it's 1000s of times slower.
end

"""
Jvec2d(v, FS, u, pdata)

Analytic Jacobian-vector product for PDE example
"""
function Jvec2d(v, FS, u, pdata)
    D2 = pdata.D2
    CV = pdata.CV
    CT = pdata.CT
    jvec = pdata.jvec
    p1 = pdata.jvect1
    #    jvec .= D2 * v
    #    p1 .= CV * u
    mul!(jvec, D2, v)
    mul!(p1, CV, u)
    p1 .*= 20.0
    p1 .*= v
    jvec .+= p1
    #    p1 .= CV * v
    mul!(p1, CV, v)
    p1 .*= 20.0
    p1 .*= u
    jvec .+= p1
    return jvec
end


"""
pdeinit(n)

collects the precomputed data for the elliptic pde example. This 
includes 

- the sparse matrix representation of the operators, 
- the right side of the equation,
- the exact solution,
- the data that the fft-based fast Poisson solver (fish2d) needs
"""
function pdeinit(n)
    # Make the grids
    n2 = n * n
    h = 1.0 / (n + 1.0)
    x = collect(h:h:1.0-h)
    # collect the operators
    D2 = Lap2d(n)
    DX = Dx2d(n)
    DY = Dy2d(n)
    CV = (DX + DY)
    # I need a spare sparse matrix to save allocations in the Jacobian computation
    CT = copy(CV)
    # Exact solution and its derivatives
    uexact = solexact(x)
    dxe = dxexact(x)
    dye = dyexact(x)
    d2e = l2dexact(x)
    dxv = reshape(dxe, (n2,))
    dyv = reshape(dye, (n2,))
    d2v = reshape(d2e, (n2,))
    uv = reshape(uexact, (n2,))
    fdata = fishinit(n)
    # The right side of the equation
    RHS = d2v + 20.0 * uv .* (dxv + dyv)
    # preallocate a few vectors
    jvec = zeros(n2)
    jvect1 = zeros(n2)
    # Pack it and ship it.
    pdedata =
        (D2 = D2, CV = CV, CT = CT, RHS = RHS, jvec, jvect1, fdata = fdata, uexact = uexact)
end

"""
This collection of functions 
builds u, u_x, u_y, and the negative Laplacian for the 
example problem in the book. Here
u(x,y) = 10 x y (1-x)(1-y) exp(x^4.5)

which is the example from FA01.
"""

function w(x)
    w = 10.0 * x .* (1.0 .- x) .* exp.(x .^ (4.5))
end

function wx(x)
    wx = 4.5 * (x .^ (3.5)) .* w(x) + 10.0 * exp.(x .^ (4.5)) .* (1.0 .- 2.0 * x)
end

function wxx(x)
    wxx =
        (4.5 * 3.5) * (x .^ (2.5)) .* w(x) +
        4.5 * (x .^ (3.5)) .* wx(x) +
        +10.0 * 4.5 * (x .^ (3.5)) .* exp.(x .^ (4.5)) .* (1.0 .- 2.0 * x) +
        -20.0 * exp.(x .^ (4.5))
end

function v(x)
    v = x .* (1.0 .- x)
end

function vx(x)
    vx = 1.0 .- 2.0 * x
end

function vxx(x)
    vxx = -2.0 * ones(size(x))
end

function solexact(x)
    solexact = w(x) * v(x)'
end

function l2dexact(x)
    l2dexact = -(w(x) * vxx(x)') - (wxx(x) * v(x)')
end

function dxexact(x)
    dxexact = wx(x) * v(x)'
end

function dyexact(x)
    dxexact = w(x) * vx(x)'
end
