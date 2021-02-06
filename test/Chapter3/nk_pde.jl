"""
nk_pde(n)

Solve the Elliptic PDE using nsol.jl on an n x n grid. 
"""
function nk_pde()
    # Get some room for the residual
    rtol = 1.e-7
    atol = 1.e-10
    n = 15
    u0 = zeros(n * n)
    FV = copy(u0)
    # Get the precomputed data from pdeinit
    pdata = pdeinit(n)
    # Storage for the Jacobian-vector products
    JV = zeros(n * n, 100)
    # Call the solver with a finite-difference Jac-Vec
    hout = nsoli(
        pdeF!,
        u0,
        FV,
        JV;
        rtol = rtol,
        atol = atol,
        pdata = pdata,
        eta = .1,
        fixedeta = false,
        maxit = 20,
    )
    # Call the solver twice with an analytic Jac-Vec
    hout2 = NsoliPDE(n; fixedeta=false)
    hout3 = NsoliPDE(n; fixedeta=true)
    soldiff = (
        norm(hout3.solution - hout.solution, Inf) +
        norm(hout3.solution - hout.solution, Inf)
    )
    solpass = (soldiff < 1.e-6)
    solpass || println("solution compare fails in nk_pde, ",soldiff)
    histdiffv = (hout.history - hout2.history) ./ hout.history[1]
    histdiff = norm(histdiffv, Inf) 
    histpass = (histdiff < .1)
    histpass || println("history compare fails in nk_pde, ",histdiff)
    cost1 = sum(hout.stats.ijac)
    cost2 = sum(hout2.stats.ijac)
    cost3 = sum(hout3.stats.ijac)
    costpass = (cost1 > 80) && (cost2 > 30) && (cost3 > cost2)
    costpass || println(cost1," ", cost2, "  ", cost3)
    return costpass
end
