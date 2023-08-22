"""
knowsdt_test()

Test the jknowsdt keyword with a simple problem

u' = u (mu - u^2)
v' = u (mu - u^2) - 7 v

Remember that PTC thinks in terms of x' = - F(x)

Here mu=4 and the initial data are (u0, v0)=(.1, 10) so the correct
stable solution is (u,v)=(2,0).

"""
function knowsdt_test()
    u0 = [0.1; 10.0]
    FU = zeros(2)
    JV = zeros(2, 2)
    # Jacobian does not know about dt + finite difference
    pout = ptcsol(Fode!, u0, FU, JV; delta0 = 0.01, maxit = 100)
    # Analytc Jacobian does not know about dt
    pout2 = ptcsol(Fode!, u0, FU, JV, Jval!; delta0 = 0.01, maxit = 100)
    # Analytc Jacobian knows about dt
    pout3 = ptcsol(Fode!, u0, FU, JV, Jval2!; delta0 = 0.01, maxit = 100, jknowsdt = true)
    #
    # Collect the output and figure out if you did things right.
    #
    hist = pout.history
    hist2 = pout2.history
    hist3 = pout3.history
    sol = pout.solution
    sol2 = pout2.solution
    sol3 = pout3.solution
    ustar = [2.0; 0.0]
    dtdiff = norm(sol2 - sol3, Inf) + norm(hist2 - hist3, Inf)
    soldiff = norm(sol - sol2, Inf)
    stardiff = norm(sol3 - ustar, Inf)
    hdiff = norm(hist - hist2, Inf)
    dtpass = (soldiff < 1.e-9) && (hdiff < 1.e-6) && (dtdiff < 1.e-15)
    dtpass || println("knowsdt_test fails")
    return dtpass
end

function Fode!(FS, x)
    #
    # Dynamics are x' = -FS(x), so the zero solution is unstable for mu > 0
    #
    mu = 4.0
    FS[1] = -x[1] * (mu - x[1] * x[1])
    FS[2] = FS[1] + 7.0 * x[2]
    return FS
end

function Jval2!(JV, FU, x, dt)
    JV .= Jval!(JV, FU, x)
    JV .= JV + (1.0 / dt) * I
end

function Jval!(JV, FU, x)
    mu = 4.0
    JV[1, 1] = -(mu - 3.0 * x[1] * x[1])
    JV[1, 2] = 0.0
    JV[2, 1] = JV[1, 1]
    JV[2, 2] = 7.0
    return JV
end
