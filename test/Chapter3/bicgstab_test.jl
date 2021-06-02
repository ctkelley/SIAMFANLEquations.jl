"""
bicgstab_test.jl

Tests the linear BiCGSTAB code, kl_bicgstab.

This is for CI only. Nothing to see here. Move along.
"""

function bicgstab_test()
    pass3 = test3x3()
    passint = test_integop()
    passr1 = testR1()
    passbicgs = pass3 && passint && passr1 
    return passbicgs
end

function test3x3()
    A = [0.001 0 0; 0 0.0011 0; 0 0 1.e4];
    V = zeros(3); b = [1.0; 1.0; 1.0]; x0 = zeros(3);
    eta = 1.e-10
    gout = kl_bicgstab(x0, b, atv, V, 1.e-10; pdata = A)
    pass = (length(gout.reshist)== 5) && 
            (norm(A*gout.sol - b,Inf) < 1.e-12) && gout.idid && (gout.lits==4)
#    return (gout=gout, pass=pass)
    return pass
end

function testR1()
    A = Float64.([1 2 3 4 5])
    E = A' * A
    A = I + E
    b = ones(5)
    x0 = zeros(5)
    V = zeros(5)
    gout = kl_bicgstab(x0, b, atv, V, 1.e-7; pdata = A)
    pass = (length(gout.reshist)== 3) && (norm(A*gout.sol - b,Inf) < 1.e-12)
#    return (gout=gout, pass=pass)
    return pass
end

function test_integop(n=100)
    pdata = integopinit(n)
    f = pdata.f
    ue = pdata.xe
    u0 = zeros(size(f)); V = zeros(size(f))
    gout = kl_bicgstab(u0, f, integop, V, 1.e-10; pdata = pdata)
    realres = (I - pdata.K)*gout.sol - f
pass = ( (norm(realres,Inf) < 1.e-15) && (length(gout.reshist)==4))
#    return (gout=gout, pass=pass)
     return pass
end


function atv(x, A)
    return A * x
end


function integop(u, pdata)
    K = pdata.K
#    f = pdata.f
    return u - K * u
end

function integopinit(n)
    h = 1 / n
    X = collect(0.5*h:h:1.0-0.5*h)
    K = [ker(x, y) for x in X, y in X]
    #    K = zeros(n, n)
    #    for j = 1:n
    #        for i = 1:n
    #            K[i, j] = ker(x[i], x[j])
    #        end
    #    end
    K .*= h
    #    sol = exp.(x) .* log.(2.0 * x .+ 1.0)
    #    sol = usol.(X)
    sol = [usol(x) for x in X]
    f = sol - K * sol
    pdata = (K = K, xe = sol, f = f)
    return pdata
end

function usol(x)
    return exp.(x) .* log.(2.0 * x .+ 1.0)
end

function ker(x, y)
    ker = 0.1 * sin(x + exp(y))
end
