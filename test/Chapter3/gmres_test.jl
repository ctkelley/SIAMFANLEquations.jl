"""
gmres_test.jl

Tests the linear GMRES code, klgmres. 

This is for CI only. Nothing to see here. Move along.
"""

function gmres_test()
pass3 = test3x3()
passint = test_integop(40)
passr1 = testR1()
passgm = pass3 && passint && passr1
return passgm
end

function test3x3()
    A = [0.001 0 0; 0 0.0011 0; 0 0 1.e4]
    V = zeros(3, 10)
    b = [1.0; 1.0; 1.0]
    x0 = zeros(3,)
    eta = 1.e-10
    passgm = true
    R = []
    rightsize = [10, 6, 5, 4]
    Methods = ("cgs1", "mgs1", "mgs2", "cgs2")
    i = 1
    for orth in Methods
        gout = kl_gmres(x0, b, atv, V, 1.e-10; pdata=A, orth = orth)
        ithist = gout.reshist
        lhist = length(ithist)
        push!(R, ithist)
        resnorm = norm(A * gout.sol - b)
        locpass = (resnorm < 1.e-8) && (lhist == rightsize[i])
        locpass || println(
            "failure at orth = ",
            orth,
            ", lhist = ",
            lhist,
            ", Res norm = ",
            resnorm,
        )
        passgm = passgm && locpass
        i += 1
    end
#    return (pass = passgm, RH = R)
    return passgm
end

function testR1()
    A = Float64.([1 2 3 4 5])
    E = A' * A
    A = I + E
    b = ones(5,)
    x0 = zeros(5,)
    V = zeros(5, 4)
    gout = kl_gmres(x0, b, atv, V, 1.e-7; pdata=A)
    lhist = length(gout.reshist)
    nerr = norm(A * gout.sol - b, Inf)
    pass = (lhist == 3) && (nerr < 1.e-14)
    pass || println("Rank one test fails")
    return pass
end

function test_integop(n)
    pdata = integopinit(n)
    f = pdata.f
    ue = pdata.xe
    u0 = zeros(size(f))
    V = zeros(n, 20)
    Methods = ("cgs1", "mgs1", "mgs2", "cgs2")
    pass = true
    for orth in Methods
        goutinteg = kl_gmres(u0, f, integop, V, 1.e-10; pdata=pdata, orth = orth)
        errn = norm(goutinteg.sol - ue, Inf)
        rhist = goutinteg.reshist
        lhist = length(rhist)
        rred = rhist[4] ./ rhist[1]
        lpass = (errn < 1.e-14) && (rred < 1.e-14) && (lhist == 4)
        lpass || println("Failure with orth = ", orth)
        pass = pass && lpass
    end
    pass || println("Integral operator test fails.")
    return pass
end

function atv(x, A)
    return A * x
end


function integop(u, pdata)
    K = pdata.K
    f = pdata.f
    return u - K * u
end

function integopinit(n)
    h = 1 / n
    x = collect(0.5*h:h:1.0-0.5*h)
    K = zeros(n, n)
    for j = 1:n
        for i = 1:n
            K[i, j] = ker(x[i], x[j])
        end
    end
    K .*= h
    sol = exp.(x) .* log.(2.0 * x .+ 1.0)
    f = sol - K * sol
    pdata = (K = K, xe = sol, f = f)
    return pdata
end


function ker(x, y)
    ker = 0.1 * sin(x + exp(y))
end
