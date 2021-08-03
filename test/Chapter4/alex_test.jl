"""
alex_test()

Test for duplication of Table 
"""
function alex_test()
    (historye, condhiste, alphanorme) = vtst()
    u0 = ones(2)
    maxit = 20
    maxm = 2
    vdim = 2 * (maxm + 1)
    Vstore = zeros(2, vdim)
    m = 2
    aout = aasol(alexfp!, u0, m, Vstore; rtol = 1.e-10)
    alexerr = (
        reldiff(aout.history, historye) +
        reldiff(aout.stats.condhist, condhiste) +
        reldiff(aout.stats.alphanorm, alphanorme)
    )
    solerr = reldiff(aout.history, historye)
    solok = (solerr < 1.e-5)
    solok || println("alex solution error")
    conderr = reldiff(aout.stats.condhist, condhiste)
    condok = (conderr < 1.e-5)
    condok || println("alex condition error")
    aerr = reldiff(aout.stats.alphanorm, alphanorme)
    aok = (aerr < 1.e-5)
    aok || println("alex coefficient error")
    aout.idid || println("idid is wrong for m=2")
    alexok2 = solok && condok && aok
    aout = aasol(alexfp!, u0, 0, Vstore; rtol = 1.e-10)
    aout.idid && println("idid is wrong for m=0")
    alexok0 = ~aout.idid && (aout.errcode == 10)
    alexok = alexok2 && alexok0
    return alexok
end


function reldiff(x, y)
    p = (x - y) ./ abs.(x)
    return norm(p, Inf)
end


function vtst()
    historye = [
        6.50111e-01
        4.48661e-01
        2.61480e-02
        7.25389e-02
        1.53107e-04
        1.18512e-05
        1.82476e-08
        1.04804e-13
    ]
    condhiste = [
        1.00000e+00
        2.01556e+10
        1.37776e+09
        3.61344e+10
        2.54947e+11
        3.67672e+10
    ]
    alphanorme = [
        1.00000e+00
        4.61720e+00
        2.15749e+00
        1.18377e+00
        1.00000e+00
        1.00171e+00
    ]
    return (historye, condhiste, alphanorme)
end

function alexfp!(G, u)
G[1]=cos(.5*(u[1]+u[2]))
G[2]=G[1]+ 1.e-8 * sin(u[1]*u[1])
return G
end
