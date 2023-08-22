"""
simple!(FV,x)
This is the function for Figure 2.1
It also shows up in CI
"""
function simple!(FV, x)
    FV[1] = x[1] * x[1] + x[2] * x[2] - 2.0
    FV[2] = exp(x[1] - 1) + x[2] * x[2] - 2.0
    #
    # The return FV is important
    #
    return FV
end

function jsimple!(JacV, FV, x)
    JacV[1, 1] = 2.0 * x[1]
    JacV[1, 2] = 2.0 * x[2]
    JacV[2, 1] = exp(x[1] - 1)
    JacV[2, 2] = 2 * x[2]
    #
    # The return JacV is important
    #
    return JacV
end

"""
JVsimple(v, FV, x)

Jacobian-vector product for simple!. There is, of course, no reason 
to use Newton-Krylov for this problem other than CI or demonstrating 
how to call nsoli.jl.
"""
function JVsimple(v, FV, x)
    jvec = zeros(2)
    jvec[1] = 2.0 * x' * v
    jvec[2] = v[1] * exp(x[1] - 1.0) + 2.0 * v[2] * x[2]
    #
    # The return jvec is important
    #
    return jvec
end
