"""
fpeval_newton

Evaluates f' by differences or the user's code.

"""
function fpeval_newton(x, f, fc, fp, h)
    if fp == difffp
        df = difffp(x, f, fc, h)
    else
        df = fp(x) 
    end
    return df 
end 


"""
difffp

forward differencing for scalar equations
"""
function difffp(x, f, fc, h)
    df = (f(x + h) - fc) / h
    return df
end
