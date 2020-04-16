function fpeval_newton(x, f, fc, fp, h)
    if fp == difffp
        df = difffp(x, f, fc, h)
    else
        df = fp(x) 
    end
    return df 
end 
