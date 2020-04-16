function difffp(x, f, fc, h)
    df = (f(x + h) - fc) / h
    return df
end
