function PTCUpdate(FPS::AbstractArray, FS, x, ItRules, step, residm, dt)
    F! = ItRules.f
    pdata = ItRules.pdata
    #
    FPF = PrepareJac!(FPS, FS, x, ItRules, dt)
    #
    step .= -(FPF \ FS)
    #
    # update solution/function value
    #
    x .= x + step
    EvalF!(F!, FS, x, pdata)
    resnorm = norm(FS)
    #
    # Update dt
    #
    dt *= (residm / resnorm)
    return (x, dt, resnorm)
end

function PTCUpdate(df::Real, fval, x, ItRules, step, residm, dt)
        f = ItRules.f
        dx = ItRules.dx
        fp = ItRules.fp
        df = fpeval_newton(x, f, fval, fp, dx)
        idt = 1.0 / dt
        step = -fval / (idt + df)
        x = x + step
        fval = f(x)
        # SER
        resnorm = abs(fval)
        dt *=  (residm/resnorm)
        return (x, dt, fval, resnorm)
end
