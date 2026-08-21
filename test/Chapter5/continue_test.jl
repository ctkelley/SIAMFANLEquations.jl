function continue_test(n = 128)
    v1ok = test_v1(n)
    v1ok || println("test continue v1 fails")
    v2ok = test_PAC(n + 1)
    v2ok || println("test continue v2fails")
    return continueok = v1ok && v2ok
end

function test_v1(n = 129)
    version = "orig"
    (pval, nval, x, lambda) = heq_continue(n; version = version)
    dpath = path_test.(pval, nval)
    del1 = dpath[1:(end - 1)]
    del2 = dpath[end]
    v1_pass = (norm(del1, Inf) < 4.0e-9) && (del2 < 1.5e-4)
    return v1_pass
end

function test_PAC(n = 128)
    version = "pac"
    (pval, nval, x, lambda) = heq_continue(n; version = version)
    dpath = path_test.(pval, nval)
    nsingular = argmax(dpath)
    del2 = dpath[nsingular]
    del1 = [dpath[1:(nsingular - 1)]; dpath[(nsingular + 1):end]]
    v2_pass = (norm(del1, Inf) < 1.0e-5) && (del2 < 1.0e-4)
    return v2_pass
end

function path_test(pval, nval)
    if pval > 0
        rp = (1.0 + sqrt.(1.0 .- pval)) / (0.5 .* pval)
        rm = (1.0 - sqrt.(1.0 .- pval)) / (0.5 .* pval)
    else
        rm = 1.0
    end
    (nval .<= 2) ? dp = abs.(nval - rm) : dp = abs.(nval - rp)
    return dp
end
