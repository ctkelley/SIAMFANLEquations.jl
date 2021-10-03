function reldiff(x, y)
    p = (x - y) ./ abs.(x)
    return norm(p, Inf)
end
