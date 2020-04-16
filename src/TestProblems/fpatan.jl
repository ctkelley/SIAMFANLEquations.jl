"""
fpatan is the derivative of atan
used in the demo above
"""
function fpatan(x)
    return 1 / (1 + x * x)
end
