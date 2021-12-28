"""
fpatan is the derivative of atan
used in the demo above
"""
function fpatan(x)
    return 1.0 / (1.0 + x * x)
end
