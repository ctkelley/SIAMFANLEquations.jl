"""
Test problem to give the chord method problems.
The line search will fail in the middle of the iteration
and demand a recompute of the derivative.
"""
function linatan(x)
alpha=.01
return (1.0+alpha*x)*atan(x)
end
