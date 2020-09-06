"""
simple!(FV,x)
This is the function for Figure 2.1
It also shows up in CI
"""
function simple!(FV,x)
FV[1]=x[1]*x[1] + x[2]*x[2] -2.0;
FV[2]=exp(x[1]-1) + x[2]*x[2] - 2.0;
end

function jsimple!(JacV,FV,x)
JacV[1,1]=2.0*x[1]
JacV[1,2]=2.0*x[2]
JacV[2,1]=exp(x[1]-1) 
JacV[2,2]=2*x[2]
end
