function simple!(FV,x)
FV[1]=x[1]*x[1] + x[2]*x[2] -2.0;
FV[2]=exp(x[1]-1) + x[2]*x[2] - 2.0;
end

function jsimple!(FV,JV,x)
JV[1,1]=2.0*x[1]
JV[1,2]=2.0*x[2]
JV[2,1]=exp(x[1]-1) 
JV[2,2]=2*x[2]
end
