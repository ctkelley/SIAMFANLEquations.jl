"""
fig1dot2()

This is the Julia code for the tan(x) = x example.

This code makes Figure 1.2

See the notebook.

"""
function fig1dot2()
    kwnewt=(maxit=14, rtol=1.e-17, atol=1.e-17, printerr=false,  
            stagnationok=true)
    nnout=nsolsc(ftanx,4.5; kwnewt...)
    lnn=length(nnout.history)
    nncounter=0:lnn-1
semilogy(0:5, nnout.history[1:6],"k-")
ylabel("Log Absolute Nonlinear Residual")
xlabel("Nonlinear Iterations")
plt.text(4.1, 2.e-15,"Stagnation")
ylabel("Log Absolute Nonlinear Residual")
xlabel("Nonlinear Iterations")
title("Figure 1.2 from print book")
#return nnout
end

#function ftanx(x)
#   return tan(x) - x
#end

#function ftanxp(x)
#   return sec(x)^2 - 1
#end
