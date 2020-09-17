"""
fig1dot3()

This is the Julia code for the tan(x) = x example.

This code makes Figure 1.3 in the print book.

See the notebook.

"""
function fig1dot3()
    kwnewt=(maxit=14, rtol=1.e-17, atol=1.e-17, printerr=false,  
            stagnationok=true)
    kwchord=(maxit=14, rtol=1.e-17, atol=1.e-17, printerr=false,  
            stagnationok=true,solver="chord")
    kwsec=(maxit=6, rtol=1.e-17, atol=1.e-17, printerr=false,  
            stagnationok=true,solver="secant")
    nnout=nsolsc(ftanx,4.5; kwnewt...)
    lnn=length(nnout.history)
    nncounter=0:lnn-1
    chout=nsolsc(ftanx,4.5; kwchord...)
    lnc=length(chout.history)
    nccounter=0:lnc-1
    scout=nsolsc(ftanx,4.5; kwsec...)
    lns=length(scout.history)
    sccounter=0:lns-1
semilogy(nncounter, nnout.history,"k-",
nccounter, chout.history,"k--",
sccounter, scout.history,"k-."
)
ylabel("Log Absolute Nonlinear Residual")
xlabel("Nonlinear Iterations")
legend(("Newton","Chord","Secant"))
title("Figure 1.3 from print book")
#return nnout
end

#function ftanx(x)
#   return tan(x) - x
#end

#function ftanxp(x)
#   return sec(x)^2 - 1
#end
