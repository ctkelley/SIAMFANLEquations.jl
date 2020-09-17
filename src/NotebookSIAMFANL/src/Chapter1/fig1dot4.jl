"""
fig1dot4()

This is the Julia code for Figures 1.1 and 1.4 of the 
Green book. 

I documented the solver in nsolsc.jl. The challenge here
is to get the plots to do what I want. I do not recommend your
studying this code.

The call to the solver is pretty simple

global_hist = nsolsc(1.0, atan, fpatan; rtol = 1.e-8, maxit = 20)

and fpatan(x) = 1.0/(1.0+x^2)

is in this .jl file. 

"""
function fig1dot4()
#
# Figure 1.4. 
#
        global_hist = nsolsc(
            atan,
            10.0,
            fpatan; rtol = 1.e-8, atol = 1.e-10, maxit = 20,
            armfix=true
        )
figure(3)
        rval = global_hist.history
        itc=length(rval)
        ival = 0:itc-1
        semilogy(ival, abs.(rval), "k-")
        plt.text(ival[2], abs(rval[2] * 5), "3")
        plot(ival[2], abs(rval[2]), "ko")
        plt.text(ival[3], abs(rval[3] * 5), "3")
        plot(ival[3], abs(rval[3]), "ko")
        plt.text(ival[4], abs(rval[4] * 5), "2")
        plot(ival[4], abs(rval[4]), "ko")
        plt.text(ival[5], abs(rval[5] * 5), "2")
        plot(ival[5], abs(rval[5]), "ko")
        ylim(1.e-11, 1000)
        xlim(0, 12)
        ylabel("Absolute Nonlinear Residual")
        xlabel("Nonlinear iterations")
        title("Figure 1.4 from print book")
    end

