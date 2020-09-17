"""
function fig1dot1()

This is the Julia code for Figure 1.1 in the book.

I documented the solver in nsolsc.jl. The challenge here
is to get the plots to do what I want. I do not recommend your
studying this code.

The call to the solver is pretty simple

local_hist = nsolsc(1.0, atan, fpatan; rtol = 1.e-8, maxit = 20)

fp_atan(x) = 1.0/(1.0+x^2) comes from SIAMFANLTestProblems

"""
function fig1dot1()

    #
    # Figure 1.1 from the print book.
    #

    local_hist = nsolsc(atan,1.0,fpatan; rtol = 1.e-8, maxit = 20)
    iplot = true
    if iplot
        figure(1)
        xval = local_hist.solhist
        yval=atan.(xval)
        xtval = -1.2:.01:1.2
        ytval = atan.(xtval)
        ztval = zeros(size(xtval))
        plot(xval, yval, "ko", xtval, ytval, "k-", xtval, ztval, "k-")
        t = 0:.1:1
        xv = xval[1] * (1.0 .- t) + t .* xval[2]
        yv = yval[1] * (1.0 .- t)
        plot(xv, yv, "k-")
        yv = yval[2] * t
        xv = xval[2] .* ones(length(t), 1)
        plot(xv, yv, "k-")
        xv = xval[2] * (1.0 .- t) + t .* xval[3]
        yv = yval[2] * (1.0 .- t)
        plot(xv, yv, "k-")
        yv = yval[3] * t
        xv = xval[3] .* ones(length(t), 1)
        plot(xv, yv, "k-")
        plt.text(1, 0.64, L"(x_0,y_0)")
        plt.text(-.2, .45, L"y=m_0(x)")
        plt.text(xval[2] - .1, .125, L"(x_1,0)")
        plt.text(xval[2] + .05, yval[2] - .1, L"(x_1,y_1)")
        plt.text(xval[3] + .1, -.1, L"(x_2,0)")
        plt.text(xval[3] + .1, yval[3], L"(x_2,y_2)")
        plt.text(-.1, -.25, L"y=m_1(x)")
        plt.text(-.07, .05, L"x^*")
        plot(xval[2], 0, "ko")
        plot(xval[3], 0, "ko")
        ylabel("atan(x)")
        xlabel("x")
        title("Figure 1.1 from print book")
     end
end

