"""
fig1dot6()

Make Figure 1.6 in the print book.

f(x) = x^3 - lambda x = 0 with lambda = 1/2

The answer(s!!) are x=0, -sqrt(2)/2, and sqrt(2)/2.

The initial iterate is x0=.1. You'll see Newton converge to the 
unstable soluution x=0 and ptc do it's job and get x=sqrt(2).

You might want to experiment with different values of dt0. I picked
a big one so that I would not have to spend forever waiting for the
transients to go away. 

"""
function fig1dot6()
x0=.1
ustable=.5*sqrt(2.0)
uunstable=0.0
#
# Convergence to two different solutions.
#
ptcdata1=ptcsolsc(sptest,x0; dt0=1.0, rtol=1.e-12)
ptcerr=ptcdata1.solhist.-ustable
ptcfun=ptcdata1.history
newtdata=nsolsc(sptest,x0)
nival=length(newtdata.history)
nival=0:nival-1
itlim=length(ptcfun)
pival=0:itlim-1
newterr=abs.(newtdata.solhist)
newtresid=newtdata.history
#
# Newton errors vs PTC errors
#
figure(1)
subplot(1,2,1)
title("Figure 1.6 from print book")
semilogy(pival,abs.(ptcerr),"k--",nival,abs.(newterr),"k-")
PTClabel=L"$\Psi$TC"
legend((PTClabel,"Newton"))
xlabel("Iterations",fontsize="12")
gylabele=L"$|x-x^*|$"
ylabel(gylabele,fontsize="12")
axis([0, itlim, 1.e-17, 1.0])
#
# Look at how |f| and dt almost mirror each other, as the forumla
# says they should.
#
subplot(1,2,2)
semilogy(pival,abs.(ptcfun),"k--",nival,newtresid,"k-")
xlabel("Iterations",fontsize="12")
gylabelfdt=L"$|f|$"
ylabel(gylabelfdt,fontsize="12")
legend((PTClabel,"Newton"))
axis([0, itlim, 1.e-15, 1.0])
PyPlot.tight_layout()

return (ptcdata=ptcdata1, newtdata=newtdata)
end
