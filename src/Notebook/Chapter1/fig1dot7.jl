"""
fig1dot7()

Make Figure 1.7 in the print book. Examine how the wrong
dt0 will cost you.

f(x) = x^3 - lambda x = 0 with lambda = 1/2

The answer(s!!) are x=0, -sqrt(2)/2, and sqrt(2)/2.

The initial iterate is x0=.1. You'll see Newton converge to the 
unstable soluution x=0 and ptc do it's job and get x=sqrt(2).

You might want to experiment with different values of dt0. I picked
a big one so that I would not have to spend forever waiting for the
transients to go away. 

"""
function fig1dot7()
x0=.1
ustable=.5*sqrt(2.0)
uunstable=0.0
#
# Look at how dt0 changes things.
#
dtlist=[1.0, 1.e-1, 1.e-2, 1.e-3, 1.e-4]
labels=[ 
L"$dt_0 = 1.0$", 
L"$dt_0 = 10^{-1}$",
L"$dt_0 = 10^{-2}$",
L"$dt_0 = 10^{-3}$"
]
outdata=[]
for id=1:4
   dti=10.0^(1-id)
   ptcdata1=ptcsolsc(sptest,x0; dt0=dti, rtol=1.e-6, maxit=1000000)
   push!(outdata,ptcdata1.history)
end
plothist(outdata,labels; semilogflag=false, ptitle="Figure 1.7 in print book")
#plothist(outdata,labels; semilogflag=false, ptitle="")
return outdata
end
