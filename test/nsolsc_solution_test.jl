"""
nsolsc_solution_test

Test nsolsc with the atan function. Check answers and iteration stats.
"""
function nsolsc_solution_test()
# 
# Local convergence with forward difference derivative
#
sdatal=nsolsc(1.0, atan)
solok=(abs(sdatal.solution) < 1.e-8)
funok=(abs(sdatal.functionval) < 1.e-8)
hs=size(sdatal.history)
histok=(hs[1]==5)
locok = funok && solok && histok
if locok 
   println("local FD ok")
end
#
# Local convergence with analytic derivative
#
sdataa=nsolsc(1.0, atan; fp=fpatan)
solok=(abs(sdataa.solution) < 1.e-8)
funok=(abs(sdataa.functionval) < 1.e-8)
hs=size(sdataa.history)
histok=(hs[1]==5)
analyticok = funok && solok && histok
if analyticok
   println("analytic derivative ok")
else
   println(sdatta)
end
#
# Global convergence
#
sdatag=nsolsc(10.0, atan)
solok=(abs(sdatag.solution) < 1.e-8)
funok=(abs(sdatag.functionval) < 1.e-8)
hs=size(sdatag.history)
histok=(hs[1]==12)
globok = funok && solok && histok
if globok 
   println("global FD ok")
end
#
# Local convergence with secant method
#
sdatas=nsolsc(1.0, atan; solver="secant")
solok=(abs(sdatas.solution) < 1.e-10)
funok=(abs(sdatas.functionval) < 1.e-10)
hs=size(sdatas.history)
histok=(hs[1]==6)
secantok = funok && solok && histok
if secantok
   println("secant ok")
end
#
# Initialize secant method when x0=0
#
zedata=nsolsc(0.0,fcos;solver="secant",rtol=1.e-9)
solution=7.390851333858823e-01
solok=(abs(zedata.solution-solution) < 1.e-9)
funok=(abs(zedata.functionval) < 1.e-9)
hs=size(zedata.history)
histok=(hs[1]==7)
zecok = funok && solok && histok
if zecok
   println("local FD at zero ok")
end
sdatal=nsolsc(200.0, linatan; sham=5,maxit=20,armmax=10,rtol=1.e-9)
solution=-100
solok=(abs(sdatal.solution-solution) < 1.e-8)
funok=(abs(sdatal.functionval) < 1.e-8)
hs=size(sdatal.history)
histok=(hs[1]==4)
shamfastok = funok && solok && histok
if shamfastok
   println("Fast Shamanskii response ok")
end
#
# Test linesearch failure complaints.
#
armfail=nsolsc(10.0,atan; armmax=1)
afok=false
if armfail.idid==false
   afok=true
   println("Armijo failure test passed.")
end
#
# Test residual failure mode and no history.
#
resok=false
resfail=nsolsc(10.0, atan; maxit=2, keepsolhist=false)
if resfail.idid==false
   resok=true
end
#
return locok && globok && secantok && analyticok && zecok && 
       shamfastok && afok && resok
end
