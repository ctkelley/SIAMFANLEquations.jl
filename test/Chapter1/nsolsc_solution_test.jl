"""
nsolsc_solution_test

Test nsolsc with the atan function. Check answers and iteration stats.
"""
function nsolsc_solution_test()
# 
# Local convergence with forward difference derivative
#
sdatal=nsolsc(atan,1.0)
solok=(abs(sdatal.solution) < 1.e-8)
funok=(abs(sdatal.functionval) < 1.e-8)
hs=size(sdatal.history)
histok=(hs[1]==5)
locok = funok && solok && histok
if ~locok 
   println("local FD fails")
end
#
# Local convergence with analytic derivative
#
sdataa=nsolsc(atan, 1.0, fpatan)
solok=(abs(sdataa.solution) < 1.e-8)
funok=(abs(sdataa.functionval) < 1.e-8)
hs=size(sdataa.history)
histok=(hs[1]==5)
analyticok = funok && solok && histok
if ~analyticok
   println("failure with analytic derivative ")
   println(sdataa)
end
#
# Global convergence
#
sdatag=nsolsc(atan,10.0; maxit=11, armfix=true)
solok=(abs(sdatag.solution) < 1.e-8)
funok=(abs(sdatag.functionval) < 1.e-8)
hs=size(sdatag.history)
histok=(hs[1]==12)
globok = funok && solok && histok
if ~globok 
   println("global FD fails")
end
#
# Global convergence with parab3p
#
sdatap3p=nsolsc(atan,30.0; rtol=1.e-10,maxit=11)
solok=(abs(sdatap3p.solution) < 1.e-8)
funok=(abs(sdatap3p.functionval) < 1.e-8)
hs=size(sdatap3p.history)
histok=(hs[1]==12)
p3pok = funok && solok && histok
if ~p3pok
   println("parab3p fails")
end

#
# Local convergence with secant method
#
sdatas=nsolsc(atan,1.0; solver="secant")
solok=(abs(sdatas.solution) < 1.e-10)
funok=(abs(sdatas.functionval) < 1.e-10)
hs=size(sdatas.history)
histok=(hs[1]==6)
secantok = funok && solok && histok
if ~secantok
   println("secant failure")
end
#
# Initialize secant method when x0=0
#
zedata=nsolsc(fcos,0.0;solver="secant",rtol=1.e-9)
solution=7.390851333858823e-01
solok=(abs(zedata.solution-solution) < 1.e-9)
funok=(abs(zedata.functionval) < 1.e-9)
hs=size(zedata.history)
histok=(hs[1]==7)
zecok = funok && solok && histok
if ~zecok
   println("local FD fixup at zero fails")
end
#
# Tricky line search problem
#
sdatal=nsolsc(linatan,200.0; sham=5,maxit=20,armmax=10,armfix=true, 
              rtol=1.e-10)
solution=-100.0
solok=(abs(sdatal.solution-solution) < 1.e-8)
funok=(abs(sdatal.functionval) < 1.e-8)
hs=size(sdatal.history)
histok=(hs[1]==6)
shamfastok = funok && solok && histok
if ~shamfastok
   println("Fast Shamanskii response FAILURE")
end
#
# Test linesearch failure complaints.
#
armfail=nsolsc(atan,10.0; armmax=1, armfix=true, printerr=false)
afok=false
if armfail.idid==false && armfail.errcode == 1
   afok=true
else
   println("Armijo failure test FAILED.")
end
#
# Test residual failure mode and no history.
#
resok=false
resfail=nsolsc(atan,10.0; maxit=3, armfix=true, keepsolhist=false)
if resfail.idid==false && resfail.errcode == 10
   resok=true
else
   println("Residual failure test FAILED.")
end
#
# Test stagnation mode
#
stagdatan=nsolsc(ftanx,4.5, ftanxp; rtol=1.e-17, atol=1.e-17, 
         armfix=true, maxit=14)
fvals=stagdatan.history
avals=stagdatan.stats.iarm
ifvals=stagdatan.stats.ifun
jvals=stagdatan.stats.ijac
stagl=(length(fvals)==6)
stagf=(fvals[5] < 1.e-15)
stags=(avals[6]==5) && (ifvals[6]==6) &&  (jvals[6]==1)
stagok=stagl && stags && stagf
if ~stagok
   println("Stagnation test FAILED")
end
#
#
# Test chord method
#
lttest=nsolsc(atan,.5;solver="chord");
fvals=lttest.history
chordl=(length(fvals)==11)
ratl=fvals[11]/fvals[10]
chordr=(abs(ratl-.25) < 1.e-7)
solok=(fvals[11] < 1.e-6)
chordok=chordl && chordr && solok
if ~chordok
   println("Chord test failed")
end
#
# Make sure nsolsc knows if the solution and the intial iterate are
# the same
#
lotout=nsolsc(flot,0.0)
lotok= (lotout.errcode == -1)
if ~lotok
   println("Lottery test failed")
end
return locok && globok && secantok && analyticok && zecok && 
       shamfastok && afok && resok && p3pok && chordok && lotok
end

function flot(x)
flot = x*exp(x)
end
