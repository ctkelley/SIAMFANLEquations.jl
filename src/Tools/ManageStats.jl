function InitStats(x, fc, resid, keepsolhist)
solution=x
functionval=fc
history=[resid]
ifun=[1]
ijac=[0]
iarm=[0]
idid=true
solhist=[]
if keepsolhist
   solhist=[x]
end
NewStats=ItStats(solution, functionval, history, 
                 ifun, ijac, iarm, idid, solhist)
return NewStats
end
