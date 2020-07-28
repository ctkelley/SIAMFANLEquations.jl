function InitStats(resid)
history=[resid]
ifun=[1]
ijac=[0]
iarm=[0]
history=[resid]
NewStats=ItStats(ifun, ijac, iarm, history)
return NewStats
end

function updateStats!(ItData, newf, newj, newiarm, resid)
append!(ItData.ifun,newf)
append!(ItData.ijac,newj)
append!(ItData.iarm,newiarm)
append!(ItData.history,resid)
end

