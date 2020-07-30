Base.@kwdef mutable struct ItStats{T <: Real}
ifun::Array{Int64,1} = [1]
ijac::Array{Int64,1} = [0]
iarm::Array{Int64,1} = [0]
history::Array{T,1}  
end

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

function updateStats!(ItData, newfun, newjac, AOUT)
newjac = newjac + AOUT.newjac
newiarm = AOUT.aiarm
newfun = newfun + newiarm + 1
resid=AOUT.resid
append!(ItData.ifun,newfun)
append!(ItData.ijac,newjac)
append!(ItData.iarm,newiarm)
append!(ItData.history,resid)
end
