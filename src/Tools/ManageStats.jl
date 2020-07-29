mutable struct ItStats{T <: Real}
ifun::Array{Int64,1}
ijac::Array{Int64,1}
iarm::Array{Int64,1}
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

