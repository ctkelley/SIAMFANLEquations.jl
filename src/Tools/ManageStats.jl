mutable struct ItStats{T <: Real}
ifun::Array{Int64,1}
ijac::Array{Int64,1}
iarm::Array{Int64,1}
history::Array{T,1}
end

function ItStats(hist)
   ItStats([1],[0],[0],[hist])
end

function updateStats!(ItData, newfun, newjac, AOUT)
newiarm = AOUT.aiarm
newfun = newfun + newiarm + 1
resnorm=AOUT.resnorm
append!(ItData.ifun,newfun)
append!(ItData.ijac,newjac)
append!(ItData.iarm,newiarm)
append!(ItData.history,resnorm)
end
