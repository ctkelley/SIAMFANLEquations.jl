mutable struct ItStats{T<:Real}
    ifun::Array{Int64,1}
    ijac::Array{Int64,1}
    iarm::Array{Int64,1}
    history::Array{T,1}
end

mutable struct ItStatsK{T<:Real}
    ifun::Array{Int64,1}
    ijac::Array{Int64,1}
    iarm::Array{Int64,1}
    ikfail::Array{Int64,1}
    history::Array{T,1}
end

mutable struct ItStatsPTC{T<:Real}
    history::Array{T,1}
end

#
# initfun = 1 unless it's the scalar secant method
#             then it's 2
#
function ItStats(hist, initfun = 1)
    ItStats([initfun], [0], [0], [hist])
end

function updateStats!(ItData::ItStats, newfun, newjac, AOUT)
    newiarm = AOUT.aiarm
    newfun = newfun + newiarm + 1
    resnorm = AOUT.resnorm
    append!(ItData.ifun, newfun)
    append!(ItData.ijac, newjac)
    append!(ItData.iarm, newiarm)
    append!(ItData.history, resnorm)
end

#
# initfun = 1 unless it's the scalar secant method
#             then it's 2
#
function ItStatsK(hist, initfun = 1)
    ItStatsK([initfun], [0], [0], [0], [hist])
end

function updateStats!(ItData::ItStatsK, newfun, newjac, AOUT, newikfail)
    newiarm = AOUT.aiarm
    newfun = newfun + newiarm + 1
    resnorm = AOUT.resnorm
    append!(ItData.ifun, newfun)
    append!(ItData.ijac, newjac)
    append!(ItData.iarm, newiarm)
    append!(ItData.ikfail, newikfail)
    append!(ItData.history, resnorm)
end

#
# Keep stats for PTC
#
function ItStatsPTC(hist)
    ItStatsPTC([hist])
end

function updateStats!(ItData::ItStatsPTC, resnorm)
    append!(ItData.history, resnorm)
end


