#
# The functions and data structures in this file organize the
# iteration statistics and report the results after the iteration
# is complete
#

#
# Keep the books for nsol and secant
#
mutable struct ItStats{T<:Real}
    ifun::Array{Int64,1}
    ijac::Array{Int64,1}
    iarm::Array{Int64,1}
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

function CollectStats(ItData::ItStats)
    stats = (ifun = ItData.ifun, ijac = ItData.ijac, iarm = ItData.iarm)
    return stats
end

#
# Keep the books for nsoli
#
mutable struct ItStatsK{T<:Real}
    ifun::Array{Int64,1}
    ijac::Array{Int64,1}
    iarm::Array{Int64,1}
    ikfail::Array{Int64,1}
    history::Array{T,1}
end

function ItStatsK(hist)
    ItStatsK([1], [0], [0], [0], [hist])
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

function CollectStats(ItData::ItStatsK)
    stats =
        (ifun = ItData.ifun, ijac = ItData.ijac, iarm = ItData.iarm, ikfail = ItData.ikfail)
    return stats
end

#
# Keep stats for PTC
#
mutable struct ItStatsPTC{T<:Real}
    history::Array{T,1}
end

function ItStatsPTC(hist)
    ItStatsPTC([hist])
end

function updateStats!(ItData::ItStatsPTC, resnorm)
    append!(ItData.history, resnorm)
end

function CollectStats(ItData::ItStatsPTC)
    stats = []
    return stats
end

"""
CloseIteration(x, FS, ItData, idid, errcode, keepsolhist, solhist = [])

Collect the solution, function value, iteration stats and send them back.
"""
function CloseIteration(x, FS, ItData, idid, errcode, keepsolhist, solhist = [])
    stats = CollectStats(ItData)
    ithist = ItData.history
    if keepsolhist
        sizehist = length(ithist)
        return (
            solution = x,
            functionval = FS,
            history = ithist,
            stats = stats,
            idid = idid,
            errcode = errcode,
            solhist = solhist[:, 1:sizehist],
        )
    else
        return (
            solution = x,
            functionval = FS,
            history = ithist,
            stats = stats,
            idid = idid,
            errcode = errcode,
        )
    end
end

