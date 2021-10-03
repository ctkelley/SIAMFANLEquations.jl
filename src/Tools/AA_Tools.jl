#
# Keep the books for aasol
#
mutable struct ItStatsA{T<:Real}
    condhist::Array{T,1}
    alphanorm::Array{T,1}
    history::Array{T,1}
end

function ItStatsA(rnorm)
    ItStatsA([1.0],[1.0],[rnorm])
end

#
# Collect the stats at the end of the iteration
# 
function CollectStats(ItData::ItStatsA)
    stats = (condhist=ItData.condhist[2:end], 
                    alphanorm=ItData.alphanorm[2:end])
return stats
end

function updateStats!(ItData::ItStatsA, condhist, alphanorm)
      append!(ItData.condhist, condhist)
      append!(ItData.alphanorm, alphanorm)
end

function updateHist!(ItData::ItStatsA, rnorm)
      append!(ItData.history, rnorm)
end

#
# Initialize Anderson iteration
#
function Anderson_Init(x0, Vstore, m, maxit, beta, keepsolhist)
blocksize=1024
(0.0 < abs(beta) <= 1) || error("abs(beta) must be in (0,1]")
sol=copy(x0)
n=length(x0)
(mv, nv) = size(Vstore)
mv == n || error("Vstore needs ", n, " rows")
(nv >= 2 * (m + 1)) || error("Vstore needs ", 2 * m + 4, " columns")
#
# Just in case you are reusing Vstore for several problems, I will
# reinitialize it to zero.
#
Vstore .= 0.0
if m==0 
   Qd=[]; QP=[]; DG=[];
   nvblock=1
else
QP = @views Vstore[:, 1:m]
DG = @views Vstore[:, m+1:2*m]
if (nv >= 3*m+3) 
    (Qd = @views Vstore[:,2*m+1:3*m-1])
    nvblock=3*m
else
    @warn "Low storage mode"
    Qd=zeros(blocksize,m-1)
    nvblock=2*m+1
end
end
gx = Anderson_vector_Init(Vstore,nvblock)
df = Anderson_vector_Init(Vstore,nvblock+1)
dg = Anderson_vector_Init(Vstore,nvblock+2)
res = Anderson_vector_Init(Vstore,nvblock+3)
keepsolhist ? (solhist = solhistinit(n, maxit, sol)) : (solhist = [])
return(sol, gx, df, dg, res, DG, QP, Qd, solhist)
end

function Anderson_vector_Init(Vstore,nvblock)
gx = @views Vstore[:,nvblock]
return gx
end


#
# Figure out what idid and errcode are. Boring but must be done.
# This is a lot simpler than for Newton. There are no linesearches
# or Krylov iterations to keep track of.
#
function AndersonOK(resnorm, tol, k, toosoon)
idid=(resnorm <= tol)
idid ? (errcode = 0) : (errcode = 10)
toosoon && (errcode = -1)
idid || println("Failure to converge in aasol.jl")
toosoon && println("Iteration terminates on entry to aasol.jl")
return (idid, errcode)
end

"""
falpha(alpha,theta,mk)
Map thetas to alphas for stats
"""
function falpha(alpha, theta, mk)
    alpha[1] = theta[1]
    for ia = 2:mk
        alpha[ia] = theta[ia] - theta[ia-1]
    end
    alpha[mk+1] = 1.0 - theta[mk]
    return norm(alpha, 1)
end
