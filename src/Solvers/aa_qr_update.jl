"""
aa_qr_update(Q, R, vnew, m, k, Qd)

Update the QR factorization for the Anderson acceleration optimization
problem.

Still need to make the allocation for Qd go away.
"""
function aa_qr_update!(Q, R, vnew, m, k, Qd)
(n,m)=size(Q)
    aaqr_dim_check(Q, R, vnew, m, k)
    if k == 0
        R[1, 1] = norm(vnew)
        @views Q[:, 1] .= vnew / norm(vnew)
    else
        if k > m - 1
            downdate_aaqr!(Q, R, m, Qd)
        end # inner if block
        kq = min(k, m - 1)
        update_aaqr!(Q, R, vnew, m, kq)
    end # outer if block
    return (Q, R)
end

function update_aaqr!(Q, R, vnew, m, k)
    (nq, mq) = size(Q)
    (k > m - 1) && error("Dimension error in Anderson QR")
    @views Qkm=Q[:,1:k]
    @views hv = vec(R[1:k+1, k+1])
    Orthogonalize!(Qkm, hv, vnew, "cgs2")
    @views R[1:k+1, k+1] .= hv
    @views Q[:, k+1] .= vnew
#    return (Q = Q, R = R)
end

function downdate_aaqr!(Q, R, m, Qd)
(nq,mq)=size(Q)
(pd,md)=size(Qd)
(md == m-1) || @error("dimension error in downdate")
    @views Rp = R[:, 2:m]
    G = qr!(Rp)
    Rd = Matrix(G.R)
    Qx = Matrix(G.Q)
    @views R[1:m-1, 1:m-1] .= Rd
    @views R[:,m].=0.0
if (pd==nq)
    mul!(Qd,Q,Qx)
    @views Q[:,1:m-1] .= Qd
else
    blocksize=pd
    (dlow,dhigh)=blockdim(nq,blocksize)
    blen=length(dlow)
    for il=1:blen
    asize=dhigh[il]-dlow[il]+1
    @views QZ=Qd[1:asize,:]
    @views Qsec=Q[dlow[il]:dhigh[il],:]
    @views mul!(QZ,Qsec,Qx)
    @views Qsec[:,1:m-1] .= QZ
    end
end
    @views Q[:, m] .= 0.0
    return (Q, R)
end

function aaqr_dim_check(Q, R, vnew, m, k)
    (mq, nq) = size(Q)
    (mr, nr) = size(R)
    n = length(vnew)
    dimqok = ((mq == n) && (nq == m))
    dimrok = ((mr == m) && (nr == m))
    dimok = (dimqok && dimrok)
    dimok || error("array size error in AA update")
end

function blockdim(n, block)
    p = Int(floor(n / block))
    res = n - p * block
    ilow = Int64[]
    ihigh = Int64[]
    for jb = 1:p
        lowval = (jb - 1) * block + 1
        push!(ilow, lowval)
        highval = ilow[jb] + block - 1
        push!(ihigh, highval)
    end
    if res > 0
        lowval = p * block + 1
        push!(ilow, lowval)
        push!(ihigh, n)
    end
    return (ilow, ihigh)
end
