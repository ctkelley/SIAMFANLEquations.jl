"""
Orthogonalize!(V, hv, vv, orth)

Orthogonalize the Krylov vectors using your (my) choice of
methods. Anything other than classical Gram-Schmidt twice (cgs2) is
likely to become an undocumented option. Methods other than cgs2 are
for CI for the linear solver.
"""
function Orthogonalize!(V, hv, vv, orth = "mgs1")
    orthopts = ["mgs1", "mgs2", "cgs1", "cgs2"]
    orth in orthopts || error("Impossible orth spec in Orthogonalize!")
    if orth == "mgs1"
        mgs!(V, hv, vv)
    elseif orth == "mgs2"
        mgs!(V, hv, vv, "twice")
    elseif orth == "cgs1"
        cgs!(V, hv, vv, "once")
    else
        cgs!(V, hv, vv, "twice")
    end
end

"""
cgs!(V, hv, vv, orth="twice")

Classical Gram-Schmidt.
"""
function cgs!(V, hv, vv, orth)
    T = eltype(V)
    k = length(hv)
    @views rk = hv[1:k-1]
    onep = T(1.0)
    zerop = T(0.0)
    pk = zeros(T, k - 1)
    qk = vv
    Qkm = V
    #
    # As the inventors of the BLAS said, "RTFM"!
    #
    # Orthogonalize!
    BLAS.gemv!('T', onep, Qkm, qk, onep, rk)
    BLAS.gemv!('N', -onep, Qkm, rk, onep, qk)
    if orth == "twice"
        #
        # Orthogonalize again!
        BLAS.gemv!('T', onep, Qkm, qk, zerop, pk)
        BLAS.gemv!('N', -onep, Qkm, pk, onep, qk)
        rk .+= pk
    end
    #
    # Keep track of what you did.
    nqk = norm(qk)
    qk ./= nqk
    hv[k] = nqk
end


#
function cgs!(V::SubArray{Float16,2}, hv, vv, orth)
    #
    #   no BLAS
    #
    k = length(hv)
    T = eltype(V)
    onep = T(1.0)
    zerop = T(0.0)
    @views rk = hv[1:k-1]
    pk = zeros(T, size(rk))
    qk = vv
    Qkm = V
    # Orthogonalize
    rk .+= Qkm' * qk
    qk .-= Qkm * rk
    if orth == "twice"
        # Orthogonalize again
        pk .= Qkm' * qk
        qk .-= Qkm * pk
        rk .+= pk
    end
    # Keep track of what you did.
    nqk = norm(qk)
    nqk != 0.0 || println("breakdown")
    nqk == 0.0 || qk ./= nqk
    hv[k] = nqk
end


"""
mgs!(V, hv, vv, orth)
"""
function mgs!(V, hv, vv, orth = "once")
    k = length(hv) - 1
    normin = norm(vv)
    #p=copy(vv)
    @views for j = 1:k
        p = vec(V[:, j])
        hv[j] = p' * vv
        vv .-= hv[j] * p
    end
    hv[k+1] = norm(vv)
    if (normin + 0.001 * hv[k+1] == normin) #&& (orth == "twice")
        @views for j = 1:k
            p = vec(V[:, j])
            hr = p' * vv
            hv[j] += hr
            vv .-= hr * p
        end
        hv[k+1] = norm(vv)
    end
    nv = hv[k+1]
    #
    # Watch out for happy breakdown
    #
    #if hv[k+1] != 0
    #@views vv .= vv/hv[k+1]
    nv != 0 || println("breakdown")
    if nv != 0
        vv ./= nv
    end
end
