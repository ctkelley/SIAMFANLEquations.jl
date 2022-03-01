"""
Orthogonalize!(V, hv, vv, orth; verbose=false)

Orthogonalize the Krylov vectors using your (my) choice of
methods. Anything other than classical Gram-Schmidt twice (cgs2) is
likely to become an undocumented and UNSUPPORTED option. Methods other 
than cgs2 are mostly for CI for the linear solver.

DO NOT use anything other than "cgs2" with Anderson acceleration.
"""
function Orthogonalize!(V, hv, vv, orth = "cgs2"; verbose = false)
    orthopts = ["mgs1", "mgs2", "cgs1", "cgs2"]
    orth in orthopts || error("Impossible orth spec in Orthogonalize!")
    if orth == "mgs1"
        mgs!(V, hv, vv; verbose = verbose)
    elseif orth == "mgs2"
        mgs!(V, hv, vv, "twice"; verbose = verbose)
    elseif orth == "cgs1"
        cgs!(V, hv, vv, "once"; verbose = verbose)
    else
        cgs!(V, hv, vv, "twice"; verbose = verbose)
    end
end

"""
mgs!(V, hv, vv, orth; verbose=false)
"""
function mgs!(V, hv, vv, orth = "once"; verbose = false)
    k = length(hv) - 1
    normin = norm(vv)
    #p=copy(vv)
    @views for j = 1:k
        p = vec(V[:, j])
        hv[j] = p' * vv
        vv .-= hv[j] * p
    end
    hv[k+1] = norm(vv)
    if (normin + 0.001 * hv[k+1] == normin) && (orth == "twice")
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
    (nv != 0) || (verbose && (println("breakdown in mgs1")))
    if nv != 0
        vv ./= nv
    end
end

"""
cgs!(V, hv, vv, orth="twice"; verbose=false)

Classical Gram-Schmidt.
"""
function cgs!(V, hv, vv, orth = "twice"; verbose = false)
    #
    #   no BLAS. Seems faster than BLAS since 1.6 and allocates
    #   far less memory.
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
    # New low allocation stuff
    mul!(rk, Qkm', qk, 1.0, 1.0)
    ###    mul!(pk, Qkm', qk)
    ###    rk .+= pk
    ##    rk .+= Qkm' * qk
    #    qk .-= Qkm * rk
    mul!(qk, Qkm, rk, -1.0, 1.0)
    if orth == "twice"
        # Orthogonalize again
        # New low allocation stuff
        mul!(pk, Qkm', qk)
        ##        pk .= Qkm' * qk
        #        qk .-= Qkm * pk
        mul!(qk, Qkm, pk, -1.0, 1.0)
        rk .+= pk
    end
    # Keep track of what you did.
    nqk = norm(qk)
    (nqk != 0) || (verbose && (println("breakdown in cgs")))
    (nqk > 0.0) && (qk ./= nqk)
    hv[k] = nqk
end
