"""
function beam_test()
Test the time-dependent and steady state beam problem.

"""
function beam_test()
    #
    dt = 0.02
    n = 20
    stepnum = 5
    (t, se, xe, fhist, fhistt) = ivpBeam(n, dt, stepnum)
    beamtdout = (length(fhist) == 6) && (norm(fhistt, Inf) < 5.e-5)
    (pout, nout) = ptcBeam(10, 100)
    nsolp = norm(pout.solution)
    nsoln = norm(nout.solution)
    itp = length(pout.history)
    pnormok = (nsolp > 5.0) && (nsoln < 1.e-15)
    presok = (itp < 100) && (pout.history[itp] < 1.e-10)
    return beamtdout && pnormok && presok
end
