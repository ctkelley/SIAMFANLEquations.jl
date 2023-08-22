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
    beamtdout || println("error in beam_test beamtdout")
    (pout, nout) = ptcBeam(10, 100)
    pout2 = ptcBeam(10, 100; jknowsdt = true)
    kdtdiff =
        norm(pout.solution - pout2.solution, Inf) + norm(pout2.history - pout.history, Inf)
    kdtok = (kdtdiff < 1.e-14)
    kdtok || println("error is knowsdt test")
    nsolp = norm(pout.solution)
    nsoln = norm(nout.solution)
    itp = length(pout.history)
    pnormok = (nsolp > 5.0) && (nsoln < 1.e-15)
    pnormok || println("error in beam_test pnromok")
    presok = (itp < 100) && (pout.history[itp] < 1.e-10)
    presok || println("error in beam_test presok")
    return beamtdout && pnormok && presok && kdtok
end
