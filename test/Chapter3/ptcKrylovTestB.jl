function ptcKrylovTestB(n = 63)
    delta0 = 0.01
    pout1 = ptciBeam(; lsolver = "bicgstab")
    pout2 = ptciBeam(n, delta0, false; lsolver = "bicgstab")
    pout3 = ptciBeam(n, delta0, false, "left"; lsolver = "bicgstab")
    sol1 = pout1.solution
    sol2 = pout2.solution
    sol3 = pout3.solution
    #
    # sol3 is the wrong stable branch. Left preconditioning bites you!
    #
    solpass1a = (norm(sol1 - sol2, Inf) < 1.e-9)
    solpass1b = (norm(sol1 - sol3, Inf) < 1.e-9)
    solpass1 = solpass1a && solpass1b
    solpass1 || println("solpass1 fails for ptcsoli-bicgstab")
    histpass = (length(pout1.history) == 25)
    histpass || println("histpass fails for ptcsoli-bicgstab")
    solpass2 = (abs(norm(sol1, Inf) - 2.191) < 1.e-3)
    solpass2 || println("solpass2 fails for ptcsoli-bicgstab")
    ptcipass = solpass1 && histpass && solpass2

end
