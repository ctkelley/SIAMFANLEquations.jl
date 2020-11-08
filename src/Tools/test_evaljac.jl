function test_evaljac(ItRules, itc, newiarm, residratio)
    solver = ItRules.solver
    sham = ItRules.sham
    resdec = ItRules.resdec
    evaljacit = (itc % sham == 0 || newiarm > 0 || residratio > resdec)
    chordinit = (solver == "chord") && itc == 0
    evaljac = (evaljacit && solver == "newton") || chordinit || solver == "secant"
end
