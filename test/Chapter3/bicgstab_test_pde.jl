"""
bicgstab_test_pde(n)

PDE test from FR16. Test of kl_bicgstab with all kinds of preconditioning.
"""
function bicgstab_test_pde(n; write = false, eta = 9.8 * 1.e-4)
    pdata = pdegminit(n)
    RHS = pdata.RHS
    ue = pdata.ue
    u0 = zeros(n * n)
#
    # Solve with left preconditioning hard-wired in
#
    fdata = pdata.fdata
    b = Pfish2d(RHS, fdata)
    goutp = kl_bicgstab(u0, b, pdelpatv, eta; pdata = pdata, lmaxit=200)
    pcres = goutp.reshist
    pcres /= pcres[1]
    sollhw = goutp.sol
#
    # Solve with right preconditioning hard-wired in
#
    goutrp = kl_bicgstab(u0, RHS, pderatv, eta; pdata = pdata, lmaxit=200)
    pcresr = goutrp.reshist
    pcresr /= pcresr[1]
    solrhw = copy(u0)
    solrhw .= pdeptv(goutrp.sol,pdata)
    solldiff=norm(solrhw-sollhw,Inf)
#
    # Solve with right preconditioning 
#
    goutrp1 = kl_bicgstab(u0, RHS, pdeatv, eta, pdeptv;
             pdata = pdata, lmaxit=200)
    pcresr1 = goutrp1.reshist
    pcresr1 /= pcresr1[1]
    solr = goutrp1.sol
#    solldiff += norm(solr-sollhw,Inf)
    solldiff = max(solldiff,norm(solr-sollhw,Inf))
#
    # Solve with left preconditioning 
#
    goutl1 = kl_bicgstab(u0, RHS, pdeatv, eta, pdeptv;
             pdata = pdata, lmaxit=200, side="left")
    pcresl1 = goutl1.reshist
    pcresl1 /= pcresl1[1]
    soll = goutl1.sol
#    solldiff += norm(soll-sollhw,Inf)
    solldiff = max(solldiff,norm(soll-sollhw,Inf))
#
# Hardwired and normal give same results?
#
    leftdel= norm(soll-sollhw,Inf) + norm(pcres-pcresl1,Inf)
    leftpass=(leftdel < 1.e-15)
    rightdel = norm(solr-solrhw,Inf) + norm(pcresr1-pcresr,Inf)
    rightpass=(rightdel < 1.e-15)
#
    # Solve with no preconditioning to duplicate fig 3.4 in red book
#
    goutnp = kl_bicgstab(u0, RHS, pdeatv, eta; lmaxit=200, pdata = pdata)
    pcresnp = goutnp.reshist
    pcresnp /= pcresnp[1]
    solnone=goutnp.sol
    solldiff = max(solldiff,norm(solnone-sollhw,Inf))
#
# Solve with no preconditioning to get a failure
#
    goutnf = kl_bicgstab(u0, RHS, pdeatv, V, eta; lmaxit=20, pdata = pdata)
    failpass = ~goutnf.idid && (goutnf.lits==20)
#
# Are the answers close enough?
#
    println(solldiff/eta)
    sollpass=(solldiff < 2.0*eta)
#
    pass=sollpass && rightpass && leftpass && failpass
    return pass
end
