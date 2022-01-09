#
# Test the transport solve with s=infty against the data 
# from Tables 1 and 2 of
#
# author="R.D.M. Garcia and C.E. Siewert",
# title = "Radiative transfer in finite inhomogeneous plane-parallel 
#          atmospheres",
# journal="J. Quant. Spectrosc. Radiat. Transfer",
# year = 1982,
# volume=27,
# pages="141--148"
# 
function transport_test()
nx=2^8; na2=40; s=Inf;
vleft=1.0
vright=0.0
sn_data=sn_init(nx, na2, x -> exp(-x/s), 5.0, vleft, vright)
source=zeros(nx)
#
tol=1.e-5
kout=find_flux(source, sn_data, tol)
#
(sn_left, sn_right) = sn_tabulate(s, nx, kout.sol, source; maketab=false)
(out_left, out_right) = ces_data();
diff=norm(out_left-sn_left,Inf) + norm(out_right-sn_right, Inf)
kynum=length(kout.reshist)
transok =  (diff < 1.e-4) && (kynum <= 13)
transok || println("Transport test fails: dataerr = $diff; itcount = $kynum")
return transok
end

function find_flux(source, sn_data, tol)
b=getrhs(source, sn_data)
kout=kl_gmres(sn_data.phi0,b,AxB,sn_data.V,tol; pdata=sn_data)
return kout
end


function ces_data()
out_left=[8.97797e-01, 8.87836e-01, 8.69581e-01, 8.52299e-01, 8.35503e-01, 
8.18996e-01, 8.02676e-01, 7.86493e-01, 7.70429e-01, 7.54496e-01, 7.38721e-01];
out_right=[1.02202e-01, 1.12164e-01, 1.30419e-01, 1.47701e-01, 1.64497e-01, 
1.81004e-01, 1.97324e-01, 2.13507e-01, 2.29571e-01, 2.45504e-01, 2.61279e-01];
return (out_left, out_right)
end


"""
sn_tabulate(s, nx, flux, psi_left, psi_right, source ; maketab=true)

Make the tables to compare with Garcia/Siewert

Uses the converged flux from the solve.
"""
function sn_tabulate(s, nx, flux, source; maketab=true)
    angleout = [-.05; collect(-.1:-.1:-1.0); 0.05; collect(0.1:0.1:1.0)]
    #
    # I don't really need the weights, but sn_init expects some
    weights = angleout
    #
    na2 = length(angleout)
    na = floor(Int, na2 / 2)
    vleft=1.0
    vright=0.0
    np=nx 
    tsn_data = sn_init(nx, na2, x -> exp(-x/s), 5.0, vleft, vright; siewert=true)
    psi_right=tsn_data.psi_right
    psi_left=tsn_data.psi_left
    psi = tsn_data.psi
    psi = transport_sweep!(psi, flux, psi_left, psi_right, source, tsn_data)
    if maketab
    header = " mu         I(0,-mu)        I(tau,mu)"
    @printf("%s \n", header)
    for it=1:na
    @printf("%5.2f %15.5e %15.5e \n", 
           angleout[it+na], psi[it,1], psi[na+it,np])
    end
    end
    return (left = psi[1:na, 1], right = psi[na+1:na2, np])
end
