function Newton_Krylov_Init(
    x0,
    dx,
    F!,
    Jvec,
    Pvec,
    lsolver,
    eta,
    fixedeta,
    armmax,
    armfix,
    maxit,
    lmaxit,
    printerr,
    pdata,
)
    #
    #   Initialize the iteration.
    #
    eta > 0 || error("eta must be positive")
    n = length(x0)
    x = copy(x0)
    ItRules = (
        dx = dx,
        f = F!,
        Jvec=Jvec,
        Pvec=Pvec,
        lsolver = lsolver,
        eta=eta,
        fixedeta=fixedeta,
        lmaxit = lmaxit,
        armmax = armmax,
        armfix = armfix,
        maxit = maxit,
        printerr = printerr,
        pdata = pdata,
    )
    return (ItRules, x, n)
end

function Krylov_Step!(step, x, FS, FPS, ItRules, etag)
#
# Test for too much, too soon.
#
lsolver=ItRules.lsolver
lmaxit=ItRules.lmaxit
T=eltype(FPS)
(nk,mk)=size(FPS)
n=length(step)
n == nk || error("Krylov vectors wrong length")
lmaxit < mk || error("Restarts not enabled yet")
lsolver == "gmres" || error(lsolver," ","not supported")
#T == Float64 || error("Mixed precision not supported yet. Must figure
#out RHS scaling.")
Jvec=ItRules.Jvec
Pvec=ItRules.Pvec
pdata=ItRules.pdata
dx=ItRules.dx
f=ItRules.f
fixedeta=ItRules.fixedeta
s0=zeros(size(step))
kdata=(pdata=pdata, dx=dx, xc=x, f=f, FS=FS, Jvec=Jvec, Pvec=Pvec)
#
# map the Jacobian-vector project from nsoli format to what
# kl_gmres wants to see
#
Jvecg = Jvec2
Jvec == dirder && (Jvecg=Jvec)
#
#RHS=FS
#T == Float64 || (RHS=T.(FS))
kout=kl_gmres(s0, FS, Jvecg, FPS, etag; pdata=kdata)
step = -kout.sol
reshist=kout.reshist
lits=kout.lits
idid=kout.idid
Lstats= (reshist=reshist, lits=lits, idid=idid)
return (step=step, Lstats=Lstats)
end

function Jvec2(v,kdata)
F=kdata.f
FS=kdata.FS
xc=kdata.xc
JV=kdata.Jvec
atv=JV(v, FS, xc)
return atv
end

function dirder(v,kdata)
pdata=kdata.pdata
dx=kdata.dx
F=kdata.f
FS=kdata.FS
xc=kdata.xc
delx=copy(xc)
delx.=xc + dx*v
FPP=copy(xc)
EvalF!(F, FPP, delx, pdata)
atv=(FPP-FS)/dx
return atv
end

function forcing(itc, residratio, etag, ItRules)
gamma=.9
etamax=ItRules.eta
fixedeta=ItRules.fixedeta
if fixedeta || (itc==0)
   etag=etamax
else
   safeguard=gamma*(etag^2)
   etaA=gamma*(residratio^2)
   if safeguard <= .1
      etag=min(etamax,etaA)
   else
      etag=min(etamax, max(etaA, safeguard))
   end 
end
   return etag
end


