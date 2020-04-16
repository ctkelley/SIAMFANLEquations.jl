"""
spitchfork(u,lambda)

The nonlinearity f(u) = u^3 - lamba u. The dynamics for du/du = -f(u)
have a pitchfork bifurcation at lambda=0. The steady-state solution 
u=0 is unique for lambda < 0 and there are three steady-state solutions
if lambda > 0. This is a simple-minded version of the buckling beam problem.

The function sptest(u) = spitchfork(u,.5) is the one I call in the testing.

"""
function spitchfork(u,lambda)
fu=u^3 - lambda*u
return fu
end

function sptest(u)
lambda=.5
return spitchfork(u,lambda)
end

function spitchp(u,lambda)
fp=3*u^2 - lambda
return fp
end

function sptestp(u)
lambda=.5
return spitchp(u,lambda)
end
