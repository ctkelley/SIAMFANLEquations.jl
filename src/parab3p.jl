"""
parab3p(lambdac, lambdam, ff0, ffc, ffm)

Three point parabolic line search.

input:
       lambdac = current steplength
       lambdam = previous steplength
       ff0 = value of || F(x_c) ||^2
       ffc = value of || F(x_c + lambdac d) ||^2
       ffm = value of || F(x_c + lambdam d) ||^2

output:
       lambdap = new value of lambda

internal parameters:
       sigma0 = .1, sigma1=.5, safeguarding bounds for the linesearch

You get here if cutting the steplength in half doesn't get you
sufficient decrease. Now you have three points and can build a parabolic
model. I do not like cubic models because they either need four points
or a derivative. 

So let's think about how this works. I cheat a bit and check the model
for negative curvature, which I don't want to see.

 The polynomial is

 p(lambda) = ff0 + (c1 lambda + c2 lambda^2)/d1

 d1 = (lambdac - lambdam)*lambdac*lambdam < 0
 So if c2 > 0 we have negative curvature and default to
      lambdap = sigma0 * lambda
 The logic is that negative curvature is telling us that
 the polynomial model is not helping much, so it looks better
 to take the smallest possible step. This is not what I did in the
 matlab code because I did it wrong. I have sinced fixed it.

 So (Students, listen up!) if c2 < 0 then all we gotta do is minimize
 c1 lambda + c2 lambda^2 over [.1* lambdac, .5*lambdac]
 So I find the zero of the derivative and check the endpoints.

"""

function parab3p(lambdac, lambdam, ff0, ffc, ffm)
#
# internal parameters
#
sigma0=.1; sigma1=.5;
#
c2 = lambdam*(ffc-ff0)-lambdac*(ffm-ff0)
if c2 >=0
#
# Sanity check for negative curvature
#
   lambdap = sigma0*lambdac
else
#
# It's a convex parabola, so use calculus!
#
   c1=lambdac*lambdac*(ffm-ff0)-lambdam*lambdam*(ffc-ff0);
   lambdap=-c1*.5/c2
#
   lambdaup=sigma1*lambdac
   lambdadown=sigma0*lambdac
   lambdap = max(lambdadown,min(lambdaup,lambdap))
end
return lambdap
end

