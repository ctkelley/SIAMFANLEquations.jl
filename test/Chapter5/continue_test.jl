function continue_test()
v1ok = test_v1()
v2ok = test_v2()
continueok= v1ok && v2ok
end

function test_v1()
n=100; version= 1
(pval, nval, x, lambda) = heq_continue(n; version=version);
dpath = path_test.(pval,nval);
del1 = dpath[1:end-1];
del2 = dpath[end];
v1_pass = (norm(del1, Inf) < 4.e-9) && (del2 < 1.5e-4)
end

function test_PAC()
n=100; version = 3;
(pval, nval, x, lambda) = heq_continue(n; version=version);
dpath = path_test.(pval,nval);
nsingular=argmax(dpath);
del2=dpath[nsingular];
del1=[dpath[1:nsingular-1]; dpath[nsingular+1:end]];
v2_pass = (norm(del1,Inf) < 1.e-5) && (del2 < 1.e-4)
end

function path_test(pval,nval)
if pval > 0
rp=(1.0 + sqrt.(1.0 .- pval))/(.5 .* pval)
rm=(1.0 - sqrt.(1.0 .- pval))/(.5 .* pval)
else
rm = 1.0
end
(nval .<= 2) ? dp=abs.(nval-rm) : dp=abs.(nval-rp)
return dp
end

