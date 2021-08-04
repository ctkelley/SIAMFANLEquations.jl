#
# Get into the MGS reorthogonalization loop and see if it
# does its job.
#
function mgs_test(cond=1.e6)
(A, x0, b)=data_cook(cond)
V=zeros(3,20)
gout=kl_gmres(x0, b, matvec, V, 1.e-9; orth="mgs1",pdata=A)
gout2=kl_gmres(x0, b, matvec, V, 1.e-9; orth="mgs2",pdata=A)
del=gout.reshist-gout2.reshist
mgs2ok = (norm(del,Inf) > 1.e-12) && gout.idid && gout2.idid
mgs2ok || println("mgs_test fails")
return mgs2ok
#return(gout, gout2, mgs2ok)
end

function matvec(x,A)
return A*x
end

function data_cook(cond)
u1=[1,-2, 0]/sqrt(5.0)
u2=[0,0,1]
u3=[2,1,0]/sqrt(5.0)
U=[u1 u2 u3]
V=[u3 u1 u2]
D=diagm([1, cond, sqrt(cond)])
A=U * D * V'
xstar=ones(3,)
b=A*xstar
x0=[10.0,10.0,10.0]
return(A, x0, b)
end

