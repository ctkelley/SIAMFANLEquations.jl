function PitchFork1()
lambdam=collect(-1:.01:0)
lambdap=-lambdam
zm=zeros(size(lambdam))
stp=sqrt.(lambdap)
stm=-stp
plot(lambdam,zm,"k-")
plot(lambdap,zm,"k--")
plot(lambdap,stp,"k-")
plot(lambdap,stm,"k-")
xlabel(L"$\lambda$")
ylabel(L"$x^*(\lambda)$")
legend(["Stable","Unstable"])
end
