function plothist(args...)
   fmtplot=("k-","k--","k-.","k-.","k>:")
   na=length(args)
if na > 10
   error("Too many lines for the graph. Use fewer.")
end
   nplot=na/2
   labelarray=()
   figure(1)
   for iq=1:nplot
   ip=2*iq-1
   labelarray=(labelarray, args[ip+1])
   itmax=length(args[ip].history)
   itc=0:itmax-1
   semilogy(itc,args[ip].history,fmtplot[ip])
end
legend(labelarray)
ylabel("Nonlinear Residual Norm")
xlabel("Nonlinear Iterations")
end
