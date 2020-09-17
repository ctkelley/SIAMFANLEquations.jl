function plothist(avals,labelarray; semilogflag=true, ptitle="")
   fmtplot=["k-","k--","k-.","k:"]
   nplot=length(labelarray)
if nplot > 4
   error("Too many lines for the graph. Use fewer.")
end
   figure(1)
   for ip=1:nplot
   itmax=length(avals[ip])
   itc=0:itmax-1
   if semilogflag 
      semilogy(itc,avals[ip],fmtplot[ip])
   else
      loglog(itc,avals[ip],fmtplot[ip])
   end
end
legend(labelarray)
if ptitle != ""
   title(ptitle)
end
ylabel("Nonlinear Residual Norm")
xlabel("Nonlinear Iterations")
end
