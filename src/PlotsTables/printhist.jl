"""
printhist(tablein,headers;TeX=false,figures=5)

You are not supposed to look at this code. It is mostly repulsive
bookkeeping.

Inputs:
tablein: columns of equal length with the data
         Pad the short ones with NaNs. My codes which all this
         do that. The formatting turns the NaNs into spaces of
         the correct length.
headers: the titles for the columns
         example: headers=["foo", "bar"]
         Do not add the header for the iteration counter. I do that.
TeX:     Format the tabel in LaTeX or not?
figures: The printf format will be (figures+7).figures e
         So if figures = 5 the data will be formatted in %12e5.
         This makes everything line up.


Compare iteration histories from SIAMFANLEquations family of solvers.
This is usually called by something else.
Nothing to see here. Move along.
"""
function printhist(tablein,headers;TeX=false, figures=5)
ntab=length(tablein[1,:])
bighead=Array{String,2}(undef,1,ntab+1)
bighead[1]="n"
for ih=2:ntab+1
    bighead[ih]=headers[ih-1]
end
if ntab > 5
   error("Too many columns for the table. Use fewer.")
end
fmtout=buildformat(ntab, TeX, figures)
tabfmt=fmtout.format
nanspace=fmtout.nanspace
headerfmt=fmtout.headerfmt
printf(fmt::String,args...) = @eval @printf($fmt,$(args...))
sprintf(fmt::String,args...) = @eval @sprintf($fmt,$(args...))
itmax=length(tablein[:,1])
itc=0:itmax-1
if TeX
@printf("\\begin{tabular}{");
for i=1:ntab+1
    @printf("l");
end
@printf("} \n")
printf(headerfmt,bighead...)
else
printf(headerfmt,bighead...)
end
   for it=1:itmax
      st=sprintf(tabfmt,it-1,tablein[it,:]...)
      snan=findfirst(isequal('N'),st)
      lt=length(st)
      while typeof(snan) != Nothing
           st=string(st[1:snan-1],nanspace,st[snan+3:lt])
           snan=findfirst(isequal('N'),st)
      end
      printf("%s", st)
   end
if TeX
@printf("\\hline \n")
@printf("\\end{tabular} \n")
end
end

function buildformat(ntab, TeX, figures)
format="%2d "
nanspace="   "
fwide=figures+7
fm1=string(fwide)
fm2=string(figures)
basefmt=string(" %",fm1,".",fm2,"e ")
hfmt=string(" %",fm1,"s")
headerfmt=" %s"
if TeX
for i=1:ntab
    format=string(format,"&",basefmt)
    headerfmt=string(headerfmt,"&",hfmt)
end
format=string(format," \\\\ \n")
headerfmt=string(headerfmt," \\\\ \\hline \n")
else
for i=1:ntab
    format=string(format,basefmt)
    headerfmt=string(headerfmt,hfmt)
end
format=string(format," \n")
headerfmt=string(headerfmt," \n \n")
end
return (format=format, nanspace=nanspace, headerfmt=headerfmt)
end
