"""
ftanx()

Function and derivative for Figures 1.1 and 1.3 in the print book.
Also used as simple test problems.
"""
function ftanx(x)
   return tan(x) - x
end

function ftanxp(x)
   return sec(x)^2 - 1
end

