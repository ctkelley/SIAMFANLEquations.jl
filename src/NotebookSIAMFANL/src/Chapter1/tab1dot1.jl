"""
tab1dot1()

This is the Julia code for the tan(x) = x example.

This code makes Table 1.1 in the print book.

See the notebook.

"""
function tab1dot1()
    ftout=nsolsc(ftanx,4.5; maxit=14, rtol=1.e-17, atol=1.e-17, 
          printerr=false);
    printhist(ftout.history[1:6],["|f(x)|"])
#return ftout
end

