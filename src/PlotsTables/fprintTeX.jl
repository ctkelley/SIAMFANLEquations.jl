"""
fprintTeX(headers,formats,data)

Print a LaTeX table from a Julia array.

Inputs:
headers: the titles for the columns
          example: headers=["foo", "bar"]
formats: c-style formatting for the columns.
          fprintTeX will add the carriage returns for you.
          example: formats="%d & %7.2e";
"""
function fprintTeX(headers, formats, data)
    (mr, mc) = size(data)
    @printf("\\begin{tabular}{")
    for i = 1:mc
        @printf("l")
    end
    @printf("} \n")
    for i = 1:mc-1
        @printf("%9s &", headers[i])
    end
    @printf("%9s \\\\ \n" , headers[mc])
    @printf("\\hline \n")
    #
    # I am not sure why @printf needs this, but it does.
    # See https://github.com/JuliaLang/julia/issues/4248
    #
    printf(fmt::String, args...) = @eval @printf($fmt, $(args...))
    #
    bigform = string(formats, "   \\\\ \n")
    for i = 1:mr
        printf(bigform, data[i, :]...)
    end
    @printf("\\hline \n")
    @printf("\\end{tabular} \n")
end
