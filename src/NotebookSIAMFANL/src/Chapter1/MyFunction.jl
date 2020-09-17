module MyFunction

export(myfun)
export(myfunp)

function myfun(x)
    fun=cos(x)-x
    return fun
end

function myfunp(x)
    funp=-sin(x)-1
    return funp
end

end
