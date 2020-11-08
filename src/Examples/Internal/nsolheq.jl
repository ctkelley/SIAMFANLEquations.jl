"""
nsolheq(x0, FS, FPS, hdata; diff=:fd)
Internal function to run with CI. Nothing to see here, move along.
"""
function nsolheq(x0, FS, FPS, hdata; diff = :fd)
    if diff == :fd
        heqout = nsol(heqf!, x0, FS, FPS; pdata = hdata, sham = 1)
    else
        heqout = nsol(heqf!, x0, FS, FPS, heqJ!; pdata = hdata, sham = 1)
    end
    return heqout
end
