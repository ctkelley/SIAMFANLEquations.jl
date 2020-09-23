"""
nsolheq(x0, FS, FPS, hdata; diff=:fd)
Internal function to run with CI. Nothing to see here, move along.
"""
function nsolheq(x0, FS, FPS, hdata; diff=:fd)
if diff == :fd
heqout=nsold(heqf!, x0, FS, FPS; pdata = hdata);
else
heqout=nsold(heqf!, x0, FS, FPS, heqJ!; pdata = hdata);
end
return heqout
end
