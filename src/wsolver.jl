struct WSolver
    solver::WSymbol
end

macro WS_str(s::AbstractString)
    WSolver(WSymbol(s))
end

function wjaggedlist(x::AbstractArray{WElem,N}) where {N}
    if N == 1
        WExpr(W"List", map(y -> y.val, x))
    else
        s = size(x)
        x = reshape(x, s[1], :)
        WExpr(W"List", wjaggedlist.(reshape(x[i, :], s[2:N]) for i in 1:s[1]))
    end
end
function parsewlist(x::WExpr, reccall::Bool = false)
    iswlist(x) || error("")
    if all(y -> iswlist(y), x.args)
        x = parsewlist.(x.args, true)
        x = cat(x...; dims = ndims(first(x)) + 1)
        reccall ? x : permutedims(x, ndims(x):-1:1)
    else
        WElem.(x.args)
    end
end

function (s::WSolver)(a::AbstractArray{WElem}, rissingle::Bool = false)
    r = weval((s.solver)(wjaggedlist(a)))
    if rissingle
        iswlist(r) ? parsewlist(r) : WElem(r)
    else
        iswlist(r) || error("")
        map(x -> iswlist(x) ? parsewlist(x) : WElem(x), r.args)
    end
end
