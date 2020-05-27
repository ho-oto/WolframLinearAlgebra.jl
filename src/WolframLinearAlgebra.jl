module WolframLinearAlgebra

import Base: promote_rule, zero, one, +, -, *, /, ^, sqrt, conj, adjoint, copy
import Base: sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
import Base: log, log2, log10, exp

using LinearAlgebra
import LinearAlgebra: Eigen, eigen, SVD, svd, Schur, schur
import LinearAlgebra: qr, lq
import LinearAlgebra: pinv, nullspace

import MathLink: WSymbol, WExpr, WInteger, WReal, weval, parseexpr, @W_str, @W_cmd

export WElem, WSolver, @WE_str, @WE_cmd, @WS_str, @W_str, @W_cmd
export weval, simplify, n
export eigen, svd, qr, lq, lu, schur, pinv, nullspace

iswlist(x) = false
iswlist(x::WExpr) = x.head == W"List"

struct WElem
    val::T where {T<:Union{Integer,AbstractFloat,WSymbol,WExpr,WInteger,WReal}}
    WElem(val) = iswlist(val) ? error("") : new(val)
end

promote_rule(::Type{WElem}, ::Type{T}) where {T<:Number} = WElem
copy(x::WElem) = WElem(deepcopy(x.val))
weval(x::WElem) = WElem(weval(x.val))

macro WE_str(s::AbstractString)
    WElem(WSymbol(s))
end
macro WE_cmd(c::Cmd)
    WElem(parseexpr(c))
end

function (s::WElem)(args...)
    s.val isa WSymbol || error("")
    all(x -> x isa Union{WElem,Number,WSymbol,WExpr,WInteger,WReal}, args) || error("")
    args = map(x -> x isa WElem ? x : WElem(x), args)
    WElem((s.val)(map(x -> x.val, args)...))
end

(+)(a::WElem, b::WElem) = WE"Plus"(a, b)
(-)(a::WElem, b::WElem) = WE"Subtract"(a, b)
(*)(a::WElem, b::WElem) = WE"Times"(a, b)
(/)(a::WElem, b::WElem) = WE"Divide"(a, b)
(^)(a::WElem, b::WElem) = WE"Power"(a, b)

(+)(a::WElem, b::Number) = a + WElem(b)
(-)(a::WElem, b::Number) = a - WElem(b)
(*)(a::WElem, b::Number) = a * WElem(b)
(/)(a::WElem, b::Number) = a / WElem(b)
(^)(a::WElem, b::Number) = a^WElem(b)

(+)(a::Number, b::WElem) = WElem(a) + b
(-)(a::Number, b::WElem) = WElem(a) - b
(*)(a::Number, b::WElem) = WElem(a) * b
(/)(a::Number, b::WElem) = WElem(a) / b
(^)(a::Number, b::WElem) = WElem(a)^b

WElem(val::Irrational{:π}) = WE"Pi"
WElem(val::Irrational{:ℯ}) = WE"E"
WElem(val::Irrational{:φ}) = WE"GoldenRatio"
WElem(val::Irrational{:γ}) = WE"EulerGamma"
WElem(val::Irrational{:catalan}) = WE"Catalan"
WElem(val::Complex) = WElem(real(val)) + WElem(imag(val)) * WE"I"
WElem(val::Rational) = WElem(numerator(val)) / WElem(denominator(val))

zero(::Type{WElem}) = WElem(0)
zero(x::WElem) = WElem(0)
one(::Type{WElem}) = WElem(1)
one(x::WElem) = WElem(1)

simplify(x::WElem) = weval(WE"Simplify"(x)) # TODO: keyword arguments
n(x::WElem) = weval(WE"N"(x)) # TODO: keyword arguments

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

include("basicfunctions.jl")
include("linargfunctions.jl")

end # module
