module WolframLinearAlgebra

import Base: promote_rule, zero, one, +, -, *, /, ^, sqrt, conj, adjoint
import Base: sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
import Base: log, log2, log10, exp

import MathLink: WSymbol, WExpr, WInteger, WReal, weval, parseexpr, @W_str

export WElem, WSolver, @WE_str, @WE_cmd, @WS_str
export weval, simplify
export weigen, wsvd, wqr

_iswlist(x) = false
_iswlist(x::WExpr) = x.head == W"List"
struct WElem
    val::T where {T<:Union{Integer,AbstractFloat,WSymbol,WExpr,WInteger,WReal}}
    WElem(val) = _iswlist(val) ? error("") : new(val)
end
promote_rule(::Type{WElem}, ::Type{T}) where {T<:Number} = WElem
macro WE_str(s::AbstractString)
    WElem(WSymbol(s))
end
macro WE_cmd(c::Cmd)
    WElem(parseexpr(c))
end

function (s::WElem)(args...)
    s.val isa WSymbol || error("")
    all(x -> x isa Union{WElem,Number,WSymbol,WExpr,WInteger,WReal}, args) || error("")
    any(x -> x isa AbstractIrrational, args) && error("Irrational is unsupported")
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
(^)(a::WElem, b::Number) = a ^ WElem(b)

(+)(a::Number, b::WElem) = WElem(a) + b
(-)(a::Number, b::WElem) = WElem(a) - b
(*)(a::Number, b::WElem) = WElem(a) * b
(/)(a::Number, b::WElem) = WElem(a) / b
(^)(a::Number, b::WElem) = WElem(a) ^ b

WElem(val::Complex) = WElem(real(val)) + WElem(imag(val)) * WE"I"
WElem(val::Rational) = WElem(numerator(val)) / WElem(denominator(val))

zero(::Type{WElem}) = WElem(0)
one(::Type{WElem}) = WElem(1)

weval(x::WElem) = WElem(weval(x.val))
simplify(x::WElem) = weval(WE"Simplify"(x))

struct WSolver
    solver::WSymbol
end
macro WS_str(s::AbstractString)
    WSolver(WSymbol(s))
end

function _wjaggedarray(x::AbstractArray{WElem,N}) where {N}
    if N == 1
        map(y -> y.val, x)
    else
        s = size(x)
        x = reshape(x, s[1], :)
        _wjaggedarray.(reshape(x[i, :], s[2:N]) for i in 1:s[1])
    end
end
function _parsearray(x::WExpr)
    x.head == W"List" || error("")
    if all(y -> y isa WExpr && y.head == W"List", x.args)
        x = _parsearray.(x.args)
        cat(x...; dims = ndims(first(x)) + 1)
    else
        WElem.(x.args)
    end
end
function (s::WSolver)(a::AbstractMatrix{WElem}, rnum::Int = 1)
    r = weval((s.solver)(_wjaggedarray(a)))
    if rnum == 1
        _iswlist(r) ? _parsearray(r) : r
    else
        r.head == W"List" || error("")
        r = r.args
        map(x -> _iswlist(x) ? _parsearray(x) : x, r)
    end
end

# -----

sqrt(a::WElem) = WE"Sqrt"(a)

sin(a::WElem) = WE"Sin"(a)
cos(a::WElem) = WE"Cos"(a)
tan(a::WElem) = WE"Tan"(a)

sinh(a::WElem) = WE"Sinh"(a)
cosh(a::WElem) = WE"Cosh"(a)
tanh(a::WElem) = WE"Tanh"(a)

asin(a::WElem) = WE"ArcSin"(a)
acos(a::WElem) = WE"ArcCos"(a)
atan(a::WElem) = WE"ArcTan"(a)

asinh(a::WElem) = WE"ArcSinh"(a)
acosh(a::WElem) = WE"ArcCosh"(a)
atanh(a::WElem) = WE"ArcTanh"(a)

log(a::WElem) = WE"Log"(a)
log2(a::WElem) = WE"Log2"(a)
log10(a::WElem) = WE"Log10"(a)

exp(a::WElem) = WE"Exp"(a)

conj(a::WElem) = WE"Conjugate"(a)
adjoint(a::WElem) = WE"Conjugate"(a)

# -----

weigen(x::AbstractMatrix{WElem}) = WS"Eigensystem"(x, 2)
wsvd(x::AbstractMatrix{WElem}) = WS"SingularValueDecomposition"(x, 3)
wqr(x::AbstractMatrix{WElem}) = WS"QRDecomposition"(x, 2)

end # module
