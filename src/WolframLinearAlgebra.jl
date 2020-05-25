module WolframLinearAlgebra

import Base: zero, one, +, -, *, /, ^, sqrt, conj, adjoint
import Base: sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
import Base: log, log2, log10, exp

import MathLink: WSymbol, WExpr, WInteger, WReal, weval, parseexpr, @W_str

export WElements, WOperator, @WE_str, @WE_cmd, @WO_str
export weval, simplify
export weigen, wsvd, wqr

_iswlist(x) = false
_iswlist(x::WExpr) = x.head == W"List"
struct WElements
    val::T where {T<:Union{Integer,AbstractFloat,WSymbol,WExpr,WInteger,WReal}}
    WElements(val) = _iswlist(val) ? error("") : new(val)
end
macro WE_str(s::AbstractString)
    WElements(WSymbol(s))
end
macro WE_cmd(c::Cmd)
    WElements(parseexpr(c))
end

function (s::WElements)(args...)
    s.val isa WSymbol || error("")
    all(x -> x isa Union{WElements,Number,WSymbol,WExpr,WInteger,WReal}, args) || error("")
    any(x -> x isa AbstractIrrational) && error("")
    args = map(x -> x isa WElements ? x : WElements(x), args)
    WElements((s.val)(map(x -> x.val, args)...))
end
(+)(a::WElements, b::WElements) = WE"Plus"(a, b)
(-)(a::WElements, b::WElements) = WE"Subtract"(a, b)
(*)(a::WElements, b::WElements) = WE"Times"(a, b)
(/)(a::WElements, b::WElements) = WE"Divide"(a, b)
(^)(a::WElements, b::WElements) = WE"Power"(a, b)

WElements(val::Complex) = WElements(real(val)) + WElements(imag(val)) * WE"I"
WElements(val::Rational) = WElements(numerator(val)) / WElements(denominator(val))

zero(::Type{WElements}) = WElements(0)
one(::Type{WElements}) = WElements(1)

weval(x::WElements) = WElements(weval(x.val))
simplify(x::WElements) = weval(WE"Simplify"(x))

struct WOperator
    op::WSymbol
end
macro WO_str(s::AbstractString)
    WOperator(WSymbol(s))
end

function _wjaggedarray(x::AbstractArray{WElements,N}) where {N}
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
        WElements.(x.args)
    end
end
function (s::WOperator)(a::AbstractMatrix{WElements}, rnum::Int = 1)
    r = weval((s.op)(_wjaggedarray(a)))
    if rnum == 1
        _iswlist(r) ? _parsearray(r) : r
    else
        r.head == W"List" || error("")
        r = r.args
        map(x -> _iswlist(x) ? _parsearray(x) : x, r)
    end
end

# -----

sqrt(a::WElements) = WE"Sqrt"(a)

sin(a::WElements) = WE"Sin"(a)
cos(a::WElements) = WE"Cos"(a)
tan(a::WElements) = WE"Tan"(a)

sinh(a::WElements) = WE"Sinh"(a)
cosh(a::WElements) = WE"Cosh"(a)
tanh(a::WElements) = WE"Tanh"(a)

asin(a::WElements) = WE"ArcSin"(a)
acos(a::WElements) = WE"ArcCos"(a)
atan(a::WElements) = WE"ArcTan"(a)

asinh(a::WElements) = WE"ArcSinh"(a)
acosh(a::WElements) = WE"ArcCosh"(a)
atanh(a::WElements) = WE"ArcTanh"(a)

log(a::WElements) = WE"Log"(a)
log2(a::WElements) = WE"Log2"(a)
log10(a::WElements) = WE"Log10"(a)

exp(a::WElements) = WE"Exp"(a)

conj(a::WElements) = WE"Conjugate"(a)
adjoint(a::WElements) = WE"Conjugate"(a)

# -----

weigen(x::AbstractMatrix{WElements}) = WO"Eigensystem"(x, 2)
wsvd(x::AbstractMatrix{WElements}) = WO"SingularValueDecomposition"(x, 3)
wqr(x::AbstractMatrix{WElements}) = WO"QRDecomposition"(x, 2)

end # module
