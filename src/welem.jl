struct WElem
    val::T where {T<:Union{Integer,AbstractFloat,WSymbol,WExpr,WInteger,WReal}}
    WElem(val) = iswlist(val) ? error("WExpr of List cannot become a WElem") : new(val)
end

function (s::WElem)(args...)
    s.val isa WSymbol || error("")
    all(x -> x isa Union{WElem,Number,WSymbol,WExpr,WInteger,WReal}, args) || error("")
    args = map(x -> x isa WElem ? x : WElem(x), args)
    WElem((s.val)(map(x -> x.val, args)...))
end

macro WE_str(s::AbstractString)
    WElem(WSymbol(s))
end

macro WE_cmd(c::Cmd)
    WElem(parseexpr(c))
end

WElem(val::Irrational{:π}) = WE"Pi"
WElem(val::Irrational{:ℯ}) = WE"E"
WElem(val::Irrational{:φ}) = WE"GoldenRatio"
WElem(val::Irrational{:γ}) = WE"EulerGamma"
WElem(val::Irrational{:catalan}) = WE"Catalan"
WElem(val::Complex) = WE"Complex"(real(val), imag(val))
WElem(val::Rational) = WE"Rational"(numerator(val), denominator(val))

promote_rule(::Type{WElem}, ::Type{T}) where {T<:Number} = WElem

convert(::Type{WElem}, x::Number) = WElem(x)
# isinteger(x::WElem)
zero(::Type{WElem}) = WElem(0)
zero(x::WElem) = WElem(0)
one(::Type{WElem}) = WElem(1)
one(x::WElem) = WElem(1)
# iszero(x::WElem)
# isone(x::WElem)
size(x::WElem) = ()
size(x::WElem, d::Integer) = d < 1 ? throw(BoundsError()) : 1
axes(x::WElem) = ()
axes(x::WElem, d::Integer) = d < 1 ? throw(BoundsError()) : OneTo(1)
# eltype(::Type{T}) where {T<:Number} = T
ndims(x::WElem) = 0
ndims(::Type{WElem}) = 0
length(x::WElem) = 1
firstindex(x::WElem) = 1
lastindex(x::WElem) = 1
IteratorSize(::Type{WElem}) = HasShape{0}()
keys(::WElem) = OneTo(1)
first(x::WElem) = x
last(x::WElem) = x
copy(x::WElem) = WElem(deepcopy(x.val))

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

inv(a::WElem) = 1 / a
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
transpose(a::AbstractMatrix{WElem}) = permutedims(a, (2,1))
adjoint(a::AbstractMatrix{WElem}) = transpose(conj(a))

weval(x::WElem; kargs...) = WElem(weval(x.val; kargs...))

simplify(x::WElem; kargs...) = weval(WE"Simplify"(x); kargs...)

function n(x::WElem; kargs...)
    x = weval(WE"N"(x); kargs...)
    if x.val isa Real
        x.val
    elseif x.val isa WExpr && x.val.head == W"Complex"
        r, i = x.val.args
        Complex(r, i)
    end
end
n(x::AbstractArray{WElem}; kargs...) = n.(x; kargs...)