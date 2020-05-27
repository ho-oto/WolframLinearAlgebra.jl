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
