module WolframLinearAlgebra

import Base: OneTo, HasShape, promote_rule, convert, zero, one, size, axes, ndims, length
import Base: firstindex, lastindex, IteratorSize, keys, first, last, copy, promote_op
import Base: +, -, *, /, ^, sqrt, conj, adjoint, transpose
import Base: sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
import Base: log, log2, log10, exp

import LinearAlgebra: Eigen, eigen, SVD, svd, Schur, schur
import LinearAlgebra: qr, lq
import LinearAlgebra: pinv, nullspace, diag

import MathLink: WSymbol, WExpr, WInteger, WReal, weval, parseexpr, @W_str, @W_cmd

export WElem, WSolver, @WE_str, @WE_cmd, @WS_str, @W_str, @W_cmd
export weval, simplify, n
export eigen, svd, qr, lq, schur, pinv, nullspace

iswlist(x) = false
iswlist(x::WExpr) = x.head == W"List"

include("welem.jl")
include("wsolver.jl")
include("linargfunctions.jl")

end # module
