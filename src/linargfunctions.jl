function eigen(x::Matrix{WElem})
    @info "use wolframscript"
    u, v = WS"Eigensystem"(x)
    Eigen(u, transpose(v))
end

function svd(x::Matrix{WElem})
    @info "use wolframscript"
    u, v, w = WS"SingularValueDecomposition"(x)
    SVD(u, diag(v), w')
end

function qr(x::Matrix{WElem})
    @info "use wolframscript"
    q, r = WS"QRDecomposition"(x)
    (Q = q', R = r)
end

function lq(x::Matrix{WElem})
    @info "use wolframscript"
    q, r = WS"QRDecomposition"(permutedims(x, (2,1)))
    (L = transpose(r), Q = conj(q))
end

function schur(x::Matrix{WElem})
    @info "use wolframscript"
    z, t = WS"SchurDecomposition"(x)
    Schur(t, z, diag(t))
end

inv(x::Matrix{WElem}) = WS"Inverse"(x, true)
pinv(x::Matrix{WElem}) = WS"PseudoInverse"(x, true)
nullspace(x::Matrix{WElem}) = permutedims(WS"NullSpace"(x, true), (2,1))
