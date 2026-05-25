# SPDX-License-Identifier: LGPL-2.1-or-later

# ─────────────────────────────────────────────────────────────────────────────
# covariance_squeeze.jl — MnCovarianceSqueeze.
#
# Mirrors reference/Minuit2_cpp/src/MnCovarianceSqueeze.cxx.
#
# Removes a row+column from a symmetric matrix. Three flavors:
#  1. Pure matrix:        H_squeezed = H[setdiff, setdiff]
#  2. MinimumError:       invert → squeeze H → invert back, with fall-back
#                         diagonal if either inversion fails.
#  3. User covariance:    invert covariance → squeeze H → invert back, with
#                         the same fall-back. Used when a parameter is
#                         fixed-then-released or vice versa.
# ─────────────────────────────────────────────────────────────────────────────

"""
    squeeze_symmetric(M::Symmetric, n::Integer) -> Symmetric

Return a new symmetric matrix with row + column `n` removed. Index is
1-based (Julia). Mirrors
`reference/Minuit2_cpp/src/MnCovarianceSqueeze.cxx:89-109` (the
3-argument overload) up to base-index conversion.
"""
function squeeze_symmetric(M::Symmetric{Float64,Matrix{Float64}}, n::Integer)
    p = parent(M)
    nrow = LinearAlgebra.checksquare(p)
    1 <= n <= nrow ||
        throw(ArgumentError("squeeze index $n out of bounds for $(nrow)x$(nrow) matrix"))
    nrow > 1 ||
        throw(ArgumentError("cannot squeeze a 1x1 matrix"))

    out_n = nrow - 1
    out = zeros(Float64, out_n, out_n)
    @inbounds for i in 1:nrow
        i == n && continue
        oi = i < n ? i : i - 1
        for k in i:nrow
            k == n && continue
            ok = k < n ? k : k - 1
            out[oi, ok] = M[i, k]
        end
    end
    return Symmetric(out, :U)
end

"""
    squeeze_error(err::MinimumError, n::Integer; prec=MachinePrecision())
        -> MinimumError

Remove parameter `n` from a `MinimumError`. The error matrix is
inverted to get the Hessian, the Hessian is squeezed, then re-inverted.
On inversion failure returns a diagonal-only `MinimumError` tagged
`MnInvertFailed` (matches C++ behavior at
`reference/Minuit2_cpp/src/MnCovarianceSqueeze.cxx:76-84`).

The returned `MinimumError`'s `dcovar` is preserved from the input
(matches C++ line 86: `MinimumError(squeezed, err.Dcovar())`).
"""
function squeeze_error(err::MinimumError, n::Integer;
                        prec::MachinePrecision = MachinePrecision())
    # Step 1: invert the inverse-Hessian (V) to get the Hessian (H)
    H = copy(parent(err.inv_hessian))
    Hsym = Symmetric(H, :U)
    inv_ok_1 = true
    try
        sym_invert!(Hsym)
    catch
        inv_ok_1 = false
    end

    if !inv_ok_1
        # Diagonal fallback per C++ MnCovarianceSqueeze.cxx:76-84.
        # Mark MnInvertFailed; diagonal = inverted-diagonal of the
        # original error.
        n_in = LinearAlgebra.checksquare(parent(err.inv_hessian))
        diag = zeros(Float64, n_in - 1, n_in - 1)
        oi = 0
        for i in 1:n_in
            i == n && continue
            oi += 1
            diag[oi, oi] = err.inv_hessian[i, i]  # keep original diagonal entry
        end
        return MinimumError(Symmetric(diag, :U), MnInvertFailed)
    end

    # Step 2: squeeze the Hessian
    Hs = squeeze_symmetric(Hsym, n)

    # Step 3: invert squeezed Hessian back to error matrix
    Hs_mat = copy(parent(Hs))
    Hs_sym = Symmetric(Hs_mat, :U)
    inv_ok_2 = true
    try
        sym_invert!(Hs_sym)
    catch
        inv_ok_2 = false
    end

    if !inv_ok_2
        # Inversion of squeezed Hessian failed: diagonal of 1/H_squeezed[i,i]
        n_out = LinearAlgebra.checksquare(parent(Hs))
        diag = zeros(Float64, n_out, n_out)
        for i in 1:n_out
            diag[i, i] = 1.0 / Hs[i, i]
        end
        return MinimumError(Symmetric(diag, :U), MnInvertFailed)
    end

    # Preserve dcovar from the input error (C++ line 86).
    return MinimumError(Hs_sym, err.dcovar)
end
