# SPDX-License-Identifier: LGPL-2.1-or-later

# ─────────────────────────────────────────────────────────────────────────────
# hesse.jl — full numerical Hessian.
#
# Mirrors reference/Minuit2_cpp/src/MnHesse.cxx:93-330 (the
# operator()(MnFcn, MinimumState, MnUserTransformation, maxcalls)
# overload — the "real" MnHesse).
#
# Algorithm:
#   1. Evaluate FCN at the current point (amin).
#   2. Diagonal pass: for each parameter i, multi-cycle central-difference
#      refinement to determine g2[i] = ∂²f/∂x_i². Each cycle:
#        - Find a step d such that sag = (f(x+d) + f(x-d) - 2·f) ≠ 0
#          (multiplier loop up to 5× growth, bounded at 0.5 if param
#          has limits).
#        - Update d using `d = sqrt(2·aimsag / |g2|)`; clamp to
#          [dmin, 10·dlast] / [0.1·dlast, 0.5 if has-limits].
#        - Break if d-step or g2 has converged below the strategy
#          tolerances.
#      vhmat[i, i] = g2[i].
#   3. (Strategy > 0): refine gradient via HessianGradientCalculator.
#      Phase 1 first cut SKIPS this and uses the gradient as-is. The
#      refined-gradient path lands in a Phase 1 follow-up.
#   4. Off-diagonal pass: for each pair (i, j) with i < j, compute
#      `(f(x + d_i + d_j) + f(x) - f(x + d_i) - f(x + d_j)) / (d_i d_j)`.
#      Uses cached single-direction values `yy[i]` from the diagonal
#      pass.
#   5. MnPosDef enforcement.
#   6. Sym invert. If fails → MnInvertFailed diagonal matrix.
#   7. New EDM via the standard estimator.
#   8. Return new MinimumState with the updated MinimumError.
#
# This is the standalone MnHesse. Calling it from inside MIGRAD when
# Strategy ≥ 1 + Dcovar > 0.05 (the inner-HESSE refinement) is the
# Phase 1 integration step — see migrad.jl follow-up.
# ─────────────────────────────────────────────────────────────────────────────

"""
    hesse(cf, state, strategy=Strategy(1); prec=MachinePrecision(), maxcalls=0)
        -> MinimumState

Compute the full numerical Hessian at `state.parameters.x` and return
a new `MinimumState` with the refined error matrix, recomputed EDM,
and updated FCN call count.

Mirrors C++ `MnHesse::operator()(MnFcn, MinimumState,
MnUserTransformation, maxcalls)` —
`reference/Minuit2_cpp/src/MnHesse.cxx:93-330`.

# Arguments

- `cf::CostFunction` — the user FCN (operates on the parameter
  vector that `state.parameters.x` reports — Phase 1 first cut: no
  bounds, so internal == external; bounded HESSE is a follow-up).
- `state::MinimumState` — current state. The gradient field provides
  initial step sizes (`gst[i] = state.gradient.gstep[i]`) and the
  algorithm refines `g2[i]`.
- `strategy::Strategy` — controls `hessian_ncycles` (cycles per
  parameter), `hessian_step_tolerance`, `hessian_g2_tolerance`,
  and (Strategy ≥ 1) gradient refinement.
- `prec::MachinePrecision` — floor for step sizes and pos-def gate.
- `maxcalls::Integer` — FCN call cap; `0` means use the default
  `200 + 100n + 5n²` per `MnApplication.cxx:43`.

# Return statuses

The returned `MinimumState`'s `error.status` may be:
- `MnHesseValid` — full success.
- `MnMadePosDef` — pos-def perturbation was applied.
- `MnHesseFailed` — sag stayed zero or maxcalls hit; matrix is
  diagonal `1/g2[i]` (or 1 where g2 is too small).
- `MnInvertFailed` — inversion failed; matrix is the same diagonal.
"""
function hesse(
    cf::CostFunction,
    state::MinimumState,
    strategy::Strategy = Strategy(1);
    prec::MachinePrecision = MachinePrecision(),
    maxcalls::Integer = 0,
)
    n = length(state)
    if maxcalls == 0
        maxcalls = 200 + 100 * n + 5 * n * n
    end

    x = copy(state.parameters.x)
    amin = cf(x)
    aimsag = sqrt(prec.eps2) * (abs(amin) + cf.up)

    # Scratch — independent vectors so we don't mutate state.gradient
    g2 = copy(state.gradient.g2)
    gst = copy(state.gradient.gstep)
    grd = copy(state.gradient.grad)
    dirin = copy(gst)
    yy = zeros(Float64, n)

    vhmat = zeros(Float64, n, n)

    # No `has_limits` info per parameter in Phase 1 first-cut (bounds
    # integration is the follow-up). Treat all params as unbounded for
    # the d-clamp branches; this matches the migrad.jl Phase 0 caller.
    has_limits = false

    # ── Diagonal pass ─────────────────────────────────────────────
    for i in 1:n
        xtf = x[i]
        dmin = 8.0 * prec.eps2 * (abs(xtf) + prec.eps2)
        d = abs(gst[i])
        d < dmin && (d = dmin)

        fs1 = 0.0
        fs2 = 0.0
        sag = 0.0
        converged = false

        for icyc in 1:strategy.hessian_ncycles
            # Multiplier loop — grow d until sag ≠ 0
            sag = 0.0
            mlp_failed = false
            for multpy in 1:5
                x[i] = xtf + d
                fs1 = cf(x)
                x[i] = xtf - d
                fs2 = cf(x)
                x[i] = xtf
                sag = 0.5 * (fs1 + fs2 - 2.0 * amin)
                sag != 0 && break
                if has_limits
                    d > 0.5 && (mlp_failed = true; break)
                    d *= 10
                    d > 0.5 && (d = 0.51)
                else
                    d *= 10.0
                end
                multpy == 5 && (mlp_failed = true)
            end

            if sag == 0 || mlp_failed
                # Sag stayed zero → return failure with diagonal matrix
                return _hesse_diagonal_failure(state, g2, prec, ncalls(cf), MnHesseFailed)
            end

            g2bfor = g2[i]
            g2[i] = 2.0 * sag / (d * d)
            grd[i] = (fs1 - fs2) / (2.0 * d)
            gst[i] = d
            dirin[i] = d
            yy[i] = fs1
            dlast = d

            d = sqrt(2.0 * aimsag / abs(g2[i]))
            has_limits && (d = min(0.5, d))
            d < dmin && (d = dmin)

            # Convergence checks
            if abs((d - dlast) / d) < strategy.hessian_step_tolerance
                converged = true
                break
            end
            if g2[i] != 0 && abs((g2[i] - g2bfor) / g2[i]) < strategy.hessian_g2_tolerance
                converged = true
                break
            end
            d = min(d, 10.0 * dlast)
            d = max(d, 0.1 * dlast)
        end

        vhmat[i, i] = g2[i]

        if ncalls(cf) > maxcalls
            return _hesse_diagonal_failure(state, g2, prec, ncalls(cf), MnHesseFailed)
        end
    end

    # ── Strategy > 0: gradient refinement (Phase 1 follow-up) ──
    # The C++ HessianGradientCalculator refines `grd` and `gst` here when
    # strategy.level > 0. Phase 1 first cut keeps the current grad; the
    # full refinement lands when hessian_gradient.jl is ported.

    # ── Off-diagonal pass ─────────────────────────────────────────
    # All pairs (i, j) with i < j.
    if n > 1
        for i in 1:n
            x[i] += dirin[i]
            for j in (i + 1):n
                x[j] += dirin[j]
                fs1 = cf(x)
                vhmat[i, j] = (fs1 + amin - yy[i] - yy[j]) / (dirin[i] * dirin[j])
                x[j] -= dirin[j]
            end
            x[i] -= dirin[i]
        end
    end

    # ── Pos-def enforcement on the H matrix ───────────────────────
    # vhmat is the Hessian (second derivatives); MnPosDef ensures
    # positive-definiteness on the to-be-inverted matrix.
    err_tmp = make_posdef(MinimumError(Symmetric(vhmat, :U), 1.0), prec)
    vhmat_pd = copy(parent(err_tmp.inv_hessian))

    # ── Symmetric invert ──────────────────────────────────────────
    inv_ok = true
    try
        sym_invert!(Symmetric(vhmat_pd, :U))
    catch _
        inv_ok = false
    end

    if !inv_ok
        return _hesse_diagonal_failure(state, g2, prec, ncalls(cf), MnInvertFailed)
    end

    # ── New EDM with the refined error ────────────────────────────
    refined_grad = FunctionGradient(grd, g2, gst)
    new_err_status = is_made_pos_def(err_tmp) ? MnMadePosDef : MnHesseValid
    new_dcov = is_made_pos_def(err_tmp) ? 1.0 : 0.0
    new_err = MinimumError(Symmetric(vhmat_pd, :U), new_dcov, new_err_status, true)
    new_edm = estimate_edm(refined_grad, new_err)

    return MinimumState(state.parameters, new_err, refined_grad,
                        new_edm, ncalls(cf))
end

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    _hesse_diagonal_failure(state, g2, prec, nfcn, status) -> MinimumState

Build the failure-mode `MinimumState` with a diagonal inverse-Hessian
of `1/g2[i]` (clamped to 1 when g2 is too small). Mirrors C++
`MnHesse.cxx:177-184` (and the analogous block at lines 216-223 for
maxcalls overrun).
"""
function _hesse_diagonal_failure(state::MinimumState, g2::Vector{Float64},
                                  prec::MachinePrecision, nfcn::Integer,
                                  status::CovStatus)
    n = length(g2)
    M = zeros(Float64, n, n)
    @inbounds for j in 1:n
        tmp = g2[j] < prec.eps2 ? 1.0 : 1.0 / g2[j]
        M[j, j] = tmp < prec.eps2 ? 1.0 : tmp
    end
    err = MinimumError(Symmetric(M, :U), status)
    return MinimumState(state.parameters, err, state.gradient,
                        state.edm, nfcn)
end
