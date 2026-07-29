# Derived-target dimension verification

- Run ID: `20260729T054838Z-derived-dimension-docs-r1`
- Baseline: `aa8d4a7`
- Scope: the dimension entering the profile-likelihood $\Delta\chi^2$
  threshold for scalar and vector-valued derived targets

## Independent route 1: profile likelihood and Wilks

For a target $\psi=g(\theta)$, define the profile-likelihood-ratio statistic

$$
q(\psi_0)=-2\log\frac{\sup_{\theta:g(\theta)=\psi_0}L(\theta)}
{\sup_\theta L(\theta)}.
$$

Under local asymptotic normality, an interior identifiable truth, nonsingular
information in the identifiable directions, and a differentiable target map
of constant local rank, the null constraint has tangent space
$Dg(\theta_0)h=0$. The loss of optimized log likelihood is the squared Fisher
distance from the unconstrained estimator to that tangent space. Therefore

$$
q\left(g(\theta_0)\right)\xrightarrow{d}\chi_r^2,
\qquad
r=\operatorname{rank}Dg(\theta_0)
$$

where the rank is taken on the identifiable parameter directions.

## Independent route 2: Gaussian geometry

In the local linear-Gaussian limit, let
$\widehat\theta-\theta_0\sim\mathcal N(0,C)$ and
$J=Dg(\theta_0)$. Then the derived displacement has covariance
$\Sigma=JCJ^{\mathsf T}$ and

$$
\left[g(\widehat\theta)-g(\theta_0)\right]^{\mathsf T}
\Sigma^+
\left[g(\widehat\theta)-g(\theta_0)\right]
\sim\chi_r^2,
\qquad
r=\operatorname{rank}\!\left(JC^{1/2}\right).
$$

When $C$ is nonsingular on the identifiable space, this is the same rank as
in the Wilks derivation. The two derivations are therefore equivalent in
their common local-quadratic regime.

## Agreed consequences

| Reported target | Effective dimension | 68.27% threshold |
| --- | ---: | ---: |
| One non-degenerate scalar, including one fixed-energy spectrum value, a real scattering length, or one pole component | 1 | $\Delta\chi^2=1.00$ |
| Joint $(\operatorname{Re}E_{\mathrm{pole}},\operatorname{Im}E_{\mathrm{pole}})$ with two locally independent directions | 2 | $\Delta\chi^2=2.2957489289\simeq2.30$ |
| Real and imaginary pole parts quoted separately | 1 for each interval | $\Delta\chi^2=1.00$ for each; not a 68.27% joint region |

At each fixed abscissa, a curve value is a scalar target. A collection of
such intervals is a pointwise band. Simultaneous coverage of the entire curve
requires separate calibration or a joint region for the identifiable curve
directions; the plotting-grid length is not the number of degrees of freedom.

## Limits of the result

The standard $\chi^2_r$ threshold is asymptotic outside the exactly
linear-Gaussian case. Boundaries, non-identifiability, a zero or changing
rank, non-smooth branch changes, or small samples can invalidate its nominal
coverage. Profiling respects nuisance parameters and parameter limits but
does not, by itself, calibrate a non-regular likelihood-ratio statistic.

The two blind derivations were method-diverse but from the same model family.
An impartial comparison found all three claims mathematically equivalent with
no outliers. This verifies the abstract statistical statements above; it does
not verify a particular fitting implementation or finite-sample coverage.
