# Named parameter containers

A fit with more than a handful of parameters is easier to read — and much
harder to mis-index — when the parameters carry names inside the objective
itself:

```julia
using ComponentArrays, NativeMinuit

start = ComponentArray(mass = 3.872, width = 0.001, coupling = 0.5)

function chi2(p)
    model = breit_wigner.(xdata; m = p.mass, Γ = p.width) .* p.coupling
    sum(((model .- ydata) ./ σ).^2)
end

m = Minuit(chi2, start; error = [1e-3, 1e-4, 0.05])
migrad!(m)
minos!(m)

m.fmin.ext_values.mass      # fitted mass, by name
```

NativeMinuit has no dependency on ComponentArrays and no extension for it.
The support is structural: `x0` is used as the **allocation template** for
every coordinate-shaped workspace, so any `AbstractVector{Float64}` whose
`similar`, `copy`, broadcast and mutation preserve its structure is carried
along for free. `ComponentVector` is the common case; a `LabelledArrays`
`LVector`, or your own array type, works on the same terms.

## What is guaranteed

Your container is used for **every value that crosses back into your code**:

- the objective `f(θ)` and a user-supplied gradient `g(θ)`, on every path —
  MIGRAD, HESSE, MINOS, the contour family, `profile` / `mnprofile`, `scan`,
  SIMPLEX, `extremize` / `profile_band`, the `Optim` bridge, MCMC and the
  Bayesian posterior, and the `threaded_gradient = :auto` safety probe;
- a user prior in `bayesian` / `posterior_sample`;
- the derived-quantity function in `quantiles` and `quantile_band`, and the
  rows from `ens[i]` / iterating a `LikelihoodEnsemble`;
- the reported results: `m.fmin.ext_values`, `m.fmin.ext_errors`,
  `copy(m.values)`, `copy(m.errors)`, the MINOS crossing snapshots
  `m.merrors[name].upper_state` / `.lower_state`, and `extremize`'s endpoint
  vectors `plo` / `phi`.

This holds with bounds, with fixed parameters, and with an analytic gradient.

Your `x0` itself is never mutated — everything is allocated fresh from it.

## Internal coordinates are labelled honestly

MIGRAD does not always work in your coordinates. Two things change the space:

- a **bound** makes the stored number an arcsin- or sqrt-transformed
  coordinate, so it is not the physical parameter value;
- a **fixed** parameter is dropped, so the vector is shorter and slot `i` no
  longer corresponds to external parameter `i`.

In either case the internal workspace is a plain `Vector{Float64}`. It would
be worse than useless to label a transformed coordinate `mass`: reading
`.mass` off it would hand you a number that looks physical and is not.

```julia
mb = Minuit(chi2, start; error = errs, limits = [(3.8, 3.9), nothing, nothing])
migrad!(mb)

mb.fmin.internal.state.parameters.x    # Vector{Float64} — transformed coords
mb.fmin.ext_values.mass                # the physical fitted mass
```

Only when the internal and external spaces coincide — every parameter free
**and** unbounded — does the internal state keep the container, because there
each slot still denotes the parameter its label names.

The rule to remember: **read results from `ext_values` / `m.values`, never
from the internal state.** That was already true before named containers; the
labelling just makes it visible.

## Limits

- `load_ensemble` returns rows as plain vectors. A JSON round-trip cannot
  carry a container, so a reloaded ensemble has no template to rebuild onto.
- `names` still default to `x0`, `x1`, … — the component names are not (yet)
  used as parameter names, so `m.values["x0"]` and `m.merrors["x0"]` keep the
  positional names. Pass `name = ["mass", "width", "coupling"]` if you want
  them to agree.
- Anything you build yourself from `collect(m.values)` is, of course, a plain
  vector. Use `copy(m.values)` when you want the container.
