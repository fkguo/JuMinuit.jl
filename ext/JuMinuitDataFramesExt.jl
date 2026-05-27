# SPDX-License-Identifier: LGPL-2.1-or-later
#
# ─────────────────────────────────────────────────────────────────────────────
# JuMinuitDataFramesExt — `Data(::DataFrame)` IMinuit.jl drop-in compat.
#
# IMinuit.jl exposes `Data(df::DataFrame) = Data(df[:,1], df[:,2], df[:,3])`
# (assumes 3 columns interpreted as x, y, err). This was missing from the
# JuMinuit core API, breaking notebooks that pass DataFrames directly
# (e.g., BenchmarkExamples/IAM_2Pformfactor/iamfit.ipynb).
#
# Added as a Package Extension rather than a hard dependency: DataFrames
# pulls in PrettyTables + Compat + a stack of CSV-ecosystem deps, so
# making it mandatory would inflate JuMinuit's install footprint. Users
# who already have DataFrames loaded get the convenience method for free.
# ─────────────────────────────────────────────────────────────────────────────

module JuMinuitDataFramesExt

using JuMinuit
using DataFrames

"""
    JuMinuit.Data(df::DataFrame)

Construct a `Data` from a 3-column DataFrame; columns interpreted as
`(x, y, err)` by position (not by name). Mirrors IMinuit.jl's
`Data(::DataFrame)` for drop-in notebook compatibility.

# Examples

```julia
using DataFrames, CSV, JuMinuit
df = DataFrame(CSV.File("data.csv"; header=[:w, :y, :err]))
d = Data(df)        # → Data(df[:,1], df[:,2], df[:,3])
```

If you need explicit column selection by name, use the 3-arg form:
`Data(df.w, df.y, df.err)`.
"""
JuMinuit.Data(df::DataFrame) = JuMinuit.Data(df[:, 1], df[:, 2], df[:, 3])

end # module JuMinuitDataFramesExt
