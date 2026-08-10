# SPDX-License-Identifier: LGPL-2.1-or-later

module NativeMinuitComponentArraysExt

using ComponentArrays
using NativeMinuit

struct ComponentArrayAdapter{A} <: NativeMinuit.AbstractParameterAdapter
    axes::A
end

NativeMinuit._parameter_adapter(x0::ComponentArrays.ComponentVector) =
    ComponentArrayAdapter(ComponentArrays.getaxes(x0))

NativeMinuit._restore_parameters(adapter::ComponentArrayAdapter,
                                 x::AbstractVector) =
    ComponentArrays.ComponentVector(x, adapter.axes)

NativeMinuit._adapt_parameter_view(adapter::ComponentArrayAdapter,
                                   view::AbstractVector) =
    NativeMinuit._restore_parameters(adapter, view)

end # module NativeMinuitComponentArraysExt
