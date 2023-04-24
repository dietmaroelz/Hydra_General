"""
    Conservative components of the energy functional
"""
function EnergyFunctional(Nodes::Matrix{T}, Triangles::Vector{Vector{Int}}, Psivals::Vector{S}; Loud::Bool=false)::Real where {T <: Real, S <: Real}
    global LastNodes = Nodes
    Vnodes = SphereVolume(Nodes, Triangles)
    EVolume = η/2*(V0target - Vnodes)^2

    facekappavals = kappa.(Psivals)
    Estretch = sum(facekappavals .* StretchingEnergies(Sphere0, Nodes))


    if Loud
        println("Volume energy: $EVolume")
        println("Stretching energy: $Estretch")
        println()
        println("Total energy: $(EVolume + Estretch)")
    end
    return EVolume + Estretch
end

function EnergyFunctional(xvec::Vector, Triangles::Vector{Vector{Int}}, Psivals::Vector; kwargs...)::Real
    Nodes= xvec2Nodes(xvec)
    return EnergyFunctional(Nodes, Triangles, Psivals; kwargs...)
end
EnergyFunctional(Spher::Spheroid; args...)::Real = EnergyFunctional(Nodes2xvec(Spher.Nodes), Spher.Triangles, Spher.Psivals; args...)
EnergyFunctional(Sim::ODESimulation, t::Int; args...) = EnergyFunctional(Sim.Sphervec[t]; args...)

"""
    Convert Nodes into a single vector for use inside an energy functional, ignoring entries that aren't degrees of freedom - the first X, Y values, the final N-ring of Z values
"""
function Nodes2xvec(Nodes)
    xvec = Nodes[:, 1] #X values
    append!(xvec, Nodes[:,2]) #Y values
    append!(xvec, Nodes[:,3]) #Z values
    return xvec
end

"""
    Inputs xvec in the form used by the energy functional, converts to Nodes and HalfNodes by reendowing points that aren't degrees of freedom 
"""
function xvec2Nodes(xvec)
    X = xvec[1:noNodes]
    Y = xvec[noNodes+1:2*noNodes]
    Z = xvec[2*noNodes+1:3*noNodes]
    return hcat(X,  Y,  Z)
end

