""" Area of a triangle provided in matrix form, assume 2D coordinates for the triangle"""
function Area(Tri::Matrix{Float64})
    @assert size(Tri) == (3,2) "wrong dimensions"
    v1 = Tri[2,:] - Tri[1,:]
    v2 = Tri[3,:] - Tri[1,:]
    return norm(cross([v1; 0.], [v2; 0.]))/2
end

""" Make the D matrix for given triangle, contains the basis vectors of triangular coordinates """
function MakeD(Tri::Matrix{T}) where T <: Real
    @assert size(Tri) == (3,2) "wrong dimensions"
    return [
        Tri[2,1]-Tri[1,1] Tri[3,1]-Tri[1,1]; 
        Tri[2,2]-Tri[1,2] Tri[3,2]-Tri[1,2]
    ]
end

""" Make the Greens strain tensor for original triangle T0 and transformed triangle T1 """
function MakeG(T0::Matrix{Float64}, T1::Matrix{T}) where T <: Real
    @assert size(T0) == size(T1) == (3,2) "wrong dimensions"
    Fmat = MakeD(T1) * inv(MakeD(T0))
    return (transpose(Fmat)*Fmat - I)/2
end

""" Find the stretching energy associated with the specified transformation, FIX THE PARAMETERS HERE """
function StretchEnergy(T0::Matrix{Float64}, T1::Matrix{T})::T where {T <: Real}
    @assert size(T0) == size(T1) == (3,3) "wrong dimensions"
    T02D, T12D = TriTo2D(T0), TriTo2D(T1)
    G = MakeG(T02D, T12D)
    return Area(T02D)*(tr(G^2) + tr(G)^2)
end

""" transform a triangle in R³ to one in R² with a line along the x axis, preserve internal angles and edge lengths"""
function TriTo2D(Tri::Matrix{T})::Matrix{T} where T <: Real
    @assert size(Tri) == (3,3)
    AB, AC = Tri[2,:] - Tri[1,:], Tri[3,:] - Tri[1,:]
    temp = dot(AB, AC)/norm(AB)/norm(AC)
    θA = temp >= 1 ? acos(1.0) : acos(temp)
    return [
        0. 0.;
        norm(AB) 0.;
        norm(AC)*cos(θA) norm(AC)*sin(θA)
    ]
end

""" Find the individual stretching energies for all triangles within pre-transformation Spher0, and post-transformation Spher1 """
function StretchingEnergies(Spher0::Spheroid, Spher1::Spheroid)
    @assert Spher0.Triangles == Spher1.Triangles "I know this is annoying but it's what we're doing for now, besides it's a useful check for compatability"
    return StretchingEnergies(Spher0, Spher1.Nodes)
end

function StretchingEnergies(Spher0::Spheroid, Nodes::Matrix{T}) where T <: Real
    Evals = zeros(T, length(Spher0.Triangles))
    for (i,Tri) in enumerate(Spher0.Triangles)
        Evals[i] = StretchEnergy(Spher0.Nodes[Tri,:], Nodes[Tri,:])
    end
    return Evals
end