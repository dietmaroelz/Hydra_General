"""
    Input a matrix of nodes and convert to matrices x, y, z in the format needed for GLMakie.surface
"""
function Nodes2xyz(Nodes::Matrix; wrap::Bool=true)
    x, y, z = zeros(M, N), zeros(M, N), zeros(M, N)
    #m=1
    x[1,:] = fill(Nodes[1,1], N)
    y[1,:] = fill(Nodes[1,2], N)
    z[1,:] = fill(Nodes[1,3], N)
    for m in 2:M
        for n in 1:N
            node = Nodes[mn2i(m,n),:]
            x[m,n] = node[1]
            y[m,n] = node[2]
            z[m,n] = node[3]
        end
    end

    if wrap
        x = hcat(x, x[:,1])
        y = hcat(y, y[:,1])
        z = hcat(z, z[:,1])
    end
    return x, y, z
end

"""
    Convert a matrix of nodes to a vector of point vectors, with length noNodes
"""
function Nodes2Points(Nodes::Matrix)
    return [Nodes[i,:] for i in 1:size(Nodes)[1]]
end


"""
    Calculate the distance between two points, provide integers to find the distance between nodes corresponding to those indices
"""
Distance(n1::Vector, n2::Vector) = norm(n2 - n1)
Distance(i1::Int, i2::Int, Nodes::Matrix) = Distance(Nodes[i1,:], Nodes[i2,:])


struct Spheroid
    Nodes::Matrix
    Triangles::Vector{Vector{Int}}
    Psivals::Vector{Float64}
    Spheroid(Nodes::Matrix, Triangles::Vector{Vector{Int}}, Psivals::Vector{Float64}) = new(Nodes, Triangles, Psivals)
end
Spheroid(Spher::Spheroid) = Spheroid(Spher.Nodes, Spher.Triangles, Spher.Psivals)
Distance(i1::Int, i2::Int, spher::Spheroid) = Distance(i1, i2, spher.Nodes)
length(Spher::Spheroid) = size(Spher.Nodes)[1]
NoTriangles(Spher::Spheroid) = length(Spher.Triangles)


"""
    Find the Volume of sphere by integrating (x,y,z)ᵀ ⋅ n̂ across the discretised surface. At the moment this integration is only performed at a single point on each triangle, but it could be improved (made precise?) using a 7-point integration. 
"""
function SphereVolume(Nodes::Matrix, Triangles::Vector)
    Vtrack = 0.0
    for T in Triangles
        centre = sum(Nodes[t,:] for t in T)/3
        v1, v2 = Nodes[T[2],:] - Nodes[T[1],:], Nodes[T[3],:] - Nodes[T[1],:]
        outernormal = cross(v1, v2)

        outernormal = (dot(centre, outernormal) < 0) ? (outernormal * -1/norm(outernormal)) : outernormal / norm(outernormal)
        Vtrack += (TriArea([Nodes[t,:] for t in T]...) * dot(centre, outernormal))/3
    end
    return Vtrack
end

SphereVolume(Spher::Spheroid) = SphereVolume(Spher.Nodes, Spher.Triangles)

""" Unit normal vector associated with a specified triangle """
function NormalToTri(p1, p2, p3)
    v1, v2 = p2-p1, p3-p1
    normal = cross(v1, v2)
    return normal/norm(normal)
end
NormalToTri(Nodes::Matrix, T::Vector{Int}) = NormalToTri([Nodes[t,:] for t in T]...)

"""
Area of a triangle
"""
function TriArea(p1::Vector, p2::Vector, p3::Vector)
    @assert length(p1) == length(p2) == length(p3) == 3 "All points should be 3D"
    v1 = p2 .- p1
    v2 = p3 .- p1
    return norm(cross(v1, v2))/2
end

function TriArea(i1::Int, i2::Int, i3::Int, Nodes::Matrix)
    p1, p2, p3 = Nodes[i1,:], Nodes[i2,:], Nodes[i3,:]
    return TriArea(p1, p2, p3)    
end

"""
    Find the areas of all faces specified by Triangles
"""
function FaceAreas(Nodes, Triangles)
    return [TriArea(T..., Nodes) for T in Triangles]
end
FaceAreas(Spher::Spheroid) = FaceAreas(Spher.Nodes, Spher.Triangles)