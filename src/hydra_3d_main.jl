using GLMakie
using LinearAlgebra
using Statistics
using PolygonInbounds
using LineSearches
using Optim
using ForwardDiff
using BenchmarkTools
using DifferentialEquations
using GeometryBasics 
using FileIO
using Loess


"""
    Loads a data file created by SaveRun and appends it to the global vectors in the current workspace. 
"""
function LoadData(Name::String)
    return load("data\\"*Name*".jld2")["s"]
end

import Main.length
import Main.show

const R = 1 #sphere initial radius
const η = 10.0 #Volume preservation
const κ0 = 50.0 #modulate contractility
const σ0 = 1.0 #chemical secretion rate
const ξ = 1.0 #Nodewise friction parameter
const α = 2.0 #chemical decay rate
const β = 5.0 #Poisson ratio of medium

const V0sphere = 4π*R^3/3
const A0sphere = 4π*R^2

function kappa(Ψ)
    return κ0*exp(-β*Ψ)
end

include("hydra_3d_auxillary.jl") #contains various helper functions, as well as spheroid definition

"""
    struct to store a simulation, Spherevec is a vector of Spheroidal configurations (user-defiend object), tvec is the time corresponding to each configuration and RunTime is the computation time needed to run the simulation 
"""
struct ODESimulation
    Spherevec::Vector{Spheroid}
    tvec::Vector{Float64}
    RunTime::Float64
end

"""
    Load in one of the pre-assembled meshes of a sphere and convert it into a Spheroid compatible with this code
"""
function MeshToSpheroid(fn::String)
    mesh = load("meshes//"*fn)
    vertexcoordinates = coordinates(mesh)
    f = faces(mesh)
    alltriangles=[GeometryBasics.value.(f[i])[1:end] for i=1:length(f)]

    noNodes = length(vertexcoordinates)
    Nodes = zeros(noNodes, 3)
    for (i,v) in enumerate(vertexcoordinates)
        Nodes[i,:] = v
    end
    return Spheroid(Nodes, alltriangles, zeros(length(alltriangles)))
end

## Load in initial condition and define some convenience constants from it
# Spher0 = MeshToSpheroid("SphereMesh_Even81.off")
Sphere0 = LoadData("IC_Even81")
Nodes = Sphere0.Nodes
const noNodes = size(Sphere0.Nodes)[1]

include("hydra_3d_plotting.jl") #various plotting and visualisation helper functions
include("hydra_3d_stretch_energy.jl") #evaluating the Green's strain energy density for 2D objects
include("hydra_3d_energy.jl") #energy functional to implement EOM's
include("hydra_3d_gradientdescent_decoupled.jl") #Includes gradient-descent functionality

const V0numeric = SphereVolume(Sphere0) #numerical volume of discretised sphere, for coarse discretisations this won't well approximate 4πR^3/3, so the numerical version is used as the basis
const V0target = V0numeric*1.2 #target volume
const Areas0numeric = FaceAreas(Sphere0) #original area of each triangle, used to define strain energy
const S0numeric = sum(Areas0numeric) #total unperturbed surface area (discretised)

psivals0 = zeros(noNodes)

"""
    Return the θ, ϕ coordinate of a position in R³
"""
function Node2angle(x, y, z)
    iszero(y) && return z>0 ? (0., 0.) : (π, 0.)
    return acos(z/sqrt(x^2+y^2+z^2)), sign(y)*acos(x/sqrt(x^2+y^2))
end
Node2angle(Nodes, i) = Node2angle(Nodes[i,:]...)


"""
    Chemistry as a linear 'hat' function centred at (θ, ϕ) = (π/2, π)
"""
function psi1(theta, phi)
    # return max(1.5-abs(theta-π/2),0) * max(1.5-abs(phi-π),0)
    return max(1.5-theta, 0.0)
end

psivals1 = [mean(psi1(Node2angle(Nodes, n)...) for n in Tri) for Tri in Sphere0.Triangles];
Spher1 = Spheroid(Sphere0.Nodes, Sphere0.Triangles, psivals1)


""" Set concentration uniformly Ψ=1 everywhere"""
psivals2 = ones(length(Sphere0.Triangles));
Spher2 = Spheroid(Sphere0.Nodes, Sphere0.Triangles, psivals2)

""" A more extreme version of the hat function psi1"""
function psi3(theta, phi)
    return max(2.5-abs(theta-π/2),0) * max(2.5-abs(phi-π),0)
end
psivals3 = [mean(psi3(Node2angle(Nodes, n)...) for n in Tri) for Tri in Sphere0.Triangles];
Sphere3 = Spheroid(Sphere0.Nodes, Sphere0.Triangles, psivals3)

""" A psi function with smooth variance along θ, used to create a dumbbell """
function psi4(theta, phi)
    return min(abs(π/2-theta), 1.0)
end
psivals4 = [mean(psi4(Node2angle(Nodes, n)...) for n in Tri) for Tri in Sphere0.Triangles];
Sphere4 = Spheroid(Sphere0.Nodes, Sphere0.Triangles, psivals4)

""" A psi function with Ψ as a step function along θ, best for creating a dumbbell"""
function psi5(theta, phi)
    return abs((theta-π/2)) < π/8 ? 0.0 : 1.0
end
psivals5 = [mean(psi5(Node2angle(Nodes, n)...) for n in Tri) for Tri in Sphere0.Triangles];
Sphere5 = Spheroid(Sphere0.Nodes, Sphere0.Triangles, psivals5)

""" Initially set Ψ as a small random number to test for instantaneous pattern formation """
function psi6(theta, phi)
    return rand()/3
end
psivals6 = [mean(psi6(Node2angle(Nodes, n)...) for n in Tri) for Tri in Sphere0.Triangles];
Spher6 = Spheroid(Sphere0.Nodes, Sphere0.Triangles, psivals6)


"""
    Save a run into our data file
"""
function SaveData(data::Union{Spheroid, ODESimulation}, filename::String)
    @assert !isfile("data\\$filename.jld2") "A file already exists with that name"
    save("data\\$filename.jld2", "s", data)
end
