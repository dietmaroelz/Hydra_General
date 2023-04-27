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
const η = 100.0 #Volume preservation
const κ0 = 5.0 #modulate contractility
const κ1 = 50.0
const σ0 = 20.0 #chemical secretion rate
const ξ = 1.0 #Nodewise friction parameter
const α = 1.0 #chemical decay rate
const β = 5.0 #Poisson ratio of medium

const V0sphere = 4π*R^3/3
const A0sphere = 4π*R^2

kappa(Ψ) = κ1*exp(-β*Ψ) + κ0


include("hydra_3d_auxillary.jl") #contains various helper functions, as well as spheroid definition

"""
    struct to store a simulation, Spherevec is a vector of Spheroidal configurations (user-defined object), tvec is the time corresponding to each configuration and RunTime is the computation time needed to run the simulation 
"""
struct Simulation
    Spherevec::Vector{Spheroid}
    tvec::Vector{Float64}
    V0vec::Vector{Float64}
    RunTime::Float64
end
length(Sim::Simulation) = length(Sim.Spherevec)

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
include("hydra_3d_updatescheme.jl") #Implement the energy minisation into simulation structs
# include("hydra_3d_gradientdescent.jl") #Includes gradient-descent functionality, currently I'm working with the energy minisation scheme for simplicity

const V0numeric = SphereVolume(Sphere0) #numerical volume of discretised sphere, for coarse discretisations this won't well approximate 4πR^3/3, so the numerical version is used as the basis
const V0target = V0numeric*1.5 #target volume
const Areas0numeric = FaceAreas(Sphere0) #original area of each triangle, used to define strain energy
const S0numeric = sum(Areas0numeric) #total unperturbed surface area (discretised)

psivals0 = zeros(noNodes)

""" Find the gradient corresponding to the Phi update using the deformation calculated by comparing Sphere to Sphere0. An alternate method specifies the current configuration by a Nodes matrix rather than Sphere """
function GradientPsi(psivals::Vector{Float64}, Sphere::Spheroid, Spher0::Spheroid)
    @assert Sphere.Triangles == Spher0.Triangles
    return GradientPsi(psivals, Sphere.Nodes, Spher0)
end
GradientPsi(psivals::Vector{Float64}, xvec::Vector{Float64}, Spher0::Spheroid) = GradientPsi(psivals, xvec2Nodes(xvec), Spher0)

function GradientPsi(psivals::Vector{Float64}, Nodes::Matrix{Float64}, Spher0::Spheroid; n=10, K=2.0)
    FaceStretching = StretchingEnergies(Spher0, Nodes)
    # return zeros(length(psivals))
    # return σ0*FaceStretching.^n./(K^n .+ FaceStretching.^n) - α*psivals
    return σ0*FaceStretching - α*psivals
end

"""
    Return the θ, ϕ coordinate of a position in R³
"""
function Node2angle(x, y, z)
    iszero(y) && return z>0 ? (0., 0.) : (π, 0.)
    return acos(z/sqrt(x^2+y^2+z^2)), sign(y)*acos(x/sqrt(x^2+y^2))
end
Node2angle(Nodes, i) = Node2angle(Nodes[i,:]...)

""" Set concentration uniformly Ψ=1 everywhere"""
psiUniform = ones(length(Sphere0.Triangles));
SphereUniform = Spheroid(Sphere0.Nodes, Sphere0.Triangles, psiUniform)

""" A psi function with smooth variance along θ, used to create a dumbbell """
psiSmooth(theta, phi) = min(abs(π/2-theta), 1.0)
psivalsSmooth = [mean(psiSmooth(Node2angle(Nodes, n)...) for n in Tri) for Tri in Sphere0.Triangles];
SphereSmooth = Spheroid(Sphere0.Nodes, Sphere0.Triangles, psivalsSmooth)

""" A psi function with Ψ as a step function along θ, best for creating a dumbbell"""
psiStep(theta, phi) = abs((theta-π/2)) < π/8 ? 0.0 : 1.0
psivalsStep = [mean(psiStep(Node2angle(Nodes, n)...) for n in Tri) for Tri in Sphere0.Triangles];
SphereStep = Spheroid(Sphere0.Nodes, Sphere0.Triangles, psivalsStep)

""" Initially set Ψ as a small random number to test for instantaneous pattern formation """
psiRandom(theta, phi) = rand()
psivalsRandom = [mean(psiRandom(Node2angle(Nodes, n)...) for n in Tri) for Tri in Sphere0.Triangles];
SphereRandom = Spheroid(Sphere0.Nodes, Sphere0.Triangles, psivalsRandom)

psiMonopole(theta, phi) = theta<π/2 ? 1.0 : 0.0
psivalsMonopole = [mean(psiMonopole(Node2angle(Nodes, n)...) for n in Tri) for Tri in Sphere0.Triangles];
SphereMonopole = Spheroid(Sphere0.Nodes, Sphere0.Triangles, psivalsMonopole)

psiDimple(theta, phi) = theta<π/4 ? 1.0 : 0.0
psivalsDimple = [mean(psiDimple(Node2angle(Nodes, n)...) for n in Tri) for Tri in Sphere0.Triangles];
SphereDimple = Spheroid(Sphere0.Nodes, Sphere0.Triangles, psivalsDimple)

psiCosine(theta, phi) = 0.5 + cos(3theta)/2
psivalsCosine = [mean(psiCosine(Node2angle(Nodes, n)...) for n in Tri) for Tri in Sphere0.Triangles];
SphereCosine = Spheroid(Sphere0.Nodes, Sphere0.Triangles, psivalsCosine)

"""
    Save a run into our data file
"""
function SaveData(data::Union{Spheroid, Simulation}, filename::String)
    @assert !isfile("data\\$filename.jld2") "A file already exists with that name"
    save("data\\$filename.jld2", "s", data)
end
