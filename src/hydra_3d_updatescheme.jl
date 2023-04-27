"""
    Use an energy minimisation to find the steady-state mechanical configuration given the vector Psivals, defaulting to the value stored in Sphere, use the currnt configuration of sphere as the initial guess of the position. 
"""
function UpdateShape(Sphere::Spheroid, Psivals::Vector{Float64}=Sphere.Psivals; Loud::Bool=true, V0=V0target)
    t0 = time()
    xvec= Nodes2xvec(Sphere.Nodes)
    optimresult = optimize(
        x -> EnergyFunctional(x, Sphere.Triangles, Psivals, V0=V0, Loud=false),
        xvec, 
        BFGS(linesearch=LineSearches.BackTracking()), autodiff=:forward,
        Optim.Options(iterations=Int(1e8), time_limit=300, show_trace=false)
    )
    !Optim.converged(optimresult) && throw("Optimisation failed to converge")
    println("Optimisation took $(round(optimresult.time_run, digits=2)) seconds and $(optimresult.iterations) iterations.")

    xvec = optimresult.minimizer
    Loud && EnergyFunctional(xvec, Sphere.Triangles, Psivals, Loud=true)
    Nodes= xvec2Nodes(xvec)
    Sphere = Spheroid(Nodes, Sphere.Triangles, Psivals)
    Loud && display(ShowPsi(Sphere))
    Loud && println("Updating the shape took $(time() - t0) seconds")
    return Sphere
end

# struct MinimisationSimulation
#     Spherevec::Vector{Spheroid}
#     tvec::Vector{Float64}
#     RunTime::Float64
# end
# ODESimulation(Sim::MinimisationSimulation) = ODESimulation(Sim.Spherevec, Sim.tvec, Sim.RunTime)

function RunSimulation(InitialSphere::Spheroid, InitialPsivals::Vector{Float64}, T::Int, Δt::Float64; V0grad::Float64=0.0, V0rest::Float64=V0target, V0thresh::Float64 = V0rest+0.5)
    t0 = time()
    @assert length(InitialPsivals) == length(InitialSphere.Psivals) "Psivals is the incorrect length"
    Psivec = copy(InitialPsivals)
    Spherevec = [Spheroid(InitialSphere.Nodes, InitialSphere.Triangles, Psivec)]
    V0vec = zeros(T+1)
    V0vec[1] = V0rest

    for t in 1:T
        (mod(t, 10)==0) && println("timestep $t of $T")
        Psivec = max.(Psivec + Δt * GradientPsi(Psivec, Spherevec[end], Sphere0), 0.0)
        Psivec = min.(Psivec, 1.0)
        # @show typeof(Spherevec[end])
        Sphere = UpdateShape(Spherevec[end], Psivec, V0=V0vec[t], Loud=false)
        push!(Spherevec, Sphere)
        V0vec[t+1] = V0vec[t] > V0thresh ? V0rest : V0vec[t] + Δt*V0grad
    end
    # display(Show(Spherevec[end]))
    trun = (time()-t0)
    println("Running simulation took $trun seconds")
    return Simulation(Spherevec, Δt*collect(0:T), V0vec, trun)
end

RunSimulation(InitialSphere::Spheroid, T::Int, Δt::Float64; kwargs...) = RunSimulation(InitialSphere, InitialSphere.Psivals, T, Δt; kwargs...)


