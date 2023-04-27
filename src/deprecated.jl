
function StretchAgainstKappaVideo(Sim::Simulation; name::String="Stretch_kappa")
    Spher0 = Sim.Spherevec[1]
    Evals = Observable(StretchingEnergies(Spher0, Spher0))
    Kappavals = Observable(kappa.(Spher0.Psivals))
    fig, ax, sc = scatter(Kappavals, Evals)
    xlims!(ax, [0, κ0])
    ax.xlabel = "κ"
    ax.ylabel = "Deformation"
    ax.title = "Deformation against κ per mesh triangle over time"
    record(fig, name*".mp4", 1:length(Sim), framerate=10) do t
        Evals[] = StretchingEnergies(Sim.Spherevec[t], Spher0)
        Kappavals[] = kappa.(Sim.Spherevec[t].Psivals)
        ylims!(ax, [0,maximum(Evals[])*1.2])
    end
end

function StretchAgainstKappaTrajectories(Sim::Simulation)
    Spher0 = Sim.Spherevec[1]
    DeformationVals = zeros(length(Spher0.Triangles), length(Sim))
    Kappavals = zeros(length(Spher0.Triangles), length(Sim))

    for t in 1:length(Sim)
        Evals = StretchingEnergies(Sim.Spherevec[t], Spher0)
        kappalocal = kappa.(Sim.Spherevec[t].Psivals)
        for i in 1:length(Spher0.Triangles)
            DeformationVals[i,t] = Evals[i]>0 ? Evals[i] : NaN
            Kappavals[i,t] = kappalocal[i]
        end
    end
    fig, ax, ln = lines(Kappavals[1,:], DeformationVals[1,:])
    ax.yscale=log10
    ax.xscale=log10
    for i in 2:length(Spher0.Triangles)
        lines!(Kappavals[i,:], DeformationVals[i,:])
    end
    xlims!(ax, [0.1, κ0])
    ax.xlabel = "κ"
    ax.ylabel = "Deformation"
    ax.title = "Deformation against κ per mesh triangle over time"
    return fig
end


function StretchAgainstPsiVideo(Sim::Simulation; name::String="Stretch_psi")
    Spher0 = Sim.Spherevec[1]
    Evals = Observable(StretchingEnergies(Spher0, Spher0))
    Psivals = Observable(Spher0.Psivals)
    fig, ax, sc = scatter(Psivals, Evals)
    xlims!(ax, [0, 1])
    ax.xlabel = "Ψ"
    ax.ylabel = "Deformation"
    ax.title = "Deformation against Ψ per mesh triangle over time"
    record(fig, name*".mp4", 1:length(Sim), framerate=10) do t
        Evals[] = StretchingEnergies(Sim.Spherevec[t], Spher0)
        Psivals[] = Sim.Spherevec[t].Psivals
        ylims!(ax, [0,maximum(Evals[])*1.2])
    end
end

function StretchAgainstPsiTrajectories(Sim::Simulation)
    Spher0 = Sim.Spherevec[1]
    DeformationVals = zeros(length(Spher0.Triangles), length(Sim))
    Psivals = zeros(length(Spher0.Triangles), length(Sim))

    for t in 1:length(Sim)
        Evals = StretchingEnergies(Sim.Spherevec[t], Spher0)
        psilocal = Sim.Spherevec[t].Psivals
        for i in 1:length(Spher0.Triangles)
            DeformationVals[i,t] = Evals[i]
            Psivals[i,t] = psilocal[i]
        end
    end
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Ψ", ylabel="Deformation", title="Deformation against Ψ per mesh triangle over time")#, xscale=log10, yscale=log10)
    lines!(ax, Psivals[1,2:end], DeformationVals[1,2:end])
    for i in 2:length(Spher0.Triangles)
        lines!(ax, Psivals[i,2:end], DeformationVals[i,2:end])
    end
    return fig
end

"""
    Plot the stretch relative to Spher0 of Spher1 against the kappa values associated with Spher1
"""
function StretchAgainstKappa(Spher0::Spheroid, Spher1::Spheroid)
    Evals = StretchingEnergies(Spher0, Spher1)
    Kappavals = kappa.(Spher1.Psivals)
    fig, ax, sc = scatter(Kappavals, Evals)
    ax.xlabel = "κ"
    ax.ylabel = "Stretch energy"
    ax.title = "Stretch energy per triangle against κ"
    return fig
end
StretchAgainstKappa(Sim::Simulation, j::Int) = StretchAgainstKappa(Sim.Spherevec[1], Sim.Spherevec[j])

function StretchAgainstPsi(Spher0::Spheroid, Spher1::Spheroid)
    Evals = StretchingEnergies(Spher0, Spher1)
    fig, ax, sc = scatter(Spher1.Psivals, Evals)
    ax.xlabel = "Ψ"
    ax.ylabel = "Stretch energy"
    ax.title = "Stretch energy per triangle against Ψ"
    return fig
end
StretchAgainstPsi(Sim::Simulation, j::Int) = StretchAgainstPsi(Sim.Spherevec[1], Sim.Spherevec[j])


function RadiusAgainstTime(Sim::Simulation)
    Rvals = Float64[]
    for Spher in Sim.Spherevec
        Rvalslocal = Float64[]
        Nodes = Spher.Nodes
        for n in 1:size(Nodes)[1]
            x, y, z = Nodes[n,:]
            push!(Rvalslocal, sqrt(x^2+y^2+z^2))
        end
        push!(Rvals, mean(Rvalslocal))
    end
    # return angles
    fig, ax, sc = lines(Sim.tvec, Rvals)

    ax.xlabel = "t"
    ax.ylabel = "R"
    ax.title = "Average Radius Over Time"
    return fig
end