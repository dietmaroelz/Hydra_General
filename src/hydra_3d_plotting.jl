"""
    Plot the concentration Ψ across the spheroid. Color each face by averaging Ψ on the vertices. 
"""
function ShowPsi(spher::Spheroid; col=cgrad(:thermal))
    Points = Point3[Nodes2Points(spher.Nodes)...]
    Triangles = spher.Triangles

    facepsivals = spher.Psivals
    # facekappavals = kappa.(facepsivals)

    fig, ax, pl = poly(Points, Triangles[1], color=get(col, facepsivals[1]))
    for (i, T) in enumerate(Triangles)
        i==1 && continue
        poly!(Points, T, color=get(col, facepsivals[i]))
    end

    scatter!(spher.Nodes, color=:black)
    
    combinedNodes = spher.Nodes
    
    trix, triy, triz = Float64[], Float64[], Float64[]
    for T in Triangles
        append!(trix, [[combinedNodes[t,1] for t in T]; [combinedNodes[T[1], 1], NaN]])
        append!(triy, [[combinedNodes[t,2] for t in T]; [combinedNodes[T[1], 2], NaN]])
        append!(triz, [[combinedNodes[t,3] for t in T]; [combinedNodes[T[1], 3], NaN]])
    end
    lines!(trix,triy,triz)
    fig2, ax2, sc2 = scatter(rand(length(spher.Psivals)), color=spher.Psivals)
    Colorbar(fig[1,2], sc2, label="Ψ")
    return fig
end


"""
    Assumes that the shape is centred at the origin, also that there are no concavities with respect to the radial direction 
"""
function PlotRadius(Nodes::Matrix{Float64}; name::String = "Radius for each node")
    θvals, ϕvals, Rvals = Float64[], Float64[], Float64[]
    angles = Point2[]
    for n in 1:size(Nodes)[1]
        x, y, z = Nodes[n,:]
        push!(Rvals, sqrt(x^2+y^2+z^2))
        push!(θvals, acos(z/Rvals[end]))
        push!(angles, Point2([acos(z/Rvals[end]), sign(y)*acos(x/sqrt(x^2+y^2))]))
        push!(ϕvals, sign(y)*acos(x/sqrt(x^2+y^2)))
    end

    fig, ax, sc = scatter(angles, color=Rvals)

    ax.xlabel = "θ"
    ax.ylabel = "ϕ"
    ax.title = name
    Colorbar(fig[1,2], sc)
    return fig
end
PlotRadius(Spher::Spheroid; args...) = PlotRadius(Spher.Nodes; args...)

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
StretchAgainstKappa(Sim::ODESimulation, j::Int) = StretchAgainstKappa(Sim.Spherevec[1], Sim.Spherevec[j])

function StretchAgainstPsi(Spher0::Spheroid, Spher1::Spheroid)
    Evals = StretchingEnergies(Spher0, Spher1)
    fig, ax, sc = scatter(Spher1.Psivals, Evals)
    ax.xlabel = "Ψ"
    ax.ylabel = "Stretch energy"
    ax.title = "Stretch energy per triangle against Ψ"
    return fig
end
StretchAgainstPsi(Sim::ODESimulation, j::Int) = StretchAgainstPsi(Sim.Spherevec[1], Sim.Spherevec[j])


function LoessPlot(xvals, yvals)
    model = loess(xvals, yvals)

    uvals = range(extrema(xvals)...; length=100)
    vvals = predict(model, uvals)
    return uvals, vvals
end

function MeanPsi(Sim::ODESimulation)
    meanPsivals = mean(S.Psivals for S in Sim.Spherevec)
    return meanPsivals
    return plot(meanPsivals)
end

function StretchAgainstKappaVideo(Sim::ODESimulation; name::String="Stretch_kappa")
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

function StretchAgainstKappaTrajectories(Sim::ODESimulation)
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


function StretchAgainstPsiVideo(Sim::ODESimulation; name::String="Stretch_psi")
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

function StretchAgainstPsiTrajectories(Sim::ODESimulation)
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

function AnimateSimulation(Sim::ODESimulation; col=cgrad(:thermal), PlotEvery::Int=20, name::String="SimVideo")
    sphere = Sim.Spherevec[1]
    Nodes = sphere.Nodes
    Points = Observable(Point3[Nodes2Points(Nodes)...])
    # Psivals = Observable(spher.Psivals)
    Psivals = sphere.Psivals
    Triangles = sphere.Triangles


    
    fig, ax, plt = poly(Points, Triangles[1], color=get(col, Psivals[1]))
    pltvec = [plt]
    for (i, T) in enumerate(Triangles)
        i==1 && continue
        push!(pltvec, poly!(Points, T, color=get(col, Psivals[i])))
    end

    record(fig, name*".mp4", 1:PlotEvery:length(Sim), framerate=2) do t
        Points[] = Point3[Nodes2Points(Sim.Spherevec[t].Nodes)...]
        for i in 1:length(Triangles)
            pltvec[i].color = get(col, Sim.Spherevec[t].Psivals[i])
        end
        # Psivals[] .+= rand(length(Psivals[]))
    end
    
    return fig
end

function RadiusAgainstTime(Sim::ODESimulation)
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

function Volumes(Sim::ODESimulation; V0target::Float64=V0target)
    Vvals = SphereVolume.(Sim.Spherevec)
    fig, ax, ln =  lines(Sim.tvec, Vvals)
    lines!([Sim.tvec[1], Sim.tvec[end]], [V0target, V0target])
    ax.title = "Volume over time"
    ax.xlabel = "Time"
    ax.ylabel = "Volume"
    return fig
end