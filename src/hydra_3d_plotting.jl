"""
    Plot the concentration Ψ across the spheroid. Color each face by averaging Ψ on the vertices. 
"""
function Show(spher::Spheroid; col=cgrad(:thermal))
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
PlotRadius(Sim::Simulation, t::Int; kwargs...) = PlotRadius(Sim.Spherevec[t]; kwargs...)



function LoessPlot(xvals, yvals)
    model = loess(xvals, yvals)

    uvals = range(extrema(xvals)...; length=100)
    vvals = predict(model, uvals)
    return uvals, vvals
end

function MeanPsi(Sim::Simulation)
    meanPsivals = mean(S.Psivals for S in Sim.Spherevec)
    return meanPsivals
    return plot(meanPsivals)
end


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


xlims(ax::Axis) = @lift ($(ax.finallimits).origin[1], $(ax.finallimits).origin[1] + $(ax.finallimits).widths[1])
ylims(ax::Axis) = @lift ($(ax.finallimits).origin[2], $(ax.finallimits).origin[2] + $(ax.finallimits).widths[2])

xlims() = xlims(current_axis())
ylims() = ylims(current_axis())

function AnimateSimulation(Sim::Simulation; col=cgrad(:thermal), PlotInterval::Int=1, name::String="SimVideo")
    sphere = Sim.Spherevec[1]
    Nodes = sphere.Nodes
    Points = Observable(Point3[Nodes2Points(Nodes)...])
    # Psivals = Observable(spher.Psivals)
    Psivals = sphere.Psivals
    Triangles = sphere.Triangles

    fig = Figure()

    pltax = poly(fig[1:2,1], Points, Triangles[1], color=get(col, Psivals[1]))
    pltvec = [pltax.plot]
    for (i, T) in enumerate(Triangles)
        i==1 && continue
        push!(pltvec, poly!(Points, T, color=get(col, Psivals[i])))
    end
    scatter!(Points, color=:black)

    trix, triy, triz = Float64[], Float64[], Float64[]
    for T in Triangles
        append!(trix, [[sphere.Nodes[t,1] for t in T]; [sphere.Nodes[T[1], 1], NaN]])
        append!(triy, [[sphere.Nodes[t,2] for t in T]; [sphere.Nodes[T[1], 2], NaN]])
        append!(triz, [[sphere.Nodes[t,3] for t in T]; [sphere.Nodes[T[1], 3], NaN]])
    end
    trix, triy, triz = Observable.([trix, triy, triz])
    lines!(trix, triy, triz)


    Rvals = Float64[]
    angles = Point2[]
    for n in 1:size(Nodes)[1]
        x, y, z = Nodes[n,:]
        push!(Rvals, sqrt(x^2+y^2+z^2))
        push!(angles, Point2([acos(z/Rvals[end]), sign(y)*acos(x/sqrt(x^2+y^2))]))
    end
    angles = Observable(angles)

    sc = scatter(fig[1,2], angles, color=Rvals)
    axradius = sc.axis
    axradius.xlabel = "θ"
    axradius.ylabel = "ϕ"
    axradius.title = name
    Colorbar(fig[1,3], sc.plot)

    lines(fig[2,2], collect(0:length(Sim)-1), Sim.V0vec)
    trun = Observable([0.0, 0.0])
    ylims!(minimum(Sim.V0vec)-0.5, maximum(Sim.V0vec)+0.5)
    lines!(trun, Float64.([ylims()[]...]), color=:black, linestyle=:dash)
    axVol = fig.content[end]
    axVol.xlabel= "timestep"
    axVol.ylabel = "Preferred Volume"


    record(fig, "Animations/"*name*".mp4", 1:PlotInterval:length(Sim), framerate=2) do t
        sphere = Sim.Spherevec[t]
        Points[] = Point3[Nodes2Points(sphere.Nodes)...]
        for i in 1:length(Triangles)
            pltvec[i].color = get(col, sphere.Psivals[i])
        end

        trixtemp, triytemp, triztemp = Float64[], Float64[], Float64[]
        for T in Triangles
            append!(trixtemp, [[sphere.Nodes[t,1] for t in T]; [sphere.Nodes[T[1], 1], NaN]])
            append!(triytemp, [[sphere.Nodes[t,2] for t in T]; [sphere.Nodes[T[1], 2], NaN]])
            append!(triztemp, [[sphere.Nodes[t,3] for t in T]; [sphere.Nodes[T[1], 3], NaN]])
        end
        trix[], triy[], triz[] = trixtemp, triytemp, triztemp

        Rvals = Float64[]
        anglestemp = Point2[]
        for n in 1:size(sphere.Nodes)[1]
            x, y, z = sphere.Nodes[n,:]
            push!(Rvals, sqrt(x^2+y^2+z^2))
            push!(anglestemp, Point2([acos(z/Rvals[end]), sign(y)*acos(x/sqrt(x^2+y^2))]))
        end
        angles[] = anglestemp 
        sc.plot.color = Rvals

        trun[] = [t-2, t-2]
    end
    
    return fig
end


function Volumes(Sim::Simulation; V0target::Float64=V0target)
    Vvals = SphereVolume.(Sim.Spherevec)
    fig, ax, ln =  lines(Sim.tvec, Vvals)
    lines!([Sim.tvec[1], Sim.tvec[end]], [V0target, V0target])
    ax.title = "Volume over time"
    ax.xlabel = "Time"
    ax.ylabel = "Volume"
    return fig
end