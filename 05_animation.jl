using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom, PartialCellBottom
using Printf
using GLMakie
using Oceananigans.ImmersedBoundaries: mask_immersed_field!


arch = CPU()
#tracer_advection = CenteredSecondOrder()
#momentum_advection = CenteredSecondOrder()
momentum_advection = WENO5()
tracer_advection = WENO5()
#momentum_advection = CenteredFourthOrder()
#tracer_advection = CenteredFourthOrder()
#momentum_advection = UpwindBiasedFirstOrder()
#tracer_advection = UpwindBiasedFirstOrder()
#momentum_advection = UpwindBiasedThirdOrder()
#tracer_advection = UpwindBiasedThirdOrder()
#momentum_advection = UpwindBiasedFifthOrder()
#tracer_advection = UpwindBiasedFifthOrder()


underlying_grid = RectilinearGrid(arch,
                                  size=(160, 80), halo=(3, 3), 
                                  y = (-4, 4),
                                  z = (-1, 0),
                                  topology=(Flat, Periodic, Bounded))

# A bump
h₀ = 0.1 # bump height
L = 1 # bump width
@inline h(y) = h₀ * exp(- y^2 / L^2)
@inline seamount(x, y) = - 1 + h(y)

seamount_field = Field{Center, Center, Nothing}(underlying_grid)
set!(seamount_field, seamount)
fill_halo_regions!(seamount_field)

minimum_fractional_Δz = 0.2

immersed_boundaries = [
                       PartialCellBottom(seamount_field.data;
                                         minimum_fractional_Δz),
                       GridFittedBottom(seamount_field.data)
                      ]

b = []
v = []

function progress(sim)
    vmax = maximum(abs, sim.model.velocities.v)
    @info @sprintf("Iter: %d, time: %.2e, max|v|: %.2e",
                   iteration(sim), time(sim), vmax)

    return nothing
end

tracer_errors  = Dict()

for ib in immersed_boundaries
    grid = ImmersedBoundaryGrid(underlying_grid, ib)
    #grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(seamount_field.data))


    @show grid

    





### In this example of the paper the velocities were prescribe
### which means that they won't evolve during the simulation
### and using the function below we prescribe the velocities
#velocities = PrescribedVelocityFields( v=V, w=W)
  
  model = HydrostaticFreeSurfaceModel(; grid,
                                        tracer_advection,
                                        momentum_advection,
                                        #coriolis = FPlane(f=0.1),
                                        tracers = :b,
                                        #velocities = velocities,
                                        buoyancy = BuoyancyTracer())

B= Field{Center, Center, Center}(grid)   
N² = 1
                                  
set!(B, (x, y, z) -> (N² * z) ) 
mask_immersed_field!(B)

U, V, W = model.velocities

Ψ = Field{Center, Face, Face}(grid)

h(y)    = h₀*exp(-y^2/L^2)
ζ(y, z) = z/(h(y) - 1)
set!(Ψ, (x, y, z) -> (1 - ζ(y, z))^2)
fill_halo_regions!(Ψ, arch)


### We mask psi
mask_immersed_field!(Ψ)

### V is (Center, Face, Center)
V = YFaceField(grid)
### W is (Center, Centere, Face)
W = ZFaceField(grid)

### V and W are deravatives of psi 
### which in this code we use exact expression 
### which we computed in maple
### we are calling this method 'analytical'
V.=  ∂z(Ψ)
W.= -∂y(Ψ)

### We mask V and W        
mask_immersed_field!(V)
mask_immersed_field!(W)

### We fill the halo regions of V and W

fill_halo_regions!(V, arch)
fill_halo_regions!(W, arch)


  
    set!(model, b = B,  w=W, v=V)
    tracer_initial= sum(interior(model.tracers.b))

    simulation = Simulation(model; Δt=1e-3, stop_time=1)
    simulation.callbacks[:p] = Callback(progress, IterationInterval(10))



    filename = "flow_over_seamount_partial"

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers,
                                                      schedule = TimeInterval(0.01),
                                                      filename = filename * ".jld2",
                                                      overwrite_existing = true)

    run!(simulation)

    push!(b, Array(interior(model.tracers.b, 1, :, :)))
    push!(v, Array(interior(model.velocities.v, 1, :, :)))
    tracer_final = sum(interior(model.tracers.b))
    tracer_errors[ib] = ((tracer_initial - tracer_final)/tracer_initial)*100
    @info """
  
    Error in tracer conservation is $(tracer_errors[ib]) percent.
    """
end

Δb= []

grid_full = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(seamount_field.data))
ΔB= Field{Center, Center, Center}(grid_full)   
#set!(ΔB, (x, y, z) -> (0) ) 
mask_immersed_field!(ΔB, NaN)

Δb = ΔB[1,:,:]

b_partial = b[1]
b_full    = b[2]
Δb = b_full .- b_partial

v_partial = v[1]
v_full    = v[2]
Δv = v_full .- v_partial





fig = Figure(resolution=(1200, 1800))

partial_cell_title = @sprintf("PartialCellBottom with ϵ = %.1f", minimum_fractional_Δz)
ax_bp = Axis(fig[1, 2], title=partial_cell_title)
ax_bf = Axis(fig[2, 2], title="GridFittedBottom")
ax_bd = Axis(fig[3, 2], title="Difference (GridFitted - PartialCell)")

ax_vp = Axis(fig[1, 4], title=partial_cell_title)
ax_vf = Axis(fig[2, 4], title="GridFittedBottom")
ax_vd = Axis(fig[3, 4], title="Difference (GridFitted - PartialCell)")

color = (:black, 0.5)
linewidth = 3
levels = 15

#grid_partial = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(seamount_field.data,minimum_fractional_Δz))
#grid_full = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(seamount_field.data))

xb, yb, zb = nodes((Center, Center, Center), underlying_grid)
xv, yv, zv = nodes((Center, Face, Center), underlying_grid)


 
n = Observable(1)

ω = @lift interior(underlying_grid[$n], :, :, 1)
s = @lift interior(underlying_grid[$n], :, :, 1)


hmbp = heatmap!(ax_bp,yb, zb, b_partial, colorrange = (-1, 0))
contour!(ax_bp,yb, zb, b_partial; levels=-1:0.1:0, color, linewidth)
Colorbar(fig[1, 1], hmbp, label="Buoyancy", flipaxis=false)

hmbf = heatmap!(ax_bf, yb, zb, b_full, colorrange = (-1, 0))
contour!(ax_bf,yb, zb, b_full; levels=-1:0.1:1, color, linewidth)
Colorbar(fig[2, 1], hmbf, label="Buoyancy", flipaxis=false)

hmbd = heatmap!(ax_bd,yb, zb, Δb)
Colorbar(fig[3, 1], hmbd, label="Buoyancy", flipaxis=false)

hmvp = heatmap!(ax_vp, yv, zv, v_partial)
contour!(ax_vp, yv, zv, v_partial; levels, color, linewidth)
Colorbar(fig[1, 3], hmvp, label="Velocity", flipaxis=false)

hmvf = heatmap!(ax_vf,yv, zv, v_full,)
contour!(ax_vf, yv, zv, v_full; levels, color, linewidth)
Colorbar(fig[2, 3], hmvf, label="Velocity", flipaxis=false)

hmvd = heatmap!(ax_vd,yv, zv, Δv)
Colorbar(fig[3, 3], hmvd, label="Velocity", flipaxis=false)

display(fig)
title = @lift "t = " * string(round(times[$n], digits=2))
Label(fig[:, :], title, textsize=24, tellwidth=false)

frames = 1:length(times)

@info "Making a neat animation of vorticity and speed..."

record(fig, filename * ".mp4", frames, framerate=24) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")
    n[] = i
end

maximum(Δv)
minimum(Δv)
maximum(Δb)
minimum(Δb)
