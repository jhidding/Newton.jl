# Welcome to Newton.jl!

Newton.jl computes orbits of particles under the force of Gravity.

## Tutorial

```@eval
using Unitful
using CairoMakie
using Newton: set_still!, random_particles, position, run_simulation

function plot_n_particles_stabilized(n; seed=nothing)
    fig = Figure()
    ax = Axis3(fig[1, 1])

    orbits = run_simulation(
        set_still!(random_particles(n, seed=seed)), 1.0u"s", 5000)
    for i in 1:n
        orbit = position.(getindex.(orbits, i))
        lines!(ax, orbit / u"m")
	    scatter!(ax, orbit[end] / u"m")
    end
	save("three-particle-orbit.svg", fig)
end

plot_n_particles_stabilized(3, seed=4)

nothing
```

![](three-particle-orbit.svg)

## API

```@autodocs
Modules = [Newton]
```
