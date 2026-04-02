# Welcome to Newton.jl!

Newton.jl computes orbits of particles under the force of Gravity.

## Tutorial

The `Newton.jl` package provides a framework to compute orbits of a set of particles under the force of gravity.
A particle is a structure with `mass`, `position` and `momentum`.

```@example 1
using Unitful
using CairoMakie
using Printf
using Newton: set_still!, random_particles, position, mass, velocity, momentum, run_simulation
```

We provide a way to generate a set of random particles.

```@example 1
ps = random_particles(3, seed=4)

stat(p) = @sprintf("mass: %f,    pos: %5.3f %5.3f %5.3f,    vel: %5.3f %5.3f %5.3f",
  mass(p), position(p)..., velocity(p)...)
  
for (i, p) in enumerate(ps)
  println("Particle $(i): $(stat(p))")
end

@printf "Total momentum: %f %f %f" (momentum(ps)...)
```

We can set these still without changing their relative velocities:

```@example 1
set_still!(ps)
@printf "Total momentum: %f %f %f" (momentum(ps)...)
```

Let's compute their orbits!

```@example 1
orbits = run_simulation(ps, 1.0u"s", 5000)
nothing
```

And plot them,

```@example 1
fig = Figure()
ax = Axis3(fig[1, 1])

for i in 1:3
  orbit = position.(getindex.(orbits, i))
  lines!(ax, orbit / u"m")
  scatter!(ax, orbit[end] / u"m")
end

save("three-particle-orbit.svg", fig)
nothing
```

![](three-particle-orbit.svg)

## API

```@autodocs
Modules = [Newton]
```
