module Newton

using Unitful
using GeometryBasics
using LinearAlgebra
using Random

"""
Universal Gravitational Constant
"""
const G = 6.6743e-11u"m^3*kg^-1*s^-2"

"""
    gravitation_force(m1, m2, r)

Takes `r` to be the scalar distance between two objects of masses `m1` and `m2`.
Returns the strength of the force of gravitational attraction between the
two objects.
"""
gravitational_force(m1, m2, r) =
    G * m1 * m2 / r^2

"""
    gravitational_force(m1, m2, r::AbstractVector)

Takes `r` to be the distance vector between two objects of masses `m1` and `m2`.
Returns the gravitational force due to Newton's law in of the direction `r`.
"""
gravitational_force(m1, m2, r::AbstractVector) =
	r * (G * m1 * m2 * (r ⋅ r)^(-1.5))

"""
Type for masses in units of kilograms.
"""
const Mass = typeof(1.0u"kg")

"""
Type for a 3d momentum vector in units of Newton seconds.
"""
const MomentumVector = typeof(Vec3d(1)u"kg*m/s")

"""
Type for a 3d position in vector in units of meters.
"""
const PositionVector = typeof(Vec3d(1)u"m")

"""
Type for a 3d velocity vector in units of meters per second.
"""
const VelocityVector = typeof(Vec3d(1)u"m/s")

"""
    Particle(mass, position, momentum)

Particle structure. The `position` and `momentum` should be 3-vectors with the correct units,
and `mass` a scalar mass.
"""
mutable struct Particle
    mass::Mass
    position::PositionVector
    momentum::MomentumVector
end

mass(p::Particle) = p.mass
position(p::Particle) = p.position
momentum(p::Particle) = p.momentum

mass(p::AbstractArray{Particle}) = sum(mass, p)
momentum(p::AbstractArray{Particle}) = sum(momentum, p)
velocity(p) = momentum(p) / mass(p)

"""
    random_particle(mass=1e6"kg", spread=1.0u"m", dispersion=2.0u"mm/s")

Generate a particle with given `mass`, but random position and velocity.
The position and velocity are drawn from a normal distribution and scaled
with given `spread` and `dispersion`.

The default values are chosen to give a high probability for interesting
behaviour.
"""
random_particle(mass=1e6u"kg", spread=1.0u"m", dispersion=2.0u"mm/s") =
    Particle(mass, randn(Vec3d) * spread, randn(Vec3d) * dispersion * mass)

"""
    random_particles(n; seed=0, args...)

Generate `n` random particles with a given random seed. Extra keyword
arguments `args...` are forwarded to the `random_particle` function.
"""
function random_particles(n; seed=0, args...)
    Random.seed!(seed)
    [random_particle(args...) for _ in 1:n]
end

"""
    set_still!(particles)

Computes the net velocity of a set of particles, and changes the momentum
of each particle to match this frame of reference.

Returns the particle set.
"""
function set_still!(particles)
    v = velocity(particles)
    for p in particles
        p.momentum -= v * mass(p)
    end
    return particles
end

"""
    kick!(particles::AbstractVector{Particle}, dt)

Change the momentum of a each particle in the vector `particles`, following
direct one-to-one computation of their respective attractive forces.
"""
function kick!(particles, dt)
    for i in eachindex(particles)
        for j in 1:(i-1)
            r = particles[j].position - particles[i].position
            force = gravitational_force(particles[i].mass, particles[j].mass, r)
            particles[i].momentum += dt * force
            particles[j].momentum -= dt * force
        end
    end
    return particles
end

"""
    potential_energy(particles::AbstractVector{Particle})

Computes the potential energy of the system of particles.
"""
function potential_energy(particles)
    total = 0.0u"J"
    for i in eachindex(particles)
        for j in 1:(i-1)
            r = particles[j].position - particles[i].position
            m1 = particles[i].mass
            m2 = particles[j].mass
            total -= G * m1 * m2 / sqrt(r ⋅ r)
        end
    end
    return total
end

"""
    kinetic_energy(p::Particle)

Compute the kinetic energy of a particle.
"""
kinetic_energy(s::Particle) = let p = momentum(s)
    (p ⋅ p) / (2 * mass(s))
end

kinetic_energy(particles::AbstractVector{Particle}) =
    sum(kinetic_energy(p) for p in particles)

total_energy(particles::AbstractVector{Particle}) =
    potential_energy(particles) + kinetic_energy(particles)

"""
    drift!(p::Particle, dt)

Evolve the position of particle `p` for a time `dt` from its given momentum.
"""
function drift!(p::Particle, dt)
    p.position += dt * p.momentum / p.mass
end

"""
    drift!(particles: AbstractVector{Partcile}, dt)

Evolve the position of all particles.
"""
function drift!(particles, dt)
    for p in particles
        drift!(p, dt)
    end
    return particles
end

"""
    leap_frog!(particles, dt)

One leap-frog integration time step `dt` for `particles` under
gravity.
"""
function leap_frog!(particles, dt)
    drift!(particles, dt/2)
    kick!(particles, dt)
    drift!(particles, dt/2)
end

"""
    run_simulation(particles, dt, n_steps)

Leap-frog a set of particles `n_steps` times for a time step `dt`.
Copies the particle set after every iteration, returning a
vector containing the full state for each time step.
"""
function run_simulation(particles, dt, n_steps)
    x = deepcopy(particles)
    [deepcopy(leap_frog!(x, dt)) for _ in 1:n_steps]
end

"""
    random_orbits(n, mass; dt=1.0u"s", steps=5000, args...)

Generate random orbits of `n` particles with given `mass`.
"""
function random_orbits(n, mass; dt=1.0u"s", steps=5000, args...)
    particles = random_particles(n; args...)
    run_simulation(particles, dt, steps) |> collect
end

end  # module Newton
