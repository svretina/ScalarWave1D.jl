module Run1D

import ..Domains
import ..GridFunctions
import ..Grids
import ..InitialData
import ..ODE1D
import ..Galerkin

using LinearAlgebra
using OrdinaryDiffEq
#using StaticArrays

## Wave equation in 1D on a cartesian grid
function cartesian(L::Real, ncells::Integer, tf::Real, cfl::Real,
                   boundary_type::Symbol=:radiative)
    @assert boundary_type === :radiative || boundary_type === :reflective

    T = Float64
    dom = Domains.Domain([-L, L])
    grid = Grids.Grid(dom, ncells)
    N = grid.npoints

    domt = Domains.Domain([zero(tf), tf])
    nt = ceil(Int64, tf / (cfl * grid.spacing))
    t = Grids.Grid(domt, nt)

    # Initial Data
    Φ = InitialData.Gaussian1D.(0, grid.coords, 2)
    Π = InitialData.dtGaussian1D.(0, grid.coords, 2)
    Ψ = InitialData.dxGaussian1D.(0, grid.coords, 2)
    tspan = Tuple(T.(domt.domain))
    params = (h=grid.spacing, N=N, bc=boundary_type)
    statevector = hcat(Φ, Π, Ψ)

    save_every = 2
    ode = ODEProblem{true}(ODE1D.rhs!, statevector, tspan, params)

    sol = OrdinaryDiffEq.DiffEqBase.solve(ode,
                                          OrdinaryDiffEq.RK4();
                                          adaptive=false,
                                          dt=t.spacing)
    # VectorOfArray(sol.u)
    return (t, grid, sol)
end

## Wave equation in 1D on a cartesian grid but Initial Data is projected
## onto basis functions ( piecewise constant or linear )
function galerkin_wave1D(L::Real, ncells::Integer, tf::Real, cfl::Real,
                         boundary_type::Symbol=:radiative)
    @assert boundary_type === :radiative || boundary_type === :reflective

    T = Float64
    dom = Domains.Domain([-L, L])
    grid = Grids.Grid(dom, ncells)
    N = grid.npoints

    domt = Domains.Domain([zero(tf), tf])
    nt = ceil(Int64, tf / (cfl * grid.spacing))
    t = Grids.Grid(domt, nt)
    @show size(t.coords)
    # Initial Data
    Φ = Galerkin.coefs(x -> InitialData.Gaussian1D.(0, x, 2), grid)
    Π = Galerkin.coefs(x -> InitialData.dtGaussian1D.(0, x, 2), grid)
    Ψ = Galerkin.coefs(x -> InitialData.dxGaussian1D.(0, x, 2), grid)

    tspan = Tuple(T.(domt.domain))
    params = (h=grid.spacing, N=N, bc=boundary_type)

    statevector = hcat(Φ, Π, Ψ)

    save_every = 2
    ode = ODEProblem{true}(ODE1D.rhs!, statevector, tspan, params)

    sol = OrdinaryDiffEq.DiffEqBase.solve(ode,
                                          OrdinaryDiffEq.RK4();
                                          adaptive=false,
                                          dt=t.spacing)
    return (t, grid, sol)
end

end
