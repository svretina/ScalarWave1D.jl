module Convergence1D

import ..Domains
import ..GridFunctions
import ..Grids
import ..Run1D
import ..ODE1D

using LinearAlgebra
using OrdinaryDiffEq
using StaticArrays
using HDF5
using Grep
## Sol of ODE -> (spatial_Grid, state_vector, time_Grid)

const proj_path = pkgdir(Convergence1D)
if occursin("cn", gethostname()) || occursin("mn", gethostname())
    const output_path = "/gpfs/svretinaris/ScalarWave/runs/"
    # const error_path = "/gpfs/svretinaris/ScalarWave/errors/"
    # const conv_path = "/gpfs/svretinaris/ScalarWave/convergence/"
else
    const output_path = string(proj_path, "/output/")
end
function error_path(folder)
    return string(output_path, folder, "/errors/")
end
function sims_path(folder)
    return string(output_path, folder, "/sims/")
end
function convergence_path(folder)
    return string(output_path, folder, "/convergence/")
end

function run_and_downsample_true(depth::Integer, tf::Real)
    T = Float64
    L = 10
    nc0 = 10
    dx = 2
    cfl = 1 / 10
    dt = cfl * dx
    nt = floor(Int64, tf / dt) + 1
    nx = nc0 + 1
    paramsID = SVector{4,Int64}(2, 1, 1, L)

    # allocate stuff
    time_grids = Array{T,2}(undef, nt, depth + 1)
    grids = Array{T,2}(undef, nx, depth + 1)
    sims = Array{T,6}(undef, nx, nx, nx, 5, depth + 1, nt)
    truesol = Array{T,5}(undef, nx, nx, nx, 5, nt)
    errors = Array{T,4}(undef, 1, 5, depth + 1, nt)

    @inbounds for i in 1:(depth + 1)
        multiplier = 2^(i - 1)
        ncells = multiplier * nc0
        @show ncells
        simi = Run1D.cartesian(L, ncells, tf, cfl)
        time_grids[:, i] = simi[1].coords[1:multiplier:end]
        grids[:, i] = simi[2].x[1:multiplier:end]
        sims[:, :, :, :, i, :] = simi[3][1:multiplier:end,
                                         1:multiplier:end,
                                         1:multiplier:end,
                                         :,
                                         1:multiplier:end]
    end
    t = time_grids[:, 1]
    g = grids[:, 1]
    @inbounds for i in 1:nt
        statevector_ti = @view truesol[:, :, :, :, i]
        Run.sample_statevector!(statevector_ti, g, paramsID, t[i])
    end
    # check initial data are consistent
    @inbounds for i in 1:(depth + 1)
        @assert all(abs.(sims[:, :, :, :, i, 1] .- truesol[:, :, :, :, 1]) .== zero(T))
    end

    @inbounds for j in 1:nt, i in 1:(depth + 1), k in 1:5
        errors[1, k, i, j] = norm(sims[:, :, :, k, i, j] .- truesol[:, :, :, k, j])
    end
    return (t, g, sims, errors)
end

function galerkin_run_and_downsample_unknown(depth::Integer, tf::Real)
    @assert depth >= 2
    T = Float64
    L = 20
    nc0 = 320
    dx = 2L / nc0
    @show dx
    cfl = 1 / 10
    dt = cfl * dx
    nt = floor(Int64, tf / dt) + 1
    nx = nc0 + 1
    # paramsID = SVector{4,Int64}(2, 1, 1, L)
    @show nt
    # allocate stuff
    time_grids = Array{T,2}(undef, nt, depth + 1)
    grids = Array{T,2}(undef, nx, depth + 1)
    sims = Array{T,4}(undef, nx, 3, depth + 1, nt)
    truesol = Array{T,3}(undef, nx, 3, nt)
    errors = Array{T,4}(undef, 1, 3, depth, nt)

    @inbounds for i in 1:(depth + 1)
        multiplier = 2^(i - 1)
        ncells = multiplier * nc0
        @show ncells
        simi = Run1D.cartesian(L, ncells, tf, cfl, :reflective)
        time_grids[:, i] = simi[1].coords[1:multiplier:end]
        grids[:, i] = simi[2].coords[1:multiplier:end]
        sims[:, :, i, :] = simi[3][1:multiplier:end,
                                   :,
                                   1:multiplier:end]
    end
    t = time_grids[:, 1]
    g = grids[:, 1]
    # @inbounds for i in 1:nt
    #     statevector_ti = @view truesol[:, :, :, :, i]
    #     Run1D.sample_statevector!(statevector_ti, g, paramsID, t[i])
    # end
    # # check initial data are consistent
    # @inbounds for i in 1:(depth + 1)
    #     @assert all(abs.(sims[:, :, :, :, i, 1] .- truesol[:, :, :, :, 1]) .== zero(T))
    # end

    @inbounds for j in 1:nt, i in 1:(depth), k in 1:3
        errors[1, k, i, j] = norm(sims[:, k, i, j] .- sims[:, k, i + 1, j])
    end
    @show g
    return (t, g, sims, errors)
end

function run_and_downsample_unknown(depth::Integer, tf::Real)
    @assert depth >= 2
    T = Float64
    L = 20
    nc0 = 320
    dx = 2L / nc0
    @show dx
    cfl = 1 / 10
    dt = cfl * dx
    nt = floor(Int64, tf / dt) + 1
    nx = nc0 + 1
    # paramsID = SVector{4,Int64}(2, 1, 1, L)
    @show nt
    # allocate stuff
    time_grids = Array{T,2}(undef, nt, depth + 1)
    grids = Array{T,2}(undef, nx, depth + 1)
    sims = Array{T,4}(undef, nx, 3, depth + 1, nt)
    truesol = Array{T,3}(undef, nx, 3, nt)
    errors = Array{T,4}(undef, 1, 3, depth, nt)

    @inbounds for i in 1:(depth + 1)
        multiplier = 2^(i - 1)
        ncells = multiplier * nc0
        @show ncells
        simi = Run1D.cartesian(L, ncells, tf, cfl, :reflective)
        time_grids[:, i] = simi[1].coords[1:multiplier:end]
        grids[:, i] = simi[2].coords[1:multiplier:end]
        sims[:, :, i, :] = simi[3][1:multiplier:end,
                                   :,
                                   1:multiplier:end]
    end
    t = time_grids[:, 1]
    g = grids[:, 1]
    # @inbounds for i in 1:nt
    #     statevector_ti = @view truesol[:, :, :, :, i]
    #     Run1D.sample_statevector!(statevector_ti, g, paramsID, t[i])
    # end
    # # check initial data are consistent
    # @inbounds for i in 1:(depth + 1)
    #     @assert all(abs.(sims[:, :, :, :, i, 1] .- truesol[:, :, :, :, 1]) .== zero(T))
    # end

    @inbounds for j in 1:nt, i in 1:(depth), k in 1:3
        errors[1, k, i, j] = norm(sims[:, k, i, j] .- sims[:, k, i + 1, j])
    end
    @show g
    return (t, g, sims, errors)
end

end
