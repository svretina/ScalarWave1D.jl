module ODE1D

using OrdinaryDiffEq: throwex
using LinearAlgebra
using StaticArrays
using HDF5
using OrdinaryDiffEq

const proj_path = pkgdir(ODE1D)
if occursin("cn", gethostname()) || occursin("mn", gethostname())
    const output_path = "/gpfs/svretinaris/ScalarWave/runs/"
else
    const output_path = string(proj_path, "/output/")
end

@show output_path

function sims_path(folder)
    return string(output_path, folder, "/sims/")
end

function get_iter_str(i, N)
    i = string(i)
    while length(i) < length(string(N))
        i = string("0", i)
    end
    if length(i) <= 2
        return string("0", i)
    else
        return i
    end
end

function rhs!(du::AbstractArray{<:Real},
              U::AbstractArray{<:Real},
              params,
              t::Real)
    h = params.h
    N = params.N
    boundary_type = params.bc
    @fastmath @inbounds begin
        h2 = 2h
        N1 = N - 1
        Π = @view U[:, 2]
        Ψx = @view U[:, 3]

        dtΦ = @view du[:, 1]
        dtΠ = @view du[:, 2]
        dtΨx = @view du[:, 3]

        for i in 1:N
            #calculate RHS everywhere
            if i == 1
                dxΨx = (Ψx[2] - Ψx[1]) / h
                dxΠ = (Π[2] - Π[1]) / h
            elseif i == N
                dxΨx = (Ψx[end] - Ψx[N1]) / h
                dxΠ = (Π[end] - Π[N1]) / h
            else
                dxΠ = (Π[i + 1] - Π[i - 1]) / h2
                dxΨx = (Ψx[i + 1] - Ψx[i - 1]) / h2
            end
            dtΦ[i] = Π[i]
            # Apply boundary conditions
            # absorbing/radiative
            if i == 1
                if boundary_type === :radiative
                    ## Variable change to characteristic vars
                    a = dxΨx + dxΠ # speed = -1
                    b = dxΨx - dxΠ # speed = +1
                    ## Set incoming characteristic to 0
                    b = 0
                    ## Change back to primitive variables
                    dxΨx = (a + b) / 2
                    dxΠ = (a - b) / 2
                    ## ODE
                    dtΠ[i] = dxΨx
                    dtΨx[i] = dxΠ
                elseif boundary_type === :reflective
                    ## Set ingoing to outgoing characteristic
                    dtΠ[i] = 0 # dxΨx
                    dtΨx[i] = dxΠ
                else
                    throw("boundary_type variable can only be :radiative or :reflective, you provided
                           $boundary_type")
                end
            elseif i == N
                if boundary_type === :radiative
                    ## Variable change to characteristic vars
                    a = dxΨx + dxΠ
                    b = dxΨx - dxΠ
                    ## Set incoming characteristic to 0
                    a = 0
                    ## Change back to primitive variables
                    dxΨx = (a + b) / 2
                    dxΠ = (a - b) / 2
                    # ODE
                    dtΠ[i] = dxΨx
                    dtΨx[i] = dxΠ
                elseif boundary_type === :reflective
                    dtΠ[i] = 0 # dxΨx
                    dtΨx[i] = dxΠ
                else
                    throw("boundary_type variable can only be :radiative or :reflective, you provided
                           $boundary_type")
                end
            else
                dtΠ[i] = dxΨx
                dtΨx[i] = dxΠ
            end
        end
    end
    return nothing
end

end
