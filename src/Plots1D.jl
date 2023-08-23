module Plots1D

import ..Utils
import ..Domains
import ..Grids
import ..GridFunctions
import ..InitialData
import ..Convergence
import ..ODE

using LaTeXStrings
using Plots
using Printf
using Grep
using HDF5
using Base.Threads

###  convert -delay 10 -loop 0 `ls -v Φsurface*` surface.gif
const proj_path = pkgdir(Plots1D)
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
function figure_path(folder)
    return string(output_path, folder, "/figures/")
end

function plot_over_time(sol, t=0, save=true)
    @show sol[1].coords
    idx = findall(x -> t - sol[1].spacing < x < t + sol[1].spacing, sol[1].coords)[1]
    @show idx, sol[1].coords[idx]
    for i in idx:(sol[1].npoints)
        @show i
        p1 = plot(sol[2].coords, sol[3][:, 1, i]; lw=1,
                  marker=(:circle, 2, 0.5))
        title!(p1, "t=$(sol[1].coords[i])")
        xlabel!(p1, "x")
        ylabel!(p1, "Φ")
        if save
            figure_path = string(proj_path, "/figures/")
            savefig(p1, string(figure_path, "Φ_1D_i=$(i).png"))
        end
    end
    return nothing
end

function plot_energy(sol; save=true)
    sz = size(sol[3])
    N = sz[begin]
    nt = sz[end]
    n0 = length(sol[1].coords) - 1
    n1 = nt - 1
    step = Int64(n0 / n1)
    @views begin
        Π = sol[3][:, 2, :]
        Ψ = sol[3][:, 3, :]
    end
    val = zeros(nt)
    for t in 1:nt
        @fastmath @inbounds for ci in CartesianIndices((1:N))
            i = ci[1]
            if i == 1 || i == N
                val[t] += (Π[i, t]^2 +
                           Ψ[i, t]^2) / 2
            else
                val[t] += Π[i, t]^2 + Ψ[i, t]^2
            end
        end
    end
    @show size(sol[1].coords)
    @show size(val)
    @show step
    p1 = plot(sol[1].coords[begin:step:end],
              val ./ val[1];
              lw=1,
              marker=(:circle, 2, 0.5),
              xlabel="time",
              ylabel="energy",
              tex_output_standalone=true,
              title=L"E = \displaystyle \int \Pi^2 +\sum_i \Psi^2_i dV / E_0",
              leg=false)

    figure_path = string(proj_path, "/figures/")
    @show figure_path
    if save
        savefig(p1, string(figure_path, "energy1D.tex"))
        savefig(p1, string(figure_path, "energy1D.pdf"))
        savefig(p1, string(figure_path, "energy1D.svg"))
    end
    return nothing
end

function plot_error_unknown(time_grid, grid, errors; save=true)
    depth = size(errors)[3]
    n0 = length(grid) - 1
    p1 = plot(; title=L"\rm Error")
    for i in 1:(depth)
        mul = 2^(i - 1)
        n = mul * n0
        p1 = plot!(time_grid,
                   errors[1, 2, i, :];
                   legend=:outerbottom,
                   legendcolumns=depth,
                   tex_output_standalone=true)
    end
    fig_path = string(proj_path, "/figures/")
    if save
        savefig(p1, string(fig_path, "error.tex"))
        savefig(p1, string(fig_path, "error.pdf"))
        savefig(p1, string(fig_path, "error.svg"))
    end
end

function plot_convorder(time_grid, grid, errors; save=true)
    depth = size(errors)[3]
    n0 = length(grid) - 1
    p1 = plot(; title=L"\rm Convergence \,\, order")

    for i in 2:(depth)
        mul = 2^(i - 1)
        n = mul * n0
        p = @. log2(errors[1, 2, i - 1, :] / errors[1, 2, i, :])
        p1 = plot!(time_grid,
                   p;
                   label="Ncells = $n",
                   leg=false,
                   tex_output_standalone=true,
                   xlabel="time",
                   ylabel=L"p")
    end
    ylims!(p1, -2, 4)
    fig_path = string(proj_path, "/figures/")
    if save
        savefig(p1, string(fig_path, "p.tex"))
        savefig(p1, string(fig_path, "p.pdf"))
        savefig(p1, string(fig_path, "p.svg"))
    end
end

end
