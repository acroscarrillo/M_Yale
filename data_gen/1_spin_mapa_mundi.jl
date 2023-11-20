include("../src/src.jl") 
using Plots
using ProgressBars
using DataFrames
using CSV
using LaTeXStrings

O_basis = O_basis_spin_half(1)

zero = [1,0]
one = [0,1]

function M_one_spin(θ,ϕ)
    ψ = cos(θ*π/2)*zero + exp(im*ϕ*π)*sin(θ*π/2)*one
    ρ = ψ*ψ'
    return M(ρ, O_basis)
end


θ_array = Vector( range(0,π, length=999) )
ϕ_array = Vector( range(0,2*π, length=1000) )

contourf(θ_array./π,ϕ_array./π,M_one_spin,xlab=L"\theta/\pi", ylab=L"\phi/\pi",title=L"M(\theta,\phi) = M\left( \cos(\theta/2)|0\rangle + \exp(i\phi)\sin(θ/2)|1\rangle \right)",dpi=600)

savefig("figs/1_spin_mapa_mundi.png")