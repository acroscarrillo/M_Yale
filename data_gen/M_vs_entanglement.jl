include("../src/src.jl") 
using Plots
using ProgressBars
using DataFrames
using CSV
using LaTeXStrings

O_basis = O_basis_spin_half(4)

zero = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
one = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]

O_basis = O_basis_spin_half(2)

zero = [1,0,0,0]
one = [0,0,0,1]

θ_array = Vector( range(0,1, length=100) )
ϕ_array = Vector( range(0,2, length=100) )

M_array = zeros(length(θ_array),length(ϕ_array))
for (j,θ) in ProgressBar(enumerate(θ_array))
    for (k,ϕ) in enumerate(ϕ_array)
        ψ = cos(θ*π/2)*zero + exp(im*ϕ*π)*sin(θ*π/2)*one
        ρ = ψ*ψ'
        M_array[j,k] = M(ρ, O_basis)
    end
end
heatmap(M_array',xlab=L"\theta/\pi", ylab=L"\phi/\pi",title=L"M\left(\cos(\theta/2)|0000\rangle + \exp(i\phi)\sin(θ/2)|1111\rangle \right)")


entanglement_array = zeros(length(θ_array),length(ϕ_array))
for (j,θ) in ProgressBar(enumerate(θ_array))
    for (k,ϕ) in enumerate(ϕ_array)
        ψ = cos(θ*π/2)*zero + exp(im*ϕ*π)*sin(θ*π/2)*one
        ρ = ψ*ψ'
        entanglement_array[j,k] = entanglement_entropy(ρ, [0,1])
    end
end
heatmap(entanglement_array',xlab=L"\theta/\pi", ylab=L"\phi/\pi",title=L"\mathcal{E}\left(\cos(\theta/2)|0000\rangle + \exp(i\phi)\sin(θ/2)|1111\rangle \right)")




function M_f(θ,ϕ)
    j = findfirst(x->x==θ,θ_array)
    k = findfirst(x->x==ϕ,ϕ_array)
    return M_array[j,k]
end

M_plot = contourf(θ_array,ϕ_array,M_f,xlab=L"\theta/\pi", ylab=L"\phi/\pi",title=L"M(\theta,\phi)",dpi=600,titlefontsize=14,guidefont=14)

function entang_f(θ,ϕ)
    j = findfirst(x->x==θ,θ_array)
    k = findfirst(x->x==ϕ,ϕ_array)
    return entanglement_array[j,k]
end
E_plot = contourf(θ_array,ϕ_array,entang_f,xlab=L"\theta/\pi", ylab=L"\phi/\pi",title=L"\mathcal{E}(\theta,\phi)",dpi=600,titlefontsize=14,guidefont=14)



l = @layout [M_plot E_plot] 

plot(M_plot, E_plot, layout = l,size=(600,250),dpi=650)
savefig("figs/M_vs_E.pdf")
