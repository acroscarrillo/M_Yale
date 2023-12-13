include("../src/src.jl") 
using Plots
using ProgressBars
using DataFrames
using CSV
using LaTeXStrings


function ψ_test(L)
    zero = 1
    one = 1
    for n=1:L
        zero = zero ⊗ [1,0]
        one = one ⊗ [0,1]
    end
    return cos(π/4)*zero + exp(im*π/4)*sin(π/4)*one
end


L_array = Vector(1:1:5)

M_scaling_array = zeros(length(L_array))
for (j,L) in ProgressBar(enumerate(L_array))
    O_basis = O_basis_spin_half(L)
    ρ = ψ_test(L)*ψ_test(L)'
    M_scaling_array[j] = M(ρ, O_basis)
end
ttl=L"M( \cos(π/4)\bigotimes_{n=1}^L|0\rangle + \exp(iπ/4)\sin(π/4)\bigotimes_{n=1}^L|1\rangle )"
plot(L_array,M_scaling_array,yscale=:log2,legend=false,color=:black)
scatter!(L_array,M_scaling_array,color=:orange,ylabel=L"M",xlabel=L"L",title=ttl,dpi=600)
savefig("figs/M_scaling.png")
savefig("figs/M_scaling.pdf")




L_array = Vector(2:2:8)

E_scaling_array = zeros(length(L_array))
for (j,L) in ProgressBar(enumerate(L_array))
    ρ = ψ_test(L)*ψ_test(L)'
    E_scaling_array[j] = entanglement_entropy(ρ)
end
scatter(L_array,E_scaling_array)

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

contourf(θ_array,ϕ_array,M_f,xlab=L"\theta/\pi", ylab=L"\phi/\pi",title=L"M(\theta,\phi) = M\left( \cos(\theta/2)|00\rangle + \exp(i\phi)\sin(θ/2)|11\rangle \right)",dpi=600)

function entang_f(θ,ϕ)
    j = findfirst(x->x==θ,θ_array)
    k = findfirst(x->x==ϕ,ϕ_array)
    return entanglement_array[j,k]
end
contourf(θ_array,ϕ_array,entang_f,xlab=L"\theta/\pi", ylab=L"\phi/\pi",title=L"\mathcal{E}(\theta,\phi) = \mathcal{E}\left( \cos(\theta/2)|00\rangle + \exp(i\phi)\sin(θ/2)|11\rangle \right)",dpi=600)


savefig("figs/test.png")