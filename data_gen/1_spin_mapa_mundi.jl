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

θ_array = Vector( range(0,1, length=997) )
ϕ_array = Vector( range(0,1, length=1001) )

M_max = 0
max_θ, max_ϕ = 0, 0
counter = 1
for θ in ProgressBar(θ_array)
    for ϕ in ϕ_array 
        M_temp = M_one_spin(θ,ϕ)
        if M_temp > M_max
            M_max = M_temp
            max_θ, max_ϕ = θ, ϕ
        end
    end
end

# hline!([0.25],color=:deepskyblue,linestyle=:dash,lw=1.5)
# hline!([0.75],color=:deepskyblue,linestyle=:dash,lw=1.5)
# hline!([1.25],color=:deepskyblue,linestyle=:dash,lw=1.5)
# hline!([1.75],color=:deepskyblue,linestyle=:dash,lw=1.5)
# vline!([0.30410821643286573],color=:deepskyblue,linestyle=:dash,lw=1.5)
# vline!([0.6958917835671342],color=:deepskyblue,linestyle=:dash,legend=false,lw=1.5)


right = plot(θ_array,M_one_spin.(θ_array,0.25),ylabel=L"M",xlabel=L"\theta/\pi",legend=false,color=:orange,lw=1.5,title="Cut B")
scatter!([0.30410821643286573,0.5,0.6958917835671342],[0.5773500286109644,0.5,0.5773500286109644],ms=3,color=:blue,titlefontsize=8,guidefont=10,size=(400,200))



left = plot(θ_array,M_one_spin.(θ_array,0.5),ylabel=L"M",xlabel=L"\theta/\pi",legend=false,color=:orange,lw=1.5,title="Cut A",titlefontsize=8,guidefont=10,size=(400,200))


θ_array = Vector( range(0,1, length=999) )
ϕ_array = Vector( range(0,2, length=1000) )
top = contourf(θ_array,ϕ_array,M_one_spin,xlab=L"\theta/\pi", ylab=L"\phi/\pi",title=L"M(\theta,\phi) = M\left(\cos(\theta/2)|0\rangle + \exp(i\phi)\sin(θ/2)|1\rangle \right)",dpi=600,guidefontsize=10,titlefontsize=10)
hline!([0.25,0.5],color=:black,lw=1.5)
hline!([0.25,0.5],color=:white,linestyle=:dash,lw=1.5)
# scatter!([0.05,0.05], [0.35,0.6], series_annotations = ["b",Plots.text("a)", :white,fontsize=6)],alpha=0,legend=false,titlefontsize=10,guidefont=10,color=:white)
annotate!([(0.05,0.5, ("A", 8, :bottom, :white)),(0.05,0.03, ("B", 8, :bottom, :white))],legend=false)

l = @layout [top; left right] 



plot(top, left, right, layout = l,size=(400,400),dpi=650)
savefig("figs/M_mapa_mundi.pdf")
savefig("figs/M_mapa_mundi.png")

