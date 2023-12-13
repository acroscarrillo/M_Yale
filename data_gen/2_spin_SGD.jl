include("../src/src.jl") 
using Plots
using ProgressBars
using DataFrames
using CSV
using LaTeXStrings
using Flux

O_basis = O_basis_spin_half(2)
t_fixed = t_tensor(O_basis)

zerozero = [1,0,0,0]
zeroone = [0,1,0,0]
onezero = [0,0,1,0]
oneone = [0,0,0,1]

function ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_1 = cos(θ_1*π/2)*zerozero
    ψ_2 = exp(im*ϕ_1*π)*sin(θ_1*π/2)*cos(θ_2*π/2)*zeroone
    ψ_3 = exp(im*ϕ_2*π)*sin(θ_1*π/2)*sin(θ_2*π/2)*cos(θ_3*π/2)*onezero
    ψ_4 = exp(im*ϕ_3*π)*sin(θ_1*π/2)*sin(θ_2*π/2)*sin(θ_3*π/2)*oneone
    ψ = ψ_1 + ψ_2 + ψ_3 + ψ_4
    return ψ
end

function M_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_temp = ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ρ_temp = ψ_temp*ψ_temp'
    return M(ρ_temp,O_basis)
end


####################
# Derivatives of ψ #
####################
function ∂ϕ1ψ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_1 = 0*zerozero
    ψ_2 = im*exp(im*ϕ_1*π)*sin(θ_1*π/2)*cos(θ_2*π/2)*zeroone
    ψ_3 = 0*onezero
    ψ_4 = 0*oneone
    ψ = ψ_1 + ψ_2 + ψ_3 + ψ_4
    return ψ
end
function ∂ϕ2ψ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_1 = 0*zerozero
    ψ_2 = 0*zeroone
    ψ_3 = im*exp(im*ϕ_2*π)*sin(θ_1*π/2)*sin(θ_2*π/2)*cos(θ_3*π/2)*onezero
    ψ_4 = 0*oneone
    ψ = ψ_1 + ψ_2 + ψ_3 + ψ_4
    return ψ
end
function ∂ϕ3ψ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_1 = 0*zerozero
    ψ_2 = 0*zeroone
    ψ_3 = 0*onezero
    ψ_4 = im*exp(im*ϕ_3*π)*sin(θ_1*π/2)*sin(θ_2*π/2)*sin(θ_3*π/2)*oneone
    ψ = ψ_1 + ψ_2 + ψ_3 + ψ_4
    return ψ
end
function ∂θ1ψ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_1 = -sin(θ_1*π/2)*π/2*zerozero
    ψ_2 = exp(im*ϕ_1*π)*cos(θ_1*π/2)*π/2*cos(θ_2*π/2)*zeroone
    ψ_3 = exp(im*ϕ_2*π)*cos(θ_1*π/2)*π/2*sin(θ_2*π/2)*cos(θ_3*π/2)*onezero
    ψ_4 = exp(im*ϕ_3*π)*cos(θ_1*π/2)*π/2*sin(θ_2*π/2)*sin(θ_3*π/2)*oneone
    ψ = ψ_1 + ψ_2 + ψ_3 + ψ_4
    return ψ
end
function ∂θ2ψ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_1 = 0*zerozero
    ψ_2 = -exp(im*ϕ_1*π)*sin(θ_1*π/2)*sin(θ_2*π/2)*π/2*zeroone
    ψ_3 = exp(im*ϕ_2*π)*sin(θ_1*π/2)*cos(θ_2*π/2)*π/2*cos(θ_3*π/2)*onezero
    ψ_4 = exp(im*ϕ_3*π)*sin(θ_1*π/2)*cos(θ_2*π/2)*π/2*sin(θ_3*π/2)*oneone
    ψ = ψ_1 + ψ_2 + ψ_3 + ψ_4
    return ψ
end
function ∂θ3ψ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_1 = 0*zerozero
    ψ_2 = 0*zeroone
    ψ_3 = -exp(im*ϕ_2*π)*sin(θ_1*π/2)*sin(θ_2*π/2)*sin(θ_3*π/2)*π/2*onezero
    ψ_4 = exp(im*ϕ_3*π)*sin(θ_1*π/2)*sin(θ_2*π/2)*cos(θ_3*π/2)*π/2*oneone
    ψ = ψ_1 + ψ_2 + ψ_3 + ψ_4
    return ψ
end

####################
# Derivatives of ρ #
####################
function ∂ϕ1ρ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_temp = ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ∂ψ_temp = ∂ϕ1ψ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    temp = ψ_temp*∂ψ_temp'
    return temp + temp'
end
function ∂ϕ2ρ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_temp = ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ∂ψ_temp = ∂ϕ2ψ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    temp = ψ_temp*∂ψ_temp'
    return temp + temp'
end
function ∂ϕ3ρ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_temp = ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ∂ψ_temp = ∂ϕ3ψ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    temp = ψ_temp*∂ψ_temp'
    return temp + temp'
end
function ∂θ1ρ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_temp = ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ∂ψ_temp = ∂θ1ψ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    temp = ψ_temp*∂ψ_temp'
    return temp + temp'
end
function ∂θ2ρ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_temp = ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ∂ψ_temp = ∂θ2ψ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    temp = ψ_temp*∂ψ_temp'
    return temp + temp'
end
function ∂θ3ρ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ψ_temp = ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ∂ψ_temp = ∂θ3ψ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    temp = ψ_temp*∂ψ_temp'
    return temp + temp'
end

function newϕ1(η,θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ∂ρ_temp = ∂ϕ1ρ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)

    ρ_temp = ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)*ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)'
    ρ_vec_temp = ρ_vec(ρ_temp,O_basis)
    ρ_vec_2 = ρ_vec_temp.^2
    ρ_vec_∂ρ_vec = ρ_vec_temp.*ρ_vec(∂ρ_temp,O_basis)
    return real( ϕ_1 + η*(ρ_vec_2'*t_fixed*ρ_vec_∂ρ_vec) )
end
function newϕ2(η,θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ∂ρ_temp = ∂ϕ2ρ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)

    ρ_temp = ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)*ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)'
    ρ_vec_temp = ρ_vec(ρ_temp,O_basis)
    ρ_vec_2 = ρ_vec_temp.^2
    ρ_vec_∂ρ_vec = ρ_vec_temp.*ρ_vec(∂ρ_temp,O_basis)
    return real( ϕ_2 + η*(ρ_vec_2'*t_fixed*ρ_vec_∂ρ_vec) )
end
function newϕ3(η,θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ∂ρ_temp = ∂ϕ3ρ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)

    ρ_temp = ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)*ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)'
    ρ_vec_temp = ρ_vec(ρ_temp,O_basis)
    ρ_vec_2 = ρ_vec_temp.^2
    ρ_vec_∂ρ_vec = ρ_vec_temp.*ρ_vec(∂ρ_temp,O_basis)
    return real( ϕ_3 + η*(ρ_vec_2'*t_fixed*ρ_vec_∂ρ_vec) )
end
function newθ1(η,θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ∂ρ_temp = ∂θ1ρ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)

    ρ_temp = ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)*ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)'
    ρ_vec_temp = ρ_vec(ρ_temp,O_basis)
    ρ_vec_2 = ρ_vec_temp.^2
    ρ_vec_∂ρ_vec = ρ_vec_temp.*ρ_vec(∂ρ_temp,O_basis)
    return real( θ_1 + η*(ρ_vec_2'*t_fixed*ρ_vec_∂ρ_vec) )
end
function newθ2(η,θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ∂ρ_temp = ∂θ2ρ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)

    ρ_temp = ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)*ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)'
    ρ_vec_temp = ρ_vec(ρ_temp,O_basis)
    ρ_vec_2 = ρ_vec_temp.^2
    ρ_vec_∂ρ_vec = ρ_vec_temp.*ρ_vec(∂ρ_temp,O_basis)
    return real( θ_2 + η*(ρ_vec_2'*t_fixed*ρ_vec_∂ρ_vec) )
end
function newθ3(η,θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ∂ρ_temp = ∂θ3ρ(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)

    ρ_temp = ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)*ψ_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)'
    ρ_vec_temp = ρ_vec(ρ_temp,O_basis)
    ρ_vec_2 = ρ_vec_temp.^2
    ρ_vec_∂ρ_vec = ρ_vec_temp.*ρ_vec(∂ρ_temp,O_basis)
    return real( θ_3 + η*(ρ_vec_2'*t_fixed*ρ_vec_∂ρ_vec) )
end




loops_M = 100000
# M_max_array = zeros(loops_M)
# params_max_array = zeros(loops_M,6)

loops = 5000
η = .1
for m=ProgressBar(1:loops_M)
    θ_1,θ_2,θ_3 = rand(3)
    ϕ_1,ϕ_2,ϕ_3 = 2*rand(3)

    for n=1:loops
        θ_1_new = newθ1(η/√(n),θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
        θ_2_new = newθ2(η/√(n),θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
        θ_3_new = newθ3(η/√(n),θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
        ϕ_1_new = newϕ1(η/√(n),θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
        ϕ_2_new = newϕ2(η/√(n),θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
        ϕ_3_new = newϕ3(η/√(n),θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
        θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3 = θ_1_new,θ_2_new,θ_3_new,ϕ_1_new,ϕ_2_new,ϕ_3_new
    end
    M_max_array[m] = M_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    params_max_array[m,:] .= θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3
end


l = @layout [top{0.3h}; left right]

left = histogram(params_max_array[:,1],alpha=0.4,normalize=true,xlab=L"\theta_\max/\pi",label=L"\theta_1/\pi")
histogram!(params_max_array[:,2],alpha=0.4,normalize=true,label=L"\theta_2/\pi")
histogram!(params_max_array[:,3],alpha=0.4,normalize=true,label=L"\theta_3/\pi",guidefont=15,legendfontsize=12,foreground_color_legend = nothing,background_color_legend=nothing)

right = histogram(params_max_array[:,4],alpha=0.4,normalize=true,xlab=L"\phi_\max/\pi",label=L"\phi_1/\pi")
histogram!(params_max_array[:,5],alpha=0.4,normalize=true,label=L"\phi_2/\pi")
histogram!(params_max_array[:,6],alpha=0.4,normalize=true,label=L"\phi_3/\pi",guidefont=15,legendfontsize=12,foreground_color_legend = nothing,background_color_legend=nothing)

top = histogram(M_max_array,xlims=(0.50,0.534),normalize=true,legend=false,xlabel=L"M_{\max}",guidefont=15,size=(600,200))


plot(top,left,right, layout = l,size=(600,400),dpi=650)
savefig("figs/M_max_2_spins.pdf")
savefig("figs/M_max_2_spins.png")














############################
# plot some "trajectories" #
############################
loops = 5000
η = .1

θ_1,θ_2,θ_3 = rand(3)
ϕ_1,ϕ_2,ϕ_3 = 2*rand(3)
M_initial = M_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
display("Initial M=$M_initial at θ_1=$θ_1, θ_2=$θ_2, θ_3=$θ_3, ϕ_1=$ϕ_1, ϕ_2=$ϕ_2, ϕ_3=$ϕ_3")

M_array = zeros(loops)
for n=1:loops
    θ_1_new = newθ1(η/√(n),θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    θ_2_new = newθ2(η/√(n),θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    θ_3_new = newθ3(η/√(n),θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ϕ_1_new = newϕ1(η/√(n),θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ϕ_2_new = newϕ2(η/√(n),θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    ϕ_3_new = newϕ3(η/√(n),θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
    θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3 = θ_1_new,θ_2_new,θ_3_new,ϕ_1_new,ϕ_2_new,ϕ_3_new
    M_array[n] = M_2_spins(θ_1,θ_2,θ_3,ϕ_1,ϕ_2,ϕ_3)
end


M_final = M_array[end] 
display("Final M=$M_final at θ_1=$θ_1, θ_2=$θ_2, θ_3=$θ_3, ϕ_1=$ϕ_1, ϕ_2=$ϕ_2, ϕ_3=$ϕ_3")
plot!(M_array,legend=false,ylabel=L"M",xlab="learning steps",title="Gradient descend of M on the 2 spin 1/2 space",alpha=0.4)