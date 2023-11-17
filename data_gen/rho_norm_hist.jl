include("../src/src.jl") # import src.jl which has creation/annahilation operators defined
using Plots
using ProgressBars
using DataFrames
using CSV

nsmax = 5
samples = 1000

O_basis = O_basis_spin_half(n_spins)
hist_ρ_norm  = zeros(samples,nsmax)

for n_spins = 1:nsmax
    O_basis = O_basis_spin_half(n_spins)
    for n=ProgressBar(1:samples)
        ψ = randn(2^n_spins) + im*randn(2^n_spins)
        ψ = ψ/norm(ψ)
        ρ = ψ*ψ'

        hist_ρ_norm[n,n_spins] = norm(ρ_vec(ρ, O_basis ).^2)
    end
end

histogram(hist_ρ_norm,show =true,
ylab = "Counts",linecolor = nothing
,label = ["1Spin" "2Spins" "3Spins" "4Spins" "5Spins"],
xlab = L"||\rho_{\lambda}\odot \rho_{\lambda}||",
xscale = :lin, yscale = :lin)

# hist_df2 = DataFrame(hist_m,["1spin","2spin","3spin","4spin","5spin"])
# # hist_df = DataFrame(hist_m,["1spin","2spin","3spin"])
# CSV.write("data/M_spin_hist2.csv",hist_df2)



