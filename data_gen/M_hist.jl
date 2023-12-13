include("../src/src.jl") # import src.jl which has creation/annahilation operators defined
using Plots
using ProgressBars
using DataFrames
using CSV

nsmax = 4
samples = 1000

hist_m  = zeros(samples,nsmax)

for n_spins = 1:nsmax
    O_basis = O_basis_spin_half(n_spins)
    for n=ProgressBar(1:samples)
        ψ = randn(2^n_spins) + im*randn(2^n_spins)
        ψ = ψ/norm(ψ)
        ρ = ψ*ψ'

        hist_m[n,n_spins] = M(ρ, O_basis)
    end
end

hist_df = DataFrame(hist_m,["1spin","2spin","3spin","4spin"])
CSV.write("data/M_spin_hist.csv",hist_df)