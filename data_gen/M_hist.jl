include("../../src/src.jl") # import src.jl which has creation/annahilation operators defined

n_spins = 4

O_basis = O_basis_spin_half(n_spins)
hist_vec  = zeros(500)
for n=ProgressBar(1:500)
    ψ = randn(2^n_spins) + im*randn(2^n_spins)
    ψ = ψ/norm(ψ)
    ρ = ψ*ψ'

    ρ_vec(ρ, O_basis )

    hist_vec[n] = M(ρ, O_basis )
end
histogram(hist_vec)