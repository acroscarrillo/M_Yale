
include("../src/src.jl") # import src.jl which has creation/annahilation operators defined

N_spin = 5

O_basis = O_basis_spin_half(N_spin)
lamb, v = eigen(t_tensor(O_basis))
scatter(lamb)

for j=1:4^N_spin
    ρ = zeros(2^N_spin,2^N_spin)
    for i=1:4^N_spin
        ρ += ComplexF64.(v[i,j])*O_basis[:,:,i]
    end 
    if abs(tr(ρ))>1e-15
        display(j)
        display(tr(ρ))
    end
end

display(ρ)
display(tr(ρ))
ρ = ρ/tr(ρ)
display(tr(ρ))
display(tr(ρ^2))