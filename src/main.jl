using LinearAlgebra
using Einsum

σ_0 = I(2)
σ_1 = [[0, 1] [1, 0]]
σ_2 = [[0, -im] [im, 0]]
σ_3 = [[1, 0] [0, -1]]

σ_vec = [σ_0,σ_1,σ_2,σ_3]

function spin_combinations(i,N_spins)
    return digits(i, base=4, pad=N_spins)
end

function O_basis_spin_half(N_spins)
    O_temp = ComplexF64.( ones(2^N_spins, 2^N_spins, 4^N_spins) )
    for n=1:(4^N_spins)
        indx_array =  spin_combinations(n-1,N_spins)
        element = σ_vec[indx_array[1]+1]
        for m=2:N_spins
            element = kron(element,σ_vec[indx_array[m]+1])
        end
        O_temp[:,:,n] = element./(2^(N_spins/2))
    end
    return O_temp
end


function ρ_vec(ρ, O_basis)
    N = size(O_basis)[3]
    ρ_temp = ComplexF64.(zeros(N))
    for n=1:N
        ρ_temp[n] = tr(ρ*O_basis[:,:,n]')
    end
    return ρ_temp
end

function C_λμ(O_λ,O_μ)
    return O_λ*O_μ - O_μ*O_λ
end

function t_tensor(O_basis)
    N = size(O_basis)[3]
    t_temp = zeros(N,N)
    for λ=1:N
        for μ=1:λ
            C_temp = C_λμ(O_basis[:,:,λ],O_basis[:,:,μ])
            val = real.(  tr(C_temp*C_temp') )
            t_temp[λ,μ] = val
            t_temp[μ,λ] = val
        end
    end
    return t_temp
end

function M(ρ,O_basis)
    ρ_2_temp =  norm.(ρ_vec(ρ, O_basis )).^2
    return √( ρ_2_temp'*t_tensor(O_basis)*ρ_2_temp )
end