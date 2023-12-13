using LinearAlgebra
using Einsum

σ_0 = I(2)
σ_1 = [[0, 1] [1, 0]]
σ_2 = [[0, -im] [im, 0]]
σ_3 = [[1, 0] [0, -1]]

σ_vec = [σ_0,σ_1,σ_2,σ_3]

⊗(A,B) = kron(A,B)

function spin_combinations(i,N_spins)
    return digits(i, base=4, pad=N_spins)
end

function O_basis_spin_half(N_spins)
    O_temp = ComplexF64.( ones(2^N_spins, 2^N_spins, 4^N_spins) )
    for n=1:(4^N_spins)
        indx_array =  spin_combinations(n-1,N_spins)
        element = 1
        for m=1:N_spins
            element = element ⊗ σ_vec[indx_array[m]+1]
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



"""
    spin_partial_trace(ρ,trace_mask)

# Theory
In a Hilbert space `H = H_A⨂H_B`, any state `ρ` can be written as 

`ρ = ∑_jk ρ_jk A_j⨂B_k`

where `A_j ≐ |j_A⟩⟨j_A|` with `|j_A⟩` a basis element of `H_A` (and similarly for `B`). In this basis, the partial trace is straight forward:

`ρ_B ≐ Tr_A(∑_jk ρ_jk A_j⨂B_k) = ∑_jkn ρ_jk B_k ⟨n_A|A_j|n_A⟩ =  ∑_jk ρ_jk B_k.`

Note however, that from the above, it is clear that the partial trace over `A` is simply a map `Tr_A : H_A⨂H_B -> H_B` such that when acting on a basis element like above

`Tr_A(A_j⨂B_k) ≐ Tr_A(A_j) B_k  ≡ ∑_n ( ⟨n_A|⨂I_B ) A_j⨂B_k ( |n_A⟩⨂I_B ).`

So we can define the partial trace succintly as the operatorial map 

`ρ_B =  ∑_n ( ⟨n_A|⨂I_B ) ρ ( |n_A⟩⨂I_B ) `.

And this is what the code does: first it constructs `( |n_A⟩⨂I_B )` with a helper function `partial_trace_op(n,trace_mask)` and then performs the sum.


# Code
This code uses `trace_mask` to define what part of `H` is `H_A` and `H_B`. For instance, in a 1D chain of 3 spins, `trace_mask = [1,0,1]` indicates `H_A` englobe the first and last spin (the `1`s) whereas `H_B` is the second spin (the `0`s). It performs the trace over `B`.

# Examples
```julia-repl
julia> ω_a(1,[0.00075, 1.27*10^(-7)],0)
0.999996631
```
"""
function spin_partial_trace(ρ,trace_mask)
    #check trace_mask has correct length
    if (2^length(trace_mask)) != size(ρ)[1]
        n_spins = log2(size(ρ)[1])
        trace_mask_len = length(trace_mask)
        throw(ArgumentError("ρ is a state in a Hilbert space of $n_spins spins. However, a trace_mask of length $trace_mask_len was given."))
    end

    dim_H_A = 2^count( Bool.(trace_mask) )
    dim_H_B = 2^( length(trace_mask) - count( Bool.(trace_mask) ) )

    ρ_A = zeros(dim_H_A,dim_H_A)
    for n=1:dim_H_B
        op_temp = partial_trace_op(n,trace_mask)
        ρ_A += op_temp'*ρ*op_temp
    end
    return ρ_A/tr(ρ_A)
end

function partial_trace_op(n,trace_mask)
    spin_vec = [[1,0], [0,1]]
    op_temp = 1
    for (j,trace_n) in enumerate(trace_mask)
        if trace_n == 1
            op_temp = op_temp ⊗ I(2)
        else
            indx = digits(n, base=2, pad=length(trace_mask))[j]
            op_temp = op_temp ⊗ spin_vec[ indx + 1 ]
        end
    end
    return op_temp
end

function entanglement_entropy(ρ,trace_mask)
    rho_B = spin_partial_trace(ρ,trace_mask)
    lambs, _ = eigen(rho_B)
    entropy = 0
    for lamb in lambs
        if lamb == 0
            entropy += 0
        else
            entropy += -lamb * log2(lamb)
        end
    end
    return entropy
end

function entanglement_entropy(ρ)
    n_spins = Int(log2(size(ρ)[1]))
    trace_mask = zeros(n_spins)
    for n=1:Int(n_spins/2)
        trace_mask[n] = 1
    end
    return entanglement_entropy(ρ,trace_mask)
end