using SpecialFunctions
using LinearAlgebra
# using GR
using OffsetArrays
# using Cubature
using Plots

Float = Float64
function function_f(x, α)
    f = 1
    return f
end

function u_exact(x, α)
    α = Float(α)
    CC = 2^(-α) * gamma(1 / 2) / (gamma(1 + α / 2) * gamma((1 + α) / 2))
    ut = CC * (x * (1 - x))^(α / 2)
    return ut
end


function Gridmesh(Ω, N, r)
    grid = OffsetArray(zeros(Float, 2N + 1), 0:2N)
    a = Ω[1]
    b = Ω[2]

    for j = 0:N
        grid[j] = (Float((b - a)) / 2) * (j // N)^r .+ a
    end
    for j = N+1:2N
        grid[j] = -(Float((b - a)) / 2) * (2 - j // N)^r .+ b
    end
    # end

    return grid
end





Ω = (0, 1)


α = 1+ 1 // 10
α = Float(α)

@show 1 < α < 2

r = 4/α
N = 40


# function num_err(Ω, r, α, N)

x = Gridmesh(Ω, N, r)        # [x0, x1, ..., x2N]
h = x |> parent |> diff      # [h1, ...h2N], h_j = x_j - x_{j-1}

CR = 1 / (2 * cos(π * α / 2) * gamma(2 - α))

## K_{ij} = |x_i-x_j|^{3-α}

K = zeros(Float, 2N + 1, 2N + 1)
K = OffsetArray(K, 0:2N, 0:2N)

for i = 0:2N, j = 0:2N
    K[i, j] = abs(x[i] - x[j])^(3 - α)
end
println("Kij is OK")
## ̃a_{ij} = Cₐ (|x_{j+1} - x_i|^{3-α} / h_{j+1} - (h_j+h_{j+1})/(h_j*h_{j+1}) * |x_j - x_i|^{3-α} +  |x_{j-1} - x_i|^{3-α}/h_j)

A_1 = OffsetArray(zeros(Float, 2N + 1, 2N - 1), 0:2N, 1:2N-1)

C_a = 1 / (2 - α) / (3 - α)
Threads.@threads for j = 1:2N-1
    A_1[:, j] = C_a * (K[:, j+1] / h[j+1] - (h[j] + h[j+1]) / (h[j] * h[j+1]) * K[:, j] + K[:, j-1] / h[j])
end
println("A_1 is OK")
## a_{ij} = 2 Cᵣ ( ̃a_{i+1,j}/ h_{i+1} - (h_i+h_{i+1})/(h_i*h_{i+1}) ̃a_{ij} +  a_{i-1,j}/ h_{i})

A = zeros(Float, 2N - 1, 2N - 1)

Threads.@threads for i = 1:2N-1
    A[i, :] = 2 * CR * (A_1[i+1, :] / h[i+1] / (h[i] + h[i+1]) - A_1[i, :] / (h[i] * h[i+1]) + A_1[i-1, :] / h[i] / (h[i] + h[i+1]))
end
println("A is OK")


# ∂ₜu - Δᵅ/² u  = f    f ≡ 1

T = 1
M = 2
τ = T // M

# Uⁿ⁺¹ - Uⁿ + τ A Uⁿ⁺¹ = τ F

U_0 = u_exact.(x, α)[1:2N-1]

Un = U_0

EN = zeros(Float, M, 2N-1)

for i in 1:M
    global Un, Unp1
    Unp1 = (I + τ * A ) \ (Un .+ τ*1 )  ;   Un = Unp1
    EN[i, :] = Un - U_0
end

# Unp1
# Unp1 - U_0[1:2N-1] .|> abs |> maximum

dEN = zeros(Float, M, 2N-1)
dEN[1, :] = EN[1, :]
for i in 2:M
    dEN[i, :] = (EN[i, :] - EN[i-1, :])
end

xi = x[1:2N-1]
plot()
for i in 1:M
    plot!(x[1:2N-1], dEN[i, :], legend=false)
end
plot!(xi, EN[end,:], legend=false)
ylims!(0, 0.0004)
# plot(xi, dEN[end,:]./ maximum(abs.(dEN[end,:])) .|> abs, legend=false)

# append!(dENend, dEN[end,:])
# H = (h[1:2N-1] + h[2:2N]) / 2