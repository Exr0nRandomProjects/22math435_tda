using LinearAlgebra
using Ripserer
using Plots
using NPZ
using Chain

# println(ARGS[1])

torus = @chain npzread("torus.npy") reshape(_, :, 3)
println("torus size", torus |> size)

# function annulus(n, r1=1, r2=2, offset=(0, 0))
#     result = Tuple{Float64,Float64}[]
#     while length(result) < n
#         point = 2 * r2 * rand(2) .- r2
#         if r1 < norm(point) < r2
#             push!(result, (point[1] + offset[1], point[2] + offset[2]))
#         end
#     end
#     return result
# end
# 
# data = annulus(300)
# println("annulus", data |> typeof, data |> size)

data = torus |> eachrow |> collect 

# data = [ i for i in size(torus)[1] ]
# # println(data |> size)
# println(torus |> eachrow |> collect |> typeof)
# # println(torus |> typeof, torus |> size, torus |> eachrow |> vec |> typeof, torus |> vec |> size)

println("ripping...")
diagram_cocycles = ripserer(data; reps=true)
println("ripped.")

exit()

# plot(torus[:, 1], torus[:, 2], torus[:, 3], seriestype=:scatter, show=true)
# println("torus plotted")
# readline()

diagram_cocycles = ripserer(data; reps=true)
most_persistent_co = diagram_cocycles[2][end]
filtration = diagram_cocycles[2].filtration

midpoint = (death(most_persistent_co) - birth(most_persistent_co)) / 2
reconstructed_at_midpoint = reconstruct_cycle(filtration, most_persistent_co, midpoint)

scatter(data; label="data", markersize=2, aspect_ratio=1)
print(reconstructed_at_midpoint |> shape)
plot!(reconstructed_at_midpoint, data; label="reconstruction", show=true)
