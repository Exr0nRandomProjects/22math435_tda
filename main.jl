using LinearAlgebra
# using Plots; pythonplot() # dependency error
using Plots; plotlyjs()
using Ripserer
using NPZ
using Chain

torus = @chain npzread("torus.npy") reshape(_, :, 3)
println(torus |> size)

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
#
#

# println("data generated")
# 
# scatter(data; label="data", markersize=2, aspect_ratio=1, show=true)
# println("plotted")

plot(torus[:, 1], torus[:, 2], torus[:, 3], seriestype=:scatter, show=true)
println("torus plotted")
readline()


# 
# scatter(data; label="data", markersize=2, aspect_ratio=1)
# plot!(most_persistent_ho, data; label="cycle", show=true)


diagram_cocycles = ripserer(data; reps=true)
most_persistent_co = diagram_cocycles[2][end]
filtration = diagram_cocycles[2].filtration

midpoint = (death(most_persistent_co) - birth(most_persistent_co)) / 2
reconstructed_at_midpoint = reconstruct_cycle(filtration, most_persistent_co, midpoint)

scatter(data; label="data", markersize=2, aspect_ratio=1)
plot!(reconstructed_at_midpoint, data; label="reconstruction", show=true)

# reconstructed_at_birth = reconstruct_cycle(filtration, most_persistent_co)
# 
# scatter(data; label="data", markersize=2, aspect_ratio=1)
# plot!(reconstructed_at_birth, data; label="reconstruction", show=true)


# readline()



# diagram = ripserer(data)
# println("ripped")
# plot(diagram, show=true)
# 
# 
# most_persistent = diagram[2][end]
# 
# 
# readline()
