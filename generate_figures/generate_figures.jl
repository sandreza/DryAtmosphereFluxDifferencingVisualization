figure_directory = "figures/"

### Comment ###
#= 
fig is defined within each of the included files 
=#

figure_number = 1
# Discontinuous Galerkin Intro visualizations 
include("dg_projection.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
println("done with ", figure_number)
figure_number += 1

# Convection Mixing Layer Depth
include("convection_statistics.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
println("done with ", figure_number)
figure_number += 1

# Three Dimensional Convection Visualization
include("convection_visualization.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
println("done with ", figure_number)
figure_number += 1

# Held-Suarez Statistics
include("held_suarez_statistics.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
println("done with ", figure_number)
figure_number += 1

# Held-Suarez Statistics: Favre Flavor
include("held_suarez_statistics_favre.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
println("done with ", figure_number)
figure_number += 1

# Small-Planet Statistics 
include("held_suarez_statistics_small.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
println("done with ", figure_number)
figure_number += 1

# Held-Suarez Small vs Regular Planet
include("held_suarez_compare.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
println("done with ", figure_number)
figure_number += 1
