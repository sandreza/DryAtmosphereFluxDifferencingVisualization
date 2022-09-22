figure_directory = "figures/"

### Comment ###
#= 
fig is defined within each of the included files 
=#

figure_number = 1
# Discontinuous Galerkin Intro visualizations 
include("dg_projection.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
figure_number += 1
println("done with ", figure_number)

# Convection Mixing Layer Depth
include("convection_statistics.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
figure_number += 1
println("done with ", figure_number)

# Three Dimensional Convection Visualization
include("convection_visualization.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
figure_number += 1
println("done with ", figure_number)


# Held-Suarez Statistics
include("held_suarez_statistics.jl")
save(figure_directory * "Figure4.png", fig)
figure_number += 1
println("done with ", figure_number)

# Small-Planet Statistics 
include("held_suarez_statistics_small.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
figure_number += 1
println("done with ", figure_number)

# Held-Suarez Small vs Regular Planet
include("held_suarez_compare.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
figure_number += 1
println("done with ", figure_number)