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

# Convection Mixing Layer Depth
include("mixing_layer_depth.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
figure_number += 1

# Three Dimensional Convection Visualization
include("convection_visualization.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
figure_number += 1

# Held-Suarez Statistics
include("held_suarez_statistics.jl")
save(figure_directory * "Figure4.png", fig)

# Small-Planet Statistics 
include("held_suarez_statistics.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
figure_number += 1

# Held-Suarez Small vs Regular Planet
include("held_suarez_small_vs_regular.jl")
save(figure_directory * "Figure" * string(figure_number) * ".png", fig)
figure_number += 1