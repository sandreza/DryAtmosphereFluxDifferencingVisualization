data_path = "/Users/andresouza/Desktop/Repositories/HeldSuarezVisualizationScripts/data/"

data_path = "/Users/andresouza/Desktop/Data/FluxDifferencingPaper/"
# data_path = "/Users/andresouza/Desktop/Julia/DryAtmosphereFluxDifferencingVisualization/"
data1 = "HeldSuarezStatisticsConsistent_Nev8_Neh10_Nq1_5_Nq2_5_Nq3_5.jld2"
data2 = "TraditionalSmallHeldSuarezStatisticsConsistent_Nev8_Neh10_Nq1_5_Nq2_5_Nq3_5_X_20.0.jld2"
data3 = "SmallHeldSuarezStatisticsConsistent_Nev8_Neh10_Nq1_5_Nq2_5_Nq3_5_X_20.0.jld2"

jl_file1 = jldopen(data_path * data1, "r+")
jl_file2 = jldopen(data_path * data2, "r+")
jl_file3 = jldopen(data_path * data3, "r+")

# fig = Figure(resolution=(3400, 1000))
fig = Figure(resolution=(2700,700))
add_label = true
state_names = []

i = 1
ii = (i - 1) % 3 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")
s_string = "u"
slice_zonal = grab_state(s_string, jl_file1)
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal)
colorrange = (-40, 40) # override
λ, ϕ, r, p_coord = grab_grid(jl_file)
s_string = "standard planet, full Coriolis"
push!(state_names, s_string)
# slice_zonal = (state_val[1:end-1, :] + state_val[2:end, :]) * 0.5 / 1000 # average
ax1 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax1, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels=add_label, random_seed=1)

i = 2
ii = (i - 1) % 3 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")
s_string = "u"
slice_zonal = grab_state(s_string, jl_file2)
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal)
colorrange = (-40, 40) # override
λ, ϕ, r, p_coord = grab_grid(jl_file)
s_string = "small planet, traditional Coriolis"
push!(state_names, s_string)
ax2 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax2, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels=add_label, random_seed=6)
hideydecorations!(ax2, grid=false)

i = 3
ii = (i - 1) % 3 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 3 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")
s_string = "u"
slice_zonal = grab_state(s_string, jl_file3)
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal)
colorrange = (-40, 40) # override
λ, ϕ, r, p_coord = grab_grid(jl_file)
s_string = "small planet, full Coriolis"
push!(state_names, s_string)
ax3 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax3, ϕ, p_coord, slice_zonal, contour_levels, colorrange, add_labels=add_label, random_seed=1)
hideydecorations!(ax3, grid=false)

display(fig)
close(jl_file1)
close(jl_file2)
close(jl_file3)