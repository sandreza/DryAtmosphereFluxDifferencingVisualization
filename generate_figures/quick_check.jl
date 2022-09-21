

data_path = "/Users/andresouza/Desktop/Data/FluxDifferencingPaper/"
data = "HeldSuarezStatistics_Nev6_Neh12_Nq1_5_Nq2_5_Nq3_5.jld2"
filename = data_path * data
jl_file1 = jldopen(filename)


data = "SmallHeldSuarezStatistics_Nev6_Neh12_Nq1_5_Nq2_5_Nq3_5_X_20.0.jld2"
data_path = "/Users/andresouza/Desktop/Data/FluxDifferencingPaper/"
filename = data_path * data
jl_file2 = jldopen(filename)

##
v1 = grab_state("v", jl_file1)
v2 = grab_state("v", jl_file2)
w1 = grab_state("w", jl_file1)
w2 = grab_state("w", jl_file2)
ww1 = grab_second_moment("ww", jl_file1)
ww2 = grab_second_moment("ww", jl_file2)
p1 = mean(grab_state("p", jl_file1), dims=1)[:]
p2 = mean(grab_state("p", jl_file2), dims=1)[:]

ϕ = collect(1:180);

##
fig = Figure()
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])

heatmap!(ax1, ϕ, -p1, v1, colormap=:balance, colorrange=(-2, 2), interpolate=true)
heatmap!(ax2, ϕ, -p2, v2, colormap=:balance, colorrange=(-2, 2), interpolate=true)

display(fig)

##
fig = Figure()
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
heatmap!(ax1, ww1 .- (w1 .^ 2))
heatmap!(ax2, ww2 .- (w2 .^ 2))
# heatmap!(ax1, w1) # ww2 .- (w2 .^ 2)
# heatmap!(ax2, w2)
display(fig)

##
close(jl_file1)
close(jl_file2)

##
data_path = "/Users/andresouza/Desktop/Repositories/HeldSuarezVisualizationScripts/data/"
data = "averages_small_earth_long_hs_he_12_hp_4_ve_12_vp_4_roefanov.jld2"
jl_file = jldopen(data_path * data)

##