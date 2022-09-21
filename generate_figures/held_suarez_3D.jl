include("utils.jl")
# filename = "data/small_earth_snapshots.jld2"
#filename = "data/earth_snapshots.jld2"
# jl_file = jldopen(filename, "r+")

# anomaly
# state_val = jl_file["T"]["1"][:,:,:,1] # ./ jl_file["ρ"]["1"][:,:,:,1]
# state_val = state_val .- mean(state_val, dims = 1)
# velocity

#=
data_path = "/Users/andresouza/Desktop/Repositories/HeldSuarezVisualizationScripts/data/"
data = "small_earth_snapshots.jld2"
# data = "earth_snapshots.jld2"
jl_file = jldopen(data_path * data, "r+")
state_val = jl_file["ρw"]["1"] ./ jl_file["ρ"]["1"]
close(jl_file)
=#


data_path = "/Users/andresouza/Desktop/Data/FluxDifferencingPaper/"
data = "HeldSuarezStatistics_Nev6_Neh12_Nq1_5_Nq2_5_Nq3_5.jld2"
data = "SmallHeldSuarezStatistics_Nev6_Neh12_Nq1_5_Nq2_5_Nq3_5_X_20.0.jld2"
filename = data_path * data
jl_file = jldopen(data_path * data, "r+")
state_val = jl_file["instantaneous"]["w"]
close(jl_file)


title_string = "Vertical Velocity"
height_index_start = 1
height_index_end = size(state_val, 3)

clims = quantile(state_val[:, :, height_index_start:height_index_end][:], 0.99)
clims = (-clims, clims)

fig = Figure(resolution=(1520, 980))
ax = LScene(fig, scenekw=(camera=cam3d!, show_axis=true))
ax_text = Label(fig, title_string,
    textsize=30, color=(:black, 0.85))

cmap = :balance # :Blues_9
cmapa = RGBAf.(to_colormap(cmap))
cmap = vcat(cmapa[1:50], fill(RGBAf(0, 0, 0, 0), 50), cmapa[end-50+1:end])

v1 = volume!(ax, 0 .. 20, 0 .. 10, 0 .. 5, state_val[:, 1:180, height_index_start:height_index_end],
    colorrange=clims, algorithm=:absorption, absorption=10.0f0,
    colormap=cmap)
axis = ax.scene[OldAxis]
axis[:names, :axisnames] = ("longitude [ᵒ]", "latitude [ᵒ]", "height [km]")
tstyle = axis[:names] #  get the nested attributes and work directly with them

tstyle[:textsize] = 15
tstyle[:textcolor] = (:black, :black, :black)
tstyle[:font] = "helvetica"
tstyle[:gap] = 10
axis[:ticks][:textcolor] = :black
axis[:ticks][:textsize] = 10
cbar1 = Colorbar(fig, v1, label=L" $w$ [m/s]", width=25, ticklabelsize=30,
    labelsize=30, ticksize=25, tickalign=1, height=Relative(3 / 4)
)

axis[:ticks][:ranges] = ([0.0, 5.0, 10.0, 15.0, 20.0], [0.0, 2.5, 5.0, 7.5, 10.0], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
axis[:ticks][:labels] = (["180W", "90W", "0", "90E", "180E"], ["90S", "45S", "0", "45N", "90N"], ["0", "6", "12", "18", "24", "30"])

fig[2:10, 1:10] = ax
fig[3:8, 11] = cbar1
fig[1, 5:6] = ax_text

zoom!(ax.scene, 0.8)
# update_cam!(fig.scene)
display(fig)

close(jl_file)