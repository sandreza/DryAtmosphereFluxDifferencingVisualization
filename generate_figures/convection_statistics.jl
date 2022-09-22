using HDF5, GLMakie

data_path = "/Users/andresouza/Desktop/Data/FluxDifferencingPaper/"
data = "mixing_layer_depth_fixed.h5"
data = "mixing_layer_depth_more_rez.h5"
h5file = h5open(data_path * data, "r+")
entrainment_layer_depth = read(h5file["entrainment_layer_depth"])
maxind_t = read(h5file["maxind_t"])
zlist = read(h5file["zlist"])
tlist = read(h5file["tlist"])
close(h5file)



fig = Figure(resolution=(800, 600))
options = (; titlesize=30, ylabelsize=32,
    xlabelsize=32, xgridstyle=:dash, ygridstyle=:dash, xtickalign=1,
    xticksize=10, ytickalign=1, yticksize=10,
    xticklabelsize=30, yticklabelsize=30, xlabel="time [hours]", ylabel="entrainment layer height [km]")

ax1 = Axis(fig[1, 1]; options...)
ylims!(ax1, (0, 3.0))
sc = scatter!(ax1, tlist ./ (60 * 60), zlist[maxind_t] ./ 1000, color=:blue, label="Numerical")
ln = lines!(ax1, tlist ./ (60 * 60), entrainment_layer_depth ./ 1000, color=:red, label="Theoretical", linewidth=3)
axislegend(ax1, position=:rt)
display(fig)
