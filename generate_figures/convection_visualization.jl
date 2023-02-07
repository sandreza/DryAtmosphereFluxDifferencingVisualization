using HDF5, GLMakie, Statistics, Random

data = "convection_visualization_more_rez.h5"
h5file = h5open(data_path * data, "r+")

new_θ = read(h5file["th3d"])
thm = read(h5file["thm"])
thpthp = read(h5file["thpthp"])
wpthp = read(h5file["wpthp"])
zlist = read(h5file["z"])
close(h5file)
M = size(thm, 1)



options = (; titlesize=40, ylabelsize=40,
    xlabelsize=40, xgridstyle=:dash, ygridstyle=:dash, xtickalign=1,
    xticksize=20, ytickalign=1, yticksize=20,
    xticklabelsize=40, yticklabelsize=40, ylabel="z [km]")

fig = Figure(resolution=(2500, 1400))
volume_width = 3
height_index = floor(Int, M / 10 * 4) # 120 for M = 300
updraftquantile = 0.83 # 0.83


cbar_fig = fig[1, 1] = GridLayout()
volume_fig = fig[1, 2] = GridLayout()
stats_fig = fig[1, 3] = GridLayout()
ax = LScene(volume_fig[1, 1], scenekw=(camera=cam3d!, show_axis=true))
stat_ax1 = Axis(stats_fig[1, 1]; title="⟨θ⟩ [K]", options...)
stat_ax2 = Axis(stats_fig[2, 1]; title="⟨θ'w'⟩ [K m s⁻¹]", options...)
stat_ax3 = Axis(stats_fig[3, 1]; title="⟨θ'θ'⟩ [K²]", options...)

cmap = :linear_kryw_0_100_c71_n256

cmapa = reverse(RGBAf.(to_colormap(cmap)))
cmap = vcat(fill(RGBAf(0, 0, 0, 0), floor(Int, 256 * (1 / (1 - updraftquantile) - 1))), cmapa[1:256])

plot_state = new_θ[:, :, 1:height_index]
println(zlist[height_index])

clims = (quantile(plot_state[:], 0.00), mean(plot_state[:, :, 1]))

v1 = volume!(ax, 0 .. 1, 0 .. 1, 0 .. 1, plot_state,
    colorrange=clims, algorithm=:absorption, absorption=10.0f0,
    colormap=cmap)

axis = ax.scene[OldAxis]
axis[:names, :axisnames] = ("x [km]", "y [km]", "z [km]")
tstyle = axis[:names] #  get the nested attributes and work directly with them

tstyle[:fontsize] = 05
tstyle[:textcolor] = (:black, :black, :black)
tstyle[:font] = "helvetica"
tstyle[:gap] = 10
axis[:ticks][:textcolor] = :black
axis[:ticks][:fontsize] = 05

cbar1 = Colorbar(cbar_fig[1, 1], v1, label="θ [K]", width=25, ticklabelsize=40,
    labelsize=50, ticksize=40, tickalign=1, height=Relative(3 / 4)
)

axis[:ticks][:ranges] = ([0.0, 0.5, 1.0], [0.0, 0.5, 1.0], [0.0, 0.5, 1.0])
axis[:ticks][:labels] = (["-1.5", "0", "1.5"], ["-1.5", "0", "1.5"], ["0", "0.7", "1.2"])

line_options = (; linewidth=3, color=cmapa[end-4])
lines!(stat_ax1, thm, zlist ./ 1e3; line_options...)
xlims!(stat_ax1, (302.7, 306.3))
ylims!(stat_ax1, (0.0, 2.0))
lines!(stat_ax2, wpthp, zlist ./ 1e3; line_options...)
ylims!(stat_ax2, (0.0, 2.0))
xlims!(stat_ax3, (-0.01, 0.10))
lines!(stat_ax3, thpthp, zlist ./ 1e3; line_options...)
# xlims!(stat_ax2, (-0.03, 0.11))
ylims!(stat_ax3, (0.0, 2.0))

rowgap!(stats_fig, 40.0)
# rowsize!(volume_fig, 1, Auto(10.5))
colsize!(fig.layout, 2, Auto(3.0))

# Modify the 3D plot
translate_cam!(ax.scene, (-0.10, -0.05, 0.0))
zoom!(ax.scene, 0.85)
rotate_cam!(ax.scene, (π / 16, 0 * -π / 16, 0))
display(fig)