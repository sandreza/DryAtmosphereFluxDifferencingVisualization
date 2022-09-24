using HDF5, GLMakie, Statistics, Random

data_path = "/Users/andresouza/Desktop/Data/FluxDifferencingPaper/"
data = "convection_visualization_more_rez.h5"
h5file = h5open(data_path * data, "r+")

new_θ = read(h5file["th3d"])
thm = read(h5file["thm"])
thpthp = read(h5file["thpthp"])
wpthp = read(h5file["wpthp"])
zlist = read(h5file["z"])
close(h5file)
M = size(thm, 1)
##
begin
    options = (; titlesize=30, ylabelsize=32,
        xlabelsize=32, xgridstyle=:dash, ygridstyle=:dash, xtickalign=1,
        xticksize=10, ytickalign=1, yticksize=10,
        xticklabelsize=30, yticklabelsize=30, ylabel="z [km]")

    fig = Figure(resolution=(2500, 1400))
    volume_width = 3
    height_index = floor(Int, M / 10 * 4) # 120 for M = 300
    updraftquantile = 0.83 # 0.83

    ax = LScene(fig[1:8, 2:volume_width+1], scenekw=(camera=cam3d!, show_axis=true))
    stat_ax1 = Axis(fig[2:3, volume_width+2]; title="⟨θ⟩ [K]", options...)
    stat_ax2 = Axis(fig[4:5, volume_width+2]; title="⟨θ'w'⟩ [K m s⁻¹]", options...)
    stat_ax3 = Axis(fig[6:7, volume_width+2]; title="⟨θ'θ'⟩ [K²]", options...)

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

    tstyle[:textsize] = 05
    tstyle[:textcolor] = (:black, :black, :black)
    tstyle[:font] = "helvetica"
    tstyle[:gap] = 10
    axis[:ticks][:textcolor] = :black
    axis[:ticks][:textsize] = 05
    cbar1 = Colorbar(fig[2:7, 1], v1, label="θ [K]", width=25, ticklabelsize=30,
        labelsize=30, ticksize=35, tickalign=1, height=Relative(3 / 4)
    )

    axis[:ticks][:ranges] = ([0.0, 0.5, 1.0], [0.0, 0.5, 1.0], [0.0, 0.5, 1.0])
    axis[:ticks][:labels] = (["-1.5", "0", "1.5"], ["-1.5", "0", "1.5"], ["0", "0.7", "1.2"])

    line_options = (; linewidth=3, color=cmapa[end-4])
    lines!(stat_ax1, thm, zlist ./ 1e3; line_options...)
    xlims!(stat_ax1, (302.7, 306.3))
    ylims!(stat_ax1, (0.0, 2.0))
    lines!(stat_ax2, wpthp, zlist ./ 1e3; line_options...)
    ylims!(stat_ax2, (0.0, 2.0))
    xlims!(stat_ax2, (-0.03, 0.11))
    lines!(stat_ax3, thpthp, zlist ./ 1e3; line_options...)
    ylims!(stat_ax3, (0.0, 2.0))
    rotate_cam!(ax.scene, (π/16, -π/16, 0))
    display(fig)
end
