include("utils.jl")

figure_directory = "figures/"

data_path = "/Users/andresouza/Desktop/Data/FluxDifferencingPaper/"
data = "SmallHeldSuarezStatisticsFavre_Nev8_Neh10_Nq1_5_Nq2_5_Nq3_5_X_20.0.jld2"

filename = data_path * data
jl_file = jldopen(filename, "r+")

# Grab states
u = grab_state("v", jl_file)
ρu = grab_state("ρu", jl_file)
v = grab_state("v", jl_file)
ρv = grab_state("ρv", jl_file)
T = grab_state("T", jl_file)
ρT = grab_state("ρT", jl_file)
ρ = grab_state("ρ", jl_file)

ρuv = grab_second_moment("ρuv", jl_file)
ρvT = grab_second_moment("ρvT", jl_file)
ρTT = grab_second_moment("ρTT", jl_file)
ρuu = grab_second_moment("ρuu", jl_file)
ρvv = grab_second_moment("ρvv", jl_file)

slice_zonal1 = ρu ./ ρ
slice_zonal2 = ρT ./ ρ
slice_zonal3 = ρTT ./ ρ - ρT .* ρT ./ (ρ .^ 2)
slice_zonal4 = ρuv ./ ρ - ρu .* ρv ./ (ρ .^ 2)
slice_zonal5 = ρvT ./ ρ - ρv .* ρT ./ (ρ .^ 2)
slice_zonal6 = 0.5 * (ρuu ./ ρ - ρu .* ρu ./ (ρ .^ 2) + ρvv ./ ρ - ρv .* ρv ./ (ρ .^ 2))

##
fig = Figure(resolution=(1700 + 600, 1000 + 400))
add_label = true
println("looking at ", filename)
state_names = []

i = 1
ii = (i - 1) % 3 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 3 + 1 # +1 on why 1 based indexing is wrong
s_string = "u"
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal1)
colorrange = (-40, 40) # override
s_string = "⟨ρu⟩ / ⟨ρ⟩"
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax1 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax1, ϕ, p_coord, slice_zonal1,
    contour_levels, colorrange,
    add_labels=add_label, random_seed=1)

i = 2
ii = (i - 1) % 3 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 3 + 1 # +1 on why 1 based indexing is wrong
s_string = "T"
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal2)
colorrange = (180, 310)
λ, ϕ, r, p_coord = grab_grid(jl_file)
s_string = "⟨ρT⟩ / ⟨ρ⟩"
push!(state_names, s_string)
ax2 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax2, ϕ, p_coord, slice_zonal2, contour_levels, colorrange, add_labels=add_label, colormap=:thermometer)
hideydecorations!(ax2, grid=false)

i = 3
ii = (i - 1) % 3 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 3 + 1 # +1 on why 1 based indexing is wrong

s_string = "⟨T'T'⟩"
s_string = "⟨ρTT⟩/⟨ρ⟩ - ⟨ρT⟩⟨ρT⟩/⟨ρ⟩²"
colorrange = (-8, 70) # override
contour_levels = [0, 20, 30, 40, 50, 60]
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax3 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax3, ϕ, p_coord, slice_zonal3, contour_levels,
    colorrange, add_labels=add_label,
    colormap=:thermometer, random_seed=18)
hideydecorations!(ax3, grid=false)

i = 4
ii = (i - 1) % 3 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 3 + 1 # +1 on why 1 based indexing is wrong
s_string = "uv"
_, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal4)
colorrange = (-60, 60)
s_string = "⟨ρuv⟩/⟨ρ⟩ - ⟨ρu⟩⟨ρv⟩/⟨ρ⟩²"
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax4 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax4, ϕ, p_coord, slice_zonal4, contour_levels, colorrange, add_labels=add_label)

i = 5
ii = (i - 1) % 3 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 3 + 1 # +1 on why 1 based indexing is wrong
s_string = "vT"
_, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal5)
# contour_levels = [-21, -18, -15, -12, -9, -6, -3, 3, 6, 9, 12, 15, 18, 21]
colorrange = (-24, 24) # override
s_string = "⟨ρvT⟩/⟨ρ⟩ - ⟨ρv⟩⟨ρT⟩/⟨ρ⟩²"
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax5 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax5, ϕ, p_coord, slice_zonal5, contour_levels, colorrange, add_labels=add_label)
hideydecorations!(ax5, grid=false)

i = 6
ii = (i - 1) % 3 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 3 + 1 # +1 on why 1 based indexing is wrong

colorrange = (0, 360) # modify to extrema 

s_string = L"\langle (u' u' + v' v')/2 \rangle"
s_string = "⟨u'u' + v'v'⟩/2"
push!(state_names, s_string)

contour_levels = [0, 40, 80, 160, 240, 320]
λ, ϕ, r, p_coord = grab_grid(jl_file)
ax6 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax6, ϕ, p_coord, slice_zonal6,
    contour_levels, colorrange, add_labels=add_label,
    colormap=:thermometer, random_seed=12)
hideydecorations!(ax6, grid=false)
display(fig)
close(jl_file)