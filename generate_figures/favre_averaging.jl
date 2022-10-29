include("utils.jl")

data_path = "/Users/andresouza/Desktop/Data/FluxDifferencingPaper/"
data = "SmallHeldSuarezStatisticsFavre_Nev8_Neh10_Nq1_5_Nq2_5_Nq3_5_X_20.0.jld2"

filename = data_path * data

fig = Figure(resolution=(2496, 1387))
add_label = true
println("looking at ", filename)
state_names = []

i = 1
ii = (i - 1) % 4 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 4 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")
s_string = "u"
slice_zonal = grab_state(s_string, jl_file)
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal)
colorrange = (-40, 40) # override
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax1 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax1, ϕ, p_coord, slice_zonal,
    contour_levels, colorrange,
    add_labels=add_label, random_seed=1)
oldslice_zonal = copy(slice_zonal)

i = 2
ii = (i - 1) % 4 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 4 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")
s_string = "ρu"
slice_zonal = grab_state(s_string, jl_file)
s_string2 = "ρ"
slice_zonal2 = grab_state(s_string2, jl_file)
_, _, _ = plot_helper(s_string, slice_zonal)
colorrange = (-40, 40) # override
λ, ϕ, r, p_coord = grab_grid(jl_file)
s_string = "⟨ ρu ⟩ / ⟨ ρ ⟩"
push!(state_names, s_string)
ax2 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax2, ϕ, p_coord, slice_zonal ./ slice_zonal2,
    contour_levels, colorrange,
    add_labels=add_label, random_seed=1)
hideydecorations!(ax2)

i = 3
ii = (i - 1) % 4 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 4 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")
s_string = "ρu"
slice_zonal = grab_state(s_string, jl_file)
s_string2 = "ρ"
slice_zonal2 = grab_state(s_string2, jl_file)
_, _, _ = plot_helper(s_string, slice_zonal)
colorrange = (-40, 40) # override
λ, ϕ, r, p_coord = grab_grid(jl_file)
s_string = "Difference x 200"
push!(state_names, s_string)
ax3 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax3, ϕ, p_coord, 200 * (oldslice_zonal - slice_zonal ./ slice_zonal2),
    contour_levels, colorrange,
    add_labels=add_label, random_seed=1)
hideydecorations!(ax3)

i = 4
ii = (i - 1) % 4 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 4 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")
s_string = "ρu"
slice_zonal = grab_state(s_string, jl_file)
s_string2 = "ρ"
slice_zonal2 = grab_state(s_string2, jl_file)
_, _, _ = plot_helper(s_string, slice_zonal)
colorrange = (-40, 40) # override
λ, ϕ, r, p_coord = grab_grid(jl_file)
s_string = "⟨ρu⟩"
push!(state_names, s_string)
ax4 = fig[jj, ii] = Axis(fig, title=state_names[4], titlesize=40)
contour_heatmap!(ax4, ϕ, p_coord, slice_zonal,
    contour_levels, colorrange,
    add_labels=add_label, random_seed=1)
hideydecorations!(ax4)

i = 5
ii = (i - 1) % 4 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 4 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")
s_string = "T"
slice_zonal = grab_state(s_string, jl_file)
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal)
colorrange = (180, 310)
λ, ϕ, r, p_coord = grab_grid(jl_file)
push!(state_names, s_string)
ax5 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax5, ϕ, p_coord, slice_zonal,
    contour_levels, colorrange,
    add_labels=add_label, random_seed=1, colormap=:thermometer)
oldslice_zonal = copy(slice_zonal)


i = 6
ii = (i - 1) % 4 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 4 + 1 # +1 on why 1 based indexing is wrong
s_string = "ρT"
slice_zonal = grab_state(s_string, jl_file)
s_string = "⟨ ρT ⟩ / ⟨ ρ ⟩"
push!(state_names, s_string)
ax6 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax6, ϕ, p_coord, slice_zonal ./ slice_zonal2,
    contour_levels, colorrange,
    add_labels=add_label, random_seed=1, colormap=:thermometer)

hideydecorations!(ax6)

i = 7
ii = (i - 1) % 4 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 4 + 1 # +1 on why 1 based indexing is wrong
s_string = "Difference x 200"
push!(state_names, s_string)
colorrange_Δ = (0, 50)
ax7 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_levels_Δ = [0, 10, 20, 30, 40, 50]
contour_heatmap!(ax7, ϕ, p_coord, 200 * (oldslice_zonal - slice_zonal ./ slice_zonal2),
    contour_levels_Δ, colorrange_Δ,
    add_labels=add_label, random_seed=1, colormap=:bone_1)
hideydecorations!(ax7)

i = 8
ii = (i - 1) % 4 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 4 + 1 # +1 on why 1 based indexing is wrong
s_string = "⟨ρT⟩"
push!(state_names, s_string)
ax8 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax8, ϕ, p_coord, slice_zonal,
    contour_levels, colorrange,
    add_labels=add_label, random_seed=1, colormap=:thermometer)
hideydecorations!(ax8)


# Meriodonal Heat Transport
ρvT = grab_second_moment("ρvT", jl_file)
ρv = grab_state("ρv", jl_file)
ρT = grab_state("ρT", jl_file)
v = grab_state("v", jl_file)
T = grab_state("T", jl_file)
vT = grab_second_moment("vT", jl_file)
ρ = grab_state("ρ", jl_file)
ρvT ./ ρv

s_string = "vT"
slice_zonal, s_string = eddy_variance(s_string, jl_file)
colorrange, contour_levels, s_string = plot_helper(s_string, slice_zonal)
colorrange = (-24, 24) # override


i = 9
ii = (i - 1) % 4 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 4 + 1 # +1 on why 1 based indexing is wrong
jl_file = jldopen(filename, "r+")
s_string = "⟨v'T'⟩"

vpTp = vT - v .* T
push!(state_names, s_string)
ax9 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax9, ϕ, p_coord, vpTp,
    contour_levels, colorrange,
    add_labels=add_label, random_seed=1, colormap=:balance)
display(fig)

i = 10
ii = (i - 1) % 4 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 4 + 1 # +1 on why 1 based indexing is wrong
s_string = "⟨ ρvT ⟩ / ⟨ ρ ⟩ - ⟨ ρv ⟩ ⟨ ρT ⟩ / ⟨ ρ ⟩²"
vpTpF = ρvT ./ ρ - ρv .* ρT ./ ( ρ.^2)
push!(state_names, s_string)
ax10 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax10, ϕ, p_coord, vpTpF,
    contour_levels, colorrange,
    add_labels=add_label, random_seed=1, colormap=:balance)

hideydecorations!(ax10)
display(fig)

i = 11
ii = (i - 1) % 4 + 1 # +1 on why 1 based indexing is wrong
jj = (i - 1) ÷ 4 + 1 # +1 on why 1 based indexing is wrong
s_string = "Difference"
push!(state_names, s_string)
ax11 = fig[jj, ii] = Axis(fig, title=state_names[i], titlesize=40)
contour_heatmap!(ax11, ϕ, p_coord, vpTp - vpTpF,
    contour_levels, colorrange,
    add_labels=add_label, random_seed=1, colormap=:balance)
hideydecorations!(ax11)

close(jl_file)