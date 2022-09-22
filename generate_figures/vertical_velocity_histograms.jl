using HDF5, GLMakie

data_path = "/Users/andresouza/Desktop/Data/FluxDifferencingPaper/"
data_path = "/Users/andresouza/Desktop/Julia/DryAtmosphereFluxDifferencingVisualization/"
data1 = "zonalbandw_regular_planet.h5"
data2 = "zonalbandw_small_planet.h5"

# Load data
file1 = h5open(data_path * data1, "r")
file2 = h5open(data_path * data2, "r")

wzonalbands_regular = copy(read(file1["wzonalband"]))
wzonalbands_small = copy(read(file2["wzonalband"]))
close(file1)
close(file2)

##
histfig = Figure()
ax11 = Axis(histfig[1, 1]; title="latitude = 90")
ax12 = Axis(histfig[1, 2]; title="latitude = 45")
ax21 = Axis(histfig[2, 1]; title="latitude = 60")
ax22 = Axis(histfig[2, 2]; title="latitude = 15")

lat_index = 90
hist!(ax11, wzonalbands_small[:, lat_index, :][:], bins=1000, color=:red, alpha=0.5)
hist!(ax11, wzonalbands_regular[:, lat_index, :][:], bins=1000, color=:blue, alpha=0.5)
xlims!(ax11, (-0.5, 0.5))

lat_index = 45
hist!(ax12, wzonalbands_small[:, lat_index, :][:], bins=1000, color=:red, alpha=0.5)
hist!(ax12, wzonalbands_regular[:, lat_index, :][:], bins=1000, color=:blue, alpha=0.5)
xlims!(ax12, (-0.5, 0.5))

lat_index = 60
hist!(ax21, wzonalbands_small[:, lat_index, :][:], bins=1000, color=:red, alpha=0.5)
hist!(ax21, wzonalbands_regular[:, lat_index, :][:], bins=1000, color=:blue, alpha=0.5)
xlims!(ax21, (-0.5, 0.5))


lat_index = 15
hist!(ax22, wzonalbands_small[:, lat_index, :][:], bins=1000, color=:red, alpha=0.5)
hist!(ax22, wzonalbands_regular[:, lat_index, :][:], bins=1000, color=:blue, alpha=0.5)
xlims!(ax22, (-0.5, 0.5))

display(histfig)
##
histfig = Figure()
ax11 = Axis(histfig[1, 1]; title="HS")
ax12 = Axis(histfig[1, 2]; title="Small Planet")

lat_index = 90
hist!(ax11, wzonalbands_small[:, lat_index, :][:], bins=1000, color=:red, alpha=0.5)
hist!(ax11, 15 .* wzonalbands_regular[:, lat_index, :][:], bins=1000, color=:blue, alpha=0.5)
xlims!(ax11, (-0.5, 0.5))
xlims!(ax12, (-0.5, 0.5))

display(histfig)

##

for i in 1:5:91
    μ = mean(wzonalbands_small[:, i, :][:])
    σ = std(wzonalbands_small[:, i, :][:])
    median_μ = median(wzonalbands_small[:, i, :][:])
    skewness1 = mean((wzonalbands_small[:, i, :][:] .- μ) .^ 3) / σ^3
    skewness2 = 3 * (μ - median_μ) / σ
    println("-----")
    println("for i = $i")
    println("mean: $μ, std: $σ, median: $median_μ, skewness1: $skewness1, skewness2: $skewness2")
    println("------")
end