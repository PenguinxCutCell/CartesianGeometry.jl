using CartesianGeometry
using Test
using LinearAlgebra
using CairoMakie
using Printf
using Statistics
using BenchmarkTools

const T = Float64
const R = 0.21
const a = 0.5

# Création de maillages 3D
function create_mesh_3d(nx, ny, nz)
    dx, dy, dz = 1.0/nx, 1.0/ny, 1.0/nz
    x0, y0, z0 = 0.0, 0.0, 0.0
    return ([x0 + i*dx for i in 0:nx-1], 
            [y0 + i*dy for i in 0:ny-1], 
            [z0 + i*dz for i in 0:nz-1])
end

# Définition des fonctions levelset
levelset_sphere = (x, y, z) -> (x - a)^2 + (y - a)^2 + (z - a)^2 - R^2

function ellipsoid_levelset(x, y, z)
    a_axis = 0.3; b_axis = 0.2; c_axis = 0.15
    x -= a; y -= a; z -= a
    return ((x/a_axis)^2 + (y/b_axis)^2 + (z/c_axis)^2) - 1
end

function cube_levelset(x, y, z)
    size = 0.2
    x -= a; y -= a; z -= a
    return max(abs(x), max(abs(y), abs(z))) - size
end

function torus_levelset(x, y, z)
    x -= a; y -= a; z -= a
    r_major = 0.2  # Rayon majeur du tore
    r_minor = 0.07  # Rayon mineur du tore
    q = sqrt(x^2 + y^2) - r_major
    return sqrt(q^2 + z^2) - r_minor
end

println("===== ANALYSE DES INTERFACES ET VOLUMES 3D =====")

# 1. Étude de convergence sur la sphère
println("\n1. Convergence de l'interface et du volume pour une sphère")
resolutions = [10, 15, 20, 25, 30, 40]
interface_errors = Float64[]
volume_errors = Float64[]
analytical_surface = 4π * R^2
analytical_volume = (4/3) * π * R^3

for res in resolutions
    mesh_3d = create_mesh_3d(res, res, res)
    
    t_start = time()
    V, _, interface_area, _ = integrate(Tuple{0}, levelset_sphere, mesh_3d, T, zero)
    t_elapsed = time() - t_start
    
    numerical_volume = sum(V)
    numerical_surface = sum(interface_area)
    
    vol_error = abs(numerical_volume - analytical_volume)
    surf_error = abs(numerical_surface - analytical_surface)
    
    push!(volume_errors, vol_error)
    push!(interface_errors, surf_error)
    
    vol_rel_error = 100 * vol_error / analytical_volume
    surf_rel_error = 100 * surf_error / analytical_surface
    
    println(@sprintf("  Résolution %2d³: Temps = %.2fs, Volume = %.6f (%.2f%%), Surface = %.6f (%.2f%%)", 
                   res, t_elapsed, numerical_volume, vol_rel_error, numerical_surface, surf_rel_error))
end

fig_conv = Figure(size = (1000, 500))
ax_conv1 = Axis(fig_conv[1, 1],
               xscale = log10, yscale = log10,
               xlabel = "Résolution", ylabel = "Erreur absolue",
               title = "Convergence de la surface d'interface")
scatter!(ax_conv1, resolutions, interface_errors)
lines!(ax_conv1, resolutions, interface_errors)

# Références pour la convergence
ref_lin = interface_errors[2] * (resolutions[2] ./ resolutions)
ref_quad = interface_errors[2] * (resolutions[2] ./ resolutions).^2
lines!(ax_conv1, resolutions, ref_lin, linestyle = :dash, color = :red, label = "O(h)")
lines!(ax_conv1, resolutions, ref_quad, linestyle = :dash, color = :blue, label = "O(h²)")
axislegend(ax_conv1)

ax_conv2 = Axis(fig_conv[1, 2],
               xscale = log10, yscale = log10,
               xlabel = "Résolution", ylabel = "Erreur absolue",
               title = "Convergence du volume")
scatter!(ax_conv2, resolutions, volume_errors)
lines!(ax_conv2, resolutions, volume_errors)

# Références pour la convergence
ref_lin = volume_errors[2] * (resolutions[2] ./ resolutions)
ref_quad = volume_errors[2] * (resolutions[2] ./ resolutions).^2
lines!(ax_conv2, resolutions, ref_lin, linestyle = :dash, color = :red, label = "O(h)")
lines!(ax_conv2, resolutions, ref_quad, linestyle = :dash, color = :blue, label = "O(h²)")
axislegend(ax_conv2)

display(fig_conv)

# 2. Comparaison entre différentes géométries
println("\n2. Comparaison des interfaces et volumes entre différentes géométries 3D")
levelsets = [
    (levelset_sphere, "Sphère", 4π * R^2, (4/3) * π * R^3),
    (ellipsoid_levelset, "Ellipsoïde", nothing, nothing),
    (cube_levelset, "Cube", 6 * (2*0.2)^2, (2*0.2)^3),
    (torus_levelset, "Tore", nothing, nothing)
]

# Résolution modérée pour le benchmark des géométries
res = 30
mesh_3d = create_mesh_3d(res, res, res)

# Table des résultats
results_table = []

for (lvlset, name, analytical_surf, analytical_vol) in levelsets
    t_start = time()
    V, _, interface_area, _ = integrate(Tuple{0}, lvlset, mesh_3d, T, zero)
    t_elapsed = time() - t_start
    
    numerical_volume = sum(V)
    numerical_surface = sum(interface_area)
    
    # Statistiques
    println("  $name:")
    println("    - Volume calculé: $(numerical_volume)")
    println("    - Surface calculée: $(numerical_surface)")
    if analytical_vol !== nothing
        vol_rel_error = 100 * abs(numerical_volume - analytical_vol) / analytical_vol
        println("    - Volume théorique: $(analytical_vol)")
        println(@sprintf("    - Erreur relative volume: %.4f%%", vol_rel_error))
    end
    if analytical_surf !== nothing
        surf_rel_error = 100 * abs(numerical_surface - analytical_surf) / analytical_surf
        println("    - Surface théorique: $(analytical_surf)")
        println(@sprintf("    - Erreur relative surface: %.4f%%", surf_rel_error))
    end
    
    # Distribution des valeurs d'interface par cellule
    interface_cells = findall(interface_area .> 0)
    interface_values = interface_area[interface_cells]
    println("    - Nombre de cellules d'interface: $(length(interface_cells))")
    println(@sprintf("    - Surface min/max/moy par cellule: %.6f / %.6f / %.6f", 
                  minimum(interface_values), maximum(interface_values), mean(interface_values)))
    println(@sprintf("    - Temps de calcul: %.2f secondes", t_elapsed))
    
    push!(results_table, (name=name, volume=numerical_volume, surface=numerical_surface, 
                         cells=length(interface_cells), time=t_elapsed))
end

# 3. Analyse de la distribution des surfaces d'interface
println("\n3. Distribution des surfaces d'interface par cellule")
res = 30
mesh_3d = create_mesh_3d(res, res, res)
fig_dist = Figure(size = (1200, 800))

for (i, (lvlset, name, _, _)) in enumerate(levelsets)
    row, col = div(i-1, 2)+1, rem(i-1, 2)+1
    
    _, _, interface_area, _ = integrate(Tuple{0}, lvlset, mesh_3d, T, zero)
    
    # Extraire les valeurs non-nulles
    interface_cells = findall(interface_area .> 0)
    interface_values = interface_area[interface_cells]
    
    # Histogramme
    ax = Axis(fig_dist[row, col], 
              title = "Distribution ($name)",
              xlabel = "Surface d'interface", ylabel = "Fréquence")
    hist!(ax, interface_values, bins = 30, normalization = :probability)
    
    # Calculer la surface moyenne par cellule
    avg_area = sum(interface_area) / length(interface_cells)
    vlines!(ax, [avg_area], color = :red, linewidth = 2, 
            label = @sprintf("Moyenne: %.4f", avg_area))
    axislegend(ax)
    
    # Calcul de statistiques additionnelles
    println("  $name:")
    println(@sprintf("    - Surface moy/médiane: %.6f / %.6f", 
                   mean(interface_values), median(interface_values)))
    println(@sprintf("    - Écart-type: %.6f", std(interface_values)))
    println(@sprintf("    - Coefficient de variation: %.4f", std(interface_values)/mean(interface_values)))
end

display(fig_dist)

# 4. Influence de la résolution sur la distribution des surfaces d'interface
println("\n4. Influence de la résolution sur la distribution des surfaces d'interface")
resolutions = [15, 25, 40]
fig_res = Figure(size = (1200, 400))

for (i, res) in enumerate(resolutions)
    mesh_3d = create_mesh_3d(res, res, res)
    _, _, interface_area, _ = integrate(Tuple{0}, levelset_sphere, mesh_3d, T, zero)
    
    # Normaliser par la taille de la cellule (dx^2)
    dx = 1.0 / res
    cell_face_area = dx^2
    
    # Extraire les valeurs non-nulles
    interface_cells = findall(interface_area .> 0)
    interface_values = interface_area[interface_cells]
    normalized_values = interface_values ./ cell_face_area
    
    # Histogramme
    ax = Axis(fig_res[1, i], 
              title = @sprintf("Résolution %d×%d×%d", res, res, res),
              xlabel = "Surface d'interface / face de cellule", 
              ylabel = "Fréquence")
    hist!(ax, normalized_values, bins = 30, normalization = :probability)
    
    # Statistiques
    println(@sprintf("  Résolution %d³:", res))
    println("    - Nombre de cellules d'interface: $(length(interface_cells))")
    println(@sprintf("    - Surface normalisée min/max/moy: %.4f / %.4f / %.4f", 
                  minimum(normalized_values), maximum(normalized_values), mean(normalized_values)))
    println(@sprintf("    - Écart-type: %.4f", std(normalized_values)))
    
    # Ligne verticale pour la valeur moyenne
    vlines!(ax, [mean(normalized_values)], color = :red, linewidth = 2,
            label = @sprintf("Moyenne: %.4f", mean(normalized_values)))
    # Référence pour la diagonale d'une face (√2)
    vlines!(ax, [sqrt(2)], color = :green, linestyle = :dash, linewidth = 2,
            label = "√2")
    # Référence pour la diagonale d'un cube (√3)
    vlines!(ax, [sqrt(3)], color = :blue, linestyle = :dash, linewidth = 2,
            label = "√3")
    axislegend(ax)
end

display(fig_res)

# 5. Comparaison entre le volume calculé et le volume théorique en fonction de la résolution
println("\n5. Efficacité du calcul de volume en fonction de la résolution")
res_perf = [10, 15, 20, 25, 30, 40, 50]
times = Float64[]
interface_cells_count = Int[]
total_cells = Int[]
vol_errors = Float64[]

for res in res_perf
    mesh_3d = create_mesh_3d(res, res, res)
    total_cell_count = res^3
    
    t_start = time()
    V, _, interface_area, _ = integrate(Tuple{0}, levelset_sphere, mesh_3d, T, zero)
    t_elapsed = time() - t_start
    
    push!(times, t_elapsed)
    push!(total_cells, total_cell_count)
    
    interface_cells = findall(interface_area .> 0)
    push!(interface_cells_count, length(interface_cells))
    
    numerical_vol = sum(V)
    rel_error = abs(numerical_vol - analytical_volume) / analytical_volume
    push!(vol_errors, rel_error)
    
    println(@sprintf("  Résolution %2d³: %.2fs, Interface: %d cellules (%.2f%%), Erreur vol: %.4f%%", 
                   res, t_elapsed, length(interface_cells),
                   100*length(interface_cells)/total_cell_count, 100*rel_error))
end

fig_perf = Figure(size = (1200, 400))
ax1 = Axis(fig_perf[1, 1], 
          title = "Temps de calcul",
          xlabel = "Taille totale du maillage", ylabel = "Temps (s)")
scatter!(ax1, total_cells, times)
lines!(ax1, total_cells, times)

ax2 = Axis(fig_perf[1, 2],
          xscale = log10, yscale = log10,
          title = "Erreur du volume vs. résolution",
          xlabel = "Résolution", ylabel = "Erreur relative")
scatter!(ax2, res_perf, vol_errors)
lines!(ax2, res_perf, vol_errors)

# Référence pour la convergence quadratique
ref = vol_errors[3] * (res_perf[3] ./ res_perf).^2
lines!(ax2, res_perf, ref, linestyle = :dash, color = :red, label = "O(h²)")
axislegend(ax2)

ax3 = Axis(fig_perf[1, 3],
          title = "Nombre de cellules d'interface",
          xlabel = "Résolution", ylabel = "Cellules d'interface")
scatter!(ax3, res_perf, interface_cells_count)
lines!(ax3, res_perf, interface_cells_count)

# Référence pour croissance en O(n²)
ref_growth = interface_cells_count[3] * (res_perf ./ res_perf[3]).^2
lines!(ax3, res_perf, ref_growth, linestyle = :dash, color = :red, label = "O(n²)")
axislegend(ax3)

display(fig_perf)

println("\nRésumé des résultats pour la sphère:")
println("  - Ordre de convergence pour le volume: approximativement O(h²)")
println("  - Ordre de convergence pour la surface: approximativement O(h)")
println("  - Croissance du nombre de cellules d'interface: approximativement O(n²)")
println("  - Croissance du temps de calcul: approximativement O(n³)")