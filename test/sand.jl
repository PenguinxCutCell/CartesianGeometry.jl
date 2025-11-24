using ImplicitIntegration
using CartesianGeometry
using StaticArrays
const T=Float64
using Plots
Plots.default(show=true)

# Parameters
nx, ny = 80, 80
x0, y0 = 0.0, 0.0
Lx, Ly = 1.0, 1.0
Δx, Δy = Lx/nx, Ly/ny

mesh = (collect(x0:Δx:Lx), collect(y0:Δy:Ly))
mesh_center = (0.5Δx:Δx:Lx-0.5Δx, 0.5Δy:Δy:Ly-0.5Δy)
println("Mesh Size: ", length(mesh[1]), " x ", length(mesh[2]))
println("Mesh Center Size: ", length(mesh_center[1]), " x ", length(mesh_center[2]))
Φ(x,y,_=0) = y - 0.5 - 0.1*sin(3π*x)

V, bary, interface_length, cell_type = CartesianGeometry.integrate(Tuple{0}, Φ, mesh, T, zero)

# Plot
V = reshape(V, nx+1, ny+1)'[1:end-1,1:end-1]
heatmap(mesh[1], mesh[2], V, aspect_ratio=1, c=:viridis, color=:grays, grid=false, cbar=true)
readline()

# Calculate Height in each column
H = sum(V, dims=1)

println("Size of Height: ", size(H))
println("Height: ", H)

plot(mesh_center[1],H[:],seriestype=:scatter)
readline()

# Compute the interface position
s = y0 .+ H ./ Δx

println("Size of s: ", size(s))
println("s: ", s)

plot(mesh_center[1],s[:],seriestype=:scatter)
readline()

# Interpolate using Interpolations.jl from cell center : No volume conservation
using Interpolations
itp_lin = linear_interpolation((mesh_center[1]), vec(s), extrapolation_bc=Line())
itp_cub = cubic_spline_interpolation((mesh_center[1]), vec(s), extrapolation_bc=Line())

plot!(mesh[1], itp(mesh[1]), c=:red, lw=2)
plot!(mesh[1], itp_cub(mesh[1]), c=:green, lw=2)
readline()

Φ_l(x,y,_=0) = y - itp(x)
Φ_c(x,y,_=0) = y - itp_cub(x)

V_l, bary, interface_length, cell_type = CartesianGeometry.integrate(Tuple{0}, Φ_l, mesh, T, zero)
V_c, bary, interface_length, cell_type = CartesianGeometry.integrate(Tuple{0}, Φ_c, mesh, T, zero)

# Plot
V_l = reshape(V_l, nx+1, ny+1)'[1:end-1,1:end-1]
V_c = reshape(V_c, nx+1, ny+1)'[1:end-1,1:end-1]

p1 = heatmap(mesh[1], mesh[2], V_l, aspect_ratio=1, c=:viridis, color=:grays, grid=false, cbar=true, title="Linear Interpolation")
p2 = heatmap(mesh[1], mesh[2], V_c, aspect_ratio=1, c=:viridis, color=:grays, grid=false, cbar=true, title="Cubic Interpolation")

plot(p1, p2, layout=(1,2), size=(800,400))

readline()

# Vérifier la conservation du volume
new_H_l = sum(V_l, dims=1)
new_H_c = sum(V_c, dims=1)  
relative_error_l = maximum(abs.(new_H_l .- H) ./ H)
relative_error_c = maximum(abs.(new_H_c .- H) ./ H)
println("Erreur relative maximale no volume conservation (Interpolation Linéaire): ", relative_error_l)
println("Erreur relative maximale no volume conservation (Interpolation Cubique): ", relative_error_c)

# Interpolate with volume conservation
function height_function_volume_preserving(mesh_center, H, Δx)
    nx = length(mesh_center[1])
    
    # Système de 2*nx équations linéaires
    # Variables: a_i et b_i pour chaque cellule i
    A = zeros(2*nx, 2*nx)
    b = zeros(2*nx)
    
    # Première série d'équations: conservation du volume
    # a_i + 0.5*Δx*b_i = H_i
    for i in 1:nx
        A[i, i] = 1.0              # Coefficient de a_i
        A[i, nx+i] = 0.5*Δx        # Coefficient de b_i
        b[i] = H[i]                # Terme de droite: H_i
    end
    
    # Deuxième série d'équations: continuité aux interfaces
    # a_i + Δx*b_i = a_{i+1}
    for i in 1:nx
        A[nx+i, i] = 1.0           # Coefficient de a_i
        A[nx+i, nx+i] = Δx         # Coefficient de b_i
        # Pour assurer la périodicité
        next_i = i % nx + 1
        A[nx+i, next_i] = -1.0     # Coefficient de a_{i+1}
    end

    # BC pour la périodicité : a_1 = a_n+1
    A[2*nx, 1] = 1.0
    A[2*nx, nx] = -1.0


    
    # Résolution du système linéaire
    coeffs = A\b
    
    # Extraction des coefficients
    a_coeffs = coeffs[1:nx]
    b_coeffs = coeffs[nx+1:2*nx]
    
    # Fonction d'interpolation
    function interp(x)
        # Trouver la cellule contenant x
        # Gérer la périodicité
        x_periodic = mod(x, nx*Δx)
        i = min(floor(Int, x_periodic/Δx) + 1, nx)
        
        # Position relative dans la cellule
        x_i = (i-1)*Δx
        s = x_periodic - x_i
        
        # Évaluer le polynôme P_i(x) = a_i + b_i*(x - x_i)
        return a_coeffs[i] + b_coeffs[i]*s
    end
    
    return interp, a_coeffs, b_coeffs
end

# Après avoir calculé H (les hauteurs)
# Remplacer l'interpolation actuelle par celle préservant le volume
s = y0 .+ H ./ Δx
interp_fn, a_coeffs, b_coeffs = height_function_volume_preserving(mesh_center, vec(s), Δx)

# Visualiser l'interpolation
plot(mesh_center[1], vec(s), seriestype=:scatter, label="Points de hauteur")
x_fine = range(0, Lx, length=500)
plot!(mesh_center[1], interp_fn.(mesh_center[1]), c=:blue, lw=2, label="Interpolation Volume Preserving")
plot!(x_fine, interp_fn.(x_fine), c=:blue, lw=2, label="Interpolation Volume Preserving")
readline()

# Utiliser cette fonction pour définir l'interface
Φ(x,y,_=0) = y - interp_fn(x)

# Calculer le volume avec la nouvelle définition d'interface
V, bary, interface_length, cell_type = CartesianGeometry.integrate(Tuple{0}, Φ, mesh, T, zero)

# Afficher le résultat
V = reshape(V, nx+1, ny+1)'[1:end-1,1:end-1]
heatmap(mesh[1], mesh[2], V, aspect_ratio=1, c=:viridis, color=:grays, 
        grid=false, cbar=true, title="Volume Preserving Height Function")
readline()

# Vérifier la conservation du volume
new_H = sum(V, dims=1)
relative_error = maximum(abs.(new_H .- H) ./ H)
println("Erreur relative maximale volume conserving (Linear Interpolation): ", relative_error)

# Interpolation quadratique qui conserve le volume
function height_function_quadratic_volume_preserving(mesh_center, H, Δx)
    nx = length(mesh_center[1])
    
    # Système de 3*nx équations linéaires
    # Variables: a_i, b_i, c_i pour chaque cellule i
    A = zeros(3*nx, 3*nx)
    b = zeros(3*nx)
    
    for i = 1:nx
        # Indice de la première équation pour la cellule i
        row = 3*(i-1) + 1
        
        # Colonne de base pour les coefficients a_i, b_i, c_i
        col = 3*(i-1) + 1
        
        # Équation 1: Conservation du volume (V_i)
        # a_i + (Δx/2)*b_i + (Δx²/3)*c_i = H_i
        A[row, col] = 1.0                    # Coefficient de a_i
        A[row, col+1] = Δx/2                 # Coefficient de b_i
        A[row, col+2] = (Δx^2)/3             # Coefficient de c_i
        b[row] = H[i]                        # RHS = H_i
        
        # Équation 2: Continuité des valeurs (C1_i)
        # a_i + Δx*b_i + Δx²*c_i - a_{i+1} = 0
        A[row+1, col] = 1.0                  # Coefficient de a_i
        A[row+1, col+1] = Δx                 # Coefficient de b_i
        A[row+1, col+2] = Δx^2               # Coefficient de c_i
        
        # Périodicité pour a_{i+1}
        next_i = i % nx + 1
        next_col = 3*(next_i-1) + 1
        A[row+1, next_col] = -1.0            # Coefficient de a_{i+1}
        
        # Équation 3: Continuité des dérivées (C2_i)
        # b_i + 2*Δx*c_i - b_{i+1} = 0
        A[row+2, col+1] = 1.0                # Coefficient de b_i
        A[row+2, col+2] = 2*Δx               # Coefficient de c_i
        A[row+2, next_col+1] = -1.0          # Coefficient de b_{i+1}
    end
    
    # Résolution du système linéaire
    coeffs = A\b
    
    # Extraction des coefficients
    a_coeffs = coeffs[1:3:end]
    b_coeffs = coeffs[2:3:end]
    c_coeffs = coeffs[3:3:end]
    
    # Fonction d'interpolation avec extrapolation linéaire
    function interp(x)
        # Déterminer les bornes du domaine
        x_min = 0.0
        x_max = nx * Δx
        
        # Si on est à l'intérieur du domaine, interpolation normale avec périodicité
        if x_min <= x <= x_max
            # Gérer la périodicité
            x_periodic = mod(x, nx*Δx)
            i = min(floor(Int, x_periodic/Δx) + 1, nx)
            
            # Position relative dans la cellule
            x_i = (i-1)*Δx
            s = x_periodic - x_i
            
            # Évaluer le polynôme quadratique P_i(x) = a_i + b_i*(x - x_i) + c_i*(x - x_i)^2
            return a_coeffs[i] + b_coeffs[i]*s + c_coeffs[i]*s^2
        else
            # Extrapolation linéaire
            if x < x_min
                # Extrapolation à gauche en utilisant la valeur et la dérivée au bord gauche
                edge_value = a_coeffs[1]  # Valeur à x=0
                edge_slope = b_coeffs[1]  # Dérivée à x=0
                return edge_value + edge_slope * (x - x_min)
            else  # x > x_max
                # Extrapolation à droite en utilisant la valeur et la dérivée au bord droit
                # Valeur au bord x=x_max (fin de la dernière cellule)
                edge_value = a_coeffs[nx] + b_coeffs[nx]*Δx + c_coeffs[nx]*(Δx^2)
                # Dérivée au bord x=x_max
                edge_slope = b_coeffs[nx] + 2*c_coeffs[nx]*Δx
                return edge_value + edge_slope * (x - x_max)
            end
        end
    end
    
    # Fonction pour calculer la dérivée (avec extrapolation constante de la dérivée)
    function deriv(x)
        # Déterminer les bornes du domaine
        x_min = 0.0
        x_max = nx * Δx
        
        # Si on est à l'intérieur du domaine, dérivée normale
        if x_min <= x <= x_max
            # Gérer la périodicité
            x_periodic = mod(x, nx*Δx)
            i = min(floor(Int, x_periodic/Δx) + 1, nx)
            
            # Position relative dans la cellule
            x_i = (i-1)*Δx
            s = x_periodic - x_i
            
            # Dérivée du polynôme: P_i'(x) = b_i + 2*c_i*(x - x_i)
            return b_coeffs[i] + 2*c_coeffs[i]*s
        else
            # Extrapolation constante de la dérivée
            if x < x_min
                # Dérivée constante à gauche = dérivée au bord gauche
                return b_coeffs[1]
            else  # x > x_max
                # Dérivée constante à droite = dérivée au bord droit
                return b_coeffs[nx] + 2*c_coeffs[nx]*Δx
            end
        end
    end
    
    return interp, deriv, a_coeffs, b_coeffs, c_coeffs
end

# Appliquer l'interpolation quadratique qui préserve le volume
interp_fn, deriv_fn, a_coeffs, b_coeffs, c_coeffs = 
    height_function_quadratic_volume_preserving(mesh_center, vec(s), Δx)

# Visualiser l'interpolation et sa dérivée
p1 = plot(mesh_center[1], vec(s), seriestype=:scatter, label="Points de hauteur")
x_fine = range(0, Lx, length=500)
plot!(p1, x_fine, [interp_fn(x) for x in x_fine], c=:blue, lw=2, 
      label="Interpolation quadratique C¹")

# Visualiser la dérivée
p2 = plot(x_fine, [deriv_fn(x) for x in x_fine], c=:red, lw=2, 
          label="Dérivée de l'interpolation", legend=:topleft)

# Afficher les deux graphiques
plot(p1, p2, layout=(2,1), size=(800,600))
readline()

# Utiliser cette fonction pour définir l'interface
Φ(x,y,_=0) = y - interp_fn(x)

# Calculer le volume avec la nouvelle définition d'interface
V, bary, interface_length, cell_type = CartesianGeometry.integrate(Tuple{0}, Φ, mesh, T, zero)

# Afficher le résultat
V = reshape(V, nx+1, ny+1)'[1:end-1,1:end-1]
heatmap(mesh[1], mesh[2], V, aspect_ratio=1, c=:viridis, color=:grays, 
        grid=false, cbar=true, title="Volume Preserving C¹ Quadratic Interpolation")
readline()

# Vérifier la conservation du volume
new_H = sum(V, dims=1)
relative_error = maximum(abs.(new_H .- H) ./ H)
println("Erreur relative maximale volume conserving (Quadratic Interpolation): ", relative_error)

