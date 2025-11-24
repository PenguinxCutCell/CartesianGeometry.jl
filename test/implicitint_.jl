using ImplicitIntegration
using CartesianGeometry
using CairoMakie
using StaticArrays
const T=Float64

nx, ny = 20, 20
Lx, Ly = 2.0, 2.0
x0, y0 = 0.0, 0.0
dx, dy = Lx/nx, Ly/ny

x,y = collect(0:Lx/nx:Lx), collect(0:Ly/ny:Ly)
mesh = (x,y)

P = (x) -> sqrt((x[1]-1.01)^2 + (x[2]-1.01)^2) - 0.5
Φ(x,y,_=0) = P([x,y])

# VOFI reference
V, bary, interface_length, cell_type = CartesianGeometry.integrate(Tuple{0}, Φ, mesh, T, zero)
As = CartesianGeometry.integrate(Tuple{1}, Φ, mesh, T, zero)
Ws = CartesianGeometry.integrate(Tuple{0}, Φ, mesh, T, zero, bary)
Bs = CartesianGeometry.integrate(Tuple{1}, Φ, mesh, T, zero, bary)

# Utilities
isfull(val, full_val) = isapprox(val, full_val; atol=1e-8)
isempty(val) = isapprox(val, 0.0; atol=1e-8)

# Parameters
tol = 1e-8
function implicit_integration_nd(mesh, P; tol=1e-8)
    # Déterminer la dimension du problème
    dim = length(mesh)
    
    # Nombre de cellules dans chaque dimension
    n_dims = [length(mesh[d])-1 for d in 1:dim]
    
    # Nombre total de cellules
    total_cells = prod(n_dims .+ 1)  # Ajout du +1 pour corriger les dimensions
    
    # Initialiser les tableaux pour stocker les résultats
    volumes = zeros(total_cells)
    cell_types = zeros(Int, total_cells)
    interface_areas = zeros(total_cells)
    interface_centroids = [SVector{dim,Float64}(zeros(dim)) for _ in 1:total_cells]
    cell_centroids = [SVector{dim,Float64}(zeros(dim)) for _ in 1:total_cells]
    
    # Fonction d'indexation linéaire pour le stockage
    function linear_index(indices)
        idx = indices[1]
        stride = 1
        for d in 2:dim
            stride *= (n_dims[d-1] + 1)
            idx += (indices[d] - 1) * stride
        end
        return idx
    end
    
    # Fonction récursive pour parcourir les cellules en N dimensions
    function process_cell(indices, current_dim)
        if current_dim > dim
            # Tous les indices sont définis, traiter la cellule
            idx = linear_index(indices)
            
            # Points de la cellule (min et max en ND)
            a = ntuple(d -> mesh[d][indices[d]], dim)
            b = ntuple(d -> mesh[d][indices[d] + 1], dim)
            
            # Calculer le volume de la cellule
            vol_result = ImplicitIntegration.integrate(x -> 1.0, P, a, b; tol=tol)
            volumes[idx] = vol_result.val
            
            # Volume complet d'une cellule
            full_cell_volume = prod(mesh[d][indices[d] + 1] - mesh[d][indices[d]] for d in 1:dim)
            
            # Déterminer le type de cellule
            if isfull(volumes[idx], full_cell_volume)
                cell_types[idx] = 1
                
                # Barycentre de la cellule pleine
                moments = ntuple(d -> ImplicitIntegration.integrate(x -> x[d], P, a, b; tol=tol).val, dim)
                cell_centroids[idx] = SVector{dim,Float64}(moments ./ volumes[idx])
            elseif isempty(volumes[idx])
                cell_types[idx] = 0
                
                # Centre géométrique pour cellule vide
                center = ntuple(d -> 0.5*(mesh[d][indices[d]] + mesh[d][indices[d] + 1]), dim)
                cell_centroids[idx] = SVector{dim,Float64}(center)
            else
                cell_types[idx] = -1
                
                # Barycentre de la cellule coupée
                moments = ntuple(d -> ImplicitIntegration.integrate(x -> x[d], P, a, b; tol=tol).val, dim)
                cell_centroids[idx] = SVector{dim,Float64}(moments ./ volumes[idx])
            end
            
            # Calculer l'aire d'interface
            interface_result = ImplicitIntegration.integrate(x -> 1.0, P, a, b; surface=true, tol=tol)
            interface_areas[idx] = interface_result.val
            
            # Calculer le barycentre de l'interface
            if !isempty(interface_areas[idx])
                interface_moments = ntuple(d -> ImplicitIntegration.integrate(x -> x[d], P, a, b; surface=true, tol=tol).val, dim)
                interface_centroids[idx] = SVector{dim,Float64}(interface_moments ./ interface_areas[idx])
            else
                interface_centroids[idx] = SVector{dim,Float64}(zeros(dim))
            end
            
            return
        end
        
        # Parcourir les indices pour la dimension courante
        for i in 1:n_dims[current_dim]
            indices[current_dim] = i
            process_cell(indices, current_dim + 1)
        end
    end
    
    # Lancer le traitement récursif avec des indices initiaux
    indices = ones(Int, dim)
    process_cell(indices, 1)
    
    return volumes, cell_types, interface_areas, interface_centroids, cell_centroids
end

volumes, cell_types, interface_areas, interface_centroids, cell_centroids = implicit_integration_nd(mesh, P; tol=tol)

linear_idx(i,j) = i + (j-1)*(nx+1)
Ax = zeros((nx+1)*(ny+1))
Ay = zeros((nx+1)*(ny+1))
As_i = (Ax,Ay)
# For Ax (line x fixed, integrate over y):
for i in 1:nx+1
    xi = x[i]
    for j in 1:ny
        idx = linear_idx(i, j)
        yj = y[j]
        yj1 = y[j+1]
        P_1d = (y) -> P((xi,y[1]))
        Ax_result = ImplicitIntegration.integrate((x)->1.0, P_1d, (yj,), (yj1,); tol=tol)
        As_i[1][idx] = Ax_result.val
    end
end

# For Ay (line y fixed, integrate over x):
for i in 1:nx
    xi = x[i]
    xi1 = x[i+1]
    for j in 1:ny+1
        idx = linear_idx(i, j)

        yj = y[j]
        P_1d = (x) -> P((x[1],yj))
        Ay_result = ImplicitIntegration.integrate((x)->1.0, P_1d, (xi,), (xi1,); tol=tol)
        As_i[2][idx] = Ay_result.val
    end
end

# Compute the Second Capacity geometric information using ImplicitIntegration
Ws_i = (zeros((nx+1)*(ny+1)), zeros((nx+1)*(ny+1)))
Bs_i = (zeros((nx+1)*(ny+1)), zeros((nx+1)*(ny+1)))

# For Bx (line x_centroid fixed, integrate over y):
for i in 1:nx
    for j in 1:ny
        idx = linear_idx(i, j)
        xi = cell_centroids[idx][1]
        yj = y[j]
        yj1 = y[j+1]
        P_1d = (y) -> P((xi,y[1]))
        Bx_result = ImplicitIntegration.integrate((x)->1.0, P_1d, (yj,), (yj1,); tol=tol)
        Bs_i[1][idx] = Bx_result.val
    end
end

# For By (line y_centroid fixed, integrate over x):
for i in 1:nx
    for j in 1:ny
        idx = linear_idx(i, j)
        xi = x[i]
        xi1 = x[i+1]
        yj = cell_centroids[idx][2]
        P_1d = (x) -> P((x[1],yj))
        By_result = ImplicitIntegration.integrate((x)->1.0, P_1d, (xi,), (xi1,); tol=tol)
        Bs_i[2][idx] = By_result.val
    end
end

# For Wx (staggered volume based on x_centroid)
for i in 1:nx
    for j in 1:ny
        idx = linear_idx(i, j)
        a = (cell_centroids[linear_idx(max(i,1), 1)][1], y[j])
        b = (cell_centroids[linear_idx(min(i+1,nx), 1)][1], y[j+1])
        Wx_result = ImplicitIntegration.integrate((x)->1.0, P, a,b; tol=tol)
        Ws_i[1][idx] = Wx_result.val

    end
end

# For Wy (staggered volume based on y_centroid)
for i in 1:nx
    for j in 1:ny
        idx = linear_idx(i, j)
        a = (x[i], cell_centroids[linear_idx(1, max(j,1))][2])
        b = (x[i+1], cell_centroids[linear_idx(1, min(j+1,ny))][2])
        Wy_result = ImplicitIntegration.integrate((x)->1.0, P, a,b; tol=tol)
        Ws_i[2][idx] = Wy_result.val
    end
end

